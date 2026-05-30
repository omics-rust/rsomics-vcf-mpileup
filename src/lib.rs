//! bcftools mpileup SNP mode — per-position genotype-likelihood VCF from a
//! single coordinate-sorted BAM with a required reference.
//!
//! Implements the revised MAQ error model (`errmod_cal` / `htslib errmod.c`)
//! to compute per-genotype Phred-scaled likelihoods (PL), together with the
//! AD (allelic depth) and DP FORMAT fields that bcftools mpileup emits by
//! default for single-sample SNP positions.
//!
//! PL vector ordering follows VCF 4.2 / bcftools convention: for a diploid
//! with `n` alleles, entries are lexicographic (i,j) pairs where i≤j —
//! for biallelic (REF,ALT): PL[RR], PL[RA], PL[AA].
//!
//! Source references (MIT-licensed):
//!   bcftools 1.21 `mpileup.c`, `bam2bcf.c`
//!   htslib 1.21 `errmod.c` (MAQ error model)

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_precision_loss
)]

mod bases;
mod errmod;
mod gl;
mod vcf;

use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::Path;

use noodles::fasta;
use rsomics_bamio::raw::{RawRecord, read_record};
use rsomics_common::{Result, RsomicsError};
use rsomics_pileup::{Column, PileupEngine, PileupOpts};

use bases::{INT_TO_BASE, nt16_of_ascii};
use errmod::ErrMod;
use gl::process_column;
use vcf::write_header;

/// Default minimum mapping quality (`-q`).
pub const DEFAULT_MIN_MAPQ: u8 = 0;
/// Default minimum base quality (`-Q`).
pub const DEFAULT_MIN_BASEQ: u8 = 1;
/// Default maximum depth (bcftools `max_depth = 250`).
pub const DEFAULT_MAX_DEPTH: u32 = 250;

/// FLAG bits filtered by default: UNMAP|SECONDARY|QCFAIL|DUP.
const DEFAULT_RFLAG_FILTER: u16 = 0x004 | 0x100 | 0x200 | 0x400;

#[derive(Debug, Clone)]
pub struct VcfMpileupOpts {
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub max_depth: u32,
    pub region: Option<String>,
    pub threads: NonZeroUsize,
}

impl Default for VcfMpileupOpts {
    fn default() -> Self {
        Self {
            min_mapq: DEFAULT_MIN_MAPQ,
            min_baseq: DEFAULT_MIN_BASEQ,
            max_depth: DEFAULT_MAX_DEPTH,
            region: None,
            threads: NonZeroUsize::new(1).unwrap(),
        }
    }
}

/// Run vcf-mpileup, writing VCF (uncompressed) to `out`.
///
/// `bam`: coordinate-sorted BAM (indexed only when `opts.region` is set).
/// `fasta_ref`: faidx-indexed FASTA. Required.
pub fn run(bam: &Path, fasta_ref: &Path, out: impl Write, opts: &VcfMpileupOpts) -> Result<()> {
    let ref_seqs = load_reference(fasta_ref)?;

    let mut reader = rsomics_bamio::open_with_workers(bam, opts.threads)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    let sq: Vec<(String, u64)> = header
        .reference_sequences()
        .iter()
        .map(|(name, map)| {
            let len: u64 = usize::from(map.length()) as u64;
            (name.to_string(), len)
        })
        .collect();

    let mut out = BufWriter::with_capacity(256 * 1024, out);
    let sample_name = bam.file_stem().and_then(|s| s.to_str()).unwrap_or("SAMPLE");
    write_header(&mut out, &sq, sample_name).map_err(RsomicsError::Io)?;

    let tid_to_name: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|n| n.to_string())
        .collect();

    // Region filter: pre-parse so we can skip columns outside the requested window.
    let region_filter: Option<(i32, i64, i64)> = opts.region.as_deref().map(|reg| {
        let (chrom, beg, end) = parse_region(reg);
        let tid = tid_to_name
            .iter()
            .position(|n| n == &chrom)
            .map(|i| i as i32)
            .unwrap_or(-1);
        (tid, beg as i64, end as i64)
    });

    let pileup_opts = PileupOpts {
        min_mapq: opts.min_mapq,
        no_orphan: true,
        rflag_filter: DEFAULT_RFLAG_FILTER,
        ..PileupOpts::default()
    };

    let errmod = ErrMod::new();
    let mut bases_buf: Vec<u16> = Vec::with_capacity(512);
    let mut gl_buf: Vec<f32> = vec![0.0; 25];
    let mut engine = PileupEngine::new(pileup_opts);

    let emit = |col: &Column| -> Result<()> {
        if let Some((req_tid, beg, end)) = region_filter
            && (col.tid != req_tid || col.pos < beg || col.pos >= end)
        {
            return Ok(());
        }

        let tid = col.tid as usize;
        let pos0 = col.pos as usize;

        let ref_byte = tid_to_name
            .get(tid)
            .and_then(|name| ref_seqs.get(name))
            .and_then(|seq| seq.get(pos0).copied())
            .unwrap_or(b'N');

        let ref_nt16 = nt16_of_ascii(ref_byte);

        let gl = process_column(
            col,
            ref_nt16,
            opts.min_baseq,
            opts.max_depth,
            &errmod,
            &mut bases_buf,
            &mut gl_buf,
        );

        if gl.dp == 0 || !gl.has_alt {
            return Ok(());
        }

        let contig_name = tid_to_name.get(tid).map_or(".", String::as_str);
        let ref_char = INT_TO_BASE[gl.ref_base.min(4) as usize];
        let alt_char = INT_TO_BASE[gl.alt_base as usize];
        let pos1 = pos0 + 1;

        // QUAL: PL[AA] is a proxy for variant quality (bcftools convention).
        let qual_str = gl
            .pl
            .get(2)
            .map_or_else(|| ".".to_owned(), |v| v.to_string());

        write!(out, "{contig_name}\t{pos1}\t.\t")?;
        out.write_all(&[ref_char]).map_err(RsomicsError::Io)?;
        write!(out, "\t")?;
        out.write_all(&[alt_char]).map_err(RsomicsError::Io)?;
        write!(out, "\t{qual_str}\t.\tDP={dp}\tPL:DP:AD\t", dp = gl.dp)?;

        let pl_str = gl
            .pl
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let ad_str = gl
            .ad
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(out, "{pl_str}:{dp}:{ad_str}", dp = gl.dp).map_err(RsomicsError::Io)
    };

    let inner = reader.get_mut();
    let mut rec = RawRecord::default();

    rsomics_pileup::run(
        &mut engine,
        || -> Result<Option<RawRecord>> {
            let n = read_record(inner, &mut rec)?;
            if n == 0 {
                Ok(None)
            } else {
                Ok(Some(rec.clone()))
            }
        },
        emit,
    )?;

    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

fn load_reference(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    let file = std::fs::File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut reader = fasta::io::Reader::new(std::io::BufReader::new(file));
    let mut map = HashMap::new();
    for result in reader.records() {
        let rec =
            result.map_err(|e| RsomicsError::InvalidInput(format!("FASTA read error: {e}")))?;
        let name = String::from_utf8_lossy(rec.name()).into_owned();
        let seq: Vec<u8> = rec.sequence().as_ref().to_vec();
        map.insert(name, seq);
    }
    Ok(map)
}

/// Parse "chrom", "chrom:start-end" (1-based inclusive) → (chrom, 0-based beg, 0-based excl end).
fn parse_region(region: &str) -> (String, u64, u64) {
    if let Some((chrom, range)) = region.split_once(':') {
        if let Some((s, e)) = range.split_once('-') {
            let beg = s
                .replace(',', "")
                .parse::<u64>()
                .unwrap_or(1)
                .saturating_sub(1);
            let end = e.replace(',', "").parse::<u64>().unwrap_or(u64::MAX);
            return (chrom.to_owned(), beg, end);
        }
        let beg = range
            .replace(',', "")
            .parse::<u64>()
            .unwrap_or(1)
            .saturating_sub(1);
        (chrom.to_owned(), beg, u64::MAX)
    } else {
        (region.to_owned(), 0, u64::MAX)
    }
}
