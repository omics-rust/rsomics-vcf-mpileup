// Per-position GL/PL computation (bcftools bam2bcf.c SNP path).

use rsomics_pileup::Column;

use crate::bases::NT16_TO_INT;
use crate::errmod::ErrMod;

// Defaults matching bcftools mpileup.
pub(crate) const MAX_BASEQ: u8 = 60;
pub(crate) const CAP_Q: u8 = 60;
pub(crate) const DEF_MAPQ: u8 = 20;

pub(crate) struct GlResult {
    pub(crate) ref_base: u8, // ACGTN index (0-4)
    pub(crate) alt_base: u8, // ACGTN index; only meaningful when has_alt
    pub(crate) has_alt: bool,
    pub(crate) ad: Vec<u32>,
    pub(crate) dp: u32,
    pub(crate) pl: Vec<u32>, // PL[RR], PL[RA], PL[AA] or PL[RR] only
}

pub(crate) fn process_column(
    col: &Column,
    ref_nt16: u8,
    min_baseq: u8,
    max_depth: u32,
    errmod: &ErrMod,
    bases_buf: &mut Vec<u16>,
    gl_buf: &mut Vec<f32>,
) -> GlResult {
    let ref_int = NT16_TO_INT[(ref_nt16 & 0xf) as usize];

    let mut counts = [0u32; 5];
    bases_buf.clear();

    let mut seen = 0u32;
    for (read, rec) in col.reads.iter().zip(col.records.iter()) {
        if read.is_del || read.is_refskip {
            continue;
        }
        if seen >= max_depth {
            break;
        }
        let qual_scores = rec.quality_scores();
        let qpos = read.qpos;

        let bq_raw = if qpos < qual_scores.len() {
            qual_scores[qpos]
        } else {
            0
        };
        if bq_raw < min_baseq {
            continue;
        }

        let nib = rec.seq_nibble(qpos);
        let base_int = NT16_TO_INT[(nib & 0xf) as usize];

        let mapq: u8 = {
            let mq = rec.mapping_quality();
            if mq == 255 { DEF_MAPQ } else { mq }
        };

        // Effective quality: min(bq, max_baseQ, mapQ, capQ, 63), floor 4.
        let q = bq_raw.min(MAX_BASEQ).min(mapq).min(CAP_Q).clamp(4, 63);

        let strand: u8 = u8::from(rec.flags() & 0x10 != 0);
        bases_buf.push((q as u16) << 5 | (strand as u16) << 4 | (base_int as u16 & 0xf));
        counts[base_int as usize] += 1;
        seen += 1;
    }

    let dp = bases_buf.len() as u32;
    let m = 5;

    gl_buf.resize(m * m, 0.0);
    errmod.cal(bases_buf, m, gl_buf);

    // Best ALT by count, excluding REF and N.
    let mut best_alt = 5usize;
    let mut best_alt_count = 0u32;
    for (i, &cnt) in counts[..4].iter().enumerate() {
        if ref_int < 4 && i == ref_int as usize {
            continue;
        }
        if cnt > best_alt_count {
            best_alt_count = cnt;
            best_alt = i;
        }
    }
    let has_alt = best_alt < 5 && best_alt_count > 0;

    let ref_ad = if ref_int < 4 {
        counts[ref_int as usize]
    } else {
        0
    };
    let alt_ad = if has_alt { counts[best_alt] } else { 0 };
    let ad = if has_alt {
        vec![ref_ad, alt_ad]
    } else {
        vec![ref_ad]
    };

    let pl = if dp == 0 || !has_alt {
        vec![0u32; if has_alt { 3 } else { 1 }]
    } else {
        let r = if ref_int < 4 { ref_int as usize } else { 4 };
        let a = best_alt;
        let pl_rr = gl_buf[r * m + r];
        let pl_ra = gl_buf[r * m + a];
        let pl_aa = gl_buf[a * m + a];
        let min_pl = pl_rr.min(pl_ra).min(pl_aa);
        vec![
            (pl_rr - min_pl).round() as u32,
            (pl_ra - min_pl).round() as u32,
            (pl_aa - min_pl).round() as u32,
        ]
    };

    GlResult {
        ref_base: ref_int,
        alt_base: if has_alt { best_alt as u8 } else { 4 },
        has_alt,
        ad,
        dp,
        pl,
    }
}
