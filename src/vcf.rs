// VCF header and per-record emit.

use std::io::Write;

pub(crate) fn write_header(
    w: &mut impl Write,
    sq: &[(String, u64)],
    sample: &str,
) -> std::io::Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    for (name, len) in sq {
        writeln!(w, "##contig=<ID={name},length={len}>")?;
    }
    writeln!(
        w,
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for ref and alt alleles\">"
    )?;
    writeln!(
        w,
        "##source=rsomics-vcf-mpileup {}",
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(
        w,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}"
    )?;
    Ok(())
}
