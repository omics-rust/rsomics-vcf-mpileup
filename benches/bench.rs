use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_vcf_mpileup(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-vcf-mpileup");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/small.bam");
    let fa = manifest.join("tests/golden/small.fa");
    c.bench_function("rsomics-vcf-mpileup golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([bam.to_str().unwrap(), "-f", fa.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_vcf_mpileup);
criterion_main!(benches);
