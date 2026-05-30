// 4-bit BAM base encoding and FASTA NT16 helpers.

/// BAM seq_nt16_int: 4-bit NT16 code → ACGTN index (N=4).
pub(crate) const NT16_TO_INT: [u8; 16] = [4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4];

/// ACGTN index → uppercase ASCII.
pub(crate) const INT_TO_BASE: [u8; 5] = *b"ACGTN";

/// ASCII FASTA base → 4-bit NT16 (htslib seq_nt16_table).
pub(crate) fn nt16_of_ascii(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 15,
    }
}
