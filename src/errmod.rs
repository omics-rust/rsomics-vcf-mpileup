// MAQ error model — htslib errmod.c (bcftools 1.21 / htslib 1.21, MIT-licensed).
//
// Precomputed tables allocated once, reused per pileup column:
//   fk[n]       = pow(1-depcorr, n) * (1-eta) + eta
//   beta[q,n,k] = Phred cost for k matches out of n at error-rate 10^(-q/10)
//   lhet[n,k]   = log P(k minor allele out of n in a het genotype)

pub(crate) struct ErrMod {
    pub(crate) fk: Vec<f64>,
    /// beta[q<<16 | n<<8 | k]
    pub(crate) beta: Vec<f64>,
    /// lhet[n<<8 | k]
    pub(crate) lhet: Vec<f64>,
}

impl ErrMod {
    pub(crate) fn new() -> Self {
        const DEPCORR: f64 = 0.17;
        const ETA: f64 = 0.03;

        let mut fk = vec![0f64; 256];
        fk[0] = 1.0;
        for (n, slot) in fk.iter_mut().enumerate().skip(1) {
            *slot = (1.0 - DEPCORR).powi(n as i32) * (1.0 - ETA) + ETA;
        }

        let n_size: usize = 256;
        let mut lc = vec![0f64; n_size * n_size];
        for n in 1..n_size {
            let lf_n = lgamma(n as f64 + 1.0);
            for k in 1..=n {
                lc[n << 8 | k] = lf_n - lgamma(k as f64 + 1.0) - lgamma((n - k) as f64 + 1.0);
            }
        }

        let mut beta = vec![0f64; 64 * 256 * 256];
        for q in 1usize..64 {
            let e = 10f64.powf(-(q as f64) / 10.0);
            let le = e.ln();
            let le1 = (1.0 - e).ln();
            for n in 1usize..=255 {
                let base_idx = q << 16 | n << 8;
                let mut sum1 = lc[n << 8 | n] + n as f64 * le;
                beta[base_idx | n] = f64::INFINITY;
                let mut k = n as isize - 1;
                while k >= 0 {
                    let uk = k as usize;
                    let prev = sum1;
                    let delta = lc[n << 8 | uk] + uk as f64 * le + (n - uk) as f64 * le1 - prev;
                    sum1 = prev + delta.exp().ln_1p();
                    // -10/ln(10) converts ln-scale to Phred.
                    beta[base_idx | uk] = -4.342_944_819_032_518 * (prev - sum1);
                    k -= 1;
                }
            }
        }

        let mut lhet = vec![0f64; 256 * 256];
        for n in 0..256usize {
            for k in 0..256usize {
                lhet[n << 8 | k] = lc[n << 8 | k] - std::f64::consts::LN_2 * n as f64;
            }
        }

        Self { fk, beta, lhet }
    }

    /// `errmod_cal`: fill `q_out[j*m+k]` with phred-scale likelihood of diploid genotype (j,k).
    ///
    /// `bases_in`: packed `q<<5 | strand<<4 | base` (4-bit base index: 0=A 1=C 2=G 3=T 4=N).
    pub(crate) fn cal(&self, bases_in: &mut [u16], m: usize, q_out: &mut [f32]) {
        let n_orig = bases_in.len();
        q_out.iter_mut().for_each(|v| *v = 0.0);
        if n_orig == 0 {
            return;
        }
        let n = n_orig.min(255);
        bases_in[..n].sort_unstable();

        let mut bsum = [0f64; 16];
        let mut c = [0usize; 16];
        let mut w = [0usize; 32];

        for j in (0..n).rev() {
            let b = bases_in[j];
            let q = {
                let q_raw = (b >> 5) as usize;
                q_raw.clamp(4, 63)
            };
            let basestrand = (b & 0x1f) as usize;
            let base = (b & 0x0f) as usize;
            bsum[base] += self.fk[w[basestrand]] * self.beta[q << 16 | n << 8 | c[base]];
            c[base] += 1;
            w[basestrand] += 1;
        }

        for j in 0..m {
            let (mut tmp1, mut tmp2) = (0f64, 0usize);
            for k in 0..m {
                if k == j {
                    continue;
                }
                tmp1 += bsum[k];
                tmp2 += c[k];
            }
            if tmp2 > 0 {
                q_out[j * m + j] = tmp1 as f32;
            }
            for k in (j + 1)..m {
                let cjk = c[j] + c[k];
                let (mut t1, mut t2) = (0f64, 0usize);
                for i in 0..m {
                    if i == j || i == k {
                        continue;
                    }
                    t1 += bsum[i];
                    t2 += c[i];
                }
                let lh_idx = (cjk << 8 | c[k]).min(256 * 256 - 1);
                let v =
                    (-4.342_944_819_032_518 * self.lhet[lh_idx]) + if t2 > 0 { t1 } else { 0.0 };
                q_out[j * m + k] = v as f32;
                q_out[k * m + j] = v as f32;
            }
            for k in 0..m {
                if q_out[j * m + k] < 0.0 {
                    q_out[j * m + k] = 0.0;
                }
            }
        }
    }
}

/// Lanczos lgamma, ~15-digit accuracy for x > 0 (Numerical Recipes §6.1).
pub(crate) fn lgamma(x: f64) -> f64 {
    const C: [f64; 6] = [
        76.180_091_729_471_46,
        -86.505_320_329_416_77,
        24.014_098_240_830_91,
        -1.231_739_572_450_155,
        1.208_650_973_866_179e-3,
        -5.395_239_384_953_e-6,
    ];
    let x = x - 1.0;
    let y = x + 5.5;
    let tmp = (x + 0.5) * y.ln() - y;
    let mut ser = 1.000_000_000_190_015;
    let mut xx = x + 1.0;
    for ci in &C {
        ser += ci / xx;
        xx += 1.0;
    }
    tmp + (2.506_628_274_631_001 * ser).ln()
}
