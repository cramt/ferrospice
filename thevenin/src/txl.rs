//! TXL (Single Lossy Transmission Line) model.
//!
//! Uses Padé approximation of the characteristic admittance Y(s) and
//! propagation function exp(-gamma*l) to obtain poles and residues for
//! efficient transient convolution.
//!
//! Reference: ngspice `src/spicelib/devices/txl/`

use thevenin_types::ModelDef;

// ---------------------------------------------------------------------------
// Padé term: pole + residue + convolution accumulator
// ---------------------------------------------------------------------------

/// A single pole-residue term for impulse response approximation.
#[derive(Debug, Clone)]
pub struct Term {
    /// Residue coefficient.
    pub c: f64,
    /// Pole location (negative real part for stability).
    pub x: f64,
    /// Accumulated convolution value at the input port.
    pub cnv_i: f64,
    /// Accumulated convolution value at the output port.
    pub cnv_o: f64,
}

impl Term {
    pub fn zero() -> Self {
        Self {
            c: 0.0,
            x: 0.0,
            cnv_i: 0.0,
            cnv_o: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// TXLine: internal transmission line state (mirrors ngspice TXLine struct)
// ---------------------------------------------------------------------------

/// Internal TXL line state: poles, residues, and convolution history.
#[derive(Debug, Clone)]
pub struct TxLine {
    /// True if the line is lossless (R/L < 5e5 and G < 1e-2).
    pub lsl: bool,
    /// True if timestep > tau (extended time step).
    pub ext: bool,
    /// Ratio for extended time step interpolation.
    pub ratio: f64,
    /// Propagation delay in picoseconds.
    pub taul: f64,
    /// sqrt(C/L) — characteristic admittance scaling.
    pub sqt_cdl: f64,
    /// Attenuation factor for h2 (propagation).
    pub h2_aten: f64,
    /// Attenuation factor for h3 (combined Y*exp).
    pub h3_aten: f64,
    /// Sum of h1 residues (h1C = sum of h1_term[i].c after scaling).
    pub h1c: f64,
    /// Cached exp(h1_term[i].x * h) values.
    pub h1e: [f64; 3],
    /// Whether the propagation function has complex conjugate poles.
    pub if_img: bool,
    /// h1 terms: 3 poles for Y(s) Padé approximation.
    pub h1_term: [Term; 3],
    /// h2 terms: 3 poles for exp(-gamma*l) Padé approximation.
    pub h2_term: [Term; 3],
    /// h3 terms: 6 poles for Y(s)*exp(-gamma*l) combined.
    pub h3_term: [Term; 6],
    /// Voltage/current history list.
    pub vi_history: Vec<ViEntry>,
    /// DC operating point at input port.
    pub dc1: f64,
    /// DC operating point at output port.
    pub dc2: f64,
}

/// A single history entry: time (in picoseconds), voltages and currents.
#[derive(Debug, Clone)]
pub struct ViEntry {
    pub time: i64,
    pub v_i: f64,
    pub v_o: f64,
    pub i_i: f64,
    pub i_o: f64,
}

impl TxLine {
    fn new() -> Self {
        Self {
            lsl: false,
            ext: false,
            ratio: 0.0,
            taul: 0.0,
            sqt_cdl: 0.0,
            h2_aten: 0.0,
            h3_aten: 0.0,
            h1c: 0.0,
            h1e: [0.0; 3],
            if_img: false,
            h1_term: [Term::zero(), Term::zero(), Term::zero()],
            h2_term: [Term::zero(), Term::zero(), Term::zero()],
            h3_term: [
                Term::zero(),
                Term::zero(),
                Term::zero(),
                Term::zero(),
                Term::zero(),
                Term::zero(),
            ],
            vi_history: Vec::new(),
            dc1: 0.0,
            dc2: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// TXL Model
// ---------------------------------------------------------------------------

/// TXL model parameters: per-unit-length R, L, G, C and length.
#[derive(Debug, Clone)]
pub struct TxlModel {
    pub r: f64,
    pub l: f64,
    pub g: f64,
    pub c: f64,
    pub length: f64,
}

impl TxlModel {
    pub fn from_model_def(def: &ModelDef) -> Self {
        let mut r = 0.0;
        let mut l = 0.0;
        let mut g = 0.0;
        let mut c = 0.0;
        let mut length = 0.0;

        for p in &def.params {
            let val = crate::expr_val_or(&p.value, 0.0);
            match p.name.to_uppercase().as_str() {
                "R" => r = val,
                "L" => l = val,
                "G" => g = val,
                "C" => c = val,
                "LENGTH" | "LEN" => length = val,
                _ => {}
            }
        }

        Self { r, l, g, c, length }
    }
}

// ---------------------------------------------------------------------------
// TXL Instance
// ---------------------------------------------------------------------------

/// A resolved TXL instance in the MNA system.
#[derive(Debug, Clone)]
pub struct TxlInstance {
    pub name: String,
    pub pos_idx: Option<usize>,
    pub neg_idx: Option<usize>,
    pub ibr1: usize,
    pub ibr2: usize,
    pub model: TxlModel,
    /// Primary line state.
    pub txline: TxLine,
    /// Backup line state (for rollback during NR).
    pub txline2: TxLine,
    /// Instance-specific length override.
    pub length: f64,
    /// DC initialization done flag.
    pub dc_given: bool,
}

impl TxlInstance {
    /// Get the effective line resistance for DC.
    pub fn dc_resistance(&self) -> f64 {
        self.model.r * self.length
    }
}

// ---------------------------------------------------------------------------
// Setup: Padé approximation
// ---------------------------------------------------------------------------

const EPSI: f64 = 1.0e-16;
const EPSI2: f64 = 1.0e-28;

/// Maclaurin series coefficients for F(z) = sqrt((1+az)/(1+bz)).
fn mac(at: f64, bt: f64) -> (f64, f64, f64, f64, f64) {
    let a = at;
    let b = bt;

    let y1 = 0.5 * (a - b);
    let y2 = 0.5 * (3.0 * b * b - 2.0 * a * b - a * a) * y1 / (a - b);
    let y3 = ((3.0 * b * b + a * a) * y1 * y1 + 0.5 * (3.0 * b * b - 2.0 * a * b - a * a) * y2)
        / (a - b);
    let y4 = ((3.0 * b * b - 3.0 * a * a) * y1 * y1 * y1
        + (9.0 * b * b + 3.0 * a * a) * y1 * y2
        + 0.5 * (3.0 * b * b - 2.0 * a * b - a * a) * y3)
        / (a - b);
    let y5 = (12.0 * a * a * y1 * y1 * y1 * y1
        + y1 * y1 * y2 * (18.0 * b * b - 18.0 * a * a)
        + (9.0 * b * b + 3.0 * a * a) * (y2 * y2 + y1 * y3)
        + (3.0 * b * b + a * a) * y1 * y3
        + 0.5 * (3.0 * b * b - 2.0 * a * b - a * a) * y4)
        / (a - b);

    (y1, y2 / 2.0, y3 / 6.0, y4 / 24.0, y5 / 120.0)
}

fn eval2(a: f64, b: f64, c: f64, x: f64) -> f64 {
    a * x * x + b * x + c
}

/// Newton-Raphson step for cubic root finding.
fn root3(a1: f64, a2: f64, a3: f64, x: f64) -> f64 {
    let t1 = x * x * x + a1 * x * x + a2 * x + a3;
    let t2 = 3.0 * x * x + 2.0 * a1 * x + a2;
    x - t1 / t2
}

/// Polynomial division: (x^3 + a1*x^2 + a2*x + a3) / (x - root) = x^2 + p1*x + p2.
fn div3(a1: f64, _a2: f64, a3: f64, x: f64) -> (f64, f64) {
    let p1 = a1 + x;
    let p2 = -a3 / x;
    (p1, p2)
}

/// Gaussian elimination on a 3×4 augmented matrix.
fn gaussian_elimination(a: &mut [[f64; 4]; 3], epsi_val: f64) -> bool {
    for i in 0..3 {
        let mut imax = i;
        let mut max = a[i][i].abs();
        for (j, row) in a.iter().enumerate().skip(i + 1) {
            if row[i].abs() > max {
                imax = j;
                max = row[i].abs();
            }
        }
        if max < epsi_val {
            return false;
        }
        if imax != i {
            a.swap(i, imax);
        }
        let f = 1.0 / a[i][i];
        a[i][i] = 1.0;
        for val in a[i].iter_mut().skip(i + 1) {
            *val *= f;
        }
        for j in 0..3 {
            if i == j {
                continue;
            }
            let f = a[j][i];
            a[j][i] = 0.0;
            // Split borrow: read from row i, write to row j
            let (row_i, row_j) = if j > i {
                let (lo, hi) = a.split_at_mut(j);
                (&lo[i], &mut hi[0])
            } else {
                let (lo, hi) = a.split_at_mut(i);
                (&hi[0], &mut lo[j])
            };
            for (rj, ri) in row_j[(i + 1)..].iter_mut().zip(row_i[(i + 1)..].iter()) {
                *rj -= f * ri;
            }
        }
    }
    true
}

/// Find roots of cubic x^3 + a1*x^2 + a2*x + a3 = 0.
/// Returns (x1, x2, x3) where all roots are real.
fn find_roots_real(a1: f64, a2: f64, a3: f64) -> (f64, f64, f64) {
    let q = (a1 * a1 - 3.0 * a2) / 9.0;
    let p = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
    let t = q * q * q - p * p;
    let mut x = if t >= 0.0 {
        let t = (p / (q * q.sqrt())).acos();
        -2.0 * q.sqrt() * (t / 3.0).cos() - a1 / 3.0
    } else if p > 0.0 {
        let t = ((-t).sqrt() + p).cbrt();
        -(t + q / t) - a1 / 3.0
    } else if p == 0.0 {
        -a1 / 3.0
    } else {
        let t = ((-t).sqrt() - p).cbrt();
        (t + q / t) - a1 / 3.0
    };

    // Refine with Newton-Raphson
    let x_backup = x;
    for i in 0..32 {
        let t = root3(a1, a2, a3, x);
        if (t - x).abs() <= 5.0e-4 {
            break;
        }
        if i == 31 {
            x = x_backup;
            break;
        }
        x = t;
    }

    let x1 = x;
    let (qa1, qa2) = div3(a1, a2, a3, x);

    let disc = qa1 * qa1 - 4.0 * qa2;
    // Scale to avoid precision loss (matching ngspice)
    let disc_scaled = disc * 1.0e-18;
    let t = disc_scaled.abs().sqrt() * 1.0e9;

    let (x2, x3) = if qa1 >= 0.0 {
        let x2 = -0.5 * (qa1 + t);
        (x2, qa2 / x2)
    } else {
        let x2 = -0.5 * (qa1 - t);
        (x2, qa2 / x2)
    };

    (x1, x2, x3)
}

/// Find roots of cubic for the propagation function.
/// May return complex conjugate pair.
fn find_roots_exp(a1: f64, a2: f64, a3: f64) -> (f64, f64, f64, bool) {
    let q = (a1 * a1 - 3.0 * a2) / 9.0;
    let p = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
    let t = q * q * q - p * p;
    let mut x = if t >= 0.0 {
        let t = (p / (q * q.sqrt())).acos();
        -2.0 * q.sqrt() * (t / 3.0).cos() - a1 / 3.0
    } else if p > 0.0 {
        let t = ((-t).sqrt() + p).cbrt();
        -(t + q / t) - a1 / 3.0
    } else if p == 0.0 {
        -a1 / 3.0
    } else {
        let t = ((-t).sqrt() - p).cbrt();
        (t + q / t) - a1 / 3.0
    };

    // Refine
    let x_backup = x;
    for i in 0..32 {
        let t = root3(a1, a2, a3, x);
        if (t - x).abs() <= 5.0e-4 {
            break;
        }
        if i == 31 {
            x = x_backup;
            break;
        }
        x = t;
    }

    let ex1 = x;
    let (qa1, qa2) = div3(a1, a2, a3, x);

    let disc = qa1 * qa1 - 4.0 * qa2;
    if disc < 0.0 {
        // Complex conjugate pair
        let ex3 = 0.5 * (-disc).sqrt(); // imaginary part
        let ex2 = -0.5 * qa1; // real part
        (ex1, ex2, ex3, true)
    } else {
        let disc_scaled = disc * 1.0e-16;
        let t = disc_scaled.abs().sqrt() * 1.0e8;
        let (ex2, ex3) = if qa1 >= 0.0 {
            let ex2 = -0.5 * (qa1 + t);
            (ex2, qa2 / ex2)
        } else {
            let ex2 = -0.5 * (qa1 - t);
            (ex2, qa2 / ex2)
        };
        (ex1, ex2, ex3, false)
    }
}

/// Complex residue calculation for complex conjugate pole pair.
fn get_c_complex(eq1: f64, eq2: f64, eq3: f64, ep1: f64, ep2: f64, a: f64, b: f64) -> (f64, f64) {
    let d = (3.0 * (a * a - b * b) + 2.0 * ep1 * a + ep2).powi(2)
        + (6.0 * a * b + 2.0 * ep1 * b).powi(2);
    let n_i = -(eq1 * (a * a - b * b) + eq2 * a + eq3) * (6.0 * a * b + 2.0 * ep1 * b)
        + (2.0 * eq1 * a * b + eq2 * b) * (3.0 * (a * a - b * b) + 2.0 * ep1 * a + ep2);
    let ci = n_i / d;
    let n_r = (3.0 * (a * a - b * b) + 2.0 * ep1 * a + ep2)
        * (eq1 * (a * a - b * b) + eq2 * a + eq3)
        + (6.0 * a * b + 2.0 * ep1 * b) * (2.0 * eq1 * a * b + eq2 * b);
    let cr = n_r / d;
    (cr, ci)
}

/// Complex division: (ar + ai*i) / (br + bi*i).
fn div_c(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    let d = br * br + bi * bi;
    let cr = (ar * br + ai * bi) / d;
    let ci = (-ar * bi + ai * br) / d;
    (cr, ci)
}

/// Y(s) Padé approximation: compute h1 terms (3 poles for admittance).
fn y_pade(r: f64, l: f64, g: f64, c: f64, h: &mut TxLine) {
    let sqt_cdl = (c / l).sqrt();
    let rdl = r / l;
    let gdc = g / c;

    let (b1, b2, b3, b4, b5) = mac(gdc, rdl);

    let mut a = [[0.0f64; 4]; 3];
    a[0][0] = 1.0 - (gdc / rdl).sqrt();
    a[0][1] = b1;
    a[0][2] = b2;
    a[0][3] = -b3;
    a[1][0] = b1;
    a[1][1] = b2;
    a[1][2] = b3;
    a[1][3] = -b4;
    a[2][0] = b2;
    a[2][1] = b3;
    a[2][2] = b4;
    a[2][3] = -b5;

    gaussian_elimination(&mut a, EPSI);

    let p3 = a[0][3];
    let p2 = a[1][3];
    let p1 = a[2][3];

    let q1 = p1 + b1;
    let q2 = b1 * p1 + p2 + b2;
    let q3 = p3 * (gdc / rdl).sqrt();

    let (x1, x2, x3) = find_roots_real(p1, p2, p3);
    let c1 = eval2(q1 - p1, q2 - p2, q3 - p3, x1) / eval2(3.0, 2.0 * p1, p2, x1);
    let c2 = eval2(q1 - p1, q2 - p2, q3 - p3, x2) / eval2(3.0, 2.0 * p1, p2, x2);
    let c3 = eval2(q1 - p1, q2 - p2, q3 - p3, x3) / eval2(3.0, 2.0 * p1, p2, x3);

    h.sqt_cdl = sqt_cdl;
    h.h1_term[0].c = c1;
    h.h1_term[1].c = c2;
    h.h1_term[2].c = c3;
    h.h1_term[0].x = x1;
    h.h1_term[1].x = x2;
    h.h1_term[2].x = x3;
}

/// exp(-gamma*l) Padé approximation: compute h2 terms.
fn exp_pade(r: f64, l: f64, g: f64, c: f64, len: f64, h: &mut TxLine) {
    let tau = (l * c).sqrt();
    let rdl = r / l;
    let gdc = g / c;
    let rg = r * g;

    // Maclaurin series of the exponent
    let a_val = rdl;
    let b_val = gdc;

    let y1 = 0.5 * (a_val + b_val);
    let y2 = a_val * b_val - y1 * y1;
    let y3 = -3.0 * y1 * y2;
    let y4 = -3.0 * y2 * y2 - 4.0 * y1 * y3;
    let y5 = -5.0 * y1 * y4 - 10.0 * y2 * y3;
    let y6 = -10.0 * y3 * y3 - 15.0 * y2 * y4 - 6.0 * y1 * y5;

    let t = tau;
    let a0 = y1 * t * len;
    let a1 = y2 * t * t / 2.0 * len;
    let a2 = y3 * t * t * t / 6.0 * len;
    let a3 = y4 * t.powi(4) / 24.0 * len;
    let a4 = y5 * t.powi(5) / 120.0 * len;
    let a5 = y6 * t.powi(6) / 720.0 * len;

    // Padé approximation of the modified exponential
    let mut b = [0.0f64; 6];
    let aa = [-a1, -a2, -a3, -a4, -a5];
    b[0] = 1.0;
    b[1] = aa[0];
    for i in 2..=5 {
        b[i] = 0.0;
        for j in 1..=i {
            b[i] += j as f64 * aa[j - 1] * b[i - j];
        }
        b[i] /= i as f64;
    }

    let rg_sqrt = rg.sqrt();
    let mut aa_mat = [[0.0f64; 4]; 3];
    aa_mat[0][0] = 1.0 - (a0 - len * rg_sqrt).exp();
    aa_mat[0][1] = b[1];
    aa_mat[0][2] = b[2];
    aa_mat[0][3] = -b[3];
    aa_mat[1][0] = b[1];
    aa_mat[1][1] = b[2];
    aa_mat[1][2] = b[3];
    aa_mat[1][3] = -b[4];
    aa_mat[2][0] = b[2];
    aa_mat[2][1] = b[3];
    aa_mat[2][2] = b[4];
    aa_mat[2][3] = -b[5];

    gaussian_elimination(&mut aa_mat, EPSI2);

    let ep3 = aa_mat[0][3];
    let ep2 = aa_mat[1][3];
    let ep1 = aa_mat[2][3];

    let eq1 = ep1 + b[1];
    let eq2 = b[1] * ep1 + ep2 + b[2];
    let eq3 = ep3 * (a0 - len * rg_sqrt).exp();

    // Scale by tau
    let ep3s = ep3 / (tau * tau * tau);
    let ep2s = ep2 / (tau * tau);
    let ep1s = ep1 / tau;
    let eq3s = eq3 / (tau * tau * tau);
    let eq2s = eq2 / (tau * tau);
    let eq1s = eq1 / tau;

    let (ex1, ex2, ex3, if_img) = find_roots_exp(ep1s, ep2s, ep3s);

    let ec1 = eval2(eq1s - ep1s, eq2s - ep2s, eq3s - ep3s, ex1) / eval2(3.0, 2.0 * ep1s, ep2s, ex1);
    let (ec2, ec3) = if if_img {
        get_c_complex(eq1s - ep1s, eq2s - ep2s, eq3s - ep3s, ep1s, ep2s, ex2, ex3)
    } else {
        let ec2 =
            eval2(eq1s - ep1s, eq2s - ep2s, eq3s - ep3s, ex2) / eval2(3.0, 2.0 * ep1s, ep2s, ex2);
        let ec3 =
            eval2(eq1s - ep1s, eq2s - ep2s, eq3s - ep3s, ex3) / eval2(3.0, 2.0 * ep1s, ep2s, ex3);
        (ec2, ec3)
    };

    h.taul = tau * len;
    h.h2_aten = (-a0).exp();
    h.h2_term[0].c = ec1;
    h.h2_term[1].c = ec2;
    h.h2_term[2].c = ec3;
    h.h2_term[0].x = ex1;
    h.h2_term[1].x = ex2;
    h.h2_term[2].x = ex3;
    h.if_img = if_img;
}

/// Compute h3 = h1 * h2 combined terms.
fn get_h3(h: &mut TxLine) {
    h.h3_aten = h.h2_aten * h.sqt_cdl;

    let (cc1, cc2, cc3) = (h.h1_term[0].c, h.h1_term[1].c, h.h1_term[2].c);
    let (xx1, xx2, xx3) = (h.h1_term[0].x, h.h1_term[1].x, h.h1_term[2].x);
    let (cc4, cc5, cc6) = (h.h2_term[0].c, h.h2_term[1].c, h.h2_term[2].c);
    let (xx4, xx5, xx6) = (h.h2_term[0].x, h.h2_term[1].x, h.h2_term[2].x);

    // Copy pole locations
    h.h3_term[0].x = xx1;
    h.h3_term[1].x = xx2;
    h.h3_term[2].x = xx3;
    h.h3_term[3].x = xx4;
    h.h3_term[4].x = xx5;
    h.h3_term[5].x = xx6;

    if h.if_img {
        h.h3_term[0].c = cc1
            + cc1
                * (cc4 / (xx1 - xx4)
                    + 2.0 * (cc5 * xx1 - xx6 * cc6 - xx5 * cc5)
                        / (xx1 * xx1 - 2.0 * xx5 * xx1 + xx5 * xx5 + xx6 * xx6));
        h.h3_term[1].c = cc2
            + cc2
                * (cc4 / (xx2 - xx4)
                    + 2.0 * (cc5 * xx2 - xx6 * cc6 - xx5 * cc5)
                        / (xx2 * xx2 - 2.0 * xx5 * xx2 + xx5 * xx5 + xx6 * xx6));
        h.h3_term[2].c = cc3
            + cc3
                * (cc4 / (xx3 - xx4)
                    + 2.0 * (cc5 * xx3 - xx6 * cc6 - xx5 * cc5)
                        / (xx3 * xx3 - 2.0 * xx5 * xx3 + xx5 * xx5 + xx6 * xx6));

        h.h3_term[3].c = cc4 + cc4 * (cc1 / (xx4 - xx1) + cc2 / (xx4 - xx2) + cc3 / (xx4 - xx3));

        h.h3_term[4].c = cc5;
        h.h3_term[5].c = cc6;

        let (r, im) = div_c(cc5, cc6, xx5 - xx1, xx6);
        h.h3_term[4].c += r * cc1;
        h.h3_term[5].c += im * cc1;
        let (r, im) = div_c(cc5, cc6, xx5 - xx2, xx6);
        h.h3_term[4].c += r * cc2;
        h.h3_term[5].c += im * cc2;
        let (r, im) = div_c(cc5, cc6, xx5 - xx3, xx6);
        h.h3_term[4].c += r * cc3;
        h.h3_term[5].c += im * cc3;
    } else {
        h.h3_term[0].c = cc1 + cc1 * (cc4 / (xx1 - xx4) + cc5 / (xx1 - xx5) + cc6 / (xx1 - xx6));
        h.h3_term[1].c = cc2 + cc2 * (cc4 / (xx2 - xx4) + cc5 / (xx2 - xx5) + cc6 / (xx2 - xx6));
        h.h3_term[2].c = cc3 + cc3 * (cc4 / (xx3 - xx4) + cc5 / (xx3 - xx5) + cc6 / (xx3 - xx6));

        h.h3_term[3].c = cc4 + cc4 * (cc1 / (xx4 - xx1) + cc2 / (xx4 - xx2) + cc3 / (xx4 - xx3));
        h.h3_term[4].c = cc5 + cc5 * (cc1 / (xx5 - xx1) + cc2 / (xx5 - xx2) + cc3 / (xx5 - xx3));
        h.h3_term[5].c = cc6 + cc6 * (cc1 / (xx6 - xx1) + cc2 / (xx6 - xx2) + cc3 / (xx6 - xx3));
    }
}

/// Final scaling of h1, h2, h3 coefficients and compute h1C.
fn update_h1c(h: &mut TxLine) {
    let mut d = 0.0;
    for i in 0..3 {
        h.h1_term[i].c *= h.sqt_cdl;
        d += h.h1_term[i].c;
    }
    h.h1c = d;

    for i in 0..3 {
        h.h2_term[i].c *= h.h2_aten;
    }
    for i in 0..6 {
        h.h3_term[i].c *= h.h3_aten;
    }
}

/// Full Padé approximation setup for a lossy line.
fn main_pade(r: f64, l: f64, g: f64, c: f64, len: f64, h: &mut TxLine) {
    y_pade(r, l, g, c, h);
    exp_pade(r, l, g, c, len, h);
    get_h3(h);
    h.taul *= 1.0e12; // Convert to picoseconds
    update_h1c(h);
}

/// Initialize TXL line from model parameters.
pub fn setup_txline(model: &TxlModel, length: f64) -> TxLine {
    let mut t = TxLine::new();
    let r = model.r;
    let l = model.l.max(1e-12); // minimum inductance
    let c = model.c;
    let g = model.g;

    if r / l < 5.0e5 && g < 1.0e-2 {
        t.lsl = true;
        t.taul = (c * l).sqrt() * length * 1.0e12;
        t.h3_aten = (c / l).sqrt();
        t.sqt_cdl = t.h3_aten;
        t.h2_aten = 1.0;
        t.h1c = 0.0;
    }

    if !t.lsl {
        main_pade(r, l, g, c, length, &mut t);
    }

    t
}

// ---------------------------------------------------------------------------
// Transient load
// ---------------------------------------------------------------------------

/// Complex exponential: e^((ar + ai*i)*h).
fn exp_c(ar: f64, ai: f64, h: f64) -> (f64, f64) {
    let e = (ar * h).exp();
    let cs = (ai * h).cos();
    let si = (ai * h).sin();
    (e * cs, e * si)
}

/// Complex multiplication: (ar+ai*i)*(br+bi*i).
fn mult_c(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    (ar * br - ai * bi, ar * bi + ai * br)
}

/// Update h1 convolution accumulators with new voltage data.
pub fn update_cnv_txl(tx: &mut TxLine, h: f64, v_i: f64, v_o: f64, dv_i: f64, dv_o: f64) {
    for i in 0..3 {
        let e = tx.h1e[i];
        let t = tx.h1_term[i].c / tx.h1_term[i].x;
        let bi = dv_i * t;
        let bo = dv_o * t;
        tx.h1_term[i].cnv_i = (tx.h1_term[i].cnv_i - bi * h) * e
            + (e - 1.0) * (v_i * t + 1.0e12 * bi / tx.h1_term[i].x);
        tx.h1_term[i].cnv_o = (tx.h1_term[i].cnv_o - bo * h) * e
            + (e - 1.0) * (v_o * t + 1.0e12 * bo / tx.h1_term[i].x);
    }
}

/// Get previous voltage/current values at delayed time (time - tau).
/// Returns (v1_i, v2_i, i1_i, i2_i, v1_o, v2_o, i1_o, i2_o, ext_flag, ratio).
fn get_pvs_vi(
    t1: i64,
    t2: i64,
    history: &[ViEntry],
    dc1: f64,
    dc2: f64,
    taul: f64,
) -> (f64, f64, f64, f64, f64, f64, f64, f64, bool, f64) {
    let ta = t1 as f64 - taul;
    let tb = t2 as f64 - taul;

    if tb <= 0.0 {
        return (dc1, dc1, 0.0, 0.0, dc2, dc2, 0.0, 0.0, false, 0.0);
    }

    let (v1_i, v1_o, i1_i, i1_o);
    let mut start_idx = 0;

    if ta <= 0.0 {
        i1_i = 0.0;
        i1_o = 0.0;
        v1_i = dc1;
        v1_o = dc2;
    } else {
        // Find history entries bracketing ta
        let mut idx = 0;
        while idx + 1 < history.len() && (history[idx + 1].time as f64) < ta {
            idx += 1;
        }
        let vi1 = &history[idx];
        let vi2 = if idx + 1 < history.len() {
            &history[idx + 1]
        } else {
            vi1
        };
        let f = if vi2.time != vi1.time {
            (ta - vi1.time as f64) / (vi2.time - vi1.time) as f64
        } else {
            0.0
        };
        v1_i = vi1.v_i + f * (vi2.v_i - vi1.v_i);
        v1_o = vi1.v_o + f * (vi2.v_o - vi1.v_o);
        i1_i = vi1.i_i + f * (vi2.i_i - vi1.i_i);
        i1_o = vi1.i_o + f * (vi2.i_o - vi1.i_o);
        start_idx = idx;
    }

    let (v2_i, v2_o, i2_i, i2_o, ext, ratio);

    if tb > t1 as f64 {
        // Extended time step: timestep > tau
        ext = true;
        ratio = (tb - t1 as f64) / (t2 - t1) as f64;
        // Use the current time point with scaling
        let vi = &history[history.len() - 1];
        let f = 1.0 - ratio;
        v2_i = vi.v_i * f;
        v2_o = vi.v_o * f;
        i2_i = vi.i_i * f;
        i2_o = vi.i_o * f;
    } else {
        ext = false;
        ratio = 0.0;
        // Find history entries bracketing tb
        let mut idx = start_idx;
        while idx + 1 < history.len() && (history[idx + 1].time as f64) < tb {
            idx += 1;
        }
        let vi1 = &history[idx];
        let vi2 = if idx + 1 < history.len() {
            &history[idx + 1]
        } else {
            vi1
        };
        let f = if vi2.time != vi1.time {
            (tb - vi1.time as f64) / (vi2.time - vi1.time) as f64
        } else {
            0.0
        };
        v2_i = vi1.v_i + f * (vi2.v_i - vi1.v_i);
        v2_o = vi1.v_o + f * (vi2.v_o - vi1.v_o);
        i2_i = vi1.i_i + f * (vi2.i_i - vi1.i_i);
        i2_o = vi1.i_o + f * (vi2.i_o - vi1.i_o);
    }

    (v1_i, v2_i, i1_i, i2_i, v1_o, v2_o, i1_o, i2_o, ext, ratio)
}

/// Compute the RHS constants for the TXL branch equations.
/// Returns (ff, gg, ext, ratio) where ff goes to ibr1 RHS, gg to ibr2 RHS.
pub fn right_consts_txl(
    tx: &mut TxLine,
    t: i64,
    time: i64,
    h: f64,
    h1: f64,
    v_in: f64,
    v_out: f64,
) -> (f64, f64, bool, f64) {
    let mut ff = 0.0;
    let mut gg = 0.0;

    if !tx.lsl {
        let mut ff1 = 0.0;
        for i in 0..3 {
            tx.h1e[i] = (tx.h1_term[i].x * h).exp();
            ff1 -= tx.h1_term[i].c * tx.h1e[i];
            ff -= tx.h1_term[i].cnv_i * tx.h1e[i];
            gg -= tx.h1_term[i].cnv_o * tx.h1e[i];
        }
        ff += ff1 * h1 * v_in;
        gg += ff1 * h1 * v_out;
    }

    let (v1_i, v2_i, i1_i, i2_i, v1_o, v2_o, i1_o, i2_o, ext, ratio) =
        get_pvs_vi(t, time, &tx.vi_history, tx.dc1, tx.dc2, tx.taul);

    if tx.lsl {
        ff = tx.h3_aten * v2_o + tx.h2_aten * i2_o;
        gg = tx.h3_aten * v2_i + tx.h2_aten * i2_i;
    } else if tx.if_img {
        // Complex pole pair handling
        for i in 0..4 {
            let tm = &mut tx.h3_term[i];
            let e = (tm.x * h).exp();
            tm.cnv_i = tm.cnv_i * e + h1 * tm.c * (v1_i * e + v2_i);
            tm.cnv_o = tm.cnv_o * e + h1 * tm.c * (v1_o * e + v2_o);
        }
        let (er, ei) = exp_c(tx.h3_term[4].x, tx.h3_term[5].x, h);
        let a2 = h1 * tx.h3_term[4].c;
        let b2 = h1 * tx.h3_term[5].c;

        let (a, b) = mult_c(tx.h3_term[4].cnv_i, tx.h3_term[5].cnv_i, er, ei);
        let (a1, b1) = mult_c(a2, b2, v1_i * er + v2_i, v1_i * ei);
        tx.h3_term[4].cnv_i = a + a1;
        tx.h3_term[5].cnv_i = b + b1;

        let (a, b) = mult_c(tx.h3_term[4].cnv_o, tx.h3_term[5].cnv_o, er, ei);
        let (a1, b1) = mult_c(a2, b2, v1_o * er + v2_o, v1_o * ei);
        tx.h3_term[4].cnv_o = a + a1;
        tx.h3_term[5].cnv_o = b + b1;

        ff += tx.h3_aten * v2_o;
        gg += tx.h3_aten * v2_i;

        for i in 0..5 {
            ff += tx.h3_term[i].cnv_o;
            gg += tx.h3_term[i].cnv_i;
        }
        ff += tx.h3_term[4].cnv_o;
        gg += tx.h3_term[4].cnv_i;

        // h2 terms with complex poles
        {
            let tm = &mut tx.h2_term[0];
            let e = (tm.x * h).exp();
            tm.cnv_i = tm.cnv_i * e + h1 * tm.c * (i1_i * e + i2_i);
            tm.cnv_o = tm.cnv_o * e + h1 * tm.c * (i1_o * e + i2_o);
        }
        let (er, ei) = exp_c(tx.h2_term[1].x, tx.h2_term[2].x, h);
        let a2 = h1 * tx.h2_term[1].c;
        let b2 = h1 * tx.h2_term[2].c;

        let (a, b) = mult_c(tx.h2_term[1].cnv_i, tx.h2_term[2].cnv_i, er, ei);
        let (a1, b1) = mult_c(a2, b2, i1_i * er + i2_i, i1_i * ei);
        tx.h2_term[1].cnv_i = a + a1;
        tx.h2_term[2].cnv_i = b + b1;

        let (a, b) = mult_c(tx.h2_term[1].cnv_o, tx.h2_term[2].cnv_o, er, ei);
        let (a1, b1) = mult_c(a2, b2, i1_o * er + i2_o, i1_o * ei);
        tx.h2_term[1].cnv_o = a + a1;
        tx.h2_term[2].cnv_o = b + b1;

        ff += tx.h2_aten * i2_o + tx.h2_term[0].cnv_o + 2.0 * tx.h2_term[1].cnv_o;
        gg += tx.h2_aten * i2_i + tx.h2_term[0].cnv_i + 2.0 * tx.h2_term[1].cnv_i;
    } else {
        // All real poles
        for i in 0..6 {
            let tm = &mut tx.h3_term[i];
            let e = (tm.x * h).exp();
            tm.cnv_i = tm.cnv_i * e + h1 * tm.c * (v1_i * e + v2_i);
            tm.cnv_o = tm.cnv_o * e + h1 * tm.c * (v1_o * e + v2_o);
        }

        ff += tx.h3_aten * v2_o;
        gg += tx.h3_aten * v2_i;

        for i in 0..6 {
            ff += tx.h3_term[i].cnv_o;
            gg += tx.h3_term[i].cnv_i;
        }

        for i in 0..3 {
            let tm = &mut tx.h2_term[i];
            let e = (tm.x * h).exp();
            tm.cnv_i = tm.cnv_i * e + h1 * tm.c * (i1_i * e + i2_i);
            tm.cnv_o = tm.cnv_o * e + h1 * tm.c * (i1_o * e + i2_o);
        }

        ff += tx.h2_aten * i2_o;
        gg += tx.h2_aten * i2_i;

        for i in 0..3 {
            ff += tx.h2_term[i].cnv_o;
            gg += tx.h2_term[i].cnv_i;
        }
    }

    (ff, gg, ext, ratio)
}

/// Update delayed convolution values after extended time step.
pub fn update_delayed_cnv(tx: &mut TxLine, h: f64, ratio: f64) {
    let h_scaled = 0.5e-12 * h;
    let vi = tx.vi_history.last().unwrap();

    if ratio > 0.0 {
        let f_v_i = h_scaled * ratio * vi.v_i;
        for i in 0..6 {
            tx.h3_term[i].cnv_i += f_v_i * tx.h3_term[i].c;
        }
        let f_v_o = h_scaled * ratio * vi.v_o;
        for i in 0..6 {
            tx.h3_term[i].cnv_o += f_v_o * tx.h3_term[i].c;
        }

        let f_i_i = h_scaled * ratio * vi.i_i;
        for i in 0..3 {
            tx.h2_term[i].cnv_i += f_i_i * tx.h2_term[i].c;
        }
        let f_i_o = h_scaled * ratio * vi.i_o;
        for i in 0..3 {
            tx.h2_term[i].cnv_o += f_i_o * tx.h2_term[i].c;
        }
    }
}

/// Initialize DC state for a TXL instance.
pub fn init_dc_state(tx: &mut TxLine, dc_v1: f64, dc_v2: f64, dc_i1: f64, dc_i2: f64) {
    tx.dc1 = dc_v1;
    tx.dc2 = dc_v2;

    // Initialize h1 convolution accumulators
    for i in 0..3 {
        tx.h1_term[i].cnv_i = -dc_v1 * tx.h1_term[i].c / tx.h1_term[i].x;
        tx.h1_term[i].cnv_o = -dc_v2 * tx.h1_term[i].c / tx.h1_term[i].x;
    }
    for i in 0..3 {
        tx.h2_term[i].cnv_i = 0.0;
        tx.h2_term[i].cnv_o = 0.0;
    }
    for i in 0..6 {
        tx.h3_term[i].cnv_i = -dc_v1 * tx.h3_term[i].c / tx.h3_term[i].x;
        tx.h3_term[i].cnv_o = -dc_v2 * tx.h3_term[i].c / tx.h3_term[i].x;
    }

    // Initial history entry
    tx.vi_history.push(ViEntry {
        time: 0,
        v_i: dc_v1,
        v_o: dc_v2,
        i_i: dc_i1,
        i_o: dc_i2,
    });
}

/// Stamp the TXL DC equations into a linear system.
/// DC model: series resistance R*length with two branch currents.
///
/// Equations:
///   KCL at posNode: += ibr1
///   KCL at negNode: += ibr2
///   ibr1 equation: ibr1 + ibr2 = 0  (current conservation)
///   ibr2 equation: V_pos - V_neg - g*ibr1 = 0  (Ohm's law)
pub fn stamp_txl_dc(inst: &TxlInstance, system: &mut crate::LinearSystem, num_nodes: usize) {
    let g = inst.dc_resistance();
    let ibr1 = num_nodes + inst.ibr1;
    let ibr2 = num_nodes + inst.ibr2;

    // KCL at posNode: current ibr1 enters
    if let Some(p) = inst.pos_idx {
        system.matrix.add(p, ibr1, 1.0);
    }

    // KCL at negNode: current ibr2 enters
    if let Some(n) = inst.neg_idx {
        system.matrix.add(n, ibr2, 1.0);
    }

    // ibr1 equation: ibr1 + ibr2 = 0 (current conservation: I1 = -I2)
    system.matrix.add(ibr1, ibr1, 1.0);
    system.matrix.add(ibr1, ibr2, 1.0);

    // ibr2 equation: V_pos - V_neg - g*ibr1 = 0 (series resistance)
    system.matrix.add(ibr2, ibr1, -g);
    if let Some(p) = inst.pos_idx {
        system.matrix.add(ibr2, p, 1.0);
    }
    if let Some(n) = inst.neg_idx {
        system.matrix.add(ibr2, n, -1.0);
    }
}

/// Pre-computed TXL transient stamp data.
/// Computed once per timestep, applied in the load function (possibly multiple
/// times during NR iterations).
#[derive(Debug, Clone)]
pub struct TxlTransientStamp {
    pub ibr1: usize,
    pub ibr2: usize,
    pub pos_idx: Option<usize>,
    pub neg_idx: Option<usize>,
    /// Admittance: sqt_cdl + h1 * h1c.
    pub admittance: f64,
    /// RHS for ibr1 branch equation.
    pub ff: f64,
    /// RHS for ibr2 branch equation.
    pub gg: f64,
    /// Extra matrix entries for coupling (row, col, value).
    pub coupling: Vec<(usize, usize, f64)>,
}

/// Prepare TXL transient stamp data. This updates the TxLine convolution
/// state and must be called exactly once per timestep (before NR iterations).
///
/// `h_seconds` is the timestep in seconds (poles are per-second).
/// `time_ps` and `prev_time_ps` are in picoseconds (for history lookup).
pub fn prepare_txl_transient(
    inst: &mut TxlInstance,
    num_nodes: usize,
    solution: &[f64],
    time_ps: i64,
    prev_time_ps: i64,
    h_seconds: f64,
) -> TxlTransientStamp {
    let tx = &mut inst.txline;
    let ibr1 = num_nodes + inst.ibr1;
    let ibr2 = num_nodes + inst.ibr2;
    let h = h_seconds;
    let h1 = 0.5 * h;

    let v_in = inst.pos_idx.map_or(0.0, |i| solution[i]);
    let v_out = inst.neg_idx.map_or(0.0, |i| solution[i]);

    let admittance = tx.sqt_cdl + h1 * tx.h1c;

    // Compute RHS from convolution (updates internal state).
    // h/h1 are in seconds (matching pole units); times are in picoseconds.
    let (ff, gg, ext, ratio) = right_consts_txl(tx, prev_time_ps, time_ps, h, h1, v_in, v_out);

    // Coupling terms (from delayed signals)
    let mut coupling = Vec::new();
    if ext || tx.lsl {
        let r = if tx.lsl { 1.0 } else { ratio };
        if r > 0.0 || tx.lsl {
            if tx.lsl {
                let f_h3 = r * tx.h3_aten;
                if let Some(n) = inst.neg_idx {
                    coupling.push((ibr1, n, -f_h3));
                }
                if let Some(p) = inst.pos_idx {
                    coupling.push((ibr2, p, -f_h3));
                }
                let f_h2 = r * tx.h2_aten;
                coupling.push((ibr1, ibr2, -f_h2));
                coupling.push((ibr2, ibr1, -f_h2));
            } else {
                tx.ext = true;
                tx.ratio = ratio;
                let sum_h3: f64 = tx.h3_term.iter().map(|t| t.c).sum();
                let f = ratio * (h1 * sum_h3 + tx.h3_aten);
                if let Some(n) = inst.neg_idx {
                    coupling.push((ibr1, n, -f));
                }
                if let Some(p) = inst.pos_idx {
                    coupling.push((ibr2, p, -f));
                }
                let sum_h2: f64 = tx.h2_term.iter().map(|t| t.c).sum();
                let f = ratio * (h1 * sum_h2 + tx.h2_aten);
                coupling.push((ibr1, ibr2, -f));
                coupling.push((ibr2, ibr1, -f));
            }
        }
    } else {
        tx.ext = false;
    }

    TxlTransientStamp {
        ibr1,
        ibr2,
        pos_idx: inst.pos_idx,
        neg_idx: inst.neg_idx,
        admittance,
        ff,
        gg,
        coupling,
    }
}

/// Apply pre-computed TXL transient stamps into a linear system.
/// Safe to call multiple times (e.g. during NR iterations).
pub fn apply_txl_transient(stamp: &TxlTransientStamp, system: &mut crate::LinearSystem) {
    let ibr1 = stamp.ibr1;
    let ibr2 = stamp.ibr2;

    // Branch equation diagonals
    system.matrix.add(ibr1, ibr1, -1.0);
    system.matrix.add(ibr2, ibr2, -1.0);

    // Admittance stamps
    if let Some(p) = stamp.pos_idx {
        system.matrix.add(p, ibr1, 1.0);
        system.matrix.add(ibr1, p, stamp.admittance);
    }
    if let Some(n) = stamp.neg_idx {
        system.matrix.add(n, ibr2, 1.0);
        system.matrix.add(ibr2, n, stamp.admittance);
    }

    // Coupling terms
    for &(row, col, val) in &stamp.coupling {
        system.matrix.add(row, col, val);
    }

    // RHS
    system.rhs[ibr1] += stamp.ff;
    system.rhs[ibr2] += stamp.gg;
}
