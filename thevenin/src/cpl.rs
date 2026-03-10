//! CPL (Coupled Multiconductor Lossy Transmission Lines) model.
//!
//! Ports the ngspice `cpl/` device: Jacobi eigendecomposition of the
//! per-unit-length R, L, G, C matrices followed by polynomial fitting
//! across frequency and Padé approximation to obtain poles/residues for
//! each modal line.
//!
//! Reference: ngspice `src/spicelib/devices/cpl/`

use crate::txl::Term;
use thevenin_types::ModelDef;

// ---------------------------------------------------------------------------
// Constants (matching ngspice multi_line.h / cplsetup.c)
// ---------------------------------------------------------------------------

/// Matrix type alias for 2D Vec used in delayed value results.
type Mat2D = Vec<Vec<f64>>;

/// Return type for `get_pvs_vi`: (ext, v1_i, v2_i, i1_i, i2_i, v1_o, v2_o, i1_o, i2_o, ratio).
type DelayedValues = (
    bool,
    Mat2D,
    Mat2D,
    Mat2D,
    Mat2D,
    Mat2D,
    Mat2D,
    Mat2D,
    Mat2D,
    Vec<f64>,
);

const MAX_DIM: usize = 16;
const LEFT_DEG: usize = 7;
const MAX_DEG: usize = 8;
const EPSILON: f64 = 1.0e-88;
const EPSI_MULT: f64 = 1e-28;
const EPSI2_DIAG: f64 = 1e-8;

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------

/// Pole/residue term set with 3 terms (matching ngspice TMS).
#[derive(Debug, Clone)]
pub struct Tms {
    pub if_img: bool,
    pub aten: f64,
    pub tm: [Term; 3],
}

impl Tms {
    fn new() -> Self {
        Self {
            if_img: false,
            aten: 0.0,
            tm: [Term::zero(), Term::zero(), Term::zero()],
        }
    }
}

/// VI history entry for CPL (vectors per line).
#[derive(Debug, Clone)]
pub struct CplViEntry {
    pub time: i64,
    pub v_i: Vec<f64>,
    pub v_o: Vec<f64>,
    pub i_i: Vec<f64>,
    pub i_o: Vec<f64>,
}

/// CPL line internal state (mirrors ngspice CPLine struct).
#[derive(Debug, Clone)]
pub struct CpLine {
    pub no_l: usize,
    pub ext: bool,
    pub ratio: Vec<f64>,
    pub taul: Vec<f64>,
    /// h1t[i][j]: admittance (Y) terms.
    pub h1t: Vec<Vec<Option<Tms>>>,
    /// h2t[i][j][k]: current coupling through mode k.
    pub h2t: Vec<Vec<Vec<Option<Tms>>>>,
    /// h3t[i][j][k]: voltage coupling through mode k.
    pub h3t: Vec<Vec<Vec<Option<Tms>>>>,
    /// h1C coefficients.
    pub h1c_coeff: Vec<Vec<f64>>,
    /// h2C coefficients.
    pub h2c_coeff: Vec<Vec<Vec<f64>>>,
    /// h3C coefficients.
    pub h3c_coeff: Vec<Vec<Vec<f64>>>,
    /// Cached exp values for h1 terms.
    pub h1e: Vec<Vec<[f64; 3]>>,
    /// DC voltages at input ports.
    pub dc1: Vec<f64>,
    /// DC voltages at output ports.
    pub dc2: Vec<f64>,
    /// Voltage/current history.
    pub vi_history: Vec<CplViEntry>,
}

impl CpLine {
    fn new(no_l: usize) -> Self {
        let mk2 = |n| vec![vec![0.0; n]; n];
        let mk3 = |n| vec![vec![vec![0.0; n]; n]; n];
        Self {
            no_l,
            ext: false,
            ratio: vec![0.0; no_l],
            taul: vec![0.0; no_l],
            h1t: vec![vec![None; no_l]; no_l],
            h2t: vec![vec![vec![None; no_l]; no_l]; no_l],
            h3t: vec![vec![vec![None; no_l]; no_l]; no_l],
            h1c_coeff: mk2(no_l),
            h2c_coeff: mk3(no_l),
            h3c_coeff: mk3(no_l),
            h1e: vec![vec![[0.0; 3]; no_l]; no_l],
            dc1: vec![0.0; no_l],
            dc2: vec![0.0; no_l],
            vi_history: Vec::new(),
        }
    }
}

/// CPL model: per-unit-length R, L, G, C matrices and length.
#[derive(Debug, Clone)]
pub struct CplModel {
    pub r_m: Vec<Vec<f64>>,
    pub l_m: Vec<Vec<f64>>,
    pub g_m: Vec<Vec<f64>>,
    pub c_m: Vec<Vec<f64>>,
    pub length: f64,
    pub dim: usize,
}

impl CplModel {
    /// Parse CPL model from a ModelDef.
    /// R, L, G, C are upper-triangular packed matrices.
    pub fn from_model_def(def: &ModelDef, dim: usize) -> Self {
        let mut r_vals = Vec::new();
        let mut l_vals = Vec::new();
        let mut g_vals = Vec::new();
        let mut c_vals = Vec::new();
        let mut length = 0.0;

        for p in &def.params {
            let val = crate::expr_val_or(&p.value, 0.0);
            match p.name.to_uppercase().as_str() {
                "R" => r_vals.push(val),
                "L" => l_vals.push(val),
                "G" => g_vals.push(val),
                "C" => c_vals.push(val),
                "LENGTH" | "LEN" => length = val,
                _ => {}
            }
        }

        #[allow(clippy::needless_range_loop)]
        let unpack = |vals: &[f64], dim: usize| -> Vec<Vec<f64>> {
            let mut m = vec![vec![0.0; dim]; dim];
            let mut idx = 0;
            for i in 0..dim {
                for j in i..dim {
                    let v = if idx < vals.len() { vals[idx] } else { 0.0 };
                    m[i][j] = v;
                    m[j][i] = v;
                    idx += 1;
                }
            }
            // Ensure R diagonal >= 1e-4 (ngspice clamp)
            m
        };

        let mut r_m = unpack(&r_vals, dim);
        // Clamp R diagonal (ngspice: MAX(f, 1.0e-4))
        for (i, row) in r_m.iter_mut().enumerate().take(dim) {
            if row[i] < 1.0e-4 {
                row[i] = 1.0e-4;
            }
        }

        Self {
            r_m,
            l_m: unpack(&l_vals, dim),
            g_m: unpack(&g_vals, dim),
            c_m: unpack(&c_vals, dim),
            length,
            dim,
        }
    }
}

/// A resolved CPL instance in the MNA system.
#[derive(Debug, Clone)]
pub struct CplInstance {
    pub name: String,
    pub no_l: usize,
    /// Node indices for input ports.
    pub pos_nodes: Vec<Option<usize>>,
    /// Node indices for output ports.
    pub neg_nodes: Vec<Option<usize>>,
    /// Branch equation indices (input side), one per line.
    pub ibr1: Vec<usize>,
    /// Branch equation indices (output side), one per line.
    pub ibr2: Vec<usize>,
    pub model: CplModel,
    pub cpline: CpLine,
    pub cpline2: CpLine,
    pub dc_given: bool,
    pub length: f64,
}

// ---------------------------------------------------------------------------
// Internal setup structures (matching ngspice static globals)
// ---------------------------------------------------------------------------

/// Workspace for coupled-line setup computations.
struct CplSetup {
    dim: usize,
    zy: [[f64; MAX_DIM]; MAX_DIM],
    sv: [[f64; MAX_DIM]; MAX_DIM],
    d: [f64; MAX_DIM],
    y5: [[f64; MAX_DIM]; MAX_DIM],
    y5_1: [[f64; MAX_DIM]; MAX_DIM],
    sv_1: [[f64; MAX_DIM]; MAX_DIM],
    si: [[f64; MAX_DIM]; MAX_DIM],
    si_1: [[f64; MAX_DIM]; MAX_DIM],
    a: [[f64; 2 * MAX_DIM]; MAX_DIM],
    frequency: [f64; MAX_DEG],
    tau: [f64; MAX_DIM],
    scaling_f: f64,
    scaling_f2: f64,
    r_m: [[f64; MAX_DIM]; MAX_DIM],
    g_m: [[f64; MAX_DIM]; MAX_DIM],
    l_m: [[f64; MAX_DIM]; MAX_DIM],
    c_m: [[f64; MAX_DIM]; MAX_DIM],
    length: f64,
    // Polynomial data across frequency samples
    si_sv_1: Vec<Vec<Vec<f64>>>,
    sip: Vec<Vec<Vec<f64>>>,
    si_1p: Vec<Vec<Vec<f64>>>,
    sv_1p: Vec<Vec<Vec<f64>>>,
    w: Vec<Vec<f64>>,
}

/// Output of Padé for a single matrix element.
struct SingleOut {
    c_0: f64,
    poly: Option<Vec<f64>>,
}

/// Output of matrix-polynomial multiplication for one element.
struct MultOut {
    c_0: Vec<f64>,
    poly: Vec<Option<Vec<f64>>>,
}

impl CplSetup {
    fn new(dim: usize, deg_o: usize) -> Self {
        let mk3 = |d, n| vec![vec![vec![0.0; n + 1]; d]; d];
        Self {
            dim,
            zy: [[0.0; MAX_DIM]; MAX_DIM],
            sv: [[0.0; MAX_DIM]; MAX_DIM],
            d: [0.0; MAX_DIM],
            y5: [[0.0; MAX_DIM]; MAX_DIM],
            y5_1: [[0.0; MAX_DIM]; MAX_DIM],
            sv_1: [[0.0; MAX_DIM]; MAX_DIM],
            si: [[0.0; MAX_DIM]; MAX_DIM],
            si_1: [[0.0; MAX_DIM]; MAX_DIM],
            a: [[0.0; 2 * MAX_DIM]; MAX_DIM],
            frequency: [0.0; MAX_DEG],
            tau: [0.0; MAX_DIM],
            scaling_f: 1.0,
            scaling_f2: 1.0,
            r_m: [[0.0; MAX_DIM]; MAX_DIM],
            g_m: [[0.0; MAX_DIM]; MAX_DIM],
            l_m: [[0.0; MAX_DIM]; MAX_DIM],
            c_m: [[0.0; MAX_DIM]; MAX_DIM],
            length: 0.0,
            si_sv_1: mk3(dim, deg_o),
            sip: mk3(dim, deg_o),
            si_1p: mk3(dim, deg_o),
            sv_1p: mk3(dim, deg_o),
            w: vec![vec![0.0; MAX_DEG]; dim],
        }
    }
}

// ---------------------------------------------------------------------------
// Setup: Jacobi diagonalization
// ---------------------------------------------------------------------------

/// Jacobi rotation entry.
struct MaxEntry {
    row: usize,
    col: usize,
    value: f64,
}

#[allow(clippy::needless_range_loop)]
fn diag(
    zy: &mut [[f64; MAX_DIM]; MAX_DIM],
    sv: &mut [[f64; MAX_DIM]; MAX_DIM],
    d: &mut [f64; MAX_DIM],
    dim: usize,
) {
    // Scale matrix
    let mut fmin = zy[0][0].abs();
    let mut fmax = fmin;
    for i in 0..dim {
        for j in i..dim {
            let v = zy[i][j].abs();
            if v > fmax {
                fmax = v;
            }
            if v < fmin {
                fmin = v;
            }
        }
    }
    let scale = 2.0 / (fmin + fmax);
    for i in 0..dim {
        for j in i..dim {
            zy[i][j] *= scale;
        }
    }

    // Initialize eigenvectors to identity
    for i in 0..dim {
        for j in 0..dim {
            sv[i][j] = if i == j { 1.0 } else { 0.0 };
        }
    }

    // Build sorted list of off-diagonal elements
    let mut entries: Vec<MaxEntry> = Vec::new();
    for i in 0..dim.saturating_sub(1) {
        let mut m = i + 1;
        let mut mv = zy[i][m].abs();
        for j in (m + 1)..dim {
            if (zy[i][j].abs() * 1e7) as i64 > (mv * 1e7) as i64 {
                mv = zy[i][j].abs();
                m = j;
            }
        }
        entries.push(MaxEntry {
            row: i,
            col: m,
            value: mv,
        });
    }
    entries.sort_by(|a, b| {
        b.value
            .partial_cmp(&a.value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Iterate Jacobi rotations
    while !entries.is_empty() && entries[0].value > EPSI2_DIAG {
        let p = entries[0].row;
        let q = entries[0].col;

        // Rotate
        rotate(zy, sv, dim, p, q);

        // Reorder entries for rows p and q
        reorder_entries(&mut entries, zy, dim, p);
        if q < dim - 1 {
            reorder_entries(&mut entries, zy, dim, q);
        }
        entries.sort_by(|a, b| {
            b.value
                .partial_cmp(&a.value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    // Extract eigenvalues
    for i in 0..dim {
        d[i] = zy[i][i] / scale;
    }
}

fn reorder_entries(
    entries: &mut [MaxEntry],
    zy: &[[f64; MAX_DIM]; MAX_DIM],
    dim: usize,
    row: usize,
) {
    if row >= dim - 1 {
        return;
    }
    let mut m = row + 1;
    let mut mv = zy[row][m].abs();
    for (j, &val) in zy[row].iter().enumerate().skip(m + 1).take(dim - m - 1) {
        if (val.abs() * 1e7) as i64 > (mv * 1e7) as i64 {
            mv = val.abs();
            m = j;
        }
    }
    if let Some(e) = entries.iter_mut().find(|e| e.row == row) {
        e.col = m;
        e.value = mv;
    }
}

fn rotate(
    zy: &mut [[f64; MAX_DIM]; MAX_DIM],
    sv: &mut [[f64; MAX_DIM]; MAX_DIM],
    dim: usize,
    p: usize,
    q: usize,
) {
    let ld = -zy[p][q];
    let mu = 0.5 * (zy[p][p] - zy[q][q]);
    let ve = (ld * ld + mu * mu).sqrt();
    let co = ((ve + mu.abs()) / (2.0 * ve)).sqrt();
    let si = mu.signum() * ld / (2.0 * ve * co);

    // Save row p
    let mut t = [0.0f64; MAX_DIM];
    t[(p + 1)..dim].copy_from_slice(&zy[p][(p + 1)..dim]);
    for j in 0..p {
        t[j] = zy[j][p];
    }

    for j in (p + 1)..dim {
        if j == q {
            continue;
        }
        if j > q {
            zy[p][j] = t[j] * co - zy[q][j] * si;
        } else {
            zy[p][j] = t[j] * co - zy[j][q] * si;
        }
    }
    for j in (q + 1)..dim {
        if j == p {
            continue;
        }
        zy[q][j] = t[j] * si + zy[q][j] * co;
    }
    for j in 0..p {
        if j == q {
            continue;
        }
        zy[j][p] = t[j] * co - zy[j][q] * si;
    }
    for j in 0..q {
        if j == p {
            continue;
        }
        zy[j][q] = t[j] * si + zy[j][q] * co;
    }

    let old_pp = zy[p][p];
    let old_qq = zy[q][q];
    let old_pq = zy[p][q];
    zy[p][p] = old_pp * co * co + old_qq * si * si - 2.0 * old_pq * si * co;
    zy[q][q] = old_pp * si * si + old_qq * co * co + 2.0 * old_pq * si * co;
    zy[p][q] = 0.0;

    // Update eigenvectors
    let mut tp = [0.0f64; MAX_DIM];
    let mut rr = [0.0f64; MAX_DIM];
    for j in 0..dim {
        tp[j] = sv[j][p];
        rr[j] = sv[j][q];
    }
    for j in 0..dim {
        sv[j][p] = tp[j] * co - rr[j] * si;
        sv[j][q] = tp[j] * si + rr[j] * co;
    }
}

// ---------------------------------------------------------------------------
// Setup: Gaussian elimination
// ---------------------------------------------------------------------------

#[allow(clippy::needless_range_loop)]
fn gaussian_elimination_large(a: &mut [[f64; 2 * MAX_DIM]; MAX_DIM], dims: usize, inv: bool) {
    let dim_cols = if inv { 2 * dims } else { dims };

    for i in 0..dims {
        let mut imax = i;
        let mut max_val = a[i][i].abs();
        for j in (i + 1)..dim_cols {
            if j >= dims && !inv {
                break;
            }
            if j < dims && a[j][i].abs() > max_val {
                imax = j;
                max_val = a[j][i].abs();
            }
        }
        // Also check rows for pivot
        for j in (i + 1)..dims {
            if a[j][i].abs() > max_val {
                imax = j;
                max_val = a[j][i].abs();
            }
        }
        if max_val < EPSILON {
            // Singular — skip
            continue;
        }
        if imax != i {
            for k in i..=dim_cols {
                let tmp = a[i][k];
                a[i][k] = a[imax][k];
                a[imax][k] = tmp;
            }
        }

        let f = 1.0 / a[i][i];
        a[i][i] = 1.0;
        for j in (i + 1)..=dim_cols {
            a[i][j] *= f;
        }
        for j in 0..dims {
            if i == j {
                continue;
            }
            let f = a[j][i];
            a[j][i] = 0.0;
            for k in (i + 1)..=dim_cols {
                a[j][k] -= f * a[i][k];
            }
        }
    }
}

#[allow(clippy::needless_range_loop)]
fn gaussian_elimination_small(at: &mut [[f64; 4]; 4], dims: usize) {
    for i in 0..dims {
        let mut imax = i;
        let mut max_val = at[i][i].abs();
        for j in (i + 1)..dims {
            if at[j][i].abs() > max_val {
                imax = j;
                max_val = at[j][i].abs();
            }
        }
        if max_val < EPSI_MULT {
            continue;
        }
        if imax != i {
            for k in i..=dims {
                let tmp = at[i][k];
                at[i][k] = at[imax][k];
                at[imax][k] = tmp;
            }
        }

        let f = 1.0 / at[i][i];
        at[i][i] = 1.0;
        for j in (i + 1)..=dims {
            at[i][j] *= f;
        }
        for j in 0..dims {
            if i == j {
                continue;
            }
            let f = at[j][i];
            at[j][i] = 0.0;
            for k in (i + 1)..=dims {
                at[j][k] -= f * at[i][k];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Setup: Polynomial fitting
// ---------------------------------------------------------------------------

/// Neville interpolation: given xa[0..n-1], ya[0..n-1] and x, return P(x).
fn polint(xa: &[f64], ya: &[f64], n: usize, x: f64) -> f64 {
    let mut c: Vec<f64> = ya[..n].to_vec();
    let mut d: Vec<f64> = ya[..n].to_vec();

    let mut ns = 0;
    let mut dif = (x - xa[0]).abs();
    for (i, &xa_val) in xa.iter().enumerate().skip(1).take(n - 1) {
        let dift = (x - xa_val).abs();
        if dift < dif {
            ns = i;
            dif = dift;
        }
    }

    let mut y = ya[ns];
    ns += 1; // 1-based for Neville descent
    for m in 1..n {
        for i in 0..(n - m) {
            let ho = xa[i] - x;
            let hp = xa[i + m] - x;
            let w = c[i + 1] - d[i];
            let den = ho - hp;
            if den.abs() < 1e-300 {
                continue;
            }
            let den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        let dy = if 2 * ns < n - m {
            c[ns]
        } else {
            ns -= 1;
            d[ns]
        };
        y += dy;
    }
    y
}

/// Fit polynomial coefficients cof[0..n] such that ya[i] = sum_j cof[j]*xa[i]^j.
fn poly_match(n: usize, xa: &[f64], ya: &[f64]) -> Vec<f64> {
    let mut cof = vec![0.0; n + 1];
    let mut x: Vec<f64> = xa[..=n].to_vec();
    let mut y: Vec<f64> = ya[..=n].to_vec();

    for (j, cof_j) in cof.iter_mut().enumerate().take(n + 1) {
        *cof_j = polint(&x, &y, n + 1 - j, 0.0);
        let mut k_min = 0usize;
        let mut xmin = 1.0e38;
        for i in 0..=(n - j) {
            if x[i].abs() < xmin {
                xmin = x[i].abs();
                k_min = i;
            }
            if x[i] != 0.0 {
                y[i] = (y[i] - *cof_j) / x[i];
            }
        }
        for i in (k_min + 1)..=(n - j) {
            y[i - 1] = y[i];
            x[i - 1] = x[i];
        }
    }
    cof
}

// ---------------------------------------------------------------------------
// Setup: Padé approximation (matching ngspice CPL)
// ---------------------------------------------------------------------------

fn eval2(a: f64, b: f64, c: f64, x: f64) -> f64 {
    a * x * x + b * x + c
}

fn root3(a1: f64, a2: f64, a3: f64, x: f64) -> f64 {
    let t1 = x * (x * (x + a1) + a2) + a3;
    let t2 = x * (2.0 * a1 + 3.0 * x) + a2;
    x - t1 / t2
}

/// Find roots of x^3 + a1*x^2 + a2*x + a3 = 0.
/// Returns (x1, x2, x3, has_complex) where:
///   if has_complex: x1 is real, x2 ± i*x3 are complex conjugate pair
///   else: x1, x2, x3 are all real
fn find_roots(a1: f64, a2: f64, a3: f64) -> (f64, f64, f64, bool) {
    let q = (a1 * a1 - 3.0 * a2) / 9.0;
    let p = (2.0 * a1 * a1 * a1 - 9.0 * a1 * a2 + 27.0 * a3) / 54.0;
    let t = q * q * q - p * p;

    let mut x = if t >= 0.0 {
        let arg = p / (q * q.sqrt());
        let arg = arg.clamp(-1.0, 1.0);
        let t_val = arg.acos();
        -2.0 * q.sqrt() * (t_val / 3.0).cos() - a1 / 3.0
    } else if p > 0.0 {
        let t_val = ((-t).sqrt() + p).cbrt();
        -(t_val + q / t_val) - a1 / 3.0
    } else if p == 0.0 {
        -a1 / 3.0
    } else {
        let t_val = ((-t).sqrt() - p).cbrt();
        (t_val + q / t_val) - a1 / 3.0
    };

    // Newton polish
    let x_backup = x;
    for iter in 0..32 {
        let t = root3(a1, a2, a3, x);
        if (t - x).abs() <= 5.0e-4 {
            break;
        }
        if iter == 31 {
            x = x_backup;
            break;
        }
        x = t;
    }

    let x1 = x;
    // Deflate: divide out (x - x1) to get quadratic
    let new_a1 = a1 + x1;
    let new_a2 = -a3 / x1;

    let disc = new_a1 * new_a1 - 4.0 * new_a2;
    if disc < 0.0 {
        let x3 = 0.5 * (-disc).sqrt();
        let x2 = -0.5 * new_a1;
        (x1, x2, x3, true)
    } else {
        let sq = disc.sqrt();
        let x2 = if new_a1 >= 0.0 {
            -0.5 * (new_a1 + sq)
        } else {
            -0.5 * (new_a1 - sq)
        };
        let x3 = new_a2 / x2;
        (x1, x2, x3, false)
    }
}

fn get_c(q1: f64, q2: f64, q3: f64, p1: f64, p2: f64, a: f64, b: f64) -> (f64, f64) {
    let d1 = 3.0 * (a * a - b * b) + 2.0 * p1 * a + p2;
    let d2 = 6.0 * a * b + 2.0 * p1 * b;
    let d = d1 * d1 + d2 * d2;

    let n_i = -(q1 * (a * a - b * b) + q2 * a + q3) * d2 + (2.0 * q1 * a * b + q2 * b) * d1;
    let ci = n_i / d;

    let n_r = d1 * (q1 * (a * a - b * b) + q2 * a + q3) + d2 * (2.0 * q1 * a * b + q2 * b);
    let cr = n_r / d;

    (cr, ci)
}

/// Padé approximation for CPL.
/// Given b[0..5] (with b[0]=1) and a_b parameter,
/// returns (c1, c2, c3, x1, x2, x3, rtv) where rtv=2 means complex, 1 means real, 0 means failed.
fn pade_apx(a_b: f64, b: &[f64]) -> Option<(f64, f64, f64, f64, f64, f64, i32)> {
    let mut at = [[0.0f64; 4]; 4];
    at[0][0] = 1.0 - a_b;
    at[0][1] = b[1];
    at[0][2] = b[2];
    at[0][3] = -b[3];
    at[1][0] = b[1];
    at[1][1] = b[2];
    at[1][2] = b[3];
    at[1][3] = -b[4];
    at[2][0] = b[2];
    at[2][1] = b[3];
    at[2][2] = b[4];
    at[2][3] = -b[5];

    gaussian_elimination_small(&mut at, 3);

    let p3 = at[0][3];
    let p2 = at[1][3];
    let p1 = at[2][3];

    let q1 = p1 + b[1];
    let q2 = b[1] * p1 + p2 + b[2];
    let q3 = p3 * a_b;

    let (x1, x2, x3, has_complex) = find_roots(p1, p2, p3);

    let c1 = eval2(q1 - p1, q2 - p2, q3 - p3, x1) / eval2(3.0, 2.0 * p1, p2, x1);

    if has_complex {
        let (c2, c3) = get_c(q1 - p1, q2 - p2, q3 - p3, p1, p2, x2, x3);
        Some((c1, c2, c3, x1, x2, x3, 2))
    } else {
        let c2 = eval2(q1 - p1, q2 - p2, q3 - p3, x2) / eval2(3.0, 2.0 * p1, p2, x2);
        let c3 = eval2(q1 - p1, q2 - p2, q3 - p3, x3) / eval2(3.0, 2.0 * p1, p2, x3);
        Some((c1, c2, c3, x1, x2, x3, 1))
    }
}

// ---------------------------------------------------------------------------
// Setup: Mode approximation
// ---------------------------------------------------------------------------

fn approx_mode(x: &[f64], length: f64, scaling_f: f64) -> (f64, Vec<f64>) {
    let w0 = x[0];
    if w0.abs() < 1e-300 {
        return (0.0, vec![0.0; 6]);
    }
    let w1 = x[1] / w0;
    let w2 = x[2] / w0;
    let w3 = x[3] / w0;
    let w4 = x[4] / w0;
    let w5 = x[5] / w0;

    let y1 = 0.5 * w1;
    let y2 = w2 - y1 * y1;
    let y3 = 3.0 * w3 - 3.0 * y1 * y2;
    let y4 = 12.0 * w4 - 3.0 * y2 * y2 - 4.0 * y1 * y3;
    let y5 = 60.0 * w5 - 5.0 * y1 * y4 - 10.0 * y2 * y3;
    let y6 = -10.0 * y3 * y3 - 15.0 * y2 * y4 - 6.0 * y1 * y5;

    let delay = w0.sqrt() * length / scaling_f;
    let atten = (-delay * y1).exp();

    let a = [
        0.0, // unused
        y2 / 2.0 * (-delay),
        y3 / 6.0 * (-delay),
        y4 / 24.0 * (-delay),
        y5 / 120.0 * (-delay),
        y6 / 720.0 * (-delay),
    ];

    let mut b = vec![0.0; 6];
    b[0] = 1.0;
    b[1] = a[1];
    for i in 2..=5 {
        b[i] = 0.0;
        for j in 1..=i {
            b[i] += j as f64 * a[j] * b[i - j];
        }
        b[i] /= i as f64;
    }

    for v in b.iter_mut() {
        *v *= atten;
    }

    (delay, b)
}

// ---------------------------------------------------------------------------
// Setup: Polynomial multiplication
// ---------------------------------------------------------------------------

fn mult_p(p1: &[f64], p2: &[f64], d1: usize, d2: usize, d3: usize) -> Vec<f64> {
    let mut p3 = vec![0.0; d3 + 1];
    for i in 0..=d1 {
        for (k, j) in (i..).zip(0..=d2) {
            if j > d3 {
                break;
            }
            if k > d3 {
                break;
            }
            p3[k] += p1[i] * p2[k - i];
        }
    }
    p3
}

// ---------------------------------------------------------------------------
// Setup: Main coupled() function
// ---------------------------------------------------------------------------

/// Run the full CPL setup: eigendecomposition + polynomial fitting + Padé.
/// Returns the CpLine structure with pre-computed h1t/h2t/h3t terms.
#[allow(clippy::needless_range_loop)]
pub fn setup_cpline(model: &CplModel) -> CpLine {
    let dim = model.dim;
    let deg_o = LEFT_DEG;
    let mut ws = CplSetup::new(dim, deg_o);
    ws.dim = dim;
    ws.length = model.length;

    // Copy model matrices into workspace
    for i in 0..dim {
        for j in 0..dim {
            ws.r_m[i][j] = model.r_m[i][j];
            ws.l_m[i][j] = model.l_m[i][j];
            ws.g_m[i][j] = model.g_m[i][j];
            ws.c_m[i][j] = model.c_m[i][j];
        }
    }

    // Step 0: y=0, ZY = LC
    ws.scaling_f = 1.0;
    ws.scaling_f2 = 1.0;
    loop_zy(&mut ws, 0.0);
    eval_frequency(&mut ws, deg_o);
    eval_si_si_1(&mut ws, 0.0);
    store_si_sv_1(&mut ws, 0);
    store_vals(&mut ws, 0);

    // Steps 1..deg_o
    for idx in 1..=deg_o {
        let freq = ws.frequency[idx];
        loop_zy(&mut ws, freq);
        eval_si_si_1(&mut ws, freq);
        store_si_sv_1(&mut ws, idx);
        store_vals(&mut ws, idx);
    }

    // Fit polynomials to frequency samples
    let freq: Vec<f64> = ws.frequency[..=deg_o].to_vec();
    poly_matrix_fit(&mut ws.sip, dim, deg_o, &freq);
    poly_matrix_fit(&mut ws.si_1p, dim, deg_o, &freq);
    poly_matrix_fit(&mut ws.sv_1p, dim, deg_o, &freq);

    // Fit W polynomials and compute TAU
    for i in 0..dim {
        let w_vals: Vec<f64> = (0..=deg_o).map(|k| ws.w[i][k]).collect();
        let cof = poly_match(deg_o, &freq, &w_vals);
        ws.w[i][..=deg_o].copy_from_slice(&cof[..=deg_o]);
        let (delay, b) = approx_mode(&ws.w[i], ws.length, ws.scaling_f);
        ws.tau[i] = delay;
        for k in 0..6.min(ws.w[i].len()) {
            ws.w[i][k] = if k < b.len() { b[k] } else { 0.0 };
        }
    }

    // Matrix-polynomial multiplications: IWI = Sip * diag(W) * Si_1p
    let iwi = matrix_p_mult_fn(&ws.sip, &ws.w, &ws.si_1p, dim, deg_o);
    // IWV = Sip * diag(W) * Sv_1p
    let iwv = matrix_p_mult_fn(&ws.sip, &ws.w, &ws.sv_1p, dim, deg_o);

    // Fit SiSv_1 polynomials
    poly_matrix_fit(&mut ws.si_sv_1, dim, deg_o, &freq);

    // Generate output: Padé for SIV, IWI, IWV
    let siv = generate_siv(&ws, dim, deg_o);
    let iwi_out = generate_iwi_iwv(&ws, &iwi, dim, true);
    let iwv_out = generate_iwi_iwv(&ws, &iwv, dim, false);

    // Build CpLine
    let mut cp = CpLine::new(dim);
    for i in 0..dim {
        cp.taul[i] = ws.tau[i] * 1.0e12;
    }

    // Fill h1t from SIV
    for i in 0..dim {
        for j in 0..dim {
            let s = &siv[i][j];
            if s.c_0 == 0.0 || s.poly.is_none() {
                cp.h1t[i][j] = None;
                continue;
            }
            let p = s.poly.as_ref().unwrap();
            let d = s.c_0;
            let mut tms = Tms::new();
            tms.aten = d;
            tms.if_img = p[6] as i32 == 2;
            tms.tm[0].c = p[0] * d;
            tms.tm[1].c = p[1] * d;
            tms.tm[2].c = p[2] * d;
            tms.tm[0].x = p[3];
            tms.tm[1].x = p[4];
            tms.tm[2].x = p[5];
            if tms.if_img {
                cp.h1c_coeff[i][j] = tms.tm[0].c + 2.0 * tms.tm[1].c;
            } else {
                cp.h1c_coeff[i][j] = tms.tm[0].c + tms.tm[1].c + tms.tm[2].c;
            }
            cp.h1t[i][j] = Some(tms);
        }
    }

    // Fill h2t from IWI
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                let m = &iwi_out[i][j];
                if m.c_0[k] == 0.0 || m.poly[k].is_none() {
                    cp.h2t[i][j][k] = None;
                    continue;
                }
                let p = m.poly[k].as_ref().unwrap();
                let d = m.c_0[k];
                let mut tms = Tms::new();
                tms.aten = d;
                tms.if_img = p[6] as i32 == 2;
                tms.tm[0].c = p[0] * d;
                tms.tm[1].c = p[1] * d;
                tms.tm[2].c = p[2] * d;
                tms.tm[0].x = p[3];
                tms.tm[1].x = p[4];
                tms.tm[2].x = p[5];
                if tms.if_img {
                    cp.h2c_coeff[i][j][k] = tms.tm[0].c + 2.0 * tms.tm[1].c;
                } else {
                    cp.h2c_coeff[i][j][k] = tms.tm[0].c + tms.tm[1].c + tms.tm[2].c;
                }
                cp.h2t[i][j][k] = Some(tms);
            }
        }
    }

    // Fill h3t from IWV
    for i in 0..dim {
        for j in 0..dim {
            for k in 0..dim {
                let m = &iwv_out[i][j];
                if m.c_0[k] == 0.0 || m.poly[k].is_none() {
                    cp.h3t[i][j][k] = None;
                    continue;
                }
                let p = m.poly[k].as_ref().unwrap();
                let d = m.c_0[k];
                let mut tms = Tms::new();
                tms.aten = d;
                tms.if_img = p[6] as i32 == 2;
                tms.tm[0].c = p[0] * d;
                tms.tm[1].c = p[1] * d;
                tms.tm[2].c = p[2] * d;
                tms.tm[0].x = p[3];
                tms.tm[1].x = p[4];
                tms.tm[2].x = p[5];
                if tms.if_img {
                    cp.h3c_coeff[i][j][k] = tms.tm[0].c + 2.0 * tms.tm[1].c;
                } else {
                    cp.h3c_coeff[i][j][k] = tms.tm[0].c + tms.tm[1].c + tms.tm[2].c;
                }
                cp.h3t[i][j][k] = Some(tms);
            }
        }
    }

    cp
}

// ---------------------------------------------------------------------------
// Setup helpers
// ---------------------------------------------------------------------------

#[allow(clippy::needless_range_loop)]
fn loop_zy(ws: &mut CplSetup, y: f64) {
    let dim = ws.dim;

    // ZY = Scaling_F * C + G * y
    for i in 0..dim {
        for j in 0..dim {
            ws.zy[i][j] = ws.scaling_f * ws.c_m[i][j] + ws.g_m[i][j] * y;
        }
    }

    // First diagonalization (of C+G*y)
    diag(&mut ws.zy, &mut ws.sv, &mut ws.d, dim);

    let mut fmin = ws.d[0];
    for i in 1..dim {
        if ws.d[i] < fmin {
            fmin = ws.d[i];
        }
    }
    if fmin < 0.0 {
        // Not positive definite — clamp
        fmin = 1e-30;
    }
    let fmin_sqrt = fmin.sqrt();
    let fmin1 = 1.0 / fmin_sqrt;

    for i in 0..dim {
        ws.d[i] = ws.d[i].abs().sqrt();
    }

    // Y5 = D^{1/2} * Sv^T, Y5_1 = Sv^T / D^{1/2}
    for i in 0..dim {
        for j in 0..dim {
            ws.y5[i][j] = ws.d[i] * ws.sv[j][i];
            ws.y5_1[i][j] = ws.sv[j][i] / ws.d[i].max(1e-300);
        }
    }

    // Sv_1 = Sv * Y5
    for i in 0..dim {
        for j in 0..dim {
            ws.sv_1[i][j] = 0.0;
            for k in 0..dim {
                ws.sv_1[i][j] += ws.sv[i][k] * ws.y5[k][j];
            }
        }
    }
    for i in 0..dim {
        for j in 0..dim {
            ws.y5[i][j] = ws.sv_1[i][j];
        }
    }

    // Y5_1 = Sv * Y5_1_old
    for i in 0..dim {
        for j in 0..dim {
            ws.sv_1[i][j] = 0.0;
            for k in 0..dim {
                ws.sv_1[i][j] += ws.sv[i][k] * ws.y5_1[k][j];
            }
        }
    }
    for i in 0..dim {
        for j in 0..dim {
            ws.y5_1[i][j] = ws.sv_1[i][j];
        }
    }

    // ZY = (Scaling_F * L + R * y) * Y5
    for i in 0..dim {
        for j in 0..dim {
            ws.zy[i][j] = 0.0;
            for k in 0..dim {
                ws.zy[i][j] += (ws.scaling_f * ws.l_m[i][k] + ws.r_m[i][k] * y) * ws.y5[k][j];
            }
        }
    }

    // Sv_1 = Y5 * ZY
    for i in 0..dim {
        for j in 0..dim {
            ws.sv_1[i][j] = 0.0;
            for k in 0..dim {
                ws.sv_1[i][j] += ws.y5[i][k] * ws.zy[k][j];
            }
        }
    }
    for i in 0..dim {
        for j in 0..dim {
            ws.zy[i][j] = ws.sv_1[i][j];
        }
    }

    // Second diagonalization (of (L+R*y)*Y5^2)
    diag(&mut ws.zy, &mut ws.sv, &mut ws.d, dim);

    // Transform eigenvectors
    for i in 0..dim {
        for j in 0..dim {
            ws.sv_1[i][j] = 0.0;
            for k in 0..dim {
                ws.sv_1[i][j] += ws.sv[k][i] * ws.y5[k][j];
            }
            ws.sv_1[i][j] *= fmin1;
        }
    }
    for i in 0..dim {
        for j in 0..dim {
            ws.zy[i][j] = 0.0;
            for k in 0..dim {
                ws.zy[i][j] += ws.y5_1[i][k] * ws.sv[k][j];
            }
            ws.zy[i][j] *= fmin_sqrt;
        }
    }
    for i in 0..dim {
        for j in 0..dim {
            ws.sv[i][j] = ws.zy[i][j];
        }
    }
}

fn eval_frequency(ws: &mut CplSetup, deg_o: usize) {
    let dim = ws.dim;
    let mut min_d = ws.d[0];
    for i in 1..dim {
        if ws.d[i] < min_d {
            min_d = ws.d[i];
        }
    }
    if min_d <= 0.0 {
        min_d = 1e-30;
    }

    ws.scaling_f2 = 1.0 / min_d;
    ws.scaling_f = ws.scaling_f2.sqrt();
    let min_val = ws.length * 8.0;

    ws.frequency[0] = 0.0;
    for i in 1..=deg_o {
        ws.frequency[i] = ws.frequency[i - 1] + min_val;
    }

    for i in 0..dim {
        ws.d[i] *= ws.scaling_f2;
    }
}

#[allow(clippy::needless_range_loop)]
fn eval_si_si_1(ws: &mut CplSetup, y: f64) {
    let dim = ws.dim;

    // Si_1 = Sv_1 * (y*R + Scaling_F*L)
    for i in 0..dim {
        for j in 0..dim {
            ws.si_1[i][j] = 0.0;
            for k in 0..dim {
                ws.si_1[i][j] += ws.sv_1[i][k] * (y * ws.r_m[k][j] + ws.scaling_f * ws.l_m[k][j]);
            }
        }
    }

    // Divide by sqrt(D)
    for i in 0..dim {
        let sd = ws.d[i].abs().sqrt().max(1e-300);
        for j in 0..dim {
            ws.si_1[i][j] /= sd;
        }
    }

    // Si = inv(Si_1) via Gauss-Jordan
    for i in 0..dim {
        for j in 0..dim {
            ws.a[i][j] = ws.si_1[i][j];
        }
        for j in dim..(2 * dim) {
            ws.a[i][j] = 0.0;
        }
        ws.a[i][i + dim] = 1.0;
    }
    gaussian_elimination_large(&mut ws.a, dim, true);

    for i in 0..dim {
        for j in 0..dim {
            ws.si[i][j] = ws.a[i][j + dim];
        }
    }
}

#[allow(clippy::needless_range_loop)]
fn store_si_sv_1(ws: &mut CplSetup, ind: usize) {
    let dim = ws.dim;
    for i in 0..dim {
        for j in 0..dim {
            let mut temp = 0.0;
            for k in 0..dim {
                temp += ws.si[i][k] * ws.sv_1[k][j];
            }
            ws.si_sv_1[i][j][ind] = temp;
        }
    }
}

#[allow(clippy::needless_range_loop)]
fn store_vals(ws: &mut CplSetup, ind: usize) {
    let dim = ws.dim;
    for i in 0..dim {
        for j in 0..dim {
            ws.sip[i][j][ind] = ws.si[i][j];
            ws.si_1p[i][j][ind] = ws.si_1[i][j];
            ws.sv_1p[i][j][ind] = ws.sv_1[i][j];
        }
        ws.w[i][ind] = ws.d[i];
    }
}

fn poly_matrix_fit(data: &mut [Vec<Vec<f64>>], dim: usize, deg: usize, freq: &[f64]) {
    for row in data.iter_mut().take(dim) {
        for col in row.iter_mut().take(dim) {
            let vals: Vec<f64> = (0..=deg).map(|k| col[k]).collect();
            let cof = poly_match(deg, freq, &vals);
            col[..=deg].copy_from_slice(&cof[..=deg]);
        }
    }
}

#[allow(clippy::needless_range_loop)]
fn matrix_p_mult_fn(
    a_in: &[Vec<Vec<f64>>],
    d1: &[Vec<f64>],
    b: &[Vec<Vec<f64>>],
    dim: usize,
    deg_o: usize,
) -> Vec<Vec<MultOut>> {
    // T[i][j] = B[i][j] * D1[i]
    let mut t_poly: Vec<Vec<Vec<f64>>> = vec![vec![vec![0.0; deg_o + 1]; dim]; dim];
    for i in 0..dim {
        for (tp_elem, b_elem) in t_poly[i].iter_mut().zip(b[i].iter()).take(dim) {
            let d_slice = &d1[i][..=deg_o.min(d1[i].len() - 1)];
            let b_slice = &b_elem[..=deg_o.min(b_elem.len() - 1)];
            *tp_elem = mult_p(
                b_slice,
                d_slice,
                b_slice.len() - 1,
                d_slice.len() - 1,
                deg_o,
            );
        }
    }

    let mut result: Vec<Vec<MultOut>> = Vec::with_capacity(dim);
    for i in 0..dim {
        let mut row = Vec::with_capacity(dim);
        for j in 0..dim {
            let mut c_0 = vec![0.0; dim];
            let mut poly = vec![None; dim];
            for k in 0..dim {
                let a_slice = &a_in[i][k][..=deg_o.min(a_in[i][k].len() - 1)];
                let t_slice = &t_poly[k][j];
                let mut p = mult_p(
                    a_slice,
                    t_slice,
                    a_slice.len() - 1,
                    t_slice.len() - 1,
                    deg_o,
                );
                // Truncate to deg_o+1
                p.truncate(deg_o + 1);
                let t1 = p[0];
                c_0[k] = t1;
                if t1 != 0.0 {
                    p[0] = 1.0;
                    for val in p.iter_mut().skip(1) {
                        *val /= t1;
                    }
                }
                poly[k] = Some(p);
            }
            row.push(MultOut { c_0, poly });
        }
        result.push(row);
    }
    result
}

fn generate_siv(ws: &CplSetup, dim: usize, deg_o: usize) -> Vec<Vec<SingleOut>> {
    let mut siv: Vec<Vec<SingleOut>> = Vec::with_capacity(dim);
    for i in 0..dim {
        let mut row = Vec::with_capacity(dim);
        for j in 0..dim {
            let mut p: Vec<f64> = (0..=deg_o).map(|k| ws.si_sv_1[i][j][k]).collect();
            let c = p[0];
            if c == 0.0 {
                row.push(SingleOut {
                    c_0: 0.0,
                    poly: None,
                });
                continue;
            }
            for v in p.iter_mut() {
                *v /= c;
            }
            let a_b = if i == j {
                (ws.g_m[i][i] / ws.r_m[i][i]).sqrt() / c
            } else {
                0.0
            };
            if let Some((c1, c2, c3, x1, x2, x3, rtv)) = pade_apx(a_b, &p) {
                row.push(SingleOut {
                    c_0: c,
                    poly: Some(vec![c1, c2, c3, x1, x2, x3, rtv as f64]),
                });
            } else {
                row.push(SingleOut {
                    c_0: 0.0,
                    poly: None,
                });
            }
        }
        siv.push(row);
    }
    siv
}

#[allow(clippy::needless_range_loop)]
fn generate_iwi_iwv(
    ws: &CplSetup,
    data: &[Vec<MultOut>],
    dim: usize,
    is_iwi: bool,
) -> Vec<Vec<MultOut>> {
    let mut result: Vec<Vec<MultOut>> = Vec::with_capacity(dim);
    for i in 0..dim {
        let mut row = Vec::with_capacity(dim);
        for j in 0..dim {
            let mut c_0 = vec![0.0; dim];
            let mut poly = vec![None; dim];
            for k in 0..dim {
                let c = data[i][j].c_0[k];
                if c == 0.0 {
                    continue;
                }
                let p = data[i][j].poly[k].as_ref().unwrap();
                let a_b = if i == j && k == i {
                    if is_iwi {
                        (-((ws.g_m[i][i] * ws.r_m[i][i]).sqrt()) * ws.length).exp() / c
                    } else {
                        (ws.g_m[i][i] / ws.r_m[i][i]).sqrt()
                            * (-((ws.g_m[i][i] * ws.r_m[i][i]).sqrt()) * ws.length).exp()
                            / c
                    }
                } else {
                    0.0
                };
                if let Some((c1, c2, c3, x1, x2, x3, rtv)) = pade_apx(a_b, p) {
                    c_0[k] = c;
                    poly[k] = Some(vec![c1, c2, c3, x1, x2, x3, rtv as f64]);
                }
            }
            row.push(MultOut { c_0, poly });
        }
        result.push(row);
    }
    result
}

// ---------------------------------------------------------------------------
// DC stamp
// ---------------------------------------------------------------------------

/// Stamp CPL DC operating point: each line acts as R*length resistance.
pub fn stamp_cpl_dc(inst: &CplInstance, system: &mut crate::LinearSystem, num_nodes: usize) {
    let no_l = inst.no_l;
    let mut resindex = 0usize;

    for m in 0..no_l {
        let ibr1 = num_nodes + inst.ibr1[m];
        let ibr2 = num_nodes + inst.ibr2[m];

        let g = inst.model.r_m[m][m] * inst.length;

        // KCL: posIbr1, negIbr2
        if let Some(p) = inst.pos_nodes[m] {
            system.matrix.add(p, ibr1, 1.0);
        }
        if let Some(n) = inst.neg_nodes[m] {
            system.matrix.add(n, ibr2, 1.0);
        }

        // ibr1 equation: ibr1 + ibr2 = 0
        system.matrix.add(ibr1, ibr1, 1.0);
        system.matrix.add(ibr1, ibr2, 1.0);

        // ibr2 equation: V_pos - V_neg - g*ibr1 = 0
        if let Some(p) = inst.pos_nodes[m] {
            system.matrix.add(ibr2, p, 1.0);
        }
        if let Some(n) = inst.neg_nodes[m] {
            system.matrix.add(ibr2, n, -1.0);
        }
        system.matrix.add(ibr2, ibr1, -g);

        resindex += no_l - m;
    }
    let _ = resindex; // suppress unused warning
}

// ---------------------------------------------------------------------------
// Transient stamp
// ---------------------------------------------------------------------------

/// Pre-computed CPL transient stamp (computed once per timestep, applied each NR iteration).
#[derive(Debug, Clone)]
pub struct CplTransientStamp {
    pub no_l: usize,
    pub ibr1: Vec<usize>,
    pub ibr2: Vec<usize>,
    pub pos_nodes: Vec<Option<usize>>,
    pub neg_nodes: Vec<Option<usize>>,
    /// Admittance matrix entries: ibr1[m]->pos[p] and ibr2[m]->neg[p].
    pub admittance: Vec<Vec<f64>>,
    /// RHS contributions for ibr1 equations.
    pub ff: Vec<f64>,
    /// RHS contributions for ibr2 equations.
    pub gg: Vec<f64>,
    /// Coupling entries: ibr1[m]->neg[p] = -f, ibr2[m]->pos[p] = -f.
    pub coupling_v: Vec<Vec<Vec<f64>>>,
    /// Coupling entries: ibr1[m]->ibr2[p] = -f, ibr2[m]->ibr1[p] = -f.
    pub coupling_i: Vec<Vec<Vec<f64>>>,
    /// Whether extended timestep coupling is active.
    pub ext: bool,
}

/// Complex exponential: exp((ar + i*ai) * h) = (e*cos(ai*h), e*sin(ai*h)).
fn exp_c(ar: f64, ai: f64, h: f64) -> (f64, f64) {
    let e = (ar * h).exp();
    let cs = (ai * h).cos();
    let si = (ai * h).sin();
    (e * cs, e * si)
}

/// Complex multiplication: (ar + i*ai) * (br + i*bi).
fn mult_c(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    (ar * br - ai * bi, ar * bi + ai * br)
}

/// Complex division: (ar + i*ai) / (br + i*bi).
fn div_c(ar: f64, ai: f64, br: f64, bi: f64) -> (f64, f64) {
    let t = br * br + bi * bi;
    ((ar * br + ai * bi) / t, (ai * br - ar * bi) / t)
}

/// Initialize CPL DC state from operating point solution.
pub fn init_dc_state(cp: &mut CpLine, dc_v_in: &[f64], dc_v_out: &[f64]) {
    let no_l = cp.no_l;

    cp.dc1[..no_l].copy_from_slice(&dc_v_in[..no_l]);
    cp.dc2[..no_l].copy_from_slice(&dc_v_out[..no_l]);

    // Initialize h1t convolutions
    for i in 0..no_l {
        for j in 0..no_l {
            if let Some(tms) = &mut cp.h1t[i][j] {
                if tms.if_img {
                    tms.tm[0].cnv_i = -dc_v_in[j] * tms.tm[0].c / tms.tm[0].x;
                    tms.tm[0].cnv_o = -dc_v_out[j] * tms.tm[0].c / tms.tm[0].x;
                    let (a, b) = div_c(tms.tm[1].c, tms.tm[2].c, tms.tm[1].x, tms.tm[2].x);
                    tms.tm[1].cnv_i = -dc_v_in[j] * a;
                    tms.tm[1].cnv_o = -dc_v_out[j] * a;
                    tms.tm[2].cnv_i = -dc_v_in[j] * b;
                    tms.tm[2].cnv_o = -dc_v_out[j] * b;
                } else {
                    for k in 0..3 {
                        if tms.tm[k].x.abs() > 1e-300 {
                            tms.tm[k].cnv_i = -dc_v_in[j] * tms.tm[k].c / tms.tm[k].x;
                            tms.tm[k].cnv_o = -dc_v_out[j] * tms.tm[k].c / tms.tm[k].x;
                        }
                    }
                }
            }

            // h2t: initialize to zero
            for l in 0..no_l {
                if let Some(tms) = &mut cp.h2t[i][j][l] {
                    for k in 0..3 {
                        tms.tm[k].cnv_i = 0.0;
                        tms.tm[k].cnv_o = 0.0;
                    }
                }
            }

            // h3t: initialize like h1t
            for l in 0..no_l {
                if let Some(tms) = &mut cp.h3t[i][j][l] {
                    if tms.if_img {
                        tms.tm[0].cnv_i = -dc_v_in[j] * tms.tm[0].c / tms.tm[0].x;
                        tms.tm[0].cnv_o = -dc_v_out[j] * tms.tm[0].c / tms.tm[0].x;
                        let (a, b) = div_c(tms.tm[1].c, tms.tm[2].c, tms.tm[1].x, tms.tm[2].x);
                        tms.tm[1].cnv_i = -dc_v_in[j] * a;
                        tms.tm[1].cnv_o = -dc_v_out[j] * a;
                        tms.tm[2].cnv_i = -dc_v_in[j] * b;
                        tms.tm[2].cnv_o = -dc_v_out[j] * b;
                    } else {
                        for k in 0..3 {
                            if tms.tm[k].x.abs() > 1e-300 {
                                tms.tm[k].cnv_i = -dc_v_in[j] * tms.tm[k].c / tms.tm[k].x;
                                tms.tm[k].cnv_o = -dc_v_out[j] * tms.tm[k].c / tms.tm[k].x;
                            }
                        }
                    }
                }
            }
        }
    }

    // Initial VI entry
    let vi = CplViEntry {
        time: 0,
        v_i: dc_v_in.to_vec(),
        v_o: dc_v_out.to_vec(),
        i_i: vec![0.0; no_l],
        i_o: vec![0.0; no_l],
    };
    cp.vi_history.clear();
    cp.vi_history.push(vi);
}

/// Get delayed voltage/current values for all modes.
#[allow(clippy::needless_range_loop)]
fn get_pvs_vi(t1: i64, t2: i64, cp: &mut CpLine) -> DelayedValues {
    let no_l = cp.no_l;
    let mk = || vec![vec![0.0; no_l]; no_l];
    let (mut v1_i, mut v2_i, mut i1_i, mut i2_i) = (mk(), mk(), mk(), mk());
    let (mut v1_o, mut v2_o, mut i1_o, mut i2_o) = (mk(), mk(), mk(), mk());
    let mut ratio = vec![0.0; no_l];
    let mut ext = false;

    let mut ta = vec![0i64; no_l];
    let mut tb = vec![0i64; no_l];
    let mut _mini = 0usize;
    let mut minta = i64::MAX;

    for i in 0..no_l {
        ta[i] = t1 - cp.taul[i] as i64;
        tb[i] = t2 - cp.taul[i] as i64;
        if ta[i] < minta {
            minta = ta[i];
            _mini = i;
        }
    }

    let hist = &cp.vi_history;

    for i in 0..no_l {
        ratio[i] = 0.0;

        if tb[i] <= 0 {
            for j in 0..no_l {
                i1_i[i][j] = 0.0;
                i2_i[i][j] = 0.0;
                i1_o[i][j] = 0.0;
                i2_o[i][j] = 0.0;
                v1_i[i][j] = cp.dc1[j];
                v2_i[i][j] = cp.dc1[j];
                v1_o[i][j] = cp.dc2[j];
                v2_o[i][j] = cp.dc2[j];
            }
        } else {
            if ta[i] <= 0 {
                for j in 0..no_l {
                    i1_i[i][j] = 0.0;
                    i1_o[i][j] = 0.0;
                    v1_i[i][j] = cp.dc1[j];
                    v1_o[i][j] = cp.dc2[j];
                }
            } else {
                // Interpolate at ta[i]
                let (idx, frac) = find_interp(hist, ta[i]);
                for j in 0..no_l {
                    v1_i[i][j] =
                        hist[idx].v_i[j] + frac * (hist[idx + 1].v_i[j] - hist[idx].v_i[j]);
                    v1_o[i][j] =
                        hist[idx].v_o[j] + frac * (hist[idx + 1].v_o[j] - hist[idx].v_o[j]);
                    i1_i[i][j] =
                        hist[idx].i_i[j] + frac * (hist[idx + 1].i_i[j] - hist[idx].i_i[j]);
                    i1_o[i][j] =
                        hist[idx].i_o[j] + frac * (hist[idx + 1].i_o[j] - hist[idx].i_o[j]);
                }
            }

            if tb[i] > t1 {
                ext = true;
                let f = (tb[i] - t1) as f64 / (t2 - t1) as f64;
                ratio[i] = f;
                let last = hist.last().unwrap();
                let f1 = 1.0 - f;
                for j in 0..no_l {
                    v2_i[i][j] = last.v_i[j] * f1;
                    v2_o[i][j] = last.v_o[j] * f1;
                    i2_i[i][j] = last.i_i[j] * f1;
                    i2_o[i][j] = last.i_o[j] * f1;
                }
            } else {
                let (idx, frac) = find_interp(hist, tb[i]);
                for j in 0..no_l {
                    v2_i[i][j] =
                        hist[idx].v_i[j] + frac * (hist[idx + 1].v_i[j] - hist[idx].v_i[j]);
                    v2_o[i][j] =
                        hist[idx].v_o[j] + frac * (hist[idx + 1].v_o[j] - hist[idx].v_o[j]);
                    i2_i[i][j] =
                        hist[idx].i_i[j] + frac * (hist[idx + 1].i_i[j] - hist[idx].i_i[j]);
                    i2_o[i][j] =
                        hist[idx].i_o[j] + frac * (hist[idx + 1].i_o[j] - hist[idx].i_o[j]);
                }
            }
        }
    }

    (ext, v1_i, v2_i, i1_i, i2_i, v1_o, v2_o, i1_o, i2_o, ratio)
}

fn find_interp(hist: &[CplViEntry], target: i64) -> (usize, f64) {
    if hist.len() < 2 {
        return (0, 0.0);
    }
    for idx in 0..(hist.len() - 1) {
        if hist[idx + 1].time >= target {
            let dt = (hist[idx + 1].time - hist[idx].time) as f64;
            if dt.abs() < 1e-300 {
                return (idx, 0.0);
            }
            let f = (target - hist[idx].time) as f64 / dt;
            return (idx, f);
        }
    }
    (hist.len() - 2, 1.0)
}

/// Update h1 convolutions for a CPL after recording new VI.
pub fn update_cnv_cpl(cp: &mut CpLine, h: f64) {
    let no_l = cp.no_l;
    let tail = cp.vi_history.last().unwrap();
    let prev = if cp.vi_history.len() >= 2 {
        &cp.vi_history[cp.vi_history.len() - 2]
    } else {
        tail
    };

    for j in 0..no_l {
        for k in 0..no_l {
            let ai = tail.v_i[k];
            let ao = tail.v_o[k];
            let delta = if cp.vi_history.len() >= 2 {
                (tail.time - prev.time) as f64
            } else {
                1.0
            };
            let bi = if delta.abs() > 1e-300 {
                (tail.v_i[k] - prev.v_i[k]) / delta
            } else {
                0.0
            };
            let bo = if delta.abs() > 1e-300 {
                (tail.v_o[k] - prev.v_o[k]) / delta
            } else {
                0.0
            };

            if let Some(tms) = &mut cp.h1t[j][k] {
                if tms.if_img {
                    let e = cp.h1e[j][k][0];
                    let er = cp.h1e[j][k][1];
                    let ei = cp.h1e[j][k][2];

                    // Update complex pair
                    update_cnv_a_cpl(tms, h, ai, ao, ai - bi * h, ao - bo * h, er, ei);

                    // Update real term
                    let tm = &mut tms.tm[0];
                    let t = tm.c / tm.x;
                    let bic = bi * t;
                    let boc = bo * t;
                    tm.cnv_i =
                        (tm.cnv_i - bic * h) * e + (e - 1.0) * (ai * t + 1.0e12 * bic / tm.x);
                    tm.cnv_o =
                        (tm.cnv_o - boc * h) * e + (e - 1.0) * (ao * t + 1.0e12 * boc / tm.x);
                } else {
                    for idx in 0..3 {
                        let e = cp.h1e[j][k][idx];
                        let tm = &mut tms.tm[idx];
                        let t = tm.c / tm.x;
                        let bic = bi * t;
                        let boc = bo * t;
                        tm.cnv_i =
                            (tm.cnv_i - bic * h) * e + (e - 1.0) * (ai * t + 1.0e12 * bic / tm.x);
                        tm.cnv_o =
                            (tm.cnv_o - boc * h) * e + (e - 1.0) * (ao * t + 1.0e12 * boc / tm.x);
                    }
                }
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn update_cnv_a_cpl(tms: &mut Tms, h: f64, ai: f64, ao: f64, bi: f64, bo: f64, er: f64, ei: f64) {
    let h2 = h * 0.5e-12;
    let (a1, b1) = mult_c(tms.tm[1].c, tms.tm[2].c, er, ei);
    let (a, b) = mult_c(tms.tm[1].cnv_i, tms.tm[2].cnv_i, er, ei);
    tms.tm[1].cnv_i = a + h2 * (a1 * bi + ai * tms.tm[1].c);
    tms.tm[2].cnv_i = b + h2 * (b1 * bi + ai * tms.tm[2].c);

    let (a, b) = mult_c(tms.tm[1].cnv_o, tms.tm[2].cnv_o, er, ei);
    tms.tm[1].cnv_o = a + h2 * (a1 * bo + ao * tms.tm[1].c);
    tms.tm[2].cnv_o = b + h2 * (b1 * bo + ao * tms.tm[2].c);
}

/// Update delayed convolutions for h2t/h3t when extended timestep is active.
pub fn update_delayed_cnv_cpl(cp: &mut CpLine, h: f64) {
    let no_l = cp.no_l;
    let h2 = h * 0.5e-12;
    let vi = cp.vi_history.last().unwrap().clone();

    for k in 0..no_l {
        if cp.ratio[k] <= 0.0 {
            continue;
        }
        let r = cp.ratio[k];
        for i in 0..no_l {
            for j in 0..no_l {
                if let Some(tms) = &mut cp.h3t[i][j][k] {
                    let f_i = h2 * r * vi.v_i[j];
                    let f_o = h2 * r * vi.v_o[j];
                    for t in 0..3 {
                        tms.tm[t].cnv_i += f_i * tms.tm[t].c;
                        tms.tm[t].cnv_o += f_o * tms.tm[t].c;
                    }
                }

                if let Some(tms) = &mut cp.h2t[i][j][k] {
                    let f_i = h2 * r * vi.i_i[j];
                    let f_o = h2 * r * vi.i_o[j];
                    for t in 0..3 {
                        tms.tm[t].cnv_i += f_i * tms.tm[t].c;
                        tms.tm[t].cnv_o += f_o * tms.tm[t].c;
                    }
                }
            }
        }
    }
}

/// Prepare CPL transient stamp (once per timestep).
#[allow(clippy::needless_range_loop)]
pub fn prepare_cpl_transient(
    inst: &mut CplInstance,
    num_nodes: usize,
    solution: &[f64],
    time_ps: i64,
    prev_time_ps: i64,
    h_seconds: f64,
) -> CplTransientStamp {
    let no_l = inst.no_l;
    let cp = &mut inst.cpline;
    let h = h_seconds;
    let h1 = 0.5 * h;

    // Record new VI entry
    let time = prev_time_ps;
    if cp.vi_history.is_empty() || time > cp.vi_history.last().unwrap().time {
        let mut vi = CplViEntry {
            time,
            v_i: vec![0.0; no_l],
            v_o: vec![0.0; no_l],
            i_i: vec![0.0; no_l],
            i_o: vec![0.0; no_l],
        };
        for m in 0..no_l {
            vi.v_i[m] = inst.pos_nodes[m].map_or(0.0, |i| solution[i]);
            vi.v_o[m] = inst.neg_nodes[m].map_or(0.0, |i| solution[i]);
            vi.i_i[m] = solution[num_nodes + inst.ibr1[m]];
            vi.i_o[m] = solution[num_nodes + inst.ibr2[m]];
        }

        // Update h1 cnv before adding new point
        if cp.vi_history.len() >= 2 {
            let prev_tail = cp.vi_history[cp.vi_history.len() - 1].clone();
            let prev_prev = cp.vi_history[cp.vi_history.len() - 2].clone();
            let delta = (prev_tail.time - prev_prev.time) as f64;
            if delta > 0.0 {
                // Compute dv for each node
                // update_cnv is called with delta (not h)
                update_cnv_cpl(cp, delta);
                if cp.ext {
                    update_delayed_cnv_cpl(cp, delta);
                }
            }
        }

        cp.vi_history.push(vi);
    }

    // Compute admittance matrix: h1 terms
    let mut ff = vec![0.0; no_l];
    let mut gg = vec![0.0; no_l];

    let v_in: Vec<f64> = (0..no_l)
        .map(|m| inst.pos_nodes[m].map_or(0.0, |i| solution[i]))
        .collect();
    let v_out: Vec<f64> = (0..no_l)
        .map(|m| inst.neg_nodes[m].map_or(0.0, |i| solution[i]))
        .collect();

    for j in 0..no_l {
        for k in 0..no_l {
            if let Some(tms) = &mut cp.h1t[j][k] {
                if tms.if_img {
                    let e = (tms.tm[0].x * h).exp();
                    cp.h1e[j][k][0] = e;
                    let (er, ei) = exp_c(tms.tm[1].x, tms.tm[2].x, h);
                    cp.h1e[j][k][1] = er;
                    cp.h1e[j][k][2] = ei;

                    let ff1 = tms.tm[0].c * e * h1;
                    ff[j] -= tms.tm[0].cnv_i * e;
                    gg[j] -= tms.tm[0].cnv_o * e;
                    ff[j] -= ff1 * v_in[k];
                    gg[j] -= ff1 * v_out[k];

                    let (a1, _b1) = mult_c(tms.tm[1].c, tms.tm[2].c, er, ei);
                    let (a, _b) = mult_c(tms.tm[1].cnv_i, tms.tm[2].cnv_i, er, ei);
                    ff[j] -= 2.0 * (a1 * h1 * v_in[k] + a);
                    let (a, _b) = mult_c(tms.tm[1].cnv_o, tms.tm[2].cnv_o, er, ei);
                    gg[j] -= 2.0 * (a1 * h1 * v_out[k] + a);
                } else {
                    let mut ff1 = 0.0;
                    for idx in 0..3 {
                        let e = (tms.tm[idx].x * h).exp();
                        cp.h1e[j][k][idx] = e;
                        ff1 -= tms.tm[idx].c * e;
                        ff[j] -= tms.tm[idx].cnv_i * e;
                        gg[j] -= tms.tm[idx].cnv_o * e;
                    }
                    ff[j] += ff1 * h1 * v_in[k];
                    gg[j] += ff1 * h1 * v_out[k];
                }
            }
        }
    }

    // Get delayed values
    let (ext, v1_i, v2_i, i1_i, i2_i, v1_o, v2_o, i1_o, i2_o, ratio_vals) =
        get_pvs_vi(prev_time_ps, time_ps, cp);

    // h3t contributions (voltage coupling)
    for j in 0..no_l {
        for k in 0..no_l {
            for l in 0..no_l {
                if let Some(tms) = &mut cp.h3t[j][k][l] {
                    if tms.if_img {
                        let (er, ei) = exp_c(tms.tm[1].x, tms.tm[2].x, h);
                        let a2 = h1 * tms.tm[1].c;
                        let b2 = h1 * tms.tm[2].c;

                        let (a, b) = mult_c(tms.tm[1].cnv_i, tms.tm[2].cnv_i, er, ei);
                        let (a1, b1) =
                            mult_c(a2, b2, v1_i[l][k] * er + v2_i[l][k], v1_i[l][k] * ei);
                        tms.tm[1].cnv_i = a + a1;
                        tms.tm[2].cnv_i = b + b1;

                        let (a, b) = mult_c(tms.tm[1].cnv_o, tms.tm[2].cnv_o, er, ei);
                        let (a1, b1) =
                            mult_c(a2, b2, v1_o[l][k] * er + v2_o[l][k], v1_o[l][k] * ei);
                        tms.tm[1].cnv_o = a + a1;
                        tms.tm[2].cnv_o = b + b1;

                        let tm = &mut tms.tm[0];
                        let e = (tm.x * h).exp();
                        tm.cnv_i = tm.cnv_i * e + h1 * tm.c * (v1_i[l][k] * e + v2_i[l][k]);
                        tm.cnv_o = tm.cnv_o * e + h1 * tm.c * (v1_o[l][k] * e + v2_o[l][k]);

                        ff[j] += tms.aten * v2_o[l][k] + tms.tm[0].cnv_o + 2.0 * tms.tm[1].cnv_o;
                        gg[j] += tms.aten * v2_i[l][k] + tms.tm[0].cnv_i + 2.0 * tms.tm[1].cnv_i;
                    } else {
                        for idx in 0..3 {
                            let tm = &mut tms.tm[idx];
                            let e = (tm.x * h).exp();
                            tm.cnv_i = tm.cnv_i * e + h1 * tm.c * (v1_i[l][k] * e + v2_i[l][k]);
                            tm.cnv_o = tm.cnv_o * e + h1 * tm.c * (v1_o[l][k] * e + v2_o[l][k]);
                            ff[j] += tm.cnv_o;
                            gg[j] += tm.cnv_i;
                        }
                        ff[j] += tms.aten * v2_o[l][k];
                        gg[j] += tms.aten * v2_i[l][k];
                    }
                }
            }
        }

        // h2t contributions (current coupling)
        for k in 0..no_l {
            for l in 0..no_l {
                if let Some(tms) = &mut cp.h2t[j][k][l] {
                    if tms.if_img {
                        let (er, ei) = exp_c(tms.tm[1].x, tms.tm[2].x, h);
                        let a2 = h1 * tms.tm[1].c;
                        let b2 = h1 * tms.tm[2].c;

                        let (a, b) = mult_c(tms.tm[1].cnv_i, tms.tm[2].cnv_i, er, ei);
                        let (a1, b1) =
                            mult_c(a2, b2, i1_i[l][k] * er + i2_i[l][k], i1_i[l][k] * ei);
                        tms.tm[1].cnv_i = a + a1;
                        tms.tm[2].cnv_i = b + b1;

                        let (a, b) = mult_c(tms.tm[1].cnv_o, tms.tm[2].cnv_o, er, ei);
                        let (a1, b1) =
                            mult_c(a2, b2, i1_o[l][k] * er + i2_o[l][k], i1_o[l][k] * ei);
                        tms.tm[1].cnv_o = a + a1;
                        tms.tm[2].cnv_o = b + b1;

                        let tm = &mut tms.tm[0];
                        let e = (tm.x * h).exp();
                        tm.cnv_i = tm.cnv_i * e + h1 * tm.c * (i1_i[l][k] * e + i2_i[l][k]);
                        tm.cnv_o = tm.cnv_o * e + h1 * tm.c * (i1_o[l][k] * e + i2_o[l][k]);

                        ff[j] += tms.aten * i2_o[l][k] + tms.tm[0].cnv_o + 2.0 * tms.tm[1].cnv_o;
                        gg[j] += tms.aten * i2_i[l][k] + tms.tm[0].cnv_i + 2.0 * tms.tm[1].cnv_i;
                    } else {
                        for idx in 0..3 {
                            let tm = &mut tms.tm[idx];
                            let e = (tm.x * h).exp();
                            tm.cnv_i = tm.cnv_i * e + h1 * tm.c * (i1_i[l][k] * e + i2_i[l][k]);
                            tm.cnv_o = tm.cnv_o * e + h1 * tm.c * (i1_o[l][k] * e + i2_o[l][k]);
                            ff[j] += tm.cnv_o;
                            gg[j] += tm.cnv_i;
                        }
                        ff[j] += tms.aten * i2_o[l][k];
                        gg[j] += tms.aten * i2_i[l][k];
                    }
                }
            }
        }
    }

    // Build admittance matrix
    let mut admittance = vec![vec![0.0; no_l]; no_l];
    for (m, adm_row) in admittance.iter_mut().enumerate().take(no_l) {
        for (p, adm_val) in adm_row.iter_mut().enumerate().take(no_l) {
            if let Some(tms) = &cp.h1t[m][p] {
                *adm_val = tms.aten + h1 * cp.h1c_coeff[m][p];
            }
        }
    }

    // Build coupling matrices for extended timestep
    let mut coupling_v = vec![vec![vec![0.0; no_l]; no_l]; no_l];
    let mut coupling_i = vec![vec![vec![0.0; no_l]; no_l]; no_l];

    if ext {
        for q in 0..no_l {
            cp.ratio[q] = ratio_vals[q];
            if ratio_vals[q] > 0.0 {
                for m in 0..no_l {
                    for p in 0..no_l {
                        if let Some(tms) = &cp.h3t[m][p][q] {
                            coupling_v[m][p][q] =
                                ratio_vals[q] * (h1 * cp.h3c_coeff[m][p][q] + tms.aten);
                        }
                        if let Some(tms) = &cp.h2t[m][p][q] {
                            coupling_i[m][p][q] =
                                ratio_vals[q] * (h1 * cp.h2c_coeff[m][p][q] + tms.aten);
                        }
                    }
                }
            }
        }
    }
    cp.ext = ext;

    CplTransientStamp {
        no_l,
        ibr1: inst.ibr1.clone(),
        ibr2: inst.ibr2.clone(),
        pos_nodes: inst.pos_nodes.clone(),
        neg_nodes: inst.neg_nodes.clone(),
        admittance,
        ff,
        gg,
        coupling_v,
        coupling_i,
        ext,
    }
}

/// Apply CPL transient stamp to the linear system.
pub fn apply_cpl_transient(
    stamp: &CplTransientStamp,
    system: &mut crate::LinearSystem,
    num_nodes: usize,
) {
    let no_l = stamp.no_l;

    for m in 0..no_l {
        let ibr1 = num_nodes + stamp.ibr1[m];
        let ibr2 = num_nodes + stamp.ibr2[m];

        // ibr1 = -1, ibr2 = -1 (transient: opposite sign from DC)
        system.matrix.add(ibr1, ibr1, -1.0);
        system.matrix.add(ibr2, ibr2, -1.0);

        // KCL connections
        if let Some(p) = stamp.pos_nodes[m] {
            system.matrix.add(p, ibr1, 1.0);
        }
        if let Some(n) = stamp.neg_nodes[m] {
            system.matrix.add(n, ibr2, 1.0);
        }

        // Admittance matrix entries
        for p in 0..no_l {
            if stamp.admittance[m][p] != 0.0 {
                let pos_p = stamp.pos_nodes[p];
                let neg_p = stamp.neg_nodes[p];
                if let Some(pp) = pos_p {
                    system.matrix.add(ibr1, pp, stamp.admittance[m][p]);
                }
                if let Some(np) = neg_p {
                    system.matrix.add(ibr2, np, stamp.admittance[m][p]);
                }
            }
        }

        // Extended timestep coupling
        if stamp.ext {
            for q in 0..no_l {
                for p in 0..no_l {
                    let f_v = stamp.coupling_v[m][p][q];
                    if f_v != 0.0 {
                        if let Some(np) = stamp.neg_nodes[p] {
                            system.matrix.add(ibr1, np, -f_v);
                        }
                        if let Some(pp) = stamp.pos_nodes[p] {
                            system.matrix.add(ibr2, pp, -f_v);
                        }
                    }
                    let f_i = stamp.coupling_i[m][p][q];
                    if f_i != 0.0 {
                        let ibr2_p = num_nodes + stamp.ibr2[p];
                        let ibr1_p = num_nodes + stamp.ibr1[p];
                        system.matrix.add(ibr1, ibr2_p, -f_i);
                        system.matrix.add(ibr2, ibr1_p, -f_i);
                    }
                }
            }
        }

        // RHS
        system.rhs[ibr1] += stamp.ff[m];
        system.rhs[ibr2] += stamp.gg[m];
    }
}

/// Copy CPL line state (for rollback during rejected timesteps).
pub fn copy_cpline(dst: &mut CpLine, src: &CpLine) {
    dst.no_l = src.no_l;
    dst.ext = src.ext;
    dst.ratio = src.ratio.clone();
    dst.taul = src.taul.clone();

    for i in 0..src.no_l {
        for j in 0..src.no_l {
            dst.h1c_coeff[i][j] = src.h1c_coeff[i][j];
            dst.h1e[i][j] = src.h1e[i][j];
            if let Some(ref s) = src.h1t[i][j] {
                if dst.h1t[i][j].is_none() {
                    dst.h1t[i][j] = Some(Tms::new());
                }
                let d = dst.h1t[i][j].as_mut().unwrap();
                d.if_img = s.if_img;
                d.aten = s.aten;
                for k in 0..3 {
                    d.tm[k] = s.tm[k].clone();
                }
            }
            for l in 0..src.no_l {
                dst.h2c_coeff[i][j][l] = src.h2c_coeff[i][j][l];
                dst.h3c_coeff[i][j][l] = src.h3c_coeff[i][j][l];
                if let Some(ref s) = src.h2t[i][j][l] {
                    if dst.h2t[i][j][l].is_none() {
                        dst.h2t[i][j][l] = Some(Tms::new());
                    }
                    let d = dst.h2t[i][j][l].as_mut().unwrap();
                    d.if_img = s.if_img;
                    d.aten = s.aten;
                    for k in 0..3 {
                        d.tm[k] = s.tm[k].clone();
                    }
                }
                if let Some(ref s) = src.h3t[i][j][l] {
                    if dst.h3t[i][j][l].is_none() {
                        dst.h3t[i][j][l] = Some(Tms::new());
                    }
                    let d = dst.h3t[i][j][l].as_mut().unwrap();
                    d.if_img = s.if_img;
                    d.aten = s.aten;
                    for k in 0..3 {
                        d.tm[k] = s.tm[k].clone();
                    }
                }
            }
        }
    }

    dst.dc1 = src.dc1.clone();
    dst.dc2 = src.dc2.clone();
    dst.vi_history = src.vi_history.clone();
}
