//! Lossy Transmission Line (LTRA) model.
//!
//! Implements the LTRA device from ngspice: a 2-port distributed RLC/RC/LC
//! transmission line using convolution-based transient analysis.
//!
//! Reference: ngspice `src/spicelib/devices/ltra/`

use std::f64::consts::PI;
use thevenin_types::ModelDef;

// ---------------------------------------------------------------------------
// Special case types
// ---------------------------------------------------------------------------

/// Classification of the transmission line based on which per-unit-length
/// parameters are nonzero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LtraCase {
    /// Lossless: R=0, G=0, L≠0, C≠0
    Lc,
    /// Lossy with both R and L: R≠0, G=0, L≠0, C≠0
    Rlc,
    /// RC (diffusive): R≠0, G=0, L=0, C≠0
    Rc,
    /// RG (DC only): R≠0, G≠0, L=0, C=0
    Rg,
}

// ---------------------------------------------------------------------------
// Model
// ---------------------------------------------------------------------------

/// LTRA model parameters and precomputed values.
#[derive(Debug, Clone)]
pub struct LtraModel {
    // Per-unit-length input parameters
    pub r: f64,
    pub l: f64,
    pub g: f64,
    pub c: f64,
    pub length: f64,

    // Calculated parameters
    pub td: f64,
    pub imped: f64,
    pub admit: f64,
    pub alpha: f64,
    pub beta: f64,
    pub attenuation: f64,

    // RC-specific
    pub c_by_r: f64,
    pub rclsqr: f64,

    // Integral values of impulse responses (for initial conditions)
    pub int_h1dash: f64,
    pub int_h2: f64,
    pub int_h3dash: f64,

    // RG DC parameters
    pub coshl_root_gr: f64,
    pub r_rs_l_rgr_or_g: f64,
    pub r_gs_l_rgr_or_r: f64,

    // Control parameters
    pub reltol: f64,
    pub abstol: f64,
    pub st_line_reltol: f64,
    pub st_line_abstol: f64,
    pub chop_reltol: f64,
    pub chop_abstol: f64,
    pub step_limit: bool,
    pub special_case: LtraCase,
}

impl LtraModel {
    pub fn from_model_def(def: &ModelDef) -> Self {
        let mut r = 0.0;
        let mut l = 0.0;
        let mut g = 0.0;
        let mut c = 0.0;
        let mut length = 0.0;
        let mut reltol = 1.0;
        let mut abstol = 1.0;
        let mut step_limit = false;
        let mut chop_reltol = 0.0;
        let mut chop_abstol = 0.0;

        let mut r_given = false;
        let mut l_given = false;
        let mut g_given = false;
        let mut c_given = false;
        let mut length_given = false;

        for p in &def.params {
            let val = crate::expr_val_or(&p.value, 0.0);
            match p.name.to_uppercase().as_str() {
                "R" => {
                    r = val;
                    r_given = true;
                }
                "L" => {
                    l = val;
                    l_given = true;
                }
                "G" => {
                    g = val;
                    g_given = true;
                }
                "C" => {
                    c = val;
                    c_given = true;
                }
                "LEN" => {
                    length = val;
                    length_given = true;
                }
                "REL" | "RELTOL" => reltol = val,
                "ABS" | "ABSTOL" => abstol = val,
                "COMPACTREL" => chop_reltol = val,
                "COMPACTABS" => chop_abstol = val,
                "STEPLIMIT" => step_limit = true,
                "NOSTEPLIMIT" => step_limit = false,
                _ => {}
            }
        }

        if !r_given {
            r = 0.0;
        }
        if !g_given {
            g = 0.0;
        }
        if !l_given {
            l = 0.0;
        }
        if !c_given {
            c = 0.0;
        }
        assert!(length_given, "LTRA model requires LEN parameter");

        // Determine special case
        let special_case = classify_line(r, l, g, c);

        let mut model = LtraModel {
            r,
            l,
            g,
            c,
            length,
            td: 0.0,
            imped: 0.0,
            admit: 0.0,
            alpha: 0.0,
            beta: 0.0,
            attenuation: 1.0,
            c_by_r: 0.0,
            rclsqr: 0.0,
            int_h1dash: 0.0,
            int_h2: 0.0,
            int_h3dash: 0.0,
            coshl_root_gr: 0.0,
            r_rs_l_rgr_or_g: 0.0,
            r_gs_l_rgr_or_r: 0.0,
            reltol,
            abstol,
            st_line_reltol: 0.0,
            st_line_abstol: 0.0,
            chop_reltol,
            chop_abstol,
            step_limit,
            special_case,
        };

        model.precompute();
        model
    }

    /// Precompute derived values (equivalent to LTRAtemp).
    fn precompute(&mut self) {
        match self.special_case {
            LtraCase::Lc => {
                self.imped = (self.l / self.c).sqrt();
                self.admit = 1.0 / self.imped;
                self.td = (self.l * self.c).sqrt() * self.length;
                self.attenuation = 1.0;
            }
            LtraCase::Rlc => {
                self.imped = (self.l / self.c).sqrt();
                self.admit = 1.0 / self.imped;
                self.td = (self.l * self.c).sqrt() * self.length;
                self.alpha = 0.5 * (self.r / self.l);
                self.beta = self.alpha;
                self.attenuation = (-self.beta * self.td).exp();

                if self.alpha > 0.0 {
                    self.int_h1dash = -1.0;
                    self.int_h2 = 1.0 - self.attenuation;
                    self.int_h3dash = -self.attenuation;
                }
            }
            LtraCase::Rc => {
                self.c_by_r = self.c / self.r;
                self.rclsqr = self.r * self.c * self.length * self.length;
                self.int_h1dash = 0.0;
                self.int_h2 = 1.0;
                self.int_h3dash = 0.0;
            }
            LtraCase::Rg => {
                // DC only
            }
        }
    }
}

fn classify_line(r: f64, l: f64, g: f64, c: f64) -> LtraCase {
    if r == 0.0 && g == 0.0 && l != 0.0 && c != 0.0 {
        LtraCase::Lc
    } else if r != 0.0 && g == 0.0 && l != 0.0 && c != 0.0 {
        LtraCase::Rlc
    } else if r != 0.0 && g == 0.0 && l == 0.0 && c != 0.0 {
        LtraCase::Rc
    } else if r != 0.0 && g != 0.0 && l == 0.0 && c == 0.0 {
        LtraCase::Rg
    } else {
        panic!("LTRA: unsupported line type (need at least 2 of R,L,G,C nonzero, RL not supported)")
    }
}

// ---------------------------------------------------------------------------
// Instance
// ---------------------------------------------------------------------------

/// A resolved LTRA instance with matrix indices.
#[derive(Debug, Clone)]
pub struct LtraInstance {
    pub name: String,
    pub pos1_idx: Option<usize>,
    pub neg1_idx: Option<usize>,
    pub pos2_idx: Option<usize>,
    pub neg2_idx: Option<usize>,
    /// Branch equation index for port 1 (index into vsource array).
    pub br_eq1: usize,
    /// Branch equation index for port 2.
    pub br_eq2: usize,
    pub model: LtraModel,
}

// ---------------------------------------------------------------------------
// Transient state per instance
// ---------------------------------------------------------------------------

/// History storage for LTRA transient simulation.
#[derive(Debug, Clone)]
pub struct LtraState {
    /// Past voltage values at port 1.
    pub v1: Vec<f64>,
    /// Past current values at port 1.
    pub i1: Vec<f64>,
    /// Past voltage values at port 2.
    pub v2: Vec<f64>,
    /// Past current values at port 2.
    pub i2: Vec<f64>,
    /// Initial voltage at port 1 (from DC OP).
    pub init_v1: f64,
    /// Initial current at port 1.
    pub init_i1: f64,
    /// Initial voltage at port 2.
    pub init_v2: f64,
    /// Initial current at port 2.
    pub init_i2: f64,
}

impl Default for LtraState {
    fn default() -> Self {
        Self {
            v1: Vec::new(),
            i1: Vec::new(),
            v2: Vec::new(),
            i2: Vec::new(),
            init_v1: 0.0,
            init_i1: 0.0,
            init_v2: 0.0,
            init_i2: 0.0,
        }
    }
}

impl LtraState {
    pub fn new() -> Self {
        Self::default()
    }

    /// Initialize from DC operating point solution.
    pub fn init_from_dc(&mut self, v1: f64, i1: f64, v2: f64, i2: f64) {
        self.init_v1 = v1;
        self.init_i1 = i1;
        self.init_v2 = v2;
        self.init_i2 = i2;
    }

    /// Record values at current timepoint.
    pub fn accept(&mut self, v1: f64, i1: f64, v2: f64, i2: f64) {
        self.v1.push(v1);
        self.i1.push(i1);
        self.v2.push(v2);
        self.i2.push(i2);
    }
}

// ---------------------------------------------------------------------------
// Coefficient state per model (shared across instances of same model)
// ---------------------------------------------------------------------------

/// Convolution coefficients for impulse response functions.
#[derive(Debug, Clone)]
pub struct LtraCoeffs {
    pub h1dash_first_coeff: f64,
    pub h2_first_coeff: f64,
    pub h3dash_first_coeff: f64,
    pub h1dash_coeffs: Vec<f64>,
    pub h2_coeffs: Vec<f64>,
    pub h3dash_coeffs: Vec<f64>,
    /// Auxiliary index for h2/h3dash (marks where delayed signal starts).
    pub aux_index: usize,
}

impl Default for LtraCoeffs {
    fn default() -> Self {
        Self {
            h1dash_first_coeff: 0.0,
            h2_first_coeff: 0.0,
            h3dash_first_coeff: 0.0,
            h1dash_coeffs: Vec::new(),
            h2_coeffs: Vec::new(),
            h3dash_coeffs: Vec::new(),
            aux_index: 0,
        }
    }
}

impl LtraCoeffs {
    pub fn new() -> Self {
        Self::default()
    }

    /// Ensure coefficient arrays are large enough.
    pub fn ensure_size(&mut self, size: usize) {
        if self.h1dash_coeffs.len() < size {
            self.h1dash_coeffs.resize(size, 0.0);
            self.h2_coeffs.resize(size, 0.0);
            self.h3dash_coeffs.resize(size, 0.0);
        }
    }
}

// ---------------------------------------------------------------------------
// Bessel functions (from Numerical Recipes)
// ---------------------------------------------------------------------------

fn bess_i0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 3.75 {
        let y = (x / 3.75).powi(2);
        1.0 + y
            * (3.5156229
                + y * (3.0899424
                    + y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))))
    } else {
        let y = 3.75 / ax;
        (ax.exp() / ax.sqrt())
            * (0.39894228
                + y * (0.1328592e-1
                    + y * (0.225319e-2
                        + y * (-0.157565e-2
                            + y * (0.916281e-2
                                + y * (-0.2057706e-1
                                    + y * (0.2635537e-1
                                        + y * (-0.1647633e-1 + y * 0.392377e-2))))))))
    }
}

fn bess_i1(x: f64) -> f64 {
    let ax = x.abs();
    let ans = if ax < 3.75 {
        let y = (x / 3.75).powi(2);
        ax * (0.5
            + y * (0.87890594
                + y * (0.51498869
                    + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))))
    } else {
        let y = 3.75 / ax;
        let mut ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
        ans = 0.39894228
            + y * (-0.3988024e-1
                + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
        ans * (ax.exp() / ax.sqrt())
    };
    if x < 0.0 { -ans } else { ans }
}

/// I1(x)/x — avoids division by zero when x → 0.
fn bess_i1x_over_x(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 3.75 {
        let y = (x / 3.75).powi(2);
        0.5 + y
            * (0.87890594
                + y * (0.51498869
                    + y * (0.15084934 + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3)))))
    } else {
        let y = 3.75 / ax;
        let mut ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
        ans = 0.39894228
            + y * (-0.3988024e-1
                + y * (-0.362018e-2 + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
        ans * (ax.exp() / (ax * ax.sqrt()))
    }
}

// ---------------------------------------------------------------------------
// Impulse response functions (RLC case)
// ---------------------------------------------------------------------------

/// h1'(t) = α·e^{-β·t}·(I₁(α·t) - I₀(α·t))
#[allow(dead_code)]
fn rlc_h1dash_func(time: f64, alpha: f64, beta: f64) -> f64 {
    if alpha == 0.0 {
        return 0.0;
    }
    let exp_arg = -beta * time;
    let bessel_arg = alpha * time;
    (bess_i1(bessel_arg) - bess_i0(bessel_arg)) * alpha * exp_arg.exp()
}

/// ∫∫ h1'(τ) dτ dτ (twice-integrated h1dash)
#[allow(dead_code)]
fn rlc_h1dash_twice_int_func(time: f64, beta: f64) -> f64 {
    if beta == 0.0 {
        return time;
    }
    let arg = beta * time;
    if arg == 0.0 {
        return 0.0;
    }
    (bess_i1(arg) + bess_i0(arg)) * time * (-arg).exp() - time
}

/// h2(t) for RLC: involves I₁(α·√(t²-T²))/√(t²-T²)
fn rlc_h2_func(time: f64, td: f64, alpha: f64, beta: f64) -> f64 {
    if alpha == 0.0 || time < td {
        return 0.0;
    }
    let bessel_arg = if time != td {
        alpha * (time * time - td * td).sqrt()
    } else {
        0.0
    };
    alpha * alpha * td * (-beta * time).exp() * bess_i1x_over_x(bessel_arg)
}

/// ∫ h3'(τ) dτ (integrated h3dash)
#[allow(dead_code)]
fn rlc_h3dash_int_func(time: f64, td: f64, beta: f64) -> f64 {
    if time <= td || beta == 0.0 {
        return 0.0;
    }
    let bessel_arg = beta * (time * time - td * td).sqrt();
    (-beta * time).exp() * bess_i0(bessel_arg) - (-beta * td).exp()
}

// ---------------------------------------------------------------------------
// RC impulse response twice-integrated functions
// ---------------------------------------------------------------------------

fn rc_h1dash_twice_int_func(time: f64, cbyr: f64) -> f64 {
    (4.0 * cbyr * time / PI).sqrt()
}

#[allow(dead_code)]
fn rc_h2_twice_int_func(time: f64, rclsqr: f64) -> f64 {
    if time == 0.0 {
        return 0.0;
    }
    let temp = rclsqr / (4.0 * time);
    (time + rclsqr * 0.5) * erfc(temp.sqrt()) - (time * rclsqr / PI).sqrt() * (-temp).exp()
}

#[allow(dead_code)]
fn rc_h3dash_twice_int_func(time: f64, cbyr: f64, rclsqr: f64) -> f64 {
    if time == 0.0 {
        return 0.0;
    }
    let temp = rclsqr / (4.0 * time);
    let val = 2.0 * (time / PI).sqrt() * (-temp).exp() - rclsqr.sqrt() * erfc(temp.sqrt());
    cbyr.sqrt() * val
}

// ---------------------------------------------------------------------------
// Integration helper functions (from ltramisc.c)
// ---------------------------------------------------------------------------

/// ∫_{lo}^{hi} h(τ) dτ where h is linear between (t1,lo_val) and (t2,hi_val).
fn int_lin_func(lo: f64, hi: f64, lo_val: f64, hi_val: f64, t1: f64, t2: f64) -> f64 {
    let width = t2 - t1;
    if width == 0.0 {
        return 0.0;
    }
    let m = (hi_val - lo_val) / width;
    (hi - lo) * lo_val + 0.5 * m * ((hi - t1).powi(2) - (lo - t1).powi(2))
}

/// Double integral of linear function.
fn twice_int_lin_func(
    lo: f64,
    hi: f64,
    other_lo: f64,
    lo_val: f64,
    hi_val: f64,
    t1: f64,
    t2: f64,
) -> f64 {
    let width = t2 - t1;
    if width == 0.0 {
        return 0.0;
    }
    let m = (hi_val - lo_val) / width;
    let temp1 = hi - t1;
    let temp2 = lo - t1;
    let temp3 = other_lo - t1;
    let mut dummy = lo_val * ((hi - other_lo).powi(2) - (lo - other_lo).powi(2));
    dummy += m * ((temp1.powi(3) - temp2.powi(3)) / 3.0 - temp3.powi(2) * (hi - lo));
    dummy * 0.5
}

// ---------------------------------------------------------------------------
// Coefficient setup
// ---------------------------------------------------------------------------

/// Set up convolution coefficients for RC line.
#[expect(unused_assignments)]
pub fn rc_coeffs_setup(
    coeffs: &mut LtraCoeffs,
    cbyr: f64,
    rclsqr: f64,
    cur_time: f64,
    time_points: &[f64],
    time_index: usize,
    reltol: f64,
) {
    coeffs.ensure_size(time_index + 1);
    let aux_index = time_index;

    let delta1_init = cur_time - time_points[aux_index];
    let hilimit1_init = delta1_init;

    let h1hi_init = rc_h1dash_twice_int_func(hilimit1_init, cbyr);
    let h1dummy1_init = if delta1_init != 0.0 {
        h1hi_init / delta1_init
    } else {
        0.0
    };
    coeffs.h1dash_first_coeff = h1dummy1_init;
    let h1relval = (h1dummy1_init * reltol).abs();

    let (h2hi_init, temp2_init, temp3_init) = if hilimit1_init != 0.0 {
        let temp = rclsqr / (4.0 * hilimit1_init);
        let t2 = if temp >= 100.0 {
            0.0
        } else {
            erfc(temp.sqrt())
        };
        let t3 = (-temp).exp();
        let h2 = (hilimit1_init + rclsqr * 0.5) * t2 - (hilimit1_init * rclsqr / PI).sqrt() * t3;
        (h2, t2, t3)
    } else {
        (0.0, 0.0, 1.0)
    };

    let h2dummy1_init = if delta1_init != 0.0 {
        h2hi_init / delta1_init
    } else {
        0.0
    };
    coeffs.h2_first_coeff = h2dummy1_init;
    let h2relval = (h2dummy1_init * reltol).abs();

    let temp4 = rclsqr.sqrt();
    let temp5 = cbyr.sqrt();
    let h3hi_init = if hilimit1_init != 0.0 {
        let v = 2.0 * (hilimit1_init / PI).sqrt() * temp3_init - temp4 * temp2_init;
        temp5 * v
    } else {
        0.0
    };
    let h3dummy1_init = if delta1_init != 0.0 {
        h3hi_init / delta1_init
    } else {
        0.0
    };
    coeffs.h3dash_first_coeff = h3dummy1_init;
    let h3relval = (h3dummy1_init * reltol).abs();

    let mut h1lo = 0.0_f64;
    let mut h1hi = h1hi_init;
    let mut h1dummy = h1dummy1_init;
    let mut h2lo = 0.0_f64;
    let mut h2hi = h2hi_init;
    let mut h2dummy = h2dummy1_init;
    let mut h3lo = 0.0_f64;
    let mut h3hi = h3hi_init;
    let mut h3dummy = h3dummy1_init;
    let mut prev_hilimit = hilimit1_init;
    let mut doh1 = true;
    let mut doh2 = true;
    let mut doh3 = true;

    for i in (1..=aux_index).rev() {
        let delta = time_points[i] - time_points[i - 1];
        let hilimit = cur_time - time_points[i - 1];

        if doh1 {
            let h1lo_prev = h1lo;
            let _ = h1lo_prev;
            h1lo = h1hi;
            h1hi = rc_h1dash_twice_int_func(hilimit, cbyr);
            let prev = h1dummy;
            h1dummy = if delta != 0.0 {
                (h1hi - h1lo) / delta
            } else {
                0.0
            };
            coeffs.h1dash_coeffs[i] = h1dummy - prev;
            if coeffs.h1dash_coeffs[i].abs() < h1relval {
                doh1 = false;
            }
        } else {
            coeffs.h1dash_coeffs[i] = 0.0;
        }

        let (t2, t3) = if (doh2 || doh3) && hilimit != 0.0 {
            let temp = rclsqr / (4.0 * hilimit);
            let t2 = if temp >= 100.0 {
                0.0
            } else {
                erfc(temp.sqrt())
            };
            let t3 = (-temp).exp();
            (t2, t3)
        } else {
            (0.0, 1.0)
        };

        if doh2 {
            h2lo = h2hi;
            h2hi = if hilimit != 0.0 {
                (hilimit + rclsqr * 0.5) * t2 - (hilimit * rclsqr / PI).sqrt() * t3
            } else {
                0.0
            };
            let prev = h2dummy;
            h2dummy = if delta != 0.0 {
                (h2hi - h2lo) / delta
            } else {
                0.0
            };
            coeffs.h2_coeffs[i] = h2dummy - prev;
            if coeffs.h2_coeffs[i].abs() < h2relval {
                doh2 = false;
            }
        } else {
            coeffs.h2_coeffs[i] = 0.0;
        }

        if doh3 {
            h3lo = h3hi;
            h3hi = if hilimit != 0.0 {
                let v = 2.0 * (hilimit / PI).sqrt() * t3 - temp4 * t2;
                temp5 * v
            } else {
                0.0
            };
            let prev = h3dummy;
            h3dummy = if delta != 0.0 {
                (h3hi - h3lo) / delta
            } else {
                0.0
            };
            coeffs.h3dash_coeffs[i] = h3dummy - prev;
            if coeffs.h3dash_coeffs[i].abs() < h3relval {
                doh3 = false;
            }
        } else {
            coeffs.h3dash_coeffs[i] = 0.0;
        }

        prev_hilimit = hilimit;
    }
    let _ = prev_hilimit;
}

/// Set up convolution coefficients for RLC line.
#[expect(unused_assignments)]
#[expect(clippy::too_many_arguments)]
pub fn rlc_coeffs_setup(
    coeffs: &mut LtraCoeffs,
    td: f64,
    alpha: f64,
    beta: f64,
    cur_time: f64,
    time_points: &[f64],
    time_index: usize,
    reltol: f64,
) {
    coeffs.ensure_size(time_index + 1);

    // Find aux_index: timepoint at or just before cur_time - td
    let aux_index = if td == 0.0 {
        time_index
    } else if cur_time - td <= 0.0 {
        0
    } else {
        let mut exact = false;
        let mut idx = 0;
        for i in (0..=time_index).rev() {
            if cur_time - time_points[i] == td {
                exact = true;
                idx = i;
                break;
            }
            if cur_time - time_points[i] > td {
                idx = i;
                break;
            }
        }
        if exact { idx.saturating_sub(1) } else { idx }
    };
    coeffs.aux_index = aux_index;

    let exp_beta_t = (-beta * td).exp();
    let alpha_sq_t = alpha * alpha * td;

    // h2 and h3dash first coefficients (only if aux_index > 0)
    let mut h2lo = 0.0_f64;
    let mut h2hi = 0.0_f64;
    let mut h2dummy = 0.0_f64;
    let h2relval;
    let mut h3lo = 0.0_f64;
    let mut h3hi = 0.0_f64;
    let mut h3dummy = 0.0_f64;
    let h3relval;

    if aux_index != 0 {
        let lo = td;
        let hi = cur_time - time_points[aux_index];
        let delta = hi - lo;

        h2lo = rlc_h2_func(td, td, alpha, beta);
        let bessel_arg = if hi > td {
            alpha * (hi * hi - td * td).sqrt()
        } else {
            0.0
        };
        let exp_term = (-beta * hi).exp();
        let bi1ox = bess_i1x_over_x(bessel_arg);
        h2hi = if alpha == 0.0 || hi < td {
            0.0
        } else {
            alpha_sq_t * exp_term * bi1ox
        };
        h2dummy = if delta != 0.0 {
            twice_int_lin_func(lo, hi, lo, h2lo, h2hi, lo, hi) / delta
        } else {
            0.0
        };
        coeffs.h2_first_coeff = h2dummy;
        h2relval = (reltol * h2dummy).abs();

        h3lo = 0.0;
        let bi0 = bess_i0(bessel_arg);
        h3hi = if hi <= td || beta == 0.0 {
            0.0
        } else {
            exp_term * bi0 - exp_beta_t
        };
        h3dummy = if delta != 0.0 {
            int_lin_func(lo, hi, h3lo, h3hi, lo, hi) / delta
        } else {
            0.0
        };
        coeffs.h3dash_first_coeff = h3dummy;
        h3relval = (h3dummy * reltol).abs();
    } else {
        coeffs.h2_first_coeff = 0.0;
        coeffs.h3dash_first_coeff = 0.0;
        h2relval = 0.0;
        h3relval = 0.0;
    }

    // h1dash first coefficient (always computed, starts from t=0)
    let hi0 = cur_time - time_points[time_index];
    let delta0 = hi0;
    let exp_term0 = (-beta * hi0).exp();

    let h1hi_init = if beta == 0.0 {
        hi0
    } else if hi0 == 0.0 {
        0.0
    } else {
        (bess_i1(beta * hi0) + bess_i0(beta * hi0)) * hi0 * exp_term0 - hi0
    };
    let h1dummy_init = if delta0 != 0.0 {
        h1hi_init / delta0
    } else {
        0.0
    };
    coeffs.h1dash_first_coeff = h1dummy_init;
    let h1relval = (h1dummy_init * reltol).abs();

    let mut h1lo = 0.0_f64;
    let mut h1hi = h1hi_init;
    let mut h1dummy = h1dummy_init;
    let mut doh1 = true;
    let mut doh2 = true;
    let mut doh3 = true;

    let mut prev_lo = 0.0_f64;
    let mut prev_hi = hi0;
    let mut prev_delta = delta0;
    let mut prev_h2lo = h2lo;
    let mut prev_h2hi = h2hi;

    for i in (1..=time_index).rev() {
        let cur_lo = prev_hi;
        let cur_hi = cur_time - time_points[i - 1];
        let delta = time_points[i] - time_points[i - 1];

        if doh1 || doh2 || doh3 {
            let exp_arg = -beta * cur_hi;
            let exp_term = exp_arg.exp();

            if doh1 {
                h1lo = h1hi;
                h1hi = if beta == 0.0 {
                    cur_hi
                } else if cur_hi == 0.0 {
                    0.0
                } else {
                    (bess_i1(-exp_arg) + bess_i0(-exp_arg)) * cur_hi * exp_term - cur_hi
                };
                let prev = h1dummy;
                h1dummy = if delta != 0.0 {
                    (h1hi - h1lo) / delta
                } else {
                    0.0
                };
                coeffs.h1dash_coeffs[i] = h1dummy - prev;
                if coeffs.h1dash_coeffs[i].abs() <= h1relval {
                    doh1 = false;
                }
            } else {
                coeffs.h1dash_coeffs[i] = 0.0;
            }

            if i <= aux_index {
                let bessel_arg = if cur_hi > td {
                    alpha * (cur_hi * cur_hi - td * td).sqrt()
                } else {
                    0.0
                };

                if doh2 {
                    prev_h2lo = h2lo;
                    prev_h2hi = h2hi;
                    let prev = h2dummy;
                    h2lo = h2hi;
                    let bi1ox = bess_i1x_over_x(bessel_arg);
                    h2hi = if alpha == 0.0 || cur_hi < td {
                        0.0
                    } else {
                        alpha_sq_t * exp_term * bi1ox
                    };
                    h2dummy = if delta != 0.0 {
                        twice_int_lin_func(cur_lo, cur_hi, cur_lo, h2lo, h2hi, cur_lo, cur_hi)
                            / delta
                    } else {
                        0.0
                    };
                    coeffs.h2_coeffs[i] = h2dummy - prev
                        + int_lin_func(prev_lo, prev_hi, prev_h2lo, prev_h2hi, prev_lo, prev_hi);
                    if coeffs.h2_coeffs[i].abs() <= h2relval {
                        doh2 = false;
                    }
                } else {
                    coeffs.h2_coeffs[i] = 0.0;
                }

                if doh3 {
                    h3lo = h3hi;
                    let bi0 = bess_i0(bessel_arg);
                    h3hi = if cur_hi <= td || beta == 0.0 {
                        0.0
                    } else {
                        exp_term * bi0 - exp_beta_t
                    };
                    let prev = h3dummy;
                    h3dummy = if delta != 0.0 {
                        int_lin_func(cur_lo, cur_hi, h3lo, h3hi, cur_lo, cur_hi) / delta
                    } else {
                        0.0
                    };
                    coeffs.h3dash_coeffs[i] = h3dummy - prev;
                    if coeffs.h3dash_coeffs[i].abs() <= h3relval {
                        doh3 = false;
                    }
                } else {
                    coeffs.h3dash_coeffs[i] = 0.0;
                }
            }
        } else {
            coeffs.h1dash_coeffs[i] = 0.0;
            if i <= aux_index {
                coeffs.h2_coeffs[i] = 0.0;
                coeffs.h3dash_coeffs[i] = 0.0;
            }
        }

        prev_lo = cur_lo;
        prev_hi = cur_hi;
        prev_delta = delta;
    }
    let _ = prev_delta;
}

// ---------------------------------------------------------------------------
// Linear interpolation helpers
// ---------------------------------------------------------------------------

/// Linear interpolation coefficients for time t between t1 and t2.
pub fn lin_interp(t: f64, t1: f64, t2: f64) -> (f64, f64) {
    if t1 == t2 {
        return (0.5, 0.5);
    }
    let f = (t - t1) / (t2 - t1);
    (1.0 - f, f)
}

// ---------------------------------------------------------------------------
// DC stamp
// ---------------------------------------------------------------------------

/// Stamp LTRA in DC mode (simple resistor-like behavior).
pub fn stamp_ltra_dc(
    inst: &LtraInstance,
    system: &mut crate::sparse::LinearSystem,
    num_nodes: usize,
) {
    let br1 = num_nodes + inst.br_eq1;
    let br2 = num_nodes + inst.br_eq2;

    match inst.model.special_case {
        LtraCase::Rg => {
            let model = &inst.model;
            // Compute RG DC parameters
            let arg = model.length * (model.r * model.g).sqrt();
            let exp_pos = arg.exp();
            let exp_neg = (-arg).exp();
            let cosh_val = 0.5 * (exp_pos + exp_neg);

            let r_rs = if model.g <= 1.0e-10 {
                model.length * model.r
            } else {
                0.5 * (exp_pos - exp_neg) * (model.r / model.g).sqrt()
            };

            let r_gs = if model.r <= 1.0e-10 {
                model.length * model.g
            } else {
                0.5 * (exp_pos - exp_neg) * (model.g / model.r).sqrt()
            };

            // Branch eq 1: V(pos1) - V(neg1) - cosh*V(pos2) + cosh*V(neg2) + rRs*I2 = 0
            add_mna(system, br1, inst.pos1_idx, 1.0);
            add_mna(system, br1, inst.neg1_idx, -1.0);
            add_mna(system, br1, inst.pos2_idx, -cosh_val);
            add_mna(system, br1, inst.neg2_idx, cosh_val);
            system.matrix.add(br1, br2, r_rs);

            // Branch eq 2: cosh*I2 - rGs*(V(pos2)-V(neg2)) + I1 = 0
            system.matrix.add(br2, br2, cosh_val);
            add_mna(system, br2, inst.pos2_idx, -r_gs);
            add_mna(system, br2, inst.neg2_idx, r_gs);
            system.matrix.add(br2, br1, 1.0);

            // Topology: I1 enters pos1, exits neg1; I2 enters pos2, exits neg2
            add_mna_opt(system, inst.pos1_idx, br1, 1.0);
            add_mna_opt(system, inst.neg1_idx, br1, -1.0);
            add_mna_opt(system, inst.pos2_idx, br2, 1.0);
            add_mna_opt(system, inst.neg2_idx, br2, -1.0);
        }
        LtraCase::Rc | LtraCase::Lc | LtraCase::Rlc => {
            // Simple resistor-like behavior for DC
            // Branch eq 1: I1 + I2 = 0  (series connection)
            // Branch eq 2: V(pos1) - V(pos2) - R*l*I1 = 0
            add_mna_opt(system, inst.pos1_idx, br1, 1.0);
            add_mna_opt(system, inst.neg1_idx, br1, -1.0);
            add_mna_opt(system, inst.pos2_idx, br2, 1.0);
            add_mna_opt(system, inst.neg2_idx, br2, -1.0);

            system.matrix.add(br1, br1, 1.0);
            system.matrix.add(br1, br2, 1.0);
            add_mna(system, br2, inst.pos1_idx, 1.0);
            add_mna(system, br2, inst.pos2_idx, -1.0);
            system
                .matrix
                .add(br2, br1, -inst.model.r * inst.model.length);
        }
    }
}

/// Helper to add entry to sparse matrix, handling Option<usize> ground nodes.
fn add_mna(system: &mut crate::sparse::LinearSystem, row: usize, col: Option<usize>, val: f64) {
    if let Some(c) = col {
        system.matrix.add(row, c, val);
    }
}

fn add_mna_opt(system: &mut crate::sparse::LinearSystem, row: Option<usize>, col: usize, val: f64) {
    if let Some(r) = row {
        system.matrix.add(r, col, val);
    }
}

// ---------------------------------------------------------------------------
// Transient stamp
// ---------------------------------------------------------------------------

/// Result of computing LTRA transient excitation.
pub struct LtraExcitation {
    pub input1: f64,
    pub input2: f64,
}

/// Stamp the LTRA matrix for transient analysis and compute RHS excitation.
///
/// Returns the matrix coefficients added. The RHS excitation is added to rhs[br1] and rhs[br2].
#[expect(clippy::too_many_arguments)]
pub fn stamp_ltra_transient(
    inst: &LtraInstance,
    state: &LtraState,
    coeffs: &LtraCoeffs,
    system: &mut crate::sparse::LinearSystem,
    num_nodes: usize,
    cur_time: f64,
    time_points: &[f64],
    time_index: usize,
) {
    let br1 = num_nodes + inst.br_eq1;
    let br2 = num_nodes + inst.br_eq2;
    let model = &inst.model;

    match model.special_case {
        LtraCase::Rlc => {
            // h1dash first coeff contribution
            let d1 = model.admit * coeffs.h1dash_first_coeff;
            add_mna(system, br1, inst.pos1_idx, d1);
            add_mna(system, br1, inst.neg1_idx, -d1);
            add_mna(system, br2, inst.pos2_idx, d1);
            add_mna(system, br2, inst.neg2_idx, -d1);

            // Lossless-like parts (fall through from RLC to LC)
            add_mna(system, br1, inst.pos1_idx, model.admit);
            add_mna(system, br1, inst.neg1_idx, -model.admit);
            system.matrix.add(br1, br1, -1.0);
            add_mna_opt(system, inst.pos1_idx, br1, 1.0);
            add_mna_opt(system, inst.neg1_idx, br1, -1.0);

            add_mna(system, br2, inst.pos2_idx, model.admit);
            add_mna(system, br2, inst.neg2_idx, -model.admit);
            system.matrix.add(br2, br2, -1.0);
            add_mna_opt(system, inst.pos2_idx, br2, 1.0);
            add_mna_opt(system, inst.neg2_idx, br2, -1.0);

            // Compute RHS excitation
            let excitation =
                compute_rlc_excitation(inst, state, coeffs, time_points, time_index, cur_time);
            system.rhs[br1] += excitation.input1;
            system.rhs[br2] += excitation.input2;
        }
        LtraCase::Lc => {
            // Lossless line
            add_mna(system, br1, inst.pos1_idx, model.admit);
            add_mna(system, br1, inst.neg1_idx, -model.admit);
            system.matrix.add(br1, br1, -1.0);
            add_mna_opt(system, inst.pos1_idx, br1, 1.0);
            add_mna_opt(system, inst.neg1_idx, br1, -1.0);

            add_mna(system, br2, inst.pos2_idx, model.admit);
            add_mna(system, br2, inst.neg2_idx, -model.admit);
            system.matrix.add(br2, br2, -1.0);
            add_mna_opt(system, inst.pos2_idx, br2, 1.0);
            add_mna_opt(system, inst.neg2_idx, br2, -1.0);

            // Compute LC excitation (lossless-like parts only)
            let excitation = compute_lc_excitation(inst, state, time_points, time_index, cur_time);
            system.rhs[br1] += excitation.input1;
            system.rhs[br2] += excitation.input2;
        }
        LtraCase::Rc => {
            // Non-convolution parts
            system.matrix.add(br1, br1, -1.0);
            add_mna_opt(system, inst.pos1_idx, br1, 1.0);
            add_mna_opt(system, inst.neg1_idx, br1, -1.0);

            system.matrix.add(br2, br2, -1.0);
            add_mna_opt(system, inst.pos2_idx, br2, 1.0);
            add_mna_opt(system, inst.neg2_idx, br2, -1.0);

            // Convolution first terms
            let d1 = coeffs.h1dash_first_coeff;
            add_mna(system, br1, inst.pos1_idx, d1);
            add_mna(system, br1, inst.neg1_idx, -d1);
            add_mna(system, br2, inst.pos2_idx, d1);
            add_mna(system, br2, inst.neg2_idx, -d1);

            let d2 = coeffs.h2_first_coeff;
            system.matrix.add(br1, br2, -d2);
            system.matrix.add(br2, br1, -d2);

            let d3 = coeffs.h3dash_first_coeff;
            add_mna(system, br1, inst.pos2_idx, -d3);
            add_mna(system, br1, inst.neg2_idx, d3);
            add_mna(system, br2, inst.pos1_idx, -d3);
            add_mna(system, br2, inst.neg1_idx, d3);

            // Compute RC excitation
            let excitation = compute_rc_excitation(inst, state, coeffs, time_index);
            system.rhs[br1] += excitation.input1;
            system.rhs[br2] += excitation.input2;
        }
        LtraCase::Rg => {
            // RG case uses DC stamp in transient
            stamp_ltra_dc(inst, system, num_nodes);
        }
    }
}

// ---------------------------------------------------------------------------
// Excitation computation
// ---------------------------------------------------------------------------

fn compute_lc_excitation(
    inst: &LtraInstance,
    state: &LtraState,
    time_points: &[f64],
    time_index: usize,
    cur_time: f64,
) -> LtraExcitation {
    let model = &inst.model;
    let mut input1 = 0.0;
    let mut input2 = 0.0;

    let td_over = cur_time > model.td;

    if td_over {
        // Find delayed timepoint
        let (v1d, i1d, v2d, i2d) =
            interpolate_delayed(state, time_points, time_index, cur_time, model.td);
        input1 += model.attenuation * (v2d * model.admit + i2d);
        input2 += model.attenuation * (v1d * model.admit + i1d);
    } else {
        input1 += model.attenuation * (state.init_v2 * model.admit + state.init_i2);
        input2 += model.attenuation * (state.init_v1 * model.admit + state.init_i1);
    }

    LtraExcitation { input1, input2 }
}

fn compute_rlc_excitation(
    inst: &LtraInstance,
    state: &LtraState,
    coeffs: &LtraCoeffs,
    time_points: &[f64],
    time_index: usize,
    cur_time: f64,
) -> LtraExcitation {
    let model = &inst.model;
    let mut input1 = 0.0;
    let mut input2 = 0.0;

    let td_over = cur_time > model.td;

    // Convolution of h1dash with v1 and v2
    let mut d1 = 0.0;
    let mut d2 = 0.0;
    for i in (1..=time_index).rev() {
        if coeffs.h1dash_coeffs[i] != 0.0 {
            d1 += coeffs.h1dash_coeffs[i] * (state.v1[i] - state.init_v1);
            d2 += coeffs.h1dash_coeffs[i] * (state.v2[i] - state.init_v2);
        }
    }
    d1 += state.init_v1 * model.int_h1dash;
    d2 += state.init_v2 * model.int_h1dash;
    d1 -= state.init_v1 * coeffs.h1dash_first_coeff;
    d2 -= state.init_v2 * coeffs.h1dash_first_coeff;
    input1 -= d1 * model.admit;
    input2 -= d2 * model.admit;

    // Convolution of h2 with i2 and i1
    d1 = 0.0;
    d2 = 0.0;
    if td_over {
        let (v1d, i1d, v2d, i2d) =
            interpolate_delayed(state, time_points, time_index, cur_time, model.td);

        d1 = (i2d - state.init_i2) * coeffs.h2_first_coeff;
        d2 = (i1d - state.init_i1) * coeffs.h2_first_coeff;

        for i in (1..=coeffs.aux_index).rev() {
            if coeffs.h2_coeffs[i] != 0.0 {
                d1 += coeffs.h2_coeffs[i] * (state.i2[i] - state.init_i2);
                d2 += coeffs.h2_coeffs[i] * (state.i1[i] - state.init_i1);
            }
        }

        // Initial condition terms
        d1 += state.init_i2 * model.int_h2;
        d2 += state.init_i1 * model.int_h2;
        input1 += d1;
        input2 += d2;

        // Convolution of h3dash with v2 and v1
        d1 = (v2d - state.init_v2) * coeffs.h3dash_first_coeff;
        d2 = (v1d - state.init_v1) * coeffs.h3dash_first_coeff;

        for i in (1..=coeffs.aux_index).rev() {
            if coeffs.h3dash_coeffs[i] != 0.0 {
                d1 += coeffs.h3dash_coeffs[i] * (state.v2[i] - state.init_v2);
                d2 += coeffs.h3dash_coeffs[i] * (state.v1[i] - state.init_v1);
            }
        }

        d1 += state.init_v2 * model.int_h3dash;
        d2 += state.init_v1 * model.int_h3dash;
        input1 += model.admit * d1;
        input2 += model.admit * d2;

        // Lossless-like parts
        input1 += model.attenuation * (v2d * model.admit + i2d);
        input2 += model.attenuation * (v1d * model.admit + i1d);
    } else {
        d1 += state.init_i2 * model.int_h2;
        d2 += state.init_i1 * model.int_h2;
        input1 += d1;
        input2 += d2;

        d1 = state.init_v2 * model.int_h3dash;
        d2 = state.init_v1 * model.int_h3dash;
        input1 += model.admit * d1;
        input2 += model.admit * d2;

        // Lossless-like parts (before td)
        input1 += model.attenuation * (state.init_v2 * model.admit + state.init_i2);
        input2 += model.attenuation * (state.init_v1 * model.admit + state.init_i1);
    }

    LtraExcitation { input1, input2 }
}

fn compute_rc_excitation(
    inst: &LtraInstance,
    state: &LtraState,
    coeffs: &LtraCoeffs,
    time_index: usize,
) -> LtraExcitation {
    let model = &inst.model;
    let mut input1 = 0.0;
    let mut input2 = 0.0;

    // Convolution of h1dash with v1 and v2
    let mut d1 = 0.0;
    let mut d2 = 0.0;
    for i in (1..=time_index).rev() {
        if coeffs.h1dash_coeffs[i] != 0.0 {
            d1 += coeffs.h1dash_coeffs[i] * (state.v1[i] - state.init_v1);
            d2 += coeffs.h1dash_coeffs[i] * (state.v2[i] - state.init_v2);
        }
    }
    d1 += state.init_v1 * model.int_h1dash;
    d2 += state.init_v2 * model.int_h1dash;
    d1 -= state.init_v1 * coeffs.h1dash_first_coeff;
    d2 -= state.init_v2 * coeffs.h1dash_first_coeff;
    input1 -= d1;
    input2 -= d2;

    // Convolution of h2 with i2 and i1
    d1 = 0.0;
    d2 = 0.0;
    for i in (1..=time_index).rev() {
        if coeffs.h2_coeffs[i] != 0.0 {
            d1 += coeffs.h2_coeffs[i] * (state.i2[i] - state.init_i2);
            d2 += coeffs.h2_coeffs[i] * (state.i1[i] - state.init_i1);
        }
    }
    d1 += state.init_i2 * model.int_h2;
    d2 += state.init_i1 * model.int_h2;
    d1 -= state.init_i2 * coeffs.h2_first_coeff;
    d2 -= state.init_i1 * coeffs.h2_first_coeff;
    input1 += d1;
    input2 += d2;

    // Convolution of h3dash with v2 and v1
    d1 = 0.0;
    d2 = 0.0;
    for i in (1..=time_index).rev() {
        if coeffs.h3dash_coeffs[i] != 0.0 {
            d1 += coeffs.h3dash_coeffs[i] * (state.v2[i] - state.init_v2);
            d2 += coeffs.h3dash_coeffs[i] * (state.v1[i] - state.init_v1);
        }
    }
    d1 += state.init_v2 * model.int_h3dash;
    d2 += state.init_v1 * model.int_h3dash;
    d1 -= state.init_v2 * coeffs.h3dash_first_coeff;
    d2 -= state.init_v1 * coeffs.h3dash_first_coeff;
    input1 += d1;
    input2 += d2;

    LtraExcitation { input1, input2 }
}

/// Interpolate delayed signal values (at time cur_time - td).
fn interpolate_delayed(
    state: &LtraState,
    time_points: &[f64],
    time_index: usize,
    cur_time: f64,
    td: f64,
) -> (f64, f64, f64, f64) {
    let delayed_time = cur_time - td;

    // Find index i such that time_points[i] < delayed_time <= time_points[i+1]
    let mut isaved = 0;
    for i in (0..=time_index).rev() {
        if time_points[i] < delayed_time {
            isaved = i;
            break;
        }
    }

    if isaved == time_index {
        isaved = isaved.saturating_sub(1);
    }

    // Linear interpolation
    let t2 = time_points[isaved];
    let t3 = time_points[isaved + 1];
    let (lf2, lf3) = lin_interp(delayed_time, t2, t3);

    let v1d = state.v1[isaved] * lf2 + state.v1[isaved + 1] * lf3;
    let i1d = state.i1[isaved] * lf2 + state.i1[isaved + 1] * lf3;
    let v2d = state.v2[isaved] * lf2 + state.v2[isaved + 1] * lf3;
    let i2d = state.i2[isaved] * lf2 + state.i2[isaved + 1] * lf3;

    (v1d, i1d, v2d, i2d)
}

// ---------------------------------------------------------------------------
// erfc implementation (complementary error function)
// ---------------------------------------------------------------------------

/// Complementary error function: erfc(x) = 1 - erf(x).
/// Uses rational approximation for good accuracy.
fn erfc(x: f64) -> f64 {
    // Use libm's erfc if available, or a polynomial approximation
    libm::erfc(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bessel_i0_small() {
        let val = bess_i0(0.0);
        assert!((val - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_bessel_i0_large() {
        let val = bess_i0(5.0);
        // I0(5) ≈ 27.2399
        assert!((val - 27.2399).abs() < 0.01);
    }

    #[test]
    fn test_bessel_i1_zero() {
        let val = bess_i1(0.0);
        assert!(val.abs() < 1e-10);
    }

    #[test]
    fn test_bessel_i1x_over_x_zero() {
        let val = bess_i1x_over_x(0.0);
        assert!((val - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_ltra_classify_lc() {
        assert_eq!(classify_line(0.0, 1e-9, 0.0, 1e-12), LtraCase::Lc);
    }

    #[test]
    fn test_ltra_classify_rlc() {
        assert_eq!(
            classify_line(12.45, 8.972e-9, 0.0, 0.468e-12),
            LtraCase::Rlc
        );
    }

    #[test]
    fn test_ltra_classify_rc() {
        assert_eq!(classify_line(10.0, 0.0, 0.0, 1e-12), LtraCase::Rc);
    }

    #[test]
    fn test_ltra_classify_rg() {
        assert_eq!(classify_line(100.0, 0.0, 0.01, 0.0), LtraCase::Rg);
    }

    #[test]
    fn test_ltra_model_rlc_precompute() {
        let mut model = LtraModel {
            r: 12.45,
            l: 8.972e-9,
            g: 0.0,
            c: 0.468e-12,
            length: 16.0,
            td: 0.0,
            imped: 0.0,
            admit: 0.0,
            alpha: 0.0,
            beta: 0.0,
            attenuation: 1.0,
            c_by_r: 0.0,
            rclsqr: 0.0,
            int_h1dash: 0.0,
            int_h2: 0.0,
            int_h3dash: 0.0,
            coshl_root_gr: 0.0,
            r_rs_l_rgr_or_g: 0.0,
            r_gs_l_rgr_or_r: 0.0,
            reltol: 1.0,
            abstol: 1.0,
            st_line_reltol: 0.0,
            st_line_abstol: 0.0,
            chop_reltol: 0.0,
            chop_abstol: 0.0,
            step_limit: false,
            special_case: LtraCase::Rlc,
        };
        model.precompute();

        // Z = sqrt(L/C) = sqrt(8.972e-9 / 0.468e-12) ≈ 138.4
        assert!((model.imped - 138.4).abs() < 1.0);
        // td = sqrt(LC)*len
        assert!(model.td > 0.0);
        assert!(model.attenuation > 0.0 && model.attenuation < 1.0);
    }

    #[test]
    fn test_lin_interp() {
        let (c1, c2) = lin_interp(0.5, 0.0, 1.0);
        assert!((c1 - 0.5).abs() < 1e-10);
        assert!((c2 - 0.5).abs() < 1e-10);
    }
}
