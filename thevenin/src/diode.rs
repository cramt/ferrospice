//! Diode device model — Shockley equation with Newton-Raphson companion model.
//!
//! Implements the standard SPICE diode: I = IS * (exp(V / (N * Vt)) - 1)
//! with optional series resistance RS and breakdown voltage BV.

use thevenin_types::{Expr, ModelDef, Param};

use crate::physics::{KBOQ, safe_exp};

/// Thermal voltage at nominal temperature (27°C = 300.15K).
/// VT = k/q * T = 8.617087e-5 * 300.15 = 0.025864187 V
/// Matches ngspice: VT_NOM = KOverQ * (273.15 + 27.0).
pub const VT_NOM: f64 = KBOQ * 300.15;

/// Diode model parameters matching ngspice defaults.
#[derive(Debug, Clone)]
pub struct DiodeModel {
    /// Saturation current (default 1e-14 A).
    pub is: f64,
    /// Emission coefficient (default 1.0).
    pub n: f64,
    /// Series resistance (default 0.0 Ω).
    pub rs: f64,
    /// Breakdown voltage (default: None = no breakdown).
    pub bv: Option<f64>,
    /// Current at breakdown voltage (default 1e-3 A).
    pub ibv: f64,
    /// Zero-bias junction capacitance (default 0.0 F) — for future AC/transient.
    pub cjo: f64,
    /// Junction potential (default 1.0 V).
    pub vj: f64,
    /// Grading coefficient (default 0.5).
    pub m: f64,
    /// Transit time (default 0.0 s).
    pub tt: f64,
    /// Flicker noise coefficient (default 0).
    pub kf: f64,
    /// Flicker noise exponent (default 1).
    pub af: f64,
}

impl Default for DiodeModel {
    fn default() -> Self {
        Self {
            is: 1e-14,
            n: 1.0,
            rs: 0.0,
            bv: None,
            ibv: 1e-3,
            cjo: 0.0,
            vj: 1.0,
            m: 0.5,
            tt: 0.0,
            kf: 0.0,
            af: 1.0,
        }
    }
}

impl DiodeModel {
    /// Create a `DiodeModel` from a netlist `.model` definition.
    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let mut m = Self::default();
        for p in &model_def.params {
            if let Expr::Num(v) = &p.value {
                match p.name.to_uppercase().as_str() {
                    "IS" => m.is = *v,
                    "N" => m.n = *v,
                    "RS" => m.rs = *v,
                    "BV" => m.bv = Some(*v),
                    "IBV" => m.ibv = *v,
                    "CJO" | "CJ0" => m.cjo = *v,
                    "VJ" => m.vj = *v,
                    "M" => m.m = *v,
                    "TT" => m.tt = *v,
                    "KF" => m.kf = *v,
                    "AF" => m.af = *v,
                    _ => {} // ignore unknown params
                }
            }
        }
        m
    }

    /// Create a `DiodeModel` from instance parameters (on the D element line).
    pub fn with_instance_params(mut self, params: &[Param]) -> Self {
        for p in params {
            if let Expr::Num(v) = &p.value {
                match p.name.to_uppercase().as_str() {
                    "IS" => self.is = *v,
                    "N" => self.n = *v,
                    "RS" => self.rs = *v,
                    _ => {}
                }
            }
        }
        self
    }

    /// Effective thermal voltage: N * Vt.
    fn vte(&self) -> f64 {
        self.n * VT_NOM
    }

    /// Diode junction current: I = IS * (exp(V / Vte) - 1).
    pub fn current(&self, v: f64) -> f64 {
        let vte = self.vte();
        if v >= -3.0 * vte {
            // Forward / moderate reverse bias.
            self.is * (safe_exp(v / vte) - 1.0)
        } else if let Some(bv) = self.bv {
            // Reverse breakdown region.
            let vbr = -bv;
            if v <= vbr {
                -self.ibv * safe_exp(-(v + bv) / vte)
            } else {
                -self.is
            }
        } else {
            // Deep reverse, no breakdown — just -IS.
            -self.is
        }
    }

    /// Diode junction conductance: dI/dV.
    pub fn conductance(&self, v: f64) -> f64 {
        let vte = self.vte();
        if v >= -3.0 * vte {
            self.is / vte * safe_exp(v / vte)
        } else if let Some(bv) = self.bv {
            let vbr = -bv;
            if v <= vbr {
                self.ibv / vte * safe_exp(-(v + bv) / vte)
            } else {
                // Deep reverse, negligible conductance.
                self.is / vte
            }
        } else {
            self.is / vte
        }
    }

    /// Compute the NR companion model at a given junction voltage.
    ///
    /// Returns `(g_eq, i_eq)` where:
    /// - `g_eq` is the equivalent conductance to stamp into the matrix
    /// - `i_eq` is the equivalent current source to stamp into the RHS
    ///
    /// The linearization is: I ≈ g_eq * V + i_eq
    /// So: i_eq = I(V_op) - g_eq * V_op
    pub fn companion(&self, v: f64) -> (f64, f64) {
        let g = self.conductance(v);
        let i = self.current(v);
        let i_eq = i - g * v;
        (g, i_eq)
    }

    /// Junction (depletion) capacitance at a given voltage.
    ///
    /// Cj(V) = CJO / (1 - V/VJ)^M for V < FC*VJ (forward coefficient FC=0.5)
    /// Uses linear extrapolation beyond FC*VJ to avoid singularity.
    pub fn junction_capacitance(&self, v: f64) -> f64 {
        if self.cjo <= 0.0 {
            return 0.0;
        }
        let fc = 0.5; // forward coefficient (ngspice default)
        let fc_vj = fc * self.vj;
        if v < fc_vj {
            self.cjo / (1.0 - v / self.vj).powf(self.m)
        } else {
            // Linear extrapolation beyond FC*VJ
            let base = self.cjo / (1.0 - fc).powf(self.m);
            let slope = base * self.m / (self.vj * (1.0 - fc));
            base + slope * (v - fc_vj)
        }
    }

    /// Whether this diode has series resistance.
    pub fn has_series_resistance(&self) -> bool {
        self.rs > 0.0
    }
}

/// Limit voltage change between NR iterations to aid convergence.
///
/// Matches ngspice `DEVpnjlim` — limits the voltage step to prevent
/// the exponential from exploding.
pub fn pnjlim(v_new: f64, v_old: f64, vt: f64, vcrit: f64) -> f64 {
    if v_new > vcrit && (v_new - v_old).abs() > 2.0 * vt {
        if v_old > 0.0 {
            let arg = (v_new - v_old) / vt;
            if arg > 0.0 {
                v_old + vt * (2.0 + arg.ln())
            } else {
                v_old - vt * (2.0 + (-arg).ln())
            }
        } else {
            vt * (v_new / vt).ln()
        }
    } else {
        v_new
    }
}

/// Critical voltage for junction voltage limiting.
/// Vcrit = Vt * ln(Vt / (sqrt(2) * IS))
pub fn vcrit(vt: f64, is: f64) -> f64 {
    vt * (vt / (std::f64::consts::SQRT_2 * is)).ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    #[test]
    fn test_default_diode_model() {
        let d = DiodeModel::default();
        assert_eq!(d.is, 1e-14);
        assert_eq!(d.n, 1.0);
        assert_eq!(d.rs, 0.0);
        assert!(d.bv.is_none());
    }

    #[test]
    fn test_diode_forward_current() {
        let d = DiodeModel::default();
        // At V=0, I should be ~0
        assert_abs_diff_eq!(d.current(0.0), 0.0, epsilon = 1e-14);
        // At V=0.7V, should have significant forward current
        let i = d.current(0.7);
        assert!(i > 1e-3, "forward current at 0.7V should be > 1mA, got {i}");
    }

    #[test]
    fn test_diode_reverse_current() {
        let d = DiodeModel::default();
        // Deep reverse: I ≈ -IS
        let i = d.current(-1.0);
        assert!(i < 0.0);
        assert_abs_diff_eq!(i, -1e-14, epsilon = 1e-15);
    }

    #[test]
    fn test_companion_model() {
        let d = DiodeModel::default();
        let v = 0.6;
        let (g, i_eq) = d.companion(v);
        // g should be positive
        assert!(g > 0.0);
        // Verify linearization: I(v) = g*v + i_eq
        let i_check = g * v + i_eq;
        assert_abs_diff_eq!(i_check, d.current(v), epsilon = 1e-15);
    }

    #[test]
    fn test_from_model_def() {
        let model = ModelDef {
            name: "D1N4148".to_string(),
            kind: "D".to_string(),
            params: vec![
                Param {
                    name: "IS".to_string(),
                    value: Expr::Num(2.52e-9),
                },
                Param {
                    name: "RS".to_string(),
                    value: Expr::Num(0.568),
                },
                Param {
                    name: "N".to_string(),
                    value: Expr::Num(1.752),
                },
                Param {
                    name: "BV".to_string(),
                    value: Expr::Num(100.0),
                },
            ],
        };
        let d = DiodeModel::from_model_def(&model);
        assert_abs_diff_eq!(d.is, 2.52e-9, epsilon = 1e-15);
        assert_abs_diff_eq!(d.rs, 0.568, epsilon = 1e-15);
        assert_abs_diff_eq!(d.n, 1.752, epsilon = 1e-15);
        assert_eq!(d.bv, Some(100.0));
    }

    #[test]
    fn test_pnjlim() {
        let vt = VT_NOM;
        let is = 1e-14;
        let vc = vcrit(vt, is);

        // Small change — no limiting
        let v = pnjlim(0.65, 0.6, vt, vc);
        assert_abs_diff_eq!(v, 0.65, epsilon = 1e-15);

        // Large jump from 0.6 to 10.0 — should be limited
        let v = pnjlim(10.0, 0.6, vt, vc);
        assert!(v < 10.0, "should be limited, got {v}");
        assert!(v > 0.6, "should increase, got {v}");
    }
}
