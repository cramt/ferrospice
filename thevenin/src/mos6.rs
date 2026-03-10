//! MOSFET Level 6 (Sakurai-Newton n-th power) device model.
//!
//! Implements the MOS level 6 model from ngspice, which uses an exponential
//! (n-th power) I-V characteristic rather than the quadratic Shichman-Hodges
//! model of Level 1. This provides a better fit to measured device
//! characteristics, particularly in the subthreshold region.
//!
//! Reference: ngspice mos6load.c, mos6temp.c, mos6mpar.c

use thevenin_types::{Expr, ModelDef};

use crate::diode::VT_NOM;
use crate::mosfet::{MosfetCompanion, MosfetType};
use crate::physics::{EXP_LIMIT, safe_exp};

/// MOS Level 6 model parameters.
#[derive(Debug, Clone)]
pub struct Mos6Model {
    /// NMOS or PMOS.
    pub mos_type: MosfetType,
    /// Threshold voltage (V). Default 0.
    pub vto: f64,
    /// Saturation current coefficient (A/V^NC). Default 0.
    pub kc: f64,
    /// Saturation current exponent. Default 2.
    pub nc: f64,
    /// Saturation voltage coefficient (V^(1-NV)). Default 1.
    pub kv: f64,
    /// Saturation voltage exponent. Default 1.
    pub nv: f64,
    /// Threshold voltage coefficient. Default 0.
    pub nvth: f64,
    /// PS model parameter. Default 0.
    pub ps: f64,
    /// Body effect coefficient. Default 0.
    pub gamma: f64,
    /// Body effect linear coefficient. Default 0.
    pub gamma1: f64,
    /// DIBL coefficient (drain-induced barrier lowering). Default 0.
    pub sigma: f64,
    /// Surface potential (V). Default 0.6.
    pub phi: f64,
    /// Channel length modulation (legacy, used if lambda0 not given). Default 0.
    pub lambda: f64,
    /// Channel length modulation coefficient. Default 0.
    pub lambda0: f64,
    /// Body-bias dependent CLM coefficient. Default 0.
    pub lambda1: f64,
    /// Drain resistance (Ω). Default 0.
    pub rd: f64,
    /// Source resistance (Ω). Default 0.
    pub rs: f64,
    /// Bulk-drain zero-bias junction capacitance (F). Default 0.
    pub cbd: f64,
    /// Bulk-source zero-bias junction capacitance (F). Default 0.
    pub cbs: f64,
    /// Bulk junction saturation current (A). Default 1e-14.
    pub is: f64,
    /// Bulk junction potential (V). Default 0.8.
    pub pb: f64,
    /// Gate-source overlap capacitance per unit width (F/m). Default 0.
    pub cgso: f64,
    /// Gate-drain overlap capacitance per unit width (F/m). Default 0.
    pub cgdo: f64,
    /// Gate-bulk overlap capacitance per unit length (F/m). Default 0.
    pub cgbo: f64,
    /// Bottom junction capacitance per unit area (F/m²). Default 0.
    pub cj: f64,
    /// Bottom junction grading coefficient. Default 0.5.
    pub mj: f64,
    /// Sidewall junction capacitance per unit length (F/m). Default 0.
    pub cjsw: f64,
    /// Sidewall junction grading coefficient. Default 0.5.
    pub mjsw: f64,
    /// Oxide thickness (m). Default 1e-7.
    pub tox: f64,
    /// Lateral diffusion (m). Default 0.
    pub ld: f64,
    /// Substrate doping (1/cm³). Default 0.
    pub nsub: f64,
    /// Surface mobility (cm²/V·s). Default 600.
    pub u0: f64,
    /// Forward bias depletion cap coefficient. Default 0.5.
    pub fc: f64,
    /// Flicker noise coefficient. Default 0.
    pub kf: f64,
    /// Flicker noise exponent. Default 1.
    pub af: f64,
    /// Whether lambda0 was explicitly given.
    lambda0_given: bool,
}

impl Mos6Model {
    /// Create a new `Mos6Model` with default parameters.
    pub fn new(mos_type: MosfetType) -> Self {
        Self {
            mos_type,
            vto: 0.0,
            kc: 0.0,
            nc: 2.0,
            kv: 1.0,
            nv: 1.0,
            nvth: 0.0,
            ps: 0.0,
            gamma: 0.0,
            gamma1: 0.0,
            sigma: 0.0,
            phi: 0.6,
            lambda: 0.0,
            lambda0: 0.0,
            lambda1: 0.0,
            rd: 0.0,
            rs: 0.0,
            cbd: 0.0,
            cbs: 0.0,
            is: 1e-14,
            pb: 0.8,
            cgso: 0.0,
            cgdo: 0.0,
            cgbo: 0.0,
            cj: 0.0,
            mj: 0.5,
            cjsw: 0.0,
            mjsw: 0.5,
            tox: 1e-7,
            ld: 0.0,
            nsub: 0.0,
            u0: 600.0,
            fc: 0.5,
            kf: 0.0,
            af: 1.0,
            lambda0_given: false,
        }
    }

    /// Create a `Mos6Model` from a netlist `.model` definition.
    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let mos_type = if model_def.kind.to_uppercase().contains("PMOS") {
            MosfetType::Pmos
        } else {
            MosfetType::Nmos
        };
        let mut m = Self::new(mos_type);
        for p in &model_def.params {
            if let Expr::Num(v) = &p.value {
                match p.name.to_uppercase().as_str() {
                    "VTO" | "VT0" => m.vto = *v,
                    "KC" => m.kc = *v,
                    "NC" => m.nc = *v,
                    "KV" => m.kv = *v,
                    "NV" => m.nv = *v,
                    "NVTH" => m.nvth = *v,
                    "PS" => m.ps = *v,
                    "GAMMA" => m.gamma = *v,
                    "GAMMA1" => m.gamma1 = *v,
                    "SIGMA" => m.sigma = *v,
                    "PHI" => m.phi = *v,
                    "LAMBDA" => m.lambda = *v,
                    "LAMBDA0" | "LAMDA0" => {
                        m.lambda0 = *v;
                        m.lambda0_given = true;
                    }
                    "LAMBDA1" | "LAMDA1" => m.lambda1 = *v,
                    "RD" => m.rd = *v,
                    "RS" => m.rs = *v,
                    "CBD" => m.cbd = *v,
                    "CBS" => m.cbs = *v,
                    "IS" => m.is = *v,
                    "PB" => m.pb = *v,
                    "CGSO" => m.cgso = *v,
                    "CGDO" => m.cgdo = *v,
                    "CGBO" => m.cgbo = *v,
                    "CJ" => m.cj = *v,
                    "MJ" => m.mj = *v,
                    "CJSW" => m.cjsw = *v,
                    "MJSW" => m.mjsw = *v,
                    "TOX" => m.tox = *v,
                    "LD" => m.ld = *v,
                    "NSUB" => m.nsub = *v,
                    "U0" | "UO" => m.u0 = *v,
                    "FC" => m.fc = *v,
                    "KF" => m.kf = *v,
                    "AF" => m.af = *v,
                    _ => {} // ignore unknown params (LEVEL, TPG, NSS, XJ, RSH, etc.)
                }
            }
        }
        // If lambda0 was not given, use lambda as fallback
        if !m.lambda0_given {
            m.lambda0 = m.lambda;
        }
        m
    }

    /// Number of internal nodes needed (for series resistances RD, RS).
    pub fn internal_node_count(&self) -> usize {
        let mut count = 0;
        if self.rd > 0.0 {
            count += 1;
        }
        if self.rs > 0.0 {
            count += 1;
        }
        count
    }

    /// Compute the MOS6 companion model (n-th power I-V characteristic).
    ///
    /// `vgs`, `vds`, `vbs` are the terminal voltages (already type-sign adjusted).
    /// `betac` is KC * W / L_eff (pre-scaled by W/L).
    /// Returns `MosfetCompanion` with conductances and equivalent currents.
    pub fn companion(&self, vgs: f64, vds: f64, vbs: f64, betac: f64) -> MosfetCompanion {
        let vt = VT_NOM;
        let sign = self.mos_type.sign();

        // Determine mode: if Vds >= 0, normal; if Vds < 0, reversed (swap S/D).
        let (vgs_eff, vds_eff, vbs_eff, mode) = if vds >= 0.0 {
            (vgs, vds, vbs, 1)
        } else {
            (vgs - vds, -vds, vbs - vds, -1)
        };

        // Body effect: compute sqrt(phi - vbsvbd) for threshold
        let vbsvbd = if mode == 1 {
            vbs_eff
        } else {
            vbs_eff - vds_eff
        };
        let sarg1 = if vbsvbd <= 0.0 {
            (self.phi - vbsvbd).sqrt()
        } else {
            let s = self.phi.sqrt();
            (s - vbsvbd / (2.0 * s)).max(0.0)
        };

        // Threshold voltage (MOS6 formulation):
        // tVbi = VT0 - type * gamma * sqrt(phi) (at TNOM, temperature adjustment is identity)
        // von = tVbi * type + gamma * sarg1 - gamma1 * vbsvbd - sigma * vds_eff
        let vbi = self.vto - sign * self.gamma * self.phi.sqrt();
        let von = vbi * sign + self.gamma * sarg1 - self.gamma1 * vbsvbd - self.sigma * vds_eff;
        let vgon = vgs_eff - von;

        // Bulk-source and bulk-drain junction diode currents
        let vbd = vbs - vds;
        let (gbs, cbs_current) = bulk_diode_current(vbs, self.is, vt);
        let (gbd, cbd_current) = bulk_diode_current(vbd, self.is, vt);

        if vgon <= 0.0 {
            // Cutoff region
            let ceq_bs = cbs_current - gbs * vbs;
            let ceq_bd = cbd_current - gbd * vbd;
            return MosfetCompanion {
                gm: 0.0,
                gds: 0.0,
                gmbs: 0.0,
                gbd,
                gbs,
                cdrain: 0.0,
                ceq_d: 0.0,
                ceq_bs,
                ceq_bd,
                mode,
                vdsat: 0.0,
            };
        }

        // Body effect derivative: d(von)/d(vbs)
        let vonbm = if sarg1 > 0.0 {
            if vbsvbd <= 0.0 {
                self.gamma1 + self.gamma / (2.0 * sarg1)
            } else {
                self.gamma1 + self.gamma / (2.0 * self.phi.sqrt())
            }
        } else {
            0.0
        };

        // N-th power model: saturation voltage and current
        // vdsat = KV * vgon^NV
        // idsat = betac * vgon^NC
        let ln_vgon = vgon.ln();
        let vdsat = self.kv * (ln_vgon * self.nv).exp();
        let idsat = betac * (ln_vgon * self.nc).exp();
        let lambda_eff = self.lambda0 - self.lambda1 * vbsvbd;

        // Saturation region
        let mut cdrain = idsat * (1.0 + lambda_eff * vds_eff);
        let mut gm = cdrain * self.nc / vgon;
        let mut gds = gm * self.sigma + idsat * lambda_eff;
        let mut gmbs = gm * vonbm - idsat * self.lambda1 * vds_eff;

        if vdsat > vds_eff {
            // Linear region: scale by (2 - vdst) * vdst
            let vdst = vds_eff / vdsat;
            let vdst2 = (2.0 - vdst) * vdst;
            let vdstg = -vdst * self.nv / vgon;
            let ivdst1 = cdrain * (2.0 - 2.0 * vdst);

            cdrain *= vdst2;
            gm = gm * vdst2 + ivdst1 * vdstg;
            gds = gds * vdst2 + ivdst1 * (1.0 / vdsat + vdstg * self.sigma);
            gmbs = gmbs * vdst2 + ivdst1 * vdstg * vonbm;
        }

        // Equivalent current sources for NR linearization
        let ceq_d = cdrain - gm * vgs_eff - gds * vds_eff - gmbs * vbs_eff;
        let ceq_bs = cbs_current - gbs * vbs;
        let ceq_bd = cbd_current - gbd * vbd;

        MosfetCompanion {
            gm,
            gds,
            gmbs,
            gbd,
            gbs,
            cdrain,
            ceq_d,
            ceq_bs,
            ceq_bd,
            mode,
            vdsat,
        }
    }
}

/// Bulk junction diode current and conductance.
fn bulk_diode_current(v: f64, is: f64, vt: f64) -> (f64, f64) {
    let gmin = 1e-12;
    if v <= -3.0 * vt {
        let g = gmin;
        let i = g * v - is;
        (g, i)
    } else {
        let ev = safe_exp((v / vt).min(EXP_LIMIT));
        let g = is * ev / vt + gmin;
        let i = is * (ev - 1.0) + gmin * v;
        (g, i)
    }
}

/// Resolved node indices for a MOS6 instance in the MNA system.
#[derive(Debug, Clone)]
pub struct Mos6Instance {
    /// Element name.
    pub name: String,
    /// External drain node index (None = ground).
    pub drain_idx: Option<usize>,
    /// Gate node index (None = ground).
    pub gate_idx: Option<usize>,
    /// External source node index (None = ground).
    pub source_idx: Option<usize>,
    /// Bulk/substrate node index (None = ground).
    pub bulk_idx: Option<usize>,
    /// Internal drain prime node (when RD > 0), else same as drain_idx.
    pub drain_prime_idx: Option<usize>,
    /// Internal source prime node (when RS > 0), else same as source_idx.
    pub source_prime_idx: Option<usize>,
    /// Resolved model parameters.
    pub model: Mos6Model,
    /// Channel width (m).
    pub w: f64,
    /// Channel length (m).
    pub l: f64,
    /// Drain area for junction cap (m²).
    pub ad: f64,
    /// Source area for junction cap (m²).
    pub as_: f64,
    /// Drain perimeter for junction cap (m).
    pub pd: f64,
    /// Source perimeter for junction cap (m).
    pub ps: f64,
    /// Parallel multiplier.
    pub m: f64,
}

impl Mos6Instance {
    /// Get terminal voltages from the solution vector, handling PMOS sign.
    pub fn terminal_voltages(&self, solution: &[f64]) -> (f64, f64, f64) {
        let v_dp = self.drain_prime_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_g = self.gate_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_sp = self.source_prime_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_b = self.bulk_idx.map(|i| solution[i]).unwrap_or(0.0);

        let sign = self.model.mos_type.sign();
        let vgs = sign * (v_g - v_sp);
        let vds = sign * (v_dp - v_sp);
        let vbs = sign * (v_b - v_sp);
        (vgs, vds, vbs)
    }

    /// Effective beta: KC * W / L_eff (temperature-adjusted KC).
    pub fn betac(&self) -> f64 {
        let l_eff = (self.l - 2.0 * self.model.ld).max(1e-12);
        self.model.kc * self.w / l_eff
    }
}

/// Stamp the MOS6 companion model into the MNA matrix and RHS.
///
/// Same stamp pattern as Level 1 MOSFET: gds, gm VCCS, gmbs, junction
/// conductances, series resistances, and equivalent current sources.
pub fn stamp_mos6(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &Mos6Instance,
    comp: &MosfetCompanion,
) {
    let dp = inst.drain_prime_idx;
    let g = inst.gate_idx;
    let sp = inst.source_prime_idx;
    let b = inst.bulk_idx;

    let sign = inst.model.mos_type.sign();
    let m = inst.m;

    let (xnrm, xrev) = if comp.mode > 0 {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };

    // gds conductance between d' and s'
    crate::stamp_conductance(matrix, dp, sp, m * comp.gds);

    // gm VCCS: Vgs controls current from s' to d'
    let gm_scaled = m * comp.gm;
    if let Some(d) = dp {
        if let Some(gate) = g {
            matrix.add(d, gate, (xnrm - xrev) * gm_scaled);
        }
        if let Some(s) = sp {
            matrix.add(d, s, -(xnrm - xrev) * gm_scaled);
        }
    }
    if let Some(s) = sp {
        if let Some(gate) = g {
            matrix.add(s, gate, -(xnrm - xrev) * gm_scaled);
        }
        matrix.add(s, s, (xnrm - xrev) * gm_scaled);
    }

    // gmbs: Vbs controls current from s' to d'
    let gmbs_scaled = m * comp.gmbs;
    if let Some(d) = dp {
        if let Some(bulk) = b {
            matrix.add(d, bulk, (xnrm - xrev) * gmbs_scaled);
        }
        if let Some(s) = sp {
            matrix.add(d, s, -(xnrm - xrev) * gmbs_scaled);
        }
    }
    if let Some(s) = sp {
        if let Some(bulk) = b {
            matrix.add(s, bulk, -(xnrm - xrev) * gmbs_scaled);
        }
        matrix.add(s, s, (xnrm - xrev) * gmbs_scaled);
    }

    // gbd conductance between b and d'
    crate::stamp_conductance(matrix, b, dp, m * comp.gbd);

    // gbs conductance between b and s'
    crate::stamp_conductance(matrix, b, sp, m * comp.gbs);

    // Series resistances
    if inst.model.rd > 0.0 {
        let grd = 1.0 / inst.model.rd;
        crate::stamp_conductance(matrix, inst.drain_idx, dp, m * grd);
    }
    if inst.model.rs > 0.0 {
        let grs = 1.0 / inst.model.rs;
        crate::stamp_conductance(matrix, inst.source_idx, sp, m * grs);
    }

    // Equivalent current sources on the RHS
    let ceq_d = sign * m * comp.ceq_d;
    let ceq_bs = sign * m * comp.ceq_bs;
    let ceq_bd = sign * m * comp.ceq_bd;

    if let Some(d) = dp {
        rhs[d] -= ceq_d + ceq_bd;
    }
    if let Some(s) = sp {
        rhs[s] += ceq_d + ceq_bs;
    }
    if let Some(bulk) = b {
        rhs[bulk] += ceq_bd + ceq_bs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use thevenin_types::Param;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    fn make_nmos_model() -> Mos6Model {
        let model_def = ModelDef {
            name: "N10L5".to_string(),
            kind: "NMOS".to_string(),
            params: vec![
                Param {
                    name: "LEVEL".to_string(),
                    value: Expr::Num(6.0),
                },
                Param {
                    name: "KC".to_string(),
                    value: Expr::Num(3.8921e-05),
                },
                Param {
                    name: "NC".to_string(),
                    value: Expr::Num(1.1739),
                },
                Param {
                    name: "KV".to_string(),
                    value: Expr::Num(0.91602),
                },
                Param {
                    name: "NV".to_string(),
                    value: Expr::Num(0.87225),
                },
                Param {
                    name: "LAMBDA0".to_string(),
                    value: Expr::Num(0.013333),
                },
                Param {
                    name: "LAMBDA1".to_string(),
                    value: Expr::Num(0.0046901),
                },
                Param {
                    name: "VT0".to_string(),
                    value: Expr::Num(0.69486),
                },
                Param {
                    name: "GAMMA".to_string(),
                    value: Expr::Num(0.60309),
                },
                Param {
                    name: "PHI".to_string(),
                    value: Expr::Num(1.0),
                },
            ],
        };
        Mos6Model::from_model_def(&model_def)
    }

    #[test]
    fn test_default_mos6_model() {
        let m = Mos6Model::new(MosfetType::Nmos);
        assert_eq!(m.vto, 0.0);
        assert_eq!(m.kc, 0.0);
        assert_eq!(m.nc, 2.0);
        assert_eq!(m.kv, 1.0);
        assert_eq!(m.nv, 1.0);
        assert_eq!(m.lambda0, 0.0);
        assert_eq!(m.lambda1, 0.0);
        assert_eq!(m.gamma1, 0.0);
        assert_eq!(m.sigma, 0.0);
    }

    #[test]
    fn test_from_model_def() {
        let m = make_nmos_model();
        assert_eq!(m.mos_type, MosfetType::Nmos);
        assert_abs_diff_eq!(m.kc, 3.8921e-05, epsilon = 1e-15);
        assert_abs_diff_eq!(m.nc, 1.1739, epsilon = 1e-15);
        assert_abs_diff_eq!(m.kv, 0.91602, epsilon = 1e-15);
        assert_abs_diff_eq!(m.nv, 0.87225, epsilon = 1e-15);
        assert_abs_diff_eq!(m.lambda0, 0.013333, epsilon = 1e-15);
        assert_abs_diff_eq!(m.lambda1, 0.0046901, epsilon = 1e-15);
        assert_abs_diff_eq!(m.vto, 0.69486, epsilon = 1e-15);
        assert!(m.lambda0_given);
    }

    #[test]
    fn test_lambda_fallback() {
        let model_def = ModelDef {
            name: "TEST".to_string(),
            kind: "NMOS".to_string(),
            params: vec![Param {
                name: "LAMBDA".to_string(),
                value: Expr::Num(0.05),
            }],
        };
        let m = Mos6Model::from_model_def(&model_def);
        assert_abs_diff_eq!(m.lambda0, 0.05, epsilon = 1e-15);
        assert!(!m.lambda0_given);
    }

    #[test]
    fn test_companion_cutoff() {
        let mut m = make_nmos_model();
        m.vto = 1.0;
        // Vgs = 0.5 < Vto → cutoff (vgon ≤ 0)
        let comp = m.companion(0.5, 5.0, 0.0, 1e-4);
        assert_abs_diff_eq!(comp.cdrain, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(comp.gm, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(comp.gds, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_companion_saturation() {
        let m = make_nmos_model();
        // Vgs=3, Vds=5, Vbs=0 → saturation (vgon > 0, vdsat < vds)
        let comp = m.companion(3.0, 5.0, 0.0, 1e-4);
        assert!(
            comp.cdrain > 0.0,
            "drain current should be positive in saturation"
        );
        assert!(comp.gm > 0.0, "transconductance should be positive");
        assert!(comp.gds > 0.0, "output conductance should be positive");
        assert_eq!(comp.mode, 1);
    }

    #[test]
    fn test_companion_linear() {
        let m = make_nmos_model();
        // Vgs=3, Vds=0.1, Vbs=0 → linear (vds < vdsat)
        let comp = m.companion(3.0, 0.1, 0.0, 1e-4);
        assert!(
            comp.cdrain > 0.0,
            "drain current should be positive in linear"
        );
        assert_eq!(comp.mode, 1);
    }

    #[test]
    fn test_companion_reversed_mode() {
        let m = make_nmos_model();
        let comp = m.companion(3.0, -1.0, 0.0, 1e-4);
        assert_eq!(comp.mode, -1);
        assert!(comp.cdrain > 0.0);
    }

    #[test]
    fn test_pmos_model() {
        let model_def = ModelDef {
            name: "P12L5".to_string(),
            kind: "PMOS".to_string(),
            params: vec![
                Param {
                    name: "KC".to_string(),
                    value: Expr::Num(6.42696e-06),
                },
                Param {
                    name: "NC".to_string(),
                    value: Expr::Num(1.6536),
                },
                Param {
                    name: "KV".to_string(),
                    value: Expr::Num(0.92145),
                },
                Param {
                    name: "NV".to_string(),
                    value: Expr::Num(0.88345),
                },
                Param {
                    name: "LAMBDA0".to_string(),
                    value: Expr::Num(0.018966),
                },
                Param {
                    name: "LAMBDA1".to_string(),
                    value: Expr::Num(0.0084012),
                },
                Param {
                    name: "VT0".to_string(),
                    value: Expr::Num(-0.60865),
                },
                Param {
                    name: "GAMMA".to_string(),
                    value: Expr::Num(0.89213),
                },
                Param {
                    name: "PHI".to_string(),
                    value: Expr::Num(1.0),
                },
            ],
        };
        let m = Mos6Model::from_model_def(&model_def);
        assert_eq!(m.mos_type, MosfetType::Pmos);
        assert_abs_diff_eq!(m.vto, -0.60865, epsilon = 1e-15);
    }

    #[test]
    fn test_internal_node_count() {
        let mut m = Mos6Model::new(MosfetType::Nmos);
        assert_eq!(m.internal_node_count(), 0);
        m.rd = 10.0;
        assert_eq!(m.internal_node_count(), 1);
        m.rs = 5.0;
        assert_eq!(m.internal_node_count(), 2);
    }
}
