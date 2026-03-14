//! BJT (Bipolar Junction Transistor) device model — Ebers-Moll / Gummel-Poon level 1.
//!
//! Implements the standard SPICE BJT model with NR companion linearization.
//! Supports NPN and PNP types with Early effect (VAF/VAR) and high-injection
//! rolloff (IKF/IKR).

use thevenin_types::{Expr, ModelDef, Param};

use crate::diode::{VT_NOM, pnjlim, vcrit};
use crate::physics::{KBOQ, safe_exp};

/// Default nominal temperature (°C) matching SPICE TNOM.
const TNOM_DEFAULT_C: f64 = 27.0;

/// BJT polarity: NPN or PNP.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BjtType {
    Npn,
    Pnp,
}

impl BjtType {
    /// Type multiplier: +1 for NPN, -1 for PNP.
    pub fn sign(self) -> f64 {
        match self {
            BjtType::Npn => 1.0,
            BjtType::Pnp => -1.0,
        }
    }
}

/// BJT model parameters matching ngspice level 1 (Gummel-Poon) defaults.
#[derive(Debug, Clone)]
pub struct BjtModel {
    /// NPN or PNP.
    pub bjt_type: BjtType,
    /// Saturation current (default 1e-16 A).
    pub is: f64,
    /// Forward current gain (default 100).
    pub bf: f64,
    /// Reverse current gain (default 1).
    pub br: f64,
    /// Forward emission coefficient (default 1.0).
    pub nf: f64,
    /// Reverse emission coefficient (default 1.0).
    pub nr: f64,
    /// B-E leakage saturation current (default 0).
    pub ise: f64,
    /// B-E leakage emission coefficient (default 1.5).
    pub ne: f64,
    /// B-C leakage saturation current (default 0).
    pub isc: f64,
    /// B-C leakage emission coefficient (default 2.0).
    pub nc: f64,
    /// Forward Early voltage (default: infinity = no Early effect).
    pub vaf: f64,
    /// Reverse Early voltage (default: infinity).
    pub var: f64,
    /// Forward high-current rolloff (default: infinity).
    pub ikf: f64,
    /// Reverse high-current rolloff (default: infinity).
    pub ikr: f64,
    /// Base resistance (default 0 Ω).
    pub rb: f64,
    /// Collector resistance (default 0 Ω).
    pub rc: f64,
    /// Emitter resistance (default 0 Ω).
    pub re: f64,
    /// B-E zero-bias junction capacitance (default 0 F).
    pub cje: f64,
    /// B-E built-in potential (default 0.75 V).
    pub vje: f64,
    /// B-E grading coefficient (default 0.33).
    pub mje: f64,
    /// B-C zero-bias junction capacitance (default 0 F).
    pub cjc: f64,
    /// B-C built-in potential (default 0.75 V).
    pub vjc: f64,
    /// B-C grading coefficient (default 0.33).
    pub mjc: f64,
    /// Collector-substrate capacitance (default 0 F).
    pub cjs: f64,
    /// Forward transit time (default 0 s).
    pub tf: f64,
    /// Reverse transit time (default 0 s).
    pub tr: f64,
    /// Fraction of B-C cap connected to internal base (default 1).
    pub xcjc: f64,
    /// Flicker noise coefficient (default 0).
    pub kf: f64,
    /// Flicker noise exponent (default 1).
    pub af: f64,
    /// Energy gap (default 1.11 eV, silicon).
    pub eg: f64,
    /// Saturation current temperature exponent (default 3.0).
    pub xti: f64,
    /// Forward-bias depletion capacitance coefficient (default 0.5).
    pub fc: f64,
    /// Minimum base resistance at high injection (default = rb).
    pub rbm: f64,
    /// Nominal temperature for model parameters in °C (default 27.0).
    pub tnom: f64,
}

impl BjtModel {
    /// Create a new BjtModel with default parameters for the given type.
    pub fn new(bjt_type: BjtType) -> Self {
        Self {
            bjt_type,
            is: 1e-16,
            bf: 100.0,
            br: 1.0,
            nf: 1.0,
            nr: 1.0,
            ise: 0.0,
            ne: 1.5,
            isc: 0.0,
            nc: 2.0,
            vaf: f64::INFINITY,
            var: f64::INFINITY,
            ikf: f64::INFINITY,
            ikr: f64::INFINITY,
            rb: 0.0,
            rc: 0.0,
            re: 0.0,
            cje: 0.0,
            vje: 0.75,
            mje: 0.33,
            cjc: 0.0,
            vjc: 0.75,
            mjc: 0.33,
            cjs: 0.0,
            tf: 0.0,
            tr: 0.0,
            xcjc: 1.0,
            kf: 0.0,
            af: 1.0,
            eg: 1.11,
            xti: 3.0,
            fc: 0.5,
            rbm: 0.0, // sentinel: set to rb after parsing
            tnom: TNOM_DEFAULT_C,
        }
    }

    /// Create a `BjtModel` from a netlist `.model` definition.
    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let bjt_type = if model_def.kind.to_uppercase() == "PNP" {
            BjtType::Pnp
        } else {
            BjtType::Npn
        };
        let mut m = Self::new(bjt_type);
        let mut rbm_given = false;
        for p in &model_def.params {
            if let Expr::Num(v) = &p.value {
                match p.name.to_uppercase().as_str() {
                    "IS" => m.is = *v,
                    "BF" => m.bf = *v,
                    "BR" => m.br = *v,
                    "NF" => m.nf = *v,
                    "NR" => m.nr = *v,
                    "ISE" => m.ise = *v,
                    "NE" => m.ne = *v,
                    "ISC" => m.isc = *v,
                    "NC" => m.nc = *v,
                    "VAF" | "VA" => m.vaf = *v,
                    "VAR" => m.var = *v,
                    "IKF" | "JBF" => m.ikf = *v,
                    "IKR" | "JBR" => m.ikr = *v,
                    "RB" => m.rb = *v,
                    "RBM" => { m.rbm = *v; rbm_given = true; }
                    "RC" => m.rc = *v,
                    "RE" => m.re = *v,
                    "CJE" => m.cje = *v,
                    "VJE" | "PE" => m.vje = *v,
                    "MJE" | "ME" => m.mje = *v,
                    "CJC" => m.cjc = *v,
                    "VJC" | "PC" => m.vjc = *v,
                    "MJC" | "MC" => m.mjc = *v,
                    "CJS" | "CCS" => m.cjs = *v,
                    "TF" => m.tf = *v,
                    "TR" => m.tr = *v,
                    "XCJC" => m.xcjc = *v,
                    "KF" => m.kf = *v,
                    "AF" => m.af = *v,
                    "EG" => m.eg = *v,
                    "XTI" => m.xti = *v,
                    "FC" => m.fc = *v,
                    "TNOM" => m.tnom = *v,
                    _ => {} // ignore unknown params (LEVEL, etc.)
                }
            }
        }
        // RBM defaults to RB if not explicitly given (SPICE convention).
        if !rbm_given {
            m.rbm = m.rb;
        }
        m
    }

    /// Apply instance-level parameter overrides.
    pub fn with_instance_params(mut self, params: &[Param]) -> Self {
        for p in params {
            if let Expr::Num(v) = &p.value {
                match p.name.to_uppercase().as_str() {
                    "IS" => self.is = *v,
                    "BF" => self.bf = *v,
                    "BR" => self.br = *v,
                    _ => {}
                }
            }
        }
        self
    }

    /// Number of internal nodes needed (for series resistances RB, RC, RE).
    pub fn internal_node_count(&self) -> usize {
        let mut count = 0;
        if self.rb > 0.0 {
            count += 1;
        }
        if self.rc > 0.0 {
            count += 1;
        }
        if self.re > 0.0 {
            count += 1;
        }
        count
    }

    /// Critical voltage for B-E junction limiting.
    pub fn vcrit_be(&self) -> f64 {
        vcrit(self.nf * VT_NOM, self.is)
    }

    /// Critical voltage for B-C junction limiting.
    pub fn vcrit_bc(&self) -> f64 {
        vcrit(self.nr * VT_NOM, self.is)
    }

    /// Internal: compute companion model given explicit thermal voltage and saturation current.
    fn companion_impl(&self, vbe: f64, vbc: f64, vt: f64, is: f64) -> BjtCompanion {
        // Junction currents and conductances (Shockley diode equations)
        let vte_f = self.nf * vt;
        let vte_r = self.nr * vt;

        let exp_be = safe_exp(vbe / vte_f);
        let cbe = is * (exp_be - 1.0);
        let gbe = is / vte_f * exp_be;

        let exp_bc = safe_exp(vbc / vte_r);
        let cbc = is * (exp_bc - 1.0);
        let gbc = is / vte_r * exp_bc;

        // Leakage currents
        let (cben, gben) = if self.ise > 0.0 {
            let vte_e = self.ne * vt;
            let exp_ben = safe_exp(vbe / vte_e);
            (self.ise * (exp_ben - 1.0), self.ise / vte_e * exp_ben)
        } else {
            (0.0, 0.0)
        };

        let (cbcn, gbcn) = if self.isc > 0.0 {
            let vte_c = self.nc * vt;
            let exp_bcn = safe_exp(vbc / vte_c);
            (self.isc * (exp_bcn - 1.0), self.isc / vte_c * exp_bcn)
        } else {
            (0.0, 0.0)
        };

        // Base charge normalization (Gummel-Poon)
        let inv_vaf = if self.vaf.is_finite() { 1.0 / self.vaf } else { 0.0 };
        let inv_var = if self.var.is_finite() { 1.0 / self.var } else { 0.0 };
        let inv_ikf = if self.ikf.is_finite() { 1.0 / self.ikf } else { 0.0 };
        let inv_ikr = if self.ikr.is_finite() { 1.0 / self.ikr } else { 0.0 };

        let q1_denom = 1.0 - inv_vaf * vbc - inv_var * vbe;
        let q1 = if q1_denom > 0.0 { 1.0 / q1_denom } else { 1.0 / 1e-6 };

        let (qb, dqbdve, dqbdvc) = if inv_ikf == 0.0 && inv_ikr == 0.0 {
            let dve = q1 * q1 * inv_var;
            let dvc = q1 * q1 * inv_vaf;
            (q1, dve, dvc)
        } else {
            let q2 = inv_ikf * cbe + inv_ikr * cbc;
            let arg = (1.0 + 4.0 * q2).max(0.0);
            let sqarg = arg.sqrt();
            let qb_val = q1 * (1.0 + sqarg) / 2.0;
            let dve = q1 * (qb_val * inv_var / q1 + inv_ikf * gbe / sqarg);
            let dvc = q1 * (qb_val * inv_vaf / q1 + inv_ikr * gbc / sqarg);
            (qb_val, dve, dvc)
        };

        // Terminal currents
        let cc = (cbe - cbc) / qb - cbc / self.br - cbcn;
        let cb = cbe / self.bf + cben + cbc / self.br + cbcn;

        // Small-signal conductances (ngspice formulation)
        let gpi = gbe / self.bf + gben;
        let gmu = gbc / self.br + gbcn;
        let go = (gbc + (cbe - cbc) * dqbdvc / qb) / qb;
        let gm = (gbe - (cbe - cbc) * dqbdve / qb) / qb - go;

        // Equivalent current sources for the RHS (ngspice formulation)
        let ceqbe = cc + cb - vbe * (gm + go + gpi) + vbc * go;
        let ceqbc = -cc + vbe * (gm + go) - vbc * (gmu + go);

        BjtCompanion { gpi, gmu, gm, go, ceqbe, ceqbc, cc, cb, qb }
    }

    /// Compute the BJT operating point and NR companion model at nominal temperature.
    ///
    /// Uses `VT_NOM` thermal voltage and `self.is` without temperature scaling.
    /// Junction voltages should already be limited via `pnjlim` before calling.
    pub fn companion(&self, vbe: f64, vbc: f64) -> BjtCompanion {
        self.companion_impl(vbe, vbc, VT_NOM, self.is)
    }

    /// Compute the BJT companion model at a specified device temperature in Kelvin.
    ///
    /// Applies IS temperature scaling:
    /// `IS_t = IS * (T/Tnom)^(xti/nf) * exp(eg/(nf*KB) * (1/Tnom - 1/T))`
    /// and uses temperature-scaled thermal voltage `VT = KB * T`.
    pub fn companion_at(&self, vbe: f64, vbc: f64, temp_k: f64) -> BjtCompanion {
        let vt = KBOQ * temp_k;
        let tnom_k = self.tnom + 273.15;
        let ratio = temp_k / tnom_k;
        let is_t = if (ratio - 1.0).abs() < 1e-9 {
            self.is
        } else {
            let xti_nf = self.xti / self.nf;
            let eg_arg = (self.eg / (self.nf * KBOQ)) * (1.0 / tnom_k - 1.0 / temp_k);
            self.is * ratio.powf(xti_nf) * eg_arg.exp()
        };
        self.companion_impl(vbe, vbc, vt, is_t)
    }

    /// B-E junction capacitance at voltage v.
    pub fn cap_be(&self, v: f64) -> f64 {
        junction_cap(v, self.cje, self.vje, self.mje, self.fc)
    }

    /// B-C junction capacitance at voltage v.
    pub fn cap_bc(&self, v: f64) -> f64 {
        junction_cap(v, self.cjc, self.vjc, self.mjc, self.fc)
    }

    /// Limit B-E junction voltage.
    pub fn limit_vbe(&self, v_new: f64, v_old: f64) -> f64 {
        pnjlim(v_new, v_old, self.nf * VT_NOM, self.vcrit_be())
    }

    /// Limit B-C junction voltage.
    pub fn limit_vbc(&self, v_new: f64, v_old: f64) -> f64 {
        pnjlim(v_new, v_old, self.nr * VT_NOM, self.vcrit_bc())
    }
}

/// Junction depletion capacitance with forward-bias linear extrapolation.
fn junction_cap(v: f64, cj0: f64, vj: f64, m: f64, fc: f64) -> f64 {
    if cj0 <= 0.0 {
        return 0.0;
    }
    let fc_vj = fc * vj;
    if v < fc_vj {
        cj0 / (1.0 - v / vj).powf(m)
    } else {
        let base = cj0 / (1.0 - fc).powf(m);
        let slope = base * m / (vj * (1.0 - fc));
        base + slope * (v - fc_vj)
    }
}

/// NR companion model result for a BJT at an operating point.
#[derive(Debug, Clone)]
pub struct BjtCompanion {
    /// Base-emitter input conductance (dI_b/dV_be).
    pub gpi: f64,
    /// Base-collector junction conductance (dI_b/dV_bc).
    pub gmu: f64,
    /// Transconductance.
    pub gm: f64,
    /// Output conductance (c'-e').
    pub go: f64,
    /// Equivalent current source at emitter prime (ngspice ceqbe).
    pub ceqbe: f64,
    /// Equivalent current source at collector prime (ngspice ceqbc).
    pub ceqbc: f64,
    /// Collector current (for output).
    pub cc: f64,
    /// Base current (for output).
    pub cb: f64,
    /// Normalized base charge QB (used for base-resistance modulation).
    pub qb: f64,
}

/// Resolved node indices for a BJT instance in the MNA system.
#[derive(Debug, Clone)]
pub struct BjtInstance {
    /// BJT element name.
    pub name: String,
    /// External collector node index (None = ground).
    pub col_idx: Option<usize>,
    /// External base node index (None = ground).
    pub base_idx: Option<usize>,
    /// External emitter node index (None = ground).
    pub emit_idx: Option<usize>,
    /// Internal base prime node (when RB > 0), else same as base_idx.
    pub base_prime_idx: Option<usize>,
    /// Internal collector prime node (when RC > 0), else same as col_idx.
    pub col_prime_idx: Option<usize>,
    /// Internal emitter prime node (when RE > 0), else same as emit_idx.
    pub emit_prime_idx: Option<usize>,
    /// Resolved BJT model parameters.
    pub model: BjtModel,
    /// Area scaling factor (default 1.0).
    pub area: f64,
    /// Base area factor for sensitivity (default 1.0).
    pub areab: f64,
    /// Collector area factor for sensitivity (default 1.0).
    pub areac: f64,
    /// Parallel multiplier (default 1.0).
    pub m: f64,
    /// Instance temperature override in °C (NaN = use circuit temperature).
    pub temp: f64,
}

impl BjtInstance {
    /// Get junction voltages from the solution vector, handling PNP sign.
    pub fn junction_voltages(&self, solution: &[f64]) -> (f64, f64) {
        let v_bp = self.base_prime_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_cp = self.col_prime_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_ep = self.emit_prime_idx.map(|i| solution[i]).unwrap_or(0.0);

        let sign = self.model.bjt_type.sign();
        let vbe = sign * (v_bp - v_ep);
        let vbc = sign * (v_bp - v_cp);
        (vbe, vbc)
    }
}

/// Stamp the BJT companion model into the MNA matrix and RHS.
///
/// Follows the ngspice decomposition:
/// 1. gpi conductance between b' and e'
/// 2. gmu conductance between b' and c'
/// 3. gm VCCS: V_be controls current from e' to c'
/// 4. go conductance between c' and e'
/// 5. Series resistances (RB, RC, RE) as conductances between external and internal nodes
/// 6. Equivalent current sources on the RHS
pub fn stamp_bjt(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &BjtInstance,
    comp: &BjtCompanion,
) {
    let bp = inst.base_prime_idx;
    let cp = inst.col_prime_idx;
    let ep = inst.emit_prime_idx;

    let sign = inst.model.bjt_type.sign();
    let m = inst.m * inst.area;

    // 1. gpi conductance between b' and e'
    crate::stamp_conductance(matrix, bp, ep, m * comp.gpi);

    // 2. gmu conductance between b' and c'
    crate::stamp_conductance(matrix, bp, cp, m * comp.gmu);

    // 3. gm VCCS: V_be controls current from e' to c'
    let gm_scaled = m * comp.gm;
    if let Some(c) = cp {
        if let Some(b) = bp {
            matrix.add(c, b, gm_scaled);
        }
        if let Some(e) = ep {
            matrix.add(c, e, -gm_scaled);
        }
    }
    if let Some(e) = ep {
        if let Some(b) = bp {
            matrix.add(e, b, -gm_scaled);
        }
        matrix.add(e, e, gm_scaled);
    }

    // 4. go conductance between c' and e'
    crate::stamp_conductance(matrix, cp, ep, m * comp.go);

    // 5. Series resistances (stamped as conductances between external and internal)
    if inst.model.rb > 0.0 {
        // QB-modulated base resistance: rbb = rbm + (rb - rbm) / QB
        // When rbm == rb (default), this reduces to rbb = rb regardless of QB.
        let qb_safe = comp.qb.max(0.01);
        let rbb = inst.model.rbm + (inst.model.rb - inst.model.rbm) / qb_safe;
        let gx = 1.0 / rbb.max(1e-12);
        crate::stamp_conductance(matrix, inst.base_idx, bp, m * gx);
    }
    if inst.model.rc > 0.0 {
        let gcpr = 1.0 / inst.model.rc;
        crate::stamp_conductance(matrix, inst.col_idx, cp, m * gcpr);
    }
    if inst.model.re > 0.0 {
        let gepr = 1.0 / inst.model.re;
        crate::stamp_conductance(matrix, inst.emit_idx, ep, m * gepr);
    }

    // 6. Equivalent current sources on the RHS
    let ceqbe = sign * m * comp.ceqbe;
    let ceqbc = sign * m * comp.ceqbc;

    if let Some(b) = bp {
        rhs[b] -= ceqbe + ceqbc;
    }
    if let Some(c) = cp {
        rhs[c] += ceqbc;
    }
    if let Some(e) = ep {
        rhs[e] += ceqbe;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    #[test]
    fn test_default_bjt_model() {
        let m = BjtModel::new(BjtType::Npn);
        assert_eq!(m.is, 1e-16);
        assert_eq!(m.bf, 100.0);
        assert_eq!(m.br, 1.0);
        assert_eq!(m.nf, 1.0);
        assert_eq!(m.nr, 1.0);
        assert_eq!(m.rb, 0.0);
    }

    #[test]
    fn test_from_model_def_npn() {
        let model_def = ModelDef {
            name: "QND".to_string(),
            kind: "NPN".to_string(),
            params: vec![
                Param {
                    name: "BF".to_string(),
                    value: Expr::Num(50.0),
                },
                Param {
                    name: "RB".to_string(),
                    value: Expr::Num(70.0),
                },
                Param {
                    name: "VA".to_string(),
                    value: Expr::Num(50.0),
                },
                Param {
                    name: "IS".to_string(),
                    value: Expr::Num(1e-15),
                },
            ],
        };
        let m = BjtModel::from_model_def(&model_def);
        assert_eq!(m.bjt_type, BjtType::Npn);
        assert_abs_diff_eq!(m.bf, 50.0, epsilon = 1e-15);
        assert_abs_diff_eq!(m.rb, 70.0, epsilon = 1e-15);
        assert_abs_diff_eq!(m.vaf, 50.0, epsilon = 1e-15);
        assert_abs_diff_eq!(m.is, 1e-15, epsilon = 1e-30);
    }

    #[test]
    fn test_from_model_def_pnp() {
        let model_def = ModelDef {
            name: "P1".to_string(),
            kind: "PNP".to_string(),
            params: vec![Param {
                name: "BF".to_string(),
                value: Expr::Num(80.0),
            }],
        };
        let m = BjtModel::from_model_def(&model_def);
        assert_eq!(m.bjt_type, BjtType::Pnp);
        assert_abs_diff_eq!(m.bf, 80.0, epsilon = 1e-15);
    }

    #[test]
    fn test_companion_forward_active() {
        let m = BjtModel::new(BjtType::Npn);
        // Forward active: Vbe=0.7V, Vbc=-5V (collector higher than base)
        let comp = m.companion(0.7, -5.0);

        // In forward active, gm should be significant (transconductance)
        assert!(comp.gm > 0.0, "gm should be positive, got {}", comp.gm);

        // gpi = gbe/BF, should be much smaller than gm
        assert!(
            comp.gpi < comp.gm,
            "gpi={} should be < gm={}",
            comp.gpi,
            comp.gm
        );

        // go should be small (reverse junction has negligible current)
        assert!(comp.go >= 0.0, "go should be non-negative, got {}", comp.go);

        // Collector current should be positive and significant
        assert!(comp.cc > 0.0, "Ic should be positive, got {}", comp.cc);

        // Base current should be positive and Ic/BF approximately
        assert!(comp.cb > 0.0, "Ib should be positive, got {}", comp.cb);
        let beta = comp.cc / comp.cb;
        assert!(
            beta > 50.0 && beta < 150.0,
            "beta={} should be near BF=100",
            beta
        );
    }

    #[test]
    fn test_companion_cutoff() {
        let m = BjtModel::new(BjtType::Npn);
        // Cutoff: both junctions reverse biased
        let comp = m.companion(-1.0, -1.0);
        assert!(
            comp.cc.abs() < 1e-10,
            "Ic in cutoff should be ~0, got {}",
            comp.cc
        );
        assert!(
            comp.cb.abs() < 1e-10,
            "Ib in cutoff should be ~0, got {}",
            comp.cb
        );
    }

    #[test]
    fn test_companion_with_early_effect() {
        let mut m = BjtModel::new(BjtType::Npn);
        m.vaf = 50.0; // Forward Early voltage
        let comp1 = m.companion(0.7, -5.0);
        let comp2 = m.companion(0.7, -10.0);
        // With Early effect, larger |Vce| should give larger Ic
        // Vce = Vcb + Vbe, so Vcb = -Vbc
        // comp2 has larger Vce (10+0.7 vs 5+0.7)
        assert!(
            comp2.cc > comp1.cc,
            "Ic should increase with Vce (Early effect): {} vs {}",
            comp2.cc,
            comp1.cc
        );
    }

    #[test]
    fn test_pnp_type_sign() {
        assert_eq!(BjtType::Npn.sign(), 1.0);
        assert_eq!(BjtType::Pnp.sign(), -1.0);
    }

    #[test]
    fn test_internal_node_count() {
        let mut m = BjtModel::new(BjtType::Npn);
        assert_eq!(m.internal_node_count(), 0);
        m.rb = 100.0;
        assert_eq!(m.internal_node_count(), 1);
        m.rc = 40.0;
        m.re = 1.0;
        assert_eq!(m.internal_node_count(), 3);
    }
}
