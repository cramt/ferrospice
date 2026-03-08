//! MOSFET Level 1 (Shichman-Hodges) device model.
//!
//! Implements the standard SPICE Level 1 MOSFET model with NR companion
//! linearization. Supports NMOS and PMOS types with body effect (GAMMA),
//! channel length modulation (LAMBDA), and bulk junction diodes.

use ferrospice_netlist::{Expr, ModelDef};

use crate::diode::VT_NOM;

/// Maximum exponent argument to prevent overflow.
const EXP_LIMIT: f64 = 500.0;

/// Safe exponential that clamps the argument to avoid overflow.
fn safe_exp(x: f64) -> f64 {
    if x > EXP_LIMIT {
        EXP_LIMIT.exp()
    } else {
        x.exp()
    }
}

/// MOSFET polarity: NMOS or PMOS.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MosfetType {
    Nmos,
    Pmos,
}

impl MosfetType {
    /// Type multiplier: +1 for NMOS, -1 for PMOS.
    pub fn sign(self) -> f64 {
        match self {
            MosfetType::Nmos => 1.0,
            MosfetType::Pmos => -1.0,
        }
    }
}

/// MOSFET Level 1 model parameters matching ngspice defaults.
#[derive(Debug, Clone)]
pub struct MosfetModel {
    /// NMOS or PMOS.
    pub mos_type: MosfetType,
    /// Threshold voltage (default 0 V).
    pub vto: f64,
    /// Transconductance parameter (default 2e-5 A/V²).
    pub kp: f64,
    /// Body effect coefficient (default 0).
    pub gamma: f64,
    /// Surface potential / 2*phi_f (default 0.6 V).
    pub phi: f64,
    /// Channel length modulation (default 0 1/V).
    pub lambda: f64,
    /// Drain resistance (default 0 Ω).
    pub rd: f64,
    /// Source resistance (default 0 Ω).
    pub rs: f64,
    /// Bulk-drain zero-bias junction capacitance (default 0 F).
    pub cbd: f64,
    /// Bulk-source zero-bias junction capacitance (default 0 F).
    pub cbs: f64,
    /// Bulk junction saturation current (default 1e-14 A).
    pub is: f64,
    /// Bulk junction potential (default 0.8 V).
    pub pb: f64,
    /// Gate-source overlap capacitance per unit width (default 0 F/m).
    pub cgso: f64,
    /// Gate-drain overlap capacitance per unit width (default 0 F/m).
    pub cgdo: f64,
    /// Gate-bulk overlap capacitance per unit length (default 0 F/m).
    pub cgbo: f64,
    /// Bottom junction capacitance per unit area (default 0 F/m²).
    pub cj: f64,
    /// Bottom junction grading coefficient (default 0.5).
    pub mj: f64,
    /// Sidewall junction capacitance per unit length (default 0 F/m).
    pub cjsw: f64,
    /// Sidewall junction grading coefficient (default 0.5).
    pub mjsw: f64,
    /// Oxide thickness (default 1e-7 m).
    pub tox: f64,
    /// Lateral diffusion (default 0 m).
    pub ld: f64,
    /// Substrate doping (default 0 1/cm³).
    pub nsub: f64,
    /// Surface mobility (default 600 cm²/V·s).
    pub u0: f64,
    /// Forward bias depletion cap coefficient (default 0.5).
    pub fc: f64,
    /// Flicker noise coefficient (default 0).
    pub kf: f64,
    /// Flicker noise exponent (default 1).
    pub af: f64,
}

impl MosfetModel {
    /// Create a new MosfetModel with default parameters for the given type.
    pub fn new(mos_type: MosfetType) -> Self {
        Self {
            mos_type,
            vto: 0.0,
            kp: 2e-5,
            gamma: 0.0,
            phi: 0.6,
            lambda: 0.0,
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
        }
    }

    /// Create a `MosfetModel` from a netlist `.model` definition.
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
                    "KP" => m.kp = *v,
                    "GAMMA" => m.gamma = *v,
                    "PHI" => m.phi = *v,
                    "LAMBDA" => m.lambda = *v,
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
                    _ => {} // ignore unknown params (LEVEL, TPG, etc.)
                }
            }
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

    /// Compute the MOSFET operating point and NR companion model.
    ///
    /// `vgs`, `vds`, `vbs` are the terminal voltages (already adjusted for type sign
    /// and mode). Returns `MosfetCompanion` with conductances and equivalent currents.
    pub fn companion(&self, vgs: f64, vds: f64, vbs: f64) -> MosfetCompanion {
        let vt = VT_NOM;

        // Determine mode: if Vds >= 0, normal; if Vds < 0, reversed (swap S/D).
        let (vgs_eff, vds_eff, vbs_eff, mode) = if vds >= 0.0 {
            (vgs, vds, vbs, 1)
        } else {
            // Reverse mode: swap source and drain
            // Vgd = Vgs - Vds, Vds_eff = -Vds, Vbd = Vbs - Vds
            (vgs - vds, -vds, vbs - vds, -1)
        };

        // Body effect: compute threshold voltage shift.
        let sarg = if vbs_eff <= 0.0 {
            (self.phi - vbs_eff).sqrt()
        } else {
            let s = self.phi.sqrt();
            (s - vbs_eff / (2.0 * s)).max(0.0)
        };

        let von = self.vto + self.gamma * sarg;
        let vgst = vgs_eff - von;
        let vdsat = vgst.max(0.0);

        // Body effect derivative: d(Von)/d(Vbs).
        let arg = if sarg > 0.0 {
            self.gamma / (2.0 * sarg)
        } else {
            0.0
        };

        // Effective beta: KP * W/L (W and L are applied through the instance).
        let beta = self.kp;

        // Drain current and small-signal conductances.
        let (cdrain, gm, gds, gmbs);
        if vgst <= 0.0 {
            // Cutoff region
            cdrain = 0.0;
            gm = 0.0;
            gds = 0.0;
            gmbs = 0.0;
        } else if vgst <= vds_eff {
            // Saturation region
            let betap = beta * (1.0 + self.lambda * vds_eff);
            cdrain = betap * vgst * vgst * 0.5;
            gm = betap * vgst;
            gds = self.lambda * beta * vgst * vgst * 0.5;
            gmbs = gm * arg;
        } else {
            // Linear region
            let betap = beta * (1.0 + self.lambda * vds_eff);
            cdrain = betap * vds_eff * (vgst - 0.5 * vds_eff);
            gm = betap * vds_eff;
            gds = betap * (vgst - vds_eff) + self.lambda * beta * vds_eff * (vgst - 0.5 * vds_eff);
            gmbs = gm * arg;
        }

        // Bulk-source and bulk-drain junction diode currents (ngspice convention).
        let vbd = vbs - vds;
        let (gbs, cbs_current) = bulk_diode_current(vbs, self.is, vt);
        let (gbd, cbd_current) = bulk_diode_current(vbd, self.is, vt);

        // Equivalent current source for NR:
        // ceq_d = cdrain - gm*vgs_eff - gds*vds_eff - gmbs*vbs_eff
        // (this is for the mode-adjusted voltages)
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
///
/// Returns (conductance, current) for a junction at voltage v.
fn bulk_diode_current(v: f64, is: f64, vt: f64) -> (f64, f64) {
    let gmin = 1e-12; // ngspice GMIN default
    if v <= -3.0 * vt {
        // Reverse bias: linear approximation
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

/// NR companion model result for a MOSFET at an operating point.
#[derive(Debug, Clone)]
pub struct MosfetCompanion {
    /// Transconductance dId/dVgs.
    pub gm: f64,
    /// Output conductance dId/dVds.
    pub gds: f64,
    /// Body effect transconductance dId/dVbs.
    pub gmbs: f64,
    /// Bulk-drain junction conductance.
    pub gbd: f64,
    /// Bulk-source junction conductance.
    pub gbs: f64,
    /// Drain current.
    pub cdrain: f64,
    /// Equivalent current source for drain (NR linearization residual).
    pub ceq_d: f64,
    /// Equivalent current source for bulk-source junction.
    pub ceq_bs: f64,
    /// Equivalent current source for bulk-drain junction.
    pub ceq_bd: f64,
    /// Operating mode: +1 normal, -1 reversed (source/drain swapped).
    pub mode: i32,
    /// Saturation voltage.
    pub vdsat: f64,
}

/// Resolved node indices for a MOSFET instance in the MNA system.
#[derive(Debug, Clone)]
pub struct MosfetInstance {
    /// MOSFET element name.
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
    /// Resolved MOSFET model parameters (with W/L scaling applied to KP).
    pub model: MosfetModel,
    /// Channel width (default 1e-4 m = 100um).
    pub w: f64,
    /// Channel length (default 1e-4 m = 100um).
    pub l: f64,
    /// Drain area for junction cap (default 0).
    pub ad: f64,
    /// Source area for junction cap (default 0).
    pub as_: f64,
    /// Drain perimeter for junction cap (default 0).
    pub pd: f64,
    /// Source perimeter for junction cap (default 0).
    pub ps: f64,
    /// Parallel multiplier (default 1.0).
    pub m: f64,
}

impl MosfetInstance {
    /// Get terminal voltages from the solution vector, handling PMOS sign.
    ///
    /// Returns (vgs, vds, vbs) with type sign already applied.
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

    /// Effective beta with W/L scaling: KP * W / L_eff.
    pub fn beta(&self) -> f64 {
        let l_eff = self.l - 2.0 * self.model.ld;
        let l_eff = l_eff.max(1e-12); // Prevent division by zero
        self.model.kp * self.w / l_eff
    }
}

/// Stamp the MOSFET companion model into the MNA matrix and RHS.
///
/// Follows the ngspice decomposition:
/// 1. gds conductance between d' and s'
/// 2. gm VCCS: Vgs controls current from s' to d'
/// 3. gmbs: Vbs controls current from s' to d'
/// 4. gbd conductance between b and d'
/// 5. gbs conductance between b and s'
/// 6. Series resistances (RD, RS)
/// 7. Equivalent current sources on the RHS
pub fn stamp_mosfet(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &MosfetInstance,
    comp: &MosfetCompanion,
) {
    let dp = inst.drain_prime_idx;
    let g = inst.gate_idx;
    let sp = inst.source_prime_idx;
    let b = inst.bulk_idx;

    let sign = inst.model.mos_type.sign();
    let m = inst.m;

    // Determine effective source/drain based on mode.
    // mode=+1: normal (dp=drain', sp=source')
    // mode=-1: reversed (dp=source', sp=drain' effectively)
    let (xnrm, xrev) = if comp.mode > 0 {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };

    // 1. gds conductance between d' and s'
    crate::stamp_conductance(matrix, dp, sp, m * comp.gds);

    // 2. gm VCCS: Vgs controls current from s' to d'
    // In normal mode: d' gets +gm from g, -gm from s'
    // s' gets -gm from g, +gm from s'
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

    // 3. gmbs: Vbs controls current from s' to d'
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

    // 4. gbd conductance between b and d'
    crate::stamp_conductance(matrix, b, dp, m * comp.gbd);

    // 5. gbs conductance between b and s'
    crate::stamp_conductance(matrix, b, sp, m * comp.gbs);

    // 6. Series resistances
    if inst.model.rd > 0.0 {
        let grd = 1.0 / inst.model.rd;
        crate::stamp_conductance(matrix, inst.drain_idx, dp, m * grd);
    }
    if inst.model.rs > 0.0 {
        let grs = 1.0 / inst.model.rs;
        crate::stamp_conductance(matrix, inst.source_idx, sp, m * grs);
    }

    // 7. Equivalent current sources on the RHS.
    // For NR: the linearized drain current model is:
    //   Id = gm*Vgs + gds*Vds + gmbs*Vbs + ceq_d
    // ceq_d = Id_prev - gm*Vgs_prev - gds*Vds_prev - gmbs*Vbs_prev
    let ceq_d = sign * m * comp.ceq_d;
    let ceq_bs = sign * m * comp.ceq_bs;
    let ceq_bd = sign * m * comp.ceq_bd;

    // Drain current enters d' and exits s' (or reversed for mode=-1,
    // but ceq_d already accounts for mode through the companion computation).
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

/// MOSFET voltage limiting for NR convergence.
///
/// Limits Vgs and Vds changes to prevent Newton-Raphson divergence.
/// Returns limited (vgs, vds).
pub fn mos_limit(vgs_new: f64, vds_new: f64, vgs_old: f64, vds_old: f64, vto: f64) -> (f64, f64) {
    let mut vgs = vgs_new;
    let mut vds = vds_new;

    // Limit Vgs change
    let vgs_diff = vgs - vgs_old;
    if vgs_diff.abs() > 0.5 {
        if vgs_diff > 0.0 {
            vgs = vgs_old + 0.5;
        } else {
            vgs = vgs_old - 0.5;
        }
    }

    // Limit Vds change near transition regions
    let vgst = vgs - vto;
    if vgst > 0.0 {
        // In active region, limit Vds more carefully
        let vds_diff = vds - vds_old;
        if vds_diff.abs() > 2.0 * vgst.abs().max(0.5) {
            if vds_diff > 0.0 {
                vds = vds_old + 2.0 * vgst.abs().max(0.5);
            } else {
                vds = vds_old - 2.0 * vgst.abs().max(0.5);
            }
        }
    }

    (vgs, vds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ferrospice_netlist::Param;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    #[test]
    fn test_default_mosfet_model() {
        let m = MosfetModel::new(MosfetType::Nmos);
        assert_eq!(m.vto, 0.0);
        assert_eq!(m.kp, 2e-5);
        assert_eq!(m.gamma, 0.0);
        assert_eq!(m.phi, 0.6);
        assert_eq!(m.lambda, 0.0);
        assert_eq!(m.rd, 0.0);
        assert_eq!(m.rs, 0.0);
    }

    #[test]
    fn test_from_model_def_nmos() {
        let model_def = ModelDef {
            name: "N1".to_string(),
            kind: "NMOS".to_string(),
            params: vec![
                Param {
                    name: "VTO".to_string(),
                    value: Expr::Num(0.7),
                },
                Param {
                    name: "KP".to_string(),
                    value: Expr::Num(1.1e-4),
                },
                Param {
                    name: "GAMMA".to_string(),
                    value: Expr::Num(0.4),
                },
                Param {
                    name: "LAMBDA".to_string(),
                    value: Expr::Num(0.04),
                },
            ],
        };
        let m = MosfetModel::from_model_def(&model_def);
        assert_eq!(m.mos_type, MosfetType::Nmos);
        assert_abs_diff_eq!(m.vto, 0.7, epsilon = 1e-15);
        assert_abs_diff_eq!(m.kp, 1.1e-4, epsilon = 1e-15);
        assert_abs_diff_eq!(m.gamma, 0.4, epsilon = 1e-15);
        assert_abs_diff_eq!(m.lambda, 0.04, epsilon = 1e-15);
    }

    #[test]
    fn test_from_model_def_pmos() {
        let model_def = ModelDef {
            name: "P1".to_string(),
            kind: "PMOS".to_string(),
            params: vec![Param {
                name: "VTO".to_string(),
                value: Expr::Num(-0.7),
            }],
        };
        let m = MosfetModel::from_model_def(&model_def);
        assert_eq!(m.mos_type, MosfetType::Pmos);
        assert_abs_diff_eq!(m.vto, -0.7, epsilon = 1e-15);
    }

    #[test]
    fn test_companion_cutoff() {
        let mut m = MosfetModel::new(MosfetType::Nmos);
        m.vto = 1.0;
        m.kp = 1e-4;
        // Vgs = 0.5 < Vto = 1.0 → cutoff
        let comp = m.companion(0.5, 5.0, 0.0);
        assert_abs_diff_eq!(comp.cdrain, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(comp.gm, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(comp.gds, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(comp.gmbs, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_companion_saturation() {
        let mut m = MosfetModel::new(MosfetType::Nmos);
        m.vto = 1.0;
        m.kp = 1e-4; // This is already KP (will be scaled by W/L via instance)
        m.lambda = 0.0;
        // Vgs = 3V, Vds = 5V, Vbs = 0
        // Vgst = 3 - 1 = 2 > 0, Vds = 5 > Vgst = 2 → saturation
        // Id = KP/2 * (Vgs-Vt)² = 1e-4/2 * 4 = 2e-4
        let comp = m.companion(3.0, 5.0, 0.0);
        assert_abs_diff_eq!(comp.cdrain, 2e-4, epsilon = 1e-10);
        assert_eq!(comp.mode, 1);
        // gm = KP * (Vgs-Vt) = 1e-4 * 2 = 2e-4
        assert_abs_diff_eq!(comp.gm, 2e-4, epsilon = 1e-10);
        // gds = 0 (lambda=0)
        assert_abs_diff_eq!(comp.gds, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_companion_linear() {
        let mut m = MosfetModel::new(MosfetType::Nmos);
        m.vto = 1.0;
        m.kp = 1e-4;
        m.lambda = 0.0;
        // Vgs = 3V, Vds = 1V, Vbs = 0
        // Vgst = 2, Vds = 1 < Vgst → linear
        // Id = KP * Vds * (Vgst - Vds/2) = 1e-4 * 1 * (2 - 0.5) = 1.5e-4
        let comp = m.companion(3.0, 1.0, 0.0);
        assert_abs_diff_eq!(comp.cdrain, 1.5e-4, epsilon = 1e-10);
        assert_eq!(comp.mode, 1);
        // gm = KP * Vds = 1e-4 * 1 = 1e-4
        assert_abs_diff_eq!(comp.gm, 1e-4, epsilon = 1e-10);
    }

    #[test]
    fn test_companion_lambda() {
        let mut m = MosfetModel::new(MosfetType::Nmos);
        m.vto = 1.0;
        m.kp = 1e-4;
        m.lambda = 0.02;
        // Saturation: Vgs=3, Vds=5, Vgst=2
        // Id = KP/2 * Vgst² * (1 + λ*Vds) = 5e-5 * 4 * 1.1 = 2.2e-4
        let comp = m.companion(3.0, 5.0, 0.0);
        assert_abs_diff_eq!(comp.cdrain, 2.2e-4, epsilon = 1e-10);
        assert!(comp.gds > 0.0, "lambda should give nonzero gds");
    }

    #[test]
    fn test_companion_reversed_mode() {
        let mut m = MosfetModel::new(MosfetType::Nmos);
        m.vto = 1.0;
        m.kp = 1e-4;
        // Negative Vds → reversed mode
        let comp = m.companion(3.0, -1.0, 0.0);
        assert_eq!(comp.mode, -1);
        // Should compute with Vgd = Vgs - Vds = 3 - (-1) = 4, Vds_eff = 1
        // Vgst = 4 - 1 = 3 > Vds_eff=1 → linear
        assert!(comp.cdrain > 0.0);
    }

    #[test]
    fn test_internal_node_count() {
        let mut m = MosfetModel::new(MosfetType::Nmos);
        assert_eq!(m.internal_node_count(), 0);
        m.rd = 10.0;
        assert_eq!(m.internal_node_count(), 1);
        m.rs = 5.0;
        assert_eq!(m.internal_node_count(), 2);
    }

    #[test]
    fn test_pmos_type_sign() {
        assert_eq!(MosfetType::Nmos.sign(), 1.0);
        assert_eq!(MosfetType::Pmos.sign(), -1.0);
    }

    #[test]
    fn test_body_effect() {
        let mut m = MosfetModel::new(MosfetType::Nmos);
        m.vto = 1.0;
        m.kp = 1e-4;
        m.gamma = 0.5;
        m.phi = 0.6;

        // With Vbs = 0 and gamma > 0, effective Vt = Vto + gamma * sqrt(phi)
        // Vt_eff = 1.0 + 0.5 * sqrt(0.6) ≈ 1.0 + 0.387 = 1.387
        let comp_zero = m.companion(3.0, 5.0, 0.0);

        // With Vbs = -2V, Vt increases more
        // Vt_eff = 1.0 + 0.5 * sqrt(0.6+2) = 1.0 + 0.5*1.612 = 1.806
        let comp_neg = m.companion(3.0, 5.0, -2.0);

        // More body bias → larger Vt → less current
        assert!(
            comp_neg.cdrain < comp_zero.cdrain,
            "negative Vbs should reduce current: {} vs {}",
            comp_neg.cdrain,
            comp_zero.cdrain
        );

        // gmbs should be nonzero with body effect
        assert!(comp_zero.gmbs > 0.0, "gmbs should be > 0 with gamma > 0");
    }
}
