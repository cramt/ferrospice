//! JFET (Junction Field-Effect Transistor) device model.
//!
//! Implements the standard SPICE Level 1 JFET model matching ngspice.
//! Supports N-channel (NJF) and P-channel (PJF) devices.
//! Gate junctions are modeled as PN diodes with leakage current.

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

/// JFET polarity: N-channel or P-channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JfetType {
    Njf,
    Pjf,
}

impl JfetType {
    /// Type multiplier: +1 for NJF, -1 for PJF.
    pub fn sign(self) -> f64 {
        match self {
            JfetType::Njf => 1.0,
            JfetType::Pjf => -1.0,
        }
    }
}

/// JFET Level 1 model parameters matching ngspice defaults.
#[derive(Debug, Clone)]
pub struct JfetModel {
    /// N-channel or P-channel.
    pub jfet_type: JfetType,
    /// Pinch-off voltage (default -2.0 V for NJF).
    pub vto: f64,
    /// Transconductance parameter (default 1e-4 A/V²).
    pub beta: f64,
    /// Channel length modulation (default 0 1/V).
    pub lambda: f64,
    /// Drain ohmic resistance (default 0 Ω).
    pub rd: f64,
    /// Source ohmic resistance (default 0 Ω).
    pub rs: f64,
    /// Gate-source zero-bias junction capacitance (default 0 F).
    pub cgs: f64,
    /// Gate-drain zero-bias junction capacitance (default 0 F).
    pub cgd: f64,
    /// Gate junction potential (default 1.0 V).
    pub pb: f64,
    /// Gate junction saturation current (default 1e-14 A).
    pub is: f64,
    /// Emission coefficient (default 1.0).
    pub n: f64,
    /// Forward bias depletion cap coefficient (default 0.5).
    pub fc: f64,
    /// Doping tail parameter (default 1.0, Sydney University mod).
    pub b: f64,
    /// Flicker noise coefficient (default 0).
    pub kf: f64,
    /// Flicker noise exponent (default 1).
    pub af: f64,
}

impl JfetModel {
    /// Create a new JfetModel with default parameters for the given type.
    pub fn new(jfet_type: JfetType) -> Self {
        Self {
            jfet_type,
            vto: -2.0,
            beta: 1e-4,
            lambda: 0.0,
            rd: 0.0,
            rs: 0.0,
            cgs: 0.0,
            cgd: 0.0,
            pb: 1.0,
            is: 1e-14,
            n: 1.0,
            fc: 0.5,
            b: 1.0,
            kf: 0.0,
            af: 1.0,
        }
    }

    /// Create a `JfetModel` from a netlist `.model` definition.
    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let jfet_type = if model_def.kind.to_uppercase().contains("PJF") {
            JfetType::Pjf
        } else {
            JfetType::Njf
        };
        let mut m = Self::new(jfet_type);
        for p in &model_def.params {
            if let Expr::Num(v) = &p.value {
                match p.name.to_uppercase().as_str() {
                    "VTO" | "VT0" => m.vto = *v,
                    "BETA" => m.beta = *v,
                    "LAMBDA" => m.lambda = *v,
                    "RD" => m.rd = *v,
                    "RS" => m.rs = *v,
                    "CGS" => m.cgs = *v,
                    "CGD" => m.cgd = *v,
                    "PB" => m.pb = *v,
                    "IS" => m.is = *v,
                    "N" => m.n = *v,
                    "FC" => m.fc = *v,
                    "B" => m.b = *v,
                    "KF" => m.kf = *v,
                    "AF" => m.af = *v,
                    _ => {} // ignore unknown params (LEVEL, etc.)
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

    /// Compute the JFET operating point and NR companion model.
    ///
    /// `vgs` and `vgd` are gate-source and gate-drain junction voltages
    /// (already adjusted for type sign). Returns `JfetCompanion` with
    /// conductances and equivalent currents.
    pub fn companion(&self, vgs: f64, vgd: f64) -> JfetCompanion {
        let vt = self.n * VT_NOM;
        let gmin = 1e-12; // ngspice GMIN default

        let vds = vgs - vgd;

        // --- Gate junction diode currents ---
        // Gate-source diode
        let (ggs, cg_gs) = gate_junction_current(vgs, vt, self.is, gmin);
        // Gate-drain diode
        let (ggd, cg_gd) = gate_junction_current(vgd, vt, self.is, gmin);

        // --- Channel current (Sydney University JFET model) ---
        // Doping tail parameter: Bfac = (1-b)/(PB-VTO)
        let b_fac = if (self.pb - self.vto).abs() > 1e-20 {
            (1.0 - self.b) / (self.pb - self.vto)
        } else {
            0.0
        };

        let (cdrain, gm, gds) = if vds >= 0.0 {
            // Normal mode (jfetload.c lines 256-288)
            let vgst = vgs - self.vto;
            if vgst <= 0.0 {
                (0.0, 0.0, 0.0)
            } else {
                let betap = self.beta * (1.0 + self.lambda * vds);
                if vgst >= vds {
                    // Linear region (jfetload.c lines 273-278)
                    let apart = 2.0 * self.b + 3.0 * b_fac * (vgst - vds);
                    let cpart = vds * (vds * (b_fac * vds - self.b) + vgst * apart);
                    let cdrain = betap * cpart;
                    let gm = betap * vds * (apart + 3.0 * b_fac * vgst);
                    let gds = betap * (vgst - vds) * apart + self.beta * self.lambda * cpart;
                    (cdrain, gm, gds)
                } else {
                    // Saturation region (jfetload.c lines 280-287)
                    let b_fac_sat = vgst * b_fac;
                    let gm = betap * vgst * (2.0 * self.b + 3.0 * b_fac_sat);
                    let cpart = vgst * vgst * (self.b + b_fac_sat);
                    let cdrain = betap * cpart;
                    let gds = self.lambda * self.beta * cpart;
                    (cdrain, gm, gds)
                }
            }
        } else {
            // Inverse mode (jfetload.c lines 291-324)
            let vgdt = vgd - self.vto;
            if vgdt <= 0.0 {
                (0.0, 0.0, 0.0)
            } else {
                let betap = self.beta * (1.0 - self.lambda * vds);
                if vgdt + vds >= 0.0 {
                    // Inverse linear (jfetload.c lines 309-314)
                    let apart = 2.0 * self.b + 3.0 * b_fac * (vgdt + vds);
                    let cpart = vds * (-vds * (-b_fac * vds - self.b) + vgdt * apart);
                    let cdrain = betap * cpart;
                    let gm = betap * vds * (apart + 3.0 * b_fac * vgdt);
                    let gds = betap * (vgdt + vds) * apart - self.beta * self.lambda * cpart - gm;
                    (cdrain, gm, gds)
                } else {
                    // Inverse saturation (jfetload.c lines 316-323)
                    let b_fac_sat = vgdt * b_fac;
                    let gm = -betap * vgdt * (2.0 * self.b + 3.0 * b_fac_sat);
                    let cpart = vgdt * vgdt * (self.b + b_fac_sat);
                    let cdrain = -betap * cpart;
                    let gds = self.lambda * self.beta * cpart - gm;
                    (cdrain, gm, gds)
                }
            }
        };

        // Total gate current = gate-source + gate-drain diode
        let _cg = cg_gs + cg_gd;
        // Drain current = channel current - gate-drain leakage
        let cd = cdrain - cg_gd;
        // Source current = -(channel current + gate-source leakage)
        // (KCL: Ig + Id + Is = 0)

        // Equivalent current sources for NR linearization:
        // ceq_gs = cg_gs - ggs * vgs
        // ceq_gd = cg_gd - ggd * vgd
        // ceq_d = cdrain - gm * vgs - gds * vds
        let ceq_gs = cg_gs - ggs * vgs;
        let ceq_gd = cg_gd - ggd * vgd;
        let ceq_d = cdrain - gm * vgs - gds * vds;

        JfetCompanion {
            gm,
            gds,
            ggs,
            ggd,
            cdrain,
            cd,
            ceq_gs,
            ceq_gd,
            ceq_d,
        }
    }

    /// Compute gate junction capacitance for AC analysis.
    ///
    /// Standard depletion capacitance with forward bias correction.
    pub fn junction_cap(&self, vj: f64, cj0: f64) -> f64 {
        if cj0 <= 0.0 {
            return 0.0;
        }
        let fc_pb = self.fc * self.pb;
        if vj < fc_pb {
            // Reverse/moderate bias: depletion model
            cj0 / (1.0 - vj / self.pb).sqrt()
        } else {
            // Forward bias: linear extrapolation
            let f1 = cj0 / (1.0 - self.fc).sqrt();
            let f2 = 0.5 / (self.pb * (1.0 - self.fc).powf(1.5));
            f1 + f2 * cj0 * (vj - fc_pb)
        }
    }
}

/// Gate junction diode current and conductance.
///
/// Returns (conductance, current) matching ngspice JFETload.c.
fn gate_junction_current(v: f64, vt: f64, is: f64, gmin: f64) -> (f64, f64) {
    if v < -3.0 * vt {
        // Deep reverse bias: avoid exp underflow
        let arg = 3.0 * vt / (v * std::f64::consts::E);
        let arg = arg * arg * arg;
        let g = is * 3.0 * arg / (-v) + gmin;
        let i = -is * (1.0 + arg) + gmin * v;
        (g, i)
    } else {
        let ev = safe_exp(v / vt);
        let g = is * ev / vt + gmin;
        let i = is * (ev - 1.0) + gmin * v;
        (g, i)
    }
}

/// NR companion model result for a JFET at an operating point.
#[derive(Debug, Clone)]
pub struct JfetCompanion {
    /// Transconductance dId/dVgs.
    pub gm: f64,
    /// Output conductance dId/dVds.
    pub gds: f64,
    /// Gate-source junction conductance.
    pub ggs: f64,
    /// Gate-drain junction conductance.
    pub ggd: f64,
    /// Channel drain current.
    pub cdrain: f64,
    /// Total drain current (channel - gate-drain leakage).
    pub cd: f64,
    /// Equivalent current source for gate-source junction.
    pub ceq_gs: f64,
    /// Equivalent current source for gate-drain junction.
    pub ceq_gd: f64,
    /// Equivalent current source for drain (channel current linearization).
    pub ceq_d: f64,
}

/// Resolved node indices for a JFET instance in the MNA system.
#[derive(Debug, Clone)]
pub struct JfetInstance {
    /// JFET element name.
    pub name: String,
    /// External drain node index (None = ground).
    pub drain_idx: Option<usize>,
    /// Gate node index (None = ground).
    pub gate_idx: Option<usize>,
    /// External source node index (None = ground).
    pub source_idx: Option<usize>,
    /// Internal drain prime node (when RD > 0), else same as drain_idx.
    pub drain_prime_idx: Option<usize>,
    /// Internal source prime node (when RS > 0), else same as source_idx.
    pub source_prime_idx: Option<usize>,
    /// Resolved JFET model parameters.
    pub model: JfetModel,
    /// Device area factor (default 1.0).
    pub area: f64,
    /// Parallel multiplier (default 1.0).
    pub m: f64,
}

impl JfetInstance {
    /// Get junction voltages from the solution vector, handling type sign.
    ///
    /// Returns (vgs, vgd) with type sign already applied.
    pub fn junction_voltages(&self, solution: &[f64]) -> (f64, f64) {
        let v_dp = self.drain_prime_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_g = self.gate_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_sp = self.source_prime_idx.map(|i| solution[i]).unwrap_or(0.0);

        let sign = self.model.jfet_type.sign();
        let vgs = sign * (v_g - v_sp);
        let vgd = sign * (v_g - v_dp);
        (vgs, vgd)
    }
}

/// Stamp the JFET companion model into the MNA matrix and RHS.
///
/// The JFET is decomposed into:
/// 1. ggs conductance g-s' (gate-source junction)
/// 2. ggd conductance g-d' (gate-drain junction)
/// 3. gds conductance d'-s' (output conductance)
/// 4. gm VCCS: Vgs controls current from s' to d'
/// 5. Series resistances (RD, RS)
/// 6. Equivalent current sources on the RHS
pub fn stamp_jfet(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &JfetInstance,
    comp: &JfetCompanion,
) {
    let dp = inst.drain_prime_idx;
    let g = inst.gate_idx;
    let sp = inst.source_prime_idx;
    let sign = inst.model.jfet_type.sign();
    let m = inst.m * inst.area;

    // 1. ggs conductance g-s' (gate-source junction)
    crate::stamp_conductance(matrix, g, sp, m * comp.ggs);

    // 2. ggd conductance g-d' (gate-drain junction)
    crate::stamp_conductance(matrix, g, dp, m * comp.ggd);

    // 3. gds conductance d'-s' (output conductance)
    crate::stamp_conductance(matrix, dp, sp, m * comp.gds);

    // 4. gm VCCS: Vgs controls current from s' to d'
    // d' gets +gm from g, -gm from s'
    // s' gets -gm from g, +gm from s'
    let gm_scaled = m * comp.gm;
    if let Some(d) = dp {
        if let Some(gate) = g {
            matrix.add(d, gate, gm_scaled);
        }
        if let Some(s) = sp {
            matrix.add(d, s, -gm_scaled);
        }
    }
    if let Some(s) = sp {
        if let Some(gate) = g {
            matrix.add(s, gate, -gm_scaled);
        }
        matrix.add(s, s, gm_scaled);
    }

    // 5. Series resistances
    if inst.model.rd > 0.0 {
        let grd = 1.0 / inst.model.rd;
        crate::stamp_conductance(matrix, inst.drain_idx, dp, m * grd);
    }
    if inst.model.rs > 0.0 {
        let grs = 1.0 / inst.model.rs;
        crate::stamp_conductance(matrix, inst.source_idx, sp, m * grs);
    }

    // 6. Equivalent current sources on the RHS.
    // Gate-source junction: ceq_gs flows from gate to source'
    // Gate-drain junction: ceq_gd flows from gate to drain'
    // Channel current linearization: ceq_d flows from source' to drain'
    let ceq_gs = sign * m * comp.ceq_gs;
    let ceq_gd = sign * m * comp.ceq_gd;
    let ceq_d = sign * m * comp.ceq_d;

    // Gate: current exits into both junctions
    if let Some(gate) = g {
        rhs[gate] -= ceq_gs + ceq_gd;
    }
    // Drain': receives gate-drain junction current + channel current
    if let Some(d) = dp {
        rhs[d] += ceq_gd - ceq_d;
    }
    // Source': receives gate-source junction current - channel current
    if let Some(s) = sp {
        rhs[s] += ceq_gs + ceq_d;
    }
}

/// JFET voltage limiting for NR convergence.
///
/// Limits Vgs and Vgd junction voltage changes using pnjlim.
pub fn jfet_limit(
    vgs_new: f64,
    vgd_new: f64,
    vgs_old: f64,
    vgd_old: f64,
    vt: f64,
    vcrit_val: f64,
) -> (f64, f64) {
    let vgs = crate::diode::pnjlim(vgs_new, vgs_old, vt, vcrit_val);
    let vgd = crate::diode::pnjlim(vgd_new, vgd_old, vt, vcrit_val);
    (vgs, vgd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ferrospice_netlist::Param;

    #[test]
    fn test_default_jfet_model() {
        let m = JfetModel::new(JfetType::Njf);
        assert_eq!(m.vto, -2.0);
        assert_eq!(m.beta, 1e-4);
        assert_eq!(m.lambda, 0.0);
        assert_eq!(m.rd, 0.0);
        assert_eq!(m.rs, 0.0);
        assert_eq!(m.cgs, 0.0);
        assert_eq!(m.cgd, 0.0);
        assert_eq!(m.pb, 1.0);
        assert_eq!(m.is, 1e-14);
    }

    #[test]
    fn test_from_model_def() {
        let model_def = ModelDef {
            name: "MODJ".to_string(),
            kind: "NJF".to_string(),
            params: vec![
                Param {
                    name: "VTO".to_string(),
                    value: Expr::Num(-3.5),
                },
                Param {
                    name: "BETA".to_string(),
                    value: Expr::Num(4.1e-4),
                },
                Param {
                    name: "LAMBDA".to_string(),
                    value: Expr::Num(0.002),
                },
                Param {
                    name: "RD".to_string(),
                    value: Expr::Num(200.0),
                },
            ],
        };
        let m = JfetModel::from_model_def(&model_def);
        assert_eq!(m.jfet_type, JfetType::Njf);
        assert_abs_diff_eq!(m.vto, -3.5, epsilon = 1e-15);
        assert_abs_diff_eq!(m.beta, 4.1e-4, epsilon = 1e-15);
        assert_abs_diff_eq!(m.lambda, 0.002, epsilon = 1e-15);
        assert_abs_diff_eq!(m.rd, 200.0, epsilon = 1e-15);
    }

    #[test]
    fn test_from_model_def_pjf() {
        let model_def = ModelDef {
            name: "PMOD".to_string(),
            kind: "PJF".to_string(),
            params: vec![Param {
                name: "VTO".to_string(),
                value: Expr::Num(2.0),
            }],
        };
        let m = JfetModel::from_model_def(&model_def);
        assert_eq!(m.jfet_type, JfetType::Pjf);
        assert_abs_diff_eq!(m.vto, 2.0, epsilon = 1e-15);
    }

    #[test]
    fn test_companion_cutoff() {
        let mut m = JfetModel::new(JfetType::Njf);
        m.vto = -2.0;
        m.beta = 1e-4;
        // VGS = -3V, VTO = -2V → Vgst = VGS - VTO = -3 - (-2) = -1 ≤ 0 → cutoff
        let comp = m.companion(-3.0, -8.0);
        assert_abs_diff_eq!(comp.cdrain, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(comp.gm, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(comp.gds, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_companion_saturation() {
        let mut m = JfetModel::new(JfetType::Njf);
        m.vto = -2.0;
        m.beta = 1e-4;
        m.lambda = 0.0;
        m.b = 1.0; // b=1 simplifies to standard quadratic model
        // VGS = 0V, VTO = -2V → Vgst = 0 - (-2) = 2 > 0
        // VDS = VGS - VGD = 0 - (-5) = 5 > Vgst → saturation
        // With b=1: b_fac2 = 0, cdrain = beta * Vgst² * 1 = 1e-4 * 4 = 4e-4
        let comp = m.companion(0.0, -5.0);
        assert_abs_diff_eq!(comp.cdrain, 4e-4, epsilon = 1e-10);
        // gm = beta * Vgst * 2 = 1e-4 * 2 * 2 = 4e-4
        assert_abs_diff_eq!(comp.gm, 4e-4, epsilon = 1e-10);
        assert_abs_diff_eq!(comp.gds, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_companion_lambda() {
        let mut m = JfetModel::new(JfetType::Njf);
        m.vto = -2.0;
        m.beta = 1e-4;
        m.lambda = 0.01;
        m.b = 1.0;
        // Saturation: VGS=0, VGD=-5, VDS=5, Vgst=2
        // cdrain = beta * (1 + lambda*vds) * Vgst² = 1e-4 * 1.05 * 4 = 4.2e-4
        let comp = m.companion(0.0, -5.0);
        assert_abs_diff_eq!(comp.cdrain, 4.2e-4, epsilon = 1e-10);
        assert!(comp.gds > 0.0, "lambda should give nonzero gds");
    }

    #[test]
    fn test_companion_2n4221_op() {
        // Match the ngspice test: 2N4221 JFET at VGS=-2, VDS=25
        // VTO=-3.5, BETA=4.1e-4, LAMBDA=0.002, RD=200
        let mut m = JfetModel::new(JfetType::Njf);
        m.vto = -3.5;
        m.beta = 4.1e-4;
        m.lambda = 0.002;
        m.b = 1.0;

        // VGS=-2, VGD=VGS-VDS=-2-25=-27
        let comp = m.companion(-2.0, -27.0);

        // Expected from ngspice output: id ≈ 0.000968268
        // Vgst = -2 - (-3.5) = 1.5
        // VDS = 25 > Vgst = 1.5 → saturation
        // cdrain = beta * (1 + lambda*vds) * Vgst²
        //        = 4.1e-4 * (1 + 0.002*25) * 1.5² = 4.1e-4 * 1.05 * 2.25
        //        = 9.6862...e-4
        assert_abs_diff_eq!(comp.cdrain, 9.68625e-4, epsilon = 1e-7);

        // gm from ngspice: 0.00129102
        assert_abs_diff_eq!(comp.gm, 1.29150e-3, epsilon = 1e-5);

        // gds from ngspice: 1.845e-06
        assert_abs_diff_eq!(comp.gds, 1.845e-6, epsilon = 1e-8);
    }

    #[test]
    fn test_internal_node_count() {
        let mut m = JfetModel::new(JfetType::Njf);
        assert_eq!(m.internal_node_count(), 0);
        m.rd = 10.0;
        assert_eq!(m.internal_node_count(), 1);
        m.rs = 5.0;
        assert_eq!(m.internal_node_count(), 2);
    }

    #[test]
    fn test_type_sign() {
        assert_eq!(JfetType::Njf.sign(), 1.0);
        assert_eq!(JfetType::Pjf.sign(), -1.0);
    }

    #[test]
    fn test_junction_cap() {
        let m = JfetModel::new(JfetType::Njf);
        // At V=0: C = CJ0
        let m2 = JfetModel { cgs: 5e-12, ..m };
        assert_abs_diff_eq!(m2.junction_cap(0.0, m2.cgs), 5e-12, epsilon = 1e-15);
        // Reverse bias: C < CJ0
        let c_rev = m2.junction_cap(-1.0, m2.cgs);
        assert!(c_rev < 5e-12, "reverse bias cap should be less");
    }
}
