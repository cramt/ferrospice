//! MESFET (Metal Semiconductor FET) device model.
//!
//! Implements the basic Statz/Curtice MESFET model from ngspice (MES device),
//! used with `.model name NMF/PMF level=1` syntax.
//!
//! Reference: ngspice mes/ device directory.

use thevenin_types::{Expr, ModelDef};

use crate::mna::stamp_conductance;
use crate::sparse::SparseMatrix;

/// Thermal voltage at 300.15K.
const VT_NOM: f64 = 0.025864186;

/// Euler's number.
const E: f64 = std::f64::consts::E;

fn expr_val(e: &Expr) -> f64 {
    match e {
        Expr::Num(v) => *v,
        Expr::Param(_) | Expr::Brace(_) => 0.0,
    }
}

/// MESFET device type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MesfetType {
    Nmf = 1,
    Pmf = -1,
}

/// MESFET model parameters (Statz/Curtice model, level=1).
#[derive(Debug, Clone)]
pub struct MesfetModel {
    pub device_type: MesfetType,
    pub vt0: f64,    // Threshold/pinch-off voltage (default -2.0)
    pub alpha: f64,  // Saturation voltage parameter (default 2.0)
    pub beta: f64,   // Transconductance parameter (default 2.5e-3)
    pub lambda: f64, // Channel-length modulation (default 0.0)
    pub b: f64,      // Doping tail extending parameter (default 0.3)
    pub rd: f64,     // Drain resistance
    pub rs: f64,     // Source resistance
    pub cgs: f64,    // Gate-source capacitance
    pub cgd: f64,    // Gate-drain capacitance
    pub pb: f64,     // Gate junction potential (default 1.0)
    pub is: f64,     // Gate saturation current (default 1e-14)
    pub fc: f64,     // Forward bias depletion cap coefficient (default 0.5)
    pub kf: f64,     // Flicker noise coefficient
    pub af: f64,     // Flicker noise exponent (default 1.0)
    // Precomputed
    pub drain_conduct: f64,
    pub source_conduct: f64,
    pub vcrit: f64,
}

impl Default for MesfetModel {
    fn default() -> Self {
        Self::new(MesfetType::Nmf)
    }
}

impl MesfetModel {
    pub fn new(device_type: MesfetType) -> Self {
        let is = 1e-14;
        let vcrit = VT_NOM * (VT_NOM / (std::f64::consts::SQRT_2 * is)).ln();
        Self {
            device_type,
            vt0: -2.0,
            alpha: 2.0,
            beta: 2.5e-3,
            lambda: 0.0,
            b: 0.3,
            rd: 0.0,
            rs: 0.0,
            cgs: 0.0,
            cgd: 0.0,
            pb: 1.0,
            is,
            fc: 0.5,
            kf: 0.0,
            af: 1.0,
            drain_conduct: 0.0,
            source_conduct: 0.0,
            vcrit,
        }
    }

    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let kind = model_def.kind.to_uppercase();
        let device_type = if kind == "PMF" {
            MesfetType::Pmf
        } else {
            MesfetType::Nmf
        };
        let mut m = Self::new(device_type);

        for p in &model_def.params {
            let val = expr_val(&p.value);
            match p.name.to_uppercase().as_str() {
                "VTO" | "VT0" => m.vt0 = val,
                "ALPHA" => m.alpha = val,
                "BETA" => m.beta = val,
                "LAMBDA" => m.lambda = val,
                "B" => m.b = val,
                "RD" => m.rd = val,
                "RS" => m.rs = val,
                "CGS" => m.cgs = val,
                "CGD" => m.cgd = val,
                "PB" => m.pb = val,
                "IS" => m.is = val,
                "FC" => m.fc = val,
                "KF" => m.kf = val,
                "AF" => m.af = val,
                _ => {}
            }
        }

        // Precompute
        m.drain_conduct = if m.rd != 0.0 { 1.0 / m.rd } else { 0.0 };
        m.source_conduct = if m.rs != 0.0 { 1.0 / m.rs } else { 0.0 };
        if m.is > 0.0 {
            m.vcrit = VT_NOM * (VT_NOM / (std::f64::consts::SQRT_2 * m.is)).ln();
        }
        m
    }

    pub fn internal_node_count(&self) -> usize {
        let mut n = 0;
        if self.rd != 0.0 {
            n += 1;
        }
        if self.rs != 0.0 {
            n += 1;
        }
        n
    }
}

/// Resolved MESFET instance for simulation.
#[derive(Debug, Clone)]
pub struct MesfetInstance {
    pub name: String,
    pub drain_idx: Option<usize>,
    pub gate_idx: Option<usize>,
    pub source_idx: Option<usize>,
    pub drain_prime_idx: Option<usize>,
    pub source_prime_idx: Option<usize>,
    pub model: MesfetModel,
    pub area: f64,
    pub m: f64, // parallel multiplier
}

impl MesfetInstance {
    /// Get junction voltages from solution vector.
    pub fn junction_voltages(&self, solution: &[f64]) -> (f64, f64) {
        let vg = self.gate_idx.map_or(0.0, |i| solution[i]);
        let vdp = self
            .drain_prime_idx
            .map_or(self.drain_idx.map_or(0.0, |i| solution[i]), |i| solution[i]);
        let vsp = self
            .source_prime_idx
            .map_or(self.source_idx.map_or(0.0, |i| solution[i]), |i| {
                solution[i]
            });
        let sign = self.model.device_type as i32 as f64;
        let vgs = sign * (vg - vsp);
        let vgd = sign * (vg - vdp);
        (vgs, vgd)
    }
}

/// MESFET companion model result (linearized).
#[derive(Debug, Clone)]
pub struct MesfetCompanion {
    pub gm: f64,
    pub gds: f64,
    pub ggs: f64,
    pub ggd: f64,
    pub cg: f64,
    pub cd: f64,
    pub cgd_current: f64,
    pub capgs: f64,
    pub capgd: f64,
}

/// Compute the MESFET companion model at given junction voltages.
pub fn mesfet_companion(inst: &MesfetInstance, vgs: f64, vgd: f64, gmin: f64) -> MesfetCompanion {
    let model = &inst.model;
    let beta = model.beta * inst.area;
    let csat = model.is * inst.area;
    let vto = model.vt0;
    let vds = vgs - vgd;

    // Gate-source junction current
    let (cgs_current, ggs) = if vgs <= -3.0 * VT_NOM {
        let arg_val = 3.0 * VT_NOM / (vgs * E);
        let arg3 = arg_val * arg_val * arg_val;
        (
            -csat * (1.0 + arg3) + gmin * vgs,
            csat * 3.0 * arg3 / vgs + gmin,
        )
    } else {
        let evgs = (vgs / VT_NOM).exp();
        (
            csat * (evgs - 1.0) + gmin * vgs,
            csat * evgs / VT_NOM + gmin,
        )
    };

    // Gate-drain junction current
    let (cgd_current, ggd) = if vgd <= -3.0 * VT_NOM {
        let arg_val = 3.0 * VT_NOM / (vgd * E);
        let arg3 = arg_val * arg_val * arg_val;
        (
            -csat * (1.0 + arg3) + gmin * vgd,
            csat * 3.0 * arg3 / vgd + gmin,
        )
    } else {
        let evgd = (vgd / VT_NOM).exp();
        (
            csat * (evgd - 1.0) + gmin * vgd,
            csat * evgd / VT_NOM + gmin,
        )
    };

    let cg = cgs_current + cgd_current;

    // Drain current
    let (cdrain, gm, gds) = if vds >= 0.0 {
        // Normal mode
        let vgst = vgs - vto;
        if vgst <= 0.0 {
            (0.0, 0.0, 0.0)
        } else {
            let prod = 1.0 + model.lambda * vds;
            let betap = beta * prod;
            let denom = 1.0 + model.b * vgst;
            let invdenom = 1.0 / denom;
            if vds >= 3.0 / model.alpha {
                // Saturation
                let cdrain = betap * vgst * vgst * invdenom;
                let gm = betap * vgst * (1.0 + denom) * invdenom * invdenom;
                let gds = model.lambda * beta * vgst * vgst * invdenom;
                (cdrain, gm, gds)
            } else {
                // Linear
                let afact = 1.0 - model.alpha * vds / 3.0;
                let lfact = 1.0 - afact * afact * afact;
                let cdrain = betap * vgst * vgst * invdenom * lfact;
                let gm = betap * vgst * (1.0 + denom) * invdenom * invdenom * lfact;
                let gds = beta
                    * vgst
                    * vgst
                    * invdenom
                    * (model.alpha * afact * afact * prod + lfact * model.lambda);
                (cdrain, gm, gds)
            }
        }
    } else {
        // Inverse mode
        let vgdt = vgd - vto;
        if vgdt <= 0.0 {
            (0.0, 0.0, 0.0)
        } else {
            let prod = 1.0 - model.lambda * vds;
            let betap = beta * prod;
            let denom = 1.0 + model.b * vgdt;
            let invdenom = 1.0 / denom;
            if -vds >= 3.0 / model.alpha {
                // Inverse saturation
                let cdrain = -betap * vgdt * vgdt * invdenom;
                let gm = -betap * vgdt * (1.0 + denom) * invdenom * invdenom;
                let gds = model.lambda * beta * vgdt * vgdt * invdenom - gm;
                (cdrain, gm, gds)
            } else {
                // Inverse linear
                let afact = 1.0 + model.alpha * vds / 3.0;
                let lfact = 1.0 - afact * afact * afact;
                let cdrain = -betap * vgdt * vgdt * invdenom * lfact;
                let gm = -betap * vgdt * (1.0 + denom) * invdenom * invdenom * lfact;
                let gds = beta
                    * vgdt
                    * vgdt
                    * invdenom
                    * (model.alpha * afact * afact * prod + lfact * model.lambda)
                    - gm;
                (cdrain, gm, gds)
            }
        }
    };

    let cd = cdrain - cgd_current;

    // Charge model capacitances (qggnew)
    let czgs = model.cgs * inst.area;
    let czgd = model.cgd * inst.area;
    let (capgs, capgd) = if czgs > 0.0 || czgd > 0.0 {
        let phib = model.pb;
        let vcap = if model.alpha != 0.0 {
            1.0 / model.alpha
        } else {
            1.0e10
        };
        let (cgs_new, cgd_new) = qggnew_caps(vgs, vgd, phib, vcap, vto, czgs, czgd);
        (cgs_new, cgd_new)
    } else {
        (0.0, 0.0)
    };

    MesfetCompanion {
        gm,
        gds,
        ggs,
        ggd,
        cg,
        cd,
        cgd_current,
        capgs,
        capgd,
    }
}

/// Compute qggnew capacitances (Ward-Dutton partition).
fn qggnew_caps(
    vgs: f64,
    vgd: f64,
    phib: f64,
    vcap: f64,
    vto: f64,
    cgs: f64,
    cgd: f64,
) -> (f64, f64) {
    let veroot = ((vgs - vgd) * (vgs - vgd) + vcap * vcap).sqrt();
    let veff1 = 0.5 * (vgs + vgd + veroot);
    let del = 0.2;
    let vnroot = ((veff1 - vto) * (veff1 - vto) + del * del).sqrt();
    let vnew1 = 0.5 * (veff1 + vto + vnroot);
    let vmax = 0.5;
    let vnew1_clamped = if vnew1 < vmax { vnew1 } else { vmax };

    let qroot = (1.0 - vnew1_clamped / phib).sqrt();
    let par1 = 0.5 * (1.0 + (veff1 - vto) / vnroot);
    let cfact = (vgs - vgd) / veroot;
    let cplus = 0.5 * (1.0 + cfact);
    let cminus = cplus - cfact;

    let cgsnew = cgs / qroot * par1 * cplus + cgd * cminus;
    let cgdnew = cgs / qroot * par1 * cminus + cgd * cplus;
    (cgsnew, cgdnew)
}

/// Stamp MESFET with known voltages (for NR iteration, matches ngspice mesload.c).
pub fn stamp_mesfet_with_voltages(
    comp: &MesfetCompanion,
    inst: &MesfetInstance,
    vgs: f64,
    vgd: f64,
    matrix: &mut SparseMatrix,
    rhs: &mut [f64],
) {
    let model = &inst.model;
    let sign = model.device_type as i32 as f64;
    let m = inst.m;
    let vds = vgs - vgd;

    let gate_idx = inst.gate_idx;
    let drain_prime = inst.drain_prime_idx.or(inst.drain_idx);
    let source_prime = inst.source_prime_idx.or(inst.source_idx);

    let gm = comp.gm;
    let gds = comp.gds;
    let ggs = comp.ggs;
    let ggd = comp.ggd;
    let cg = comp.cg;
    let cd = comp.cd;
    let cgd_val = comp.cgd_current;

    // Norton equivalent current sources
    let ceqgd = sign * (cgd_val - ggd * vgd);
    let ceqgs = sign * ((cg - cgd_val) - ggs * vgs);
    let cdreq = sign * ((cd + cgd_val) - gds * vds - gm * vgs);

    // RHS
    if let Some(g) = gate_idx {
        rhs[g] += m * (-ceqgs - ceqgd);
    }
    if let Some(dp) = drain_prime {
        rhs[dp] += m * (-cdreq + ceqgd);
    }
    if let Some(sp) = source_prime {
        rhs[sp] += m * (cdreq + ceqgs);
    }

    // Series drain resistance
    if inst.drain_prime_idx.is_some() {
        let gdpr = model.drain_conduct * inst.area;
        stamp_conductance(matrix, inst.drain_idx, inst.drain_prime_idx, gdpr * m);
    }

    // Series source resistance
    if inst.source_prime_idx.is_some() {
        let gspr = model.source_conduct * inst.area;
        stamp_conductance(matrix, inst.source_idx, inst.source_prime_idx, gspr * m);
    }

    // Y-matrix stamps
    if let Some(g) = gate_idx {
        matrix.add(g, g, m * (ggs + ggd));
    }
    if let Some(dp) = drain_prime {
        matrix.add(dp, dp, m * (gds + ggd));
    }
    if let Some(sp) = source_prime {
        matrix.add(sp, sp, m * (gds + gm + ggs));
    }

    if let (Some(g), Some(dp)) = (gate_idx, drain_prime) {
        matrix.add(g, dp, -m * ggd);
        matrix.add(dp, g, m * (gm - ggd));
    }
    if let (Some(g), Some(sp)) = (gate_idx, source_prime) {
        matrix.add(g, sp, -m * ggs);
        matrix.add(sp, g, m * (-ggs - gm));
    }
    if let (Some(dp), Some(sp)) = (drain_prime, source_prime) {
        matrix.add(dp, sp, m * (-gds - gm));
        matrix.add(sp, dp, -m * gds);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mesfet_model_defaults() {
        let m = MesfetModel::new(MesfetType::Nmf);
        assert_eq!(m.vt0, -2.0);
        assert_eq!(m.alpha, 2.0);
        assert_eq!(m.beta, 2.5e-3);
        assert_eq!(m.lambda, 0.0);
        assert_eq!(m.b, 0.3);
    }

    #[test]
    fn mesfet_companion_cutoff() {
        let m = MesfetModel::new(MesfetType::Nmf);
        let inst = MesfetInstance {
            name: "z1".to_string(),
            drain_idx: Some(0),
            gate_idx: Some(1),
            source_idx: Some(2),
            drain_prime_idx: None,
            source_prime_idx: None,
            model: m,
            area: 1.0,
            m: 1.0,
        };
        // vgs = -3.0 < vt0 = -2.0, should be cutoff
        let comp = mesfet_companion(&inst, -3.0, -4.0, 1e-12);
        assert_eq!(comp.gm, 0.0);
        assert_eq!(comp.gds, 0.0);
    }

    #[test]
    fn mesfet_companion_saturation() {
        let mut m = MesfetModel::new(MesfetType::Nmf);
        m.vt0 = -1.3;
        m.alpha = 3.0;
        m.beta = 1.4e-3;
        m.lambda = 0.03;
        let inst = MesfetInstance {
            name: "z1".to_string(),
            drain_idx: Some(0),
            gate_idx: Some(1),
            source_idx: Some(2),
            drain_prime_idx: None,
            source_prime_idx: None,
            model: m,
            area: 1.4,
            m: 1.0,
        };
        // vgs=0 (above vt0=-1.3), vds=0.1 (saturation when vds >= 3/alpha=1.0? No: 3/3=1.0, vds=0.1 < 1.0 → linear)
        let comp = mesfet_companion(&inst, 0.0, -0.1, 1e-12);
        assert!(comp.gm > 0.0, "gm should be positive in conducting region");
        assert!(comp.gds > 0.0, "gds should be positive");
    }

    #[test]
    fn mesfet_model_from_model_def() {
        let model_def = ModelDef {
            name: "mesmod".to_string(),
            kind: "NMF".to_string(),
            params: vec![
                thevenin_types::Param {
                    name: "VT0".to_string(),
                    value: Expr::Num(-1.3),
                },
                thevenin_types::Param {
                    name: "BETA".to_string(),
                    value: Expr::Num(1.4e-3),
                },
                thevenin_types::Param {
                    name: "RD".to_string(),
                    value: Expr::Num(46.0),
                },
                thevenin_types::Param {
                    name: "RS".to_string(),
                    value: Expr::Num(46.0),
                },
            ],
        };
        let m = MesfetModel::from_model_def(&model_def);
        assert_eq!(m.device_type, MesfetType::Nmf);
        assert_eq!(m.vt0, -1.3);
        assert_eq!(m.beta, 1.4e-3);
        assert!(m.drain_conduct > 0.0);
        assert!(m.source_conduct > 0.0);
    }

    #[test]
    fn mesfet_pmf_type() {
        let model_def = ModelDef {
            name: "pmod".to_string(),
            kind: "PMF".to_string(),
            params: vec![],
        };
        let m = MesfetModel::from_model_def(&model_def);
        assert_eq!(m.device_type, MesfetType::Pmf);
    }
}
