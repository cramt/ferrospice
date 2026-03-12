//! MESA (Metal Semiconductor) FET device model.
//!
//! Implements the Ytterdal/Lee/Shur/Fjeldly GaAs MESFET model from ngspice,
//! supporting three levels:
//!   - Level 2 (mesa1): basic model with mu mobility modulation
//!   - Level 3 (mesa2): dual-doped heterostructure model
//!   - Level 4 (mesa3): charge-based model with NMAX saturation
//!
//! Reference: ngspice mesa/ device directory.

use thevenin_types::{Expr, ModelDef};

/// Physical constants matching ngspice (const.h).
const CHARGE: f64 = 1.602_176_620_8e-19;
const BOLTZMANN: f64 = 1.380_648_52e-23;
pub const K_OVER_Q: f64 = BOLTZMANN / CHARGE;
const C_TO_K: f64 = 273.15;
const EPSI_GAAS: f64 = 12.244 * 8.85418e-12;
const ROOT2: f64 = std::f64::consts::SQRT_2;

fn expr_val(e: &Expr) -> f64 {
    match e {
        Expr::Num(v) => *v,
        Expr::Param(_) | Expr::Brace(_) => 0.0,
    }
}

/// MESA FET model parameters.
#[derive(Debug, Clone)]
pub struct MesaModel {
    pub level: i32,
    pub vto: f64,
    pub lambda: f64,
    pub lambdahf: f64,
    pub beta: f64,
    pub vs: f64,
    pub eta: f64,
    pub m: f64,
    pub mc: f64,
    pub alpha: f64,
    pub sigma0: f64,
    pub vsigmat: f64,
    pub vsigma: f64,
    pub mu: f64,
    pub theta: f64,
    pub mu1: f64,
    pub mu2: f64,
    pub d: f64,
    pub nd: f64,
    pub du: f64,
    pub ndu: f64,
    pub th: f64,
    pub ndelta: f64,
    pub delta: f64,
    pub tc: f64,
    pub n: f64,
    pub rd: f64,
    pub rs: f64,
    pub rg: f64,
    pub ri: f64,
    pub rf: f64,
    pub rdi: f64,
    pub rsi: f64,
    pub phib: f64,
    pub phib1: f64,
    pub astar: f64,
    pub ggr: f64,
    pub del: f64,
    pub xchi: f64,
    pub tvto: f64,
    pub tlambda: f64,
    pub teta0: f64,
    pub teta1: f64,
    pub tmu: f64,
    pub xtm0: f64,
    pub xtm1: f64,
    pub xtm2: f64,
    pub ks: f64,
    pub vsg: f64,
    pub tf: f64,
    pub flo: f64,
    pub delfo: f64,
    pub ag: f64,
    pub tc1: f64,
    pub tc2: f64,
    pub zeta: f64,
    pub nmax: f64,
    pub gamma: f64,
    pub epsi: f64,
    pub cas: f64,
    pub cbs: f64,
    pub kf: f64,
    pub af: f64,
    // Precomputed
    pub vpo: f64,
    pub vpou: f64,
    pub vpod: f64,
    pub delta_sqr: f64,
}

impl Default for MesaModel {
    fn default() -> Self {
        Self::new()
    }
}

impl MesaModel {
    pub fn new() -> Self {
        let mut m = Self {
            level: 2,
            vto: -1.26,
            lambda: 0.045,
            lambdahf: 0.045,
            beta: 0.0085,
            vs: 1.5e5,
            eta: 1.73,
            m: 2.5,
            mc: 3.0,
            alpha: 0.0,
            sigma0: 0.081,
            vsigmat: 1.01,
            vsigma: 0.1,
            mu: 0.23,
            theta: 0.0,
            mu1: 0.0,
            mu2: 0.0,
            d: 0.12e-6,
            nd: 2.0e23,
            du: 0.035e-6,
            ndu: 1.0e22,
            th: 0.01e-6,
            ndelta: 6.0e24,
            delta: 5.0,
            tc: 0.0,
            n: 1.0,
            rd: 0.0,
            rs: 0.0,
            rg: 0.0,
            ri: 0.0,
            rf: 0.0,
            rdi: 0.0,
            rsi: 0.0,
            phib: 0.5 * CHARGE,
            phib1: 0.0,
            astar: 4.0e4,
            ggr: 40.0,
            del: 0.04,
            xchi: 0.033,
            tvto: 0.0,
            tlambda: f64::MAX,
            teta0: f64::MAX,
            teta1: 0.0,
            tmu: 300.15,
            xtm0: 0.0,
            xtm1: 0.0,
            xtm2: 0.0,
            ks: 0.0,
            vsg: 0.0,
            tf: 300.15,
            flo: 0.0,
            delfo: 0.0,
            ag: 0.0,
            tc1: 0.0,
            tc2: 0.0,
            zeta: 1.0,
            nmax: 2.0e16,
            gamma: 3.0,
            epsi: EPSI_GAAS,
            cas: 1.0,
            cbs: 1.0,
            kf: 0.0,
            af: 1.0,
            vpo: 0.0,
            vpou: 0.0,
            vpod: 0.0,
            delta_sqr: 25.0,
        };
        m.precompute();
        m
    }

    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let mut m = Self::new();
        let mut lambdahf_given = false;
        let mut tf_given = false;

        for p in &model_def.params {
            let val = expr_val(&p.value);
            match p.name.to_uppercase().as_str() {
                "VTO" => m.vto = val,
                "LAMBDA" => m.lambda = val,
                "LAMBDAHF" => {
                    m.lambdahf = val;
                    lambdahf_given = true;
                }
                "BETA" => m.beta = val,
                "VS" => m.vs = val,
                "ETA" => m.eta = val,
                "M" => m.m = val,
                "MC" => m.mc = val,
                "ALPHA" => m.alpha = val,
                "SIGMA0" => m.sigma0 = val,
                "VSIGMAT" => m.vsigmat = val,
                "VSIGMA" => m.vsigma = val,
                "MU" => m.mu = val,
                "THETA" => m.theta = val,
                "MU1" => m.mu1 = val,
                "MU2" => m.mu2 = val,
                "D" => m.d = val,
                "ND" => m.nd = val,
                "DU" => m.du = val,
                "NDU" => m.ndu = val,
                "TH" => m.th = val,
                "NDELTA" => m.ndelta = val,
                "DELTA" => m.delta = val,
                "TC" => m.tc = val,
                "N" => m.n = val,
                "RD" => m.rd = val,
                "RS" => m.rs = val,
                "RG" => m.rg = val,
                "RI" => m.ri = val,
                "RF" => m.rf = val,
                "RDI" => m.rdi = val,
                "RSI" => m.rsi = val,
                "PHIB" => m.phib = val * CHARGE,
                "PHIB1" => m.phib1 = val * CHARGE,
                "ASTAR" => m.astar = val,
                "GGR" => m.ggr = val,
                "DEL" => m.del = val,
                "XCHI" => m.xchi = val,
                "TVTO" => m.tvto = val,
                "TLAMBDA" => m.tlambda = val + C_TO_K,
                "TETA0" => m.teta0 = val + C_TO_K,
                "TETA1" => m.teta1 = val + C_TO_K,
                "TMU" => m.tmu = val + C_TO_K,
                "XTM0" => m.xtm0 = val,
                "XTM1" => m.xtm1 = val,
                "XTM2" => m.xtm2 = val,
                "KS" => m.ks = val,
                "VSG" => m.vsg = val,
                "TF" => {
                    m.tf = val + C_TO_K;
                    tf_given = true;
                }
                "FLO" => m.flo = val,
                "DELFO" => m.delfo = val,
                "AG" => m.ag = val,
                "TC1" | "RTC1" => m.tc1 = val,
                "TC2" | "RTC2" => m.tc2 = val,
                "ZETA" => m.zeta = val,
                "LEVEL" => m.level = val as i32,
                "NMAX" => m.nmax = val,
                "GAMMA" => m.gamma = val,
                "EPSI" => m.epsi = val,
                "CAS" => m.cas = val,
                "CBS" => m.cbs = val,
                "KF" => m.kf = val,
                "AF" => m.af = val,
                _ => {}
            }
        }
        if !lambdahf_given {
            m.lambdahf = m.lambda;
        }
        if !tf_given {
            m.tf = 300.15;
        }
        m.precompute();
        m
    }

    pub fn precompute(&mut self) {
        if self.level == 2 {
            self.vpo = CHARGE * self.nd * self.d * self.d / (2.0 * EPSI_GAAS);
        } else {
            self.vpou = CHARGE * self.ndu * self.du * self.du / (2.0 * EPSI_GAAS);
            self.vpod =
                CHARGE * self.ndelta * self.th * (2.0 * self.du + self.th) / (2.0 * EPSI_GAAS);
            self.vpo = self.vpou + self.vpod;
        }
        self.delta_sqr = self.delta * self.delta;
    }

    pub fn internal_node_count(&self) -> usize {
        let mut count = 0;
        if self.rd != 0.0 {
            count += 1;
        }
        if self.rs != 0.0 {
            count += 1;
        }
        if self.rg != 0.0 {
            count += 1;
        }
        if self.ri != 0.0 {
            count += 1;
        }
        if self.rf != 0.0 {
            count += 1;
        }
        count
    }
}

/// Precomputed per-instance temperature-dependent parameters.
#[derive(Debug, Clone)]
pub struct MesaPrecomp {
    pub t_mu: f64,
    pub t_vto: f64,
    pub t_eta: f64,
    pub t_lambda: f64,
    pub t_lambdahf: f64,
    pub t_rsi: f64,
    pub t_rdi: f64,
    pub drain_conduct: f64,
    pub source_conduct: f64,
    pub gate_conduct: f64,
    pub t_gi: f64,
    pub t_gf: f64,
    pub gchi0: f64,
    pub beta_inst: f64,
    pub n0: f64,
    pub nsb0: f64,
    pub isatb0: f64,
    pub cf: f64,
    pub csatfs: f64,
    pub csatfd: f64,
    pub ggrwl: f64,
    pub vcrits: f64,
    pub vcritd: f64,
    pub imax: f64,
    pub fl: f64,
    pub delf: f64,
    pub ts: f64,
    pub td: f64,
}

impl MesaPrecomp {
    pub fn compute(model: &MesaModel, ts: f64, td: f64, tnom: f64, w: f64, l: f64) -> Self {
        let vt = K_OVER_Q * ts;

        let t_mu = if model.mu1 == 0.0 && model.mu2 == 0.0 {
            model.mu * (ts / model.tmu).powf(model.xtm0)
        } else {
            let muimp = model.mu * (ts / model.tmu).powf(model.xtm0);
            let mupo = model.mu1 * (model.tmu / ts).powf(model.xtm1)
                + model.mu2 * (model.tmu / ts).powf(model.xtm2);
            1.0 / (1.0 / muimp + 1.0 / mupo)
        };

        let t_vto = model.vto - model.tvto * (ts - tnom);

        let gchi0 = if model.level == 2 {
            CHARGE * w / l
        } else {
            CHARGE * w / l * t_mu
        };

        let beta_inst = 2.0 * EPSI_GAAS * model.vs * model.zeta * w / model.d;

        let t_eta = model.eta * (1.0 + ts / model.teta0) + model.teta1 / ts;
        let t_lambda = model.lambda * (1.0 - ts / model.tlambda);
        let t_lambdahf = model.lambdahf * (1.0 - ts / model.tlambda);

        let d_for_n0 = if model.level == 3 { model.du } else { model.d };
        let n0 = if model.level == 4 {
            model.epsi * t_eta * vt / (2.0 * CHARGE * d_for_n0)
        } else {
            EPSI_GAAS * t_eta * vt / (CHARGE * d_for_n0)
        };

        let nsb0 = EPSI_GAAS * t_eta * vt / (CHARGE * (model.du + model.th));
        let isatb0 = CHARGE * n0 * vt * w / l;
        let cf = if model.level == 4 {
            0.5 * model.epsi * w
        } else {
            0.5 * EPSI_GAAS * w
        };

        // Temperature-adjust phib
        let phib_t = model.phib - model.phib1 * (ts - tnom);
        let csatfs = 0.5 * model.astar * ts * ts * (-phib_t / (BOLTZMANN * ts)).exp() * l * w;
        let csatfd = 0.5 * model.astar * td * td * (-phib_t / (BOLTZMANN * td)).exp() * l * w;

        let ggrwl = model.ggr * l * w * (model.xchi * (ts - tnom)).exp();

        let vcrits = if csatfs != 0.0 {
            vt * (vt / (ROOT2 * csatfs)).ln()
        } else {
            f64::MAX
        };
        let vtd = K_OVER_Q * td;
        let vcritd = if csatfd != 0.0 {
            vtd * (vtd / (ROOT2 * csatfd)).ln()
        } else {
            f64::MAX
        };

        let temp_exp = (ts / model.tf).exp();
        let fl = model.flo * temp_exp;
        let delf = model.delfo * temp_exp;

        let tc_temp = |r: f64, t: f64| -> f64 {
            if r != 0.0 {
                let dt = t - tnom;
                r * (1.0 + model.tc1 * dt + model.tc2 * dt * dt)
            } else {
                0.0
            }
        };

        let t_rdi = tc_temp(model.rdi, td);
        let t_rsi = tc_temp(model.rsi, ts);
        let t_rg = tc_temp(model.rg, ts);
        let t_rs = tc_temp(model.rs, ts);
        let t_rd = tc_temp(model.rd, td);
        let t_ri = tc_temp(model.ri, ts);
        let t_rf = tc_temp(model.rf, td);

        let drain_conduct = if t_rd != 0.0 { 1.0 / t_rd } else { 0.0 };
        let source_conduct = if t_rs != 0.0 { 1.0 / t_rs } else { 0.0 };
        let gate_conduct = if t_rg != 0.0 { 1.0 / t_rg } else { 0.0 };
        let t_gi = if t_ri != 0.0 { 1.0 / t_ri } else { 0.0 };
        let t_gf = if t_rf != 0.0 { 1.0 / t_rf } else { 0.0 };

        Self {
            t_mu,
            t_vto,
            t_eta,
            t_lambda,
            t_lambdahf,
            t_rsi,
            t_rdi,
            drain_conduct,
            source_conduct,
            gate_conduct,
            t_gi,
            t_gf,
            gchi0,
            beta_inst,
            n0,
            nsb0,
            isatb0,
            cf,
            csatfs,
            csatfd,
            ggrwl,
            vcrits,
            vcritd,
            imax: CHARGE * model.nmax * model.vs * w,
            fl,
            delf,
            ts,
            td,
        }
    }
}

/// Resolved MESA instance in MNA system.
#[derive(Debug, Clone)]
pub struct MesaInstance {
    pub name: String,
    pub model: MesaModel,
    pub precomp: MesaPrecomp,
    pub w: f64,
    pub l: f64,
    pub drain_idx: Option<usize>,
    pub gate_idx: Option<usize>,
    pub source_idx: Option<usize>,
    pub drain_prime_idx: Option<usize>,
    pub gate_prime_idx: Option<usize>,
    pub source_prime_idx: Option<usize>,
    pub source_prm_prm_idx: Option<usize>,
    pub drain_prm_prm_idx: Option<usize>,
}

/// Result of MESA DC companion model evaluation.
#[derive(Debug, Clone)]
pub struct MesaCompanion {
    pub cdrain: f64,
    pub cg: f64,
    pub cd: f64,
    pub cgs: f64,
    pub cgd: f64,
    pub gm: f64,
    pub gds: f64,
    pub ggs: f64,
    pub ggd: f64,
    pub capgs: f64,
    pub capgd: f64,
    // For AC frequency-dependent lambda
    pub delidgch0: f64,
    pub delidvds0: f64,
    pub delidvds1: f64,
    pub gm0: f64,
    pub gm1: f64,
    pub gm2: f64,
    pub gds0: f64,
    pub vgs: f64,
    pub vgd: f64,
}

/// Internal result from level-specific DC model.
struct DcResult {
    cdrain: f64,
    gm: f64,
    gds: f64,
    capgs: f64,
    capgd: f64,
    delidgch0: f64,
    delidvds0: f64,
    delidvds1: f64,
    gm0: f64,
    gm1: f64,
    gm2: f64,
    gds0: f64,
}

/// Compute MESA companion model.
pub fn mesa_companion(inst: &MesaInstance, vgs: f64, vgd: f64, gmin: f64) -> MesaCompanion {
    let model = &inst.model;
    let pre = &inst.precomp;

    let vds = vgs - vgd;
    let vts = K_OVER_Q * pre.ts;
    let vtd = K_OVER_Q * pre.td;
    let vtes = model.n * vts;
    let vted = model.n * vtd;

    // Gate-source junction
    let arg_s = -vgs * model.del / vts;
    let earg_s = arg_s.exp();
    let evgs = (vgs / vtes).exp();
    let ggs = pre.csatfs * evgs / vtes + pre.ggrwl * earg_s * (1.0 - arg_s) + gmin;
    let cgs = pre.csatfs * (evgs - 1.0) + pre.ggrwl * vgs * earg_s + gmin * vgs;

    // Gate-drain junction
    let arg_d = -vgd * model.del / vtd;
    let earg_d = arg_d.exp();
    let evgd = (vgd / vted).exp();
    let ggd = pre.csatfd * evgd / vted + pre.ggrwl * earg_d * (1.0 - arg_d) + gmin;
    let cgd = pre.csatfd * (evgd - 1.0) + pre.ggrwl * vgd * earg_d + gmin * vgd;

    // Channel current (handle inverse mode)
    let (vds_abs, inverse) = if vds < 0.0 {
        (-vds, true)
    } else {
        (vds, false)
    };
    let vgs_ch = if inverse { vgd } else { vgs };
    let von = pre.t_vto;

    let mut dc = match model.level {
        3 => mesa2(model, pre, inst.w, inst.l, vgs_ch, vds_abs, von),
        4 => mesa3(model, pre, inst.w, inst.l, vgs_ch, vds_abs, von),
        _ => mesa1(model, pre, inst.w, inst.l, vgs_ch, vds_abs, von),
    };

    if inverse {
        dc.cdrain = -dc.cdrain;
        std::mem::swap(&mut dc.capgs, &mut dc.capgd);
    }

    let cd = dc.cdrain - cgd;

    MesaCompanion {
        cdrain: dc.cdrain,
        cg: cgs + cgd,
        cd,
        cgs,
        cgd,
        gm: dc.gm,
        gds: dc.gds,
        ggs,
        ggd,
        capgs: dc.capgs,
        capgd: dc.capgd,
        delidgch0: dc.delidgch0,
        delidvds0: dc.delidvds0,
        delidvds1: dc.delidvds1,
        gm0: dc.gm0,
        gm1: dc.gm1,
        gm2: dc.gm2,
        gds0: dc.gds0,
        vgs,
        vgd,
    }
}

/// MESA Level 2 DC model (mesa1 in ngspice).
fn mesa1(
    model: &MesaModel,
    pre: &MesaPrecomp,
    w: f64,
    l: f64,
    vgs: f64,
    vds: f64,
    von: f64,
) -> DcResult {
    let vt = K_OVER_Q * pre.ts;
    let etavth = pre.t_eta * vt;
    let rt = pre.t_rsi + pre.t_rdi;
    let vgt0 = vgs - von;
    let s = ((vgt0 - model.vsigmat) / model.vsigma).exp();
    let sigma = model.sigma0 / (1.0 + s);
    let vgt = vgt0 + sigma * vds;
    let mu = pre.t_mu + model.theta * vgt;
    let vl = model.vs / mu * l;
    let beta_eff = pre.beta_inst / (model.vpo + 3.0 * vl);
    let u = vgt / vt - 1.0;
    let t = (model.delta_sqr + u * u).sqrt();
    let vgte = 0.5 * vt * (2.0 + u + t);
    let a = 2.0 * beta_eff * vgte;
    let b = (-vgt / etavth).exp();

    let sqrt1 = if vgte > model.vpo {
        0.0
    } else {
        (1.0 - vgte / model.vpo).sqrt()
    };

    // ns = 1/(1/(ND*D*(1-sqrt1)) + b/n0)  — ngspice line: 1.0/ND/D/(1-sqrt1) + 1.0/n0*b
    let nsm = model.nd * model.d * (1.0 - sqrt1);
    let inv_nsm = if nsm.abs() < 1e-50 { 1e50 } else { 1.0 / nsm };
    let ns = 1.0 / (inv_nsm + b / pre.n0);

    if ns < 1.0e-38 {
        return DcResult {
            cdrain: 0.0,
            gm: 0.0,
            gds: 0.0,
            capgs: pre.cf,
            capgd: pre.cf,
            delidgch0: 0.0,
            delidvds0: 0.0,
            delidvds1: 0.0,
            gm0: 0.0,
            gm1: 0.0,
            gm2: 0.0,
            gds0: 0.0,
        };
    }

    let gchi = pre.gchi0 * mu * ns;
    let gch = gchi / (1.0 + gchi * rt);
    let f = (1.0 + 2.0 * a * pre.t_rsi).sqrt();
    let d_denom = 1.0 + a * pre.t_rsi + f;
    let e = 1.0 + model.tc * vgte;
    let isata = a * vgte / (d_denom * e);
    let isatb = pre.isatb0 * mu * (vgt / etavth).exp();
    let isat = isata * isatb / (isata + isatb);
    let vsate = isat / gch;
    let vdse = vds * (1.0 + (vds / vsate).powf(model.mc)).powf(-1.0 / model.mc);
    let m0 = model.m + model.alpha * vgte;
    let g = (vds / vsate).powf(m0);
    let h = (1.0 + g).powf(1.0 / m0);
    let delidgch0 = vds / h;
    let delidgch = delidgch0 * (1.0 + pre.t_lambda * vds);
    let cdrain = gch * delidgch;

    let temp_sqrt = if vgt > model.vpo {
        0.0
    } else {
        (1.0 - vgt / model.vpo).sqrt()
    };
    let cgc = w * l * EPSI_GAAS / (temp_sqrt + b) / model.d;
    let c1 = {
        let r = (vsate - vdse) / (2.0 * vsate - vdse);
        r * r
    };
    let capgs = pre.cf + 2.0 / 3.0 * cgc * (1.0 - c1);
    let c2 = {
        let r = vsate / (2.0 * vsate - vdse);
        r * r
    };
    let capgd = pre.cf + 2.0 / 3.0 * cgc * (1.0 - c2);

    // Derivatives
    let temp_gchi_rt = 1.0 + gchi * rt;
    let delgchgchi = 1.0 / (temp_gchi_rt * temp_gchi_rt);
    let delgchins = pre.gchi0 * mu;
    let delnsvgt = ns * ns / (pre.n0 * etavth) * b;
    let q = 1.0 - sqrt1;
    let delnsvgte = if sqrt1 == 0.0 {
        0.0
    } else {
        0.5 * ns * ns / (model.vpo * model.nd * model.d * sqrt1 * q * q)
    };
    let delvgtevgt = 0.5 * (1.0 + u / t);
    let delidvds0 = gch / h;
    let delidvds1 = if vds != 0.0 {
        cdrain * (vds / vsate).powf(m0 - 1.0) / (vsate * (1.0 + g))
    } else {
        0.0
    };
    let delidvds = delidvds0 * (1.0 + 2.0 * pre.t_lambda * vds) - delidvds1;
    let delidvsate = cdrain * g / (vsate * (1.0 + g));
    let delvsateisat = 1.0 / gch;
    let r_sum = isata + isatb;
    let r_sq = r_sum * r_sum;
    let delisatisata = isatb * isatb / r_sq;
    let v_term = 1.0 + 1.0 / f;
    let ddevgte = 2.0 * beta_eff * pre.t_rsi * v_term * e + d_denom * model.tc;
    let temp_d_sq = d_denom * d_denom * e * e;
    let delisatavgte = (2.0 * a * d_denom * e - a * vgte * ddevgte) / temp_d_sq;
    let delisatabeta = 2.0 * vgte * vgte * (d_denom * e - a * e * pre.t_rsi * v_term) / temp_d_sq;
    let delisatisatb = isata * isata / r_sq;
    let delvsategch = -vsate / gch;
    let dvgtvgs = 1.0 - model.sigma0 * vds * s / model.vsigma / ((1.0 + s) * (1.0 + s));
    let temp_gchi0_ns_theta = pre.gchi0 * ns * model.theta;
    let dgchivgt = delgchins * (delnsvgte * delvgtevgt + delnsvgt) + temp_gchi0_ns_theta;
    let dvgtevds = delvgtevgt * sigma;
    let dgchivds =
        delgchins * (delnsvgte * dvgtevds + delnsvgt * sigma) + temp_gchi0_ns_theta * sigma;
    let temp_beta_vl =
        delisatabeta * 3.0 * beta_eff * vl * model.theta / (mu * (model.vpo + 3.0 * vl));
    let disatavgt = delisatavgte * delvgtevgt + temp_beta_vl;
    let disatavds = delisatavgte * dvgtevds + temp_beta_vl * sigma;
    let disatbvgt = isatb / etavth + isatb / mu * model.theta;
    let p = delgchgchi * dgchivgt;
    let w_deriv = delgchgchi * dgchivds;
    let dvsatevgt =
        delvsateisat * (delisatisata * disatavgt + delisatisatb * disatbvgt) + delvsategch * p;
    let _dvsatevds = delvsateisat * (delisatisata * disatavds + delisatisatb * disatbvgt * sigma)
        + delvsategch * w_deriv;
    let (gmmadd, gdsmadd) = if model.alpha != 0.0 && vds != 0.0 {
        let gma = cdrain
            * ((1.0 + g).ln() / (m0 * m0) - g * (vds / vsate).ln() / (m0 * (1.0 + g)))
            * model.alpha
            * delvgtevgt;
        (gma, gma * sigma)
    } else {
        (0.0, 0.0)
    };
    let gm0_val = p;
    let gm1_val = delidvsate * dvsatevgt;
    let gm2_val = dvgtvgs;
    let gm_chain = delidgch * p + gm1_val;
    let gm = (gm_chain + gmmadd) * dvgtvgs;
    let gds0_val = delidvsate * _dvsatevds + delidgch * w_deriv + gdsmadd;
    let gds = delidvds + gds0_val;

    DcResult {
        cdrain,
        gm,
        gds,
        capgs,
        capgd,
        delidgch0,
        delidvds0,
        delidvds1,
        gm0: gm0_val,
        gm1: gm1_val,
        gm2: gm2_val,
        gds0: gds0_val,
    }
}

/// MESA Level 3 DC model (mesa2 in ngspice).
fn mesa2(
    model: &MesaModel,
    pre: &MesaPrecomp,
    w: f64,
    l: f64,
    vgs: f64,
    vds: f64,
    von: f64,
) -> DcResult {
    let vt = K_OVER_Q * pre.ts;
    let etavth = pre.t_eta * vt;
    let rt = pre.t_rsi + pre.t_rdi;
    let vgt0 = vgs - von;
    let s = ((vgt0 - model.vsigmat) / model.vsigma).exp();
    let sigma = model.sigma0 / (1.0 + s);
    let vgt = vgt0 + sigma * vds;
    let t = vgt / vt - 1.0;
    let q = (model.delta_sqr + t * t).sqrt();
    let vgte = 0.5 * vt * (2.0 + t + q);
    let a = 2.0 * model.beta * vgte;

    let (nsa, ca, delnsavgte) = if vgt > model.vpod {
        if vgte > model.vpo {
            (
                model.ndelta * model.th + model.ndu * model.du,
                EPSI_GAAS / model.du,
                0.0,
            )
        } else {
            let r = ((model.vpo - vgte) / model.vpou).sqrt();
            (
                model.ndelta * model.th + model.ndu * model.du * (1.0 - r),
                EPSI_GAAS / model.du / r,
                model.ndu * model.du / model.vpou / 2.0 / r,
            )
        }
    } else if model.vpod - vgte < 0.0 {
        (
            model.ndelta * model.th * (1.0 - model.du / model.th),
            EPSI_GAAS / model.du,
            0.0,
        )
    } else {
        let r = (1.0 + model.ndu / model.ndelta * (model.vpod - vgte) / model.vpou).sqrt();
        (
            model.ndelta * model.th * (1.0 - model.du / model.th * (r - 1.0)),
            EPSI_GAAS / model.du / r,
            model.du * model.ndu / 2.0 / model.vpou / r,
        )
    };

    let b = (vgt / etavth).exp();
    let cb = EPSI_GAAS / (model.du + model.th) * b;
    let nsb = pre.nsb0 * b;
    let delnsbvgt = nsb / etavth;
    let ns = nsa * nsb / (nsa + nsb);

    if ns < 1.0e-38 {
        return DcResult {
            cdrain: 0.0,
            gm: 0.0,
            gds: 0.0,
            capgs: pre.cf,
            capgd: pre.cf,
            delidgch0: 0.0,
            delidvds0: 0.0,
            delidvds1: 0.0,
            gm0: 0.0,
            gm1: 0.0,
            gm2: 0.0,
            gds0: 0.0,
        };
    }

    let gchi = pre.gchi0 * ns;
    let gch = gchi / (1.0 + gchi * rt);
    let f = (1.0 + 2.0 * a * pre.t_rsi).sqrt();
    let d_denom = 1.0 + a * pre.t_rsi + f;
    let e = 1.0 + model.tc * vgte;
    let isata = a * vgte / d_denom / e;
    let isatb = pre.isatb0 * b;
    let isat = isata * isatb / (isata + isatb);
    let vsate = isat / gch;
    let vdse = vds * (1.0 + (vds / vsate).powf(model.mc)).powf(-1.0 / model.mc);
    let g = (vds / vsate).powf(model.m);
    let h = (1.0 + g).powf(1.0 / model.m);
    let delidgch0 = vds / h;
    let delidgch = delidgch0 * (1.0 + pre.t_lambda * vds);
    let cdrain = gch * delidgch;

    let cgc = w * l * ca * cb / (ca + cb);
    let c1 = {
        let r = (vsate - vdse) / (2.0 * vsate - vdse);
        r * r
    };
    let capgs = pre.cf + 2.0 / 3.0 * cgc * (1.0 - c1);
    let c2 = {
        let r = vsate / (2.0 * vsate - vdse);
        r * r
    };
    let capgd = pre.cf + 2.0 / 3.0 * cgc * (1.0 - c2);

    let delvgtevgt = 0.5 * (1.0 + t / q);
    let delidvds0 = gch / h;
    let delidvds1 = if vds != 0.0 {
        cdrain * (vds / vsate).powf(model.m - 1.0) / vsate / (1.0 + g)
    } else {
        0.0
    };
    let delidvds = delidvds0 * (1.0 + 2.0 * pre.t_lambda * vds) - delidvds1;
    let delgchgchi = {
        let t = 1.0 + gchi * rt;
        1.0 / (t * t)
    };
    let delgchins = pre.gchi0;
    let r_ns = nsa + nsb;
    let r_ns_sq = r_ns * r_ns;
    let delnsvgt = (nsb * nsb * delvgtevgt * delnsavgte + nsa * nsa * delnsbvgt) / r_ns_sq;
    let delidvsate = cdrain * g / vsate / (1.0 + g);
    let delvsateisat = 1.0 / gch;
    let r_isat = isata + isatb;
    let r_isat_sq = r_isat * r_isat;
    let delisatisata = isatb * isatb / r_isat_sq;
    let ddevgte = 2.0 * model.beta * pre.t_rsi * (1.0 + 1.0 / f) * e + d_denom * model.tc;
    let delisatavgte = (2.0 * a * d_denom * e - a * vgte * ddevgte) / d_denom / d_denom / e / e;
    let delisatisatb = isata * isata / r_isat_sq;
    let delisatbvgt = isatb / etavth;
    let delvsategch = -vsate / gch;
    let delvgtvgs = 1.0 - model.sigma0 * vds * s / model.vsigma / (1.0 + s) / (1.0 + s);
    let p = delgchgchi * delgchins * delnsvgt;
    let delvsatevgt = delvsateisat
        * (delisatisata * delisatavgte * delvgtevgt + delisatisatb * delisatbvgt)
        + delvsategch * p;
    let gm0_val = p;
    let gm1_val = delidvsate * delvsatevgt;
    let gm2_val = delvgtvgs;
    let gm_chain = delidgch * p + gm1_val;
    let gm = gm_chain * delvgtvgs;
    let gds0_val = gm_chain * sigma;
    let gds = delidvds + gds0_val;

    DcResult {
        cdrain,
        gm,
        gds,
        capgs,
        capgd,
        delidgch0,
        delidvds0,
        delidvds1,
        gm0: gm0_val,
        gm1: gm1_val,
        gm2: gm2_val,
        gds0: gds0_val,
    }
}

/// MESA Level 4 DC model (mesa3 in ngspice).
fn mesa3(
    model: &MesaModel,
    pre: &MesaPrecomp,
    _w: f64,
    l: f64,
    vgs: f64,
    vds: f64,
    von: f64,
) -> DcResult {
    let vt = K_OVER_Q * pre.ts;
    let etavth = pre.t_eta * vt;
    let vl = model.vs / pre.t_mu * l;
    let rt = pre.t_rsi + pre.t_rdi;
    let vgt0 = vgs - von;
    let s = ((vgt0 - model.vsigmat) / model.vsigma).exp();
    let sigma = model.sigma0 / (1.0 + s);
    let vgt = vgt0 + sigma * vds;
    let u = 0.5 * vgt / vt - 1.0;
    let t = (model.delta_sqr + u * u).sqrt();
    let vgte = vt * (2.0 + u + t);
    let b = (vgt / etavth).exp();
    let nsm = 2.0 * pre.n0 * (1.0 + 0.5 * b).ln();

    if nsm < 1.0e-38 {
        return DcResult {
            cdrain: 0.0,
            gm: 0.0,
            gds: 0.0,
            capgs: pre.cf,
            capgd: pre.cf,
            delidgch0: 0.0,
            delidvds0: 0.0,
            delidvds1: 0.0,
            gm0: 0.0,
            gm1: 0.0,
            gm2: 0.0,
            gds0: 0.0,
        };
    }

    let c_nsm = (nsm / model.nmax).powf(model.gamma);
    let q_nsm = (1.0 + c_nsm).powf(1.0 / model.gamma);
    let ns = nsm / q_nsm;
    let gchi = pre.gchi0 * ns;
    let gch = gchi / (1.0 + gchi * rt);
    let gchim = pre.gchi0 * nsm;
    let h_term = (1.0 + 2.0 * gchim * model.rsi + vgte * vgte / (vl * vl)).sqrt();
    let p_denom = 1.0 + gchim * pre.t_rsi + h_term;
    let isatm = gchim * vgte / p_denom;
    let g_imax = (isatm / pre.imax).powf(model.gamma);
    let isat = isatm / (1.0 + g_imax).powf(1.0 / model.gamma);
    let vsate = isat / gch;
    let vdse = vds * (1.0 + (vds / vsate).powf(model.mc)).powf(-1.0 / model.mc);
    let d_pow = (vds / vsate).powf(model.m);
    let e_pow = (1.0 + d_pow).powf(1.0 / model.m);
    let delidgch_val = vds * (1.0 + pre.t_lambda * vds) / e_pow;
    let cdrain = gch * delidgch_val;

    let cgcm = 1.0
        / (1.0 / model.cas * model.d / model.epsi
            + 1.0 / model.cbs * etavth / CHARGE / pre.n0 * (-vgt / etavth).exp());
    let w = _w;
    let cgc = w * l * cgcm / (1.0 + c_nsm).powf(1.0 + 1.0 / model.gamma);
    let a1 = {
        let r = (vsate - vdse) / (2.0 * vsate - vdse);
        r * r
    };
    let capgs = pre.cf + 2.0 / 3.0 * cgc * (1.0 - a1);
    let a2 = {
        let r = vsate / (2.0 * vsate - vdse);
        r * r
    };
    let capgd = pre.cf + 2.0 / 3.0 * cgc * (1.0 - a2);

    let delidvsate = cdrain * d_pow / vsate / (1.0 + d_pow);
    let delidvds = gch * (1.0 + 2.0 * pre.t_lambda * vds) / e_pow
        - cdrain * (vds / vsate).powf(model.m - 1.0) / (vsate * (1.0 + d_pow));
    let a_gchi = 1.0 + gchi * rt;
    let delgchgchi = 1.0 / (a_gchi * a_gchi);
    let delgchins = pre.gchi0;
    let delnsnsm = ns / nsm * (1.0 - c_nsm / (1.0 + c_nsm));
    let delnsmvgt = pre.n0 / etavth / (1.0 / b + 0.5);
    let delvgtevgt = 0.5 * (1.0 + u / t);
    let delvsateisat = 1.0 / gch;
    let delisatisatm = isat / isatm * (1.0 - g_imax / (1.0 + g_imax));
    let delisatmvgte = gchim * (p_denom - vgte * vgte / (vl * vl * h_term)) / (p_denom * p_denom);
    let delvsategch = -vsate / gch;
    let delisatmgchim =
        vgte * (p_denom - gchim * pre.t_rsi * (1.0 + 1.0 / h_term)) / (p_denom * p_denom);
    let delvgtvgs = 1.0 - vds * model.sigma0 / model.vsigma * s / ((1.0 + s) * (1.0 + s));
    let p = delgchgchi * delgchins * delnsnsm * delnsmvgt;
    let delvsatevgt = delvsateisat
        * delisatisatm
        * (delisatmvgte * delvgtevgt + delisatmgchim * pre.gchi0 * delnsmvgt)
        + delvsategch * p;
    let g_deriv = delidgch_val * p + delidvsate * delvsatevgt;
    let gm = g_deriv * delvgtvgs;
    let gds = delidvds + g_deriv * sigma;

    let delidgch0 = vds * (1.0 + pre.t_lambda * vds) / e_pow;
    let delidvds0 = gch / e_pow;
    let delidvds1 = if vds != 0.0 {
        cdrain * (vds / vsate).powf(model.m - 1.0) / (vsate * (1.0 + d_pow))
    } else {
        0.0
    };

    DcResult {
        cdrain,
        gm,
        gds,
        capgs,
        capgd,
        delidgch0,
        delidvds0,
        delidvds1,
        gm0: p,
        gm1: delidvsate * delvsatevgt,
        gm2: delvgtvgs,
        gds0: g_deriv * sigma,
    }
}

impl MesaInstance {
    /// Extract gate-source and gate-drain junction voltages from solution.
    pub fn junction_voltages(&self, solution: &[f64]) -> (f64, f64) {
        let v_gp = self.gate_prime_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_sp = self.source_prime_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_dp = self.drain_prime_idx.map(|i| solution[i]).unwrap_or(0.0);
        let vgs = v_gp - v_sp;
        let vgd = v_gp - v_dp;
        (vgs, vgd)
    }
}

/// Stamp the MESA companion model into the MNA system.
pub fn stamp_mesa(
    comp: &MesaCompanion,
    inst: &MesaInstance,
    matrix: &mut crate::sparse::SparseMatrix,
    rhs: &mut [f64],
) {
    let pre = &inst.precomp;
    let dp = inst.drain_prime_idx;
    let gp = inst.gate_prime_idx;
    let sp = inst.source_prime_idx;
    let spp = inst.source_prm_prm_idx;
    let dpp = inst.drain_prm_prm_idx;
    let d_ext = inst.drain_idx;
    let g_ext = inst.gate_idx;
    let s_ext = inst.source_idx;

    let vgs = comp.vgs;
    let vgd = comp.vgd;
    let vds = vgs - vgd;

    // In DC: ggspp=ggdpp=cgspp=cgdpp=0
    let ccorr = inst.model.ag * (comp.cgs - comp.cgd);
    let ceqgd_full = comp.cgd - comp.ggd * vgd;
    let ceqgs_full = comp.cgs - comp.ggs * vgs;
    let cdreq = (comp.cd + comp.cgd) - comp.gds * vds - comp.gm * vgs;

    // RHS
    if let Some(gp_i) = gp {
        rhs[gp_i] += -ceqgs_full - ceqgd_full;
    }
    if let Some(dp_i) = dp {
        rhs[dp_i] += -cdreq + ceqgd_full + ccorr;
    }
    if let Some(sp_i) = sp {
        rhs[sp_i] += cdreq + ceqgs_full - ccorr;
    }
    // spp, dpp RHS are zero in DC (ggspp=ggdpp=0)

    // Series resistance stamps: each stamps individual matrix entries matching ngspice
    // d,d += drain_conduct; dp,dp += drain_conduct; d,dp -= drain_conduct; dp,d -= drain_conduct
    if let Some(d) = d_ext {
        matrix.add(d, d, pre.drain_conduct);
    }
    if let Some(s) = s_ext {
        matrix.add(s, s, pre.source_conduct);
    }
    if let Some(g) = g_ext {
        matrix.add(g, g, pre.gate_conduct);
    }
    if let (Some(d), Some(di)) = (d_ext, dp) {
        matrix.add(d, di, -pre.drain_conduct);
        matrix.add(di, d, -pre.drain_conduct);
    }
    if let (Some(s), Some(si)) = (s_ext, sp) {
        matrix.add(s, si, -pre.source_conduct);
        matrix.add(si, s, -pre.source_conduct);
    }
    if let (Some(g), Some(gi)) = (g_ext, gp) {
        matrix.add(g, gi, -pre.gate_conduct);
        matrix.add(gi, g, -pre.gate_conduct);
    }

    // Feedback resistance stamps (ri, rf create spp, dpp)
    // spp diagonal: t_gi (+ ggspp in transient, 0 in DC)
    if let Some(spp_i) = spp {
        matrix.add(spp_i, spp_i, pre.t_gi);
    }
    // dpp diagonal: t_gf (+ ggdpp in transient, 0 in DC)
    if let Some(dpp_i) = dpp {
        matrix.add(dpp_i, dpp_i, pre.t_gf);
    }

    // Gate' diagonal (ggd + ggs + gate_conduct + ggspp + ggdpp; spp/dpp terms 0 in DC)
    if let Some(gp_i) = gp {
        matrix.add(gp_i, gp_i, comp.ggd + comp.ggs + pre.gate_conduct);
    }
    // Drain' diagonal (gds + ggd + drain_conduct + t_gf)
    if let Some(dp_i) = dp {
        matrix.add(
            dp_i,
            dp_i,
            comp.gds + comp.ggd + pre.drain_conduct + pre.t_gf,
        );
    }
    // Source' diagonal (gds + gm + ggs + source_conduct + t_gi)
    if let Some(sp_i) = sp {
        matrix.add(
            sp_i,
            sp_i,
            comp.gds + comp.gm + comp.ggs + pre.source_conduct + pre.t_gi,
        );
    }

    // Off-diagonal entries (single matrix elements, not 2x2 conductances)
    // gp-dp: -ggd
    if let (Some(gi), Some(di)) = (gp, dp) {
        matrix.add(gi, di, -comp.ggd);
    }
    // gp-sp: -ggs
    if let (Some(gi), Some(si)) = (gp, sp) {
        matrix.add(gi, si, -comp.ggs);
    }
    // dp-gp: gm - ggd
    if let (Some(di), Some(gi)) = (dp, gp) {
        matrix.add(di, gi, comp.gm - comp.ggd);
    }
    // dp-sp: -gds - gm
    if let (Some(di), Some(si)) = (dp, sp) {
        matrix.add(di, si, -comp.gds - comp.gm);
    }
    // sp-gp: -ggs - gm
    if let (Some(si), Some(gi)) = (sp, gp) {
        matrix.add(si, gi, -comp.ggs - comp.gm);
    }
    // sp-dp: -gds
    if let (Some(si), Some(di)) = (sp, dp) {
        matrix.add(si, di, -comp.gds);
    }

    // sp <-> spp coupling (Ri)
    if let (Some(si), Some(spi)) = (sp, spp) {
        matrix.add(si, spi, -pre.t_gi);
        matrix.add(spi, si, -pre.t_gi);
    }
    // dp <-> dpp coupling (Rf)
    if let (Some(di), Some(dpi)) = (dp, dpp) {
        matrix.add(di, dpi, -pre.t_gf);
        matrix.add(dpi, di, -pre.t_gf);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mesa_model_defaults() {
        let m = MesaModel::new();
        assert_eq!(m.level, 2);
        assert!((m.vto - (-1.26)).abs() < 1e-10);
        assert!((m.lambda - 0.045).abs() < 1e-10);
        assert!((m.eta - 1.73).abs() < 1e-10);
        assert!((m.m - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_mesa_model_parsing() {
        use thevenin_types::Netlist;
        let netlist = Netlist::parse(
            "test\n\
             z1 2 3 0 mesmod l=1u w=20u\n\
             .model mesmod nmf level=2\n\
             .end",
        )
        .unwrap();
        let mesa_count = netlist
            .elements()
            .filter(|e| matches!(e.kind, thevenin_types::ElementKind::Mesa { .. }))
            .count();
        assert_eq!(mesa_count, 1);
    }

    #[test]
    fn test_mesa_precomp() {
        let mut m = MesaModel::new();
        m.astar = 0.0;
        m.precompute();
        let pre = MesaPrecomp::compute(&m, 300.15, 300.15, 300.15, 20e-6, 1e-6);
        assert!(pre.gchi0 > 0.0);
        assert!(pre.beta_inst > 0.0);
        assert!(pre.cf > 0.0);
    }

    #[test]
    fn test_mesa_level2_cutoff() {
        let mut m = MesaModel::new();
        m.astar = 0.0;
        m.precompute();
        let pre = MesaPrecomp::compute(&m, 300.15, 300.15, 300.15, 20e-6, 1e-6);
        let dc = mesa1(&m, &pre, 20e-6, 1e-6, -2.0, 1.0, m.vto);
        assert!(
            dc.cdrain.abs() < 1e-6,
            "Should be near cutoff, got {}",
            dc.cdrain
        );
    }

    #[test]
    fn test_mesa_level2_on() {
        let mut m = MesaModel::new();
        m.astar = 0.0;
        m.precompute();
        let pre = MesaPrecomp::compute(&m, 300.15, 300.15, 300.15, 20e-6, 1e-6);
        let dc = mesa1(&m, &pre, 20e-6, 1e-6, 0.0, 1.0, m.vto);
        assert!(dc.cdrain > 0.0, "Should conduct, got {}", dc.cdrain);
        assert!(dc.gm > 0.0, "gm should be positive");
        assert!(dc.gds > 0.0, "gds should be positive");
    }
}
