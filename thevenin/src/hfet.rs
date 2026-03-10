//! HFET (Heterojunction FET) device model.
//!
//! Implements the HFET1 (level=5) model from ngspice,
//! used with `.model name NHFET/PHFET level=5` syntax.
//!
//! Reference: ngspice hfet1/ device directory (hfetload.c, hfetdefs.h, hfettemp.c).

use thevenin_types::{Expr, ModelDef};

use crate::mna::stamp_conductance;
use crate::sparse::SparseMatrix;

const CHARGE: f64 = 1.602_176_634e-19;
const BOLTZMANN: f64 = 1.380_649e-23;
const K_OVER_Q: f64 = BOLTZMANN / CHARGE;

fn expr_val(e: &Expr) -> f64 {
    match e {
        Expr::Num(v) => *v,
        Expr::Param(_) | Expr::Brace(_) => 0.0,
    }
}

/// HFET device type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HfetType {
    Nhfet = 1,
    Phfet = -1,
}

/// HFET model parameters.
#[derive(Debug, Clone)]
pub struct HfetModel {
    pub device_type: HfetType,
    pub vt0: f64,     // Threshold voltage
    pub lambda: f64,  // Output conductance parameter
    pub eta: f64,     // Subthreshold ideality factor
    pub m: f64,       // Knee shape parameter
    pub mc: f64,      // Capacitance knee shape
    pub gamma: f64,   // Carrier saturation parameter
    pub sigma0: f64,  // Threshold voltage modulation coefficient
    pub vsigmat: f64, // Threshold of sigma modulation
    pub vsigma: f64,  // Sigma smoothing parameter
    pub mu: f64,      // Channel mobility (cm²/V·s)
    pub di: f64,      // Channel depth
    pub delta: f64,   // Smoothing parameter
    pub vs: f64,      // Saturation velocity (cm/s)
    pub nmax: f64,    // Max channel carrier concentration
    pub deltad: f64,  // Depletion thickness correction
    pub rd: f64,      // Drain resistance
    pub rs: f64,      // Source resistance
    pub rg: f64,      // Gate resistance
    pub rdi: f64,     // Internal drain resistance
    pub rsi: f64,     // Internal source resistance
    pub rgs: f64,     // Gate-source resistance
    pub rgd: f64,     // Gate-drain resistance
    pub ri: f64,      // Drain-source feedback resistance
    pub rf: f64,      // Drain-source feedforward resistance
    pub epsi: f64,    // Relative permittivity (F/m)
    pub a1: f64,      // Correction current coefficient 1
    pub a2: f64,      // Correction current coefficient 2
    pub mv1: f64,     // Correction saturation exponent
    pub p: f64,       // PM - capacitance charge partition
    pub kappa: f64,
    pub delf: f64,
    pub fgds: f64,
    pub tf: f64,      // Frequency parameter temperature
    pub cds: f64,     // Drain-source capacitance
    pub phib: f64,    // Barrier height (eV)
    pub talpha: f64,  // Temperature rise coefficient
    pub mt1: f64,     // Gate-drain temperature exponent 1
    pub mt2: f64,     // Gate-drain temperature exponent 2
    pub ck1: f64,     // Knee voltage coefficient 1
    pub ck2: f64,     // Knee voltage coefficient 2
    pub cm1: f64,     // Max voltage coefficient 1
    pub cm2: f64,     // Max voltage coefficient 2
    pub cm3: f64,     // Correction max coefficient
    pub astar: f64,   // Richardson constant
    pub eta1: f64,    // Parasitic charge ideality 1
    pub d1: f64,      // Parasitic channel depth 1
    pub vt1: f64,     // Parasitic threshold 1
    pub eta2: f64,    // Second channel ideality
    pub d2: f64,      // Second channel depth
    pub vt2: f64,     // Second channel threshold
    pub ggr: f64,     // Gate-source recombination
    pub del: f64,     // Temperature-dependent gate recombination
    pub klambda: f64, // Lambda temperature coefficient
    pub kmu: f64,     // Mu temperature coefficient
    pub kvto: f64,    // Vto temperature coefficient
    pub js1d: f64,    // Drain junction saturation current density 1
    pub js2d: f64,    // Drain junction saturation current density 2
    pub js1s: f64,    // Source junction saturation current density 1
    pub js2s: f64,    // Source junction saturation current density 2
    pub m1d: f64,     // Drain diode ideality 1
    pub m2d: f64,     // Drain diode ideality 2
    pub m1s: f64,     // Source diode ideality 1
    pub m2s: f64,     // Source diode ideality 2
    pub gatemod: i32, // Gate model selector

    // Temperature-dependent flags
    pub eta2_given: bool,
    pub d2_given: bool,
    pub vt1_given: bool,
    pub vt2_given: bool,

    // Precomputed
    pub drain_conduct: f64,
    pub source_conduct: f64,
    pub gate_conduct: f64,
    pub gi: f64,
    pub gf: f64,
    pub delta_sqr: f64,
}

impl Default for HfetModel {
    fn default() -> Self {
        Self::new(HfetType::Nhfet)
    }
}

impl HfetModel {
    pub fn new(device_type: HfetType) -> Self {
        let sign = device_type as i32 as f64;
        Self {
            device_type,
            vt0: if sign > 0.0 { 0.15 } else { -0.15 },
            lambda: 0.15,
            eta: if sign > 0.0 { 1.28 } else { 1.4 },
            m: 3.0,
            mc: 3.0,
            gamma: 3.0,
            sigma0: 0.057,
            vsigmat: 0.3,
            vsigma: 0.1,
            mu: if sign > 0.0 { 0.4 } else { 0.03 },
            di: 0.04e-6,
            delta: 3.0,
            vs: if sign > 0.0 { 1.5e5 } else { 0.8e5 },
            nmax: 2e16,
            deltad: 4.5e-9,
            rd: 0.0,
            rs: 0.0,
            rg: 0.0,
            rdi: 0.0,
            rsi: 0.0,
            rgs: 90.0,
            rgd: 90.0,
            ri: 0.0,
            rf: 0.0,
            epsi: 12.244 * 8.85418e-12, // GaAs permittivity
            a1: 0.0,
            a2: 0.0,
            mv1: 3.0,
            p: 1.0,
            kappa: 0.0,
            delf: 0.0,
            fgds: 0.0,
            tf: 300.15,
            cds: 0.0,
            phib: 0.5 * CHARGE,
            talpha: 1200.0,
            mt1: 3.5,
            mt2: 9.9,
            ck1: 1.0,
            ck2: 0.0,
            cm1: 3.0,
            cm2: 0.0,
            cm3: 0.17,
            astar: 4.0e4,
            eta1: 2.0,
            d1: 0.03e-6,
            vt1: f64::MAX,
            eta2: 2.0,
            d2: 0.2e-6,
            vt2: f64::MAX,
            ggr: 40.0,
            del: 0.04,
            klambda: 0.0,
            kmu: 0.0,
            kvto: 0.0,
            js1d: 1.0,
            js2d: 1.15e6,
            js1s: 1.0,
            js2s: 1.15e6,
            m1d: 1.32,
            m2d: 6.9,
            m1s: 1.32,
            m2s: 6.9,
            gatemod: 0,
            eta2_given: false,
            d2_given: false,
            vt1_given: false,
            vt2_given: false,
            drain_conduct: 0.0,
            source_conduct: 0.0,
            gate_conduct: 0.0,
            gi: 0.0,
            gf: 0.0,
            delta_sqr: 9.0,
        }
    }

    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let kind = model_def.kind.to_uppercase();
        let device_type = if kind == "PHFET" {
            HfetType::Phfet
        } else {
            HfetType::Nhfet
        };
        let mut m = Self::new(device_type);

        for p in &model_def.params {
            let val = expr_val(&p.value);
            match p.name.to_uppercase().as_str() {
                "VTO" | "VT0" => m.vt0 = val,
                "LAMBDA" => m.lambda = val,
                "ETA" => m.eta = val,
                "M" => m.m = val,
                "MC" => m.mc = val,
                "GAMMA" => m.gamma = val,
                "SIGMA0" => m.sigma0 = val,
                "VSIGMAT" => m.vsigmat = val,
                "VSIGMA" => m.vsigma = val,
                "MU" => m.mu = val,
                "DI" => m.di = val,
                "DELTA" => m.delta = val,
                "VS" => m.vs = val,
                "NMAX" => m.nmax = val,
                "DELTAD" => m.deltad = val,
                "RD" => m.rd = val,
                "RS" => m.rs = val,
                "RG" => m.rg = val,
                "RDI" => m.rdi = val,
                "RSI" => m.rsi = val,
                "RGS" => m.rgs = val,
                "RGD" => m.rgd = val,
                "RI" => m.ri = val,
                "RF" => m.rf = val,
                "EPSI" => m.epsi = val,
                "A1" => m.a1 = val,
                "A2" => m.a2 = val,
                "MV1" => m.mv1 = val,
                "P" => m.p = val,
                "KAPPA" => m.kappa = val,
                "DELF" => m.delf = val,
                "FGDS" => m.fgds = val,
                "TF" => m.tf = val,
                "CDS" => m.cds = val,
                "PHIB" => m.phib = val,
                "TALPHA" => m.talpha = val,
                "MT1" => m.mt1 = val,
                "MT2" => m.mt2 = val,
                "CK1" => m.ck1 = val,
                "CK2" => m.ck2 = val,
                "CM1" => m.cm1 = val,
                "CM2" => m.cm2 = val,
                "CM3" => m.cm3 = val,
                "ASTAR" => m.astar = val,
                "ETA1" => m.eta1 = val,
                "D1" => m.d1 = val,
                "VT1" => {
                    m.vt1 = val;
                    m.vt1_given = true;
                }
                "ETA2" => {
                    m.eta2 = val;
                    m.eta2_given = true;
                }
                "D2" => {
                    m.d2 = val;
                    m.d2_given = true;
                }
                "VT2" => {
                    m.vt2 = val;
                    m.vt2_given = true;
                }
                "GGR" => m.ggr = val,
                "DEL" => m.del = val,
                "KLAMBDA" => m.klambda = val,
                "KMU" => m.kmu = val,
                "KVTO" => m.kvto = val,
                "JS1D" => m.js1d = val,
                "JS2D" => m.js2d = val,
                "JS1S" => m.js1s = val,
                "JS2S" => m.js2s = val,
                "M1D" => m.m1d = val,
                "M2D" => m.m2d = val,
                "M1S" => m.m1s = val,
                "M2S" => m.m2s = val,
                "GATEMOD" => m.gatemod = val as i32,
                _ => {}
            }
        }

        m.precompute();
        m
    }

    fn precompute(&mut self) {
        self.drain_conduct = if self.rd != 0.0 { 1.0 / self.rd } else { 0.0 };
        self.source_conduct = if self.rs != 0.0 { 1.0 / self.rs } else { 0.0 };
        self.gate_conduct = if self.rg != 0.0 { 1.0 / self.rg } else { 0.0 };
        self.gi = if self.ri != 0.0 { 1.0 / self.ri } else { 0.0 };
        self.gf = if self.rf != 0.0 { 1.0 / self.rf } else { 0.0 };
        self.delta_sqr = self.delta * self.delta;
        // Sign the threshold voltage by device type
        self.vt0 *= self.device_type as i32 as f64;
        // Set defaults for vt2, vt1 if not given
        if !self.vt2_given {
            self.vt2 = self.vt0;
        }
        if !self.vt1_given {
            self.vt1 = self.vt0 + CHARGE * self.nmax * self.di / self.epsi;
        }
    }

    /// Number of internal nodes needed by this model.
    pub fn internal_node_count(&self) -> usize {
        // From hfetsetup.c: create internal nodes only when resistance > 0
        let mut n = 0;
        if self.rd != 0.0 {
            n += 1; // drainPrime
        }
        if self.rs != 0.0 {
            n += 1; // sourcePrime
        }
        if self.rg != 0.0 {
            n += 1; // gatePrime
        }
        if self.rf != 0.0 {
            n += 1; // drainPrmPrm
        }
        if self.ri != 0.0 {
            n += 1; // sourcePrmPrm
        }
        n
    }
}

/// Precomputed instance-level parameters.
#[derive(Debug, Clone)]
pub struct HfetPrecomp {
    pub n0: f64,
    pub n01: f64,
    pub n02: f64,
    pub gchi0: f64,
    pub cf: f64,
    pub imax: f64,
    pub is1d: f64,
    pub is2d: f64,
    pub is1s: f64,
    pub is2s: f64,
    pub iso: f64,
    pub ggrwl: f64,
    pub vcrit: f64,
    pub t_lambda: f64,
    pub t_mu: f64,
    pub t_vto: f64,
    pub temp: f64,
}

impl HfetPrecomp {
    pub fn compute(model: &HfetModel, temp: f64, tnom: f64, w: f64, l: f64) -> Self {
        let vt = K_OVER_Q * temp;
        let t_lambda = model.lambda + model.klambda * (temp - tnom);
        let t_mu = model.mu - model.kmu * (temp - tnom);
        let t_vto = model.vt0 - model.kvto * (temp - tnom);
        let n0 = model.epsi * model.eta * vt / (2.0 * CHARGE * (model.di + model.deltad));
        let n01 = model.epsi * model.eta1 * vt / (2.0 * CHARGE * model.d1);
        let n02 = if model.eta2_given && model.d2 != 0.0 {
            model.epsi * model.eta2 * vt / (2.0 * CHARGE * model.d2)
        } else {
            0.0
        };
        let gchi0 = CHARGE * w * t_mu / l;
        let cf = 0.5 * model.epsi * w;
        let imax = CHARGE * model.nmax * model.vs * w;
        let is1d = model.js1d * w * l / 2.0;
        let is2d = model.js2d * w * l / 2.0;
        let is1s = model.js1s * w * l / 2.0;
        let is2s = model.js2s * w * l / 2.0;
        let iso = model.astar * w * l / 2.0;
        let ggrwl = model.ggr * l * w / 2.0;

        let vcrit = if model.gatemod == 0 {
            if is1s != 0.0 {
                vt * (vt / (std::f64::consts::SQRT_2 * is1s)).ln()
            } else {
                f64::MAX
            }
        } else if iso != 0.0 {
            vt * (vt / (std::f64::consts::SQRT_2 * iso)).ln()
        } else {
            f64::MAX
        };

        Self {
            n0,
            n01,
            n02,
            gchi0,
            cf,
            imax,
            is1d,
            is2d,
            is1s,
            is2s,
            iso,
            ggrwl,
            vcrit,
            t_lambda,
            t_mu,
            t_vto,
            temp,
        }
    }
}

/// Resolved HFET instance for simulation.
#[derive(Debug, Clone)]
pub struct HfetInstance {
    pub name: String,
    pub drain_idx: Option<usize>,
    pub gate_idx: Option<usize>,
    pub source_idx: Option<usize>,
    pub gate_prime_idx: Option<usize>,
    pub drain_prime_idx: Option<usize>,
    pub source_prime_idx: Option<usize>,
    pub drain_prm_prm_idx: Option<usize>,
    pub source_prm_prm_idx: Option<usize>,
    pub model: HfetModel,
    pub precomp: HfetPrecomp,
    pub w: f64,
    pub l: f64,
}

impl HfetInstance {
    /// Get junction voltages from solution.
    pub fn junction_voltages(&self, solution: &[f64]) -> (f64, f64) {
        let gp = self
            .gate_prime_idx
            .or(self.gate_idx)
            .map_or(0.0, |i| solution[i]);
        let dp = self
            .drain_prime_idx
            .map_or(self.drain_idx.map_or(0.0, |i| solution[i]), |i| solution[i]);
        let sp = self
            .source_prime_idx
            .map_or(self.source_idx.map_or(0.0, |i| solution[i]), |i| {
                solution[i]
            });
        let sign = self.model.device_type as i32 as f64;
        (sign * (gp - sp), sign * (gp - dp))
    }
}

/// HFET companion model result.
#[derive(Debug, Clone)]
pub struct HfetCompanion {
    pub gm: f64,
    pub gds: f64,
    pub ggs: f64,
    pub ggd: f64,
    pub cg: f64,
    pub cd: f64,
    pub cgd_current: f64,
    pub cgs_current: f64,
    pub capgs: f64,
    pub capgd: f64,
    // Gate-drain model outputs
    pub gmg: f64,
    pub gmd: f64,
}

/// Gate leakage current model (from hfetload.c `leak` function).
#[expect(clippy::too_many_arguments)]
fn leak(gmin: f64, vt: f64, v: f64, rs: f64, is1: f64, is2: f64, m1: f64, m2: f64) -> (f64, f64) {
    let vt1 = vt * m1;
    let vt2 = vt * m2;

    if v > -10.0 * vt1 {
        let vteff = vt1 + vt2;
        let iseff = is2 * (is1 / is2).powf(m1 / (m1 + m2));
        let (iaprox1, iaprox2) = if rs > 0.0 {
            let unorm1 = (v + rs * is1) / vt1 + (rs * is1 / vt1).ln();
            let ia1 = vt1 * diode_fn(unorm1) / rs - is1;
            let unorm2 = (v + rs * iseff) / vteff + (rs * iseff / vteff).ln();
            let ia2 = vteff * diode_fn(unorm2) / rs - iseff;
            (ia1, ia2)
        } else {
            let ia1 = is1 * ((v / vt1).exp() - 1.0);
            let ia2 = iseff * ((v / vteff).exp() - 1.0);
            (ia1, ia2)
        };

        let iaprox = if (iaprox1 * iaprox2) != 0.0 {
            1.0 / (1.0 / iaprox1 + 1.0 / iaprox2)
        } else {
            0.5 * (iaprox1 + iaprox2)
        };

        let dvdi0 = rs + vt1 / (iaprox + is1) + vt2 / (iaprox + is2);
        let v0 = rs * iaprox + vt1 * (iaprox / is1 + 1.0).ln() + vt2 * (iaprox / is2 + 1.0).ln();
        let il = (iaprox + (v - v0) / dvdi0).max(-is1) * 0.99999;
        let gl = 1.0 / (rs + vt1 / (il + is1) + vt2 / (il + is2));
        (il, gl)
    } else {
        let gl = gmin;
        let il = gl * v - is1;
        (il, gl)
    }
}

/// Diode function from hfetload.c.
fn diode_fn(u: f64) -> f64 {
    const U0: f64 = -2.303;
    const A: f64 = 2.221;
    const B: f64 = 6.804;
    const C: f64 = 1.685;

    let expu = u.exp();
    let it = if u <= U0 {
        expu * (1.0 - expu)
    } else {
        let b_val = 0.5 * (u - U0);
        u + A * ((U0 - u) / B).exp() - (b_val + (b_val * b_val + 0.25 * C * C).sqrt()).ln()
    };

    let ut = it + it.ln();
    let bv = u - ut;
    let c_val = 1.0 + it;
    it * (1.0 + bv / c_val + 0.5 * bv * bv / (c_val * c_val * c_val))
}

/// Core HFET computation with explicit L and W.
fn hfeta_full(
    model: &HfetModel,
    pre: &HfetPrecomp,
    vgs: f64,
    vds: f64,
    l: f64,
    _w: f64,
) -> (f64, f64, f64, f64, f64, f64, f64) {
    let vt = K_OVER_Q * pre.temp;
    let etavth = model.eta * vt;
    let vl = model.vs / pre.t_mu * l;
    let rt = pre.t_mu.recip() * 0.0 + model.rsi + model.rdi; // RSI + RDI
    let vgt0 = vgs - pre.t_vto;
    let s = ((vgt0 - model.vsigmat) / model.vsigma).exp();
    let sigma = model.sigma0 / (1.0 + s);
    let vgt = vgt0 + sigma * vds;
    let u = 0.5 * vgt / vt - 1.0;
    let t = (model.delta_sqr + u * u).sqrt();
    let vgte = vt * (2.0 + u + t);

    let b = (vgt / etavth).exp();
    let nsm = if model.eta2_given && model.d2_given {
        let nsc = pre.n02 * ((vgt + pre.t_vto - model.vt2) / (model.eta2 * vt)).exp();
        let nsn = 2.0 * pre.n0 * (1.0 + 0.5 * b).ln();
        nsn * nsc / (nsn + nsc)
    } else {
        2.0 * pre.n0 * (1.0 + 0.5 * b).ln()
    };

    if nsm < 1.0e-38 {
        return (0.0, 0.0, 0.0, pre.cf, pre.cf, 0.0, 0.0);
    }

    let c = (nsm / model.nmax).powf(model.gamma);
    let q = (1.0 + c).powf(1.0 / model.gamma);
    let ns = nsm / q;
    let gchi = pre.gchi0 * ns;
    let gch = gchi / (1.0 + gchi * rt);
    let gchim = pre.gchi0 * nsm;
    let h = (1.0 + 2.0 * gchim * model.rsi + vgte * vgte / (vl * vl)).sqrt();
    let p = 1.0 + gchim * model.rsi + h;
    let isatm = gchim * vgte / p;
    let g = (isatm / pre.imax).powf(model.gamma);
    let isat = isatm / (1.0 + g).powf(1.0 / model.gamma);
    let vsate = isat / gch;
    let d = (vds / vsate).powf(model.m);
    let e = (1.0 + d).powf(1.0 / model.m);

    let delidgch = vds * (1.0 + pre.t_lambda * vds) / e;
    let cdrain = gch * delidgch;
    let delidvsate = cdrain * d / vsate / (1.0 + d);
    let delidvds = gch * (1.0 + 2.0 * pre.t_lambda * vds) / e
        - cdrain * (vds / vsate).powf(model.m - 1.0) / (vsate * (1.0 + d));

    let a = 1.0 + gchi * rt;
    let delgchgchi = 1.0 / (a * a);
    let delgchins = pre.gchi0;
    let delnsnsm = ns / nsm * (1.0 - c / (1.0 + c));
    let delvgtevgt = 0.5 * (1.0 + u / t);
    let mut delnsmvgt = pre.n0 / etavth / (1.0 / b + 0.5);
    if model.eta2_given && model.d2_given {
        let nsc = pre.n02 * ((vgt + pre.t_vto - model.vt2) / (model.eta2 * vt)).exp();
        let nsn = 2.0 * pre.n0 * (1.0 + 0.5 * b).ln();
        delnsmvgt =
            nsc * (nsc * delnsmvgt + nsn * nsn / (model.eta2 * vt)) / ((nsc + nsn) * (nsc + nsn));
    }
    let delvsateisat = 1.0 / gch;
    let delisatisatm = isat / isatm * (1.0 - g / (1.0 + g));
    let delisatmvgte = gchim * (p - vgte * vgte / (vl * vl * h)) / (p * p);
    let delvsategch = -vsate / gch;
    let delisatmgchim = vgte * (p - gchim * model.rsi * (1.0 + 1.0 / h)) / (p * p);
    let delvgtvgs = 1.0 - vds * model.sigma0 / model.vsigma * s / ((1.0 + s) * (1.0 + s));

    let p_val = delgchgchi * delgchins * delnsnsm * delnsmvgt;
    let delvsatevgt = delvsateisat
        * delisatisatm
        * (delisatmvgte * delvgtevgt + delisatmgchim * pre.gchi0 * delnsmvgt)
        + delvsategch * p_val;
    let g_val = delidgch * p_val + delidvsate * delvsatevgt;
    let gm = g_val * delvgtvgs;
    let gds = delidvds + g_val * sigma;

    // Capacitance calculations
    let temp_val = model.eta1 * vt;
    let cg1 = 1.0 / (model.d1 / model.epsi + temp_val * (-(vgs - model.vt1) / temp_val).exp());
    let cgc = _w * l * (CHARGE * delnsnsm * delnsmvgt * delvgtvgs + cg1);
    let vdse = vds * (1.0 + (vds / vsate).powf(model.mc)).powf(-1.0 / model.mc);
    let a_cap = {
        let ratio = (vsate - vdse) / (2.0 * vsate - vdse);
        ratio * ratio
    };
    let two_thirds = 2.0 / 3.0;
    let pm = model.p + (1.0 - model.p) * (-vds / vsate).exp();
    let capgs = pre.cf + 2.0 * two_thirds * cgc * (1.0 - a_cap) / (1.0 + pm);
    let a_cap2 = {
        let ratio = vsate / (2.0 * vsate - vdse);
        ratio * ratio
    };
    let capgd = pre.cf + 2.0 * pm * two_thirds * cgc * (1.0 - a_cap2) / (1.0 + pm);

    let gmg = 0.0; // Only used with gatemod != 0
    let gmd = 0.0;

    (cdrain, gm, gds, capgs, capgd, gmg, gmd)
}

/// Stamp HFET companion with known voltages into MNA.
pub fn stamp_hfet_with_voltages(
    comp: &HfetCompanion,
    inst: &HfetInstance,
    vgs: f64,
    vgd: f64,
    matrix: &mut SparseMatrix,
    rhs: &mut [f64],
) {
    let model = &inst.model;
    let sign = model.device_type as i32 as f64;
    let m = 1.0; // HFET doesn't use m multiplier in basic form
    let vds = vgs - vgd;

    let gate_prime = inst.gate_prime_idx.or(inst.gate_idx);
    let drain_prime = inst.drain_prime_idx.or(inst.drain_idx);
    let source_prime = inst.source_prime_idx.or(inst.source_idx);

    let gm = comp.gm;
    let gds = comp.gds;
    let ggs = comp.ggs;
    let ggd = comp.ggd;
    let _gmg = comp.gmg;
    let _gmd = comp.gmd;
    let _cg = comp.cg;
    let cd = comp.cd;
    let cgd_val = comp.cgd_current;
    let cgs_val = comp.cgs_current;

    // Norton current sources
    // For gatemod == 0 (standard), the stamp is simpler:
    let ceqgd = sign * (cgd_val - ggd * vgd);
    let ceqgs = sign * (cgs_val - ggs * vgs);
    let cdreq = sign * ((cd + cgd_val) - gds * vds - gm * vgs);

    // RHS
    if let Some(gp) = gate_prime {
        rhs[gp] += m * (-ceqgs - ceqgd);
    }
    if let Some(dp) = drain_prime {
        rhs[dp] += m * (-cdreq + ceqgd);
    }
    if let Some(sp) = source_prime {
        rhs[sp] += m * (cdreq + ceqgs);
    }

    // Series drain resistance
    if inst.drain_prime_idx.is_some() {
        stamp_conductance(
            matrix,
            inst.drain_idx,
            inst.drain_prime_idx,
            model.drain_conduct * m,
        );
    }

    // Series source resistance
    if inst.source_prime_idx.is_some() {
        stamp_conductance(
            matrix,
            inst.source_idx,
            inst.source_prime_idx,
            model.source_conduct * m,
        );
    }

    // Gate resistance
    if inst.gate_prime_idx.is_some() {
        stamp_conductance(
            matrix,
            inst.gate_idx,
            inst.gate_prime_idx,
            model.gate_conduct * m,
        );
    }

    // Internal feedback resistances (ri, rf)
    if inst.drain_prm_prm_idx.is_some() {
        stamp_conductance(
            matrix,
            inst.drain_prime_idx,
            inst.drain_prm_prm_idx,
            model.gf * m,
        );
    }
    if inst.source_prm_prm_idx.is_some() {
        stamp_conductance(
            matrix,
            inst.source_prime_idx,
            inst.source_prm_prm_idx,
            model.gi * m,
        );
    }

    // Y-matrix entries for channel + junction conductances
    // Gate prime diagonal: ggs + ggd
    if let Some(gp) = gate_prime {
        matrix.add(gp, gp, m * (ggs + ggd));
    }
    // Drain prime diagonal: gds + ggd
    if let Some(dp) = drain_prime {
        matrix.add(dp, dp, m * (gds + ggd));
    }
    // Source prime diagonal: gds + gm + ggs
    if let Some(sp) = source_prime {
        matrix.add(sp, sp, m * (gds + gm + ggs));
    }

    // Off-diagonals
    if let (Some(gp), Some(dp)) = (gate_prime, drain_prime) {
        matrix.add(gp, dp, -m * ggd);
        matrix.add(dp, gp, m * (gm - ggd));
    }
    if let (Some(gp), Some(sp)) = (gate_prime, source_prime) {
        matrix.add(gp, sp, -m * ggs);
        matrix.add(sp, gp, m * (-ggs - gm));
    }
    if let (Some(dp), Some(sp)) = (drain_prime, source_prime) {
        matrix.add(dp, sp, m * (-gds - gm));
        matrix.add(sp, dp, -m * gds);
    }
}

/// Full HFET companion computation using instance dimensions.
pub fn hfet_companion_full(inst: &HfetInstance, vgs: f64, vgd: f64, gmin: f64) -> HfetCompanion {
    let model = &inst.model;
    let pre = &inst.precomp;
    let vt = K_OVER_Q * pre.temp;

    let mut cgs_current = 0.0;
    let mut cgd_current = 0.0;
    let mut ggs = 0.0;
    let mut ggd = 0.0;
    let mut gmg = 0.0;
    let mut gmd = 0.0;

    // Gate leakage
    if model.gatemod == 0 {
        if pre.is1s != 0.0 && pre.is2s != 0.0 {
            let (il, gl) = leak(
                gmin, vt, vgs, model.rgs, pre.is1s, pre.is2s, model.m1s, model.m2s,
            );
            cgs_current = il;
            ggs = gl;
        }
        let arg = -vgs * model.del / vt;
        let earg = arg.exp();
        cgs_current += pre.ggrwl * vgs * earg;
        ggs += pre.ggrwl * earg * (1.0 - arg);

        if pre.is1d != 0.0 && pre.is2d != 0.0 {
            let (il, gl) = leak(
                gmin, vt, vgd, model.rgd, pre.is1d, pre.is2d, model.m1d, model.m2d,
            );
            cgd_current = il;
            ggd = gl;
        }
        let arg = -vgd * model.del / vt;
        let earg = arg.exp();
        cgd_current += pre.ggrwl * vgd * earg;
        ggd += pre.ggrwl * earg * (1.0 - arg);
    }

    // Channel current
    let vds = vgs - vgd;
    let (vgs_ch, vds_ch, inverse) = if vds >= 0.0 {
        (vgs, vds, false)
    } else {
        (vgd, -vds, true)
    };

    let (mut cdrain, gm_ch, gds_ch, mut capgs, mut capgd, gmg_ch, gmd_ch) =
        hfeta_full(model, pre, vgs_ch, vds_ch, inst.l, inst.w);

    // Gate-source current for gatemod != 0
    if model.gatemod != 0 {
        let vtn = vt * model.m2s;
        let csat = pre.iso * pre.temp * pre.temp * (-model.phib / (BOLTZMANN * pre.temp)).exp();
        if vgs <= -5.0 * vt {
            ggs = -csat / vgs + gmin;
            cgs_current = ggs * vgs;
        } else {
            let evgs = (vgs / vtn).exp();
            ggs = csat * evgs / vtn + gmin;
            cgs_current = csat * (evgs - 1.0) + gmin * vgs;
        }
        gmg = gmg_ch;
        gmd = gmd_ch;
    }

    let cg = cgs_current + cgd_current;

    if inverse {
        cdrain = -cdrain;
        std::mem::swap(&mut capgs, &mut capgd);
    }

    let cd = cdrain - cgd_current;

    HfetCompanion {
        gm: gm_ch,
        gds: gds_ch,
        ggs,
        ggd,
        cg,
        cd,
        cgd_current,
        cgs_current,
        capgs,
        capgd,
        gmg,
        gmd,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hfet_model_defaults() {
        let m = HfetModel::new(HfetType::Nhfet);
        assert_eq!(m.eta, 1.28);
        assert_eq!(m.vs, 1.5e5);
        assert_eq!(m.nmax, 2e16);
        assert_eq!(m.delta, 3.0);
        assert_eq!(m.deltad, 4.5e-9);
    }

    #[test]
    fn hfet_model_parse() {
        let model_def = ModelDef {
            name: "hfet".to_string(),
            kind: "NHFET".to_string(),
            params: vec![
                thevenin_types::Param {
                    name: "VTO".to_string(),
                    value: Expr::Num(0.3),
                },
                thevenin_types::Param {
                    name: "MU".to_string(),
                    value: Expr::Num(0.385),
                },
                thevenin_types::Param {
                    name: "NMAX".to_string(),
                    value: Expr::Num(6e15),
                },
                thevenin_types::Param {
                    name: "LAMBDA".to_string(),
                    value: Expr::Num(0.17),
                },
            ],
        };
        let m = HfetModel::from_model_def(&model_def);
        assert_eq!(m.device_type, HfetType::Nhfet);
        // vt0 is signed by device type in precompute
        assert!((m.vt0 - 0.3).abs() < 1e-10);
        assert_eq!(m.mu, 0.385);
        assert_eq!(m.nmax, 6e15);
    }

    #[test]
    fn hfet_companion_zero_bias() {
        let mut model = HfetModel::new(HfetType::Nhfet);
        model.vt0 = 0.13;
        model.mu = 0.385;
        model.vs = 1.5e5;
        model.eta = 1.32;
        model.sigma0 = 0.04;
        model.vsigma = 0.1;
        model.vsigmat = 0.3;
        model.nmax = 6e15;
        model.lambda = 0.17;
        model.m = 2.57;
        model.precompute();

        let pre = HfetPrecomp::compute(&model, 300.15, 300.15, 10e-6, 1e-6);
        let inst = HfetInstance {
            name: "z1".to_string(),
            drain_idx: Some(0),
            gate_idx: Some(1),
            source_idx: Some(2),
            gate_prime_idx: None,
            drain_prime_idx: None,
            source_prime_idx: None,
            drain_prm_prm_idx: None,
            source_prm_prm_idx: None,
            model,
            precomp: pre,
            w: 10e-6,
            l: 1e-6,
        };
        let comp = hfet_companion_full(&inst, 0.0, 0.0, 1e-12);
        // At zero bias, drain current should be 0
        // gm and gds may be small but cdrain should be ~0
        assert!(
            comp.cd.abs() < 1e-3,
            "cd at zero bias should be small, got {}",
            comp.cd
        );
    }

    #[test]
    fn hfet_precomp_values() {
        let mut model = HfetModel::new(HfetType::Nhfet);
        model.vt0 = 0.13;
        model.mu = 0.385;
        model.nmax = 6e15;
        model.epsi = 12.244 * 8.85418e-12;
        model.eta = 1.32;
        model.di = 0.04e-6;
        model.vs = 1.5e5;
        model.precompute();

        let pre = HfetPrecomp::compute(&model, 300.15, 300.15, 10e-6, 1e-6);
        assert!(pre.gchi0 > 0.0, "gchi0 should be positive");
        assert!(pre.n0 > 0.0, "n0 should be positive");
        assert!(pre.imax > 0.0, "imax should be positive");
        assert!(pre.cf > 0.0, "cf should be positive");
    }
}
