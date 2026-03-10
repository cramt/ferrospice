//! BSIM4 MOSFET device model.
//!
//! Implements the BSIM4 MOSFET model matching ngspice level 14.
//! This is a physics-based short-channel MOSFET model with advanced
//! features including gate tunneling current, GIDL/GISL, and
//! improved capacitance models.

use std::collections::BTreeMap;

use thevenin_types::{Expr, ModelDef};

use crate::mosfet::MosfetType;
use crate::physics::{
    CHARGE_Q, EPSSI, EXP_THRESHOLD, KBOQ, MAX_EXP, MIN_EXP, bsim_safe_exp as safe_exp,
};

const EPS0: f64 = 8.85418e-12;
const TEMP_DEFAULT: f64 = 300.15;
const DELTA_1: f64 = 1.0e-9;

/// BSIM4 model parameters (from .model card).
#[derive(Debug, Clone)]
#[expect(clippy::struct_excessive_bools)]
pub struct Bsim4Model {
    pub mos_type: MosfetType,

    // Mode selection
    pub mob_mod: i32,
    pub cap_mod: i32,
    pub diomod: i32,
    pub rdsmod: i32,
    pub rbodymod: i32,
    pub rgatemod: i32,
    pub permod: i32,
    pub geomod: i32,
    pub igcmod: i32,
    pub igbmod: i32,
    pub tempmod: i32,
    pub fnoimod: i32,
    pub tnoimod: i32,
    pub trnqsmod: i32,
    pub acnqsmod: i32,
    pub wpemod: i32,
    pub gidlmod: i32,
    pub bin_unit: i32,

    // Oxide
    pub toxe: f64,
    pub toxp: f64,
    pub toxm: f64,
    pub dtox: f64,
    pub epsrox: f64,
    pub toxref: f64,

    // Threshold voltage
    pub vth0: f64,
    pub k1: f64,
    pub k2: f64,
    pub k3: f64,
    pub k3b: f64,
    pub w0: f64,
    pub lpe0: f64,
    pub lpeb: f64,
    pub dvtp0: f64,
    pub dvtp1: f64,
    pub dvt0: f64,
    pub dvt1: f64,
    pub dvt2: f64,
    pub dvt0w: f64,
    pub dvt1w: f64,
    pub dvt2w: f64,
    pub dsub: f64,
    pub minv: f64,
    pub voffl: f64,
    pub phin: f64,
    pub vfb: f64,

    // Subthreshold
    pub voff: f64,
    pub nfactor: f64,
    pub cdsc: f64,
    pub cdscb: f64,
    pub cdscd: f64,
    pub cit: f64,
    pub eta0: f64,
    pub etab: f64,

    // Doping
    pub ndep: f64,
    pub nsd: f64,
    pub ngate: f64,
    pub nsub: f64,
    pub xj: f64,

    // Mobility
    pub u0: f64,
    pub ua: f64,
    pub ub: f64,
    pub uc: f64,
    pub ud: f64,
    pub ute: f64,
    pub ucste: f64,
    pub ua1: f64,
    pub ub1: f64,
    pub uc1: f64,
    pub ud1: f64,
    pub up: f64,
    pub lp: f64,
    pub eu: f64,

    // Velocity saturation
    pub vsat: f64,
    pub a0: f64,
    pub ags: f64,
    pub a1: f64,
    pub a2: f64,
    pub at: f64,
    pub keta: f64,

    // Output resistance
    pub pclm: f64,
    pub pdibl1: f64,
    pub pdibl2: f64,
    pub pdiblb: f64,
    pub pscbe1: f64,
    pub pscbe2: f64,
    pub pvag: f64,
    pub delta: f64,
    pub fprout: f64,
    pub pdits: f64,
    pub pditsd: f64,
    pub pditsl: f64,
    pub drout: f64,

    // Series resistance
    pub rdsw: f64,
    pub rdw: f64,
    pub rsw: f64,
    pub rdswmin: f64,
    pub rdwmin: f64,
    pub rswmin: f64,
    pub prwg: f64,
    pub prwb: f64,
    pub prt: f64,
    pub wr: f64,
    pub rsh: f64,

    // Width/length effects
    pub dwg: f64,
    pub dwb: f64,
    pub b0: f64,
    pub b1: f64,

    // Substrate current
    pub alpha0: f64,
    pub alpha1: f64,
    pub beta0: f64,

    // GIDL/GISL
    pub agidl: f64,
    pub bgidl: f64,
    pub cgidl: f64,
    pub egidl: f64,
    pub agisl: f64,
    pub bgisl: f64,
    pub cgisl: f64,
    pub egisl: f64,

    // Gate tunneling
    pub aigc: f64,
    pub bigc: f64,
    pub cigc: f64,
    pub aigsd: f64,
    pub bigsd: f64,
    pub cigsd: f64,
    pub nigc: f64,
    pub pigcd: f64,
    pub poxedge: f64,
    pub ntox: f64,
    pub aigbacc: f64,
    pub bigbacc: f64,
    pub cigbacc: f64,
    pub nigbacc: f64,
    pub aigbinv: f64,
    pub bigbinv: f64,
    pub cigbinv: f64,
    pub eigbinv: f64,
    pub nigbinv: f64,

    // Capacitance
    pub cgso: f64,
    pub cgdo: f64,
    pub cgbo: f64,
    pub cgsl: f64,
    pub cgdl: f64,
    pub ckappas: f64,
    pub ckappad: f64,
    pub cf: f64,
    pub clc: f64,
    pub cle: f64,
    pub xpart: f64,
    pub acde: f64,
    pub moin: f64,
    pub noff: f64,
    pub voffcv: f64,
    pub vfbcv: f64,

    // Gate/body resistance
    pub xrcrg1: f64,
    pub xrcrg2: f64,
    pub rbpb: f64,
    pub rbpd: f64,
    pub rbps: f64,
    pub rbdb: f64,
    pub rbsb: f64,
    pub gbmin: f64,
    pub rshg: f64,
    pub ngcon: f64,

    // Junction diode (source side)
    pub jss: f64,
    pub jsws: f64,
    pub jswgs: f64,
    pub njs: f64,
    pub ijthsfwd: f64,
    pub ijthsrev: f64,
    pub bvs: f64,
    pub xjbvs: f64,
    pub xtis: f64,

    // Junction diode (drain side)
    pub jsd: f64,
    pub jswd: f64,
    pub jswgd: f64,
    pub njd: f64,
    pub ijthdfwd: f64,
    pub ijthdrev: f64,
    pub bvd: f64,
    pub xjbvd: f64,
    pub xtid: f64,

    // Junction capacitance (source)
    pub pbs: f64,
    pub cjs: f64,
    pub mjs: f64,
    pub pbsws: f64,
    pub cjsws: f64,
    pub mjsws: f64,
    pub pbswgs: f64,
    pub cjswgs: f64,
    pub mjswgs: f64,

    // Junction capacitance (drain)
    pub pbd: f64,
    pub cjd: f64,
    pub mjd: f64,
    pub pbswd: f64,
    pub cjswd: f64,
    pub mjswd: f64,
    pub pbswgd: f64,
    pub cjswgd: f64,
    pub mjswgd: f64,

    // Junction temperature coefficients
    pub tpb: f64,
    pub tcj: f64,
    pub tpbsw: f64,
    pub tcjsw: f64,
    pub tpbswg: f64,
    pub tcjswg: f64,

    // Temperature
    pub tnom: f64,
    pub kt1: f64,
    pub kt1l: f64,
    pub kt2: f64,

    // Geometry
    pub dmcg: f64,
    pub dmci: f64,
    pub dmdg: f64,
    pub dmcgt: f64,
    pub dwj: f64,
    pub xgw: f64,
    pub xgl: f64,

    // Binning dimensions
    pub lint: f64,
    pub wint: f64,
    pub ll: f64,
    pub wl: f64,
    pub lln: f64,
    pub wln: f64,
    pub lw: f64,
    pub ww: f64,
    pub lwn: f64,
    pub wwn: f64,
    pub lwl: f64,
    pub wwl: f64,
    pub xl: f64,
    pub xw: f64,
    pub dlc: f64,
    pub dwc: f64,

    // Noise
    pub kf: f64,
    pub af: f64,
    pub ef: f64,
    pub noia: f64,
    pub noib: f64,
    pub noic: f64,
    pub em: f64,

    // "Given" flags
    pub vth0_given: bool,
    pub vfb_given: bool,
    pub k1_given: bool,
    pub k2_given: bool,
    pub dlc_given: bool,
    pub dwc_given: bool,
    pub cf_given: bool,
    pub dsub_given: bool,
    pub toxp_given: bool,
    pub toxm_given: bool,

    // L/W/P binning parameters (stored compactly)
    pub l_params: BTreeMap<String, f64>,
    pub w_params: BTreeMap<String, f64>,
    pub p_params: BTreeMap<String, f64>,
}

/// Size-dependent parameters computed from model + instance geometry.
#[derive(Debug, Clone)]
pub struct Bsim4SizeParam {
    // Effective dimensions
    pub leff: f64,
    pub weff: f64,
    pub leff_cv: f64,
    pub weff_cv: f64,
    pub weff_cj: f64,

    // Binned core parameters
    pub vth0: f64,
    pub k1: f64,
    pub k2: f64,
    pub k3: f64,
    pub k3b: f64,
    pub w0: f64,
    pub lpe0: f64,
    pub lpeb: f64,
    pub dvtp0: f64,
    pub dvtp1: f64,
    pub dvt0: f64,
    pub dvt1: f64,
    pub dvt2: f64,
    pub dvt0w: f64,
    pub dvt1w: f64,
    pub dvt2w: f64,
    pub dsub: f64,
    pub minv: f64,
    pub eta0: f64,
    pub etab: f64,
    pub voff: f64,
    pub nfactor: f64,
    pub cdsc: f64,
    pub cdscb: f64,
    pub cdscd: f64,
    pub cit: f64,
    pub phin: f64,

    pub u0: f64,
    pub ua: f64,
    pub ub: f64,
    pub uc: f64,
    pub ud: f64,
    pub vsat: f64,
    pub a0: f64,
    pub ags: f64,
    pub a1: f64,
    pub a2: f64,
    pub keta: f64,
    pub pclm: f64,
    pub pdibl1: f64,
    pub pdibl2: f64,
    pub pdiblb: f64,
    pub pscbe1: f64,
    pub pscbe2: f64,
    pub pvag: f64,
    pub delta: f64,
    pub fprout: f64,
    pub pdits: f64,
    pub pditsd: f64,
    pub drout: f64,
    pub rdsw: f64,
    pub prwg: f64,
    pub prwb: f64,
    pub dwg: f64,
    pub dwb: f64,
    pub b0: f64,
    pub b1: f64,
    pub alpha0: f64,
    pub alpha1: f64,
    pub beta0: f64,
    pub agidl: f64,
    pub bgidl: f64,
    pub cgidl: f64,
    pub egidl: f64,
    pub agisl: f64,
    pub bgisl: f64,
    pub cgisl: f64,
    pub egisl: f64,

    pub aigc: f64,
    pub bigc: f64,
    pub cigc: f64,
    pub aigsd: f64,
    pub bigsd: f64,
    pub cigsd: f64,
    pub nigc: f64,
    pub pigcd: f64,
    pub poxedge: f64,
    pub ntox: f64,
    pub aigbacc: f64,
    pub bigbacc: f64,
    pub cigbacc: f64,
    pub nigbacc: f64,
    pub aigbinv: f64,
    pub bigbinv: f64,
    pub cigbinv: f64,
    pub eigbinv: f64,
    pub nigbinv: f64,

    pub xrcrg1: f64,
    pub xrcrg2: f64,

    pub cgsl: f64,
    pub cgdl: f64,
    pub ckappas: f64,
    pub ckappad: f64,
    pub cf: f64,
    pub clc: f64,
    pub cle: f64,
    pub acde: f64,
    pub moin: f64,
    pub noff: f64,
    pub voffcv: f64,
    pub vfbcv: f64,

    pub kt1: f64,
    pub kt1l: f64,
    pub kt2: f64,

    pub ndep: f64,
    pub ngate: f64,
    pub nsub: f64,
    pub nsd: f64,
    pub xj: f64,
    pub vfb: f64,

    // Temperature/size derived quantities
    pub u0temp: f64,
    pub vsattemp: f64,
    pub rds0: f64,
    pub rdswmin: f64,
    pub phi: f64,
    pub sqrt_phi: f64,
    pub phis3: f64,
    pub xdep0: f64,
    pub sqrt_xdep0: f64,
    pub litl: f64,
    pub vbi: f64,
    pub cdep0: f64,
    pub k1ox: f64,
    pub k2ox: f64,
    pub theta0vb0: f64,
    pub theta_rout: f64,
    pub vfbzb_factor: f64,

    // Oxide capacitance
    pub coxe: f64,
    pub coxp: f64,
    pub factor1: f64,

    // Overlap capacitances (per instance)
    pub cgso_eff: f64,
    pub cgdo_eff: f64,
    pub cgbo_eff: f64,

    // Subthreshold
    pub mstar: f64,
    pub voffcbn: f64,
    pub ldeb: f64,

    // Gate tunneling
    pub aechvb: f64,
    pub bechvb: f64,
    pub aechvb_edge_s: f64,
    pub aechvb_edge_d: f64,
    pub bechvb_edge: f64,
    pub tox_ratio: f64,
    pub tox_ratio_edge: f64,

    // Abulk factor
    pub abulk_cv_factor: f64,

    // Mobility model derived
    pub vtfbphi1: f64,
}

/// NR companion result for BSIM4.
#[derive(Debug, Clone)]
pub struct Bsim4Companion {
    pub ids: f64,
    pub gm: f64,
    pub gds: f64,
    pub gmbs: f64,
    pub gbd: f64,
    pub gbs: f64,
    pub cbd_current: f64,
    pub cbs_current: f64,
    pub isub: f64,
    pub igidl: f64,
    pub igisl: f64,
    pub ceq_d: f64,
    pub ceq_bs: f64,
    pub ceq_bd: f64,
    pub ceq_sub: f64,
    pub mode: i32,
    pub vdsat: f64,
    // Capacitances (9 intrinsic + junction)
    pub cggb: f64,
    pub cgdb: f64,
    pub cgsb: f64,
    pub cbgb: f64,
    pub cbdb: f64,
    pub cbsb: f64,
    pub cdgb: f64,
    pub cddb: f64,
    pub cdsb: f64,
    pub capbd: f64,
    pub capbs: f64,
    pub qinv: f64,
    // Intermediate values needed for noise computation (fnoimod=1)
    pub ueff: f64,
    pub vgsteff: f64,
    pub vdseff: f64,
    pub abulk: f64,
}

/// Resolved BSIM4 instance with node indices.
#[derive(Debug, Clone)]
pub struct Bsim4Instance {
    pub name: String,
    pub drain_idx: Option<usize>,
    pub gate_idx: Option<usize>,
    pub source_idx: Option<usize>,
    pub bulk_idx: Option<usize>,
    pub drain_prime_idx: Option<usize>,
    pub source_prime_idx: Option<usize>,
    pub w: f64,
    pub l: f64,
    pub nf: f64,
    pub ad: f64,
    pub as_: f64,
    pub pd: f64,
    pub ps: f64,
    pub nrd: f64,
    pub nrs: f64,
    pub m: f64,
    pub sa: f64,
    pub sb: f64,
    pub model: Bsim4Model,
    pub size_params: Bsim4SizeParam,
}

impl Bsim4Instance {
    /// Effective drain node: drain_prime if series resistance exists, else drain.
    pub fn drain_eff_idx(&self) -> Option<usize> {
        self.drain_prime_idx.or(self.drain_idx)
    }

    /// Effective source node: source_prime if series resistance exists, else source.
    pub fn source_eff_idx(&self) -> Option<usize> {
        self.source_prime_idx.or(self.source_idx)
    }

    pub fn terminal_voltages(&self, solution: &[f64]) -> (f64, f64, f64) {
        let v_dp = self.drain_eff_idx().map(|i| solution[i]).unwrap_or(0.0);
        let v_g = self.gate_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_sp = self.source_eff_idx().map(|i| solution[i]).unwrap_or(0.0);
        let v_b = self.bulk_idx.map(|i| solution[i]).unwrap_or(0.0);

        let sign = self.model.mos_type.sign();
        (
            sign * (v_g - v_sp),
            sign * (v_dp - v_sp),
            sign * (v_b - v_sp),
        )
    }

    pub fn drain_conductance(&self) -> f64 {
        if self.model.rsh > 0.0 && self.nrd > 0.0 {
            1.0 / (self.model.rsh * self.nrd)
        } else {
            0.0
        }
    }

    pub fn source_conductance(&self) -> f64 {
        if self.model.rsh > 0.0 && self.nrs > 0.0 {
            1.0 / (self.model.rsh * self.nrs)
        } else {
            0.0
        }
    }

    /// Build AC stamp data from this instance and its companion model.
    pub fn ac_stamp(&self, comp: &Bsim4Companion) -> crate::ac::BsimAcStamp {
        crate::ac::BsimAcStamp {
            dp: self.drain_eff_idx(),
            g: self.gate_idx,
            sp: self.source_eff_idx(),
            b: self.bulk_idx,
            drain_idx: self.drain_idx,
            source_idx: self.source_idx,
            m: self.m,
            gds: comp.gds,
            gm: comp.gm,
            gmbs: comp.gmbs,
            gbd: comp.gbd,
            gbs: comp.gbs,
            g_drain: self.drain_conductance(),
            g_source: self.source_conductance(),
            cggb: comp.cggb,
            cgdb: comp.cgdb,
            cgsb: comp.cgsb,
            cbgb: comp.cbgb,
            cbdb: comp.cbdb,
            cbsb: comp.cbsb,
            cdgb: comp.cdgb,
            cddb: comp.cddb,
            cdsb: comp.cdsb,
            capbd: comp.capbd,
            capbs: comp.capbs,
        }
    }
}

impl Bsim4Model {
    pub fn new(mos_type: MosfetType) -> Self {
        let type_sign = mos_type.sign();
        let is_nmos = type_sign > 0.0;
        Self {
            mos_type,
            mob_mod: 0,
            cap_mod: 2,
            diomod: 1,
            rdsmod: 0,
            rbodymod: 0,
            rgatemod: 0,
            permod: 1,
            geomod: 0,
            igcmod: 0,
            igbmod: 0,
            tempmod: 0,
            fnoimod: 1,
            tnoimod: 0,
            trnqsmod: 0,
            acnqsmod: 0,
            wpemod: 0,
            gidlmod: 0,
            bin_unit: 1,

            toxe: 30.0e-10,
            toxp: 30.0e-10,
            toxm: 30.0e-10,
            dtox: 0.0,
            epsrox: 3.9,
            toxref: 30.0e-10,

            vth0: if is_nmos { 0.7 } else { -0.7 },
            k1: 0.53,
            k2: -0.0186,
            k3: 80.0,
            k3b: 0.0,
            w0: 2.5e-6,
            lpe0: 1.74e-7,
            lpeb: 0.0,
            dvtp0: 0.0,
            dvtp1: 0.0,
            dvt0: 2.2,
            dvt1: 0.53,
            dvt2: -0.032,
            dvt0w: 0.0,
            dvt1w: 5.3e6,
            dvt2w: -0.032,
            dsub: 0.56,
            minv: 0.0,
            voffl: 0.0,
            phin: 0.0,
            vfb: -1.0,

            voff: -0.08,
            nfactor: 1.0,
            cdsc: 2.4e-4,
            cdscb: 0.0,
            cdscd: 0.0,
            cit: 0.0,
            eta0: 0.08,
            etab: -0.07,

            ndep: 1.7e17,
            nsd: 1.0e20,
            ngate: 0.0,
            nsub: 6.0e16,
            xj: 0.15e-6,

            u0: if is_nmos { 0.067 } else { 0.025 },
            ua: 1.0e-9,
            ub: 1.0e-19,
            uc: -0.0465e-9,
            ud: 0.0,
            ute: -1.5,
            ucste: -4.775e-3,
            ua1: 1.0e-9,
            ub1: -1.0e-18,
            uc1: -0.056e-9,
            ud1: 0.0,
            up: 0.0,
            lp: 1.0e-8,
            eu: if is_nmos { 1.67 } else { 1.0 },

            vsat: 8.0e4,
            a0: 1.0,
            ags: 0.0,
            a1: 0.0,
            a2: 1.0,
            at: 3.3e4,
            keta: -0.047,

            pclm: 1.3,
            pdibl1: 0.39,
            pdibl2: 0.0086,
            pdiblb: 0.0,
            pscbe1: 4.24e8,
            pscbe2: 1.0e-5,
            pvag: 0.0,
            delta: 0.01,
            fprout: 0.0,
            pdits: 0.0,
            pditsd: 0.0,
            pditsl: 0.0,
            drout: 0.56,

            rdsw: 200.0,
            rdw: 100.0,
            rsw: 100.0,
            rdswmin: 0.0,
            rdwmin: 0.0,
            rswmin: 0.0,
            prwg: 1.0,
            prwb: 0.0,
            prt: 0.0,
            wr: 1.0,
            rsh: 0.0,

            dwg: 0.0,
            dwb: 0.0,
            b0: 0.0,
            b1: 0.0,

            alpha0: 0.0,
            alpha1: 0.0,
            beta0: 0.0,

            agidl: 0.0,
            bgidl: 2.3e9,
            cgidl: 0.5,
            egidl: 0.8,
            agisl: 0.0,
            bgisl: 2.3e9,
            cgisl: 0.5,
            egisl: 0.8,

            aigc: if is_nmos { 1.36e-2 } else { 9.80e-3 },
            bigc: if is_nmos { 1.71e-3 } else { 7.59e-4 },
            cigc: if is_nmos { 0.075 } else { 0.03 },
            aigsd: if is_nmos { 1.36e-2 } else { 9.80e-3 },
            bigsd: if is_nmos { 1.71e-3 } else { 7.59e-4 },
            cigsd: if is_nmos { 0.075 } else { 0.03 },
            nigc: 1.0,
            pigcd: 1.0,
            poxedge: 1.0,
            ntox: 1.0,
            aigbacc: 1.36e-2,
            bigbacc: 1.71e-3,
            cigbacc: 0.075,
            nigbacc: 1.0,
            aigbinv: 1.11e-2,
            bigbinv: 9.49e-4,
            cigbinv: 0.006,
            eigbinv: 1.1,
            nigbinv: 3.0,

            cgso: 0.0,
            cgdo: 0.0,
            cgbo: 0.0,
            cgsl: 0.0,
            cgdl: 0.0,
            ckappas: 0.6,
            ckappad: 0.6,
            cf: 0.0,
            clc: 0.1e-6,
            cle: 0.6,
            xpart: 0.0,
            acde: 1.0,
            moin: 15.0,
            noff: 1.0,
            voffcv: 0.0,
            vfbcv: -1.0,

            xrcrg1: 12.0,
            xrcrg2: 1.0,
            rbpb: 50.0,
            rbpd: 50.0,
            rbps: 50.0,
            rbdb: 50.0,
            rbsb: 50.0,
            gbmin: 1.0e-12,
            rshg: 0.1,
            ngcon: 1.0,

            jss: 1.0e-4,
            jsws: 0.0,
            jswgs: 0.0,
            njs: 1.0,
            ijthsfwd: 0.1,
            ijthsrev: 0.1,
            bvs: 10.0,
            xjbvs: 1.0,
            xtis: 3.0,

            jsd: 1.0e-4,
            jswd: 0.0,
            jswgd: 0.0,
            njd: 1.0,
            ijthdfwd: 0.1,
            ijthdrev: 0.1,
            bvd: 10.0,
            xjbvd: 1.0,
            xtid: 3.0,

            pbs: 1.0,
            cjs: 5.0e-4,
            mjs: 0.5,
            pbsws: 1.0,
            cjsws: 5.0e-10,
            mjsws: 0.33,
            pbswgs: 1.0,
            cjswgs: 0.0,
            mjswgs: 0.33,

            pbd: 1.0,
            cjd: 5.0e-4,
            mjd: 0.5,
            pbswd: 1.0,
            cjswd: 5.0e-10,
            mjswd: 0.33,
            pbswgd: 1.0,
            cjswgd: 0.0,
            mjswgd: 0.33,

            tpb: 0.0,
            tcj: 0.0,
            tpbsw: 0.0,
            tcjsw: 0.0,
            tpbswg: 0.0,
            tcjswg: 0.0,

            tnom: TEMP_DEFAULT,
            kt1: -0.11,
            kt1l: 0.0,
            kt2: 0.022,

            dmcg: 0.0,
            dmci: 0.0,
            dmdg: 0.0,
            dmcgt: 0.0,
            dwj: 0.0,
            xgw: 0.0,
            xgl: 0.0,

            lint: 0.0,
            wint: 0.0,
            ll: 0.0,
            wl: 0.0,
            lln: 1.0,
            wln: 1.0,
            lw: 0.0,
            ww: 0.0,
            lwn: 1.0,
            wwn: 1.0,
            lwl: 0.0,
            wwl: 0.0,
            xl: 0.0,
            xw: 0.0,
            dlc: 0.0,
            dwc: 0.0,

            kf: 0.0,
            af: 1.0,
            ef: 1.0,
            noia: if is_nmos { 6.25e41 } else { 6.188e40 },
            noib: if is_nmos { 3.125e26 } else { 1.5e25 },
            noic: 8.75e9,
            em: 4.1e7,

            vth0_given: false,
            vfb_given: false,
            k1_given: false,
            k2_given: false,
            dlc_given: false,
            dwc_given: false,
            cf_given: false,
            dsub_given: false,
            toxp_given: false,
            toxm_given: false,

            l_params: BTreeMap::new(),
            w_params: BTreeMap::new(),
            p_params: BTreeMap::new(),
        }
    }

    pub fn coxe(&self) -> f64 {
        self.epsrox * EPS0 / self.toxe
    }

    pub fn coxp(&self) -> f64 {
        self.epsrox * EPS0 / self.toxp
    }

    pub fn internal_node_count(&self, nrd: f64, nrs: f64) -> usize {
        let mut count = 0;
        // drain prime if rsh*nrd > 0 or rdsmod=1
        if (self.rsh > 0.0 && nrd > 0.0) || self.rdsmod != 0 {
            count += 1;
        }
        // source prime
        if (self.rsh > 0.0 && nrs > 0.0) || self.rdsmod != 0 {
            count += 1;
        }
        count
    }

    fn get_l(&self, name: &str) -> f64 {
        self.l_params.get(name).copied().unwrap_or(0.0)
    }

    fn get_w(&self, name: &str) -> f64 {
        self.w_params.get(name).copied().unwrap_or(0.0)
    }

    fn get_p(&self, name: &str) -> f64 {
        self.p_params.get(name).copied().unwrap_or(0.0)
    }

    fn bin_param(&self, base: f64, name: &str, inv_l: f64, inv_w: f64, inv_lw: f64) -> f64 {
        base + self.get_l(name) * inv_l + self.get_w(name) * inv_w + self.get_p(name) * inv_lw
    }

    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let mos_type = if model_def.kind.to_uppercase().contains("PMOS") {
            MosfetType::Pmos
        } else {
            MosfetType::Nmos
        };
        let mut model = Self::new(mos_type);

        let mut vth0_given = false;
        let mut vfb_given = false;
        let mut k1_given = false;
        let mut k2_given = false;
        let mut dlc_given = false;
        let mut dwc_given = false;
        let mut cf_given = false;
        let mut dsub_given = false;
        let mut toxp_given = false;
        let mut toxm_given = false;

        for p in &model_def.params {
            let v = if let Expr::Num(v) = &p.value {
                *v
            } else {
                continue;
            };
            let key = p.name.to_lowercase();
            let key_str = key.as_str();

            // Check for L/W/P prefixed binning params
            if let Some(base_name) = key_str.strip_prefix('l')
                && is_binnable_param(base_name)
            {
                model.l_params.insert(base_name.to_string(), v);
                continue;
            }
            if let Some(base_name) = key_str.strip_prefix('w')
                && is_binnable_param(base_name)
            {
                model.w_params.insert(base_name.to_string(), v);
                continue;
            }
            if let Some(base_name) = key_str.strip_prefix('p')
                && is_binnable_param(base_name)
            {
                model.p_params.insert(base_name.to_string(), v);
                continue;
            }

            match key_str {
                // Mode selection
                "mobmod" => model.mob_mod = v as i32,
                "capmod" => model.cap_mod = v as i32,
                "diomod" => model.diomod = v as i32,
                "rdsmod" => model.rdsmod = v as i32,
                "rbodymod" => model.rbodymod = v as i32,
                "rgatemod" => model.rgatemod = v as i32,
                "permod" => model.permod = v as i32,
                "geomod" => model.geomod = v as i32,
                "igcmod" => model.igcmod = v as i32,
                "igbmod" => model.igbmod = v as i32,
                "tempmod" => model.tempmod = v as i32,
                "fnoimod" => model.fnoimod = v as i32,
                "tnoimod" => model.tnoimod = v as i32,
                "trnqsmod" => model.trnqsmod = v as i32,
                "acnqsmod" => model.acnqsmod = v as i32,
                "wpemod" => model.wpemod = v as i32,
                "gidlmod" => model.gidlmod = v as i32,
                "binunit" => model.bin_unit = v as i32,
                "paramchk" => {} // ignored

                // Oxide
                "toxe" => model.toxe = v,
                "toxp" => {
                    model.toxp = v;
                    toxp_given = true;
                }
                "toxm" => {
                    model.toxm = v;
                    toxm_given = true;
                }
                "dtox" => model.dtox = v,
                "epsrox" => model.epsrox = v,
                "toxref" => model.toxref = v,

                // Threshold voltage
                "vth0" | "vtho" => {
                    model.vth0 = v;
                    vth0_given = true;
                }
                "k1" => {
                    model.k1 = v;
                    k1_given = true;
                }
                "k2" => {
                    model.k2 = v;
                    k2_given = true;
                }
                "k3" => model.k3 = v,
                "k3b" => model.k3b = v,
                "w0" => model.w0 = v,
                "lpe0" => model.lpe0 = v,
                "lpeb" => model.lpeb = v,
                "dvtp0" => model.dvtp0 = v,
                "dvtp1" => model.dvtp1 = v,
                "dvt0" => model.dvt0 = v,
                "dvt1" => model.dvt1 = v,
                "dvt2" => model.dvt2 = v,
                "dvt0w" => model.dvt0w = v,
                "dvt1w" => model.dvt1w = v,
                "dvt2w" => model.dvt2w = v,
                "dsub" => {
                    model.dsub = v;
                    dsub_given = true;
                }
                "minv" => model.minv = v,
                "voffl" => model.voffl = v,
                "phin" => model.phin = v,
                "vfb" => {
                    model.vfb = v;
                    vfb_given = true;
                }

                // Subthreshold
                "voff" => model.voff = v,
                "nfactor" => model.nfactor = v,
                "cdsc" => model.cdsc = v,
                "cdscb" => model.cdscb = v,
                "cdscd" => model.cdscd = v,
                "cit" => model.cit = v,
                "eta0" => model.eta0 = v,
                "etab" => model.etab = v,

                // Doping
                "ndep" => model.ndep = v,
                "nsd" => model.nsd = v,
                "ngate" => model.ngate = v,
                "nsub" => model.nsub = v,
                "xj" => model.xj = v,

                // Mobility
                "u0" => model.u0 = v,
                "ua" => model.ua = v,
                "ub" => model.ub = v,
                "uc" => model.uc = v,
                "ud" => model.ud = v,
                "ute" => model.ute = v,
                "ucste" => model.ucste = v,
                "ua1" => model.ua1 = v,
                "ub1" => model.ub1 = v,
                "uc1" => model.uc1 = v,
                "ud1" => model.ud1 = v,
                "up" => model.up = v,
                "lp" => model.lp = v,
                "eu" => model.eu = v,

                // Velocity saturation
                "vsat" => model.vsat = v,
                "a0" => model.a0 = v,
                "ags" => model.ags = v,
                "a1" => model.a1 = v,
                "a2" => model.a2 = v,
                "at" => model.at = v,
                "keta" => model.keta = v,

                // Output resistance
                "pclm" => model.pclm = v,
                "pdiblc1" | "pdibl1" => model.pdibl1 = v,
                "pdiblc2" | "pdibl2" => model.pdibl2 = v,
                "pdiblcb" | "pdiblb" => model.pdiblb = v,
                "pscbe1" => model.pscbe1 = v,
                "pscbe2" => model.pscbe2 = v,
                "pvag" => model.pvag = v,
                "delta" => model.delta = v,
                "fprout" => model.fprout = v,
                "pdits" => model.pdits = v,
                "pditsd" => model.pditsd = v,
                "pditsl" => model.pditsl = v,
                "drout" => model.drout = v,

                // Series resistance
                "rdsw" => model.rdsw = v,
                "rdw" => model.rdw = v,
                "rsw" => model.rsw = v,
                "rdswmin" => model.rdswmin = v,
                "rdwmin" => model.rdwmin = v,
                "rswmin" => model.rswmin = v,
                "prwg" => model.prwg = v,
                "prwb" => model.prwb = v,
                "prt" => model.prt = v,
                "wr" => model.wr = v,
                "rsh" => model.rsh = v,

                // Width/length
                "dwg" => model.dwg = v,
                "dwb" => model.dwb = v,
                "b0" => model.b0 = v,
                "b1" => model.b1 = v,

                // Substrate current
                "alpha0" => model.alpha0 = v,
                "alpha1" => model.alpha1 = v,
                "beta0" => model.beta0 = v,

                // GIDL/GISL
                "agidl" => model.agidl = v,
                "bgidl" => model.bgidl = v,
                "cgidl" => model.cgidl = v,
                "egidl" => model.egidl = v,
                "agisl" => model.agisl = v,
                "bgisl" => model.bgisl = v,
                "cgisl" => model.cgisl = v,
                "egisl" => model.egisl = v,

                // Gate tunneling
                "aigc" => model.aigc = v,
                "bigc" => model.bigc = v,
                "cigc" => model.cigc = v,
                "aigsd" => model.aigsd = v,
                "bigsd" => model.bigsd = v,
                "cigsd" => model.cigsd = v,
                "nigc" => model.nigc = v,
                "pigcd" => model.pigcd = v,
                "poxedge" => model.poxedge = v,
                "ntox" => model.ntox = v,
                "aigbacc" => model.aigbacc = v,
                "bigbacc" => model.bigbacc = v,
                "cigbacc" => model.cigbacc = v,
                "nigbacc" => model.nigbacc = v,
                "aigbinv" => model.aigbinv = v,
                "bigbinv" => model.bigbinv = v,
                "cigbinv" => model.cigbinv = v,
                "eigbinv" => model.eigbinv = v,
                "nigbinv" => model.nigbinv = v,

                // Capacitance
                "cgso" => model.cgso = v,
                "cgdo" => model.cgdo = v,
                "cgbo" => model.cgbo = v,
                "cgsl" => model.cgsl = v,
                "cgdl" => model.cgdl = v,
                "ckappas" => model.ckappas = v,
                "ckappad" => model.ckappad = v,
                "cf" => {
                    model.cf = v;
                    cf_given = true;
                }
                "clc" => model.clc = v,
                "cle" => model.cle = v,
                "xpart" => model.xpart = v,
                "acde" => model.acde = v,
                "moin" => model.moin = v,
                "noff" => model.noff = v,
                "voffcv" => model.voffcv = v,
                "vfbcv" => model.vfbcv = v,

                // Gate/body resistance
                "xrcrg1" => model.xrcrg1 = v,
                "xrcrg2" => model.xrcrg2 = v,
                "rbpb" => model.rbpb = v,
                "rbpd" => model.rbpd = v,
                "rbps" => model.rbps = v,
                "rbdb" => model.rbdb = v,
                "rbsb" => model.rbsb = v,
                "gbmin" => model.gbmin = v,
                "rshg" => model.rshg = v,
                "ngcon" => model.ngcon = v,

                // Junction source
                "jss" => model.jss = v,
                "jsws" => model.jsws = v,
                "jswgs" => model.jswgs = v,
                "njs" => model.njs = v,
                "ijthsfwd" => model.ijthsfwd = v,
                "ijthsrev" => model.ijthsrev = v,
                "bvs" => model.bvs = v,
                "xjbvs" => model.xjbvs = v,
                "xtis" => model.xtis = v,

                // Junction drain
                "jsd" => model.jsd = v,
                "jswd" => model.jswd = v,
                "jswgd" => model.jswgd = v,
                "njd" => model.njd = v,
                "ijthdfwd" => model.ijthdfwd = v,
                "ijthdrev" => model.ijthdrev = v,
                "bvd" => model.bvd = v,
                "xjbvd" => model.xjbvd = v,
                "xtid" => model.xtid = v,

                // Junction caps source
                "pbs" => model.pbs = v,
                "cjs" => model.cjs = v,
                "mjs" => model.mjs = v,
                "pbsws" => model.pbsws = v,
                "cjsws" => model.cjsws = v,
                "mjsws" => model.mjsws = v,
                "pbswgs" => model.pbswgs = v,
                "cjswgs" => model.cjswgs = v,
                "mjswgs" => model.mjswgs = v,

                // Junction caps drain
                "pbd" => model.pbd = v,
                "cjd" => model.cjd = v,
                "mjd" => model.mjd = v,
                "pbswd" => model.pbswd = v,
                "cjswd" => model.cjswd = v,
                "mjswd" => model.mjswd = v,
                "pbswgd" => model.pbswgd = v,
                "cjswgd" => model.cjswgd = v,
                "mjswgd" => model.mjswgd = v,

                // Junction temp coefficients
                "tpb" => model.tpb = v,
                "tcj" => model.tcj = v,
                "tpbsw" => model.tpbsw = v,
                "tcjsw" => model.tcjsw = v,
                "tpbswg" => model.tpbswg = v,
                "tcjswg" => model.tcjswg = v,

                // Temperature
                "tnom" => model.tnom = v + 273.15,
                "kt1" => model.kt1 = v,
                "kt1l" => model.kt1l = v,
                "kt2" => model.kt2 = v,

                // Geometry
                "dmcg" => model.dmcg = v,
                "dmci" => model.dmci = v,
                "dmdg" => model.dmdg = v,
                "dmcgt" => model.dmcgt = v,
                "dwj" => model.dwj = v,
                "xgw" => model.xgw = v,
                "xgl" => model.xgl = v,

                // Binning
                "lint" => model.lint = v,
                "wint" => model.wint = v,
                "ll" => model.ll = v,
                "wl" => model.wl = v,
                "lln" => model.lln = v,
                "wln" => model.wln = v,
                "lw" => model.lw = v,
                "ww" => model.ww = v,
                "lwn" => model.lwn = v,
                "wwn" => model.wwn = v,
                "lwl" => model.lwl = v,
                "wwl" => model.wwl = v,
                "xl" => model.xl = v,
                "xw" => model.xw = v,
                "dlc" => {
                    model.dlc = v;
                    dlc_given = true;
                }
                "dwc" => {
                    model.dwc = v;
                    dwc_given = true;
                }

                // Noise
                "kf" => model.kf = v,
                "af" => model.af = v,
                "ef" => model.ef = v,
                "noia" => model.noia = v,
                "noib" => model.noib = v,
                "noic" => model.noic = v,
                "em" => model.em = v,

                "level" | "version" => {} // handled elsewhere
                _ => {}                   // ignore unknown params
            }
        }

        // Post-process conditional defaults
        if !toxp_given {
            model.toxp = model.toxe;
        }
        if !toxm_given {
            model.toxm = model.toxe;
        }
        if !dsub_given {
            model.dsub = model.drout;
        }
        if !dlc_given {
            model.dlc = model.lint;
        }
        if !dwc_given {
            model.dwc = model.wint;
        }
        if !cf_given {
            model.cf =
                2.0 * model.epsrox * EPS0 / std::f64::consts::PI * (1.0 + 0.4e-6 / model.toxe).ln();
        }

        model.vth0_given = vth0_given;
        model.vfb_given = vfb_given;
        model.k1_given = k1_given;
        model.k2_given = k2_given;
        model.dlc_given = dlc_given;
        model.dwc_given = dwc_given;
        model.cf_given = cf_given;
        model.dsub_given = dsub_given;
        model.toxp_given = toxp_given;
        model.toxm_given = toxm_given;

        model
    }

    #[expect(clippy::too_many_lines)]
    pub fn size_dep_param(&self, inst_w: f64, inst_l: f64, nf: f64, temp: f64) -> Bsim4SizeParam {
        let l_new = inst_l + self.xl;
        let w_new = inst_w / nf + self.xw;

        // Effective dimensions with binning corrections
        let t0 = l_new.powf(self.lln);
        let t1 = w_new.powf(self.lwn);
        let dl = self.lint + self.ll / t0 + self.lw / t1 + self.lwl / (t0 * t1);

        let t2 = l_new.powf(self.wln);
        let t3 = w_new.powf(self.wwn);
        let dw = self.wint + self.wl / t2 + self.ww / t3 + self.wwl / (t2 * t3);

        let leff = (l_new - 2.0 * dl).max(1.0e-9);
        let weff = (w_new - 2.0 * dw).max(1.0e-9);

        let dlc = if self.dlc_given { self.dlc } else { dl };
        let dwc = if self.dwc_given { self.dwc } else { dw };

        let leff_cv = (l_new - 2.0 * dlc).max(1.0e-9);
        let weff_cv = (w_new - 2.0 * dwc).max(1.0e-9);

        let dwj = if self.dwc_given { self.dwj } else { dwc };
        let weff_cj = (w_new - 2.0 * dwj).max(1.0e-9);

        // Binning inverse dimensions
        let (inv_l, inv_w, inv_lw) = if self.bin_unit == 1 {
            (1.0e-6 / leff, 1.0e-6 / weff, 1.0e-12 / (leff * weff))
        } else {
            (1.0 / leff, 1.0 / weff, 1.0 / (leff * weff))
        };

        // Bin all parameters
        macro_rules! bin {
            ($base:expr, $name:expr) => {
                self.bin_param($base, $name, inv_l, inv_w, inv_lw)
            };
        }

        let vth0 = bin!(self.vth0, "vth0");
        let k1 = bin!(self.k1, "k1");
        let k2 = bin!(self.k2, "k2");
        let k3 = bin!(self.k3, "k3");
        let k3b = bin!(self.k3b, "k3b");
        let w0 = bin!(self.w0, "w0");
        let lpe0 = bin!(self.lpe0, "lpe0");
        let lpeb = bin!(self.lpeb, "lpeb");
        let dvtp0 = bin!(self.dvtp0, "dvtp0");
        let dvtp1 = bin!(self.dvtp1, "dvtp1");
        let dvt0 = bin!(self.dvt0, "dvt0");
        let dvt1 = bin!(self.dvt1, "dvt1");
        let dvt2 = bin!(self.dvt2, "dvt2");
        let dvt0w = bin!(self.dvt0w, "dvt0w");
        let dvt1w = bin!(self.dvt1w, "dvt1w");
        let dvt2w = bin!(self.dvt2w, "dvt2w");
        let dsub = bin!(self.dsub, "dsub");
        let minv = bin!(self.minv, "minv");
        let eta0 = bin!(self.eta0, "eta0");
        let etab = bin!(self.etab, "etab");
        let voff = bin!(self.voff, "voff");
        let nfactor = bin!(self.nfactor, "nfactor");
        let cdsc = bin!(self.cdsc, "cdsc");
        let cdscb = bin!(self.cdscb, "cdscb");
        let cdscd = bin!(self.cdscd, "cdscd");
        let cit = bin!(self.cit, "cit");
        let phin = bin!(self.phin, "phin");
        let u0 = bin!(self.u0, "u0");
        let ua = bin!(self.ua, "ua");
        let ub = bin!(self.ub, "ub");
        let uc = bin!(self.uc, "uc");
        let ud = bin!(self.ud, "ud");
        let vsat = bin!(self.vsat, "vsat");
        let a0 = bin!(self.a0, "a0");
        let ags = bin!(self.ags, "ags");
        let a1 = bin!(self.a1, "a1");
        let a2 = bin!(self.a2, "a2");
        let keta = bin!(self.keta, "keta");
        let pclm = bin!(self.pclm, "pclm");
        let pdibl1 = bin!(self.pdibl1, "pdibl1");
        let pdibl2 = bin!(self.pdibl2, "pdibl2");
        let pdiblb = bin!(self.pdiblb, "pdiblb");
        let pscbe1 = bin!(self.pscbe1, "pscbe1");
        let pscbe2 = bin!(self.pscbe2, "pscbe2");
        let pvag = bin!(self.pvag, "pvag");
        let delta = bin!(self.delta, "delta");
        let fprout = bin!(self.fprout, "fprout");
        let pdits = bin!(self.pdits, "pdits");
        let pditsd = bin!(self.pditsd, "pditsd");
        let drout = bin!(self.drout, "drout");
        let rdsw = bin!(self.rdsw, "rdsw");
        let prwg = bin!(self.prwg, "prwg");
        let prwb = bin!(self.prwb, "prwb");
        let dwg = bin!(self.dwg, "dwg");
        let dwb = bin!(self.dwb, "dwb");
        let b0 = bin!(self.b0, "b0");
        let b1 = bin!(self.b1, "b1");
        let alpha0 = bin!(self.alpha0, "alpha0");
        let alpha1 = bin!(self.alpha1, "alpha1");
        let beta0 = bin!(self.beta0, "beta0");
        let agidl = bin!(self.agidl, "agidl");
        let bgidl = bin!(self.bgidl, "bgidl");
        let cgidl = bin!(self.cgidl, "cgidl");
        let egidl = bin!(self.egidl, "egidl");
        let agisl = bin!(self.agisl, "agisl");
        let bgisl = bin!(self.bgisl, "bgisl");
        let cgisl = bin!(self.cgisl, "cgisl");
        let egisl = bin!(self.egisl, "egisl");
        let aigc = bin!(self.aigc, "aigc");
        let bigc = bin!(self.bigc, "bigc");
        let cigc = bin!(self.cigc, "cigc");
        let aigsd = bin!(self.aigsd, "aigsd");
        let bigsd = bin!(self.bigsd, "bigsd");
        let cigsd = bin!(self.cigsd, "cigsd");
        let nigc = bin!(self.nigc, "nigc");
        let pigcd = bin!(self.pigcd, "pigcd");
        let poxedge = bin!(self.poxedge, "poxedge");
        let ntox = bin!(self.ntox, "ntox");
        let aigbacc = bin!(self.aigbacc, "aigbacc");
        let bigbacc = bin!(self.bigbacc, "bigbacc");
        let cigbacc = bin!(self.cigbacc, "cigbacc");
        let nigbacc = bin!(self.nigbacc, "nigbacc");
        let aigbinv = bin!(self.aigbinv, "aigbinv");
        let bigbinv = bin!(self.bigbinv, "bigbinv");
        let cigbinv = bin!(self.cigbinv, "cigbinv");
        let eigbinv = bin!(self.eigbinv, "eigbinv");
        let nigbinv = bin!(self.nigbinv, "nigbinv");
        let xrcrg1 = bin!(self.xrcrg1, "xrcrg1");
        let xrcrg2 = bin!(self.xrcrg2, "xrcrg2");
        let cgsl = bin!(self.cgsl, "cgsl");
        let cgdl = bin!(self.cgdl, "cgdl");
        let ckappas = bin!(self.ckappas, "ckappas");
        let ckappad = bin!(self.ckappad, "ckappad");
        let cf = bin!(self.cf, "cf");
        let clc = bin!(self.clc, "clc");
        let cle = bin!(self.cle, "cle");
        let acde = bin!(self.acde, "acde");
        let moin = bin!(self.moin, "moin");
        let noff = bin!(self.noff, "noff");
        let voffcv = bin!(self.voffcv, "voffcv");
        let vfbcv = bin!(self.vfbcv, "vfbcv");
        let kt1 = bin!(self.kt1, "kt1");
        let kt1l = bin!(self.kt1l, "kt1l");
        let kt2 = bin!(self.kt2, "kt2");
        let ndep = bin!(self.ndep, "ndep");
        let ngate = bin!(self.ngate, "ngate");
        let nsub = bin!(self.nsub, "nsub");
        let nsd = bin!(self.nsd, "nsd");
        let xj = bin!(self.xj, "xj");
        let vfb = bin!(self.vfb, "vfb");

        let abulk_cv_factor = 1.0 + (clc / leff_cv).powf(cle);

        // Temperature computations
        let tnom = self.tnom;
        let vtm0 = KBOQ * tnom;
        let _vtm = KBOQ * temp;
        let t_ratio = temp / tnom;
        let del_temp = temp - tnom;

        // Oxide capacitances
        let coxe = self.epsrox * EPS0 / self.toxe;
        let coxp = self.epsrox * EPS0 / self.toxp;
        let factor1 = (EPSSI / (self.epsrox * EPS0) * self.toxe).sqrt();

        // Bandgap
        let eg0 = 1.16 - 7.02e-4 * tnom * tnom / (tnom + 1108.0);
        let ni = 1.45e10
            * (tnom / 300.15)
            * (tnom / 300.15).sqrt()
            * (21.5565981 - eg0 / (2.0 * vtm0)).exp();
        let _eg = 1.16 - 7.02e-4 * temp * temp / (temp + 1108.0);

        // Surface potential
        let ndep_cm3 = ndep.max(1.0e15);
        let phi = vtm0 * (ndep_cm3 / ni).ln() + phin + 0.4;
        let phi = phi.max(0.1);
        let sqrt_phi = phi.sqrt();
        let phis3 = sqrt_phi * phi;

        // Depletion width
        let xdep0 = (2.0 * EPSSI / (CHARGE_Q * ndep_cm3 * 1.0e6)).sqrt() * sqrt_phi;
        let sqrt_xdep0 = xdep0.sqrt();

        // litl (lateral diffusion)
        let litl = (3.9 / self.epsrox * xj * self.toxe).sqrt();

        // Built-in voltage
        let vbi = vtm0 * (nsd * ndep_cm3 / (ni * ni)).ln();

        // cdep0
        let cdep0 = (CHARGE_Q * EPSSI * ndep_cm3 * 1.0e6 / 2.0 / phi).sqrt();

        // k1ox
        let k1ox = k1 * self.toxe / self.toxm;
        let k2ox = k2 * self.toxe / self.toxm;

        // Temperature-dependent mobility
        let mut u0_temp = u0;
        if u0_temp > 1.0 {
            u0_temp /= 1.0e4; // convert from cm²/Vs
        }
        let t5 = 1.0 - self.up * (-leff / self.lp).exp();
        u0_temp = u0_temp * t5 * t_ratio.powf(self.ute);

        // Temperature-dependent ua, ub, uc, ud (tempMod=0)
        let t0 = t_ratio - 1.0;
        let ua = ua + self.ua1 * t0;
        let ub = ub + self.ub1 * t0;
        let uc = uc + self.uc1 * t0;
        let ud = ud + self.ud1 * t0;

        // Temperature-dependent vsat
        let vsattemp = (vsat - self.at * t0).max(1.0e3);

        // Series resistance
        let pow_weff_wr = (weff_cj * 1.0e6).powf(self.wr) * nf;
        let t10 = self.prt * t0;
        let rds0 = ((rdsw + t10) * nf / pow_weff_wr).max(0.0);
        let rdswmin_val = ((self.rdswmin + t10) * nf / pow_weff_wr).max(0.0);

        // Temperature-dependent voff
        let voff = voff * (1.0 + self.get_l("tvoff") * del_temp);
        // theta0vb0 (DIBL)
        let tmp = (EPSSI / (self.epsrox * EPS0) * self.toxe * xdep0).sqrt();
        let t0_dibl = dsub * leff / tmp;
        let theta0vb0 = if t0_dibl < EXP_THRESHOLD {
            let t1 = t0_dibl.exp();
            let t2 = t1 - 1.0;
            let t3 = t2 * t2;
            let t4 = t3 + 2.0 * t1 * MIN_EXP;
            t1 / t4
        } else {
            1.0 / (MAX_EXP - 2.0)
        };

        // thetaRout
        let t0_rout = drout * leff / tmp;
        let theta_rout_raw = if t0_rout < EXP_THRESHOLD {
            let t1 = t0_rout.exp();
            let t2 = t1 - 1.0;
            let t3 = t2 * t2;
            let t4 = t3 + 2.0 * t1 * MIN_EXP;
            t1 / t4
        } else {
            1.0 / (MAX_EXP - 2.0)
        };
        let theta_rout = pdibl1 * theta_rout_raw + pdibl2;

        // vfbzb_factor - flat-band + correction
        let type_sign = self.mos_type.sign();

        // dvt0/dvt1 short-channel correction
        let t0_sc = dvt1 * leff / tmp;
        let t9 = if t0_sc < EXP_THRESHOLD {
            let t1 = t0_sc.exp();
            let t2 = t1 - 1.0;
            let t3 = t2 * t2 + 2.0 * t1 * MIN_EXP;
            dvt0 * t1 / t3
        } else {
            dvt0 / (MAX_EXP - 2.0)
        };

        // dvt0w/dvt1w narrow-width correction
        let t1_nw = dvt1w * weff * leff / tmp;
        let t8 = if t1_nw < EXP_THRESHOLD {
            let t1 = t1_nw.exp();
            let t2 = t1 - 1.0;
            let t3 = t2 * t2 + 2.0 * t1 * MIN_EXP;
            dvt0w * t1 / t3
        } else {
            dvt0w / (MAX_EXP - 2.0)
        };

        let v0 = vbi - phi;
        let lpe_vb = (1.0 + lpe0 / leff).sqrt();
        let t5_vfb = k1ox * (lpe_vb - 1.0) * sqrt_phi
            + (kt1 + kt1l / leff + kt2 * (type_sign * vbi - phi)) * (t_ratio - 1.0);

        let vth_narrow_w = self.toxe * phi / (weff + w0);
        let vfbzb_factor = -t9 * v0 - t8 * v0 + k3 * vth_narrow_w + t5_vfb - phi - k1 * sqrt_phi;

        // Overlap caps
        let cgso_eff = (self.cgso + cf) * weff_cv;
        let cgdo_eff = (self.cgdo + cf) * weff_cv;
        let cgbo_eff = self.cgbo * leff_cv * nf;

        // Subthreshold
        let mstar = 0.5 + (minv).atan() / std::f64::consts::PI;
        let voffcbn = voff + self.voffl / leff;
        let ldeb = (EPSSI * vtm0 / (CHARGE_Q * ndep_cm3 * 1.0e6)).sqrt() / 3.0;
        let acde_scaled = acde * (ndep_cm3 / 2.0e16).powf(-0.25);

        // Gate tunneling
        let tox_ratio = (ntox * (self.toxref / self.toxe).ln()).exp() / (self.toxe * self.toxe);
        let tox_ratio_edge = (ntox * (self.toxref / (self.toxe * poxedge)).ln()).exp()
            / (self.toxe * poxedge * self.toxe * poxedge);

        let (aechvb_base, bechvb_base) = if self.mos_type == MosfetType::Nmos {
            (4.97232e-7, 7.45669e11)
        } else {
            (3.42537e-7, 1.16645e12)
        };

        let dlcig = if self.dlc_given { self.dlc } else { self.lint };

        let aechvb = aechvb_base * weff * leff * tox_ratio;
        let bechvb = -bechvb_base * self.toxe;
        let aechvb_edge_s = aechvb_base * weff * dlcig * tox_ratio_edge;
        let aechvb_edge_d = aechvb_base * weff * dlcig * tox_ratio_edge;
        let bechvb_edge = -bechvb_base * self.toxe * poxedge;

        // vtfbphi1 for mobMod 2/4/5/6 (matching b4temp.c)
        let t3_vtfb = type_sign * vth0 - vfb - phi;
        let vtfbphi1 = if self.mos_type == MosfetType::Nmos {
            (2.0 * t3_vtfb).max(0.0)
        } else {
            (2.5 * t3_vtfb).max(0.0)
        };

        Bsim4SizeParam {
            leff,
            weff,
            leff_cv,
            weff_cv,
            weff_cj,
            vth0,
            k1,
            k2,
            k3,
            k3b,
            w0,
            lpe0,
            lpeb,
            dvtp0,
            dvtp1,
            dvt0,
            dvt1,
            dvt2,
            dvt0w,
            dvt1w,
            dvt2w,
            dsub,
            minv,
            eta0,
            etab,
            voff,
            nfactor,
            cdsc,
            cdscb,
            cdscd,
            cit,
            phin,
            u0: u0_temp,
            ua,
            ub,
            uc,
            ud,
            vsat: vsattemp,
            a0,
            ags,
            a1,
            a2,
            keta,
            pclm,
            pdibl1,
            pdibl2,
            pdiblb,
            pscbe1,
            pscbe2,
            pvag,
            delta,
            fprout,
            pdits,
            pditsd,
            drout,
            rdsw,
            prwg,
            prwb,
            dwg,
            dwb,
            b0,
            b1,
            alpha0,
            alpha1,
            beta0,
            agidl,
            bgidl,
            cgidl,
            egidl,
            agisl,
            bgisl,
            cgisl,
            egisl,
            aigc,
            bigc,
            cigc,
            aigsd,
            bigsd,
            cigsd,
            nigc,
            pigcd,
            poxedge,
            ntox,
            aigbacc,
            bigbacc,
            cigbacc,
            nigbacc,
            aigbinv,
            bigbinv,
            cigbinv,
            eigbinv,
            nigbinv,
            xrcrg1,
            xrcrg2,
            cgsl,
            cgdl,
            ckappas,
            ckappad,
            cf,
            clc,
            cle,
            acde: acde_scaled,
            moin,
            noff,
            voffcv,
            vfbcv,
            kt1,
            kt1l,
            kt2,
            ndep: ndep_cm3,
            ngate,
            nsub,
            nsd,
            xj,
            vfb,
            u0temp: u0_temp,
            vsattemp,
            rds0,
            rdswmin: rdswmin_val,
            phi,
            sqrt_phi,
            phis3,
            xdep0,
            sqrt_xdep0,
            litl,
            vbi,
            cdep0,
            k1ox,
            k2ox,
            theta0vb0,
            theta_rout,
            vfbzb_factor,
            coxe,
            coxp,
            factor1,
            cgso_eff,
            cgdo_eff,
            cgbo_eff,
            mstar,
            voffcbn,
            ldeb,
            aechvb,
            bechvb,
            aechvb_edge_s,
            aechvb_edge_d,
            bechvb_edge,
            tox_ratio,
            tox_ratio_edge,
            abulk_cv_factor,
            vtfbphi1,
        }
    }
}

/// Check if a parameter name (after stripping l/w/p prefix) is a binnable BSIM4 param.
fn is_binnable_param(name: &str) -> bool {
    matches!(
        name,
        "vth0"
            | "k1"
            | "k2"
            | "k3"
            | "k3b"
            | "w0"
            | "lpe0"
            | "lpeb"
            | "dvtp0"
            | "dvtp1"
            | "dvt0"
            | "dvt1"
            | "dvt2"
            | "dvt0w"
            | "dvt1w"
            | "dvt2w"
            | "dsub"
            | "minv"
            | "eta0"
            | "etab"
            | "voff"
            | "nfactor"
            | "cdsc"
            | "cdscb"
            | "cdscd"
            | "cit"
            | "phin"
            | "u0"
            | "ua"
            | "ub"
            | "uc"
            | "ud"
            | "ute"
            | "ua1"
            | "ub1"
            | "uc1"
            | "ud1"
            | "vsat"
            | "a0"
            | "ags"
            | "a1"
            | "a2"
            | "at"
            | "keta"
            | "pclm"
            | "pdibl1"
            | "pdibl2"
            | "pdiblb"
            | "pscbe1"
            | "pscbe2"
            | "pvag"
            | "delta"
            | "fprout"
            | "pdits"
            | "pditsd"
            | "drout"
            | "rdsw"
            | "prwg"
            | "prwb"
            | "dwg"
            | "dwb"
            | "b0"
            | "b1"
            | "alpha0"
            | "alpha1"
            | "beta0"
            | "agidl"
            | "bgidl"
            | "cgidl"
            | "egidl"
            | "agisl"
            | "bgisl"
            | "cgisl"
            | "egisl"
            | "aigc"
            | "bigc"
            | "cigc"
            | "aigsd"
            | "bigsd"
            | "cigsd"
            | "nigc"
            | "pigcd"
            | "poxedge"
            | "ntox"
            | "aigbacc"
            | "bigbacc"
            | "cigbacc"
            | "nigbacc"
            | "aigbinv"
            | "bigbinv"
            | "cigbinv"
            | "eigbinv"
            | "nigbinv"
            | "xrcrg1"
            | "xrcrg2"
            | "cgsl"
            | "cgdl"
            | "ckappas"
            | "ckappad"
            | "cf"
            | "clc"
            | "cle"
            | "acde"
            | "moin"
            | "noff"
            | "voffcv"
            | "vfbcv"
            | "kt1"
            | "kt1l"
            | "kt2"
            | "ndep"
            | "ngate"
            | "nsub"
            | "nsd"
            | "xj"
            | "vfb"
            | "kvth0we"
            | "k2we"
            | "ku0we"
            | "tvoff"
    )
}

/// Compute BSIM4 companion model at given operating point.
#[expect(clippy::too_many_lines)]
pub fn bsim4_companion(
    vgs: f64,
    vds: f64,
    vbs: f64,
    sp: &Bsim4SizeParam,
    model: &Bsim4Model,
) -> Bsim4Companion {
    let type_sign = model.mos_type.sign();
    let gmin = 1.0e-12;

    // Mode detection (forward/reverse)
    let (mode, vgs_eff, vds_eff_raw, vbs_eff_raw) = if vds >= 0.0 {
        (1, vgs, vds, vbs)
    } else {
        (-1, vgs - vds, -vds, vbs - vds)
    };

    let vds_mode = vds_eff_raw;
    let vgs_mode = vgs_eff;
    let vbs_mode = vbs_eff_raw;

    // Effective Vbs (smooth clamp to phi)
    let t0 = sp.phi - 0.01;
    let vbs_eff = if vbs_mode < t0 {
        vbs_mode
    } else {
        let t1 = (t0 - vbs_mode - DELTA_1).sqrt();
        let t2 = (t1 * t1 + 4.0 * DELTA_1 * t0).sqrt();
        t0 - 0.5 * (t1 + t2)
    };

    // Surface potential
    let sqrt_phis = (sp.phi - vbs_eff).max(0.01).sqrt();

    // Vtm
    let vtm = KBOQ * TEMP_DEFAULT;

    // Lpe_Vb
    let lpe_vb = (1.0 + sp.lpeb / sp.leff).sqrt();

    // Narrow-width effect
    let vth_narrow_w = model.toxe * sp.phi / (sp.weff + sp.w0);

    // DIBL shift
    let dibl_sft = (sp.eta0 + sp.etab * vbs_eff) * sp.theta0vb0 * vds_mode;

    // Short-channel effect (Delt_vth)
    let v0 = sp.vbi - sp.phi;
    let tmp_sc = (EPSSI / (model.epsrox * EPS0) * model.toxe * sp.xdep0).sqrt();

    let t0_sc = sp.dvt1 * sp.leff / tmp_sc;
    let theta0 = if t0_sc < EXP_THRESHOLD {
        let t1 = t0_sc.exp();
        let t2 = t1 - 1.0;
        let t3 = t2 * t2 + 2.0 * t1 * MIN_EXP;
        sp.dvt0 * t1 / t3
    } else {
        sp.dvt0 / (MAX_EXP - 2.0)
    };
    let delt_vth = theta0 * (v0 + vbs_eff * sp.dvt2);

    // Narrow-width SCE
    let t1_nw = sp.dvt1w * sp.weff * sp.leff / tmp_sc;
    let t2_nw = if t1_nw < EXP_THRESHOLD {
        let t1 = t1_nw.exp();
        let t2 = t1 - 1.0;
        let t3 = t2 * t2 + 2.0 * t1 * MIN_EXP;
        sp.dvt0w * t1 / t3
    } else {
        sp.dvt0w / (MAX_EXP - 2.0)
    };
    let delt_vth_nw = t2_nw * (v0 + vbs_eff * sp.dvt2w);

    // Temperature correction
    let kt_term =
        (sp.kt1 + sp.kt1l / sp.leff + sp.kt2 * vbs_eff) * (TEMP_DEFAULT / model.tnom - 1.0);

    // Threshold voltage
    let mut vth = type_sign * sp.vth0 + (sp.k1ox * sqrt_phis - sp.k1 * sp.sqrt_phi) * lpe_vb
        - sp.k2ox * vbs_eff
        - delt_vth
        - delt_vth_nw
        + (sp.k3 + sp.k3b * vbs_eff) * vth_narrow_w
        + kt_term
        - dibl_sft;

    // Pocket implant (dvtp0/dvtp1)
    if sp.dvtp0 > 0.0 {
        let t0 = -sp.dvtp1 * vds_mode;
        let t0 = if t0 < -EXP_THRESHOLD {
            MIN_EXP
        } else {
            t0.exp()
        };
        vth -= vtm * (sp.leff / (sp.leff + sp.dvtp0 * (1.0 + t0))).ln();
    }

    // Subthreshold slope factor n
    let tmp_n = sp.cdsc + sp.cdscb * vbs_eff + sp.cdscd * vds_mode + sp.cit;
    let tmp_n = tmp_n / sp.coxe;

    let mut n = 1.0 + sp.nfactor * EPSSI / sp.xdep0 / sp.coxe;
    n += tmp_n;
    n = n.max(0.1);

    // Poly depletion effect
    let vgs_eff_poly = if sp.ngate > 1.0e18 && sp.ngate < 1.0e25 {
        poly_depletion(vgs_mode, type_sign, sp)
    } else {
        vgs_mode
    };

    // Effective Vgst
    let vgst = vgs_eff_poly - vth;

    // Smooth weak-to-strong inversion transition
    let t0_vgst = sp.mstar * vgst;
    let t1_vgst = n * vtm;
    let t2_vgst = t0_vgst / t1_vgst;

    let vgsteff = if t2_vgst > EXP_THRESHOLD {
        vgst
    } else if t2_vgst < -EXP_THRESHOLD {
        t1_vgst * MIN_EXP.ln()
    } else {
        let t3 = t2_vgst.exp();
        t1_vgst * (1.0 + t3).ln()
    };
    let vgsteff = vgsteff.max(1.0e-50);

    let vgst2vtm = vgsteff + 2.0 * vtm;

    // Effective width
    let weff_local = sp.weff - 2.0 * (sp.dwg * vgsteff + sp.dwb * (sqrt_phis - sp.sqrt_phi));

    // Abulk (body effect coefficient)
    let t0_abulk = sp.ags * sp.a0;
    let tmp_abulk = 1.0 / (1.0 + t0_abulk * vgsteff);
    let abulk0 = 1.0 + sp.a1 * (1.0 + sp.a2 * sp.phi) + sp.k1ox * (lpe_vb - 1.0) * 0.5 / sqrt_phis;
    let abulk = abulk0 * tmp_abulk;
    let abulk = abulk.max(0.1);

    // Mobility degradation (mobMod dispatch)
    let t14 = 0.0; // mtrlMod=0
    let toxe = model.toxe;

    let t5_mob = match model.mob_mod {
        1 => {
            // mobMod=1: multiplicative form
            let t0 = vgsteff + vth + vth - t14;
            let t2 = 1.0 + sp.uc * vbs_eff;
            let t3 = t0 / toxe;
            let t4 = t3 * (sp.ua + sp.ub * t3);
            let t12 = (vth * vth + 0.0001).sqrt();
            let t9 = 1.0 / (vgsteff + 2.0 * t12);
            let t10 = t9 * toxe;
            let t8 = sp.ud * t10 * t10 * vth;
            let t6 = t8 * vth;
            t4 * t2 + t6
        }
        2 => {
            // mobMod=2: power-law (universal mobility)
            let t0 = (vgsteff + sp.vtfbphi1) / toxe;
            let t0 = t0.max(1.0e-20);
            let t1 = t0.powf(model.eu);
            let t2 = sp.ua + sp.uc * vbs_eff;
            let t12 = (vth * vth + 0.0001).sqrt();
            let t9 = 1.0 / (vgsteff + 2.0 * t12);
            let t10 = t9 * toxe;
            let t8 = sp.ud * t10 * t10 * vth;
            let t6 = t8 * vth;
            t1 * t2 + t6
        }
        _ => {
            // mobMod=0 (default): linear degradation
            let t0 = vgsteff + vth + vth - t14;
            let t3 = t0 / toxe;
            let t12 = (vth * vth + 0.0001).sqrt();
            let t9 = 1.0 / (vgsteff + 2.0 * t12);
            let t10 = t9 * toxe;
            let t8 = sp.ud * t10 * t10 * vth;
            let t6 = t8 * vth;
            let t2 = sp.ua + sp.uc * vbs_eff;
            t3 * (t2 + sp.ub * t3) + t6
        }
    };

    let denomi = if t5_mob >= -0.8 {
        1.0 + t5_mob
    } else {
        let t9 = 1.0 / (7.0 + 10.0 * t5_mob);
        (0.6 + t5_mob) * t9
    };

    let ueff = sp.u0 / denomi;
    let ueff = ueff.max(1.0e-15);

    // EsatL and Vdsat
    let esat_l = 2.0 * sp.vsat / ueff * sp.leff;
    let esat = 2.0 * sp.vsat / ueff;

    // Series resistance Rds
    let rds = sp.rds0 / (1.0 + sp.prwg * vgsteff + sp.prwb * (sqrt_phis - sp.sqrt_phi));
    let rds = rds.max(0.0);

    // Vdsat
    // Vdsat: simplified — rds effect on Vdsat is small for rdsmod=0
    let vdsat = esat_l * vgst2vtm / (abulk * esat_l + vgst2vtm);
    let vdsat = vdsat.max(1.0e-18);

    // Vdseff (smooth transition at Vdsat)
    let t1_vdseff = vdsat - vds_mode - sp.delta;
    let t2_vdseff = (t1_vdseff * t1_vdseff + 4.0 * sp.delta * vdsat).sqrt();
    let vdseff = vdsat - 0.5 * (t1_vdseff + t2_vdseff);
    let vdseff = vdseff.max(1.0e-18);

    let diff_vds = vds_mode - vdseff;

    // Channel current computation
    let cox_eff_wov_l = sp.coxe * weff_local / sp.leff;
    let beta = ueff * cox_eff_wov_l;
    let fgche1 = vgsteff * (1.0 - 0.5 * abulk * vdseff / vgst2vtm);
    let fgche2 = 1.0 + vdseff / esat_l;
    let gche = beta * fgche1 / fgche2;
    let idl = gche * vdseff / (1.0 + gche * rds);
    let idl = idl.max(0.0);

    // Output resistance effects
    // Vasat
    let vasat = esat_l - vdsat + vdsat * (1.0 - 0.5 * abulk * vdsat / vgst2vtm) + esat_l;
    let vasat = vasat.max(1.0);

    // VACLM (Channel Length Modulation)
    let pvag_term = 1.0 + sp.pvag * vgsteff / esat_l;
    let fp = if sp.fprout > 0.0 {
        1.0 / (1.0 + sp.fprout * (sp.phi - vbs_eff).max(0.0).sqrt())
    } else {
        1.0
    };

    let cclm = fp * pvag_term * (1.0 + rds * idl) * (sp.leff + vdsat / esat) / (sp.pclm * sp.litl);
    let cclm = cclm.max(1.0e-18);
    let vaclm = cclm * diff_vds;
    let vaclm = vaclm.max(1.0e-18);

    // VADIBL
    let vadibl = if sp.theta_rout > 0.0 {
        let _t0 = vgst2vtm - abulk * vdsat;
        let t1 = vgst2vtm + abulk * vdsat;
        let vadibl_raw = (vgst2vtm - vgst2vtm * abulk * vdsat / t1) / sp.theta_rout;
        let vadibl_raw = vadibl_raw / (1.0 + sp.pdiblb * vbs_eff);
        (vadibl_raw * pvag_term).max(1.0e-18)
    } else {
        1.0e30
    };

    // VADITS (matching b4ld.c lines 2003-2024)
    let vadits = if sp.pdits > MIN_EXP {
        let t0 = sp.pditsd * vds_mode;
        let t1 = if t0 > EXP_THRESHOLD {
            MAX_EXP
        } else {
            t0.exp()
        };
        let t2 = 1.0 + model.pditsl * sp.leff;
        let vadits_raw = (1.0 + t2 * t1) / sp.pdits;
        (vadits_raw * fp).max(1.0e-18)
    } else {
        MAX_EXP
    };

    // VASCBE
    let vascbe = if diff_vds > 1.0e-9 && sp.pscbe1 > 0.0 {
        let t0 = sp.pscbe1 * sp.litl / diff_vds;
        sp.leff * safe_exp(t0) / sp.pscbe2
    } else {
        1.0e30
    };
    let vascbe = vascbe.max(1.0e-18);

    // Combined Early voltage
    let va = vasat + vaclm;

    // Final drain current with all corrections
    let mut idsa = idl * (1.0 + diff_vds / vadibl);
    idsa *= 1.0 + diff_vds / vadits;
    if va > vasat {
        idsa *= 1.0 + (va / vasat).ln() / cclm;
    }
    let ids = idsa * (1.0 + diff_vds / vascbe);
    let ids = ids * vdseff;

    // Transconductances (simplified — using numerical-like approach as in BSIM3)
    let dvdseff_dvg = {
        let t1 = vdsat - vds_mode - sp.delta;
        let t2 = (t1 * t1 + 4.0 * sp.delta * vdsat).sqrt();
        if t2 > 1.0e-20 {
            0.5 * (1.0 - t1 / t2) // approximate dVdseff/dVdsat * dVdsat/dVgs
        } else {
            0.0
        }
    };

    // Approximate gm using chain rule
    let gm_raw =
        beta * vdseff / fgche2 * (1.0 - 0.5 * abulk * vdseff / vgst2vtm) / (1.0 + gche * rds);
    let gds_raw = beta * vgsteff / fgche2 * (1.0 - vdseff / esat_l) / (1.0 + gche * rds)
        * (1.0 - dvdseff_dvg);

    // Apply output resistance corrections to gds
    let gm = gm_raw * vdseff;
    let mut gds = gds_raw.abs() * vdseff + ids / vascbe;
    if ids > 0.0 && va > vasat {
        gds += ids / vaclm;
    }
    gds += ids / vadibl;

    // gmbs (body effect transconductance)
    let gmbs_factor = if (sp.phi - vbs_eff).abs() > 0.01 {
        0.5 * sp.k1ox * lpe_vb / sqrt_phis
    } else {
        0.0
    };
    let gmbs = gm * gmbs_factor;

    // Junction currents
    let (gbd, cbd_current) =
        junction_current(vbs - vds, 0.0, vtm * sp.nsd.max(1.0), gmin, model.ijthdfwd);
    let (gbs, cbs_current) =
        junction_current(vbs, 0.0, vtm * sp.nsd.max(1.0), gmin, model.ijthsfwd);

    // Substrate current (impact ionization)
    let isub = if diff_vds > 1.0e-9 && (sp.alpha0 + sp.alpha1 * sp.leff) > 0.0 {
        let t2 = (sp.alpha0 + sp.alpha1 * sp.leff) / sp.leff;
        let t0 = -sp.beta0 / diff_vds;
        t2 * diff_vds * safe_exp(t0) * idsa
    } else {
        0.0
    };

    // GIDL (simplified)
    let igidl = if sp.agidl > 0.0 {
        let vgs_eff_gidl = type_sign * vgs_mode;
        let t0 = 3.0 * model.toxe;
        let t1 = (vds_mode - vgs_eff_gidl - sp.egidl) / t0;
        if t1 > 0.27 {
            let t2 = sp.bgidl / t1;
            if t2 < 100.0 {
                sp.agidl * sp.weff_cj * t1 * (-t2).exp()
            } else {
                0.0
            }
        } else {
            0.0
        }
    } else {
        0.0
    };

    // GISL (mirror of GIDL using source side)
    let igisl = if sp.agisl > 0.0 {
        let vgd_eff_gisl = type_sign * (vgs_mode - vds_mode);
        let t0 = 3.0 * model.toxe;
        let t1 = (-vds_mode - vgd_eff_gisl - sp.egisl) / t0;
        if t1 > 0.27 {
            let t2 = sp.bgisl / t1;
            if t2 < 100.0 {
                sp.agisl * sp.weff_cj * t1 * (-t2).exp()
            } else {
                0.0
            }
        } else {
            0.0
        }
    } else {
        0.0
    };

    // NR equivalent currents (no type_sign here — stamp_bsim4 applies it)
    let ceq_d = ids - gds * vds_mode - gm * vgs_mode - gmbs * vbs_mode;
    let ceq_bs = cbs_current - gbs * vbs;
    let ceq_bd = cbd_current - gbd * (vbs - vds);

    let ceq_sub = isub;

    // Capacitances — capMod=2 CTM charge model with proper derivatives (port of b4ld.c)
    let vfbzb = sp.vfbzb_factor + type_sign * sp.vth0;
    let cox_wl = sp.coxe * weff_local * sp.leff;

    // Gate charge model
    let vfbeff = {
        let t0 = vgs_mode - vfbzb - 0.02;
        let t1 = (t0 * t0 + 0.08 * vfbzb.abs()).sqrt();
        vfbzb + 0.5 * (t0 + t1)
    };

    let qac0 = cox_wl * (vfbeff - vfbzb);

    // Inversion charge
    let abulk_cv = (1.0 + sp.abulk_cv_factor * sp.k1ox / (2.0 * sqrt_phis)).max(0.1);
    let vdsat_cv = vgsteff / abulk_cv;

    let t1_cv = vdsat_cv - vds_mode - sp.delta;
    let t2_cv = (t1_cv * t1_cv + 4.0 * sp.delta * vdsat_cv).sqrt();
    let vdseff_cv = (vdsat_cv - 0.5 * (t1_cv + t2_cv)).max(1.0e-18);

    // dVdseffCV/dVds — derivative of VdseffCV smoothing w.r.t. Vds
    let dvdseff_dvd = if t2_cv > 1.0e-30 {
        0.5 * (1.0 + t1_cv / t2_cv)
    } else {
        0.5
    };

    let t0_q = abulk_cv * vdseff_cv;
    let t2_q = 12.0 * (vgsteff - 0.5 * t0_q + 1.0e-20);
    let t3_q = t0_q / t2_q;

    // Gate inversion charge
    let qgate_inv = cox_wl * (vgsteff - t0_q * (0.5 - t3_q));

    // Bulk charge
    let qbulk_val = cox_wl * (1.0 - abulk_cv) * (0.5 * vdseff_cv - t0_q * vdseff_cv / t2_q);

    // Gate charge derivatives (matching b4ld.c T4/T5 terms)
    let t4 = 1.0 - 12.0 * t3_q * t3_q; // dQg/dVg (normalized by CoxWL)
    let t5 = abulk_cv * (6.0 * t0_q * (4.0 * vgsteff - t0_q) / (t2_q * t2_q) - 0.5);

    let cgg1 = cox_wl * t4;
    let cgd1 = cox_wl * t5 * dvdseff_dvd;

    // Bulk charge derivatives
    let t7 = 1.0 - abulk_cv;
    let t8 = t2_q * t2_q;
    let t9 = if t8.abs() > 1.0e-40 && abulk_cv.abs() > 1.0e-20 {
        12.0 * t7 * t0_q * t0_q / (t8 * abulk_cv)
    } else {
        0.0
    };
    let t11 = if abulk_cv.abs() > 1.0e-20 {
        -t7 * t5 / abulk_cv
    } else {
        0.0
    };

    let cbg1 = cox_wl * t9;
    let cbd1 = cox_wl * t11 * dvdseff_dvd;

    // Source charge and derivatives — xpart-dependent partitioning
    let xpart = model.xpart;
    let (csg, csd) = if xpart > 0.5 {
        // 0/100 partition: all inversion charge to source
        let t2x = 2.0 * t2_q;
        let t3x = t2x * t2x;
        let t7x = if t3x.abs() > 1.0e-40 {
            -(0.25 - 12.0 * t0_q * (4.0 * vgsteff - t0_q) / t3x)
        } else {
            -0.25
        };
        let t4x = if t3x.abs() > 1.0e-40 {
            -(0.5 + 24.0 * t0_q * t0_q / t3x)
        } else {
            -0.5
        };
        let t5x = t7x * abulk_cv;
        (cox_wl * t4x, cox_wl * t5x * dvdseff_dvd)
    } else if xpart < 0.5 {
        // 40/60 partition
        let t2x = t2_q / 12.0;
        let t3x = if t2x.abs() > 1.0e-30 {
            0.5 * cox_wl / (t2x * t2x)
        } else {
            0.0
        };
        let vg = vgsteff;
        let v = t0_q; // = Abulk * VdseffCV
        let t4x = vg * (2.0 * v * v / 3.0 + vg * (vg - 4.0 * v / 3.0)) - 2.0 * v * v * v / 15.0;
        let t8x = 4.0 / 3.0 * vg * (vg - v) + 0.4 * v * v;
        let t5x =
            -2.0 * (-t3x * t4x) / t2x - t3x * (vg * (3.0 * vg - 8.0 * v / 3.0) + 2.0 * v * v / 3.0);
        let t6x = abulk_cv * ((-t3x * t4x) / t2x + t3x * t8x);
        (t5x, t6x * dvdseff_dvd)
    } else {
        // 50/50 partition
        (-0.5 * cgg1, -0.5 * cgd1)
    };

    // Final capacitance matrix following ngspice convention (b4ld.c lines 3890-3909)
    // qgate += Qac0 + Qsub0 - qbulk (we simplify Qsub0 ≈ 0)
    let qgate = qgate_inv + qac0 - qbulk_val;
    let qbulk = qbulk_val - qac0;

    // Gate capacitances: Cgg includes bulk charge correction
    let c_bg = cbg1; // dQbulk/dVg
    let c_bd = cbd1;
    let c_gg = cgg1 - c_bg;
    let c_gd = cgd1 - c_bd;
    let c_gb = 0.0_f64; // simplified: body bias derivative small

    let cggb = c_gg;
    let cgdb = c_gd;
    let cgsb = -(c_gg + c_gd + c_gb);

    let cbgb = c_bg;
    let cbdb = c_bd;
    let cbsb = -(c_bg + c_bd);

    // Drain cap = -(gate + bulk + source)
    let cdgb = -(cggb + cbgb + csg);
    let cddb = -(cgdb + cbdb + csd);
    let cdsb = -(cgsb + cbsb + (-(csg + csd)));

    // Junction capacitances (simplified)
    let capbd = sp.coxe * sp.weff * 0.5e-4;
    let capbs = sp.coxe * sp.weff * 0.5e-4;

    let qinv = -(qgate + qbulk + qac0);

    let (gm_final, gmbs_final, gds_final) = (gm, gmbs, gds);

    Bsim4Companion {
        ids,
        gm: gm_final,
        gds: gds_final,
        gmbs: gmbs_final,
        gbd,
        gbs,
        cbd_current,
        cbs_current,
        isub,
        igidl,
        igisl,
        ceq_d,
        ceq_bs,
        ceq_bd,
        ceq_sub,
        mode,
        vdsat,
        cggb,
        cgdb,
        cgsb,
        cbgb,
        cbdb,
        cbsb,
        cdgb,
        cddb,
        cdsb,
        capbd,
        capbs,
        qinv,
        ueff,
        vgsteff,
        vdseff,
        abulk,
    }
}

/// Poly gate depletion effect.
fn poly_depletion(vgs: f64, type_sign: f64, sp: &Bsim4SizeParam) -> f64 {
    let epsgate = EPSSI; // gate permittivity = silicon
    let t1 = 1.0e6 * CHARGE_Q * epsgate * sp.ngate / (sp.coxe * sp.coxe);
    let vgs_phi = vgs - type_sign * (sp.phi - 0.001);

    if vgs_phi <= 0.0 {
        return vgs;
    }

    let t4 = (1.0 + 2.0 * vgs_phi / t1).sqrt();
    let t2 = 2.0 * vgs_phi / (t4 + 1.0);
    let t3 = 0.5 * t2 * t2 / t1; // Vpoly

    let t7 = 1.12 - t3 - 0.05;
    let t6 = (t7 * t7 + 0.224).sqrt();
    let t5 = 1.12 - 0.5 * (t7 + t6);

    vgs - type_sign * t5
}

/// Junction diode current with threshold-based linear extrapolation.
fn junction_current(v: f64, i_sat: f64, nvtm: f64, gmin: f64, ijth: f64) -> (f64, f64) {
    if i_sat <= 0.0 {
        return (gmin, gmin * v);
    }
    let vjsm = nvtm * (ijth / i_sat + 1.0).ln();
    if v < vjsm {
        let ev = safe_exp(v / nvtm);
        let g = i_sat * ev / nvtm + gmin;
        let i = i_sat * (ev - 1.0) + gmin * v;
        (g, i)
    } else {
        let is_evjsm = i_sat * safe_exp(vjsm / nvtm);
        let g = is_evjsm / nvtm + gmin;
        let i = is_evjsm - i_sat + g * (v - vjsm) + gmin * v;
        (g, i)
    }
}

/// Stamp BSIM4 companion model into MNA matrix and RHS.
pub fn stamp_bsim4(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &Bsim4Instance,
    comp: &Bsim4Companion,
) {
    let dp = inst.drain_eff_idx();
    let g = inst.gate_idx;
    let sp = inst.source_eff_idx();
    let b = inst.bulk_idx;

    let sign = inst.model.mos_type.sign();
    let m = inst.m;

    let (_xnrm, _xrev) = if comp.mode > 0 {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };
    let mode_sign = if comp.mode > 0 { 1.0 } else { -1.0 };

    // gds between d' and s'
    crate::stamp_conductance(matrix, dp, sp, m * comp.gds);

    // gm VCCS (gate controls drain current)
    let gm_scaled = m * comp.gm * mode_sign;
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

    // gmbs (bulk controls drain current)
    let gmbs_scaled = m * comp.gmbs * mode_sign;
    if let Some(d) = dp {
        if let Some(bulk) = b {
            matrix.add(d, bulk, gmbs_scaled);
        }
        if let Some(s) = sp {
            matrix.add(d, s, -gmbs_scaled);
        }
    }
    if let Some(s) = sp {
        if let Some(bulk) = b {
            matrix.add(s, bulk, -gmbs_scaled);
        }
        matrix.add(s, s, gmbs_scaled);
    }

    // gbd between b and d'
    crate::stamp_conductance(matrix, b, dp, m * comp.gbd);
    // gbs between b and s'
    crate::stamp_conductance(matrix, b, sp, m * comp.gbs);

    // Series resistances
    let g_drain = inst.drain_conductance();
    if g_drain > 0.0 {
        crate::stamp_conductance(matrix, inst.drain_idx, dp, m * g_drain);
    }
    let g_source = inst.source_conductance();
    if g_source > 0.0 {
        crate::stamp_conductance(matrix, inst.source_idx, sp, m * g_source);
    }

    // RHS
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

/// BSIM4 voltage limiting for NR convergence.
pub fn bsim4_limit(
    vgs_new: f64,
    vds_new: f64,
    vbs_new: f64,
    vgs_old: f64,
    vds_old: f64,
    vbs_old: f64,
    vth: f64,
) -> (f64, f64, f64) {
    use crate::device_stamp::{bsim_pnjlim, fetlim};
    let vgs = fetlim(vgs_new, vgs_old, vth);
    let vds = fetlim(vds_new, vds_old, vth);
    let vbs = bsim_pnjlim(vbs_new, vbs_old);
    (vgs, vds, vbs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use thevenin_types::{Expr, ModelDef, Param};

    #[test]
    fn test_bsim4_model_defaults() {
        let m = Bsim4Model::new(MosfetType::Nmos);
        assert_eq!(m.mos_type, MosfetType::Nmos);
        assert_abs_diff_eq!(m.toxe, 3.0e-9, epsilon = 1e-25);
        assert_abs_diff_eq!(m.vth0, 0.7, epsilon = 1e-15);
        assert_abs_diff_eq!(m.u0, 0.067, epsilon = 1e-10);
        assert_abs_diff_eq!(m.vsat, 8e4, epsilon = 1e-10);
        assert_eq!(m.mob_mod, 0);
        assert_eq!(m.cap_mod, 2);
    }

    #[test]
    fn test_bsim4_model_parse() {
        let model_def = ModelDef {
            name: "NMOD".to_string(),
            kind: "NMOS".to_string(),
            params: vec![
                Param {
                    name: "VTH0".to_string(),
                    value: Expr::Num(0.423),
                },
                Param {
                    name: "TOXE".to_string(),
                    value: Expr::Num(1.85e-9),
                },
                Param {
                    name: "U0".to_string(),
                    value: Expr::Num(0.0491),
                },
                Param {
                    name: "VSAT".to_string(),
                    value: Expr::Num(124340.0),
                },
                Param {
                    name: "MOBMOD".to_string(),
                    value: Expr::Num(0.0),
                },
                Param {
                    name: "CAPMOD".to_string(),
                    value: Expr::Num(2.0),
                },
            ],
        };
        let m = Bsim4Model::from_model_def(&model_def);
        assert_eq!(m.mos_type, MosfetType::Nmos);
        assert_abs_diff_eq!(m.vth0, 0.423, epsilon = 1e-15);
        assert_abs_diff_eq!(m.toxe, 1.85e-9, epsilon = 1e-25);
        assert_abs_diff_eq!(m.u0, 0.0491, epsilon = 1e-15);
        assert_abs_diff_eq!(m.vsat, 124340.0, epsilon = 1e-10);
        assert_eq!(m.mob_mod, 0);
        assert_eq!(m.cap_mod, 2);
    }

    #[test]
    fn test_bsim4_pmos() {
        let model_def = ModelDef {
            name: "PMOD".to_string(),
            kind: "PMOS".to_string(),
            params: vec![Param {
                name: "VTH0".to_string(),
                value: Expr::Num(-0.4),
            }],
        };
        let m = Bsim4Model::from_model_def(&model_def);
        assert_eq!(m.mos_type, MosfetType::Pmos);
        assert_abs_diff_eq!(m.vth0, -0.4, epsilon = 1e-15);
    }

    #[test]
    fn test_bsim4_size_dep_params() {
        let m = Bsim4Model::new(MosfetType::Nmos);
        let sp = m.size_dep_param(10e-6, 1e-6, 1.0, 300.15);
        assert!(sp.leff > 0.0, "leff should be positive: {}", sp.leff);
        assert!(sp.weff > 0.0, "weff should be positive: {}", sp.weff);
        assert!(sp.phi > 0.0, "phi should be positive: {}", sp.phi);
        assert!(sp.u0temp > 0.0, "u0temp should be positive: {}", sp.u0temp);
        assert!(
            sp.vsattemp > 0.0,
            "vsattemp should be positive: {}",
            sp.vsattemp
        );
    }

    #[test]
    fn test_bsim4_companion_cutoff() {
        let m = Bsim4Model::new(MosfetType::Nmos);
        let sp = m.size_dep_param(10e-6, 1e-6, 1.0, 300.15);
        // Vgs well below threshold
        let comp = bsim4_companion(0.0, 1.0, 0.0, &sp, &m);
        assert!(
            comp.ids < 1e-6,
            "Ids in cutoff should be very small: {}",
            comp.ids
        );
    }

    #[test]
    fn test_bsim4_companion_saturation() {
        let m = Bsim4Model::new(MosfetType::Nmos);
        let sp = m.size_dep_param(10e-6, 1e-6, 1.0, 300.15);
        // Vgs = 1.8V, Vds = 1.8V — well in saturation
        let comp = bsim4_companion(1.8, 1.8, 0.0, &sp, &m);
        assert!(
            comp.ids > 1e-6,
            "Should have significant current in saturation: {}",
            comp.ids
        );
        assert!(comp.gm > 0.0, "gm should be positive in saturation");
        assert!(comp.gds > 0.0, "gds should be positive");
    }

    #[test]
    fn test_bsim4_companion_linear() {
        let m = Bsim4Model::new(MosfetType::Nmos);
        let sp = m.size_dep_param(10e-6, 1e-6, 1.0, 300.15);
        // Vgs = 1.8V, Vds = 0.1V — linear region
        let comp = bsim4_companion(1.8, 0.1, 0.0, &sp, &m);
        assert!(comp.ids > 0.0, "Should have current in linear region");
        assert!(comp.gds > 0.0, "gds should be positive in linear region");
        assert!(comp.gm > 0.0, "gm should be positive in linear region");
        // Linear region current should be less than saturation current
        let comp_sat = bsim4_companion(1.8, 1.8, 0.0, &sp, &m);
        assert!(comp.ids < comp_sat.ids, "Linear Ids < saturation Ids");
    }

    #[test]
    fn test_bsim4_op_convergence() {
        // Simple BSIM4 NMOS circuit: test that OP converges
        let netlist_str = r#"BSIM4 NMOS test
M1 2 1 0 0 NMOD W=10u L=0.18u
VGS 1 0 DC 1.0
VDS 2 0 DC 1.0
.model NMOD NMOS LEVEL=14 VTH0=0.7 TOXE=3e-9 U0=0.067 VSAT=8e4
.op
.end
"#;
        let netlist = thevenin_types::Netlist::parse(netlist_str).expect("parse failed");
        let result = crate::simulate_op(&netlist).expect("OP should converge");
        assert!(!result.plots.is_empty(), "Should have at least one plot");
    }
}
