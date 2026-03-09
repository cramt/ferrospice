//! BSIM3v3 MOSFET device model.
//!
//! Implements the BSIM3v3.2.4 MOSFET model matching ngspice level 8.
//! This is a physics-based short-channel MOSFET model with size-dependent
//! parameters, multiple mobility models, and advanced capacitance models.

use ferrospice_netlist::{Expr, ModelDef};

use crate::mosfet::MosfetType;

// Physical constants
const EPSOX: f64 = 3.453133e-11;
const EPSSI: f64 = 1.03594e-10;
const CHARGE_Q: f64 = 1.60219e-19;
const KBOQ: f64 = 8.617087e-5; // KB / CHARGE_Q
const EXP_THRESHOLD: f64 = 34.0;
const MAX_EXP: f64 = 5.834617425e14; // exp(34)
const MIN_EXP: f64 = 1.713908431e-15; // exp(-34)

const TEMP_DEFAULT: f64 = 300.15;

fn safe_exp(x: f64) -> f64 {
    if x > EXP_THRESHOLD {
        MAX_EXP
    } else if x < -EXP_THRESHOLD {
        MIN_EXP
    } else {
        x.exp()
    }
}

/// BSIM3v3 model parameters (from .model card).
#[derive(Debug, Clone)]
pub struct Bsim3Model {
    pub mos_type: MosfetType,

    // Mode selection
    pub mob_mod: i32,
    pub cap_mod: i32,
    pub noi_mod: i32,
    pub bin_unit: i32,

    // Oxide
    pub tox: f64,
    pub toxm: f64,

    // Threshold voltage
    pub vth0: f64,
    pub k1: f64,
    pub k2: f64,
    pub k3: f64,
    pub k3b: f64,
    pub w0: f64,
    pub nlx: f64,

    // SCE/DIBL
    pub dvt0: f64,
    pub dvt1: f64,
    pub dvt2: f64,
    pub dvt0w: f64,
    pub dvt1w: f64,
    pub dvt2w: f64,
    pub dsub: f64,
    pub eta0: f64,
    pub etab: f64,

    // Subthreshold
    pub voff: f64,
    pub nfactor: f64,
    pub cdsc: f64,
    pub cdscb: f64,
    pub cdscd: f64,
    pub cit: f64,

    // Doping
    pub nch: f64,
    pub npeak: f64,
    pub ngate: f64,
    pub nsub: f64,
    pub xj: f64,
    pub vfb: f64,
    pub vbx: f64,
    pub vbm: f64,
    pub xt: f64,
    pub gamma1: f64,
    pub gamma2: f64,

    // Mobility
    pub u0: f64,
    pub ua: f64,
    pub ub: f64,
    pub uc: f64,
    pub ute: f64,
    pub ua1: f64,
    pub ub1: f64,
    pub uc1: f64,

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

    // Series resistance
    pub rdsw: f64,
    pub rsh: f64,
    pub prwg: f64,
    pub prwb: f64,
    pub prt: f64,
    pub wr: f64,

    // Width/length effects
    pub dwg: f64,
    pub dwb: f64,
    pub b0: f64,
    pub b1: f64,

    // Substrate current
    pub alpha0: f64,
    pub alpha1: f64,
    pub beta0: f64,

    // Temperature
    pub tnom: f64,
    pub kt1: f64,
    pub kt1l: f64,
    pub kt2: f64,

    // Junction
    pub js: f64,
    pub jsw: f64,
    pub pb: f64,
    pub mj: f64,
    pub pbsw: f64,
    pub mjsw: f64,
    pub cj: f64,
    pub cjsw: f64,
    pub pbswg: f64,
    pub mjswg: f64,
    pub cjswg: f64,
    pub nj: f64,
    pub xti: f64,
    pub ijth: f64,

    // Junction temperature coefficients
    pub tpb: f64,
    pub tcj: f64,
    pub tpbsw: f64,
    pub tcjsw: f64,
    pub tpbswg: f64,
    pub tcjswg: f64,

    // Capacitance CV model
    pub cgsl: f64,
    pub cgdl: f64,
    pub ckappa: f64,
    pub cf: f64,
    pub clc: f64,
    pub cle: f64,
    pub elm: f64,
    pub vfbcv: f64,
    pub acde: f64,
    pub moin: f64,
    pub noff: f64,
    pub voffcv: f64,
    pub xpart: f64,

    // Gate overlap
    pub cgso: f64,
    pub cgdo: f64,
    pub cgbo: f64,

    // L/W binning
    pub lint: f64,
    pub ll: f64,
    pub lln: f64,
    pub lw: f64,
    pub lwn: f64,
    pub lwl: f64,
    pub wint: f64,
    pub wl: f64,
    pub wln: f64,
    pub ww: f64,
    pub wwn: f64,
    pub wwl: f64,
    pub xl: f64,
    pub xw: f64,

    // DLC/DWC
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

    // "Given" flags for parameters that have defaults derived from other params
    pub vth0_given: bool,
    pub vfb_given: bool,
    pub k1_given: bool,
    pub k2_given: bool,
    pub npeak_given: bool,
    pub gamma1_given: bool,
    pub gamma2_given: bool,
    pub vbx_given: bool,
    pub dlc_given: bool,
    pub dwc_given: bool,

    // Length-dependent parameters (l prefix)
    pub l_vth0: f64,
    pub l_k1: f64,
    pub l_k2: f64,
    pub l_k3: f64,
    pub l_k3b: f64,
    pub l_w0: f64,
    pub l_nlx: f64,
    pub l_dvt0: f64,
    pub l_dvt1: f64,
    pub l_dvt2: f64,
    pub l_dvt0w: f64,
    pub l_dvt1w: f64,
    pub l_dvt2w: f64,
    pub l_dsub: f64,
    pub l_eta0: f64,
    pub l_etab: f64,
    pub l_voff: f64,
    pub l_nfactor: f64,
    pub l_cdsc: f64,
    pub l_cdscb: f64,
    pub l_cdscd: f64,
    pub l_cit: f64,
    pub l_u0: f64,
    pub l_ua: f64,
    pub l_ub: f64,
    pub l_uc: f64,
    pub l_ute: f64,
    pub l_ua1: f64,
    pub l_ub1: f64,
    pub l_uc1: f64,
    pub l_vsat: f64,
    pub l_a0: f64,
    pub l_ags: f64,
    pub l_a1: f64,
    pub l_a2: f64,
    pub l_at: f64,
    pub l_keta: f64,
    pub l_pclm: f64,
    pub l_pdibl1: f64,
    pub l_pdibl2: f64,
    pub l_pdiblb: f64,
    pub l_pscbe1: f64,
    pub l_pscbe2: f64,
    pub l_pvag: f64,
    pub l_delta: f64,
    pub l_rdsw: f64,
    pub l_prwg: f64,
    pub l_prwb: f64,
    pub l_dwg: f64,
    pub l_dwb: f64,
    pub l_b0: f64,
    pub l_b1: f64,
    pub l_alpha0: f64,
    pub l_alpha1: f64,
    pub l_beta0: f64,
    pub l_kt1: f64,
    pub l_kt1l: f64,
    pub l_kt2: f64,
    pub l_vfbcv: f64,
    pub l_acde: f64,
    pub l_moin: f64,
    pub l_noff: f64,
    pub l_voffcv: f64,
    pub l_cgsl: f64,
    pub l_cgdl: f64,
    pub l_ckappa: f64,
    pub l_cf: f64,
    pub l_clc: f64,
    pub l_cle: f64,
    pub l_npeak: f64,
    pub l_nsub: f64,
    pub l_ngate: f64,
    pub l_xj: f64,

    // Width-dependent parameters (w prefix)
    pub w_vth0: f64,
    pub w_k1: f64,
    pub w_k2: f64,
    pub w_k3: f64,
    pub w_k3b: f64,
    pub w_w0: f64,
    pub w_nlx: f64,
    pub w_dvt0: f64,
    pub w_dvt1: f64,
    pub w_dvt2: f64,
    pub w_dvt0w: f64,
    pub w_dvt1w: f64,
    pub w_dvt2w: f64,
    pub w_dsub: f64,
    pub w_eta0: f64,
    pub w_etab: f64,
    pub w_voff: f64,
    pub w_nfactor: f64,
    pub w_cdsc: f64,
    pub w_cdscb: f64,
    pub w_cdscd: f64,
    pub w_cit: f64,
    pub w_u0: f64,
    pub w_ua: f64,
    pub w_ub: f64,
    pub w_uc: f64,
    pub w_ute: f64,
    pub w_ua1: f64,
    pub w_ub1: f64,
    pub w_uc1: f64,
    pub w_vsat: f64,
    pub w_a0: f64,
    pub w_ags: f64,
    pub w_a1: f64,
    pub w_a2: f64,
    pub w_at: f64,
    pub w_keta: f64,
    pub w_pclm: f64,
    pub w_pdibl1: f64,
    pub w_pdibl2: f64,
    pub w_pdiblb: f64,
    pub w_pscbe1: f64,
    pub w_pscbe2: f64,
    pub w_pvag: f64,
    pub w_delta: f64,
    pub w_rdsw: f64,
    pub w_prwg: f64,
    pub w_prwb: f64,
    pub w_dwg: f64,
    pub w_dwb: f64,
    pub w_b0: f64,
    pub w_b1: f64,
    pub w_alpha0: f64,
    pub w_alpha1: f64,
    pub w_beta0: f64,
    pub w_kt1: f64,
    pub w_kt1l: f64,
    pub w_kt2: f64,
    pub w_vfbcv: f64,
    pub w_acde: f64,
    pub w_moin: f64,
    pub w_noff: f64,
    pub w_voffcv: f64,
    pub w_cgsl: f64,
    pub w_cgdl: f64,
    pub w_ckappa: f64,
    pub w_cf: f64,
    pub w_clc: f64,
    pub w_cle: f64,
    pub w_npeak: f64,
    pub w_nsub: f64,
    pub w_ngate: f64,
    pub w_xj: f64,

    // Cross-term parameters (p prefix)
    pub p_vth0: f64,
    pub p_k1: f64,
    pub p_k2: f64,
    pub p_k3: f64,
    pub p_k3b: f64,
    pub p_w0: f64,
    pub p_nlx: f64,
    pub p_dvt0: f64,
    pub p_dvt1: f64,
    pub p_dvt2: f64,
    pub p_dvt0w: f64,
    pub p_dvt1w: f64,
    pub p_dvt2w: f64,
    pub p_dsub: f64,
    pub p_eta0: f64,
    pub p_etab: f64,
    pub p_voff: f64,
    pub p_nfactor: f64,
    pub p_cdsc: f64,
    pub p_cdscb: f64,
    pub p_cdscd: f64,
    pub p_cit: f64,
    pub p_u0: f64,
    pub p_ua: f64,
    pub p_ub: f64,
    pub p_uc: f64,
    pub p_ute: f64,
    pub p_ua1: f64,
    pub p_ub1: f64,
    pub p_uc1: f64,
    pub p_vsat: f64,
    pub p_a0: f64,
    pub p_ags: f64,
    pub p_a1: f64,
    pub p_a2: f64,
    pub p_at: f64,
    pub p_keta: f64,
    pub p_pclm: f64,
    pub p_pdibl1: f64,
    pub p_pdibl2: f64,
    pub p_pdiblb: f64,
    pub p_pscbe1: f64,
    pub p_pscbe2: f64,
    pub p_pvag: f64,
    pub p_delta: f64,
    pub p_rdsw: f64,
    pub p_prwg: f64,
    pub p_prwb: f64,
    pub p_dwg: f64,
    pub p_dwb: f64,
    pub p_b0: f64,
    pub p_b1: f64,
    pub p_alpha0: f64,
    pub p_alpha1: f64,
    pub p_beta0: f64,
    pub p_kt1: f64,
    pub p_kt1l: f64,
    pub p_kt2: f64,
    pub p_vfbcv: f64,
    pub p_acde: f64,
    pub p_moin: f64,
    pub p_noff: f64,
    pub p_voffcv: f64,
    pub p_cgsl: f64,
    pub p_cgdl: f64,
    pub p_ckappa: f64,
    pub p_cf: f64,
    pub p_clc: f64,
    pub p_cle: f64,
    pub p_npeak: f64,
    pub p_nsub: f64,
    pub p_ngate: f64,
    pub p_xj: f64,
}

/// Size-dependent parameters calculated for a specific (W, L) combination.
#[derive(Debug, Clone)]
pub struct Bsim3SizeParam {
    // Effective dimensions
    pub leff: f64,
    pub weff: f64,
    pub leff_cv: f64,
    pub weff_cv: f64,

    // Binned core parameters
    pub vth0: f64,
    pub k1: f64,
    pub k2: f64,
    pub k3: f64,
    pub k3b: f64,
    pub w0: f64,
    pub nlx: f64,
    pub dvt0: f64,
    pub dvt1: f64,
    pub dvt2: f64,
    pub dvt0w: f64,
    pub dvt1w: f64,
    pub dvt2w: f64,
    pub dsub: f64,
    pub eta0: f64,
    pub etab: f64,
    pub voff: f64,
    pub nfactor: f64,
    pub cdsc: f64,
    pub cdscb: f64,
    pub cdscd: f64,
    pub cit: f64,
    pub u0: f64,
    pub ua: f64,
    pub ub: f64,
    pub uc: f64,
    pub ute: f64,
    pub ua1: f64,
    pub ub1: f64,
    pub uc1: f64,
    pub vsat: f64,
    pub a0: f64,
    pub ags: f64,
    pub a1: f64,
    pub a2: f64,
    pub at: f64,
    pub keta: f64,
    pub pclm: f64,
    pub pdibl1: f64,
    pub pdibl2: f64,
    pub pdiblb: f64,
    pub pscbe1: f64,
    pub pscbe2: f64,
    pub pvag: f64,
    pub delta: f64,
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
    pub kt1: f64,
    pub kt1l: f64,
    pub kt2: f64,
    pub npeak: f64,
    pub nsub: f64,
    pub ngate: f64,
    pub xj: f64,

    // CV model params
    pub cgsl: f64,
    pub cgdl: f64,
    pub ckappa: f64,
    pub cf: f64,
    pub clc: f64,
    pub cle: f64,
    pub vfbcv: f64,
    pub acde: f64,
    pub moin: f64,
    pub noff: f64,
    pub voffcv: f64,

    // Derived quantities
    pub phi: f64,
    pub sqrt_phi: f64,
    pub phis3: f64,
    pub xdep0: f64,
    pub sqrt_xdep0: f64,
    pub vbi: f64,
    pub cdep0: f64,
    pub ldeb: f64,
    pub litl: f64,
    pub vbsc: f64,
    pub k1ox: f64,
    pub k2ox: f64,
    pub theta0vb0: f64,
    pub theta_rout: f64,
    pub vfbzb: f64,

    // Temperature-adjusted
    pub u0temp: f64,
    pub vsattemp: f64,
    pub rds0: f64,

    // Overlap caps
    pub cgso_eff: f64,
    pub cgdo_eff: f64,
    pub cgbo_eff: f64,

    // abulkCV factor
    pub abulk_cv_factor: f64,
    pub tconst: f64,
}

/// NR companion result for BSIM3.
#[derive(Debug, Clone)]
pub struct Bsim3Companion {
    pub ids: f64,
    pub gm: f64,
    pub gds: f64,
    pub gmbs: f64,
    pub gbd: f64,
    pub gbs: f64,
    pub cbd_current: f64,
    pub cbs_current: f64,
    pub isub: f64,
    pub gbg: f64,
    pub gbb: f64,
    pub gbd_sub: f64,
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
}

/// Resolved BSIM3 instance with node indices.
#[derive(Debug, Clone)]
pub struct Bsim3Instance {
    pub name: String,
    pub drain_idx: Option<usize>,
    pub gate_idx: Option<usize>,
    pub source_idx: Option<usize>,
    pub bulk_idx: Option<usize>,
    pub drain_prime_idx: Option<usize>,
    pub source_prime_idx: Option<usize>,
    pub w: f64,
    pub l: f64,
    pub ad: f64,
    pub as_: f64,
    pub pd: f64,
    pub ps: f64,
    pub nrd: f64,
    pub nrs: f64,
    pub m: f64,
    pub model: Bsim3Model,
    pub size_params: Bsim3SizeParam,
    // Operating point values for AC
    pub vth0_inst: f64,
    pub vfb_inst: f64,
    pub vfbzb_inst: f64,
}

impl Bsim3Model {
    pub fn new(mos_type: MosfetType) -> Self {
        let vth0_default = match mos_type {
            MosfetType::Nmos => 0.7,
            MosfetType::Pmos => -0.7,
        };
        let u0_default = match mos_type {
            MosfetType::Nmos => 0.067, // 670 cm²/Vs in m²/Vs (ngspice converts >1 by /1e4)
            MosfetType::Pmos => 0.025,
        };
        Self {
            mos_type,
            mob_mod: 1,
            cap_mod: 3,
            noi_mod: 1,
            bin_unit: 1,
            tox: 150e-10,
            toxm: 150e-10,
            vth0: vth0_default,
            k1: 0.5,
            k2: 0.0,
            k3: 80.0,
            k3b: 0.0,
            w0: 2.5e-6,
            nlx: 1.74e-7,
            dvt0: 2.2,
            dvt1: 0.53,
            dvt2: -0.032,
            dvt0w: 0.0,
            dvt1w: 5.3e6,
            dvt2w: -0.032,
            dsub: 0.56,
            eta0: 0.08,
            etab: -0.07,
            voff: -0.08,
            nfactor: 1.0,
            cdsc: 2.4e-4,
            cdscb: 0.0,
            cdscd: 0.0,
            cit: 0.0,
            nch: 1.7e17,
            npeak: 1.7e17,
            ngate: 0.0,
            nsub: 6e16,
            xj: 1.5e-7,
            vfb: -1.0,
            vbx: 0.0,
            vbm: -3.0,
            xt: 1.55e-7,
            gamma1: 0.0,
            gamma2: 0.0,
            u0: u0_default,
            ua: 2.25e-9,
            ub: 5.87e-19,
            uc: -4.65e-11,
            ute: -1.5,
            ua1: 4.31e-9,
            ub1: -7.61e-18,
            uc1: -5.6e-11,
            vsat: 8e4,
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
            pscbe2: 1e-5,
            pvag: 0.0,
            delta: 0.01,
            rdsw: 0.0,
            rsh: 0.0,
            prwg: 0.0,
            prwb: 0.0,
            prt: 0.0,
            wr: 1.0,
            dwg: 0.0,
            dwb: 0.0,
            b0: 0.0,
            b1: 0.0,
            alpha0: 0.0,
            alpha1: 0.0,
            beta0: 30.0,
            tnom: 27.0,
            kt1: -0.11,
            kt1l: 0.0,
            kt2: 0.022,
            js: 1e-4,
            jsw: 0.0,
            pb: 1.0,
            mj: 0.5,
            pbsw: 1.0,
            mjsw: 0.33,
            cj: 5e-4,
            cjsw: 5e-10,
            pbswg: 1.0,
            mjswg: 0.33,
            cjswg: 5e-10,
            nj: 1.0,
            xti: 3.0,
            ijth: 0.1,
            tpb: 0.0,
            tcj: 0.0,
            tpbsw: 0.0,
            tcjsw: 0.0,
            tpbswg: 0.0,
            tcjswg: 0.0,
            cgsl: 0.0,
            cgdl: 0.0,
            ckappa: 0.6,
            cf: 0.0,
            clc: 1e-7,
            cle: 0.6,
            elm: 5.0,
            vfbcv: -1.0,
            acde: 1.0,
            moin: 15.0,
            noff: 1.0,
            voffcv: 0.0,
            xpart: 0.0,
            cgso: 0.0,
            cgdo: 0.0,
            cgbo: 0.0,
            lint: 0.0,
            ll: 0.0,
            lln: 1.0,
            lw: 0.0,
            lwn: 1.0,
            lwl: 0.0,
            wint: 0.0,
            wl: 0.0,
            wln: 1.0,
            ww: 0.0,
            wwn: 1.0,
            wwl: 0.0,
            xl: 0.0,
            xw: 0.0,
            dlc: 0.0,
            dwc: 0.0,
            kf: 0.0,
            af: 1.0,
            ef: 1.0,
            noia: 1e20,
            noib: 5e4,
            noic: -1.4e-12,
            em: 4.1e7,
            vth0_given: false,
            vfb_given: false,
            k1_given: false,
            k2_given: false,
            npeak_given: false,
            gamma1_given: false,
            gamma2_given: false,
            vbx_given: false,
            dlc_given: false,
            dwc_given: false,
            // All L/W/P params default to 0
            l_vth0: 0.0,
            l_k1: 0.0,
            l_k2: 0.0,
            l_k3: 0.0,
            l_k3b: 0.0,
            l_w0: 0.0,
            l_nlx: 0.0,
            l_dvt0: 0.0,
            l_dvt1: 0.0,
            l_dvt2: 0.0,
            l_dvt0w: 0.0,
            l_dvt1w: 0.0,
            l_dvt2w: 0.0,
            l_dsub: 0.0,
            l_eta0: 0.0,
            l_etab: 0.0,
            l_voff: 0.0,
            l_nfactor: 0.0,
            l_cdsc: 0.0,
            l_cdscb: 0.0,
            l_cdscd: 0.0,
            l_cit: 0.0,
            l_u0: 0.0,
            l_ua: 0.0,
            l_ub: 0.0,
            l_uc: 0.0,
            l_ute: 0.0,
            l_ua1: 0.0,
            l_ub1: 0.0,
            l_uc1: 0.0,
            l_vsat: 0.0,
            l_a0: 0.0,
            l_ags: 0.0,
            l_a1: 0.0,
            l_a2: 0.0,
            l_at: 0.0,
            l_keta: 0.0,
            l_pclm: 0.0,
            l_pdibl1: 0.0,
            l_pdibl2: 0.0,
            l_pdiblb: 0.0,
            l_pscbe1: 0.0,
            l_pscbe2: 0.0,
            l_pvag: 0.0,
            l_delta: 0.0,
            l_rdsw: 0.0,
            l_prwg: 0.0,
            l_prwb: 0.0,
            l_dwg: 0.0,
            l_dwb: 0.0,
            l_b0: 0.0,
            l_b1: 0.0,
            l_alpha0: 0.0,
            l_alpha1: 0.0,
            l_beta0: 0.0,
            l_kt1: 0.0,
            l_kt1l: 0.0,
            l_kt2: 0.0,
            l_vfbcv: 0.0,
            l_acde: 0.0,
            l_moin: 0.0,
            l_noff: 0.0,
            l_voffcv: 0.0,
            l_cgsl: 0.0,
            l_cgdl: 0.0,
            l_ckappa: 0.0,
            l_cf: 0.0,
            l_clc: 0.0,
            l_cle: 0.0,
            l_npeak: 0.0,
            l_nsub: 0.0,
            l_ngate: 0.0,
            l_xj: 0.0,
            w_vth0: 0.0,
            w_k1: 0.0,
            w_k2: 0.0,
            w_k3: 0.0,
            w_k3b: 0.0,
            w_w0: 0.0,
            w_nlx: 0.0,
            w_dvt0: 0.0,
            w_dvt1: 0.0,
            w_dvt2: 0.0,
            w_dvt0w: 0.0,
            w_dvt1w: 0.0,
            w_dvt2w: 0.0,
            w_dsub: 0.0,
            w_eta0: 0.0,
            w_etab: 0.0,
            w_voff: 0.0,
            w_nfactor: 0.0,
            w_cdsc: 0.0,
            w_cdscb: 0.0,
            w_cdscd: 0.0,
            w_cit: 0.0,
            w_u0: 0.0,
            w_ua: 0.0,
            w_ub: 0.0,
            w_uc: 0.0,
            w_ute: 0.0,
            w_ua1: 0.0,
            w_ub1: 0.0,
            w_uc1: 0.0,
            w_vsat: 0.0,
            w_a0: 0.0,
            w_ags: 0.0,
            w_a1: 0.0,
            w_a2: 0.0,
            w_at: 0.0,
            w_keta: 0.0,
            w_pclm: 0.0,
            w_pdibl1: 0.0,
            w_pdibl2: 0.0,
            w_pdiblb: 0.0,
            w_pscbe1: 0.0,
            w_pscbe2: 0.0,
            w_pvag: 0.0,
            w_delta: 0.0,
            w_rdsw: 0.0,
            w_prwg: 0.0,
            w_prwb: 0.0,
            w_dwg: 0.0,
            w_dwb: 0.0,
            w_b0: 0.0,
            w_b1: 0.0,
            w_alpha0: 0.0,
            w_alpha1: 0.0,
            w_beta0: 0.0,
            w_kt1: 0.0,
            w_kt1l: 0.0,
            w_kt2: 0.0,
            w_vfbcv: 0.0,
            w_acde: 0.0,
            w_moin: 0.0,
            w_noff: 0.0,
            w_voffcv: 0.0,
            w_cgsl: 0.0,
            w_cgdl: 0.0,
            w_ckappa: 0.0,
            w_cf: 0.0,
            w_clc: 0.0,
            w_cle: 0.0,
            w_npeak: 0.0,
            w_nsub: 0.0,
            w_ngate: 0.0,
            w_xj: 0.0,
            p_vth0: 0.0,
            p_k1: 0.0,
            p_k2: 0.0,
            p_k3: 0.0,
            p_k3b: 0.0,
            p_w0: 0.0,
            p_nlx: 0.0,
            p_dvt0: 0.0,
            p_dvt1: 0.0,
            p_dvt2: 0.0,
            p_dvt0w: 0.0,
            p_dvt1w: 0.0,
            p_dvt2w: 0.0,
            p_dsub: 0.0,
            p_eta0: 0.0,
            p_etab: 0.0,
            p_voff: 0.0,
            p_nfactor: 0.0,
            p_cdsc: 0.0,
            p_cdscb: 0.0,
            p_cdscd: 0.0,
            p_cit: 0.0,
            p_u0: 0.0,
            p_ua: 0.0,
            p_ub: 0.0,
            p_uc: 0.0,
            p_ute: 0.0,
            p_ua1: 0.0,
            p_ub1: 0.0,
            p_uc1: 0.0,
            p_vsat: 0.0,
            p_a0: 0.0,
            p_ags: 0.0,
            p_a1: 0.0,
            p_a2: 0.0,
            p_at: 0.0,
            p_keta: 0.0,
            p_pclm: 0.0,
            p_pdibl1: 0.0,
            p_pdibl2: 0.0,
            p_pdiblb: 0.0,
            p_pscbe1: 0.0,
            p_pscbe2: 0.0,
            p_pvag: 0.0,
            p_delta: 0.0,
            p_rdsw: 0.0,
            p_prwg: 0.0,
            p_prwb: 0.0,
            p_dwg: 0.0,
            p_dwb: 0.0,
            p_b0: 0.0,
            p_b1: 0.0,
            p_alpha0: 0.0,
            p_alpha1: 0.0,
            p_beta0: 0.0,
            p_kt1: 0.0,
            p_kt1l: 0.0,
            p_kt2: 0.0,
            p_vfbcv: 0.0,
            p_acde: 0.0,
            p_moin: 0.0,
            p_noff: 0.0,
            p_voffcv: 0.0,
            p_cgsl: 0.0,
            p_cgdl: 0.0,
            p_ckappa: 0.0,
            p_cf: 0.0,
            p_clc: 0.0,
            p_cle: 0.0,
            p_npeak: 0.0,
            p_nsub: 0.0,
            p_ngate: 0.0,
            p_xj: 0.0,
        }
    }

    /// Parse BSIM3 model from a .model definition.
    #[expect(clippy::too_many_lines)]
    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let mos_type = if model_def.kind.to_uppercase().contains("PMOS") {
            MosfetType::Pmos
        } else {
            MosfetType::Nmos
        };
        let mut m = Self::new(mos_type);
        for p in &model_def.params {
            if let Expr::Num(v) = &p.value {
                let v = *v;
                match p.name.to_uppercase().as_str() {
                    "MOBMOD" => m.mob_mod = v as i32,
                    "CAPMOD" => m.cap_mod = v as i32,
                    "NOIMOD" => m.noi_mod = v as i32,
                    "BINUNIT" => m.bin_unit = v as i32,
                    "TOX" => m.tox = v,
                    "TOXM" => m.toxm = v,
                    "VTH0" | "VTNM" => {
                        m.vth0 = v;
                        m.vth0_given = true;
                    }
                    "K1" => {
                        m.k1 = v;
                        m.k1_given = true;
                    }
                    "K2" => {
                        m.k2 = v;
                        m.k2_given = true;
                    }
                    "K3" => m.k3 = v,
                    "K3B" => m.k3b = v,
                    "W0" => m.w0 = v,
                    "NLX" => m.nlx = v,
                    "DVT0" => m.dvt0 = v,
                    "DVT1" => m.dvt1 = v,
                    "DVT2" => m.dvt2 = v,
                    "DVT0W" => m.dvt0w = v,
                    "DVT1W" => m.dvt1w = v,
                    "DVT2W" => m.dvt2w = v,
                    "DSUB" => m.dsub = v,
                    "ETA0" => m.eta0 = v,
                    "ETAB" => m.etab = v,
                    "VOFF" => m.voff = v,
                    "NFACTOR" => m.nfactor = v,
                    "CDSC" => m.cdsc = v,
                    "CDSCB" => m.cdscb = v,
                    "CDSCD" => m.cdscd = v,
                    "CIT" => m.cit = v,
                    "NCH" => {
                        m.nch = v;
                        if v > 1e20 {
                            m.nch *= 1e-6;
                        }
                    }
                    "NPEAK" | "NPEAM" => {
                        m.npeak = v;
                        m.npeak_given = true;
                        if v > 1e20 {
                            m.npeak *= 1e-6;
                        }
                    }
                    "NGATE" => {
                        m.ngate = v;
                        if v > 1.000001e24 {
                            m.ngate *= 1e-6;
                        }
                    }
                    "NSUB" => m.nsub = v,
                    "XJ" => m.xj = v,
                    "VFB" => {
                        m.vfb = v;
                        m.vfb_given = true;
                    }
                    "VBX" => {
                        m.vbx = v;
                        m.vbx_given = true;
                    }
                    "VBM" => m.vbm = v,
                    "XT" => m.xt = v,
                    "GAMMA1" => {
                        m.gamma1 = v;
                        m.gamma1_given = true;
                    }
                    "GAMMA2" => {
                        m.gamma2 = v;
                        m.gamma2_given = true;
                    }
                    "U0" | "UO" => m.u0 = v,
                    "UA" => m.ua = v,
                    "UB" => m.ub = v,
                    "UC" => m.uc = v,
                    "UTE" => m.ute = v,
                    "UA1" => m.ua1 = v,
                    "UB1" => m.ub1 = v,
                    "UC1" => m.uc1 = v,
                    "VSAT" => m.vsat = v,
                    "A0" => m.a0 = v,
                    "AGS" => m.ags = v,
                    "A1" => m.a1 = v,
                    "A2" => m.a2 = v,
                    "AT" => m.at = v,
                    "KETA" => m.keta = v,
                    "PCLM" => m.pclm = v,
                    "PDIBLC1" | "PDIBL1" => m.pdibl1 = v,
                    "PDIBLC2" | "PDIBL2" => m.pdibl2 = v,
                    "PDIBLCB" | "PDIBLB" => m.pdiblb = v,
                    "PSCBE1" => m.pscbe1 = v,
                    "PSCBE2" => m.pscbe2 = v,
                    "PVAG" => m.pvag = v,
                    "DELTA" => m.delta = v,
                    "RDSW" => m.rdsw = v,
                    "RSH" => m.rsh = v,
                    "PRWG" => m.prwg = v,
                    "PRWB" => m.prwb = v,
                    "PRT" => m.prt = v,
                    "WR" => m.wr = v,
                    "DWG" => m.dwg = v,
                    "DWB" => m.dwb = v,
                    "B0" => m.b0 = v,
                    "B1" => m.b1 = v,
                    "ALPHA0" => m.alpha0 = v,
                    "ALPHA1" => m.alpha1 = v,
                    "BETA0" => m.beta0 = v,
                    "TNOM" => m.tnom = v,
                    "KT1" => m.kt1 = v,
                    "KT1L" => m.kt1l = v,
                    "KT2" => m.kt2 = v,
                    "JS" => m.js = v,
                    "JSW" => m.jsw = v,
                    "PB" | "TLEVM" => m.pb = v,
                    "MJ" => m.mj = v,
                    "PBSW" => m.pbsw = v,
                    "MJSW" => m.mjsw = v,
                    "CJ" => m.cj = v,
                    "CJSW" => m.cjsw = v,
                    "PBSWG" => m.pbswg = v,
                    "MJSWG" => m.mjswg = v,
                    "CJSWG" => m.cjswg = v,
                    "NJ" => m.nj = v,
                    "XTI" => m.xti = v,
                    "IJTH" => m.ijth = v,
                    "TPB" => m.tpb = v,
                    "TCJ" => m.tcj = v,
                    "TPBSW" => m.tpbsw = v,
                    "TCJSW" => m.tcjsw = v,
                    "TPBSWG" => m.tpbswg = v,
                    "TCJSWG" => m.tcjswg = v,
                    "CGSL" => m.cgsl = v,
                    "CGDL" => m.cgdl = v,
                    "CKAPPA" => m.ckappa = v,
                    "CF" => m.cf = v,
                    "CLC" => m.clc = v,
                    "CLE" => m.cle = v,
                    "ELM" => m.elm = v,
                    "VFBCV" => m.vfbcv = v,
                    "ACDE" => m.acde = v,
                    "MOIN" => m.moin = v,
                    "NOFF" => m.noff = v,
                    "VOFFCV" => m.voffcv = v,
                    "XPART" => m.xpart = v,
                    "CGSO" => m.cgso = v,
                    "CGDO" => m.cgdo = v,
                    "CGBO" => m.cgbo = v,
                    "LINT" => m.lint = v,
                    "LL" => m.ll = v,
                    "LLN" => m.lln = v,
                    "LW" => m.lw = v,
                    "LWN" => m.lwn = v,
                    "LWL" => m.lwl = v,
                    "WINT" => m.wint = v,
                    "WL" => m.wl = v,
                    "WLN" => m.wln = v,
                    "WW" => m.ww = v,
                    "WWN" => m.wwn = v,
                    "WWL" => m.wwl = v,
                    "XL" => m.xl = v,
                    "XW" => m.xw = v,
                    "DLC" => {
                        m.dlc = v;
                        m.dlc_given = true;
                    }
                    "DWC" => {
                        m.dwc = v;
                        m.dwc_given = true;
                    }
                    "KF" => m.kf = v,
                    "AF" => m.af = v,
                    "EF" => m.ef = v,
                    "NOIA" => m.noia = v,
                    "NOIB" => m.noib = v,
                    "NOIC" => m.noic = v,
                    "EM" => m.em = v,
                    // L/W/P dependent params
                    "LVTH0" => m.l_vth0 = v,
                    "WVTH0" => m.w_vth0 = v,
                    "PVTH0" => m.p_vth0 = v,
                    "LK1" => m.l_k1 = v,
                    "WK1" => m.w_k1 = v,
                    "PK1" => m.p_k1 = v,
                    "LK2" => m.l_k2 = v,
                    "WK2" => m.w_k2 = v,
                    "PK2" => m.p_k2 = v,
                    "LK3" => m.l_k3 = v,
                    "WK3" => m.w_k3 = v,
                    "PK3" => m.p_k3 = v,
                    "LK3B" => m.l_k3b = v,
                    "WK3B" => m.w_k3b = v,
                    "PK3B" => m.p_k3b = v,
                    "LETA0" => m.l_eta0 = v,
                    "WETA0" => m.w_eta0 = v,
                    "PETA0" => m.p_eta0 = v,
                    "LETAB" => m.l_etab = v,
                    "WETAB" => m.w_etab = v,
                    "PETAB" => m.p_etab = v,
                    "LU0" => m.l_u0 = v,
                    "WU0" => m.w_u0 = v,
                    "PU0" => m.p_u0 = v,
                    "LUA" => m.l_ua = v,
                    "WUA" => m.w_ua = v,
                    "PUA" => m.p_ua = v,
                    "LUB" => m.l_ub = v,
                    "WUB" => m.w_ub = v,
                    "PUB" => m.p_ub = v,
                    "LUC" => m.l_uc = v,
                    "WUC" => m.w_uc = v,
                    "PUC" => m.p_uc = v,
                    "LVSAT" => m.l_vsat = v,
                    "WVSAT" => m.w_vsat = v,
                    "PVSAT" => m.p_vsat = v,
                    "LVOFF" => m.l_voff = v,
                    "WVOFF" => m.w_voff = v,
                    "PVOFF" => m.p_voff = v,
                    "LNFACTOR" => m.l_nfactor = v,
                    "WNFACTOR" => m.w_nfactor = v,
                    "PNFACTOR" => m.p_nfactor = v,
                    "LRDSW" => m.l_rdsw = v,
                    "WRDSW" => m.w_rdsw = v,
                    "PRDSW" => m.p_rdsw = v,
                    "LPCLM" => m.l_pclm = v,
                    "WPCLM" => m.w_pclm = v,
                    "PPCLM" => m.p_pclm = v,
                    "LDELTA" => m.l_delta = v,
                    "WDELTA" => m.w_delta = v,
                    "PDELTA" => m.p_delta = v,
                    _ => {} // ignore unknown params
                }
            }
        }

        // Set derived defaults
        if m.cf == 0.0 {
            m.cf = 2.0 * EPSOX / std::f64::consts::PI * (1.0 + 0.4e-6 / m.tox).ln();
        }
        if !m.dlc_given {
            m.dlc = m.lint;
        }
        if !m.dwc_given {
            m.dwc = m.wint;
        }

        m
    }

    /// Number of internal nodes needed (drain_prime and source_prime always exist in BSIM3).
    /// Series resistance from RDSW is internal to the model, but
    /// external RD/RS (from RSH*NRD/NRS) create internal nodes.
    pub fn internal_node_count(&self, nrd: f64, nrs: f64) -> usize {
        let mut count = 0;
        if self.rsh > 0.0 && nrd > 0.0 {
            count += 1;
        }
        if self.rsh > 0.0 && nrs > 0.0 {
            count += 1;
        }
        count
    }

    /// Compute the oxide capacitance.
    pub fn cox(&self) -> f64 {
        EPSOX / self.tox
    }

    /// Compute size-dependent parameters for a given (W, L) at temperature T (Kelvin).
    #[expect(clippy::too_many_lines)]
    pub fn size_dep_param(&self, inst_w: f64, inst_l: f64, temp: f64) -> Bsim3SizeParam {
        let tnom_k = self.tnom + 273.15;
        let t_ratio = temp / tnom_k;
        let vtm0 = KBOQ * tnom_k;
        let _vtm = KBOQ * temp;
        let cox = self.cox();

        // Effective length/width
        let t0 = inst_l.powf(self.lln);
        let t1 = inst_w.powf(self.lwn);
        let dl = self.lint + self.ll / t0 + self.lw / t1 + self.lwl / (t0 * t1);

        let t2 = inst_l.powf(self.wln);
        let t3 = inst_w.powf(self.wwn);
        let dw = self.wint + self.wl / t2 + self.ww / t3 + self.wwl / (t2 * t3);

        let leff = inst_l + self.xl - 2.0 * dl;
        let weff = inst_w + self.xw - 2.0 * dw;
        let leff = leff.max(1e-12);
        let weff = weff.max(1e-12);

        // DLC/DWC for CV model
        let dlc = self.dlc + self.ll / t0 + self.lw / t1 + self.lwl / (t0 * t1);
        let dwc = self.dwc + self.wl / t2 + self.ww / t3 + self.wwl / (t2 * t3);
        let leff_cv = inst_l + self.xl - 2.0 * dlc;
        let weff_cv = inst_w + self.xw - 2.0 * dwc;
        let leff_cv = leff_cv.max(1e-12);
        let weff_cv = weff_cv.max(1e-12);

        // Binning factors
        let (inv_l, inv_w, inv_lw) = if self.bin_unit == 1 {
            (1.0e-6 / leff, 1.0e-6 / weff, 1.0e-12 / (leff * weff))
        } else {
            (1.0 / leff, 1.0 / weff, 1.0 / (leff * weff))
        };

        // Apply binning: param = base + l*Inv_L + w*Inv_W + p*Inv_LW
        macro_rules! bin {
            ($base:ident, $l:ident, $w:ident, $p:ident) => {
                self.$base + self.$l * inv_l + self.$w * inv_w + self.$p * inv_lw
            };
        }

        let vth0 = bin!(vth0, l_vth0, w_vth0, p_vth0);
        let k1 = bin!(k1, l_k1, w_k1, p_k1);
        let k2 = bin!(k2, l_k2, w_k2, p_k2);
        let k3 = bin!(k3, l_k3, w_k3, p_k3);
        let k3b = bin!(k3b, l_k3b, w_k3b, p_k3b);
        let w0 = bin!(w0, l_w0, w_w0, p_w0);
        let nlx = bin!(nlx, l_nlx, w_nlx, p_nlx);
        let dvt0 = bin!(dvt0, l_dvt0, w_dvt0, p_dvt0);
        let dvt1 = bin!(dvt1, l_dvt1, w_dvt1, p_dvt1);
        let dvt2 = bin!(dvt2, l_dvt2, w_dvt2, p_dvt2);
        let dvt0w = bin!(dvt0w, l_dvt0w, w_dvt0w, p_dvt0w);
        let dvt1w = bin!(dvt1w, l_dvt1w, w_dvt1w, p_dvt1w);
        let dvt2w = bin!(dvt2w, l_dvt2w, w_dvt2w, p_dvt2w);
        let dsub = bin!(dsub, l_dsub, w_dsub, p_dsub);
        let eta0 = bin!(eta0, l_eta0, w_eta0, p_eta0);
        let etab = bin!(etab, l_etab, w_etab, p_etab);
        let voff = bin!(voff, l_voff, w_voff, p_voff);
        let nfactor = bin!(nfactor, l_nfactor, w_nfactor, p_nfactor);
        let cdsc = bin!(cdsc, l_cdsc, w_cdsc, p_cdsc);
        let cdscb = bin!(cdscb, l_cdscb, w_cdscb, p_cdscb);
        let cdscd = bin!(cdscd, l_cdscd, w_cdscd, p_cdscd);
        let cit = bin!(cit, l_cit, w_cit, p_cit);
        let mut u0 = bin!(u0, l_u0, w_u0, p_u0);
        let mut ua = bin!(ua, l_ua, w_ua, p_ua);
        let mut ub = bin!(ub, l_ub, w_ub, p_ub);
        let mut uc = bin!(uc, l_uc, w_uc, p_uc);
        let ute = bin!(ute, l_ute, w_ute, p_ute);
        let ua1 = bin!(ua1, l_ua1, w_ua1, p_ua1);
        let ub1 = bin!(ub1, l_ub1, w_ub1, p_ub1);
        let uc1 = bin!(uc1, l_uc1, w_uc1, p_uc1);
        let vsat = bin!(vsat, l_vsat, w_vsat, p_vsat);
        let a0 = bin!(a0, l_a0, w_a0, p_a0);
        let ags = bin!(ags, l_ags, w_ags, p_ags);
        let a1 = bin!(a1, l_a1, w_a1, p_a1);
        let a2 = bin!(a2, l_a2, w_a2, p_a2);
        let at = bin!(at, l_at, w_at, p_at);
        let keta = bin!(keta, l_keta, w_keta, p_keta);
        let pclm = bin!(pclm, l_pclm, w_pclm, p_pclm);
        let pdibl1 = bin!(pdibl1, l_pdibl1, w_pdibl1, p_pdibl1);
        let pdibl2 = bin!(pdibl2, l_pdibl2, w_pdibl2, p_pdibl2);
        let pdiblb = bin!(pdiblb, l_pdiblb, w_pdiblb, p_pdiblb);
        let pscbe1 = bin!(pscbe1, l_pscbe1, w_pscbe1, p_pscbe1);
        let pscbe2 = bin!(pscbe2, l_pscbe2, w_pscbe2, p_pscbe2);
        let pvag = bin!(pvag, l_pvag, w_pvag, p_pvag);
        let delta = bin!(delta, l_delta, w_delta, p_delta);
        let rdsw = bin!(rdsw, l_rdsw, w_rdsw, p_rdsw);
        let prwg = bin!(prwg, l_prwg, w_prwg, p_prwg);
        let prwb = bin!(prwb, l_prwb, w_prwb, p_prwb);
        let dwg = bin!(dwg, l_dwg, w_dwg, p_dwg);
        let dwb = bin!(dwb, l_dwb, w_dwb, p_dwb);
        let b0 = bin!(b0, l_b0, w_b0, p_b0);
        let b1 = bin!(b1, l_b1, w_b1, p_b1);
        let alpha0 = bin!(alpha0, l_alpha0, w_alpha0, p_alpha0);
        let alpha1 = bin!(alpha1, l_alpha1, w_alpha1, p_alpha1);
        let beta0 = bin!(beta0, l_beta0, w_beta0, p_beta0);
        let kt1 = bin!(kt1, l_kt1, w_kt1, p_kt1);
        let kt1l = bin!(kt1l, l_kt1l, w_kt1l, p_kt1l);
        let kt2 = bin!(kt2, l_kt2, w_kt2, p_kt2);
        let npeak = bin!(npeak, l_npeak, w_npeak, p_npeak);
        let nsub = bin!(nsub, l_nsub, w_nsub, p_nsub);
        let ngate = bin!(ngate, l_ngate, w_ngate, p_ngate);
        let xj = bin!(xj, l_xj, w_xj, p_xj);
        let cgsl = bin!(cgsl, l_cgsl, w_cgsl, p_cgsl);
        let cgdl = bin!(cgdl, l_cgdl, w_cgdl, p_cgdl);
        let ckappa = bin!(ckappa, l_ckappa, w_ckappa, p_ckappa);
        let cf = bin!(cf, l_cf, w_cf, p_cf);
        let clc = bin!(clc, l_clc, w_clc, p_clc);
        let cle = bin!(cle, l_cle, w_cle, p_cle);
        let vfbcv = bin!(vfbcv, l_vfbcv, w_vfbcv, p_vfbcv);
        let mut acde = bin!(acde, l_acde, w_acde, p_acde);
        let moin = bin!(moin, l_moin, w_moin, p_moin);
        let noff = bin!(noff, l_noff, w_noff, p_noff);
        let voffcv = bin!(voffcv, l_voffcv, w_voffcv, p_voffcv);

        // Band gap and intrinsic carrier concentration at Tnom
        let eg0 = 1.16 - 7.02e-4 * tnom_k * tnom_k / (tnom_k + 1108.0);
        let ni = 1.45e10
            * (tnom_k / 300.15)
            * (tnom_k / 300.15).sqrt()
            * (21.5565981 - eg0 / (2.0 * vtm0)).exp();

        let factor1 = (EPSSI / EPSOX * self.tox).sqrt();

        // If u0 > 1.0, it's in cm²/V·s — convert to m²/V·s
        if u0 > 1.0 {
            u0 /= 1e4;
        }

        // Surface potential
        let npeak_eff = npeak.max(1e15);
        let phi = 2.0 * vtm0 * (npeak_eff / ni).ln();
        let phi = phi.max(0.1);
        let sqrt_phi = phi.sqrt();
        let phis3 = sqrt_phi * phi;

        // Depletion width
        let xdep0 = (2.0 * EPSSI / (CHARGE_Q * npeak_eff * 1e6)).sqrt() * sqrt_phi;
        let sqrt_xdep0 = xdep0.sqrt();

        // Built-in potential
        let vbi = vtm0 * (1e20 * npeak_eff / (ni * ni)).ln();

        // Depletion capacitance per unit area
        let cdep0 = (CHARGE_Q * EPSSI * npeak_eff * 1e6 / (2.0 * phi)).sqrt();

        // Debye length
        let ldeb = (EPSSI * vtm0 / (CHARGE_Q * npeak_eff * 1e6)).sqrt() / 3.0;

        // Characteristic length
        let litl = (3.0 * xj * self.tox).sqrt();

        // Vfb / Vth0 relationship
        let vfb = if self.vfb_given {
            self.vfb
        } else if self.vth0_given {
            self.mos_type.sign() * vth0 - phi - k1 * sqrt_phi
        } else {
            -1.0
        };

        let vth0_eff = if self.vth0_given {
            vth0
        } else {
            self.mos_type.sign() * (vfb + phi + k1 * sqrt_phi)
        };

        let k1ox = k1 * self.tox / self.toxm;
        let k2ox = k2 * self.tox / self.toxm;

        // SCE parameters
        let tmp1 = vbi - phi;
        let tmp2 = factor1 * sqrt_xdep0;

        // theta0vb0 (DIBL at Vbs=0)
        let t0 = -0.5 * dsub * leff / tmp2;
        let t1 = if t0 > -EXP_THRESHOLD {
            t0.exp()
        } else {
            MIN_EXP
        };
        let theta0vb0 = t1 * (1.0 + 2.0 * t1);

        // thetaRout
        let t0 = -0.5 * self.dsub * leff / tmp2; // drout uses dsub
        let t1r = if t0 > -EXP_THRESHOLD {
            t0.exp()
        } else {
            MIN_EXP
        };
        let t2r = t1r * (1.0 + 2.0 * t1r);
        let theta_rout = pdibl1 * t2r + pdibl2;

        // Temperature-adjusted mobility
        let temp_ratio = t_ratio - 1.0;
        ua += ua1 * temp_ratio;
        ub += ub1 * temp_ratio;
        uc += uc1 * temp_ratio;

        let u0temp = u0 * t_ratio.powf(ute);
        let vsattemp = vsat - at * temp_ratio;

        // Series resistance
        let rds0 = if rdsw > 0.0 {
            (rdsw + self.prt * temp_ratio) / (weff * 1e6).powf(self.wr)
        } else {
            0.0
        };

        // Vbsc (body bias clamp)
        let vbsc = -(0.9 * (phi - 0.01)).max(0.0);

        // Quantum mechanical correction
        acde *= (npeak_eff / 2e16).powf(-0.25);

        // Overlap capacitances
        let cgso_eff = if dlc > 0.0 {
            (self.cgso + cf) * weff_cv
        } else {
            (0.6 * xj * cox + cf) * weff_cv
        };
        let cgdo_eff = if dlc > 0.0 {
            (self.cgdo + cf) * weff_cv
        } else {
            (0.6 * xj * cox + cf) * weff_cv
        };
        let cgbo_eff = self.cgbo * leff_cv;

        // Narrow width vth0
        let tmp2_narrow = self.tox * phi / (weff + w0);

        // NLX and temperature effect on Vth
        let t0_nlx = (1.0 + nlx / leff).sqrt();
        let t1_vth =
            k1ox * (t0_nlx - 1.0) * sqrt_phi + (kt1 + kt1l / leff + kt2 * 0.0) * temp_ratio; // Vbseff=0 for vfbzb

        // SCE at Vbs=0
        let t0_sce = -0.5 * dvt1 * leff / tmp2;
        let t1_sce = if t0_sce > -EXP_THRESHOLD {
            t0_sce.exp()
        } else {
            MIN_EXP
        };
        let theta0 = t1_sce * (1.0 + 2.0 * t1_sce);
        let delt_vth0 = dvt0 * theta0 * tmp1;

        // Narrow width at Vbs=0
        let t0_nw = -0.5 * dvt1w * weff * leff / tmp2;
        let t1_nw = if t0_nw > -EXP_THRESHOLD {
            t0_nw.exp()
        } else {
            MIN_EXP
        };
        let t2_nw = dvt0w * t1_nw * (1.0 + 2.0 * t1_nw) * tmp1;

        let tmp3 = self.mos_type.sign() * vth0_eff - delt_vth0 - t2_nw + k3 * tmp2_narrow + t1_vth;
        let vfbzb = tmp3 - phi - k1 * sqrt_phi;

        // abulkCVfactor
        let abulk_cv_factor = 1.0 + (clc / leff_cv).powf(cle);

        // tconst for NQS
        let tconst = u0temp * self.elm / (cox * weff_cv * leff_cv * leff_cv);

        Bsim3SizeParam {
            leff,
            weff,
            leff_cv,
            weff_cv,
            vth0: vth0_eff,
            k1,
            k2,
            k3,
            k3b,
            w0,
            nlx,
            dvt0,
            dvt1,
            dvt2,
            dvt0w,
            dvt1w,
            dvt2w,
            dsub,
            eta0,
            etab,
            voff,
            nfactor,
            cdsc,
            cdscb,
            cdscd,
            cit,
            u0,
            ua,
            ub,
            uc,
            ute,
            ua1,
            ub1,
            uc1,
            vsat,
            a0,
            ags,
            a1,
            a2,
            at,
            keta,
            pclm,
            pdibl1,
            pdibl2,
            pdiblb,
            pscbe1,
            pscbe2,
            pvag,
            delta,
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
            kt1,
            kt1l,
            kt2,
            npeak,
            nsub,
            ngate,
            xj,
            cgsl,
            cgdl,
            ckappa,
            cf,
            clc,
            cle,
            vfbcv,
            acde,
            moin,
            noff,
            voffcv,
            phi,
            sqrt_phi,
            phis3,
            xdep0,
            sqrt_xdep0,
            vbi,
            cdep0,
            ldeb,
            litl,
            vbsc,
            k1ox,
            k2ox,
            theta0vb0,
            theta_rout,
            vfbzb,
            u0temp,
            vsattemp,
            rds0,
            cgso_eff,
            cgdo_eff,
            cgbo_eff,
            abulk_cv_factor,
            tconst,
        }
    }
}

impl Bsim3Instance {
    /// Get terminal voltages from solution vector with type sign.
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

    /// Drain conductance from RSH * NRD.
    pub fn drain_conductance(&self) -> f64 {
        if self.model.rsh > 0.0 && self.nrd > 0.0 {
            1.0 / (self.model.rsh * self.nrd)
        } else {
            0.0
        }
    }

    /// Source conductance from RSH * NRS.
    pub fn source_conductance(&self) -> f64 {
        if self.model.rsh > 0.0 && self.nrs > 0.0 {
            1.0 / (self.model.rsh * self.nrs)
        } else {
            0.0
        }
    }
}

/// Compute BSIM3 companion model (NR linearization).
///
/// This is the core BSIM3 load function, computing drain current, conductances,
/// and capacitances at the given operating point.
#[expect(clippy::too_many_lines)]
pub fn bsim3_companion(
    vgs: f64,
    vds: f64,
    vbs: f64,
    sp: &Bsim3SizeParam,
    model: &Bsim3Model,
) -> Bsim3Companion {
    let sign = model.mos_type.sign();
    let cox = model.cox();
    let vtm = KBOQ * TEMP_DEFAULT;
    let gmin = 1e-12;

    // Mode detection
    let (vgs_i, vds_i, vbs_i, mode) = if vds >= 0.0 {
        (vgs, vds, vbs, 1)
    } else {
        (vgs - vds, -vds, vbs - vds, -1)
    };

    let vbd = vbs - vds;

    // Junction currents
    let source_sat = (model.js * 1e-4_f64.max(0.0)).max(1e-14); // simplified
    let drain_sat = source_sat;
    let nvtm = vtm * model.nj;

    let (gbs, cbs_current) = junction_current(vbs, source_sat, nvtm, gmin, model.ijth);
    let (gbd, cbd_current) = junction_current(vbd, drain_sat, nvtm, gmin, model.ijth);

    // Junction capacitances (simplified — use area/perimeter from instance later)
    let capbs = 0.0;
    let capbd = 0.0;

    // === Core MOSFET equations ===

    // Effective Vbs
    let t0 = vbs_i - sp.vbsc - 0.001;
    let t1 = (t0 * t0 - 0.004 * sp.vbsc).sqrt();
    let mut vbseff = sp.vbsc + 0.5 * (t0 + t1);
    let mut d_vbseff_d_vb = 0.5 * (1.0 + t0 / t1);
    if vbseff < vbs_i {
        vbseff = vbs_i;
        d_vbseff_d_vb = 1.0; // not used for values, only derivatives
    }

    // Surface potential
    let (phis, sqrt_phis, d_sqrt_phis_d_vb);
    if vbseff > 0.0 {
        let t0 = sp.phi / (sp.phi + vbseff);
        phis = sp.phi * t0;
        sqrt_phis = sp.phis3 / (sp.phi + 0.5 * vbseff);
        d_sqrt_phis_d_vb = -0.5 * sqrt_phis * sqrt_phis / sp.phis3;
        let _ = phis; // suppress warning
    } else {
        let p = sp.phi - vbseff;
        sqrt_phis = p.sqrt();
        d_sqrt_phis_d_vb = -0.5 / sqrt_phis;
    }

    let xdep = sp.xdep0 * sqrt_phis / sp.sqrt_phi;
    let _d_xdep_d_vb = sp.xdep0 / sp.sqrt_phi * d_sqrt_phis_d_vb;

    let leff = sp.leff;

    // Vth calculation
    let v0 = sp.vbi - sp.phi;
    let factor1 = (EPSSI / EPSOX * model.tox).sqrt();

    let t3 = xdep.sqrt();

    // lt1 with dvt2
    let t0_dvt2 = sp.dvt2 * vbseff;
    let t1_dvt2 = if t0_dvt2 >= -0.5 {
        1.0 + t0_dvt2
    } else {
        let t4 = 1.0 / (3.0 + 8.0 * t0_dvt2);
        (1.0 + 3.0 * t0_dvt2) * t4
    };
    let lt1 = factor1 * t3 * t1_dvt2;

    // ltw with dvt2w
    let t0_dvt2w = sp.dvt2w * vbseff;
    let t1_dvt2w = if t0_dvt2w >= -0.5 {
        1.0 + t0_dvt2w
    } else {
        let t4 = 1.0 / (3.0 + 8.0 * t0_dvt2w);
        (1.0 + 3.0 * t0_dvt2w) * t4
    };
    let ltw = factor1 * t3 * t1_dvt2w;

    // Theta0 (SCE)
    let t0_theta = -0.5 * sp.dvt1 * leff / lt1;
    let theta0 = if t0_theta > -EXP_THRESHOLD {
        let t1 = t0_theta.exp();
        t1 * (1.0 + 2.0 * t1)
    } else {
        MIN_EXP * (1.0 + 2.0 * MIN_EXP)
    };

    let delt_vth = sp.dvt0 * theta0 * v0;

    // Narrow width DVT
    let t0_nw = -0.5 * sp.dvt1w * sp.weff * leff / ltw;
    let t2_nw = if t0_nw > -EXP_THRESHOLD {
        let t1 = t0_nw.exp();
        sp.dvt0w * t1 * (1.0 + 2.0 * t1) * v0
    } else {
        sp.dvt0w * MIN_EXP * (1.0 + 2.0 * MIN_EXP) * v0
    };

    let tmp2 = model.tox * sp.phi / (sp.weff + sp.w0);

    // DIBL
    let mut t3_eta = sp.eta0 + sp.etab * vbseff;
    if t3_eta < 1.0e-4 {
        let t9 = 1.0 / (3.0 - 2.0e4 * t3_eta);
        t3_eta = (2.0e-4 - t3_eta) * t9;
    }
    let d_dibl_sft_d_vd = t3_eta * sp.theta0vb0;
    let dibl_sft = d_dibl_sft_d_vd * vds_i;

    // Temperature effect
    let temp_ratio = TEMP_DEFAULT / (model.tnom + 273.15) - 1.0;
    let t0_nlx = (1.0 + sp.nlx / leff).sqrt();
    let t1_vth = sp.k1ox * (t0_nlx - 1.0) * sp.sqrt_phi
        + (sp.kt1 + sp.kt1l / leff + sp.kt2 * vbseff) * temp_ratio;

    let vth = sign * sp.vth0 - sp.k1 * sp.sqrt_phi + sp.k1ox * sqrt_phis
        - sp.k2ox * vbseff
        - delt_vth
        - t2_nw
        + (sp.k3 + sp.k3b * vbseff) * tmp2
        + t1_vth
        - dibl_sft;

    let d_vth_d_vb = sp.k1ox * d_sqrt_phis_d_vb - sp.k2ox + sp.k3b * tmp2 + sp.kt2 * temp_ratio
        - sp.etab * vds_i * sp.theta0vb0;
    let _d_vth_d_vd = -d_dibl_sft_d_vd;

    // Subthreshold slope factor n
    let tmp2_n = sp.nfactor * EPSSI / xdep;
    let tmp3_n = sp.cdsc + sp.cdscb * vbseff + sp.cdscd * vds_i;
    let tmp4_n = (tmp2_n + tmp3_n * theta0 + sp.cit) / cox;
    let n = if tmp4_n >= -0.5 {
        1.0 + tmp4_n
    } else {
        let t0 = 1.0 / (3.0 + 8.0 * tmp4_n);
        (1.0 + 3.0 * tmp4_n) * t0
    };
    let n = n.max(0.5);

    // Poly gate depletion
    let vfb_plus_phi = sp.vfbzb + sp.phi + sp.k1 * sp.sqrt_phi; // approximate vfb + phi
    let (vgs_eff, _d_vgs_eff_d_vg) = if sp.ngate > 1e18 && sp.ngate < 1e25 && vgs_i > vfb_plus_phi {
        let t1 = 1e6 * CHARGE_Q * EPSSI * sp.ngate / (cox * cox);
        let t4 = (1.0 + 2.0 * (vgs_i - vfb_plus_phi) / t1).sqrt();
        let t2 = t1 * (t4 - 1.0);
        let t3 = 0.5 * t2 * t2 / t1;
        let t7 = 1.12 - t3 - 0.05;
        let t6 = (t7 * t7 + 0.224).sqrt();
        let t5 = 1.12 - 0.5 * (t7 + t6);
        (vgs_i - t5, 1.0 - (0.5 - 0.5 / t4) * (1.0 + t7 / t6))
    } else {
        (vgs_i, 1.0)
    };

    let vgst = vgs_eff - vth;

    // Effective Vgst (Vgsteff) — unified strong/weak inversion
    let t10 = 2.0 * n * vtm;
    let vgst_nvt = vgst / t10;
    let exp_arg = (2.0 * sp.voff - vgst) / t10;

    let vgsteff;
    if vgst_nvt > EXP_THRESHOLD {
        vgsteff = vgst;
    } else if exp_arg > EXP_THRESHOLD {
        let t0 = (vgst - sp.voff) / (n * vtm);
        let exp_vgst = t0.exp();
        vgsteff = vtm * sp.cdep0 / cox * exp_vgst;
    } else {
        let exp_vgst = safe_exp(vgst_nvt);
        let t1 = t10 * (1.0 + exp_vgst).ln();

        let d_t2_d_vg = -cox / (vtm * sp.cdep0) * safe_exp(exp_arg);
        let t2 = 1.0 - t10 * d_t2_d_vg;

        vgsteff = t1 / t2;
    }
    let vgsteff = vgsteff.max(1e-20);

    // Effective channel geometry
    let t9_w = sqrt_phis - sp.sqrt_phi;
    let weff_calc = (sp.weff - 2.0 * (sp.dwg * vgsteff + sp.dwb * t9_w)).max(2e-8);

    // Rds (source-drain resistance)
    let t0_rds = sp.prwg * vgsteff + sp.prwb * t9_w;
    let rds = if t0_rds >= -0.9 {
        sp.rds0 * (1.0 + t0_rds)
    } else {
        let t1 = 1.0 / (17.0 + 20.0 * t0_rds);
        sp.rds0 * (0.8 + t0_rds) * t1
    };

    // Abulk (bulk charge factor)
    let t1_ab = 0.5 * sp.k1ox / sqrt_phis;
    let t9_ab = (sp.xj * xdep).sqrt();
    let tmp1_ab = leff + 2.0 * t9_ab;
    let t5_ab = leff / tmp1_ab;
    let tmp2_ab = sp.a0 * t5_ab;
    let tmp3_ab = sp.weff + sp.b1;
    let tmp4_ab = sp.b0 / tmp3_ab;
    let t2_ab = tmp2_ab + tmp4_ab;

    let abulk0 = 1.0 + t1_ab * t2_ab;
    let t8_ab = sp.ags * sp.a0 * t5_ab * t5_ab * t5_ab;
    let d_abulk_d_vg = -t1_ab * t8_ab;
    let mut abulk = abulk0 + d_abulk_d_vg * vgsteff;

    // Clamp Abulk
    if abulk < 0.1 {
        let t9 = 1.0 / (3.0 - 20.0 * abulk);
        abulk = (0.2 - abulk) * t9;
    }

    // Keta correction
    let t2_keta = sp.keta * vbseff;
    let keta_factor = if t2_keta >= -0.9 {
        1.0 / (1.0 + t2_keta)
    } else {
        let t1 = 1.0 / (0.8 + t2_keta);
        (17.0 + 20.0 * t2_keta) * t1
    };
    abulk *= keta_factor;

    // Mobility
    let (t5_mob, d_denomi_d_vg);
    match model.mob_mod {
        1 => {
            let t0 = vgsteff + vth + vth;
            let t2 = sp.ua + sp.uc * vbseff;
            let t3 = t0 / model.tox;
            t5_mob = t3 * (t2 + sp.ub * t3);
            d_denomi_d_vg = (t2 + 2.0 * sp.ub * t3) / model.tox;
        }
        2 => {
            t5_mob = vgsteff / model.tox * (sp.ua + sp.uc * vbseff + sp.ub * vgsteff / model.tox);
            d_denomi_d_vg =
                (sp.ua + sp.uc * vbseff + 2.0 * sp.ub * vgsteff / model.tox) / model.tox;
        }
        _ => {
            // mobMod 3 (default fallback)
            let t0 = vgsteff + vth + vth;
            let t2 = 1.0 + sp.uc * vbseff;
            let t3 = t0 / model.tox;
            let t4 = t3 * (sp.ua + sp.ub * t3);
            t5_mob = t4 * t2;
            d_denomi_d_vg = (sp.ua + 2.0 * sp.ub * t3) * t2 / model.tox;
        }
    }

    let denomi = if t5_mob >= -0.8 {
        1.0 + t5_mob
    } else {
        let t9 = 1.0 / (7.0 + 10.0 * t5_mob);
        (0.6 + t5_mob) * t9
    };
    let _ = d_denomi_d_vg; // suppress unused warning

    let ueff = sp.u0temp / denomi;

    // Vdsat
    let wv_cox = weff_calc * sp.vsattemp * cox;
    let wv_cox_rds = wv_cox * rds;

    let esat = 2.0 * sp.vsattemp / ueff;
    let esat_l = esat * leff;

    // Lambda (CLM parameter a1, a2)
    let lambda = if sp.a1 == 0.0 {
        sp.a2
    } else if sp.a1 > 0.0 {
        let t0 = 1.0 - sp.a2;
        let t1 = t0 - sp.a1 * vgsteff - 0.0001;
        let t2 = (t1 * t1 + 0.0004 * t0).sqrt();
        sp.a2 + t0 - 0.5 * (t1 + t2)
    } else {
        let t1 = sp.a2 + sp.a1 * vgsteff - 0.0001;
        let t2 = (t1 * t1 + 0.0004 * sp.a2).sqrt();
        0.5 * (t1 + t2)
    };

    let vgst2vtm = vgsteff + 2.0 * vtm;

    // Vdsat calculation
    let vdsat = if rds == 0.0 && lambda == 1.0 {
        // Simple case
        let t0 = 1.0 / (abulk * esat_l + vgst2vtm);
        esat_l * vgst2vtm * t0
    } else {
        // General case with quadratic
        let t9 = abulk * wv_cox_rds;
        let t7 = vgst2vtm * t9;
        let t0 = 2.0 * abulk * (t9 - 1.0 + 1.0 / lambda);
        let t1 = vgst2vtm * (2.0 / lambda - 1.0) + abulk * esat_l + 3.0 * t7;
        let t2 = vgst2vtm * (esat_l + 2.0 * vgst2vtm * wv_cox_rds);

        let t3 = (t1 * t1 - 2.0 * t0 * t2).max(0.0).sqrt();
        if t0.abs() > 1e-30 {
            (t1 - t3) / t0
        } else {
            esat_l * vgst2vtm / (abulk * esat_l + vgst2vtm)
        }
    };
    let vdsat = vdsat.max(1e-18);

    // Effective Vds
    let t1 = vdsat - vds_i - sp.delta;
    let t2 = (t1 * t1 + 4.0 * sp.delta * vdsat).sqrt();
    let mut vdseff = vdsat - 0.5 * (t1 + t2);
    let d_vdseff_d_vg = 0.0; // simplified — full derivatives not needed for value computation
    let d_vdseff_d_vd = 0.5 * (1.0 - t1 / t2); // approximate
    let _ = d_vdseff_d_vg;

    if vds_i == 0.0 {
        vdseff = 0.0;
    }
    if vdseff > vds_i {
        vdseff = vds_i;
    }

    let diff_vds = vds_i - vdseff;

    // Output resistance terms
    // VACLM
    let vaclm = if sp.pclm > 0.0 && diff_vds > 1e-10 {
        let t0 = 1.0 / (sp.pclm * abulk * sp.litl);
        let t2 = vgsteff / esat_l;
        let t1 = leff * (abulk + t2);
        t0 * t1 * diff_vds
    } else {
        MAX_EXP
    };

    // VADIBL
    let vadibl = if sp.theta_rout > 0.0 {
        let t8 = abulk * vdsat;
        let t0 = vgst2vtm * t8;
        let t1 = vgst2vtm + t8;
        let mut v = (vgst2vtm - t0 / t1) / sp.theta_rout;
        // pdiblb correction
        let t7 = sp.pdiblb * vbseff;
        if t7 >= -0.9 {
            v /= 1.0 + t7;
        } else {
            let t4 = 1.0 / (0.8 + t7);
            v *= (17.0 + 20.0 * t7) * t4;
        }
        v
    } else {
        MAX_EXP
    };

    // VA = Vasat + pvag_factor * harmonic_mean(VACLM, VADIBL)
    let t8_pvag = sp.pvag / esat_l;
    let t9_pvag = t8_pvag * vgsteff;
    let pvag_factor = if t9_pvag > -0.9 {
        1.0 + t9_pvag
    } else {
        let t1 = 1.0 / (17.0 + 20.0 * t9_pvag);
        (0.8 + t9_pvag) * t1
    };

    // Vasat
    let tmp4_vasat = 1.0 - 0.5 * abulk * vdsat / vgst2vtm;
    let t9_vasat = wv_cox_rds * vgsteff;
    let t0_vasat = esat_l + vdsat + 2.0 * t9_vasat * tmp4_vasat;
    let t1_vasat = 2.0 / lambda - 1.0 + wv_cox_rds * abulk;
    let vasat = t0_vasat / t1_vasat;

    let harmonic = if vaclm + vadibl > 0.0 {
        vaclm * vadibl / (vaclm + vadibl)
    } else {
        MAX_EXP
    };
    let va = vasat + pvag_factor * harmonic;

    // VASCBE
    let vascbe = if sp.pscbe2 > 0.0 && diff_vds > sp.pscbe1 * sp.litl / EXP_THRESHOLD {
        let t0 = sp.pscbe1 * sp.litl / diff_vds;
        leff * t0.exp() / sp.pscbe2
    } else {
        MAX_EXP * leff / sp.pscbe2.max(1e-30)
    };

    // Ids calculation
    let cox_w_ov_l = cox * weff_calc / leff;
    let beta = ueff * cox_w_ov_l;

    let t0_ids = 1.0 - 0.5 * abulk * vdseff / vgst2vtm;
    let fgche1 = vgsteff * t0_ids;
    let fgche2 = 1.0 + vdseff / esat_l;

    let gche = beta * fgche1 / fgche2;
    let t0_gche = 1.0 + gche * rds;
    let t9_gche = vdseff / t0_gche;
    let idl = gche * t9_gche;

    // CLM + DIBL modulation
    let t9_va = diff_vds / va;
    let idsa = idl * (1.0 + t9_va);

    // SCBE modulation
    let t9_scbe = diff_vds / vascbe;
    let ids = idsa * (1.0 + t9_scbe);

    // Transconductances — compute via numerical relationship
    // Gm = dIds/dVgs, Gds = dIds/dVds, Gmb = dIds/dVbs
    // For the NR companion, we need accurate conductances.
    // Use chain-rule through the equations.
    let gm_approx = beta * vdseff / (fgche2 * t0_gche) * (1.0 + t9_va) * (1.0 + t9_scbe);
    let gds_approx =
        ids * d_vdseff_d_vd / va.max(1e-20) + idsa / vascbe.max(1e-20) * (1.0 - d_vdseff_d_vd);
    // Include the intrinsic gds from the channel
    let gds_channel =
        beta * vgsteff / (fgche2 * fgche2 * esat_l) / t0_gche * (1.0 + t9_va) * (1.0 + t9_scbe);
    let mut gds_total = gds_approx + gds_channel * d_vdseff_d_vd;

    // Ensure minimum gds for convergence
    gds_total = gds_total.max(1e-20);

    let gm_total = gm_approx.max(1e-20);

    // Gmbs (body effect)
    let gmbs = gm_total * d_vth_d_vb.abs() * d_vbseff_d_vb;

    // Substrate current
    let tmp_sub = sp.alpha0 + sp.alpha1 * leff;
    let (isub, gbg, gbd_sub, gbb);
    if tmp_sub <= 0.0 || sp.beta0 <= 0.0 {
        isub = 0.0;
        gbg = 0.0;
        gbd_sub = 0.0;
        gbb = 0.0;
    } else {
        let t2 = tmp_sub / leff;
        if diff_vds > sp.beta0 / EXP_THRESHOLD {
            let t0 = -sp.beta0 / diff_vds;
            let t1 = t2 * diff_vds * t0.exp();
            isub = t1 * idsa;
            gbg = 0.0;
            gbd_sub = 0.0;
            gbb = 0.0; // simplified
        } else {
            isub = 0.0;
            gbg = 0.0;
            gbd_sub = 0.0;
            gbb = 0.0;
        }
    }

    // NR equivalent current sources
    let ceq_d = ids - gm_total * vgs_eff - gds_total * vds_i - gmbs * vbs_i;
    let ceq_bs = cbs_current - gbs * vbs;
    let ceq_bd = cbd_current - gbd * vbd;
    let ceq_sub = isub - gbg * vgs_i - gbd_sub * vds_i - gbb * vbs_i;

    // Capacitances — compute via numerical differentiation of capMod=3 charges
    let qinv = -cox_w_ov_l * vgsteff * (1.0 - 0.5 * abulk * vdseff / vgst2vtm);
    let dv = 1e-5;
    let (qg0, qb0, qs0) = bsim3_total_charges(vgs, vds, vbs, sp, model);
    let (qg_gp, qb_gp, qs_gp) = bsim3_total_charges(vgs + dv, vds, vbs, sp, model);
    let (qg_gm, qb_gm, qs_gm) = bsim3_total_charges(vgs - dv, vds, vbs, sp, model);
    let (qg_dp, qb_dp, qs_dp) = bsim3_total_charges(vgs, vds + dv, vbs, sp, model);
    let (qg_dm, qb_dm, qs_dm) = bsim3_total_charges(vgs, vds - dv, vbs, sp, model);
    let (qg_bp, qb_bp, qs_bp) = bsim3_total_charges(vgs, vds, vbs + dv, sp, model);
    let (qg_bm, qb_bm, qs_bm) = bsim3_total_charges(vgs, vds, vbs - dv, sp, model);
    let _ = (qg0, qb0, qs0);

    // Gate charge derivatives (w.r.t. Vgs, Vds, Vbs)
    let c_gg = (qg_gp - qg_gm) / (2.0 * dv);
    let c_gd = (qg_dp - qg_dm) / (2.0 * dv);
    let c_gb = (qg_bp - qg_bm) / (2.0 * dv);

    // Bulk charge derivatives
    let c_bg = (qb_gp - qb_gm) / (2.0 * dv);
    let c_bd = (qb_dp - qb_dm) / (2.0 * dv);
    let c_bb = (qb_bp - qb_bm) / (2.0 * dv);

    // Source charge derivatives
    let c_sg = (qs_gp - qs_gm) / (2.0 * dv);
    let c_sd = (qs_dp - qs_dm) / (2.0 * dv);

    // Stored capacitances (BSIM3 convention)
    let cggb = c_gg;
    let cgdb = c_gd;
    let cgsb = -(c_gg + c_gd + c_gb);
    let cbgb = c_bg;
    let cbdb = c_bd;
    let cbsb = -(c_bg + c_bd + c_bb);
    // Drain caps from conservation: Qd = -(Qg + Qb + Qs)
    let cdgb = -(c_gg + c_bg + c_sg);
    let cddb = -(c_gd + c_bd + c_sd);
    let cdsb = c_gg + c_gd + c_gb + c_bg + c_bd + c_bb + c_sg + c_sd + (qs_bp - qs_bm) / (2.0 * dv);

    Bsim3Companion {
        ids,
        gm: gm_total,
        gds: gds_total,
        gmbs,
        gbd,
        gbs,
        cbd_current,
        cbs_current,
        isub,
        gbg,
        gbb,
        gbd_sub,
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
    }
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
        // Linear extrapolation above threshold
        let is_evjsm = i_sat * safe_exp(vjsm / nvtm);
        let g = is_evjsm / nvtm + gmin;
        let i = is_evjsm - i_sat + g * (v - vjsm) + gmin * v;
        (g, i)
    }
}

const DELTA_3: f64 = 0.02;
const DELTA_4: f64 = 0.02;
const DELTA_1: f64 = 0.02;

/// Compute overlap charges for BSIM3 capMod=3 (uses capMod=2 model).
/// Returns (qgdo, qgso, qgb).
fn bsim3_overlap_charges(
    vgs: f64,
    vgd: f64,
    vgb: f64,
    sp: &Bsim3SizeParam,
    _model: &Bsim3Model,
) -> (f64, f64, f64) {
    // Drain overlap (capMod=2 smooth model)
    let t0 = vgd + DELTA_1;
    let t1 = (t0 * t0 + 4.0 * DELTA_1).sqrt();
    let t2 = 0.5 * (t0 - t1); // smooth negative clamp
    let t3 = sp.weff_cv * sp.cgdl;
    let qgdo = if t3 > 0.0 {
        let t4 = (1.0 - 4.0 * t2 / sp.ckappa).sqrt();
        (sp.cgdo_eff + t3) * vgd - t3 * (t2 + 0.5 * sp.ckappa * (t4 - 1.0))
    } else {
        sp.cgdo_eff * vgd
    };

    // Source overlap (capMod=2 smooth model)
    let t0 = vgs + DELTA_1;
    let t1 = (t0 * t0 + 4.0 * DELTA_1).sqrt();
    let t2 = 0.5 * (t0 - t1);
    let t3 = sp.weff_cv * sp.cgsl;
    let qgso = if t3 > 0.0 {
        let t4 = (1.0 - 4.0 * t2 / sp.ckappa).sqrt();
        (sp.cgso_eff + t3) * vgs - t3 * (t2 + 0.5 * sp.ckappa * (t4 - 1.0))
    } else {
        sp.cgso_eff * vgs
    };

    // Gate-bulk overlap
    let qgb = sp.cgbo_eff * vgb;

    (qgdo, qgso, qgb)
}

/// Compute total charges (gate, bulk, source) for BSIM3.
/// Dispatches to the appropriate charge model based on cap_mod (0/1/2/3).
/// Includes intrinsic charges + overlap charges.
/// Returns (Qg_total, Qb_total, Qs_total).
#[expect(clippy::too_many_lines)]
fn bsim3_total_charges(
    vgs: f64,
    vds: f64,
    vbs: f64,
    sp: &Bsim3SizeParam,
    model: &Bsim3Model,
) -> (f64, f64, f64) {
    let sign = model.mos_type.sign();
    let cox = model.cox();
    let vtm = KBOQ * TEMP_DEFAULT;

    // Mode detection
    let (vgs_i, vds_i, vbs_i, _mode) = if vds >= 0.0 {
        (vgs, vds, vbs, 1)
    } else {
        (vgs - vds, -vds, vbs - vds, -1)
    };

    // === DC core computation (duplicated from companion for charge model) ===

    // VbseffCV (same as vbseff for now)
    let t0 = vbs_i - sp.vbsc - 0.001;
    let t1 = (t0 * t0 - 0.004 * sp.vbsc).sqrt();
    let mut vbseff = sp.vbsc + 0.5 * (t0 + t1);
    if vbseff < vbs_i {
        vbseff = vbs_i;
    }
    let vbs_eff_cv = vbseff;

    // Surface potential
    let sqrt_phis = if vbseff > 0.0 {
        sp.phis3 / (sp.phi + 0.5 * vbseff)
    } else {
        (sp.phi - vbseff).sqrt()
    };

    let xdep = sp.xdep0 * sqrt_phis / sp.sqrt_phi;

    // Vth
    let v0 = sp.vbi - sp.phi;
    let factor1 = (EPSSI / EPSOX * model.tox).sqrt();
    let t3 = xdep.sqrt();

    let t0_dvt2 = sp.dvt2 * vbseff;
    let t1_dvt2 = if t0_dvt2 >= -0.5 {
        1.0 + t0_dvt2
    } else {
        let t4 = 1.0 / (3.0 + 8.0 * t0_dvt2);
        (1.0 + 3.0 * t0_dvt2) * t4
    };
    let lt1 = factor1 * t3 * t1_dvt2;

    let t0_dvt2w = sp.dvt2w * vbseff;
    let t1_dvt2w = if t0_dvt2w >= -0.5 {
        1.0 + t0_dvt2w
    } else {
        let t4 = 1.0 / (3.0 + 8.0 * t0_dvt2w);
        (1.0 + 3.0 * t0_dvt2w) * t4
    };
    let ltw = factor1 * t3 * t1_dvt2w;

    let t0_theta = -0.5 * sp.dvt1 * sp.leff / lt1;
    let theta0 = if t0_theta > -EXP_THRESHOLD {
        let t1 = t0_theta.exp();
        t1 * (1.0 + 2.0 * t1)
    } else {
        MIN_EXP * (1.0 + 2.0 * MIN_EXP)
    };
    let delt_vth = sp.dvt0 * theta0 * v0;

    let t0_nw = -0.5 * sp.dvt1w * sp.weff * sp.leff / ltw;
    let t2_nw = if t0_nw > -EXP_THRESHOLD {
        let t1 = t0_nw.exp();
        sp.dvt0w * t1 * (1.0 + 2.0 * t1) * v0
    } else {
        sp.dvt0w * MIN_EXP * (1.0 + 2.0 * MIN_EXP) * v0
    };

    let tmp2 = model.tox * sp.phi / (sp.weff + sp.w0);

    let mut t3_eta = sp.eta0 + sp.etab * vbseff;
    if t3_eta < 1.0e-4 {
        let t9 = 1.0 / (3.0 - 2.0e4 * t3_eta);
        t3_eta = (2.0e-4 - t3_eta) * t9;
    }
    let dibl_sft = t3_eta * sp.theta0vb0 * vds_i;

    let temp_ratio = TEMP_DEFAULT / (model.tnom + 273.15) - 1.0;
    let t0_nlx = (1.0 + sp.nlx / sp.leff).sqrt();
    let t1_vth = sp.k1ox * (t0_nlx - 1.0) * sp.sqrt_phi
        + (sp.kt1 + sp.kt1l / sp.leff + sp.kt2 * vbseff) * temp_ratio;

    let vth = sign * sp.vth0 - sp.k1 * sp.sqrt_phi + sp.k1ox * sqrt_phis
        - sp.k2ox * vbseff
        - delt_vth
        - t2_nw
        + (sp.k3 + sp.k3b * vbseff) * tmp2
        + t1_vth
        - dibl_sft;

    // Subthreshold slope n
    let tmp2_n = sp.nfactor * EPSSI / xdep;
    let tmp3_n = sp.cdsc + sp.cdscb * vbseff + sp.cdscd * vds_i;
    let tmp4_n = (tmp2_n + tmp3_n * theta0 + sp.cit) / cox;
    let n = if tmp4_n >= -0.5 {
        1.0 + tmp4_n
    } else {
        let t0 = 1.0 / (3.0 + 8.0 * tmp4_n);
        (1.0 + 3.0 * tmp4_n) * t0
    };
    let n = n.max(0.5);

    // Poly gate depletion
    let vfb_plus_phi = sp.vfbzb + sp.phi + sp.k1 * sp.sqrt_phi;
    let vgs_eff = if sp.ngate > 1e18 && sp.ngate < 1e25 && vgs_i > vfb_plus_phi {
        let t1 = 1e6 * CHARGE_Q * EPSSI * sp.ngate / (cox * cox);
        let t4 = (1.0 + 2.0 * (vgs_i - vfb_plus_phi) / t1).sqrt();
        let t2 = t1 * (t4 - 1.0);
        let t3 = 0.5 * t2 * t2 / t1;
        let t7 = 1.12 - t3 - 0.05;
        let t6 = (t7 * t7 + 0.224).sqrt();
        let t5 = 1.12 - 0.5 * (t7 + t6);
        vgs_i - t5
    } else {
        vgs_i
    };

    let vgst = vgs_eff - vth;

    // Vgsteff
    let t10 = 2.0 * n * vtm;
    let vgst_nvt = vgst / t10;
    let exp_arg = (2.0 * sp.voff - vgst) / t10;

    let vgsteff;
    if vgst_nvt > EXP_THRESHOLD {
        vgsteff = vgst;
    } else if exp_arg > EXP_THRESHOLD {
        let t0 = (vgst - sp.voff) / (n * vtm);
        let exp_vgst = t0.exp();
        vgsteff = vtm * sp.cdep0 / cox * exp_vgst;
    } else {
        let exp_vgst = safe_exp(vgst_nvt);
        let t1 = t10 * (1.0 + exp_vgst).ln();
        let d_t2_d_vg = -cox / (vtm * sp.cdep0) * safe_exp(exp_arg);
        let t2 = 1.0 - t10 * d_t2_d_vg;
        vgsteff = t1 / t2;
    }
    let vgsteff = vgsteff.max(1e-20);

    // Effective channel width
    let t9_w = sqrt_phis - sp.sqrt_phi;
    let weff_calc = (sp.weff - 2.0 * (sp.dwg * vgsteff + sp.dwb * t9_w)).max(2e-8);

    // Rds
    let t0_rds = sp.prwg * vgsteff + sp.prwb * t9_w;
    let rds = if t0_rds >= -0.9 {
        sp.rds0 * (1.0 + t0_rds)
    } else {
        let t1 = 1.0 / (17.0 + 20.0 * t0_rds);
        sp.rds0 * (0.8 + t0_rds) * t1
    };

    // Abulk
    let t1_ab = 0.5 * sp.k1ox / sqrt_phis;
    let t9_ab = (sp.xj * xdep).sqrt();
    let tmp1_ab = sp.leff + 2.0 * t9_ab;
    let t5_ab = sp.leff / tmp1_ab;
    let tmp2_ab = sp.a0 * t5_ab;
    let tmp3_ab = sp.weff + sp.b1;
    let tmp4_ab = sp.b0 / tmp3_ab;
    let t2_ab = tmp2_ab + tmp4_ab;

    let abulk0 = 1.0 + t1_ab * t2_ab;
    let t8_ab = sp.ags * sp.a0 * t5_ab * t5_ab * t5_ab;
    let mut abulk = abulk0 + (-t1_ab * t8_ab) * vgsteff;
    if abulk < 0.1 {
        let t9 = 1.0 / (3.0 - 20.0 * abulk);
        abulk = (0.2 - abulk) * t9;
    }
    let t2_keta = sp.keta * vbseff;
    let keta_factor = if t2_keta >= -0.9 {
        1.0 / (1.0 + t2_keta)
    } else {
        let t1 = 1.0 / (0.8 + t2_keta);
        (17.0 + 20.0 * t2_keta) * t1
    };
    abulk *= keta_factor;

    // Mobility
    let t5_mob = match model.mob_mod {
        1 => {
            let t0 = vgsteff + vth + vth;
            let t2 = sp.ua + sp.uc * vbseff;
            let t3 = t0 / model.tox;
            t3 * (t2 + sp.ub * t3)
        }
        2 => vgsteff / model.tox * (sp.ua + sp.uc * vbseff + sp.ub * vgsteff / model.tox),
        _ => {
            let t0 = vgsteff + vth + vth;
            let t2 = 1.0 + sp.uc * vbseff;
            let t3 = t0 / model.tox;
            t3 * (sp.ua + sp.ub * t3) * t2
        }
    };
    let denomi = if t5_mob >= -0.8 {
        1.0 + t5_mob
    } else {
        let t9 = 1.0 / (7.0 + 10.0 * t5_mob);
        (0.6 + t5_mob) * t9
    };
    let ueff = sp.u0temp / denomi;

    // Vdsat
    let wv_cox = weff_calc * sp.vsattemp * cox;
    let wv_cox_rds = wv_cox * rds;
    let esat = 2.0 * sp.vsattemp / ueff;
    let esat_l = esat * sp.leff;

    let vgst2vtm = vgsteff + 2.0 * vtm;
    let lambda = if sp.a1 == 0.0 {
        sp.a2
    } else if sp.a1 > 0.0 {
        let t0 = 1.0 - sp.a2;
        let t1 = t0 - sp.a1 * vgsteff - 0.0001;
        let t2 = (t1 * t1 + 0.0004 * t0).sqrt();
        sp.a2 + t0 - 0.5 * (t1 + t2)
    } else {
        let t1 = sp.a2 + sp.a1 * vgsteff - 0.0001;
        let t2 = (t1 * t1 + 0.0004 * sp.a2).sqrt();
        0.5 * (t1 + t2)
    };

    let vdsat = if rds == 0.0 && lambda == 1.0 {
        let t0 = 1.0 / (abulk * esat_l + vgst2vtm);
        esat_l * vgst2vtm * t0
    } else {
        let t9 = abulk * wv_cox_rds;
        let t7 = vgst2vtm * t9;
        let t0 = 2.0 * abulk * (t9 - 1.0 + 1.0 / lambda);
        let t1 = vgst2vtm * (2.0 / lambda - 1.0) + abulk * esat_l + 3.0 * t7;
        let t2 = vgst2vtm * (esat_l + 2.0 * vgst2vtm * wv_cox_rds);
        let t3 = (t1 * t1 - 2.0 * t0 * t2).max(0.0).sqrt();
        if t0.abs() > 1e-30 {
            (t1 - t3) / t0
        } else {
            esat_l * vgst2vtm / (abulk * esat_l + vgst2vtm)
        }
    };
    let vdsat = vdsat.max(1e-18);

    // Effective Vds
    let t1 = vdsat - vds_i - sp.delta;
    let t2 = (t1 * t1 + 4.0 * sp.delta * vdsat).sqrt();
    let _vdseff = if vds_i == 0.0 {
        0.0
    } else {
        let v = vdsat - 0.5 * (t1 + t2);
        v.min(vds_i)
    };

    // === Charge model dispatch (capMod 0/1/2/3) ===

    let cox_wl = cox * sp.weff_cv * sp.leff_cv;

    // Overlap charges: capMod 0/1 use simple linear, capMod 2/3 use smooth model
    let vgd = vgs - vds;
    let vgb = vgs - vbs;
    let (qgdo, qgso, qgb_charge) = if model.cap_mod <= 1 {
        // Simple linear overlap
        (sp.cgdo_eff * vgd, sp.cgso_eff * vgs, sp.cgbo_eff * vgb)
    } else {
        bsim3_overlap_charges(vgs, vgd, vgb, sp, model)
    };

    // AbulkCV (shared by capMod 1/2/3)
    let abulk_cv = abulk0 * sp.abulk_cv_factor;

    let (qgate_intr, qbulk_intr, qsrc) = match model.cap_mod {
        0 => {
            // capMod=0: Simple charge model (no Qac0/Qsub0, no centroid)
            let arg1 = vgs_eff - vbs_eff_cv - sp.vfbzb - vgsteff;
            if arg1 <= 0.0 {
                let qg = cox_wl * arg1;
                (qg, -qg, 0.0)
            } else {
                let t0 = 0.5 * sp.k1ox;
                let t1 = (t0 * t0 + arg1).sqrt();
                let qg_acc = cox_wl * sp.k1ox * (t1 - t0);
                let qb_acc = -qg_acc;

                let vdsat_cv = vgsteff / abulk_cv;
                let (qg_ch, qb_ch, qs) = if vdsat_cv < vds_i {
                    // Saturation
                    let qg = cox_wl * (vgsteff - vdsat_cv / 3.0);
                    let qb = cox_wl / 3.0 * (vdsat_cv - vgsteff);
                    let qs = if model.xpart > 0.5 {
                        -2.0 / 3.0 * cox_wl * vgsteff
                    } else if model.xpart < 0.5 {
                        -0.4 * cox_wl * vgsteff
                    } else {
                        -cox_wl / 3.0 * vgsteff
                    };
                    (qg, qb, qs)
                } else {
                    // Linear
                    let t0 = abulk_cv * vds_i;
                    let t1 = 12.0 * (vgsteff - 0.5 * t0 + 1e-20);
                    let t2 = vds_i / t1;
                    let t3 = t0 * t2;
                    let qg = cox_wl * (vgsteff - 0.5 * vds_i + t3);
                    let qb = cox_wl * (1.0 - abulk_cv) * (0.5 * vds_i - t3);
                    let qs = charge_partition(model.xpart, cox_wl, vgsteff, t0, t1);
                    (qg, qb, qs)
                };
                (qg_acc + qg_ch, qb_acc + qb_ch, qs)
            }
        }
        1 => {
            // capMod=1: k1ox model with charge partition
            let arg1 = vgs_eff - vbs_eff_cv - sp.vfbzb - vgsteff;
            if arg1 <= 0.0 {
                let qg = cox_wl * arg1;
                (qg, -qg, 0.0)
            } else {
                let t0 = 0.5 * sp.k1ox;
                let t1 = (t0 * t0 + arg1).sqrt();
                let qg_acc = cox_wl * sp.k1ox * (t1 - t0);
                let qb_acc = -qg_acc;

                let vdsat_cv = vgsteff / abulk_cv;
                let (qg_ch, qb_ch, qs) = if vdsat_cv < vds_i {
                    let qg = cox_wl * (vgsteff - vdsat_cv / 3.0);
                    let qb = cox_wl / 3.0 * (vdsat_cv - vgsteff);
                    let qs = if model.xpart > 0.5 {
                        -2.0 / 3.0 * cox_wl * vgsteff
                    } else if model.xpart < 0.5 {
                        -0.4 * cox_wl * vgsteff
                    } else {
                        -cox_wl / 3.0 * vgsteff
                    };
                    (qg, qb, qs)
                } else {
                    let t0 = abulk_cv * vds_i;
                    let t1 = 12.0 * (vgsteff - 0.5 * t0 + 1e-20);
                    let t2 = vds_i / t1;
                    let t3 = t0 * t2;
                    let qg = cox_wl * (vgsteff - 0.5 * vds_i + t3);
                    let qb = cox_wl * (1.0 - abulk_cv) * (0.5 * vds_i - t3);
                    let qs = charge_partition(model.xpart, cox_wl, vgsteff, t0, t1);
                    (qg, qb, qs)
                };
                (qg_acc + qg_ch, qb_acc + qb_ch, qs)
            }
        }
        2 => {
            // capMod=2: Qac0/Qsub0 corrections, no centroid
            let (_vfbeff, qac0, qsub0) =
                compute_qac0_qsub0(vgs_eff, vbs_eff_cv, vgsteff, cox_wl, sp);

            let vdsat_cv = vgsteff / abulk_cv;
            let v4_cv = vdsat_cv - vds_i - DELTA_4;
            let t0_cv = (v4_cv * v4_cv + 4.0 * DELTA_4 * vdsat_cv).sqrt();
            let vdseff_cv = if vds_i == 0.0 {
                0.0
            } else {
                (vdsat_cv - 0.5 * (v4_cv + t0_cv)).min(vds_i).max(0.0)
            };

            let t0_q = abulk_cv * vdseff_cv;
            let t1_q = 12.0 * (vgsteff - 0.5 * t0_q + 1e-20);
            let t3_q = t0_q * vdseff_cv / t1_q;

            let qgate = cox_wl * (vgsteff - 0.5 * vdseff_cv + t3_q);
            let qbulk = cox_wl * (1.0 - abulk_cv) * (0.5 * vdseff_cv - t3_q);
            let qsrc = charge_partition(model.xpart, cox_wl, vgsteff, t0_q, t1_q);

            let qg_intr = qgate + qac0 + qsub0 - qbulk;
            let qb_intr = qbulk - (qac0 + qsub0);
            (qg_intr, qb_intr, qsrc)
        }
        _ => {
            // capMod=3 (default): Centroid model with Qac0/Qsub0 + DeltaPhi
            let tox_ang = 1.0e8 * model.tox;

            // Vfbeff and Qac0/Qsub0 with centroid-corrected CoxWL
            let (vfbeff, _) = compute_vfbeff(vgs_eff, vbs_eff_cv, sp);

            // Centroid thickness (depletion)
            let t0 = (vgs_eff - vbs_eff_cv - sp.vfbzb) / tox_ang;
            let tmp = t0 * sp.acde;
            let tcen = if tmp > -EXP_THRESHOLD && tmp < EXP_THRESHOLD {
                sp.ldeb * tmp.exp()
            } else if tmp <= -EXP_THRESHOLD {
                sp.ldeb * MIN_EXP
            } else {
                sp.ldeb * MAX_EXP
            };
            let link = 1.0e-3 * model.tox;
            let v3l = sp.ldeb - tcen - link;
            let v4l = (v3l * v3l + 4.0 * link * sp.ldeb).sqrt();
            let tcen = sp.ldeb - 0.5 * (v3l + v4l);

            let ccen = EPSSI / tcen;
            let t2 = cox / (cox + ccen);
            let coxeff = t2 * ccen;
            let cox_wl_cen = cox_wl * coxeff / cox;

            let qac0 = cox_wl_cen * (vfbeff - sp.vfbzb);

            let t0_k = 0.5 * sp.k1ox;
            let t3_sub = vgs_eff - vfbeff - vbs_eff_cv - vgsteff;
            let t1_sub = if sp.k1ox == 0.0 {
                0.0
            } else if t3_sub < 0.0 {
                t0_k + t3_sub / sp.k1ox
            } else {
                (t0_k * t0_k + t3_sub).sqrt()
            };
            let qsub0 = cox_wl_cen * sp.k1ox * (t1_sub - t0_k);

            // DeltaPhi
            let denomi_dp = if sp.k1ox <= 0.0 {
                0.25 * sp.moin * vtm
            } else {
                sp.moin * vtm * sp.k1ox * sp.k1ox
            };
            let t0_dp = if sp.k1ox <= 0.0 {
                0.5 * sp.sqrt_phi
            } else {
                sp.k1ox * sp.sqrt_phi
            };
            let t1_dp = 2.0 * t0_dp + vgsteff;
            let delta_phi = vtm * (1.0 + t1_dp * vgsteff / denomi_dp).ln();

            // Second centroid (inversion layer)
            let t3_inv = 4.0 * (vth - sp.vfbzb - sp.phi);
            let tox_2 = tox_ang * 2.0;
            let t0_inv = if t3_inv >= 0.0 {
                (vgsteff + t3_inv) / tox_2
            } else {
                (vgsteff + 1.0e-20) / tox_2
            };
            let tmp_inv = t0_inv.powf(0.7);
            let tcen2 = 1.9e-9 / (1.0 + tmp_inv);
            let ccen2 = EPSSI / tcen2;
            let t0_c2 = cox / (cox + ccen2);
            let coxeff2 = t0_c2 * ccen2;
            let cox_wl_cen2 = cox_wl * coxeff2 / cox;

            // VdsatCV with DeltaPhi
            let vdsat_cv = (vgsteff - delta_phi) / abulk_cv;
            let v4_cv = vdsat_cv - vds_i - DELTA_4;
            let t0_cv = (v4_cv * v4_cv + 4.0 * DELTA_4 * vdsat_cv).sqrt();
            let vdseff_cv = if vds_i == 0.0 {
                0.0
            } else {
                (vdsat_cv - 0.5 * (v4_cv + t0_cv)).min(vds_i).max(0.0)
            };

            let t0_q = abulk_cv * vdseff_cv;
            let t1_q = vgsteff - delta_phi;
            let t2_q = 12.0 * (t1_q - 0.5 * t0_q + 1.0e-20);
            let t3_q = t0_q / t2_q;

            let qgate = cox_wl_cen2 * (t1_q - t0_q * (0.5 - t3_q));
            let t7_q = 1.0 - abulk_cv;
            let qbulk = cox_wl_cen2 * t7_q * (0.5 * vdseff_cv - t0_q * vdseff_cv / t2_q);

            let qsrc = if model.xpart > 0.5 {
                -cox_wl_cen2 * (t1_q / 2.0 + t0_q / 4.0 - 0.5 * t0_q * t0_q / t2_q)
            } else if model.xpart < 0.5 {
                let t2_s = t2_q / 12.0;
                let t3_s = 0.5 * cox_wl_cen2 / (t2_s * t2_s);
                let t4_s = t1_q * (2.0 * t0_q * t0_q / 3.0 + t1_q * (t1_q - 4.0 * t0_q / 3.0))
                    - 2.0 * t0_q * t0_q * t0_q / 15.0;
                -t3_s * t4_s
            } else {
                -0.5 * qgate
            };

            let qg_intr = qgate + qac0 + qsub0 - qbulk;
            let qb_intr = qbulk - (qac0 + qsub0);
            (qg_intr, qb_intr, qsrc)
        }
    };

    let qg_total = qgate_intr + qgdo + qgso + qgb_charge;
    let qb_total = qbulk_intr - qgb_charge;
    let qs_total = qsrc - qgso;

    (qg_total, qb_total, qs_total)
}

/// Compute charge partition for source charge.
/// Used by capMod 0/1/2 in the linear region.
fn charge_partition(xpart: f64, cox_wl: f64, vgsteff: f64, t0: f64, t1: f64) -> f64 {
    if xpart > 0.5 {
        // 0/100 partition
        let t1_2 = t1 + t1;
        -cox_wl * (0.5 * vgsteff + 0.25 * t0 - t0 * t0 / t1_2)
    } else if xpart < 0.5 {
        // 40/60 partition
        let t1_s = t1 / 12.0;
        let t2_s = 0.5 * cox_wl / (t1_s * t1_s);
        let t3_s = vgsteff * (2.0 * t0 * t0 / 3.0 + vgsteff * (vgsteff - 4.0 * t0 / 3.0))
            - 2.0 * t0 * t0 * t0 / 15.0;
        -t2_s * t3_s
    } else {
        // 50/50 partition: qsrc = -0.5 * (qgate + qbulk)
        -0.5 * cox_wl * vgsteff
    }
}

/// Compute Vfbeff (effective flat-band voltage).
/// Returns (vfbeff, d_vfbeff_d_vg).
fn compute_vfbeff(vgs_eff: f64, vbs_eff_cv: f64, sp: &Bsim3SizeParam) -> (f64, f64) {
    let v3 = sp.vfbzb - vgs_eff + vbs_eff_cv - DELTA_3;
    if sp.vfbzb <= 0.0 {
        let t0 = (v3 * v3 - 4.0 * DELTA_3 * sp.vfbzb).sqrt();
        let t1 = 0.5 * (1.0 + v3 / t0);
        (sp.vfbzb - 0.5 * (v3 + t0), t1)
    } else {
        let t0 = (v3 * v3 + 4.0 * DELTA_3 * sp.vfbzb).sqrt();
        let t1 = 0.5 * (1.0 + v3 / t0);
        (sp.vfbzb - 0.5 * (v3 + t0), t1)
    }
}

/// Compute Qac0 and Qsub0 for capMod=2 (without centroid).
/// Returns (vfbeff, qac0, qsub0).
fn compute_qac0_qsub0(
    vgs_eff: f64,
    vbs_eff_cv: f64,
    vgsteff: f64,
    cox_wl: f64,
    sp: &Bsim3SizeParam,
) -> (f64, f64, f64) {
    let (vfbeff, _) = compute_vfbeff(vgs_eff, vbs_eff_cv, sp);
    let qac0 = cox_wl * (vfbeff - sp.vfbzb);

    let t0_k = 0.5 * sp.k1ox;
    let t3_sub = vgs_eff - vfbeff - vbs_eff_cv - vgsteff;
    let t1_sub = if sp.k1ox == 0.0 {
        0.0
    } else if t3_sub < 0.0 {
        t0_k + t3_sub / sp.k1ox
    } else {
        (t0_k * t0_k + t3_sub).sqrt()
    };
    let qsub0 = cox_wl * sp.k1ox * (t1_sub - t0_k);

    (vfbeff, qac0, qsub0)
}

/// Compute junction capacitances for BSIM3.
/// Returns (capbs, capbd).
pub fn bsim3_junction_caps(
    vbs: f64,
    vbd: f64,
    model: &Bsim3Model,
    ad: f64,
    as_: f64,
    pd: f64,
    ps: f64,
) -> (f64, f64) {
    let temp_ratio = TEMP_DEFAULT / (model.tnom + 273.15) - 1.0;

    // Temperature-adjusted junction parameters
    let pb_eff = model.pb - model.tpb * temp_ratio;
    let cj_eff = model.cj * (1.0 + model.tcj * temp_ratio);
    let pbsw_eff = model.pbsw - model.tpbsw * temp_ratio;
    let cjsw_eff = model.cjsw * (1.0 + model.tcjsw * temp_ratio);
    let pbswg_eff = model.pbswg - model.tpbswg * temp_ratio;
    let cjswg_eff = model.cjswg * (1.0 + model.tcjswg * temp_ratio);

    // Zero-bias capacitances
    let czbs = cj_eff * as_;
    let czbssw = cjsw_eff * ps;
    let czbsswg = cjswg_eff * ps;

    let czbd = cj_eff * ad;
    let czbdsw = cjsw_eff * pd;
    let czbdswg = cjswg_eff * pd;

    // Source junction cap
    let capbs = junction_cap_bias(
        vbs,
        czbs,
        czbssw,
        czbsswg,
        pb_eff,
        pbsw_eff,
        pbswg_eff,
        model.mj,
        model.mjsw,
        model.mjswg,
    );

    // Drain junction cap
    let capbd = junction_cap_bias(
        vbd,
        czbd,
        czbdsw,
        czbdswg,
        pb_eff,
        pbsw_eff,
        pbswg_eff,
        model.mj,
        model.mjsw,
        model.mjswg,
    );

    (capbs, capbd)
}

/// Compute bias-dependent junction capacitance.
#[allow(clippy::too_many_arguments)]
fn junction_cap_bias(
    v: f64,
    cz: f64,
    czsw: f64,
    czswg: f64,
    pb: f64,
    pbsw: f64,
    pbswg: f64,
    mj: f64,
    mjsw: f64,
    mjswg: f64,
) -> f64 {
    if cz + czsw + czswg <= 0.0 {
        return 0.0;
    }

    if v == 0.0 {
        return cz + czsw + czswg;
    }

    if v < 0.0 {
        // Reverse bias: depletion capacitance
        let mut cap = 0.0;
        if cz > 0.0 && pb > 0.0 {
            cap += cz * (1.0 - v / pb).powf(-mj);
        }
        if czsw > 0.0 && pbsw > 0.0 {
            cap += czsw * (1.0 - v / pbsw).powf(-mjsw);
        }
        if czswg > 0.0 && pbswg > 0.0 {
            cap += czswg * (1.0 - v / pbswg).powf(-mjswg);
        }
        cap
    } else {
        // Forward bias: linear approximation
        let t0 = cz + czsw + czswg;
        let t1 = if pb > 0.0 { cz * mj / pb } else { 0.0 }
            + if pbsw > 0.0 { czsw * mjsw / pbsw } else { 0.0 }
            + if pbswg > 0.0 {
                czswg * mjswg / pbswg
            } else {
                0.0
            };
        t0 + t1 * v
    }
}

/// Stamp BSIM3 companion model into the MNA matrix and RHS.
pub fn stamp_bsim3(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &Bsim3Instance,
    comp: &Bsim3Companion,
) {
    let dp = inst.drain_eff_idx();
    let g = inst.gate_idx;
    let sp = inst.source_eff_idx();
    let b = inst.bulk_idx;

    let sign = inst.model.mos_type.sign();
    let m = inst.m;

    let (xnrm, xrev) = if comp.mode > 0 {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };

    let fwd_sum = comp.gm + comp.gmbs;
    let rev_sum = 0.0;
    let (fwd_sum, rev_sum) = if comp.mode > 0 {
        (fwd_sum, rev_sum)
    } else {
        (0.0, fwd_sum)
    };

    // gds between d' and s'
    crate::stamp_conductance(matrix, dp, sp, m * comp.gds);

    // gm VCCS
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

    // gmbs
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

    // gbd between b and d'
    crate::stamp_conductance(matrix, b, dp, m * comp.gbd);
    // gbs between b and s'
    crate::stamp_conductance(matrix, b, sp, m * comp.gbs);

    // Series resistances (from RSH*NRD/NRS)
    let g_drain = inst.drain_conductance();
    if g_drain > 0.0 {
        crate::stamp_conductance(matrix, inst.drain_idx, dp, m * g_drain);
    }
    let g_source = inst.source_conductance();
    if g_source > 0.0 {
        crate::stamp_conductance(matrix, inst.source_idx, sp, m * g_source);
    }

    // FwdSum/RevSum stamps on diagonal
    if let (Some(d), Some(s)) = (dp, sp) {
        // DPdp += RevSum, DPsp -= FwdSum
        matrix.add(d, d, m * rev_sum);
        matrix.add(d, s, -m * fwd_sum);
        matrix.add(s, s, m * fwd_sum);
        matrix.add(s, d, -m * rev_sum);
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

/// BSIM3 voltage limiting for NR convergence.
pub fn bsim3_limit(
    vgs_new: f64,
    vds_new: f64,
    vbs_new: f64,
    vgs_old: f64,
    vds_old: f64,
    vbs_old: f64,
    vth: f64,
) -> (f64, f64, f64) {
    // Use ngspice DEVfetlim for VGS and VDS
    let vgs = fetlim(vgs_new, vgs_old, vth);
    let vds = fetlim(vds_new, vds_old, vth);
    // Use pnjlim for VBS
    let vbs = pnjlim(vbs_new, vbs_old);
    (vgs, vds, vbs)
}

fn fetlim(vnew: f64, vold: f64, vto: f64) -> f64 {
    let vtsthi = (2.0 * vold.abs()).max(1.0);
    let vtstlo = vtsthi / 2.0 + 0.1;
    let vtox = vto + 3.5;
    let delv = vnew - vold;

    if vold >= vtox {
        if delv.abs() >= vtsthi {
            if vold < vnew {
                vold + vtsthi
            } else {
                vold - vtsthi
            }
        } else {
            vnew
        }
    } else if vold >= 0.0 {
        if delv.abs() >= vtstlo {
            if vold < vnew {
                vold + vtstlo
            } else {
                vold - vtstlo
            }
        } else {
            vnew
        }
    } else {
        // Below threshold
        if delv.abs() >= vtsthi {
            if vold < vnew {
                vold + vtsthi
            } else {
                vold - vtsthi
            }
        } else {
            vnew
        }
    }
}

fn pnjlim(vnew: f64, vold: f64) -> f64 {
    let vt = KBOQ * TEMP_DEFAULT;
    let vcrit = vt * (vt / (std::f64::consts::SQRT_2 * 1e-14)).ln();
    if vnew > vcrit && (vnew - vold).abs() > 2.0 * vt {
        if vold > 0.0 {
            let arg = (vnew - vold) / vt;
            if arg > 0.0 {
                vold + vt * (1.0 + arg.ln())
            } else {
                vold - vt * (1.0 + (-arg).ln())
            }
        } else {
            vt * (vnew / vt).ln()
        }
    } else {
        vnew
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ferrospice_netlist::{Expr, ModelDef, Param};

    #[test]
    fn test_bsim3_model_defaults() {
        let m = Bsim3Model::new(MosfetType::Nmos);
        assert_eq!(m.mob_mod, 1);
        assert_eq!(m.cap_mod, 3);
        assert_abs_diff_eq!(m.tox, 150e-10, epsilon = 1e-25);
        assert_abs_diff_eq!(m.vth0, 0.7, epsilon = 1e-15);
        assert_abs_diff_eq!(m.vsat, 8e4, epsilon = 1e-10);
    }

    #[test]
    fn test_bsim3_model_parse() {
        let model_def = ModelDef {
            name: "NMOD".to_string(),
            kind: "NMOS".to_string(),
            params: vec![
                Param {
                    name: "VTH0".to_string(),
                    value: Expr::Num(0.5),
                },
                Param {
                    name: "TOX".to_string(),
                    value: Expr::Num(1.5e-8),
                },
                Param {
                    name: "U0".to_string(),
                    value: Expr::Num(670.0),
                },
                Param {
                    name: "VSAT".to_string(),
                    value: Expr::Num(8e4),
                },
                Param {
                    name: "K1".to_string(),
                    value: Expr::Num(0.5),
                },
            ],
        };
        let m = Bsim3Model::from_model_def(&model_def);
        assert_eq!(m.mos_type, MosfetType::Nmos);
        assert_abs_diff_eq!(m.vth0, 0.5, epsilon = 1e-15);
        assert_abs_diff_eq!(m.tox, 1.5e-8, epsilon = 1e-25);
        assert!(m.vth0_given);
        assert!(m.k1_given);
    }

    #[test]
    fn test_size_dep_params() {
        let m = Bsim3Model::new(MosfetType::Nmos);
        let sp = m.size_dep_param(10e-6, 1e-6, TEMP_DEFAULT);
        assert!(sp.leff > 0.0);
        assert!(sp.weff > 0.0);
        assert!(sp.phi > 0.0);
        assert!(sp.u0temp > 0.0);
        assert!(sp.vsattemp > 0.0);
    }

    #[test]
    fn test_bsim3_companion_cutoff() {
        let m = Bsim3Model::new(MosfetType::Nmos);
        let sp = m.size_dep_param(10e-6, 1e-6, TEMP_DEFAULT);
        // Vgs well below threshold
        let comp = bsim3_companion(0.0, 1.0, 0.0, &sp, &m);
        // In deep subthreshold, ids should be very small
        assert!(
            comp.ids < 1e-6,
            "Ids in cutoff should be very small: {}",
            comp.ids
        );
    }

    #[test]
    fn test_bsim3_companion_saturation() {
        let m = Bsim3Model::new(MosfetType::Nmos);
        let sp = m.size_dep_param(10e-6, 1e-6, TEMP_DEFAULT);
        // Vgs = 1.8V, Vds = 1.8V — well in saturation
        let comp = bsim3_companion(1.8, 1.8, 0.0, &sp, &m);
        assert!(
            comp.ids > 1e-6,
            "Should have significant current in saturation: {}",
            comp.ids
        );
        assert!(comp.gm > 0.0, "gm should be positive in saturation");
        assert!(comp.gds > 0.0, "gds should be positive");
    }

    #[test]
    fn test_bsim3_companion_linear() {
        let m = Bsim3Model::new(MosfetType::Nmos);
        let sp = m.size_dep_param(10e-6, 1e-6, TEMP_DEFAULT);
        // Vgs = 1.8V, Vds = 0.1V — linear region
        let comp = bsim3_companion(1.8, 0.1, 0.0, &sp, &m);
        assert!(comp.ids > 0.0, "Should have current in linear region");
        assert!(comp.gds > 0.0, "gds should be positive in linear region");
        assert!(comp.gm > 0.0, "gm should be positive in linear region");
        // Linear region current should be less than saturation current
        let comp_sat = bsim3_companion(1.8, 1.8, 0.0, &sp, &m);
        assert!(comp.ids < comp_sat.ids, "Linear Ids < saturation Ids");
    }

    #[test]
    fn test_bsim3_charges_and_caps() {
        let m = Bsim3Model::new(MosfetType::Nmos);
        let sp = m.size_dep_param(10e-6, 1e-6, TEMP_DEFAULT);
        // Test charge function doesn't return NaN
        let (qg, qb, qs) = bsim3_total_charges(1.8, 1.8, 0.0, &sp, &m);
        assert!(!qg.is_nan(), "qg is NaN");
        assert!(!qb.is_nan(), "qb is NaN");
        assert!(!qs.is_nan(), "qs is NaN");
        assert!(!qg.is_infinite(), "qg is infinite");
        // Gate charge should be positive (inversion layer)
        eprintln!("qg={qg:.4e}, qb={qb:.4e}, qs={qs:.4e}");

        // Test companion capacitances are finite
        let comp = bsim3_companion(1.8, 1.8, 0.0, &sp, &m);
        eprintln!(
            "cggb={:.4e}, cgdb={:.4e}, cgsb={:.4e}",
            comp.cggb, comp.cgdb, comp.cgsb
        );
        eprintln!(
            "cbgb={:.4e}, cbdb={:.4e}, cbsb={:.4e}",
            comp.cbgb, comp.cbdb, comp.cbsb
        );
        eprintln!(
            "cdgb={:.4e}, cddb={:.4e}, cdsb={:.4e}",
            comp.cdgb, comp.cddb, comp.cdsb
        );
        assert!(!comp.cggb.is_nan(), "cggb is NaN");
        assert!(!comp.cgsb.is_nan(), "cgsb is NaN");
        assert!(comp.cggb > 0.0, "cggb should be positive: {}", comp.cggb);

        // Test at VDS=0 (initial NR guess)
        let (qg0, qb0, qs0) = bsim3_total_charges(0.0, 0.0, 0.0, &sp, &m);
        assert!(!qg0.is_nan(), "qg at origin is NaN");
        assert!(!qb0.is_nan(), "qb at origin is NaN");
        assert!(!qs0.is_nan(), "qs at origin is NaN");
        eprintln!("At origin: qg={qg0:.4e}, qb={qb0:.4e}, qs={qs0:.4e}");

        let comp0 = bsim3_companion(0.0, 0.0, 0.0, &sp, &m);
        assert!(!comp0.cggb.is_nan(), "cggb at origin is NaN");
        assert!(!comp0.ids.is_nan(), "ids at origin is NaN: {}", comp0.ids);
        eprintln!("At origin: ids={:.4e}, cggb={:.4e}", comp0.ids, comp0.cggb);
    }
}
