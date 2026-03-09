//! BSIM3SOI-DD (Double Diffused Silicon-On-Insulator) MOSFET model.
//!
//! Implements the BSIM3SOI-DD v2.0 model matching ngspice level 56.
//! DD is a hybrid of PD and FD: it uses the FD-style self-consistent surface
//! potential chain (Vbs0t→Vbs0→Vbs0mos→Vthfd→Vbs0eff→Vbsmos→Vbseff) combined
//! with the PD-style 4-component junction diode model and GIDL currents.
//! Impact ionization uses ALPHA0/ALPHA1/BETA0 + AII/BII/CII/DII parameters.

#![allow(unused_variables, dead_code, clippy::too_many_arguments, unused_parens)]

use ferrospice_netlist::{Expr, ModelDef};

use crate::mosfet::MosfetType;

// Physical constants
const EPSOX: f64 = 3.453133e-11;
const EPSSI: f64 = 1.03594e-10;
const CHARGE_Q: f64 = 1.60219e-19;
const KBOQ: f64 = 8.617087e-5;
const EG300: f64 = 1.115;
const EXP_THRESHOLD: f64 = 34.0;
const MAX_EXP: f64 = 5.834617425e14;
const MIN_EXP: f64 = 1.713908431e-15;
const EXPL_THRESHOLD: f64 = 100.0;
const MAX_EXPL: f64 = 2.688117142e43;
const MIN_EXPL: f64 = 3.720075976e-44;
const DELTA_1: f64 = 0.02;
const DELTA_4: f64 = 0.02;
const DELT_VBS0EFF: f64 = 0.02;
const DELT_VBSMOS: f64 = 0.005;
const DELT_VBSEFF: f64 = 0.005;

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

/// Safe exponential for SOI (larger threshold than BSIM3).
fn soi_dexp(x: f64) -> (f64, f64) {
    if x > EXPL_THRESHOLD {
        (MAX_EXPL * (1.0 + x - EXPL_THRESHOLD), MAX_EXPL)
    } else if x < -EXPL_THRESHOLD {
        (MIN_EXPL, 0.0)
    } else {
        let e = x.exp();
        (e, e)
    }
}

/// BSIM3SOI-DD model parameters (from .model card, Level=56).
#[derive(Debug, Clone)]
pub struct Bsim3SoiDdModel {
    pub mos_type: MosfetType,

    // Mode selection
    pub mob_mod: i32,
    pub cap_mod: i32,
    pub sh_mod: i32,

    // Oxide / SOI geometry
    pub tox: f64,
    pub tsi: f64,
    pub tbox: f64,

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
    pub pdiblc1: f64,
    pub pdiblc2: f64,
    pub pdiblcb: f64,
    pub drout: f64,
    pub pvag: f64,
    pub delta: f64,

    // Series resistance
    pub rdsw: f64,
    pub prwg: f64,
    pub prwb: f64,
    pub prt: f64,
    pub wr: f64,

    // Width/length effects
    pub dwg: f64,
    pub dwb: f64,
    pub b0: f64,
    pub b1: f64,
    pub wint: f64,
    pub lint: f64,
    pub dlc: f64,
    pub dwc: f64,

    // Impact ionization (DD uses ALPHA0/ALPHA1/BETA0 + AII/BII/CII/DII)
    pub alpha0: f64,
    pub alpha1: f64,
    pub beta0: f64,
    pub aii: f64,
    pub bii: f64,
    pub cii: f64,
    pub dii: f64,

    // Temperature
    pub tnom: f64,
    pub kt1: f64,
    pub kt1l: f64,
    pub kt2: f64,

    // SOI junction model (same as PD: 4-component)
    pub ndiode: f64,
    pub ntun: f64,
    pub nrecf0: f64,
    pub nrecr0: f64,
    pub vrec0: f64,
    pub ntrecf: f64,
    pub ntrecr: f64,
    pub isbjt: f64,
    pub isdif: f64,
    pub istun: f64,
    pub isrec: f64,
    pub xbjt: f64,
    pub xdif: f64,
    pub xrec: f64,
    pub xtun: f64,
    pub ahli: f64,
    pub lbjt0: f64,
    pub ln: f64,
    pub nbjt: f64,
    pub ndif: f64,
    pub aely: f64,
    pub vabjt: f64,

    // GIDL (same as PD)
    pub agidl: f64,
    pub bgidl: f64,
    pub ngidl: f64,

    // DD-specific surface potential params (shared with FD)
    pub kb1: f64,
    pub kb3: f64,
    pub dvbd0: f64,
    pub dvbd1: f64,
    pub vbsa: f64,
    pub delp: f64,
    pub abp: f64,
    pub mxc: f64,
    pub adice0: f64,
    pub kbjt1: f64,
    pub edl: f64,

    // Body resistance
    pub rbody: f64,
    pub rbsh: f64,

    // Gate overlap
    pub cgso: f64,
    pub cgdo: f64,

    // CV model
    pub clc: f64,
    pub cle: f64,
    pub cf: f64,
    pub ckappa: f64,
    pub cgdl: f64,
    pub cgsl: f64,

    // Junction capacitance
    pub cjswg: f64,
    pub mjswg: f64,
    pub pbswg: f64,
    pub tt: f64,
    pub csdesw: f64,
    pub asd: f64,

    // Self-heating
    pub rth0: f64,
    pub cth0: f64,

    // Precomputed
    pub cox: f64,
    pub vtm: f64,
    pub phi: f64,
    pub sqrt_phi: f64,
    pub vbi_default: f64,
    pub factor1: f64,
    pub ni: f64,
    pub eg: f64,
    // DD-specific precomputed (from FD)
    pub cbox: f64,
    pub csi: f64,
    pub qsi: f64,
    pub csieff: f64,
    pub qsieff: f64,
    pub vfbb: f64,
}

/// Size-dependent parameters for BSIM3SOI-DD.
#[derive(Debug, Clone)]
pub struct Bsim3SoiDdSizeParam {
    pub leff: f64,
    pub weff: f64,
    pub leff_cv: f64,
    pub weff_cv: f64,

    // Core parameters
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
    pub eta0: f64,
    pub etab: f64,
    pub dsub: f64,
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
    pub vsat: f64,
    pub a0: f64,
    pub ags: f64,
    pub a1: f64,
    pub a2: f64,
    pub at: f64,
    pub keta: f64,
    pub pclm: f64,
    pub pdiblc1: f64,
    pub pdiblc2: f64,
    pub pdiblcb: f64,
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
    pub beta0: f64,
    pub kt1: f64,
    pub kt1l: f64,
    pub kt2: f64,
    pub ndiode: f64,
    pub agidl: f64,
    pub bgidl: f64,
    pub ngidl: f64,

    // Precomputed
    pub phi: f64,
    pub sqrt_phi: f64,
    pub xdep0: f64,
    pub litl: f64,
    pub theta0vb0: f64,
    pub theta_rout: f64,
    pub k1eff: f64,
    pub npeak: f64,
    pub nsub: f64,
    pub vfb: f64,
    pub u0temp: f64,
    pub vsattemp: f64,
    pub rds0: f64,
    pub cdep0: f64,
    pub vbi: f64,
    pub lratio: f64,
    pub lratiodif: f64,
    pub vearly: f64,
    pub arfabjt: f64,
    pub wdios: f64,
    pub wdiod: f64,

    // Junction precomputed
    pub jbjt: f64,
    pub jdif: f64,
    pub jrec: f64,
    pub jtun: f64,

    // Overlap caps
    pub cgso_eff: f64,
    pub cgdo_eff: f64,

    // DD-specific
    pub dvbd0: f64,
    pub dvbd1: f64,

    pub nseg: f64,
}

/// BSIM3SOI-DD instance with node indices.
#[derive(Debug, Clone)]
pub struct Bsim3SoiDdInstance {
    pub name: String,
    pub drain_idx: Option<usize>,
    pub gate_idx: Option<usize>,
    pub source_idx: Option<usize>,
    /// Back-gate (E) node.
    pub e_idx: Option<usize>,
    /// External body contact (B/P) node — optional.
    pub body_idx: Option<usize>,
    pub drain_prime_idx: Option<usize>,
    pub source_prime_idx: Option<usize>,
    /// Internal body node (always created for SOI).
    pub body_int_idx: Option<usize>,
    pub w: f64,
    pub l: f64,
    pub m: f64,
    pub nrd: f64,
    pub nrs: f64,
    pub model: Bsim3SoiDdModel,
    pub size_params: Bsim3SoiDdSizeParam,
    pub vth0_inst: f64,
    pub nbc: f64,
}

/// NR companion result for BSIM3SOI-DD.
#[derive(Debug, Clone)]
pub struct Bsim3SoiDdCompanion {
    pub ids: f64,
    pub gm: f64,
    pub gds: f64,
    pub gmbs: f64,
    pub mode: i32,
    pub vdsat: f64,

    // Junction currents and conductances (body node KCL)
    pub ibs: f64,
    pub ibd: f64,
    pub gbs_jct: f64,
    pub gbd_jct: f64,

    // Impact ionization
    pub iii: f64,
    pub gii_d: f64,
    pub gii_g: f64,
    pub gii_b: f64,

    // GIDL
    pub igidl: f64,
    pub ggidl_d: f64,
    pub ggidl_g: f64,
    pub isgidl: f64,
    pub gsgidl_g: f64,

    // Body current components for body node KCL
    pub ceq_d: f64,
    pub ceq_bs: f64,
    pub ceq_bd: f64,

    // Capacitances (intrinsic)
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

    // DD-computed body-source voltage (for body node feedback in floating body)
    pub vbs_dd: f64,
}

impl Bsim3SoiDdModel {
    pub fn new(mos_type: MosfetType) -> Self {
        let vth0_default = match mos_type {
            MosfetType::Nmos => 0.7,
            MosfetType::Pmos => -0.7,
        };
        let u0_default = match mos_type {
            MosfetType::Nmos => 0.067,
            MosfetType::Pmos => 0.025,
        };
        let mut m = Self {
            mos_type,
            mob_mod: 1,
            cap_mod: 3,
            sh_mod: 0,
            tox: 150e-10,
            tsi: 1e-7,
            tbox: 8e-8,
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
            pdiblc1: 0.39,
            pdiblc2: 0.0086,
            pdiblcb: 0.0,
            drout: 0.56,
            pvag: 0.0,
            delta: 0.01,
            rdsw: 0.0,
            prwg: 0.0,
            prwb: 0.0,
            prt: 0.0,
            wr: 1.0,
            dwg: 0.0,
            dwb: 0.0,
            b0: 0.0,
            b1: 0.0,
            wint: 0.0,
            lint: 0.0,
            dlc: 0.0,
            dwc: 0.0,
            alpha0: 0.0,
            alpha1: 0.0,
            beta0: 30.0,
            aii: 0.0,
            bii: 0.0,
            cii: 0.0,
            dii: 0.0,
            tnom: 27.0,
            kt1: -0.11,
            kt1l: 0.0,
            kt2: 0.022,
            ndiode: 1.0,
            ntun: 10.0,
            nrecf0: 2.0,
            nrecr0: 10.0,
            vrec0: 0.0,
            ntrecf: 0.0,
            ntrecr: 0.0,
            isbjt: 1e-6,
            isdif: 0.0,
            istun: 0.0,
            isrec: 0.0,
            xbjt: 1.0,
            xdif: 1.0,
            xrec: 1.0,
            xtun: 0.0,
            ahli: 0.0,
            lbjt0: 0.2e-6,
            ln: 2e-6,
            nbjt: 1.0,
            ndif: -1.0,
            aely: 0.0,
            vabjt: 10.0,
            agidl: 0.0,
            bgidl: 2.3e9,
            ngidl: 1.2,
            kb1: 1.0,
            kb3: 1.0,
            dvbd0: 0.0,
            dvbd1: 0.0,
            vbsa: 0.0,
            delp: 0.02,
            abp: 1.0,
            mxc: -0.9,
            adice0: 1.0,
            kbjt1: 0.0,
            edl: 0.0,
            rbody: 0.0,
            rbsh: 0.0,
            cgso: 0.0,
            cgdo: 0.0,
            clc: 1e-7,
            cle: 0.6,
            cf: 0.0,
            ckappa: 0.6,
            cgdl: 0.0,
            cgsl: 0.0,
            cjswg: 0.0,
            mjswg: 0.5,
            pbswg: 1.0,
            tt: 0.0,
            csdesw: 0.0,
            asd: 0.3,
            rth0: 0.0,
            cth0: 0.0,
            cox: 0.0,
            vtm: 0.0,
            phi: 0.0,
            sqrt_phi: 0.0,
            vbi_default: 0.0,
            factor1: 0.0,
            ni: 0.0,
            eg: 0.0,
            cbox: 0.0,
            csi: 0.0,
            qsi: 0.0,
            csieff: 0.0,
            qsieff: 0.0,
            vfbb: 0.0,
        };
        m.precompute();
        m
    }

    pub fn from_model_def(def: &ModelDef) -> Self {
        let mos_type = match def.kind.to_uppercase().as_str() {
            "PMOS" => MosfetType::Pmos,
            _ => MosfetType::Nmos,
        };
        let mut m = Self::new(mos_type);

        fn pf(def: &ModelDef, name: &str) -> Option<f64> {
            def.params.iter().find_map(|p| {
                if p.name.eq_ignore_ascii_case(name) {
                    if let Expr::Num(v) = &p.value {
                        Some(*v)
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
        }

        macro_rules! set {
            ($field:ident, $name:expr) => {
                if let Some(v) = pf(def, $name) {
                    m.$field = v;
                }
            };
        }
        macro_rules! seti {
            ($field:ident, $name:expr) => {
                if let Some(v) = pf(def, $name) {
                    m.$field = v as i32;
                }
            };
        }

        seti!(mob_mod, "MOBMOD");
        seti!(cap_mod, "CAPMOD");
        seti!(sh_mod, "SHMOD");
        set!(tox, "TOX");
        set!(tsi, "TSI");
        set!(tbox, "TBOX");
        set!(vth0, "VTH0");
        set!(k1, "K1");
        set!(k2, "K2");
        set!(k3, "K3");
        set!(k3b, "K3B");
        set!(w0, "W0");
        set!(nlx, "NLX");
        set!(dvt0, "DVT0");
        set!(dvt1, "DVT1");
        set!(dvt2, "DVT2");
        set!(dvt0w, "DVT0W");
        set!(dvt1w, "DVT1W");
        set!(dvt2w, "DVT2W");
        set!(dsub, "DSUB");
        set!(eta0, "ETA0");
        set!(etab, "ETAB");
        set!(voff, "VOFF");
        set!(nfactor, "NFACTOR");
        set!(cdsc, "CDSC");
        set!(cdscb, "CDSCB");
        set!(cdscd, "CDSCD");
        set!(cit, "CIT");
        set!(nch, "NCH");
        set!(npeak, "NPEAK");
        set!(ngate, "NGATE");
        set!(nsub, "NSUB");
        set!(u0, "U0");
        set!(ua, "UA");
        set!(ub, "UB");
        set!(uc, "UC");
        set!(ute, "UTE");
        set!(ua1, "UA1");
        set!(ub1, "UB1");
        set!(uc1, "UC1");
        set!(vsat, "VSAT");
        set!(a0, "A0");
        set!(ags, "AGS");
        set!(a1, "A1");
        set!(a2, "A2");
        set!(at, "AT");
        set!(keta, "KETA");
        set!(pclm, "PCLM");
        set!(pdiblc1, "PDIBLC1");
        set!(pdiblc2, "PDIBLC2");
        set!(pdiblcb, "PDIBLCB");
        set!(drout, "DROUT");
        set!(pvag, "PVAG");
        set!(delta, "DELTA");
        set!(rdsw, "RDSW");
        set!(prwg, "PRWG");
        set!(prwb, "PRWB");
        set!(prt, "PRT");
        set!(wr, "WR");
        set!(dwg, "DWG");
        set!(dwb, "DWB");
        set!(b0, "B0");
        set!(b1, "B1");
        set!(wint, "WINT");
        set!(lint, "LINT");
        set!(dlc, "DLC");
        set!(dwc, "DWC");
        set!(alpha0, "ALPHA0");
        set!(alpha1, "ALPHA1");
        set!(beta0, "BETA0");
        set!(aii, "AII");
        set!(bii, "BII");
        set!(cii, "CII");
        set!(dii, "DII");
        set!(tnom, "TNOM");
        set!(kt1, "KT1");
        set!(kt1l, "KT1L");
        set!(kt2, "KT2");
        set!(ndiode, "NDIODE");
        set!(ntun, "NTUN");
        set!(nrecf0, "NRECF0");
        set!(nrecr0, "NRECR0");
        set!(vrec0, "VREC0");
        set!(ntrecf, "NTRECF");
        set!(ntrecr, "NTRECR");
        set!(isbjt, "ISBJT");
        set!(isdif, "ISDIF");
        set!(istun, "ISTUN");
        set!(isrec, "ISREC");
        set!(xbjt, "XBJT");
        set!(xdif, "XDIF");
        set!(xrec, "XREC");
        set!(xtun, "XTUN");
        set!(ahli, "AHLI");
        set!(lbjt0, "LBJT0");
        set!(ln, "LN");
        set!(nbjt, "NBJT");
        set!(ndif, "NDIF");
        set!(aely, "AELY");
        set!(vabjt, "VABJT");
        set!(agidl, "AGIDL");
        set!(bgidl, "BGIDL");
        set!(ngidl, "NGIDL");
        set!(kb1, "KB1");
        set!(kb3, "KB3");
        set!(dvbd0, "DVBD0");
        set!(dvbd1, "DVBD1");
        set!(vbsa, "VBSA");
        set!(delp, "DELP");
        set!(abp, "ABP");
        set!(mxc, "MXC");
        set!(adice0, "ADICE0");
        set!(kbjt1, "KBJT1");
        set!(edl, "EDL");
        set!(rbody, "RBODY");
        set!(rbsh, "RBSH");
        set!(cgso, "CGSO");
        set!(cgdo, "CGDO");
        set!(clc, "CLC");
        set!(cle, "CLE");
        set!(cf, "CF");
        set!(ckappa, "CKAPPA");
        set!(cgdl, "CGDL");
        set!(cgsl, "CGSL");
        set!(cjswg, "CJSWG");
        set!(mjswg, "MJSWG");
        set!(pbswg, "PBSWG");
        set!(tt, "TT");
        set!(csdesw, "CSDESW");
        set!(asd, "ASD");
        set!(rth0, "RTH0");
        set!(cth0, "CTH0");

        // Handle u0 units: ngspice treats u0 > 1 as cm²/Vs, converts by /1e4
        if m.u0 > 1.0 {
            m.u0 /= 1e4;
        }

        m.precompute();
        m
    }

    fn precompute(&mut self) {
        let tnom_k = self.tnom + 273.15;
        self.cox = EPSOX / self.tox;
        self.vtm = KBOQ * TEMP_DEFAULT;

        let eg0 = 1.16 - 7.02e-4 * tnom_k * tnom_k / (tnom_k + 1108.0);
        self.eg = eg0;
        self.ni = 1.45e10
            * (tnom_k / 300.15)
            * (tnom_k / 300.15).sqrt()
            * (21.5565981 - eg0 / (2.0 * KBOQ * tnom_k)).exp();

        let npeak = if self.npeak > 1e20 {
            self.npeak * 1e-6
        } else {
            self.npeak
        };

        self.phi = 2.0 * self.vtm * (npeak / self.ni).ln();
        if self.phi < 0.4 {
            self.phi = 0.4;
        }
        self.sqrt_phi = self.phi.sqrt();

        self.vbi_default = self.vtm * (npeak * 1e10 / (self.ni * self.ni)).ln();
        self.factor1 = (EPSSI / (EPSOX * self.tox)).sqrt();

        // DD-specific precomputed (same as FD)
        self.cbox = EPSOX / self.tbox;
        self.csi = EPSSI / self.tsi;
        let nsub = if self.nsub > 1e20 {
            self.nsub * 1e-6
        } else {
            self.nsub
        };
        self.qsi = CHARGE_Q * npeak * self.tsi;
        // csieff and qsieff use surface potential-corrected values
        self.csieff = self.csi;
        self.qsieff = self.qsi.max(1e-20);
        // Back-gate flat-band voltage
        self.vfbb = self.vtm * (npeak * 1e6 / (nsub * 1e6)).ln();
    }

    /// Number of internal nodes this model creates.
    pub fn internal_node_count(&self, nrd: f64, nrs: f64) -> usize {
        let mut count = 1; // Always create internal body node
        if self.rdsw > 0.0 || nrd > 0.0 {
            count += 1; // drain prime
        }
        if self.rdsw > 0.0 || nrs > 0.0 {
            count += 1; // source prime
        }
        count
    }

    pub fn size_dep_param(&self, w: f64, l: f64, temp: f64) -> Bsim3SoiDdSizeParam {
        let tnom_k = self.tnom + 273.15;
        let vtm = KBOQ * temp;

        let dl = self.lint;
        let dw = self.wint;
        let leff = l - 2.0 * dl;
        let weff = w - 2.0 * dw;
        let dlc = if self.dlc != 0.0 { self.dlc } else { dl };
        let dwc = if self.dwc != 0.0 { self.dwc } else { dw };
        let leff_cv = l - 2.0 * dlc;
        let weff_cv = w - 2.0 * dwc;

        let phi = self.phi;
        let sqrt_phi = self.sqrt_phi;
        let temp_ratio = temp / tnom_k;

        let u0temp = self.u0 * temp_ratio.powf(self.ute);
        let vsattemp = self.vsat - self.at * (temp - tnom_k);

        let rds0 = if self.rdsw > 0.0 {
            (self.rdsw + self.prt * (temp - tnom_k)) / weff.powf(self.wr) * 1e-6
        } else {
            0.0
        };

        let ua = self.ua + self.ua1 * (temp - tnom_k);
        let ub = self.ub + self.ub1 * (temp - tnom_k);
        let uc = self.uc + self.uc1 * (temp - tnom_k);

        let npeak = if self.npeak > 1e20 {
            self.npeak * 1e-6
        } else {
            self.npeak
        };
        let nsub = if self.nsub > 1e20 {
            self.nsub * 1e-6
        } else {
            self.nsub
        };

        let xdep0 = (2.0 * EPSSI / (CHARGE_Q * npeak * 1e6)).sqrt() * sqrt_phi;
        let xj = 1.5e-7; // Default XJ
        let litl = (EPSSI * xj / self.cox).sqrt();

        let t0 = -0.5 * self.dvt1 * leff / litl;
        let theta0 = if t0 > -EXP_THRESHOLD {
            let t1 = t0.exp();
            t1 * (1.0 + 2.0 * t1)
        } else {
            MIN_EXP * (1.0 + 2.0 * MIN_EXP)
        };
        let theta0vb0 = self.dvt0 * theta0;

        let t0 = -0.5 * self.drout * leff / litl;
        let theta_rout = if t0 > -EXP_THRESHOLD {
            let t1 = t0.exp();
            self.pdiblc1 * (t1 * (1.0 + 2.0 * t1)) + self.pdiblc2
        } else {
            self.pdiblc2
        };

        let k1eff = self.k1;
        let vfb = -1.0;
        let cdep0 = (CHARGE_Q * EPSSI * npeak * 1e6).sqrt() / (2.0 * phi);

        let eg = 1.16 - 7.02e-4 * temp * temp / (temp + 1108.0);
        let ni_temp = 1.45e10
            * (temp / 300.15)
            * (temp / 300.15).sqrt()
            * (21.5565981 - eg / (2.0 * vtm)).exp();
        let vbi = vtm * (npeak * 1e6 / (ni_temp * ni_temp)).abs().ln();

        // SOI junction parameters
        let t0 = (eg / (2.0 * KBOQ * tnom_k)).exp();
        let t_eg = (eg / (2.0 * KBOQ * temp)).exp();
        let t0_ratio = t0 / t_eg;

        let jbjt = self.isbjt * t0_ratio;
        let jdif = self.isdif * t0_ratio;
        let jrec = self.isrec
            * ((self.nrecf0 * 0.026 * (1.0 + self.ntrecf * (temp_ratio - 1.0)))
                / (self.nrecf0 * 0.026))
                .exp();
        let jtun = self.istun;

        let wdios = weff * self.asd;
        let wdiod = weff * self.asd;

        let lratio = if self.lbjt0 > 0.0 {
            (1.0 - leff / (leff + self.lbjt0)) / (1.0 + (self.ndif * leff / (leff + self.lbjt0)))
        } else {
            0.0
        };
        let lratiodif = lratio;
        let vearly = if self.vabjt > 0.0 { self.vabjt } else { 10.0 };
        let arfabjt = self.xbjt;

        let cgso_eff = if self.cgso > 0.0 {
            self.cgso
        } else {
            0.6 * dlc * self.cox
        };
        let cgdo_eff = if self.cgdo > 0.0 {
            self.cgdo
        } else {
            0.6 * dlc * self.cox
        };

        Bsim3SoiDdSizeParam {
            leff,
            weff,
            leff_cv,
            weff_cv,
            vth0: self.vth0,
            k1: self.k1,
            k2: self.k2,
            k3: self.k3,
            k3b: self.k3b,
            w0: self.w0,
            nlx: self.nlx,
            dvt0: self.dvt0,
            dvt1: self.dvt1,
            dvt2: self.dvt2,
            dvt0w: self.dvt0w,
            dvt1w: self.dvt1w,
            dvt2w: self.dvt2w,
            eta0: self.eta0,
            etab: self.etab,
            dsub: self.dsub,
            voff: self.voff,
            nfactor: self.nfactor,
            cdsc: self.cdsc,
            cdscb: self.cdscb,
            cdscd: self.cdscd,
            cit: self.cit,
            u0: u0temp,
            ua,
            ub,
            uc,
            vsat: vsattemp,
            a0: self.a0,
            ags: self.ags,
            a1: self.a1,
            a2: self.a2,
            at: self.at,
            keta: self.keta,
            pclm: self.pclm,
            pdiblc1: self.pdiblc1,
            pdiblc2: self.pdiblc2,
            pdiblcb: self.pdiblcb,
            pvag: self.pvag,
            delta: self.delta,
            rdsw: self.rdsw,
            prwg: self.prwg,
            prwb: self.prwb,
            dwg: self.dwg,
            dwb: self.dwb,
            b0: self.b0,
            b1: self.b1,
            alpha0: self.alpha0,
            beta0: self.beta0,
            kt1: self.kt1,
            kt1l: self.kt1l,
            kt2: self.kt2,
            ndiode: self.ndiode,
            agidl: self.agidl,
            bgidl: self.bgidl,
            ngidl: self.ngidl,
            phi,
            sqrt_phi,
            xdep0,
            litl,
            theta0vb0,
            theta_rout,
            k1eff,
            npeak,
            nsub,
            vfb,
            u0temp,
            vsattemp,
            rds0,
            cdep0,
            vbi,
            lratio,
            lratiodif,
            vearly,
            arfabjt,
            wdios,
            wdiod,
            jbjt,
            jdif,
            jrec,
            jtun,
            cgso_eff,
            cgdo_eff,
            dvbd0: self.dvbd0,
            dvbd1: self.dvbd1,
            nseg: 1.0,
        }
    }
}

impl Bsim3SoiDdInstance {
    pub fn terminal_voltages(&self, solution: &[f64]) -> (f64, f64, f64, f64) {
        let vg = self.gate_idx.map_or(0.0, |i| solution[i]);
        let vd = self.drain_eff_idx().map_or(0.0, |i| solution[i]);
        let vs = self.source_eff_idx().map_or(0.0, |i| solution[i]);
        let vb = self.body_int_idx.map_or(0.0, |i| solution[i]);

        let sign = self.model.mos_type.sign();
        let vgs = sign * (vg - vs);
        let vds = sign * (vd - vs);
        let vbs = sign * (vb - vs);
        let ves = sign * (self.e_idx.map_or(0.0, |i| solution[i]) - vs);

        (vgs, vds, vbs, ves)
    }

    pub fn drain_eff_idx(&self) -> Option<usize> {
        self.drain_prime_idx.or(self.drain_idx)
    }

    pub fn source_eff_idx(&self) -> Option<usize> {
        self.source_prime_idx.or(self.source_idx)
    }

    pub fn ac_stamp(&self, comp: &Bsim3SoiDdCompanion) -> crate::ac::BsimAcStamp {
        crate::ac::BsimAcStamp {
            dp: self.drain_eff_idx(),
            g: self.gate_idx,
            sp: self.source_eff_idx(),
            b: self.body_int_idx,
            drain_idx: self.drain_idx,
            source_idx: self.source_idx,
            m: self.m,
            gm: comp.gm,
            gds: comp.gds,
            gmbs: comp.gmbs,
            gbd: comp.gbd_jct,
            gbs: comp.gbs_jct,
            g_drain: 0.0,
            g_source: 0.0,
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

/// Compute BSIM3SOI-DD companion model (NR linearization).
///
/// DD uses the FD-style self-consistent surface potential chain for body bias,
/// combined with PD-style 4-component junction currents and GIDL.
#[expect(clippy::too_many_lines)]
pub fn bsim3soi_dd_companion(
    vgs: f64,
    vds: f64,
    vbs: f64,
    ves: f64,
    sp: &Bsim3SoiDdSizeParam,
    model: &Bsim3SoiDdModel,
) -> Bsim3SoiDdCompanion {
    let sign = model.mos_type.sign();
    let cox = model.cox;
    let vtm = KBOQ * TEMP_DEFAULT;
    let gmin = 1e-12;
    let phi = sp.phi;
    let sqrt_phi = sp.sqrt_phi;

    // Mode detection (forward/reverse)
    let (vgs_i, vds_i, vbs_i, ves_i, mode) = if vds >= 0.0 {
        (vgs, vds, vbs, ves, 1)
    } else {
        (vgs - vds, -vds, vbs - vds, ves - vds, -1)
    };

    let leff = sp.leff;
    let weff = sp.weff;
    let vbi = sp.vbi;
    let v0 = vbi - phi;

    let vesfb = ves_i - model.vfbb;

    // ========== DD surface potential chain (same as FD) ==========

    // Vbs0t
    let t0 = -sp.dvbd1 * leff / sp.litl;
    let t1 = sp.dvbd0 * (safe_exp(0.5 * t0) + 2.0 * safe_exp(t0));
    let t2 = t1 * v0;
    let t3 = 0.5 * model.qsi / model.csi;
    let vbs0t = phi - t3 + model.vbsa + t2;

    // Vbs0 (with back-gate coupling)
    let t0_kb = 1.0 + model.csieff / model.cbox;
    let t1_kb = model.kb1 / t0_kb;
    let t2_kb = t1_kb * (vbs0t - vesfb);
    let t6_vbs0 = vbs0t - t2_kb;

    // Limit Vbs0 below phi - delp
    let t1_lim = phi - model.delp;
    let t2_lim = t1_lim - t6_vbs0 - DELT_VBSEFF;
    let t3_lim = (t2_lim * t2_lim + 4.0 * DELT_VBSEFF).sqrt();
    let vbs0 = t1_lim - 0.5 * (t2_lim + t3_lim);

    // Vbs0mos
    let t1_mos = vbs0t - vbs0 - DELT_VBSMOS;
    let t2_mos = (t1_mos * t1_mos + DELT_VBSMOS * DELT_VBSMOS).sqrt();
    let t3_mos = 0.5 * (t1_mos + t2_mos);
    let t4_mos = t3_mos * model.csieff / model.qsieff;
    let vbs0mos = vbs0 - 0.5 * t3_mos * t4_mos;

    // Vthfd (threshold voltage using Vbs0mos)
    let phis_fd = phi - vbs0mos;
    let sqrt_phis_fd = phis_fd.abs().sqrt();
    let xdep_fd = sp.xdep0 * sqrt_phis_fd / sqrt_phi;

    let t0_dvt = sp.dvt2 * vbs0mos;
    let t1_dvt = if t0_dvt >= -0.5 {
        1.0 + t0_dvt
    } else {
        let t4 = 1.0 / (3.0 + 8.0 * t0_dvt);
        (1.0 + 3.0 * t0_dvt) * t4
    };
    let lt1_fd = model.factor1 * xdep_fd.sqrt() * t1_dvt;

    let t0_sce = -0.5 * sp.dvt1 * leff / lt1_fd;
    let theta0_fd = if t0_sce > -EXP_THRESHOLD {
        let t1 = t0_sce.exp();
        t1 * (1.0 + 2.0 * t1)
    } else {
        MIN_EXP * (1.0 + 2.0 * MIN_EXP)
    };
    let delt_vth_fd = sp.dvt0 * theta0_fd * v0;

    let t0_dvtw = sp.dvt2w * vbs0mos;
    let t1_dvtw = if t0_dvtw >= -0.5 {
        1.0 + t0_dvtw
    } else {
        let t4 = 1.0 / (3.0 + 8.0 * t0_dvtw);
        (1.0 + 3.0 * t0_dvtw) * t4
    };
    let ltw_fd = model.factor1 * xdep_fd.sqrt() * t1_dvtw;
    let t0_w_fd = -0.5 * sp.dvt1w * weff * leff / ltw_fd;
    let t2_w_fd = if t0_w_fd > -EXP_THRESHOLD {
        let t1 = t0_w_fd.exp();
        t1 * (1.0 + 2.0 * t1)
    } else {
        MIN_EXP * (1.0 + 2.0 * MIN_EXP)
    };
    let delt_vthw_fd = sp.dvt0w * t2_w_fd * v0;

    let temp_ratio_minus1 = (TEMP_DEFAULT / (model.tnom + 273.15)) - 1.0;
    let t0_nlx = (1.0 + sp.nlx / leff).sqrt();
    let t1_kt = sp.kt1 + sp.kt1l / leff + sp.kt2 * vbs0mos;
    let delt_vth_temp_fd = sp.k1 * (t0_nlx - 1.0) * sqrt_phi + t1_kt * temp_ratio_minus1;

    let tmp2_fd = model.tox * phi / (weff + sp.w0);
    let t3_eta_fd = sp.eta0 + sp.etab * vbs0mos;
    let t3_eta_fd_eff = if t3_eta_fd < 1e-4 {
        let t9 = 1.0 / (3.0 - 2e4 * t3_eta_fd);
        (2e-4 - t3_eta_fd) * t9
    } else {
        t3_eta_fd
    };
    let dibl_sft_fd = t3_eta_fd_eff * sp.theta0vb0 * vds_i;

    let vthfd = sign * sp.vth0 + sp.k1 * (sqrt_phis_fd - sqrt_phi)
        - sp.k2 * vbs0mos
        - delt_vth_fd
        - delt_vthw_fd
        + (sp.k3 + sp.k3b * vbs0mos) * tmp2_fd
        + delt_vth_temp_fd
        - dibl_sft_fd;

    // Poly gate depletion
    let t0_poly = sp.vfb + phi;
    let (vgs_eff, dvgs_eff_dvg) = if model.ngate > 1e18 && model.ngate < 1e25 && vgs_i > t0_poly {
        let t1 = 1e6 * CHARGE_Q * EPSSI * model.ngate / (cox * cox);
        let t4 = (1.0 + 2.0 * (vgs_i - t0_poly) / t1).sqrt();
        let t2 = t1 * (t4 - 1.0);
        let t3 = 0.5 * t2 * t2 / t1;
        let t7 = 1.12 - t3 - 0.05;
        let t6 = (t7 * t7 + 0.224).sqrt();
        let t5 = 1.12 - 0.5 * (t7 + t6);
        (vgs_i - t5, 1.0 - (0.5 - 0.5 / t4) * (1.0 + t7 / t6))
    } else {
        (vgs_i, 1.0)
    };

    // ========== DD Vbs0eff and Vbsmos calculation ==========

    let t1_eff = vthfd - vgs_eff - DELT_VBS0EFF;
    let t2_eff = (t1_eff * t1_eff + DELT_VBS0EFF * DELT_VBS0EFF).sqrt();

    let vbs0teff = vbs0t - 0.5 * (t1_eff + t2_eff);

    // Nfb (feedback factor)
    let k1 = sp.k1;
    let t3_nfb = 1.0 / (k1 * k1);
    let t4_nfb = model.kb3 * model.cbox / cox;
    let t8_nfb = (phi - vbs0mos).abs().sqrt();
    let t5_nfb = (1.0 + 4.0 * t3_nfb * (phi + k1 * t8_nfb - vbs0mos))
        .abs()
        .sqrt();
    let t6_nfb = 1.0 + t4_nfb * t5_nfb;
    let nfb = 1.0 / t6_nfb;

    let vbs0eff_dd = vbs0 - nfb * 0.5 * (t1_eff + t2_eff);
    let vbsdio = vbs0eff_dd;

    // Vbsmos
    let t1_bsmos = vbs0teff - vbsdio - DELT_VBSMOS;
    let t2_bsmos = (t1_bsmos * t1_bsmos + DELT_VBSMOS * DELT_VBSMOS).sqrt();
    let t3_bsmos = 0.5 * (t1_bsmos + t2_bsmos);
    let t4_bsmos = t3_bsmos * model.csieff / model.qsieff;
    let vbsmos = vbsdio - 0.5 * t3_bsmos * t4_bsmos;

    // ========== Vbseff (final body-source effective voltage) ==========
    let t1_vbseff = phi - model.delp;
    let t2_vbseff = t1_vbseff - vbsmos - DELT_VBSEFF;
    let t3_vbseff = (t2_vbseff * t2_vbseff + 4.0 * DELT_VBSEFF * t1_vbseff).sqrt();
    let vbseff = t1_vbseff - 0.5 * (t2_vbseff + t3_vbseff);

    // The DD-computed Vbs for body node feedback
    let vbs_dd = vbsdio;

    // ========== Main MOSFET equations ==========
    let phis = phi - vbseff;
    let sqrt_phis = phis.abs().sqrt();
    let dsqrt_phis_dvb = -0.5 / sqrt_phis.max(1e-20);
    let xdep = sp.xdep0 * sqrt_phis / sqrt_phi;
    let dxdep_dvb = (sp.xdep0 / sqrt_phi) * dsqrt_phis_dvb;

    // Vth calculation
    let t3_vth = xdep.sqrt();

    let t0 = sp.dvt2 * vbseff;
    let (t1, t2_) = if t0 >= -0.5 {
        (1.0 + t0, sp.dvt2)
    } else {
        let t4 = 1.0 / (3.0 + 8.0 * t0);
        ((1.0 + 3.0 * t0) * t4, sp.dvt2 * t4 * t4)
    };
    let lt1 = model.factor1 * t3_vth * t1;
    let dlt1_dvb = model.factor1 * (0.5 / t3_vth * t1 * dxdep_dvb + t3_vth * t2_);

    let t0w = sp.dvt2w * vbseff;
    let (t1w, t2w) = if t0w >= -0.5 {
        (1.0 + t0w, sp.dvt2w)
    } else {
        let t4 = 1.0 / (3.0 + 8.0 * t0w);
        ((1.0 + 3.0 * t0w) * t4, sp.dvt2w * t4 * t4)
    };
    let ltw = model.factor1 * t3_vth * t1w;
    let dltw_dvb = model.factor1 * (0.5 / t3_vth * t1w * dxdep_dvb + t3_vth * t2w);

    let t0_sce2 = -0.5 * sp.dvt1 * leff / lt1;
    let (theta0, dtheta0_dvb) = if t0_sce2 > -EXP_THRESHOLD {
        let t1 = t0_sce2.exp();
        (
            t1 * (1.0 + 2.0 * t1),
            (-t0_sce2 / lt1 * t1 * dlt1_dvb) * (1.0 + 4.0 * t1),
        )
    } else {
        (MIN_EXP * (1.0 + 2.0 * MIN_EXP), 0.0)
    };
    let delt_vth = sp.dvt0 * theta0 * v0;
    let ddelt_vth_dvb = sp.dvt0 * dtheta0_dvb * v0;

    let t0_w = -0.5 * sp.dvt1w * weff * leff / ltw;
    let (t2_w, dt2w_dvb) = if t0_w > -EXP_THRESHOLD {
        let t1 = t0_w.exp();
        (
            t1 * (1.0 + 2.0 * t1),
            (-t0_w / ltw * t1 * dltw_dvb) * (1.0 + 4.0 * t1),
        )
    } else {
        (MIN_EXP * (1.0 + 2.0 * MIN_EXP), 0.0)
    };
    let delt_vthw = sp.dvt0w * t2_w * v0;
    let ddelt_vthw_dvb = sp.dvt0w * dt2w_dvb * v0;

    let delt_vth_temp = sp.k1 * (t0_nlx - 1.0) * sqrt_phi + t1_kt * temp_ratio_minus1;

    let tmp2 = model.tox * phi / (weff + sp.w0);
    let t3_eta = sp.eta0 + sp.etab * vbseff;
    let (t3_eta_eff, dt3_dvb) = if t3_eta < 1e-4 {
        let t9 = 1.0 / (3.0 - 2e4 * t3_eta);
        ((2e-4 - t3_eta) * t9, t9 * t9 * sp.etab)
    } else {
        (t3_eta, sp.etab)
    };
    let dibl_sft = t3_eta_eff * sp.theta0vb0 * vds_i;
    let ddibl_sft_dvd = sp.theta0vb0 * t3_eta_eff;
    let ddibl_sft_dvb = sp.theta0vb0 * vds_i * dt3_dvb;

    let vth =
        sign * sp.vth0 + sp.k1 * (sqrt_phis - sqrt_phi) - sp.k2 * vbseff - delt_vth - delt_vthw
            + (sp.k3 + sp.k3b * vbseff) * tmp2
            + delt_vth_temp
            - dibl_sft;

    let t6 = sp.k3b * tmp2 - sp.k2 + sp.kt2 * temp_ratio_minus1;
    let dvth_dvb = sp.k1 * dsqrt_phis_dvb - ddelt_vth_dvb - ddelt_vthw_dvb + t6 - ddibl_sft_dvb;
    let dvth_dvd = -ddibl_sft_dvd;

    // Calculate n (subthreshold swing factor)
    let t2_n = sp.nfactor * EPSSI / xdep;
    let dt2n_dvb = -t2_n / xdep * dxdep_dvb;
    let t3_n = sp.cdsc + sp.cdscb * vbseff + sp.cdscd * vds_i;
    let dt3n_dvb = sp.cdscb;
    let dt3n_dvd = sp.cdscd;
    let t4_n = (t2_n + t3_n * theta0 + sp.cit) / cox;
    let dt4n_dvb = (dt2n_dvb + theta0 * dt3n_dvb + dtheta0_dvb * t3_n) / cox;
    let dt4n_dvd = theta0 * dt3n_dvd / cox;
    let (n, dn_dvb, dn_dvd) = if t4_n >= -0.5 {
        (1.0 + t4_n, dt4n_dvb, dt4n_dvd)
    } else {
        let t0 = 1.0 / (3.0 + 8.0 * t4_n);
        let n = (1.0 + 3.0 * t4_n) * t0;
        let t0sq = t0 * t0;
        (n, t0sq * dt4n_dvb, t0sq * dt4n_dvd)
    };

    // Vgsteff (effective gate overdrive)
    let vgst = vgs_eff - vth;
    let dvgst_dvg = dvgs_eff_dvg;
    let dvgst_dvd = -dvth_dvd;
    let dvgst_dvb = -dvth_dvb;

    let t10 = 2.0 * n * vtm;
    let vgst_nvt = vgst / t10;
    let exp_arg = (2.0 * sp.voff - vgst) / t10;

    // Use dvbseff_dvb = 1 for DD (body voltage is self-consistent)
    let dvbseff_dvb = 1.0;

    let (vgsteff, dvgsteff_dvg, dvgsteff_dvd, dvgsteff_dvb) = if vgst_nvt > EXPL_THRESHOLD {
        (vgst, dvgs_eff_dvg, -dvth_dvd, -dvth_dvb * dvbseff_dvb)
    } else if exp_arg > EXPL_THRESHOLD {
        let t0 = (vgst - sp.voff) / (n * vtm);
        let exp_vgst = t0.exp();
        let vgsteff_val = vtm * sp.cdep0 / cox * exp_vgst;
        let t3 = vgsteff_val / (n * vtm);
        let t1 = -t3 * (dvth_dvb + t0 * vtm * dn_dvb);
        (
            vgsteff_val,
            t3 * dvgs_eff_dvg,
            -t3 * (dvth_dvd + t0 * vtm * dn_dvd),
            t1 * dvbseff_dvb,
        )
    } else {
        let exp_vgst = vgst_nvt.exp();
        let t1 = t10 * (1.0 + exp_vgst).ln();
        let dt1_dvg = exp_vgst / (1.0 + exp_vgst);
        let dt1_dvb = -dt1_dvg * (dvth_dvb + vgst / n * dn_dvb) + t1 / n * dn_dvb;
        let dt1_dvd = -dt1_dvg * (dvth_dvd + vgst / n * dn_dvd) + t1 / n * dn_dvd;

        let dt2_dvg = -cox / (vtm * sp.cdep0) * exp_arg.exp();
        let t2_val = 1.0 - t10 * dt2_dvg;
        let dt2_dvd =
            -dt2_dvg * (dvth_dvd - 2.0 * vtm * exp_arg * dn_dvd) + (t2_val - 1.0) / n * dn_dvd;
        let dt2_dvb =
            -dt2_dvg * (dvth_dvb - 2.0 * vtm * exp_arg * dn_dvb) + (t2_val - 1.0) / n * dn_dvb;

        let vgsteff_val = t1 / t2_val;
        let t3 = t2_val * t2_val;
        let t4 = (t2_val * dt1_dvb - t1 * dt2_dvb) / t3;
        (
            vgsteff_val,
            (t2_val * dt1_dvg - t1 * dt2_dvg) / t3 * dvgs_eff_dvg,
            (t2_val * dt1_dvd - t1 * dt2_dvd) / t3,
            t4 * dvbseff_dvb,
        )
    };

    let vgst2vtm = vgsteff + 2.0 * vtm;

    // Effective channel geometry
    let mut weff_ch = weff - 2.0 * (sp.dwg * vgsteff + sp.dwb * (sqrt_phis - sqrt_phi));
    let mut dweff_dvg = -2.0 * sp.dwg;
    let mut dweff_dvb = -2.0 * sp.dwb * dsqrt_phis_dvb;
    if weff_ch < 2e-8 {
        let t0 = 1.0 / (6e-8 - 2.0 * weff_ch);
        weff_ch = 2e-8 * (4e-8 - weff_ch) * t0;
        let t0sq = t0 * t0 * 4e-16;
        dweff_dvg *= t0sq;
        dweff_dvb *= t0sq;
    }

    // Series resistance Rds
    let t0_rds = sp.prwg * vgsteff + sp.prwb * (sqrt_phis - sqrt_phi);
    let (rds, drds_dvg, drds_dvb) = if t0_rds >= -0.9 {
        (
            sp.rds0 * (1.0 + t0_rds),
            sp.rds0 * sp.prwg,
            sp.rds0 * sp.prwb * dsqrt_phis_dvb,
        )
    } else {
        let t1 = 1.0 / (17.0 + 20.0 * t0_rds);
        (
            sp.rds0 * (0.8 + t0_rds) * t1,
            sp.rds0 * sp.prwg * t1 * t1,
            sp.rds0 * sp.prwb * dsqrt_phis_dvb * t1 * t1,
        )
    };

    // Abulk calculation
    let (abulk0, dabulk0_dvb, abulk, dabulk_dvg, dabulk_dvb) = if sp.a0 == 0.0 {
        (1.0, 0.0, 1.0, 0.0, 0.0)
    } else {
        let t10_k = sp.keta * vbseff;
        let (t11, dt11_dvb) = if t10_k >= -0.9 {
            let t11 = 1.0 / (1.0 + t10_k);
            (t11, -sp.keta * t11 * t11)
        } else {
            let t12 = 1.0 / (0.8 + t10_k);
            let t11 = (17.0 + 20.0 * t10_k) * t12;
            (t11, -sp.keta * t12 * t12)
        };

        let t10_phi = phi;
        let t13 = (vbseff * t11) / t10_phi;
        let dt13_dvb = (vbseff * dt11_dvb + t11) / t10_phi;

        let (t14, dt14_dvb) = if t13 < 0.96 {
            let t14 = 1.0 / (1.0 - t13).sqrt();
            let t10 = 0.5 * t14 / (1.0 - t13);
            (t14, t10 * dt13_dvb)
        } else {
            let t11 = 1.0 / (1.0 - 1.043406 * t13);
            let t14 = (6.00167 - 6.26044 * t13) * t11;
            let t10 = 0.001742 * t11 * t11;
            (t14, t10 * dt13_dvb)
        };

        let t10_k1 = 0.5 * sp.k1eff / phi.sqrt();
        let t1 = t10_k1 * t14;
        let dt1_dvb = t10_k1 * dt14_dvb;

        let t9 = (1.5e-7 * xdep).sqrt(); // XJ=1.5e-7
        let tmp1 = leff + 2.0 * t9;
        let t5 = leff / tmp1;
        let tmp2_a = sp.a0 * t5;
        let tmp3_a = weff + sp.b1;
        let tmp4_a = sp.b0 / tmp3_a;
        let t2 = tmp2_a + tmp4_a;
        let dt2_dvb = -t9 * tmp2_a / tmp1 / xdep * dxdep_dvb;
        let t6 = t5 * t5;
        let t7 = t5 * t6;

        let mut abulk0 = 1.0 + t1 * t2;
        let mut dabulk0_dvb = t1 * dt2_dvb + t2 * dt1_dvb;

        let t8 = sp.ags * sp.a0 * t7;
        let dabulk_dvg = -t1 * t8;
        let mut abulk = abulk0 + dabulk_dvg * vgsteff;
        let mut dabulk_dvb = dabulk0_dvb - t8 * vgsteff * (dt1_dvb + 3.0 * t1 * dt2_dvb / tmp2_a);

        if abulk0 < 0.01 {
            let t9 = 1.0 / (3.0 - 200.0 * abulk0);
            abulk0 = (0.02 - abulk0) * t9;
            dabulk0_dvb *= t9 * t9;
        }
        if abulk < 0.01 {
            let t9 = 1.0 / (3.0 - 200.0 * abulk);
            abulk = (0.02 - abulk) * t9;
            dabulk_dvb *= t9 * t9;
        }

        (abulk0, dabulk0_dvb, abulk, dabulk_dvg, dabulk_dvb)
    };

    // Mobility
    let t5 = if model.mob_mod == 1 {
        let t0 = vgsteff + vth + vth;
        let t2 = sp.ua + sp.uc * vbseff;
        let t3 = t0 / model.tox;
        t3 * (t2 + sp.ub * t3)
    } else if model.mob_mod == 2 {
        vgsteff / model.tox * (sp.ua + sp.uc * vbseff + sp.ub * vgsteff / model.tox)
    } else {
        let t0 = vgsteff + vth + vth;
        let t2 = 1.0 + sp.uc * vbseff;
        let t3 = t0 / model.tox;
        t3 * (sp.ua + sp.ub * t3) * t2
    };

    let denomi = if t5 >= -0.8 {
        1.0 + t5
    } else {
        let t9 = 1.0 / (7.0 + 10.0 * t5);
        (0.6 + t5) * t9
    };

    let ueff = sp.u0temp / denomi;
    let t9 = -ueff / denomi;
    let dueff_dvg = if model.mob_mod == 1 || model.mob_mod == 3 {
        t9 * (sp.ua + 2.0 * sp.ub * (vgsteff + vth + vth) / model.tox) / model.tox
    } else {
        t9 * (sp.ua + sp.uc * vbseff + 2.0 * sp.ub * vgsteff / model.tox) / model.tox
    };
    let dueff_dvd = 0.0;
    let dueff_dvb = 0.0;

    // Saturation voltage Vdsat
    let wvcox = weff_ch * sp.vsattemp * cox;
    let wvcox_rds = wvcox * rds;
    let esat = 2.0 * sp.vsattemp / ueff;
    let esat_l = esat * leff;
    let t0_esat = -esat_l / ueff;
    let desat_l_dvg = t0_esat * dueff_dvg;
    let desat_l_dvd = t0_esat * dueff_dvd;
    let desat_l_dvb = t0_esat * dueff_dvb;

    let a1_val = sp.a1;
    let (lambda, dlambda_dvg) = if a1_val == 0.0 {
        (sp.a2, 0.0)
    } else if a1_val > 0.0 {
        let t0 = 1.0 - sp.a2;
        let t1 = t0 - a1_val * vgsteff - 0.0001;
        let t2 = (t1 * t1 + 0.0004 * t0).sqrt();
        (sp.a2 + t0 - 0.5 * (t1 + t2), 0.5 * a1_val * (1.0 + t1 / t2))
    } else {
        let t1 = sp.a2 + a1_val * vgsteff - 0.0001;
        let t2 = (t1 * t1 + 0.0004 * sp.a2).sqrt();
        (0.5 * (t1 + t2), 0.5 * a1_val * (1.0 + t1 / t2))
    };

    let vdsat;
    let dvdsat_dvg;
    let dvdsat_dvd;
    let dvdsat_dvb;

    let (tmp2_rds, tmp3_rds) = if rds > 0.0 {
        (
            drds_dvg / rds + dweff_dvg / weff_ch,
            drds_dvb / rds + dweff_dvb / weff_ch,
        )
    } else {
        (dweff_dvg / weff_ch, dweff_dvb / weff_ch)
    };

    if rds == 0.0 && lambda == 1.0 {
        let t0 = 1.0 / (abulk * esat_l + vgst2vtm);
        let t1 = t0 * t0;
        let t2 = vgst2vtm * t0;
        let t3 = esat_l * vgst2vtm;
        vdsat = t3 * t0;
        let dt0_dvg = -(abulk * desat_l_dvg + esat_l * dabulk_dvg + 1.0) * t1;
        let dt0_dvd = -(abulk * desat_l_dvd) * t1;
        let dt0_dvb = -(abulk * desat_l_dvb + esat_l * dabulk_dvb) * t1;
        dvdsat_dvg = t3 * dt0_dvg + t2 * desat_l_dvg + esat_l * t0;
        dvdsat_dvd = t3 * dt0_dvd + t2 * desat_l_dvd;
        dvdsat_dvb = t3 * dt0_dvb + t2 * desat_l_dvb;
    } else {
        let t9 = abulk * wvcox_rds;
        let t8 = abulk * t9;
        let t7 = vgst2vtm * t9;
        let t6 = vgst2vtm * wvcox_rds;
        let t0 = 2.0 * abulk * (t9 - 1.0 + 1.0 / lambda);
        let dt0_dvg = 2.0
            * (t8 * tmp2_rds - abulk * dlambda_dvg / (lambda * lambda)
                + (2.0 * t9 + 1.0 / lambda - 1.0) * dabulk_dvg);
        let dt0_dvb =
            2.0 * (t8 * (2.0 / abulk * dabulk_dvb + tmp3_rds) + (1.0 / lambda - 1.0) * dabulk_dvb);
        let dt0_dvd = 0.0;

        let t1 = vgst2vtm * (2.0 / lambda - 1.0) + abulk * esat_l + 3.0 * t7;
        let dt1_dvg = (2.0 / lambda - 1.0) - 2.0 * vgst2vtm * dlambda_dvg / (lambda * lambda)
            + abulk * desat_l_dvg
            + esat_l * dabulk_dvg
            + 3.0 * (t9 + t7 * tmp2_rds + t6 * dabulk_dvg);
        let dt1_dvb =
            abulk * desat_l_dvb + esat_l * dabulk_dvb + 3.0 * (t6 * dabulk_dvb + t7 * tmp3_rds);
        let dt1_dvd = abulk * desat_l_dvd;

        let t2 = vgst2vtm * (esat_l + 2.0 * t6);
        let dt2_dvg = esat_l + vgst2vtm * desat_l_dvg + t6 * (4.0 + 2.0 * vgst2vtm * tmp2_rds);
        let dt2_dvb = vgst2vtm * (desat_l_dvb + 2.0 * t6 * tmp3_rds);
        let dt2_dvd = vgst2vtm * desat_l_dvd;

        let t3 = (t1 * t1 - 2.0 * t0 * t2).sqrt();
        vdsat = (t1 - t3) / t0;
        dvdsat_dvg =
            (dt1_dvg - (t1 * dt1_dvg - dt0_dvg * t2 - t0 * dt2_dvg) / t3 - vdsat * dt0_dvg) / t0;
        dvdsat_dvb =
            (dt1_dvb - (t1 * dt1_dvb - dt0_dvb * t2 - t0 * dt2_dvb) / t3 - vdsat * dt0_dvb) / t0;
        dvdsat_dvd = (dt1_dvd - (t1 * dt1_dvd - t0 * dt2_dvd) / t3) / t0;
    }

    // Vdseff
    let t1 = vdsat - vds_i - sp.delta;
    let t2 = (t1 * t1 + 4.0 * sp.delta * vdsat).sqrt();
    let t0 = t1 / t2;
    let t3 = 2.0 * sp.delta / t2;
    let vdseff = vdsat - 0.5 * (t1 + t2);
    let dvdseff_dvg = dvdsat_dvg - 0.5 * (dvdsat_dvg + t0 * dvdsat_dvg + t3 * dvdsat_dvg);
    let dvdseff_dvd =
        dvdsat_dvd - 0.5 * (dvdsat_dvd - 1.0 + t0 * (dvdsat_dvd - 1.0) + t3 * dvdsat_dvd);
    let dvdseff_dvb = dvdsat_dvb - 0.5 * (dvdsat_dvb + t0 * dvdsat_dvb + t3 * dvdsat_dvb);

    let vdseff = if vdseff > vds_i { vds_i } else { vdseff };
    let diff_vds = vds_i - vdseff;

    // VA (Early voltage)
    let vaclm = if sp.pclm > 0.0 && diff_vds > 1e-10 {
        let t0 = 1.0 / (sp.pclm * abulk * sp.litl);
        let t2 = vgsteff / esat_l;
        let t1 = leff * (abulk + t2);
        t0 * t1 * diff_vds
    } else {
        MAX_EXP
    };

    let vadibl = if sp.theta_rout > 0.0 {
        let t8 = abulk * vdsat;
        let t0 = vgst2vtm * t8;
        let t1 = vgst2vtm + t8;
        let vadibl = (vgst2vtm - t0 / t1) / sp.theta_rout;
        let t7 = sp.pdiblcb * vbseff;
        if t7 >= -0.9 {
            vadibl / (1.0 + t7)
        } else {
            let t4 = 1.0 / (0.8 + t7);
            vadibl * (17.0 + 20.0 * t7) * t4
        }
    } else {
        MAX_EXP
    };

    let t8 = sp.pvag / esat_l;
    let t9 = t8 * vgsteff;
    let t0 = if t9 > -0.9 {
        1.0 + t9
    } else {
        let t1 = 1.0 / (17.0 + 20.0 * t9);
        (0.8 + t9) * t1
    };

    let tmp3_va = vaclm + vadibl;
    let t1_va = vaclm * vadibl / tmp3_va;

    let tmp4 = 1.0 - 0.5 * abulk * vdsat / vgst2vtm;
    let t9_va = wvcox_rds * vgsteff;
    let t0_va = esat_l + vdsat + 2.0 * t9_va * tmp4;
    let t9_ab = wvcox_rds * abulk;
    let t1_ab = 2.0 / lambda - 1.0 + t9_ab;
    let vasat = t0_va / t1_ab;

    let va = vasat + t0 * t1_va;

    // Ids calculation
    let cox_wov_l = cox * weff_ch / leff;
    let beta = ueff * cox_wov_l;

    let t0_ids = 1.0 - 0.5 * abulk * vdseff / vgst2vtm;
    let fgche1 = vgsteff * t0_ids;
    let t9_fgche = vdseff / esat_l;
    let fgche2 = 1.0 + t9_fgche;
    let gche = beta * fgche1 / fgche2;
    let t0_gche = 1.0 + gche * rds;
    let t9_gche = vdseff / t0_gche;
    let idl = gche * t9_gche;

    let t9_ids = diff_vds / va;
    let t0_ids2 = 1.0 + t9_ids;
    let ids = idl * t0_ids2;

    // Derivatives
    let dgche_dvg = (beta * (t0_ids + vgsteff * (-0.5 * abulk / vgst2vtm))
        + fgche1 * dueff_dvg * cox_wov_l)
        / fgche2;
    let didl_dvg =
        (gche * dvdseff_dvg + t9_gche * dgche_dvg) / t0_gche - idl * gche / t0_gche * drds_dvg;
    let didl_dvd = (gche * dvdseff_dvd) / t0_gche;
    let didl_dvb = (gche * dvdseff_dvb - idl * drds_dvb * gche) / t0_gche;

    let dids_dvg = t0_ids2 * didl_dvg * dvgsteff_dvg;
    let dids_dvd = t0_ids2 * didl_dvd
        + idl * (1.0 - dvdseff_dvd) / va
        + dids_dvg / dvgs_eff_dvg * dvgsteff_dvd;
    let dids_dvb = (t0_ids2 * didl_dvb + t0_ids2 * didl_dvg * dvgsteff_dvb) * dvbseff_dvb;

    let gm = dids_dvg;
    let gds = dids_dvd + dids_dvg * dvgsteff_dvd / dvgs_eff_dvg.max(1e-20);
    let gmbs = dids_dvb;

    // GIDL current (drain side)
    let (igidl, ggidl_d, ggidl_g) = {
        let t0 = 3.0 * model.tox;
        let t1 = (vds_i - vgs_eff - sp.ngidl) / t0;
        if sp.agidl <= 0.0 || sp.bgidl <= 0.0 || t1 <= 0.0 {
            (0.0, 0.0, 0.0)
        } else {
            let dt1_dvd = 1.0 / t0;
            let dt1_dvg = -dt1_dvd * dvgs_eff_dvg;
            let t2 = sp.bgidl / t1;
            if t2 < EXPL_THRESHOLD {
                let igidl = sp.wdiod * sp.agidl * t1 * (-t2).exp();
                let t3 = igidl / t1 * (t2 + 1.0);
                (igidl, t3 * dt1_dvd, t3 * dt1_dvg)
            } else {
                let t3 = sp.wdiod * sp.agidl * MIN_EXPL;
                (t3 * t1, t3 * dt1_dvd, t3 * dt1_dvg)
            }
        }
    };

    // GIDL source side
    let (isgidl, gsgidl_g) = {
        let t0 = 3.0 * model.tox;
        let t1 = (-vgs_eff - sp.ngidl) / t0;
        if sp.agidl <= 0.0 || sp.bgidl <= 0.0 || t1 <= 0.0 {
            (0.0, 0.0)
        } else {
            let dt1_dvg = -dvgs_eff_dvg / t0;
            let t2 = sp.bgidl / t1;
            if t2 < EXPL_THRESHOLD {
                let isgidl = sp.wdios * sp.agidl * t1 * (-t2).exp();
                let t3 = isgidl / t1 * (t2 + 1.0);
                (isgidl, t3 * dt1_dvg)
            } else {
                let t3 = sp.wdios * sp.agidl * MIN_EXPL;
                (t3 * t1, t3 * dt1_dvg)
            }
        }
    };

    // Junction currents (4-component SOI model, same as PD)
    let nvtm1 = vtm * sp.ndiode;
    let vbd = vbs_i - vds_i;

    // Ibs1/Ibd1: Diffusion
    let (ibs1, dibs1_dvb, ibd1, dibd1_dvb, dibd1_dvd) = if sp.jdif == 0.0 {
        (0.0, 0.0, 0.0, 0.0, 0.0)
    } else {
        let t0 = vbs_i / nvtm1;
        let (exp_vbs, dexp) = soi_dexp(t0);
        let t0_ws = sp.wdios * model.tsi * sp.jdif;
        let ibs1 = t0_ws * (exp_vbs - 1.0);
        let dibs1_dvb = t0_ws * dexp / nvtm1;

        let t0 = vbd / nvtm1;
        let (exp_vbd, dexp) = soi_dexp(t0);
        let t0_wd = sp.wdiod * model.tsi * sp.jdif;
        let ibd1 = t0_wd * (exp_vbd - 1.0);
        let dibd1_dvb = t0_wd * dexp / nvtm1;
        (ibs1, dibs1_dvb, ibd1, dibd1_dvb, -dibd1_dvb)
    };

    // Ibs2/Ibd2: Recombination
    let (ibs2, dibs2_dvb, ibd2, dibd2_dvb, dibd2_dvd) = if sp.jrec == 0.0 {
        (0.0, 0.0, 0.0, 0.0, 0.0)
    } else {
        let nvtmf = 0.026 * model.nrecf0;
        let t0 = vbs_i / nvtmf;
        let (t10, t2) = soi_dexp(t0);
        let dt10_dvb = t2 / nvtmf;

        let t3 = sp.wdios * model.tsi * sp.jrec;
        let ibs2 = t3 * t10;
        let dibs2_dvb = t3 * dt10_dvb;

        let t0 = vbd / nvtmf;
        let (t10, t2) = soi_dexp(t0);
        let dt10_dvb = t2 / nvtmf;
        let t3 = sp.wdiod * model.tsi * sp.jrec;
        let ibd2 = t3 * t10;
        let dibd2_dvb = t3 * dt10_dvb;
        (ibs2, dibs2_dvb, ibd2, dibd2_dvb, -dibd2_dvb)
    };

    // Ibs3/Ibd3: BJT
    let (ibs3, dibs3_dvb, ibd3, dibd3_dvb, dibd3_dvd) = if sp.jbjt == 0.0 || sp.lratio == 0.0 {
        (0.0, 0.0, 0.0, 0.0, 0.0)
    } else {
        let t0_bs = vbs_i / nvtm1;
        let (exp_vbs, dexp_bs) = soi_dexp(t0_bs);
        let t0_bd = vbd / nvtm1;
        let (exp_vbd, dexp_bd) = soi_dexp(t0_bd);

        let ien = weff / sp.nseg * model.tsi * sp.jbjt * sp.lratio;
        let t0 = 1.0 - sp.arfabjt;
        if t0 < 1e-2 {
            (0.0, 0.0, 0.0, 0.0, 0.0)
        } else {
            let t1 = t0 * ien;
            let ibs3 = t1 * (exp_vbs - 1.0);
            let dibs3_dvb = t1 * dexp_bs / nvtm1;
            let ibd3 = t1 * (exp_vbd - 1.0);
            let dibd3_dvb = t1 * dexp_bd / nvtm1;
            (ibs3, dibs3_dvb, ibd3, dibd3_dvb, -dibd3_dvb)
        }
    };

    // Ibs4/Ibd4: Tunneling
    let (ibs4, dibs4_dvb, ibd4, dibd4_dvb, dibd4_dvd) = if sp.jtun == 0.0 {
        (0.0, 0.0, 0.0, 0.0, 0.0)
    } else {
        let nvtm_tun = vtm * model.ntun;
        let t0 = -vbs_i / nvtm_tun;
        let (exp_val, dexp_val) = soi_dexp(t0);
        let t3 = sp.wdios * model.tsi * sp.jtun;
        let ibs4 = -t3 * (exp_val - 1.0);
        let dibs4_dvb = t3 * dexp_val / nvtm_tun;

        let t0 = -vbd / nvtm_tun;
        let (exp_val, dexp_val) = soi_dexp(t0);
        let t3 = sp.wdiod * model.tsi * sp.jtun;
        let ibd4 = -t3 * (exp_val - 1.0);
        let dibd4_dvb = t3 * dexp_val / nvtm_tun;
        (ibs4, dibs4_dvb, ibd4, dibd4_dvb, -dibd4_dvb)
    };

    // Total junction currents
    let ibs = ibs1 + ibs2 + ibs3 + ibs4;
    let ibd = ibd1 + ibd2 + ibd3 + ibd4;
    let gbs_jct = dibs1_dvb + dibs2_dvb + dibs3_dvb + dibs4_dvb + gmin;
    let gbd_jct = dibd1_dvb + dibd2_dvb + dibd3_dvb + dibd4_dvb + gmin;

    // Impact ionization (DD uses ALPHA0/ALPHA1/BETA0, same as FD)
    let (iii, gii_d, gii_g, gii_b) = if sp.alpha0 <= 0.0 || vds_i <= sp.beta0.max(0.0) {
        (0.0, 0.0, 0.0, 0.0)
    } else {
        let t0 = vds_i - sp.beta0;
        let t1 = (-model.alpha1 / t0.max(1e-20)).exp();
        let iii = sp.alpha0 * ids * t1;
        let gii_d = sp.alpha0 * (gds * t1 + ids * t1 * model.alpha1 / (t0 * t0));
        let gii_g = sp.alpha0 * gm * t1;
        let gii_b = sp.alpha0 * gmbs * t1;
        (iii, gii_d, gii_g, gii_b)
    };

    // Body node current balance
    let ceq_d = sign * (ids - gm * vgs_i - gds * vds_i - gmbs * vbs_i);
    let ceq_bs = -(ibs - gbs_jct * vbs_i);
    let ceq_bd = -(ibd - gbd_jct * vbd);

    // Capacitances
    let cox_wl = cox * weff_ch * leff;
    let (cggb, cgdb, cgsb) = if vgsteff > 0.0 {
        let t0 = 1.0 - abulk * vdseff / (2.0 * vgst2vtm);
        (cox_wl * (1.0 - t0 * t0), -cox_wl * t0 * 0.5, 0.0)
    } else {
        (cox_wl * 0.05, 0.0, 0.0)
    };
    let cbgb = 0.0;
    let cbdb = 0.0;
    let cbsb = 0.0;
    let cdgb = -cggb - cgdb;
    let cddb = -cgdb;
    let cdsb = -(cdgb + cddb);

    let qinv = if vgsteff > 0.0 {
        cox_wl * vgsteff * (1.0 - 0.5 * abulk * vdseff / vgst2vtm)
    } else {
        0.0
    };

    Bsim3SoiDdCompanion {
        ids: ids / sp.nseg,
        gm: gm / sp.nseg,
        gds: gds / sp.nseg,
        gmbs: gmbs / sp.nseg,
        mode,
        vdsat,
        ibs,
        ibd,
        gbs_jct,
        gbd_jct,
        iii,
        gii_d,
        gii_g,
        gii_b,
        igidl,
        ggidl_d,
        ggidl_g,
        isgidl,
        gsgidl_g,
        ceq_d,
        ceq_bs,
        ceq_bd,
        cggb: cggb + sp.cgso_eff * weff + sp.cgdo_eff * weff,
        cgdb: cgdb - sp.cgdo_eff * weff,
        cgsb: cgsb - sp.cgso_eff * weff,
        cbgb,
        cbdb,
        cbsb,
        cdgb: cdgb - sp.cgdo_eff * weff,
        cddb: cddb + sp.cgdo_eff * weff,
        cdsb,
        capbd: 0.0,
        capbs: 0.0,
        qinv,
        vbs_dd,
    }
}

/// Stamp BSIM3SOI-DD companion model into the MNA matrix and RHS.
pub fn stamp_bsim3soi_dd(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &Bsim3SoiDdInstance,
    comp: &Bsim3SoiDdCompanion,
) {
    let dp = inst.drain_eff_idx();
    let g = inst.gate_idx;
    let sp_idx = inst.source_eff_idx();
    let b = inst.body_int_idx;

    let sign = inst.model.mos_type.sign();
    let m = inst.m;

    let (xnrm, xrev) = if comp.mode > 0 {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };

    let gm = m * comp.gm;
    let gds = m * comp.gds;
    let gmbs = m * comp.gmbs;
    let gbd = m * comp.gbd_jct;
    let gbs = m * comp.gbs_jct;

    let gm_eff = gm * xnrm + comp.gds * xrev * m;
    let gds_eff = gds * xnrm + comp.gm * xrev * m;
    let gmbs_eff = gmbs;

    // D' row
    crate::stamp_conductance(matrix, dp, dp, gds_eff + gbd);
    crate::stamp_conductance(matrix, dp, g, gm_eff);
    crate::stamp_conductance(matrix, dp, sp_idx, -(gm_eff + gds_eff + gmbs_eff));
    crate::stamp_conductance(matrix, dp, b, gmbs_eff - gbd);

    // S' row
    crate::stamp_conductance(matrix, sp_idx, dp, -(gm_eff + gds_eff));
    crate::stamp_conductance(matrix, sp_idx, g, -gm_eff * xnrm + gds * xrev * m);
    crate::stamp_conductance(matrix, sp_idx, sp_idx, gm_eff + gds_eff + gmbs_eff + gbs);
    crate::stamp_conductance(matrix, sp_idx, b, -gmbs_eff - gbs);

    // B row
    crate::stamp_conductance(matrix, b, dp, -gbd);
    crate::stamp_conductance(matrix, b, sp_idx, -gbs);
    crate::stamp_conductance(matrix, b, b, gbd + gbs);

    // Impact ionization
    if comp.iii != 0.0 {
        let gii_d = m * comp.gii_d;
        let gii_g = m * comp.gii_g;
        let gii_b = m * comp.gii_b;
        crate::stamp_conductance(matrix, dp, dp, -gii_d);
        crate::stamp_conductance(matrix, dp, g, -gii_g);
        crate::stamp_conductance(matrix, dp, b, -gii_b);
        crate::stamp_conductance(matrix, b, dp, gii_d);
        crate::stamp_conductance(matrix, b, g, gii_g);
        crate::stamp_conductance(matrix, b, b, gii_b);
    }

    // RHS stamps
    let ceq_d = sign * m * comp.ceq_d;
    let ceq_bs = sign * m * comp.ceq_bs;
    let ceq_bd = sign * m * comp.ceq_bd;

    if let Some(d) = dp {
        rhs[d] -= ceq_d + ceq_bd;
    }
    if let Some(s) = sp_idx {
        rhs[s] += ceq_d + ceq_bs;
    }
    if let Some(bulk) = b {
        rhs[bulk] -= ceq_bs + ceq_bd;
    }

    // DD body feedback: drive internal body node towards DD-computed Vbs.
    // DD has junction currents (unlike FD), but floating body still needs feedback
    // to establish the body potential from the self-consistent surface potential chain.
    if inst.body_idx.is_none() {
        // Floating body: use feedback conductance
        let gbody_fb = m * 1.0; // 1 S feedback conductance
        crate::stamp_conductance(matrix, b, sp_idx, gbody_fb);
        let i_fb = gbody_fb * comp.vbs_dd * sign;
        if let Some(bi) = b {
            rhs[bi] += i_fb;
        }
        if let Some(si) = sp_idx {
            rhs[si] -= i_fb;
        }
    }

    // Body resistance to external body contact (if present)
    if let (Some(b_int), Some(b_ext)) = (inst.body_int_idx, inst.body_idx) {
        let gbody = if inst.model.rbody > 0.0 {
            m / inst.model.rbody
        } else {
            m * 1e3
        };
        crate::stamp_conductance(matrix, Some(b_int), Some(b_ext), -gbody);
        crate::stamp_conductance(matrix, Some(b_ext), Some(b_int), -gbody);
        crate::stamp_conductance(matrix, Some(b_int), Some(b_int), gbody);
        crate::stamp_conductance(matrix, Some(b_ext), Some(b_ext), gbody);
    }

    // Series resistance: D<->D', S<->S'
    if inst.drain_prime_idx.is_some() && inst.drain_prime_idx != inst.drain_idx {
        let rd = if inst.nrd > 0.0 && inst.model.rbsh > 0.0 {
            inst.model.rbsh * inst.nrd
        } else {
            0.01
        };
        let grd = m / rd;
        crate::stamp_conductance(matrix, inst.drain_idx, inst.drain_idx, grd);
        crate::stamp_conductance(matrix, inst.drain_idx, inst.drain_prime_idx, -grd);
        crate::stamp_conductance(matrix, inst.drain_prime_idx, inst.drain_idx, -grd);
        crate::stamp_conductance(matrix, inst.drain_prime_idx, inst.drain_prime_idx, grd);
    }
    if inst.source_prime_idx.is_some() && inst.source_prime_idx != inst.source_idx {
        let rs = if inst.nrs > 0.0 && inst.model.rbsh > 0.0 {
            inst.model.rbsh * inst.nrs
        } else {
            0.01
        };
        let grs = m / rs;
        crate::stamp_conductance(matrix, inst.source_idx, inst.source_idx, grs);
        crate::stamp_conductance(matrix, inst.source_idx, inst.source_prime_idx, -grs);
        crate::stamp_conductance(matrix, inst.source_prime_idx, inst.source_idx, -grs);
        crate::stamp_conductance(matrix, inst.source_prime_idx, inst.source_prime_idx, grs);
    }
}

/// BSIM3SOI-DD voltage limiting for NR convergence.
pub fn bsim3soi_dd_limit(
    vgs_new: f64,
    vds_new: f64,
    vbs_new: f64,
    ves_new: f64,
    vgs_old: f64,
    vds_old: f64,
    vbs_old: f64,
    ves_old: f64,
    vth: f64,
) -> (f64, f64, f64, f64) {
    let vgs = crate::bsim3::fetlim(vgs_new, vgs_old, vth);
    let vds = crate::bsim3::fetlim(vds_new, vds_old, vth);
    let limit = 5.0;
    let vbs = if (vbs_new - vbs_old).abs() > limit {
        if vbs_new > vbs_old {
            vbs_old + limit
        } else {
            vbs_old - limit
        }
    } else {
        vbs_new
    };
    let ves = if (ves_new - ves_old).abs() > limit {
        if ves_new > ves_old {
            ves_old + limit
        } else {
            ves_old - limit
        }
    } else {
        ves_new
    };
    (vgs, vds, vbs, ves)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dd_model_defaults() {
        let model = Bsim3SoiDdModel::new(MosfetType::Nmos);
        assert_eq!(model.mos_type, MosfetType::Nmos);
        assert!((model.vth0 - 0.7).abs() < 1e-10);
        assert!(model.cox > 0.0);
        assert!(model.phi > 0.0);
        assert!(model.cbox > 0.0);
        assert!(model.csi > 0.0);
        assert!(model.qsi > 0.0);
    }

    #[test]
    fn test_dd_model_pmos_defaults() {
        let model = Bsim3SoiDdModel::new(MosfetType::Pmos);
        assert_eq!(model.mos_type, MosfetType::Pmos);
        assert!((model.vth0 - (-0.7)).abs() < 1e-10);
    }

    #[test]
    fn test_dd_size_dep_param() {
        let model = Bsim3SoiDdModel::new(MosfetType::Nmos);
        let sp = model.size_dep_param(10e-6, 0.25e-6, TEMP_DEFAULT);
        assert!(sp.leff > 0.0);
        assert!(sp.weff > 0.0);
        assert!(sp.u0 > 0.0);
        assert!(sp.vsat > 0.0);
    }

    #[test]
    fn test_dd_companion_zero_bias() {
        let model = Bsim3SoiDdModel::new(MosfetType::Nmos);
        let sp = model.size_dep_param(10e-6, 0.25e-6, TEMP_DEFAULT);
        let comp = bsim3soi_dd_companion(0.0, 0.0, 0.0, 0.0, &sp, &model);
        // At zero bias, Ids should be very small (subthreshold)
        assert!(comp.ids.abs() < 1e-3);
        assert!(comp.gm.is_finite());
        assert!(comp.gds.is_finite());
    }

    #[test]
    fn test_dd_companion_on_state() {
        let model = Bsim3SoiDdModel::new(MosfetType::Nmos);
        let sp = model.size_dep_param(10e-6, 0.25e-6, TEMP_DEFAULT);
        let comp = bsim3soi_dd_companion(1.5, 1.0, 0.0, 0.0, &sp, &model);
        // In strong inversion with positive Vds, should have meaningful current
        assert!(comp.ids > 0.0);
        assert!(comp.gm > 0.0);
        assert!(comp.gds > 0.0);
    }

    #[test]
    fn test_dd_companion_reverse() {
        let model = Bsim3SoiDdModel::new(MosfetType::Nmos);
        let sp = model.size_dep_param(10e-6, 0.25e-6, TEMP_DEFAULT);
        let comp = bsim3soi_dd_companion(1.5, -0.5, 0.0, 0.0, &sp, &model);
        assert_eq!(comp.mode, -1);
        assert!(comp.ids.is_finite());
    }

    #[test]
    fn test_dd_voltage_limiting() {
        let (vgs, vds, vbs, ves) =
            bsim3soi_dd_limit(10.0, 10.0, 10.0, 10.0, 0.5, 0.5, 0.0, 0.0, 0.7);
        // VGS and VDS should be limited
        assert!(vgs < 10.0);
        assert!(vds < 10.0);
        // VBS and VES limited to ±5V from old values
        assert!(vbs <= 5.0);
        assert!(ves <= 5.0);
    }
}
