//! BSIM3SOI-FD (Fully Depleted Silicon-On-Insulator) MOSFET model.
//!
//! Implements the BSIM3SOI-FD v2.1 model matching ngspice level 55.
//! Based on the BSIM3v3 core equations with SOI-specific extensions for
//! fully depleted devices: self-consistent surface potential calculation,
//! back-gate (E node), optional body contact (B/P node), and simplified
//! junction model (no parasitic BJT, no GIDL).

#![allow(unused_variables, dead_code, clippy::too_many_arguments, unused_parens)]

use thevenin_types::{Expr, ModelDef};

use crate::mosfet::MosfetType;
use crate::physics::{
    CHARGE_Q, EG300, EPSOX, EPSSI, EXP_THRESHOLD, KBOQ, MAX_EXP, MIN_EXP, bsim_safe_exp as safe_exp,
};

const DELTA_1: f64 = 0.02;
const DELTA_4: f64 = 0.02;
const DELT_VBS0EFF: f64 = 0.02;
const DELT_VBSMOS: f64 = 0.005;
const DELT_VBSEFF: f64 = 0.005;
const DELT_VBS0DIO: f64 = 1e-7;
/// Smoothing delta for Vbsdio clamp (from ngspice #define DELT_Vbsdio)
const DELT_VBSDIO: f64 = 0.01;
/// Offset for Vbsdio clamp floor (from ngspice #define OFF_Vbsdio)
const OFF_VBSDIO: f64 = 0.02;

const TEMP_DEFAULT: f64 = 300.15;

/// BSIM3SOI-FD model parameters (from .model card, Level=55).
#[derive(Debug, Clone)]
pub struct Bsim3SoiFdModel {
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

    // Impact ionization (FD-specific: AII/BII/CII/DII)
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

    // SOI junction model
    pub ndiode: f64,
    pub ntun: f64,
    pub isbjt: f64,
    pub isdif: f64,
    pub istun: f64,
    pub isrec: f64,
    pub xbjt: f64,
    pub xdif: f64,
    pub xrec: f64,
    pub xtun: f64,

    // FD-specific surface potential params
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
    // FD-specific precomputed
    pub cbox: f64,
    pub csi: f64,
    pub qsi: f64,
    pub csieff: f64,
    pub qsieff: f64,
    pub vfbb: f64,
}

/// Size-dependent parameters for BSIM3SOI-FD.
#[derive(Debug, Clone)]
pub struct Bsim3SoiFdSizeParam {
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

    // ua/ub/uc temperature-adjusted
    pub uatemp: f64,
    pub ubtemp: f64,
    pub uctemp: f64,

    // rds0 denominator
    pub rds0denom: f64,

    // FD-specific
    pub dvbd0: f64,
    pub dvbd1: f64,
}

/// BSIM3SOI-FD instance with node indices.
#[derive(Debug, Clone)]
pub struct Bsim3SoiFdInstance {
    pub name: String,
    pub drain_idx: Option<usize>,
    pub gate_idx: Option<usize>,
    pub source_idx: Option<usize>,
    /// Back-gate (E) node — always present for SOI.
    pub e_idx: Option<usize>,
    /// External body contact (B/P) node — optional (floating body when None).
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
    pub model: Bsim3SoiFdModel,
    pub size_params: Bsim3SoiFdSizeParam,
    pub vth0_inst: f64,
    /// Number of body contacts (0=floating, 1=single, 2=double)
    pub nbc: f64,
}

/// NR companion result for BSIM3SOI-FD.
#[derive(Debug, Clone)]
pub struct Bsim3SoiFdCompanion {
    pub ids: f64,
    pub gm: f64,
    pub gds: f64,
    pub gmbs: f64,
    pub mode: i32,
    pub vdsat: f64,

    // FD has no junction currents (Ibs=Ibd=0)
    pub gbs_jct: f64,
    pub gbd_jct: f64,

    // Impact ionization
    pub iii: f64,
    pub gii_d: f64,
    pub gii_g: f64,
    pub gii_b: f64,

    // Equivalent current sources for NR companion
    pub ceq_d: f64,
    pub ceq_bs: f64,
    pub ceq_bd: f64,
    pub ceq_iii: f64,

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

    // FD-computed body-source voltage (for body node feedback)
    pub vbs_fd: f64,
}

impl Bsim3SoiFdModel {
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
            isbjt: 1e-6,
            isdif: 0.0,
            istun: 0.0,
            isrec: 0.0,
            xbjt: 1.0,
            xdif: 1.0,
            xrec: 1.0,
            xtun: 0.0,
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
        // NCH and NPEAK are aliases in ngspice (IOP "nch" -> B3SOIFD_MOD_NPEAK).
        set!(nch, "NCH");
        set!(npeak, "NPEAK");
        if m.nch != 1.7e17 && m.npeak == 1.7e17 {
            m.npeak = m.nch;
        } else if m.npeak != 1.7e17 && m.nch == 1.7e17 {
            m.nch = m.npeak;
        }
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
        set!(isbjt, "ISBJT");
        set!(isdif, "ISDIF");
        set!(istun, "ISTUN");
        set!(isrec, "ISREC");
        set!(xbjt, "XBJT");
        set!(xdif, "XDIF");
        set!(xrec, "XREC");
        set!(xtun, "XTUN");
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

        // vbi = Vt * ln(ND * NA / ni²), ND = 1e20 /cm³ (n+ S/D), NA = npeak /cm³.
        self.vbi_default = self.vtm * (1e20 * npeak / (self.ni * self.ni)).ln();
        self.factor1 = (EPSSI / EPSOX * self.tox).sqrt();

        // FD-specific: box and silicon capacitances
        self.cbox = EPSOX / self.tbox;
        self.csi = EPSSI / self.tsi;
        self.qsi = CHARGE_Q * npeak * 1e6 * self.tsi;

        // Effective silicon capacitance (for body potential calculation)
        // csieff = csi / 2 (half of silicon film for FD)
        self.csieff = self.csi * 0.5;
        self.qsieff = self.qsi * 0.5;

        // Flat-band voltage
        let nsub = if self.nsub > 1e20 {
            self.nsub * 1e-6
        } else {
            self.nsub
        };
        if nsub > 0.0 {
            self.vfbb = -self.mos_type.sign() * self.vtm * (npeak / nsub).ln();
        } else {
            self.vfbb =
                -self.mos_type.sign() * self.vtm * (-npeak * nsub / (self.ni * self.ni)).ln();
        }
    }

    /// Number of internal nodes this model creates.
    pub fn internal_node_count(&self, nrd: f64, nrs: f64) -> usize {
        self.internal_node_count_fd(nrd, nrs, true)
    }

    /// FD-specific internal node count: body node only created when body contact exists.
    pub fn internal_node_count_fd(&self, nrd: f64, nrs: f64, has_body_contact: bool) -> usize {
        let mut count = if has_body_contact { 1 } else { 0 };
        if nrd > 0.0 && self.rbsh > 0.0 {
            count += 1; // drain prime
        }
        if nrs > 0.0 && self.rbsh > 0.0 {
            count += 1; // source prime
        }
        count
    }

    pub fn size_dep_param(&self, w: f64, l: f64, temp: f64) -> Bsim3SoiFdSizeParam {
        let tnom_k = self.tnom + 273.15;
        let vtm = KBOQ * temp;

        // Effective length/width
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

        // Temperature ratio
        let temp_ratio = temp / tnom_k;
        let temp_ratio_minus1 = temp_ratio - 1.0;

        // Mobility temperature dependence
        let u0temp = self.u0 * temp_ratio.powf(self.ute);

        // Velocity saturation temperature dependence
        let vsattemp = self.vsat - self.at * (temp - tnom_k);

        // Series resistance
        let wr2 = 1.0 / self.wr;
        let rds0denom = (weff * 1e6).powf(wr2);
        let rds0 = (self.rdsw + self.prt * temp_ratio_minus1) / rds0denom;

        // ua/ub/uc with temperature
        let uatemp = self.ua + self.ua1 * temp_ratio_minus1;
        let ubtemp = self.ub + self.ub1 * temp_ratio_minus1;
        let uctemp = self.uc + self.uc1 * temp_ratio_minus1;

        // Depletion width
        let xdep0 = (2.0 * EPSSI / (CHARGE_Q * self.npeak * 1e6)).sqrt() * sqrt_phi;

        // litl (characteristic length)
        let litl = (EPSSI * self.tox / self.cox).sqrt();

        // Theta0vb0 for DIBL
        let t0 = -0.5 * self.dsub * leff / litl;
        let theta0vb0 = if t0 > -EXP_THRESHOLD {
            let t1 = t0.exp();
            self.dvt0 * t1 * (1.0 + 2.0 * t1)
        } else {
            self.dvt0 * MIN_EXP * (1.0 + 2.0 * MIN_EXP)
        };

        // Theta for output resistance (theta_rout)
        let t0 = -0.5 * self.drout * leff / litl;
        let theta_rout = if self.pdiblc1 == 0.0 && self.pdiblc2 == 0.0 {
            0.0
        } else if t0 > -EXP_THRESHOLD {
            let t1 = t0.exp();
            let t2 = t1 * (1.0 + 2.0 * t1);
            self.pdiblc1 * t2 + self.pdiblc2
        } else {
            self.pdiblc2
        };

        // k1eff
        let k1eff = self.k1;

        // Flat-band voltage
        let vfb = self.vbi_default - phi;

        // cdep0
        let cdep0 = (2.0 * EPSSI * CHARGE_Q * self.npeak * 1e6).sqrt() / (2.0 * vtm);

        // Junction precomputed (temperature-dependent)
        let t3_temp = temp_ratio - 1.0;
        let t4_base = EG300 / self.ndiode / vtm * t3_temp;
        let t5_exp = t4_base.exp();
        let t6_sqrt = t5_exp.sqrt();

        let t7_bjt = self.xbjt / self.ndiode;
        let t0_bjt = temp_ratio.powf(t7_bjt);
        let t7_dif = self.xdif / self.ndiode;
        let t1_dif = temp_ratio.powf(t7_dif);
        let t7_rec = self.xrec / self.ndiode / 2.0;
        let t2_rec = temp_ratio.powf(t7_rec);

        let jbjt = self.isbjt * t0_bjt * t5_exp;
        let jdif = self.isdif * t1_dif * t5_exp;
        let jrec = self.isrec * t2_rec * t6_sqrt;

        let t7_tun = self.xtun / self.ntun;
        let t0_tun = temp_ratio.powf(t7_tun);
        let jtun = self.istun * t0_tun;

        // Diode widths
        let wdios = weff * self.asd;
        let wdiod = weff * self.asd;

        // Overlap caps
        let cgso_eff = self.cgso.max(0.0);
        let cgdo_eff = self.cgdo.max(0.0);

        // npeak for this size (no binning, just base values)
        let npeak = self.npeak;
        let nsub = self.nsub;

        Bsim3SoiFdSizeParam {
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
            u0: self.u0,
            ua: self.ua,
            ub: self.ub,
            uc: self.uc,
            vsat: self.vsat,
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
            vbi: self.vbi_default,
            wdios,
            wdiod,
            jbjt,
            jdif,
            jrec,
            jtun,
            cgso_eff,
            cgdo_eff,
            uatemp,
            ubtemp,
            uctemp,
            rds0denom,
            dvbd0: self.dvbd0,
            dvbd1: self.dvbd1,
        }
    }
}

impl Bsim3SoiFdInstance {
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

    pub fn ac_stamp(&self, comp: &Bsim3SoiFdCompanion) -> crate::ac::BsimAcStamp {
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

/// Compute BSIM3SOI-FD companion model (NR linearization).
///
/// Core load function computing drain current, conductances, and capacitances.
/// FD model uses self-consistent surface potential calculation (Vbs0t, Vbs0, Vbsmos, Vbseff)
/// instead of the external Vbs used in PD. Junction currents are zero in FD.
#[expect(clippy::too_many_lines)]
pub fn bsim3soi_fd_companion(
    vgs: f64,
    vds: f64,
    vbs: f64,
    ves: f64,
    sp: &Bsim3SoiFdSizeParam,
    model: &Bsim3SoiFdModel,
    floating_body: bool,
) -> Bsim3SoiFdCompanion {
    let sign = model.mos_type.sign();
    let cox = model.cox;
    let vtm = KBOQ * TEMP_DEFAULT;
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

    // ========== FD-specific surface potential calculation ==========

    // Prepare Vbs0t (surface potential at Vgs=Vth)
    let t0 = -sp.dvbd1 * leff / sp.litl;
    let t1 = sp.dvbd0 * (safe_exp(0.5 * t0) + 2.0 * safe_exp(t0));
    let t2 = t1 * v0;
    let t3 = 0.5 * model.qsi / model.csi;
    let vbs0t = phi - t3 + model.vbsa + t2;

    // Prepare Vbs0 (accounting for back-gate coupling)
    let t0_kb = 1.0 + model.csieff / model.cbox;
    let t1_kb = model.kb1 / t0_kb;
    let t2_kb = t1_kb * (vbs0t - vesfb);
    let t6_vbs0 = vbs0t - t2_kb;

    // Limit Vbs0 below phi - delp
    let t1_lim = phi - model.delp;
    let t2_lim = t1_lim - t6_vbs0 - DELT_VBSEFF;
    let t3_lim = (t2_lim * t2_lim + 4.0 * DELT_VBSEFF).sqrt();
    let vbs0 = t1_lim - 0.5 * (t2_lim + t3_lim);

    // Vbs0mos (corrected for silicon capacitance)
    let t1_mos = vbs0t - vbs0 - DELT_VBSMOS;
    let t2_mos = (t1_mos * t1_mos + DELT_VBSMOS * DELT_VBSMOS).sqrt();
    let t3_mos = 0.5 * (t1_mos + t2_mos);
    let t4_mos = t3_mos * model.csieff / model.qsieff;
    let vbs0mos = vbs0 - 0.5 * t3_mos * t4_mos;

    // Prepare Vthfd (FD threshold voltage using Vbs0mos)
    let phis_fd = phi - vbs0mos;
    let sqrt_phis_fd = phis_fd.abs().sqrt();
    let xdep_fd = sp.xdep0 * sqrt_phis_fd / sqrt_phi;

    // SCE for Vthfd
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

    // Width effect
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

    // Temperature contribution
    let temp_ratio_minus1 = (TEMP_DEFAULT / (model.tnom + 273.15)) - 1.0;
    let t0_nlx = (1.0 + sp.nlx / leff).sqrt();
    let t1_kt = sp.kt1 + sp.kt1l / leff + sp.kt2 * vbs0mos;
    let delt_vth_temp_fd = sp.k1 * (t0_nlx - 1.0) * sqrt_phi + t1_kt * temp_ratio_minus1;

    let tmp2 = model.tox * phi / (weff + sp.w0);

    // DIBL for Vthfd
    let t3_eta = sp.eta0 + sp.etab * vbs0mos;
    let t3_eta_eff = if t3_eta < 1e-4 {
        let t9 = 1.0 / (3.0 - 2e4 * t3_eta);
        (2e-4 - t3_eta) * t9
    } else {
        t3_eta
    };
    let dibl_sft_fd = t3_eta_eff * sp.theta0vb0 * vds_i;

    let vthfd = sign * sp.vth0 + sp.k1 * (sqrt_phis_fd - sqrt_phi)
        - sp.k2 * vbs0mos
        - delt_vth_fd
        - delt_vthw_fd
        + (sp.k3 + sp.k3b * vbs0mos) * tmp2
        + delt_vth_temp_fd
        - dibl_sft_fd;

    // ========== Poly Gate Si Depletion Effect ==========
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

    // ========== FD Vbs0eff and Vbsmos calculation ==========

    // Effective Vbs0 for all Vgs (depends on whether we're in strong or weak inversion)
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

    let vbs0eff_fd = vbs0 - nfb * 0.5 * (t1_eff + t2_eff);

    // Vbsdio for FD: when body is floating (bodyMod=0), Vbsdio = Vbs0eff, dVbsdio/dVb = 0
    // (matching b3soifdld.c line 1090). When body has an external contact, use smooth-max.
    let vbsdio = if floating_body {
        // FD floating body: channel current independent of body voltage
        vbs0eff_fd
    } else {
        let t1_vbsdio = vbs_i - (vbs0eff_fd + OFF_VBSDIO) - DELT_VBSDIO;
        let t2_vbsdio = (t1_vbsdio * t1_vbsdio + DELT_VBSDIO * DELT_VBSDIO).sqrt();
        vbs0eff_fd + OFF_VBSDIO + 0.5 * (t1_vbsdio + t2_vbsdio)
    };

    // Prepare Vbsmos (account for silicon capacitance correction)
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

    // ========== Main MOSFET equations (same structure as PD) ==========
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

    let t0_sce = -0.5 * sp.dvt1 * leff / lt1;
    let (theta0, dtheta0_dvb) = if t0_sce > -EXP_THRESHOLD {
        let t1 = t0_sce.exp();
        (
            t1 * (1.0 + 2.0 * t1),
            (-t0_sce / lt1 * t1 * dlt1_dvb) * (1.0 + 4.0 * t1),
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

    let t3_eta_main = sp.eta0 + sp.etab * vbseff;
    let (t3_eta_main_eff, dt3_dvb_main) = if t3_eta_main < 1e-4 {
        let t9 = 1.0 / (3.0 - 2e4 * t3_eta_main);
        ((2e-4 - t3_eta_main) * t9, t9 * t9 * sp.etab)
    } else {
        (t3_eta_main, sp.etab)
    };
    let dibl_sft = t3_eta_main_eff * sp.theta0vb0 * vds_i;
    let ddibl_sft_dvd = sp.theta0vb0 * t3_eta_main_eff;
    let ddibl_sft_dvb = sp.theta0vb0 * vds_i * dt3_dvb_main;

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
    let t10 = 2.0 * n * vtm;
    let vgst_nvt = vgst / t10;
    let exp_arg = (2.0 * sp.voff - vgst) / t10;

    let (vgsteff, dvgsteff_dvg, dvgsteff_dvd, dvgsteff_dvb) = if vgst_nvt > EXP_THRESHOLD {
        (vgst, dvgs_eff_dvg, -dvth_dvd, -dvth_dvb)
    } else if exp_arg > EXP_THRESHOLD {
        let t0 = (vgst - sp.voff) / (n * vtm);
        let exp_vgst = t0.exp();
        let vgsteff_val = vtm * sp.cdep0 / cox * exp_vgst;
        let t3 = vgsteff_val / (n * vtm);
        let t1 = -t3 * (dvth_dvb + t0 * vtm * dn_dvb);
        (
            vgsteff_val,
            t3 * dvgs_eff_dvg,
            -t3 * (dvth_dvd + t0 * vtm * dn_dvd),
            t1,
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
            t4,
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
        let rds = sp.rds0 * (0.8 + t0_rds) * t1;
        let t1sq = t1 * t1;
        (
            rds,
            sp.rds0 * sp.prwg * t1sq,
            sp.rds0 * sp.prwb * dsqrt_phis_dvb * t1sq,
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

        let t10_k1 = 0.5 * sp.k1 / sqrt_phi;
        let t1 = t10_k1 * t11;
        let dt1_dvb = t10_k1 * dt11_dvb;

        let t9 = (model.tox.max(1e-12) * xdep).sqrt();
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
        let mut dabulk_dvb =
            dabulk0_dvb - t8 * vgsteff * (dt1_dvb + 3.0 * t1 * dt2_dvb / tmp2_a.max(1e-20));

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
    let (ueff, dueff_dvg, dueff_dvd, dueff_dvb) = {
        let t5 = if model.mob_mod == 1 {
            let t0 = vgsteff + vth + vth;
            let t2 = sp.uatemp + sp.uctemp * vbseff;
            let t3 = t0 / model.tox;
            t3 * (t2 + sp.ubtemp * t3)
        } else if model.mob_mod == 2 {
            vgsteff / model.tox * (sp.uatemp + sp.uctemp * vbseff + sp.ubtemp * vgsteff / model.tox)
        } else {
            let t0 = vgsteff + vth + vth;
            let t2 = 1.0 + sp.uctemp * vbseff;
            let t3 = t0 / model.tox;
            t3 * (sp.uatemp + sp.ubtemp * t3) * t2
        };

        let denomi = if t5 >= -0.8 {
            1.0 + t5
        } else {
            let t9 = 1.0 / (7.0 + 10.0 * t5);
            (0.6 + t5) * t9
        };

        let ueff = sp.u0temp / denomi;
        let t9 = -ueff / denomi;
        let dueff_dvg_raw = if model.mob_mod == 1 || model.mob_mod == 3 {
            t9 * (sp.uatemp + 2.0 * sp.ubtemp * (vgsteff + vth + vth) / model.tox) / model.tox
        } else {
            t9 * (sp.uatemp + sp.uctemp * vbseff + 2.0 * sp.ubtemp * vgsteff / model.tox)
                / model.tox
        };
        (ueff, dueff_dvg_raw, 0.0, 0.0)
    };

    // Saturation voltage Vdsat
    let wvcox = weff_ch * sp.vsattemp * cox;
    let wvcox_rds = wvcox * rds;
    let esat = 2.0 * sp.vsattemp / ueff;
    let esat_l = esat * leff;
    let t0_esat = -esat_l / ueff;
    let desat_l_dvg = t0_esat * dueff_dvg;
    let desat_l_dvd = t0_esat * dueff_dvd;
    let desat_l_dvb = t0_esat * dueff_dvb;

    // Lambda
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

        let disc = (t1 * t1 - 2.0 * t0 * t2).max(0.0);
        let t3 = disc.sqrt();
        vdsat = (t1 - t3) / t0;
        dvdsat_dvg = (dt1_dvg
            - (t1 * dt1_dvg - dt0_dvg * t2 - t0 * dt2_dvg) / t3.max(1e-30)
            - vdsat * dt0_dvg)
            / t0;
        dvdsat_dvb = (dt1_dvb
            - (t1 * dt1_dvb - dt0_dvb * t2 - t0 * dt2_dvb) / t3.max(1e-30)
            - vdsat * dt0_dvb)
            / t0;
        dvdsat_dvd = (dt1_dvd - (t1 * dt1_dvd - t0 * dt2_dvd) / t3.max(1e-30)) / t0;
    }

    // Vdseff (effective Vds)
    let t1 = vdsat - vds_i - sp.delta;
    let t2 = (t1 * t1 + 4.0 * sp.delta * vdsat).sqrt();
    let t0 = t1 / t2;
    let t3 = 2.0 * sp.delta / t2;
    let vdseff = vdsat - 0.5 * (t1 + t2);
    let dvdseff_dvg = dvdsat_dvg - 0.5 * (dvdsat_dvg + t0 * dvdsat_dvg + t3 * dvdsat_dvg);
    let dvdseff_dvd =
        dvdsat_dvd - 0.5 * (dvdsat_dvd - 1.0 + t0 * (dvdsat_dvd - 1.0) + t3 * dvdsat_dvd);
    let dvdseff_dvb = dvdsat_dvb - 0.5 * (dvdsat_dvb + t0 * dvdsat_dvb + t3 * dvdsat_dvb);

    let (vdseff, dvdseff_dvg, dvdseff_dvd, dvdseff_dvb) = if vdseff > vds_i {
        (vds_i, 0.0_f64, 1.0_f64, 0.0_f64)
    } else {
        (vdseff, dvdseff_dvg, dvdseff_dvd, dvdseff_dvb)
    };
    let diff_vds = vds_i - vdseff;

    // VA (Early voltage: VACLM || VADIBL)
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

    // Combined VA
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

    // Derivatives: gm, gds, gmbs
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
        + dids_dvg / dvgs_eff_dvg.max(1e-20) * dvgsteff_dvd;
    let dids_dvb = t0_ids2 * didl_dvb + t0_ids2 * didl_dvg * dvgsteff_dvb;

    let gm = dids_dvg;
    let gds = dids_dvd + dids_dvg * dvgsteff_dvd / dvgs_eff_dvg.max(1e-20);
    // FD floating body: dIds/dVb = 0 because Vbsdio = Vbs0eff is independent of Vb
    // (dVbsdio/dVb = 0, see b3soifdld.c line 1095)
    let gmbs = if floating_body { 0.0 } else { dids_dvb };

    // FD: Junction currents are ALL ZERO
    let gbs_jct = 0.0;
    let gbd_jct = 0.0;

    // Impact ionization (FD uses AII/BII/CII/DII)
    let (iii, gii_d, gii_g, gii_b) = if model.alpha0 <= 0.0 || vds_i <= sp.beta0.max(0.0) {
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
    let ceq_bs = 0.0; // No junction current in FD
    let ceq_bd = 0.0;
    let ceq_iii = iii - gii_d * vds_i - gii_g * vgs_i - gii_b * vbs_i;

    // Capacitances (simplified — gate overlap + basic intrinsic)
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

    Bsim3SoiFdCompanion {
        ids,
        gm,
        gds,
        gmbs,
        mode,
        vdsat,
        gbs_jct,
        gbd_jct,
        iii,
        gii_d,
        gii_g,
        gii_b,
        ceq_d,
        ceq_bs,
        ceq_bd,
        ceq_iii,
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
        vbs_fd: vbseff,
    }
}

/// Stamp BSIM3SOI-FD companion model into the MNA matrix and RHS.
pub fn stamp_bsim3soi_fd(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &Bsim3SoiFdInstance,
    comp: &Bsim3SoiFdCompanion,
) {
    let dp = inst.drain_eff_idx();
    let g = inst.gate_idx;
    let sp = inst.source_eff_idx();
    let b = inst.body_int_idx;

    let sign = inst.model.mos_type.sign();
    let m = inst.m;

    let (xnrm, xrev) = if comp.mode > 0 {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };

    let gm_eff = m * (comp.gm * xnrm + comp.gds * xrev);
    let gds_eff = m * (comp.gds * xnrm + comp.gm * xrev);
    let gmbs_eff = m * comp.gmbs;
    let gbd = m * comp.gbd_jct;
    let gbs = m * comp.gbs_jct;

    // Channel current: asymmetric VCCS stamps (must use matrix.add, not stamp_conductance)
    if let Some(d) = dp {
        matrix.add(d, d, gds_eff);
        if let Some(gate) = g {
            matrix.add(d, gate, gm_eff);
        }
        if let Some(s) = sp {
            matrix.add(d, s, -(gm_eff + gds_eff + gmbs_eff));
        }
        if let Some(bulk) = b {
            matrix.add(d, bulk, gmbs_eff);
        }
    }
    if let Some(s) = sp {
        if let Some(d) = dp {
            matrix.add(s, d, -gds_eff);
        }
        if let Some(gate) = g {
            matrix.add(s, gate, -gm_eff);
        }
        matrix.add(s, s, gm_eff + gds_eff + gmbs_eff);
        if let Some(bulk) = b {
            matrix.add(s, bulk, -gmbs_eff);
        }
    }

    // BD/BS junctions: symmetric two-terminal stamps (correct use of stamp_conductance)
    crate::stamp_conductance(matrix, b, dp, gbd);
    crate::stamp_conductance(matrix, b, sp, gbs);

    // Impact ionization: Iii flows drain→body (OUT of drain, INTO body).
    if comp.iii != 0.0 {
        let gii_d = m * comp.gii_d;
        let gii_g = m * comp.gii_g;
        let gii_b = m * comp.gii_b;
        if let (Some(d), Some(bi)) = (dp, b) {
            matrix.add(d, d, gii_d);
            matrix.add(bi, d, -gii_d);
        }
        if let Some(gate) = g {
            if let Some(d) = dp {
                matrix.add(d, gate, gii_g);
            }
            if let Some(bi) = b {
                matrix.add(bi, gate, -gii_g);
            }
        }
        if let Some(bi) = b {
            matrix.add(bi, bi, -gii_b);
            if let Some(d) = dp {
                matrix.add(d, bi, gii_b);
            }
        }
    }

    // RHS
    let ceq_d = sign * m * comp.ceq_d;
    let ceq_iii = sign * m * comp.ceq_iii;
    if let Some(d) = dp {
        rhs[d] -= ceq_d + ceq_iii;
    }
    if let Some(s) = sp {
        rhs[s] += ceq_d;
    }
    // Impact ionization current enters body
    if let Some(bi) = b
        && comp.iii != 0.0
    {
        rhs[bi] += ceq_iii;
    }

    // Body resistance to external body contact
    if let (Some(b_int), Some(b_ext)) = (inst.body_int_idx, inst.body_idx) {
        let gbody = if inst.model.rbody > 0.0 {
            m / inst.model.rbody
        } else {
            m * 1e3
        };
        crate::stamp_conductance(matrix, Some(b_int), Some(b_ext), gbody);
    }

    // Series resistance: single stamp_conductance per resistor
    if inst.drain_prime_idx.is_some() && inst.drain_prime_idx != inst.drain_idx {
        let rd = if inst.nrd > 0.0 && inst.model.rbsh > 0.0 {
            inst.model.rbsh * inst.nrd
        } else {
            0.01
        };
        crate::stamp_conductance(matrix, inst.drain_idx, inst.drain_prime_idx, m / rd);
    }
    if inst.source_prime_idx.is_some() && inst.source_prime_idx != inst.source_idx {
        let rs = if inst.nrs > 0.0 && inst.model.rbsh > 0.0 {
            inst.model.rbsh * inst.nrs
        } else {
            0.01
        };
        crate::stamp_conductance(matrix, inst.source_idx, inst.source_prime_idx, m / rs);
    }
}

/// BSIM3SOI-FD voltage limiting for NR convergence.
///
/// FD uses simpler limiting than PD: all voltages limited ±3V per iteration
/// (matching b3soifdld.c B3SOIFDlimit with limit=3.0), and in floating body DC,
/// Vbs is clamped to non-negative.
pub fn bsim3soi_fd_limit(
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
    // FD body voltage limiting: ±0.2V (matching b3soifdld.c SmartVbs + B3SOIFDlimit)
    let limit = 0.2;
    let vbs = if (vbs_new - vbs_old).abs() > limit {
        if vbs_new > vbs_old {
            vbs_old + limit
        } else {
            vbs_old - limit
        }
    } else {
        vbs_new
    };
    // FD SmartVbs: in DC floating body, Vbs >= 0
    let vbs = vbs.max(0.0);

    let limit_e = 3.0;
    let ves = if (ves_new - ves_old).abs() > limit_e {
        if ves_new > ves_old {
            ves_old + limit_e
        } else {
            ves_old - limit_e
        }
    } else {
        ves_new
    };
    (vgs, vds, vbs, ves)
}
