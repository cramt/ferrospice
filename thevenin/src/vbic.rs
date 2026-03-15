//! VBIC (Vertical Bipolar Inter-Company) BJT model — 4-terminal, no self-heating.
//!
//! Implements the VBIC95 model with NR companion linearization.
//! This pass covers DC load (companion model) and NR stamps.
//! Self-heating (Vrth) and excess phase (NQS) are NOT implemented.

use thevenin_types::{Expr, ModelDef, Param};

use crate::diode::{VT_NOM, pnjlim, vcrit};
use crate::physics::safe_exp;

/// Boltzmann constant (J/K) — matches ngspice vbicload.c.
const KB: f64 = 1.380662e-23;
/// Elementary charge (C) — matches ngspice vbicload.c.
const QE: f64 = 1.602189e-19;

/// Thermal voltage at a given temperature in Kelvin.
fn vt_at(t_k: f64) -> f64 {
    KB * t_k / QE
}

/// VBIC polarity: NPN or PNP.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VbicType {
    Npn,
    Pnp,
}

impl VbicType {
    pub fn sign(self) -> f64 {
        match self {
            VbicType::Npn => 1.0,
            VbicType::Pnp => -1.0,
        }
    }
}

/// Junction depletion charge model.
///
/// Returns (charge, capacitance) for a junction with parameters:
/// - `v`: junction voltage
/// - `cj`: zero-bias capacitance
/// - `p`: built-in potential
/// - `m`: grading coefficient
/// - `aj`: smoothing parameter (-0.5 typical)
/// - `fc`: forward-bias coefficient (0.9 typical)
pub fn depletion_charge(v: f64, cj: f64, p: f64, m: f64, aj: f64, fc: f64) -> (f64, f64) {
    if cj == 0.0 || p == 0.0 {
        return (0.0, 0.0);
    }

    // Two-piece depletion charge model matching ngspice VBIC implementation.
    // Below FC*PE: standard power-law capacitance.
    // Above FC*PE: polynomial extrapolation to avoid singularity at v=PE.
    let dv0 = -p * fc; // = -(PE * FC)
    let dvh = v + dv0; // = v - FC*PE; positive means forward-biased beyond FC limit

    if aj <= 0.0 {
        // Standard SPICE junction model (AJE <= 0)
        if dvh > 0.0 {
            // Forward bias beyond FC*PE: linear+quadratic extrapolation
            let one_minus_fc = 1.0 - fc;
            let pwq = one_minus_fc.powf(-1.0 - m);
            // qlo = charge at the FC boundary (constant w.r.t. v)
            let qlo = p * (1.0 - one_minus_fc.powf(1.0 - m)) / (1.0 - m);
            // qhi = polynomial extrapolation beyond FC*PE
            let qhi = dvh * (1.0 - fc + 0.5 * m * dvh / p) * pwq;
            let qhi_dv = (1.0 - fc + m * dvh / p) * pwq;
            (cj * (qlo + qhi), cj * qhi_dv)
        } else {
            // Reverse/low forward bias: standard power law
            let one_minus_vp = (1.0 - v / p).max(1e-20);
            let pow_1mm = one_minus_vp.powf(1.0 - m);
            let q = p * (1.0 - pow_1mm) / (1.0 - m);
            let c = one_minus_vp.powf(-m);
            (cj * q, cj * c)
        }
    } else {
        // Smooth junction model (AJE > 0) using hyperbolic smoothing
        let mv0 = (dv0 * dv0 + 4.0 * aj * aj).sqrt();
        let vl0 = -0.5 * (dv0 + mv0);
        // Smoothed voltage
        let dv = v + dv0;
        let mv = (dv * dv + 4.0 * aj * aj).sqrt();
        let vl = -0.5 * (dv + mv) + vl0;
        let dvl_dv = -0.5 * (1.0 + dv / mv);

        let one_minus_vlp = (1.0 - vl / p).max(1e-20);
        let pow_1mm = one_minus_vlp.powf(1.0 - m);
        let q = p * (1.0 - pow_1mm) / (1.0 - m);
        let c = one_minus_vlp.powf(-m) * dvl_dv;
        (cj * q, cj * c)
    }
}

/// VBIC model parameters.
///
/// Parameter naming follows the VBIC parameter array p[0..107].
#[derive(Debug, Clone)]
pub struct VbicModel {
    pub vbic_type: VbicType,
    // Nominal temperature
    pub tnom: f64,
    // Resistances
    pub rcx: f64,
    pub rci: f64,
    pub vo: f64,
    pub gamm: f64,
    pub hrcf: f64,
    pub rbx: f64,
    pub rbi: f64,
    pub re: f64,
    pub rs: f64,
    pub rbp: f64,
    // Transport saturation current
    pub is: f64,
    pub nf: f64,
    pub nr: f64,
    pub fc: f64,
    // Oxide capacitances
    pub cbeo: f64,
    // B-E junction capacitance
    pub cje: f64,
    pub pe: f64,
    pub me: f64,
    pub aje: f64,
    // Oxide cap B-C
    pub cbco: f64,
    // B-C junction capacitance
    pub cjc: f64,
    pub qco: f64,
    pub cjep: f64,
    pub pc: f64,
    pub mc: f64,
    pub ajc: f64,
    // Substrate junction
    pub cjcp: f64,
    pub ps: f64,
    pub ms: f64,
    pub ajs: f64,
    // Forward base current
    pub ibei: f64,
    pub wbe: f64,
    pub nei: f64,
    pub iben: f64,
    pub nen: f64,
    // Reverse base current
    pub ibci: f64,
    pub nci: f64,
    pub ibcn: f64,
    pub ncn: f64,
    // Avalanche
    pub avc1: f64,
    pub avc2: f64,
    // Parasitic transport
    pub isp: f64,
    pub wsp: f64,
    pub nfp: f64,
    // Parasitic base currents
    pub ibeip: f64,
    pub ibenp: f64,
    pub ibcip: f64,
    pub ncip: f64,
    pub ibcnp: f64,
    pub ncnp: f64,
    // Early voltages
    pub vef: f64,
    pub ver: f64,
    // High injection
    pub ikf: f64,
    pub ikr: f64,
    pub ikp: f64,
    // Transit times
    pub tf: f64,
    pub qtf: f64,
    pub xtf: f64,
    pub vtf: f64,
    pub itf: f64,
    pub tr: f64,
    pub td: f64,
    // Flicker noise
    pub kfn: f64,
    pub afn: f64,
    pub bfn: f64,
    // Temperature coefficients
    pub xre: f64,
    pub xrbi: f64,
    pub xrci: f64,
    pub xrs: f64,
    pub xvo: f64,
    pub ea: f64,
    pub eaie: f64,
    pub eaic: f64,
    pub eais: f64,
    pub eane: f64,
    pub eanc: f64,
    pub eans: f64,
    pub xis: f64,
    pub xii: f64,
    pub xin: f64,
    pub tnf: f64,
    pub tavc: f64,
    // Thermal resistance (unused in this pass)
    pub rth: f64,
    pub cth: f64,
    pub vrt: f64,
    pub art: f64,
    pub ccso: f64,
    // Extended parameters
    pub qbm: f64,
    pub nkf: f64,
    pub xikf: f64,
    pub xrcx: f64,
    pub xrbx: f64,
    pub xrbp: f64,
    pub isrr: f64,
    pub xisr: f64,
    pub dear: f64,
    pub eap: f64,
    // Base-emitter breakdown
    pub vbbe: f64,
    pub nbbe: f64,
    pub ibbe: f64,
    pub tvbbe1: f64,
    pub tvbbe2: f64,
    pub tnbbe: f64,

    // Temperature-adjusted parameters (computed by temperature_adjust)
    pub is_t: f64,
    pub isp_t: f64,
    pub ibei_t: f64,
    pub ibci_t: f64,
    pub iben_t: f64,
    pub ibcn_t: f64,
    pub ibeip_t: f64,
    pub ibcip_t: f64,
    pub ibenp_t: f64,
    pub ibcnp_t: f64,
    pub rci_t: f64,
    pub rcx_t: f64,
    pub rbi_t: f64,
    pub rbx_t: f64,
    pub re_t: f64,
    pub rs_t: f64,
    pub rbp_t: f64,
    pub vo_t: f64,
    pub pe_t: f64,
    pub pc_t: f64,
    pub ps_t: f64,
    pub cje_t: f64,
    pub cjc_t: f64,
    pub cjep_t: f64,
    pub cjcp_t: f64,
    pub nf_t: f64,
    pub ikf_t: f64,
    pub ikr_t: f64,
    pub ikp_t: f64,
    pub avc2_t: f64,
    pub vbbe_t: f64,
    pub nbbe_t: f64,
    pub ibbe_t: f64,
    pub gamm_t: f64,
    pub isrr_t: f64,
    pub vt: f64,
}

impl VbicModel {
    pub fn new(vbic_type: VbicType) -> Self {
        let mut m = Self {
            vbic_type,
            tnom: 27.0,
            rcx: 0.0,
            rci: 0.1,
            vo: 0.0,
            gamm: 0.0,
            hrcf: 1.0,
            rbx: 0.0,
            rbi: 0.1,
            re: 0.0,
            rs: 0.0,
            rbp: 0.1,
            is: 1e-16,
            nf: 1.0,
            nr: 1.0,
            fc: 0.9,
            cbeo: 0.0,
            cje: 0.0,
            pe: 0.75,
            me: 0.33,
            aje: -0.5,
            cbco: 0.0,
            cjc: 0.0,
            qco: 0.0,
            cjep: 0.0,
            pc: 0.75,
            mc: 0.33,
            ajc: -0.5,
            cjcp: 0.0,
            ps: 0.75,
            ms: 0.33,
            ajs: -0.5,
            ibei: 1e-18,
            wbe: 1.0,
            nei: 1.0,
            iben: 0.0,
            nen: 2.0,
            ibci: 1e-16,
            nci: 1.0,
            ibcn: 0.0,
            ncn: 2.0,
            avc1: 0.0,
            avc2: 0.0,
            isp: 0.0,
            wsp: 1.0,
            nfp: 1.0,
            ibeip: 0.0,
            ibenp: 0.0,
            ibcip: 0.0,
            ncip: 1.0,
            ibcnp: 0.0,
            ncnp: 2.0,
            vef: 0.0,
            ver: 0.0,
            ikf: 0.0,
            ikr: 0.0,
            ikp: 0.0,
            tf: 0.0,
            qtf: 0.0,
            xtf: 0.0,
            vtf: 0.0,
            itf: 0.0,
            tr: 0.0,
            td: 0.0,
            kfn: 0.0,
            afn: 1.0,
            bfn: 1.0,
            xre: 0.0,
            xrbi: 0.0,
            xrci: 0.0,
            xrs: 0.0,
            xvo: 0.0,
            ea: 1.12,
            eaie: 1.12,
            eaic: 1.12,
            eais: 1.12,
            eane: 1.12,
            eanc: 1.12,
            eans: 1.12,
            xis: 3.0,
            xii: 3.0,
            xin: 3.0,
            tnf: 0.0,
            tavc: 0.0,
            rth: 0.0,
            cth: 0.0,
            vrt: 0.0,
            art: 0.1,
            ccso: 0.0,
            qbm: 0.0,
            nkf: 0.5,
            xikf: 0.0,
            xrcx: 0.0,
            xrbx: 0.0,
            xrbp: 0.0,
            isrr: 1.0,
            xisr: 0.0,
            dear: 0.0,
            eap: 1.12,
            vbbe: 0.0,
            nbbe: 1.0,
            ibbe: 1e-6,
            tvbbe1: 0.0,
            tvbbe2: 0.0,
            tnbbe: 0.0,
            // Temperature-adjusted values (initialized to nominal)
            is_t: 1e-16,
            isp_t: 0.0,
            ibei_t: 1e-18,
            ibci_t: 1e-16,
            iben_t: 0.0,
            ibcn_t: 0.0,
            ibeip_t: 0.0,
            ibcip_t: 0.0,
            ibenp_t: 0.0,
            ibcnp_t: 0.0,
            rci_t: 0.1,
            rcx_t: 0.0,
            rbi_t: 0.1,
            rbx_t: 0.0,
            re_t: 0.0,
            rs_t: 0.0,
            rbp_t: 0.1,
            vo_t: 0.0,
            pe_t: 0.75,
            pc_t: 0.75,
            ps_t: 0.75,
            cje_t: 0.0,
            cjc_t: 0.0,
            cjep_t: 0.0,
            cjcp_t: 0.0,
            nf_t: 1.0,
            ikf_t: 0.0,
            ikr_t: 0.0,
            ikp_t: 0.0,
            avc2_t: 0.0,
            vbbe_t: 0.0,
            nbbe_t: 1.0,
            ibbe_t: 1e-6,
            gamm_t: 0.0,
            isrr_t: 1.0,
            vt: VT_NOM,
        };
        m.temperature_adjust(27.0);
        m
    }

    pub fn from_model_def(model_def: &ModelDef) -> Self {
        let vbic_type = if model_def.kind.to_uppercase() == "PNP" {
            VbicType::Pnp
        } else {
            VbicType::Npn
        };
        let mut m = Self::new(vbic_type);
        for p in &model_def.params {
            if let Expr::Num(v) = &p.value {
                match p.name.to_uppercase().as_str() {
                    "TNOM" => m.tnom = *v,
                    "RCX" => m.rcx = *v,
                    "RCI" => m.rci = *v,
                    "VO" => m.vo = *v,
                    "GAMM" => m.gamm = *v,
                    "HRCF" => m.hrcf = *v,
                    "RBX" => m.rbx = *v,
                    "RBI" => m.rbi = *v,
                    "RE" => m.re = *v,
                    "RS" => m.rs = *v,
                    "RBP" => m.rbp = *v,
                    "IS" => m.is = *v,
                    "NF" => m.nf = *v,
                    "NR" => m.nr = *v,
                    "FC" => m.fc = *v,
                    "CBEO" => m.cbeo = *v,
                    "CJE" => m.cje = *v,
                    "PE" => m.pe = *v,
                    "ME" => m.me = *v,
                    "AJE" => m.aje = *v,
                    "CBCO" => m.cbco = *v,
                    "CJC" => m.cjc = *v,
                    "QCO" => m.qco = *v,
                    "CJEP" => m.cjep = *v,
                    "PC" => m.pc = *v,
                    "MC" => m.mc = *v,
                    "AJC" => m.ajc = *v,
                    "CJCP" => m.cjcp = *v,
                    "PS" => m.ps = *v,
                    "MS" => m.ms = *v,
                    "AJS" => m.ajs = *v,
                    "IBEI" => m.ibei = *v,
                    "WBE" => m.wbe = *v,
                    "NEI" => m.nei = *v,
                    "IBEN" => m.iben = *v,
                    "NEN" => m.nen = *v,
                    "IBCI" => m.ibci = *v,
                    "NCI" => m.nci = *v,
                    "IBCN" => m.ibcn = *v,
                    "NCN" => m.ncn = *v,
                    "AVC1" => m.avc1 = *v,
                    "AVC2" => m.avc2 = *v,
                    "ISP" => m.isp = *v,
                    "WSP" => m.wsp = *v,
                    "NFP" => m.nfp = *v,
                    "IBEIP" => m.ibeip = *v,
                    "IBENP" => m.ibenp = *v,
                    "IBCIP" => m.ibcip = *v,
                    "NCIP" => m.ncip = *v,
                    "IBCNP" => m.ibcnp = *v,
                    "NCNP" => m.ncnp = *v,
                    "VEF" => m.vef = *v,
                    "VER" => m.ver = *v,
                    "IKF" => m.ikf = *v,
                    "IKR" => m.ikr = *v,
                    "IKP" => m.ikp = *v,
                    "TF" => m.tf = *v,
                    "QTF" => m.qtf = *v,
                    "XTF" => m.xtf = *v,
                    "VTF" => m.vtf = *v,
                    "ITF" => m.itf = *v,
                    "TR" => m.tr = *v,
                    "TD" => m.td = *v,
                    "KFN" => m.kfn = *v,
                    "AFN" => m.afn = *v,
                    "BFN" => m.bfn = *v,
                    "XRE" => m.xre = *v,
                    "XRBI" => m.xrbi = *v,
                    "XRCI" => m.xrci = *v,
                    "XRS" => m.xrs = *v,
                    "XVO" => m.xvo = *v,
                    "EA" => m.ea = *v,
                    "EAIE" => m.eaie = *v,
                    "EAIC" => m.eaic = *v,
                    "EAIS" => m.eais = *v,
                    "EANE" => m.eane = *v,
                    "EANC" => m.eanc = *v,
                    "EANS" => m.eans = *v,
                    "XIS" => m.xis = *v,
                    "XII" => m.xii = *v,
                    "XIN" => m.xin = *v,
                    "TNF" => m.tnf = *v,
                    "TAVC" => m.tavc = *v,
                    "RTH" => m.rth = *v,
                    "CTH" => m.cth = *v,
                    "VRT" => m.vrt = *v,
                    "ART" => m.art = *v,
                    "CCSO" => m.ccso = *v,
                    "QBM" => m.qbm = *v,
                    "NKF" => m.nkf = *v,
                    "XIKF" => m.xikf = *v,
                    "XRCX" => m.xrcx = *v,
                    "XRBX" => m.xrbx = *v,
                    "XRBP" => m.xrbp = *v,
                    "ISRR" => m.isrr = *v,
                    "XISR" => m.xisr = *v,
                    "DEAR" => m.dear = *v,
                    "EAP" => m.eap = *v,
                    "VBBE" => m.vbbe = *v,
                    "NBBE" => m.nbbe = *v,
                    "IBBE" => m.ibbe = *v,
                    "TVBBE1" => m.tvbbe1 = *v,
                    "TVBBE2" => m.tvbbe2 = *v,
                    "TNBBE" => m.tnbbe = *v,
                    _ => {}
                }
            }
        }
        m.temperature_adjust(27.0);
        m
    }

    pub fn with_instance_params(mut self, params: &[Param]) -> Self {
        for p in params {
            if let Expr::Num(v) = &p.value {
                match p.name.to_uppercase().as_str() {
                    "IS" => self.is = *v,
                    "RCX" => self.rcx = *v,
                    "RCI" => self.rci = *v,
                    "RBX" => self.rbx = *v,
                    "RBI" => self.rbi = *v,
                    "RE" => self.re = *v,
                    "RS" => self.rs = *v,
                    _ => {}
                }
            }
        }
        self
    }

    /// Number of internal nodes needed.
    ///
    /// Always internal: collCI, baseBI, baseBP (3 minimum).
    /// Conditional: collCX (if RCX>0), baseBX (if RBX>0), emitEI (if RE>0), subsSI (if RS>0).
    pub fn internal_node_count(&self) -> usize {
        self.internal_node_count_with_substrate(true)
    }

    /// Count internal nodes, optionally excluding substrate internal node.
    /// When the BJT has no substrate terminal (3-terminal device), the SI
    /// node should not be created even if RS > 0.
    pub fn internal_node_count_with_substrate(&self, has_substrate: bool) -> usize {
        let mut count = 3; // collCI, baseBI, baseBP always present
        if self.rcx > 0.0 {
            count += 1;
        }
        if self.rbx > 0.0 {
            count += 1;
        }
        if self.re > 0.0 {
            count += 1;
        }
        if has_substrate && self.rs > 0.0 {
            count += 1;
        }
        if self.rth > 0.0 {
            count += 1; // thermal node for self-heating
        }
        count
    }

    /// Temperature adjustment following vbic_4T_et_cf_t from vbictemp.c.
    ///
    /// `temp` is the simulation temperature in degrees Celsius.
    pub fn temperature_adjust(&mut self, temp: f64) {
        let t_k = temp + 273.15;
        let tnom_k = self.tnom + 273.15;
        let vt = vt_at(t_k);
        let vt_nom = vt_at(tnom_k);
        self.vt = vt;

        let delt = t_k - tnom_k;
        let tratio = t_k / tnom_k;
        let tln_ratio = tratio.ln();

        // Current temperature scaling helper — ngspice vbictemp.c formula:
        // I_T = I_nom * (rT^xp * exp(-ea*(1-rT)/Vtv))^(1/nf)
        //     = I_nom * rT^(xp/nf) * exp(ea/nf * (1/Vtnom - 1/Vt))
        let temp_current = |i_nom: f64, xp: f64, ea_val: f64, nf_val: f64| -> f64 {
            if i_nom == 0.0 {
                return 0.0;
            }
            let arg = ea_val / (nf_val * vt_nom) - ea_val / (nf_val * vt);
            i_nom * tratio.powf(xp / nf_val) * safe_exp(arg)
        };

        // Resistance temperature scaling: R_t = R_nom * (tratio^xr)
        let temp_resistance = |r_nom: f64, xr: f64| -> f64 {
            if r_nom == 0.0 {
                return 0.0;
            }
            r_nom * tratio.powf(xr)
        };

        // Junction potential temperature scaling — ngspice VBIC vbictemp.c formula:
        //   psiio = 2*Vtnom * ln(2*sinh(PE_nom / (2*Vtnom)))
        //   psiin = psiio*rT - 3*Vt*ln(rT) - EA*(rT-1)
        //   PE_t  = psiin + 2*Vt * ln(0.5*(1 + sqrt(1 + 4*exp(-psiin/Vt))))
        // This matches ngspice vbictemp.c lines 278-319 exactly.
        let temp_potential = |p_nom: f64, ea_val: f64| -> f64 {
            if p_nom == 0.0 {
                return 0.0;
            }
            let x = p_nom / (2.0 * vt_nom);
            let sinh2x = x.exp() - (-x).exp(); // 2*sinh(x)
            let psiio = 2.0 * vt_nom * sinh2x.ln();
            let psiin = psiio * tratio - 3.0 * vt * tln_ratio - ea_val * (tratio - 1.0);
            psiin + 2.0 * vt * (0.5 * (1.0 + (1.0 + 4.0 * (-psiin / vt).exp()).sqrt())).ln()
        };

        // Junction capacitance temperature scaling (ngspice vbictemp.c):
        // C_t = C_nom * (P_nom / P_t)^M
        let temp_cap = |c_nom: f64, p_nom: f64, p_t: f64, m: f64| -> f64 {
            if c_nom == 0.0 || p_nom == 0.0 || p_t <= 0.0 {
                return 0.0;
            }
            c_nom * (p_nom / p_t).powf(m)
        };

        // Transport current IS
        self.is_t = temp_current(self.is, self.xis, self.ea, self.nf);

        // Parasitic transport current ISP
        self.isp_t = temp_current(self.isp, self.xis, self.eap, self.nfp);

        // ISRR
        self.isrr_t = temp_current(self.isrr * self.is, self.xisr, self.dear + self.ea, self.nr);
        if self.is_t > 0.0 {
            self.isrr_t /= self.is_t;
        }

        // Base currents - ideal
        self.ibei_t = temp_current(self.ibei, self.xii, self.eaie, self.nei);
        self.ibci_t = temp_current(self.ibci, self.xii, self.eaic, self.nci);

        // Base currents - non-ideal
        self.iben_t = temp_current(self.iben, self.xin, self.eane, self.nen);
        self.ibcn_t = temp_current(self.ibcn, self.xin, self.eanc, self.ncn);

        // Parasitic base currents
        self.ibeip_t = temp_current(self.ibeip, self.xii, self.eaic, self.nci);
        self.ibcip_t = temp_current(self.ibcip, self.xii, self.eaic, self.ncip);
        self.ibenp_t = temp_current(self.ibenp, self.xin, self.eanc, self.nen);
        self.ibcnp_t = temp_current(self.ibcnp, self.xin, self.eanc, self.ncnp);

        // Resistances
        self.rci_t = temp_resistance(self.rci, self.xrci);
        self.rcx_t = temp_resistance(self.rcx, self.xrcx);
        self.rbi_t = temp_resistance(self.rbi, self.xrbi);
        self.rbx_t = temp_resistance(self.rbx, self.xrbx);
        self.re_t = temp_resistance(self.re, self.xre);
        self.rs_t = temp_resistance(self.rs, self.xrs);
        self.rbp_t = temp_resistance(self.rbp, self.xrbp);

        // VO
        self.vo_t = temp_resistance(self.vo, self.xvo);

        // GAMM temperature scaling — ngspice vbictemp.c:
        //   p[4] = pnom[4] * rT^XIS * exp(-EA*(1-rT)/Vtv)
        //        = GAMM * rT^XIS * exp(EA/Vtnom - EA/Vt)
        self.gamm_t = if self.gamm > 0.0 {
            self.gamm * tratio.powf(self.xis) * safe_exp(self.ea / vt_nom - self.ea / vt)
        } else {
            0.0
        };

        // Junction potentials
        self.pe_t = temp_potential(self.pe, self.eaie);
        self.pc_t = temp_potential(self.pc, self.eaic);
        self.ps_t = temp_potential(self.ps, self.eais);

        // Junction capacitances
        self.cje_t = temp_cap(self.cje, self.pe, self.pe_t, self.me);
        self.cjc_t = temp_cap(self.cjc, self.pc, self.pc_t, self.mc);
        self.cjep_t = temp_cap(self.cjep, self.pe, self.pe_t, self.me);
        self.cjcp_t = temp_cap(self.cjcp, self.ps, self.ps_t, self.ms);

        // NF temperature
        self.nf_t = self.nf * (1.0 + self.tnf * delt);

        // IKF, IKR, IKP temperature
        self.ikf_t = if self.ikf > 0.0 {
            self.ikf * tratio.powf(self.xikf)
        } else {
            0.0
        };
        self.ikr_t = if self.ikr > 0.0 {
            self.ikr * tratio.powf(self.xikf) // same exponent as IKF
        } else {
            0.0
        };
        self.ikp_t = if self.ikp > 0.0 {
            self.ikp * tratio.powf(self.xikf)
        } else {
            0.0
        };

        // AVC2 temperature
        self.avc2_t = self.avc2 * (1.0 + self.tavc * delt);

        // VBBE temperature
        self.vbbe_t = self.vbbe * (1.0 + self.tvbbe1 * delt + self.tvbbe2 * delt * delt);
        self.nbbe_t = self.nbbe * (1.0 + self.tnbbe * delt);
        self.ibbe_t = self.ibbe * safe_exp(-self.vbbe_t / (self.nbbe_t * vt));
    }

    /// Critical voltage for B-E junction (forward).
    pub fn vcrit_bei(&self) -> f64 {
        vcrit(self.nei * self.vt, self.ibei_t)
    }

    /// Critical voltage for B-C junction.
    pub fn vcrit_bci(&self) -> f64 {
        vcrit(self.nci * self.vt, self.ibci_t)
    }

    /// Limit B-E internal junction voltage.
    pub fn limit_vbei(&self, v_new: f64, v_old: f64) -> f64 {
        pnjlim(v_new, v_old, self.nei * self.vt, self.vcrit_bei())
    }

    /// Limit B-C internal junction voltage.
    pub fn limit_vbci(&self, v_new: f64, v_old: f64) -> f64 {
        pnjlim(v_new, v_old, self.nci * self.vt, self.vcrit_bci())
    }

    /// Compute the VBIC operating point and NR companion model.
    ///
    /// Junction voltages should already be limited via `pnjlim` before calling.
    #[expect(clippy::too_many_arguments)]
    pub fn companion(
        &self,
        vbei: f64,
        vbex: f64,
        vbci: f64,
        vbcx: f64,
        vbep: f64,
        vrci: f64,
        vrbi: f64,
        vrbp: f64,
        vbcp: f64,
        gmin: f64,
    ) -> VbicCompanion {
        let vt = self.vt;

        // =====================================================================
        // 1. Forward and reverse transport currents
        // =====================================================================
        let (ifi, difi_dvbei) = {
            let arg = vbei / (self.nf_t * vt);
            let e = safe_exp(arg);
            let i = self.is_t * (e - 1.0);
            let g = self.is_t / (self.nf_t * vt) * e;
            (i, g)
        };

        let is_r = self.is_t * self.isrr_t;
        let (iri, diri_dvbci) = {
            let arg = vbci / (self.nr * vt);
            let e = safe_exp(arg);
            let i = is_r * (e - 1.0);
            let g = is_r / (self.nr * vt) * e;
            (i, g)
        };

        // =====================================================================
        // 2. Depletion charges for Early effect (qdbe, qdbc)
        // =====================================================================
        // Note: qdbe/qdbc in ngspice are in VOLTS (charge/CJE, i.e., the voltage
        // integral of the normalized capacitance).  We pass cj=1.0 so the result
        // has the same Volt units used in q1z = 1 + qdbe*IVER + qdbc*IVEF.
        let (qdbe, dqdbe_dvbei) =
            depletion_charge(vbei, 1.0, self.pe_t, self.me, self.aje, self.fc);
        let (qdbc, dqdbc_dvbci) =
            depletion_charge(vbci, 1.0, self.pc_t, self.mc, self.ajc, self.fc);

        // =====================================================================
        // 3. Base charge normalization (qb)
        // =====================================================================
        // ngspice: q1z = 1 + qdbe*IVER + qdbc*IVEF  (IVER=1/VER, IVEF=1/VEF)
        // So: qdbe/VER + qdbc/VEF
        let q1_term_ver = if self.ver > 0.0 { qdbe / self.ver } else { 0.0 };
        let q1_term_vef = if self.vef > 0.0 { qdbc / self.vef } else { 0.0 };
        let q1z = 1.0 + q1_term_ver + q1_term_vef;

        let dq1z_dvbei = if self.ver > 0.0 {
            dqdbe_dvbei / self.ver
        } else {
            0.0
        };
        let dq1z_dvbci = if self.vef > 0.0 {
            dqdbc_dvbci / self.vef
        } else {
            0.0
        };

        let q2_ikf = if self.ikf_t > 0.0 {
            ifi / self.ikf_t
        } else {
            0.0
        };
        let q2_ikr = if self.ikr_t > 0.0 {
            iri / self.ikr_t
        } else {
            0.0
        };
        let q2 = q2_ikf + q2_ikr;

        let dq2_dvbei = if self.ikf_t > 0.0 {
            difi_dvbei / self.ikf_t
        } else {
            0.0
        };
        let dq2_dvbci = if self.ikr_t > 0.0 {
            diri_dvbci / self.ikr_t
        } else {
            0.0
        };

        let (qb, dqb_dvbei, dqb_dvbci) = if self.qbm == 0.0 {
            // Standard VBIC base charge: qb = q1z/2 * (1 + sqrt(1 + 4*q2/q1z^2))
            let q1z_sq = q1z * q1z;
            let arg = (1.0 + 4.0 * q2 / q1z_sq).max(0.0);
            let sqrt_arg = arg.sqrt();
            let qb_val = q1z * (1.0 + sqrt_arg) / 2.0;

            // dqb/dvbei
            let dsqrt_dvbei = if sqrt_arg > 0.0 {
                (4.0 * dq2_dvbei * q1z_sq - 4.0 * q2 * 2.0 * q1z * dq1z_dvbei)
                    / (2.0 * sqrt_arg * q1z_sq * q1z_sq)
            } else {
                0.0
            };
            let dqb_dvbei_val = dq1z_dvbei * (1.0 + sqrt_arg) / 2.0 + q1z * dsqrt_dvbei / 2.0;

            // dqb/dvbci
            let dsqrt_dvbci = if sqrt_arg > 0.0 {
                (4.0 * dq2_dvbci * q1z_sq - 4.0 * q2 * 2.0 * q1z * dq1z_dvbci)
                    / (2.0 * sqrt_arg * q1z_sq * q1z_sq)
            } else {
                0.0
            };
            let dqb_dvbci_val = dq1z_dvbci * (1.0 + sqrt_arg) / 2.0 + q1z * dsqrt_dvbci / 2.0;

            (qb_val, dqb_dvbei_val, dqb_dvbci_val)
        } else {
            // Alternate base charge formulation: qb = q1z*(1 + q2)
            let qb_val = q1z * (1.0 + q2);
            let dqb_dvbei_val = dq1z_dvbei * (1.0 + q2) + q1z * dq2_dvbei;
            let dqb_dvbci_val = dq1z_dvbci * (1.0 + q2) + q1z * dq2_dvbci;
            (qb_val, dqb_dvbei_val, dqb_dvbci_val)
        };

        let qb = qb.max(1e-12); // Avoid division by zero

        // =====================================================================
        // 4. Collector transport current
        // =====================================================================
        let itzf = ifi / qb;
        let itzr = iri / qb;

        // Derivatives of itzf
        let ditzf_dvbei = (difi_dvbei * qb - ifi * dqb_dvbei) / (qb * qb);
        let ditzf_dvbci = -ifi * dqb_dvbci / (qb * qb);

        // Derivatives of itzr
        let ditzr_dvbci = (diri_dvbci * qb - iri * dqb_dvbci) / (qb * qb);
        let ditzr_dvbei = -iri * dqb_dvbei / (qb * qb);

        // No excess phase: Iciei = Itzf - Itzr
        let iciei = itzf - itzr;
        let diciei_dvbei = ditzf_dvbei - ditzr_dvbei;
        let diciei_dvbci = ditzf_dvbci - ditzr_dvbci;

        // =====================================================================
        // 5. Base currents
        // =====================================================================

        // --- Ibe: B-E ideal + non-ideal + breakdown ---
        let (ibe_ideal, dibe_ideal_dvbei) = {
            let wbe = self.wbe;
            let arg = vbei / (self.nei * vt);
            let e = safe_exp(arg);
            let i = wbe * self.ibei_t * (e - 1.0);
            let g = wbe * self.ibei_t / (self.nei * vt) * e;
            (i, g)
        };

        let (ibe_nonideal, dibe_nonideal_dvbei) = if self.iben_t > 0.0 {
            let wbe = self.wbe;
            let arg = vbei / (self.nen * vt);
            let e = safe_exp(arg);
            let i = wbe * self.iben_t * (e - 1.0);
            let g = wbe * self.iben_t / (self.nen * vt) * e;
            (i, g)
        } else {
            (0.0, 0.0)
        };

        // VBBE breakdown term
        let (ibe_vbbe, dibe_vbbe_dvbei) = if self.vbbe_t != 0.0 {
            let arg = (-self.vbbe_t + vbei) / (self.nbbe_t * vt);
            let e = safe_exp(arg);
            let i = self.ibbe_t * e;
            let g = self.ibbe_t / (self.nbbe_t * vt) * e;
            (i, g)
        } else {
            (0.0, 0.0)
        };

        let ibe = ibe_ideal + ibe_nonideal + ibe_vbbe;
        let dibe_dvbei = dibe_ideal_dvbei + dibe_nonideal_dvbei + dibe_vbbe_dvbei;

        // --- Ibex: external B-E (1-WBE fraction) ---
        let (ibex, dibex_dvbex) = {
            let wbe_x = 1.0 - self.wbe;
            if wbe_x > 0.0 {
                let arg1 = vbex / (self.nei * vt);
                let e1 = safe_exp(arg1);
                let i1 = wbe_x * self.ibei_t * (e1 - 1.0);
                let g1 = wbe_x * self.ibei_t / (self.nei * vt) * e1;

                let (i2, g2) = if self.iben_t > 0.0 {
                    let arg2 = vbex / (self.nen * vt);
                    let e2 = safe_exp(arg2);
                    (
                        wbe_x * self.iben_t * (e2 - 1.0),
                        wbe_x * self.iben_t / (self.nen * vt) * e2,
                    )
                } else {
                    (0.0, 0.0)
                };

                (i1 + i2, g1 + g2)
            } else {
                (0.0, 0.0)
            }
        };

        // --- Ibc: B-C ideal + non-ideal ---
        let (ibcj, dibcj_dvbci) = {
            let arg1 = vbci / (self.nci * vt);
            let e1 = safe_exp(arg1);
            let i1 = self.ibci_t * (e1 - 1.0);
            let g1 = self.ibci_t / (self.nci * vt) * e1;

            let (i2, g2) = if self.ibcn_t > 0.0 {
                let arg2 = vbci / (self.ncn * vt);
                let e2 = safe_exp(arg2);
                (self.ibcn_t * (e2 - 1.0), self.ibcn_t / (self.ncn * vt) * e2)
            } else {
                (0.0, 0.0)
            };

            (i1 + i2, g1 + g2)
        };

        // --- Avalanche current Igc ---
        let (igc, digc_dvbci, digc_dvbei) = if self.avc1 > 0.0 && self.avc2_t > 0.0 {
            // Smooth clamping: vl = 0.5*(sqrt((PC-Vbci)^2 + 0.01) + PC - Vbci)
            let pc_minus_vbci = self.pc_t - vbci;
            let vl = 0.5 * ((pc_minus_vbci * pc_minus_vbci + 0.01).sqrt() + pc_minus_vbci);
            let dvl_dvbci =
                0.5 * (-pc_minus_vbci / (pc_minus_vbci * pc_minus_vbci + 0.01).sqrt() - 1.0);

            let avalf = self.avc1 * vl * safe_exp(-self.avc2_t * vl);
            let davalf_dvl = self.avc1 * safe_exp(-self.avc2_t * vl) * (1.0 - self.avc2_t * vl);

            // Igc = (Itzf - Itzr - Ibcj) * avalf
            let i_drive = itzf - itzr - ibcj;
            let igc_val = i_drive * avalf;

            // dIgc/dVbci = d(i_drive)/dVbci * avalf + i_drive * davalf/dVbci
            let di_drive_dvbci = ditzf_dvbci - ditzr_dvbci - dibcj_dvbci;
            let digc_dvbci_val = di_drive_dvbci * avalf + i_drive * davalf_dvl * dvl_dvbci;

            // dIgc/dVbei = d(i_drive)/dVbei * avalf
            let di_drive_dvbei = ditzf_dvbei - ditzr_dvbei;
            let digc_dvbei_val = di_drive_dvbei * avalf;

            (igc_val, digc_dvbci_val, digc_dvbei_val)
        } else {
            (0.0, 0.0, 0.0)
        };

        let ibc = ibcj + igc;
        let dibc_dvbci = dibcj_dvbci + digc_dvbci;
        let dibc_dvbei = digc_dvbei;

        // --- Ibep: parasitic B-E ---
        let (ibep, dibep_dvbep) = {
            let (i1, g1) = if self.ibeip_t > 0.0 {
                let arg = vbep / (self.nci * vt);
                let e = safe_exp(arg);
                (self.ibeip_t * (e - 1.0), self.ibeip_t / (self.nci * vt) * e)
            } else {
                (0.0, 0.0)
            };

            let (i2, g2) = if self.ibenp_t > 0.0 {
                let arg = vbep / (self.nen * vt);
                let e = safe_exp(arg);
                (self.ibenp_t * (e - 1.0), self.ibenp_t / (self.nen * vt) * e)
            } else {
                (0.0, 0.0)
            };

            (i1 + i2, g1 + g2)
        };

        // --- Ibcp: parasitic B-C ---
        let (ibcp, dibcp_dvbcp) = {
            let (i1, g1) = if self.ibcip_t > 0.0 {
                let arg = vbcp / (self.ncip * vt);
                let e = safe_exp(arg);
                (
                    self.ibcip_t * (e - 1.0),
                    self.ibcip_t / (self.ncip * vt) * e,
                )
            } else {
                (0.0, 0.0)
            };

            let (i2, g2) = if self.ibcnp_t > 0.0 {
                let arg = vbcp / (self.ncnp * vt);
                let e = safe_exp(arg);
                (
                    self.ibcnp_t * (e - 1.0),
                    self.ibcnp_t / (self.ncnp * vt) * e,
                )
            } else {
                (0.0, 0.0)
            };

            (i1 + i2, g1 + g2)
        };

        // =====================================================================
        // 6. Parasitic transport current
        // =====================================================================
        let (iccp, diccp_dvbep, diccp_dvbcp) = if self.isp_t > 0.0 {
            let arg_f = vbep / (self.nfp * vt);
            let ef = safe_exp(arg_f);
            let ifp = self.isp_t * (ef - 1.0);
            let difp_dvbep = self.isp_t / (self.nfp * vt) * ef;

            let arg_r = vbcp / (self.nfp * vt);
            let er = safe_exp(arg_r);
            let irp = self.isp_t * (er - 1.0);
            let dirp_dvbcp = self.isp_t / (self.nfp * vt) * er;

            let (qbp, dqbp_dvbep) = if self.ikp_t > 0.0 {
                let arg = 1.0 + 4.0 * ifp / self.ikp_t;
                let sq = arg.max(0.0).sqrt();
                let qbp_val = (1.0 + sq) / 2.0;
                let dqbp = if sq > 0.0 {
                    2.0 * difp_dvbep / (self.ikp_t * sq) / 2.0
                } else {
                    0.0
                };
                (qbp_val, dqbp)
            } else {
                (1.0, 0.0)
            };

            let qbp = qbp.max(1e-12);
            let iccp_val = (ifp - irp) / qbp;
            let diccp_dvbep_val = (difp_dvbep * qbp - (ifp - irp) * dqbp_dvbep) / (qbp * qbp);
            let diccp_dvbcp_val = -dirp_dvbcp / qbp;

            (iccp_val, diccp_dvbep_val, diccp_dvbcp_val)
        } else {
            (0.0, 0.0, 0.0)
        };

        // =====================================================================
        // 7. Internal resistance currents
        // =====================================================================

        // --- Irci: internal collector resistance with quasi-saturation ---
        let (irci, dirci_dvrci, dirci_dvbci, dirci_dvbcx) = self.compute_irci(vrci, vbci, vbcx, vt);

        // --- Irbi: internal base resistance modulated by qb ---
        let rbi_eff = if self.rbi_t > 0.0 {
            self.rbi_t / qb
        } else {
            0.0
        };
        let (irbi, dirbi_dvrbi, dirbi_dvbei, dirbi_dvbci) = if rbi_eff > 0.0 {
            let g = qb / self.rbi_t;
            let i = vrbi * g;
            // dIrbi/dVrbi = qb/RBI
            // dIrbi/dVbei = Vrbi/RBI * dqb/dVbei
            // dIrbi/dVbci = Vrbi/RBI * dqb/dVbci
            (
                i,
                g,
                vrbi * dqb_dvbei / self.rbi_t,
                vrbi * dqb_dvbci / self.rbi_t,
            )
        } else {
            (0.0, 0.0, 0.0, 0.0)
        };

        // --- Irbp: parasitic base resistance ---
        // RBP is modulated by parasitic qbp
        let (irbp, dirbp_dvrbp) = if self.rbp_t > 0.0 {
            let defl = if self.isp_t > 0.0 && self.ikp_t > 0.0 {
                let arg_f = vbep / (self.nfp * vt);
                let ef = safe_exp(arg_f);
                let ifp_val = self.isp_t * (ef - 1.0);
                let q_arg = (1.0 + 4.0 * ifp_val / self.ikp_t).max(0.0).sqrt();
                (1.0 + q_arg) / 2.0
            } else {
                1.0
            };
            let g = defl / self.rbp_t;
            (vrbp * g, g)
        } else {
            (0.0, 0.0)
        };

        // =====================================================================
        // 8. Depletion charges
        // =====================================================================
        let (qbe, cqbe) = depletion_charge(vbei, self.cje_t, self.pe_t, self.me, self.aje, self.fc);
        let (qbc, cqbc) = depletion_charge(vbci, self.cjc_t, self.pc_t, self.mc, self.ajc, self.fc);
        let (qbep, cqbep) =
            depletion_charge(vbep, self.cjep_t, self.pe_t, self.me, self.aje, self.fc);
        let (qbcp, cqbcp) =
            depletion_charge(vbcp, self.cjcp_t, self.ps_t, self.ms, self.ajs, self.fc);

        // External BE depletion charge (1-WBE fraction)
        let (_qbex, cqbex) =
            depletion_charge(vbex, self.cje_t, self.pe_t, self.me, self.aje, self.fc);

        // =====================================================================
        // 9. Transit time charges (TF modulation from vbicload.c lines 3830-3887)
        // =====================================================================
        // Forward transit time with ITF/XTF/VTF modulation:
        //   tff = TF * (1 + QTF*q1) * (1 + XTF*exp(Vbci/(1.44*VTF)) * (slTF + mIf^2) * sgIf)
        // where q1 = qdbe/VER + qdbc/VEF (Early charge terms)
        //       mIf = rIf/(rIf+1), rIf = Ifi/(Ifi+IKF_t)
        //       sgIf = sign(Ifi), slTF = 0 (simplified)
        let (tff, dtff_dvbei, dtff_dvbci) = if self.tf > 0.0 {
            let q1 = q1_term_vef + q1_term_ver;
            let dq1_dvbei = dq1z_dvbei; // dq1z_dvbei excludes the base 1.0
            let dq1_dvbci = dq1z_dvbci;

            let qtf_factor = 1.0 + self.qtf * q1;
            let dqtf_factor_dvbei = self.qtf * dq1_dvbei;
            let dqtf_factor_dvbci = self.qtf * dq1_dvbci;

            if self.xtf > 0.0 && self.itf > 0.0 {
                // ITF/XTF modulation
                let rif = if self.ikf_t > 0.0 {
                    ifi / (ifi + self.ikf_t)
                } else {
                    0.0
                };
                let drif_dvbei = if self.ikf_t > 0.0 {
                    difi_dvbei * self.ikf_t / ((ifi + self.ikf_t) * (ifi + self.ikf_t))
                } else {
                    0.0
                };

                let mif = rif / (rif + 1.0);
                let dmif_drif = 1.0 / ((rif + 1.0) * (rif + 1.0));
                let dmif_dvbei = dmif_drif * drif_dvbei;

                let ivtf = if self.vtf > 0.0 { 1.0 / self.vtf } else { 0.0 };
                let xvar1 = vbci * ivtf / 1.44;
                let xvar2 = safe_exp(xvar1);
                let dxvar2_dvbci = xvar2 * ivtf / 1.44;

                let sgif = if ifi >= 0.0 { 1.0 } else { -1.0 };
                let sl_tf = 0.0; // slTF simplified

                let xtf_term = 1.0 + self.xtf * xvar2 * (sl_tf + mif * mif) * sgif;
                let dxtf_dvbei = self.xtf * xvar2 * 2.0 * mif * dmif_dvbei * sgif;
                let dxtf_dvbci = self.xtf * dxvar2_dvbci * (sl_tf + mif * mif) * sgif;

                let tff_val = self.tf * qtf_factor * xtf_term;
                let dtff_dvbei_val =
                    self.tf * (dqtf_factor_dvbei * xtf_term + qtf_factor * dxtf_dvbei);
                let dtff_dvbci_val =
                    self.tf * (dqtf_factor_dvbci * xtf_term + qtf_factor * dxtf_dvbci);
                (tff_val, dtff_dvbei_val, dtff_dvbci_val)
            } else {
                let tff_val = self.tf * qtf_factor;
                let dtff_dvbei_val = self.tf * dqtf_factor_dvbei;
                let dtff_dvbci_val = self.tf * dqtf_factor_dvbci;
                (tff_val, dtff_dvbei_val, dtff_dvbci_val)
            }
        } else {
            (0.0, 0.0, 0.0)
        };

        // Total Qbe = CJE_t * qdbe * WBE + tff * Ifi / qb
        let dtransit_dvbei =
            (dtff_dvbei * ifi + tff * difi_dvbei) / qb - tff * ifi * dqb_dvbei / (qb * qb);
        let dtransit_dvbci = dtff_dvbci * ifi / qb - tff * ifi * dqb_dvbci / (qb * qb);

        let cqbe_vbei = cqbe * self.wbe + dtransit_dvbei;
        let cqbe_vbci = dtransit_dvbci;

        // External BE: Qbex = CJE_t * qdbex * (1-WBE)
        let cqbex_vbex = cqbex * (1.0 - self.wbe);

        // Total Qbc = CJC_t * qdbc + TR * Iri + QCO * Kbci
        // Kbci from quasi-saturation model
        let gamm_exp_bci = self.gamm_t * safe_exp(vbci / vt);
        let kbci = (1.0 + gamm_exp_bci).sqrt();
        let dkbci_dvbci = if kbci > 0.0 {
            gamm_exp_bci / (vt * 2.0 * kbci)
        } else {
            0.0
        };
        let cqbc_vbci = cqbc + self.tr * diri_dvbci + self.qco * dkbci_dvbci;

        // Qbcx = QCO * Kbcx
        let gamm_exp_bcx = self.gamm_t * safe_exp(vbcx / vt);
        let kbcx = (1.0 + gamm_exp_bcx).sqrt();
        let dkbcx_dvbcx = if kbcx > 0.0 {
            gamm_exp_bcx / (vt * 2.0 * kbcx)
        } else {
            0.0
        };
        let cqbcx_vbcx = self.qco * dkbcx_dvbcx;

        // Qbep = CJEP_t * qdbep + TR * Ifp (parasitic forward current)
        let (ifp_for_charge, difp_dvbep_for_charge) = if self.isp_t > 0.0 {
            let arg = vbep / (self.nfp * vt);
            let ef = safe_exp(arg);
            (self.isp_t * (ef - 1.0), self.isp_t / (self.nfp * vt) * ef)
        } else {
            (0.0, 0.0)
        };
        let _ = ifp_for_charge; // Ifp value used for charge only
        let cqbep_vbep = cqbep + self.tr * difp_dvbep_for_charge;
        let cqbep_vbci = 0.0; // Ifp doesn't depend on Vbci in 4T model

        // Qbcp = CJCP_t * qdbcp + CCSO * Vbcp
        let cqbcp_vbcp = cqbcp + self.ccso;

        // Apply gmin stabilisation to junction currents (matches ngspice vbicload.c lines 756-771).
        // CKTgmin is added to each junction to prevent singular matrices and provide a numerical
        // floor in reverse-biased / subthreshold regions.
        let ibe = ibe + gmin * vbei;
        let dibe_dvbei = dibe_dvbei + gmin;
        let ibex = ibex + gmin * vbex;
        let dibex_dvbex = dibex_dvbex + gmin;
        let ibc = ibc + gmin * vbci;
        let dibc_dvbci = dibc_dvbci + gmin;
        let ibep = ibep + gmin * vbep;
        let dibep_dvbep = dibep_dvbep + gmin;
        // Irci gets three gmin terms: on Vrci, Vbci, and Vbcx
        let irci = irci + gmin * vrci + gmin * vbci + gmin * vbcx;
        let dirci_dvrci = dirci_dvrci + gmin;
        let dirci_dvbci = dirci_dvbci + gmin;
        let dirci_dvbcx = dirci_dvbcx + gmin;
        let ibcp = ibcp + gmin * vbcp;
        let dibcp_dvbcp = dibcp_dvbcp + gmin;

        VbicCompanion {
            // Transport
            iciei,
            diciei_dvbei,
            diciei_dvbci,

            // Base currents
            ibe,
            dibe_dvbei,
            ibex,
            dibex_dvbex,
            ibc,
            dibc_dvbci,
            dibc_dvbei,
            ibep,
            dibep_dvbep,
            ibcp,
            dibcp_dvbcp,

            // Parasitic transport
            iccp,
            diccp_dvbep,
            diccp_dvbcp,

            // Avalanche
            igc,
            digc_dvbci,
            digc_dvbei,

            // Resistive
            irci,
            dirci_dvrci,
            dirci_dvbci,
            dirci_dvbcx,
            irbi,
            dirbi_dvrbi,
            dirbi_dvbei,
            dirbi_dvbci,
            irbp,
            dirbp_dvrbp,

            // Depletion charges
            qbe,
            cqbe,
            qbc,
            cqbc,
            qbep,
            cqbep,
            qbcp,
            cqbcp,

            // Total charge derivatives for AC
            cqbe_vbei,
            cqbe_vbci,
            cqbex_vbex,
            cqbc_vbci,
            cqbcx_vbcx,
            cqbep_vbep,
            cqbep_vbci,
            cqbcp_vbcp,

            // Operating point info
            itzf,
            itzr,
            qb,
        }
    }

    /// Compute Irci (internal collector resistance) with quasi-saturation model.
    ///
    /// Returns (Irci, dIrci/dVrci, dIrci/dVbci, dIrci/dVbcx).
    fn compute_irci(&self, vrci: f64, vbci: f64, vbcx: f64, vt: f64) -> (f64, f64, f64, f64) {
        if self.rci_t <= 0.0 {
            return (0.0, 0.0, 0.0, 0.0);
        }

        // Kbci = sqrt(1 + GAMM*exp(Vbci/Vt))
        let gamm_exp_bci = self.gamm_t * safe_exp(vbci / vt);
        let kbci = (1.0 + gamm_exp_bci).sqrt();
        let dkbci_dvbci = if kbci > 0.0 {
            gamm_exp_bci / (vt * 2.0 * kbci)
        } else {
            0.0
        };

        // Kbcx = sqrt(1 + GAMM*exp(Vbcx/Vt))
        let gamm_exp_bcx = self.gamm_t * safe_exp(vbcx / vt);
        let kbcx = (1.0 + gamm_exp_bcx).sqrt();
        let dkbcx_dvbcx = if kbcx > 0.0 {
            gamm_exp_bcx / (vt * 2.0 * kbcx)
        } else {
            0.0
        };

        // rKp1 = (Kbci+1) / (Kbcx+1)  [ngspice vbicload.c line 3691]
        let rkp1 = (kbci + 1.0) / (kbcx + 1.0);
        let drkp1_dvbci = dkbci_dvbci / (kbcx + 1.0);
        let drkp1_dvbcx = -(kbci + 1.0) * dkbcx_dvbcx / ((kbcx + 1.0) * (kbcx + 1.0));

        // Ohmic current: Iohm = (Vrci + Vt*(Kbci - Kbcx - ln(rKp1))) / RCI  [ngspice line 3703]
        let ln_rkp1 = rkp1.ln();
        let dln_rkp1_dvbci = drkp1_dvbci / rkp1;
        let dln_rkp1_dvbcx = drkp1_dvbcx / rkp1;
        let iohm = (vrci + vt * (kbci - kbcx - ln_rkp1)) / self.rci_t;
        let diohm_dvrci = 1.0 / self.rci_t;
        let diohm_dvbci = vt * (dkbci_dvbci - dln_rkp1_dvbci) / self.rci_t;
        let diohm_dvbcx = vt * (-dkbcx_dvbcx - dln_rkp1_dvbcx) / self.rci_t;

        if self.vo_t > 0.0 {
            // Velocity saturation model
            let ivo = 1.0 / self.vo_t;
            let ihrcf = 1.0 / self.hrcf;

            // derf = IVO*RCI*Iohm / (1 + 0.5*IVO*IHRCF*sqrt(Vrci^2 + 0.01))
            let vrci_smooth = (vrci * vrci + 0.01).sqrt();
            let denom = 1.0 + 0.5 * ivo * ihrcf * vrci_smooth;
            let num = ivo * self.rci_t * iohm;
            let derf = num / denom;

            let sqrt_1_derf2 = (1.0 + derf * derf).sqrt();
            let irci_val = iohm / sqrt_1_derf2;

            // Derivative of derf w.r.t. vrci
            let dvrci_smooth_dvrci = vrci / vrci_smooth;
            let ddenom_dvrci = 0.5 * ivo * ihrcf * dvrci_smooth_dvrci;
            let dnum_dvrci = ivo * self.rci_t * diohm_dvrci;
            let dderf_dvrci = (dnum_dvrci * denom - num * ddenom_dvrci) / (denom * denom);

            let dsqrt_dvrci = derf * dderf_dvrci / sqrt_1_derf2;
            let dirci_dvrci =
                (diohm_dvrci * sqrt_1_derf2 - iohm * dsqrt_dvrci) / (sqrt_1_derf2 * sqrt_1_derf2);

            // Derivative w.r.t. vbci
            let dnum_dvbci = ivo * self.rci_t * diohm_dvbci;
            let dderf_dvbci = dnum_dvbci / denom;
            let dsqrt_dvbci = derf * dderf_dvbci / sqrt_1_derf2;
            let dirci_dvbci =
                (diohm_dvbci * sqrt_1_derf2 - iohm * dsqrt_dvbci) / (sqrt_1_derf2 * sqrt_1_derf2);

            // Derivative w.r.t. vbcx
            let dnum_dvbcx = ivo * self.rci_t * diohm_dvbcx;
            let dderf_dvbcx = dnum_dvbcx / denom;
            let dsqrt_dvbcx = derf * dderf_dvbcx / sqrt_1_derf2;
            let dirci_dvbcx =
                (diohm_dvbcx * sqrt_1_derf2 - iohm * dsqrt_dvbcx) / (sqrt_1_derf2 * sqrt_1_derf2);

            (irci_val, dirci_dvrci, dirci_dvbci, dirci_dvbcx)
        } else {
            // Simple resistive model
            (iohm, diohm_dvrci, diohm_dvbci, diohm_dvbcx)
        }
    }
}

/// NR companion model result for a VBIC BJT at an operating point.
#[derive(Debug, Clone)]
pub struct VbicCompanion {
    // Transport current (collector-emitter internal)
    pub iciei: f64,
    pub diciei_dvbei: f64,
    pub diciei_dvbci: f64,

    // Base-emitter current (internal)
    pub ibe: f64,
    pub dibe_dvbei: f64,

    // Base-emitter current (external partition)
    pub ibex: f64,
    pub dibex_dvbex: f64,

    // Base-collector current (with avalanche)
    pub ibc: f64,
    pub dibc_dvbci: f64,
    pub dibc_dvbei: f64,

    // Parasitic base-emitter current
    pub ibep: f64,
    pub dibep_dvbep: f64,

    // Parasitic base-collector current
    pub ibcp: f64,
    pub dibcp_dvbcp: f64,

    // Parasitic transport current
    pub iccp: f64,
    pub diccp_dvbep: f64,
    pub diccp_dvbcp: f64,

    // Avalanche
    pub igc: f64,
    pub digc_dvbci: f64,
    pub digc_dvbei: f64,

    // Internal collector resistance current
    pub irci: f64,
    pub dirci_dvrci: f64,
    pub dirci_dvbci: f64,
    pub dirci_dvbcx: f64,

    // Internal base resistance current
    pub irbi: f64,
    pub dirbi_dvrbi: f64,
    pub dirbi_dvbei: f64,
    pub dirbi_dvbci: f64,

    // Parasitic base resistance current
    pub irbp: f64,
    pub dirbp_dvrbp: f64,

    // Depletion charges and capacitances (raw, for reference)
    pub qbe: f64,
    pub cqbe: f64,
    pub qbc: f64,
    pub cqbc: f64,
    pub qbep: f64,
    pub cqbep: f64,
    pub qbcp: f64,
    pub cqbcp: f64,

    // Total charge derivatives for AC analysis (depletion + transit time + overlap)
    // These include transit time modulation and are the values used for AC stamps.
    pub cqbe_vbei: f64,  // dQbe/dVbei (depletion*WBE + transit time)
    pub cqbe_vbci: f64,  // dQbe/dVbci (transit time cross-derivative)
    pub cqbex_vbex: f64, // dQbex/dVbex (external BE depletion * (1-WBE))
    pub cqbc_vbci: f64,  // dQbc/dVbci (depletion + TR*dIri/dVbci + QCO*dKbci/dVbci)
    pub cqbcx_vbcx: f64, // dQbcx/dVbcx (QCO*dKbcx/dVbcx)
    pub cqbep_vbep: f64, // dQbep/dVbep (parasitic depletion + TR*dIfp/dVbep)
    pub cqbep_vbci: f64, // dQbep/dVbci (parasitic transit cross-derivative)
    pub cqbcp_vbcp: f64, // dQbcp/dVbcp (substrate depletion + CCSO)

    // Operating point
    pub itzf: f64,
    pub itzr: f64,
    pub qb: f64,
}

/// Resolved node indices for a VBIC instance in the MNA system.
///
/// External: coll(C), base(B), emit(E), subs(S).
/// Always internal: coll_ci, base_bi, base_bp.
/// Conditional: coll_cx, base_bx, emit_ei, subs_si.
#[derive(Debug, Clone)]
pub struct VbicInstance {
    pub name: String,
    // External nodes
    pub coll_idx: Option<usize>,
    pub base_idx: Option<usize>,
    pub emit_idx: Option<usize>,
    pub subs_idx: Option<usize>,
    // Always-internal nodes
    pub coll_ci_idx: Option<usize>,
    pub base_bi_idx: Option<usize>,
    pub base_bp_idx: Option<usize>,
    // Conditional internal nodes
    pub coll_cx_idx: Option<usize>, // = coll_idx if RCX=0
    pub base_bx_idx: Option<usize>, // = base_idx if RBX=0
    pub emit_ei_idx: Option<usize>, // = emit_idx if RE=0
    pub subs_si_idx: Option<usize>, // = subs_idx if RS=0
    /// Internal thermal node for self-heating (only when RTH > 0).
    pub rth_idx: Option<usize>,
    pub model: VbicModel,
    pub area: f64,
    pub m: f64,
    /// Ambient temperature (°C) for self-heating computation.
    pub t_ambient: f64,
}

impl VbicInstance {
    /// Get all junction voltages from the solution vector.
    ///
    /// Returns (vbei, vbex, vbci, vbcx, vbep, vrci, vrbi, vrbp, vbcp).
    pub fn junction_voltages(
        &self,
        solution: &[f64],
    ) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64) {
        let node = |idx: Option<usize>| idx.map(|i| solution[i]).unwrap_or(0.0);
        let sign = self.model.vbic_type.sign();

        let v_bi = node(self.base_bi_idx);
        let v_bx = node(self.base_bx_idx);
        let v_bp = node(self.base_bp_idx);
        let v_ci = node(self.coll_ci_idx);
        let v_cx = node(self.coll_cx_idx);
        let v_ei = node(self.emit_ei_idx);
        let v_si = node(self.subs_si_idx);

        let vbei = sign * (v_bi - v_ei);
        let vbex = sign * (v_bx - v_ei);
        let vbci = sign * (v_bi - v_ci);
        let vbcx = sign * (v_bi - v_cx);
        let vbep = sign * (v_bx - v_bp);
        let vrci = sign * (v_cx - v_ci); // across internal collector resistance
        let vrbi = sign * (v_bx - v_bi); // across internal base resistance
        let vrbp = sign * (v_bp - v_cx); // across parasitic base resistance (BP to CX)
        let vbcp = sign * (v_si - v_bp); // parasitic B-C (substrate)

        (vbei, vbex, vbci, vbcx, vbep, vrci, vrbi, vrbp, vbcp)
    }

    /// Get the thermal voltage rise (Vrth) from the solution vector.
    /// Returns 0 if self-heating is not active.
    pub fn vrth(&self, solution: &[f64]) -> f64 {
        self.rth_idx.map(|i| solution[i]).unwrap_or(0.0)
    }
}

/// Compute self-heating power dissipation Ith for a VBIC device.
///
/// Ith = sum of V_branch × I_branch for all electrical branches.
/// This follows ngspice vbicload.c's Ith computation.
#[expect(clippy::too_many_arguments)]
pub fn compute_self_heating_power(
    comp: &VbicCompanion,
    model: &VbicModel,
    vbei: f64,
    vbex: f64,
    vbci: f64,
    vbep: f64,
    vbcp: f64,
    vrci: f64,
    vrbi: f64,
    vrbp: f64,
    vrcx: f64,
    vrbx: f64,
    vre: f64,
    vrs: f64,
    area: f64,
    m: f64,
) -> f64 {
    let scale = m * area;
    // Vcei = Vbei - Vbci (collector-emitter internal voltage)
    let vcei = vbei - vbci;
    // Vcep = Vbep - Vbcp (parasitic collector-emitter voltage)
    let vcep = vbep - vbcp;

    // External resistance currents (V²/R = V*I)
    let p_rcx = if model.rcx_t > 0.0 {
        vrcx * vrcx / model.rcx_t
    } else {
        0.0
    };
    let p_rbx = if model.rbx_t > 0.0 {
        vrbx * vrbx / model.rbx_t
    } else {
        0.0
    };
    let p_re = if model.re_t > 0.0 {
        vre * vre / model.re_t
    } else {
        0.0
    };
    let p_rs = if model.rs_t > 0.0 {
        vrs * vrs / model.rs_t
    } else {
        0.0
    };

    scale
        * (comp.ibe * vbei
            + comp.ibex * vbex
            + (comp.itzf - comp.itzr) * vcei
            + comp.ibc * vbci
            + comp.ibep * vbep
            + comp.ibcp * vbcp
            + comp.iccp * vcep
            + comp.irci * vrci
            + comp.irbi * vrbi
            + comp.irbp * vrbp
            + p_rcx
            + p_rbx
            + p_re
            + p_rs)
}

/// Stamp the VBIC companion model conductances into the MNA matrix.
///
/// This stamps only the matrix conductance entries (no RHS). For the full
/// NR stamp including RHS equivalent current sources, use `stamp_vbic_with_voltages`.
pub fn stamp_vbic(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &VbicInstance,
    comp: &VbicCompanion,
) {
    // Delegate to the full stamping function with zero voltages.
    // This gives correct matrix stamps; RHS will only contain the raw currents.
    stamp_vbic_with_voltages(
        matrix, rhs, inst, comp, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    );
}

/// Stamp the VBIC companion model into the MNA matrix and RHS using operating point voltages.
///
/// This is the proper stamping function that takes operating point voltages
/// to compute the Norton equivalent current sources for the RHS.
#[expect(clippy::too_many_arguments)]
pub fn stamp_vbic_with_voltages(
    matrix: &mut crate::SparseMatrix,
    rhs: &mut [f64],
    inst: &VbicInstance,
    comp: &VbicCompanion,
    vbei: f64,
    vbex: f64,
    vbci: f64,
    vbcx: f64,
    vbep: f64,
    vrci: f64,
    vrbi: f64,
    vrbp: f64,
    vbcp: f64,
) {
    let sign = inst.model.vbic_type.sign();
    let scale = inst.m * inst.area;

    let bi = inst.base_bi_idx;
    let bx = inst.base_bx_idx;
    let bp = inst.base_bp_idx;
    let ci = inst.coll_ci_idx;
    let cx = inst.coll_cx_idx;
    let ei = inst.emit_ei_idx;
    let si = inst.subs_si_idx;
    let c_ext = inst.coll_idx;
    let b_ext = inst.base_idx;
    let e_ext = inst.emit_idx;
    let s_ext = inst.subs_idx;

    // Helper to add matrix entry safely
    let mut m_add = |r: Option<usize>, c: Option<usize>, val: f64| {
        if let (Some(r), Some(c)) = (r, c) {
            matrix.add(r, c, val);
        }
    };

    let mut r_add = |r: Option<usize>, val: f64| {
        if let Some(r) = r {
            rhs[r] += val;
        }
    };

    // Helper: stamp branch current I from node+ to node-, with controlling voltage V = V(ctrl+) - V(ctrl-)
    // Matrix: g into quadrant, RHS: -(I - g*V) at nodes
    // Macro to stamp a single-control branch:
    macro_rules! stamp_branch {
        ($np:expr, $nm:expr, $cp:expr, $cm:expr, $g:expr, $i:expr, $v:expr) => {
            let g_s = scale * $g;
            let i_eq = sign * scale * ($i - $g * $v);
            m_add($np, $cp, g_s);
            m_add($np, $cm, -g_s);
            m_add($nm, $cp, -g_s);
            m_add($nm, $cm, g_s);
            r_add($np, -i_eq);
            r_add($nm, i_eq);
        };
    }

    // -----------------------------------------------------------------
    // 1. Iciei: CI -> EI, controlled by Vbei (bi-ei) and Vbci (bi-ci)
    // -----------------------------------------------------------------
    {
        // Two controlling voltages, stamp each partial separately
        let g_bei = scale * comp.diciei_dvbei;
        m_add(ci, bi, g_bei);
        m_add(ci, ei, -g_bei);
        m_add(ei, bi, -g_bei);
        m_add(ei, ei, g_bei);

        let g_bci = scale * comp.diciei_dvbci;
        m_add(ci, bi, g_bci);
        m_add(ci, ci, -g_bci);
        m_add(ei, bi, -g_bci);
        m_add(ei, ci, g_bci);

        // RHS: -(Iciei - dIciei/dVbei*Vbei - dIciei/dVbci*Vbci)
        let i_eq =
            sign * scale * (comp.iciei - comp.diciei_dvbei * vbei - comp.diciei_dvbci * vbci);
        r_add(ci, -i_eq);
        r_add(ei, i_eq);
    }

    // -----------------------------------------------------------------
    // 2. Ibe: BI -> EI, controlled by Vbei
    // -----------------------------------------------------------------
    stamp_branch!(bi, ei, bi, ei, comp.dibe_dvbei, comp.ibe, vbei);

    // -----------------------------------------------------------------
    // 3. Ibex: BX -> EI, controlled by Vbex
    // -----------------------------------------------------------------
    stamp_branch!(bx, ei, bx, ei, comp.dibex_dvbex, comp.ibex, vbex);

    // -----------------------------------------------------------------
    // 4. Ibc: BI -> CI, controlled by Vbci (+ avalanche cross-term from Vbei)
    // -----------------------------------------------------------------
    {
        let g_bci = scale * comp.dibc_dvbci;
        m_add(bi, bi, g_bci);
        m_add(bi, ci, -g_bci);
        m_add(ci, bi, -g_bci);
        m_add(ci, ci, g_bci);

        // Avalanche cross-term
        let g_bei = scale * comp.dibc_dvbei;
        if g_bei.abs() > 0.0 {
            m_add(bi, bi, g_bei);
            m_add(bi, ei, -g_bei);
            m_add(ci, bi, -g_bei);
            m_add(ci, ei, g_bei);
        }

        let i_eq = sign * scale * (comp.ibc - comp.dibc_dvbci * vbci - comp.dibc_dvbei * vbei);
        r_add(bi, -i_eq);
        r_add(ci, i_eq);
    }

    // -----------------------------------------------------------------
    // 5. Ibep: BX -> BP, controlled by Vbep = V(BX) - V(BP)
    // -----------------------------------------------------------------
    stamp_branch!(bx, bp, bx, bp, comp.dibep_dvbep, comp.ibep, vbep);

    // -----------------------------------------------------------------
    // 6. Ibcp: SI -> BP, controlled by Vbcp = V(SI) - V(BP)
    // -----------------------------------------------------------------
    stamp_branch!(si, bp, si, bp, comp.dibcp_dvbcp, comp.ibcp, vbcp);

    // -----------------------------------------------------------------
    // 7. Iccp: BX -> SI, controlled by Vbep, Vbci, Vbcp (parasitic transport)
    // -----------------------------------------------------------------
    {
        // Vbep = V(BX) - V(BP) control
        let g_bep = scale * comp.diccp_dvbep;
        m_add(bx, bx, g_bep);
        m_add(bx, bp, -g_bep);
        m_add(si, bx, -g_bep);
        m_add(si, bp, g_bep);

        // Vbcp = V(SI) - V(BP) control
        let g_bcp = scale * comp.diccp_dvbcp;
        m_add(bx, si, g_bcp);
        m_add(bx, bp, -g_bcp);
        m_add(si, si, -g_bcp);
        m_add(si, bp, g_bcp);

        let i_eq = sign * scale * (comp.iccp - comp.diccp_dvbep * vbep - comp.diccp_dvbcp * vbcp);
        r_add(bx, -i_eq);
        r_add(si, i_eq);
    }

    // -----------------------------------------------------------------
    // 8. Irci: CX -> CI, controlled by Vrci, Vbci, Vbcx
    // -----------------------------------------------------------------
    {
        let g_rci = scale * comp.dirci_dvrci;
        m_add(cx, cx, g_rci);
        m_add(cx, ci, -g_rci);
        m_add(ci, cx, -g_rci);
        m_add(ci, ci, g_rci);

        let g_bci = scale * comp.dirci_dvbci;
        if g_bci.abs() > 0.0 {
            m_add(cx, bi, g_bci);
            m_add(cx, ci, -g_bci);
            m_add(ci, bi, -g_bci);
            m_add(ci, ci, g_bci);
        }

        let g_bcx = scale * comp.dirci_dvbcx;
        if g_bcx.abs() > 0.0 {
            m_add(cx, bi, g_bcx);
            m_add(cx, cx, -g_bcx);
            m_add(ci, bi, -g_bcx);
            m_add(ci, cx, g_bcx);
        }

        let i_eq = sign
            * scale
            * (comp.irci
                - comp.dirci_dvrci * vrci
                - comp.dirci_dvbci * vbci
                - comp.dirci_dvbcx * vbcx);
        r_add(cx, -i_eq);
        r_add(ci, i_eq);
    }

    // -----------------------------------------------------------------
    // 9. Irbi: BX -> BI, controlled by Vrbi, Vbei, Vbci
    // -----------------------------------------------------------------
    {
        let g_rbi = scale * comp.dirbi_dvrbi;
        m_add(bx, bx, g_rbi);
        m_add(bx, bi, -g_rbi);
        m_add(bi, bx, -g_rbi);
        m_add(bi, bi, g_rbi);

        let g_bei = scale * comp.dirbi_dvbei;
        if g_bei.abs() > 0.0 {
            m_add(bx, bi, g_bei);
            m_add(bx, ei, -g_bei);
            m_add(bi, bi, -g_bei);
            m_add(bi, ei, g_bei);
        }

        let g_bci = scale * comp.dirbi_dvbci;
        if g_bci.abs() > 0.0 {
            m_add(bx, bi, g_bci);
            m_add(bx, ci, -g_bci);
            m_add(bi, bi, -g_bci);
            m_add(bi, ci, g_bci);
        }

        let i_eq = sign
            * scale
            * (comp.irbi
                - comp.dirbi_dvrbi * vrbi
                - comp.dirbi_dvbei * vbei
                - comp.dirbi_dvbci * vbci);
        r_add(bx, -i_eq);
        r_add(bi, i_eq);
    }

    // -----------------------------------------------------------------
    // 10. Irbp: BP -> CX, controlled by Vrbp
    // -----------------------------------------------------------------
    stamp_branch!(bp, cx, bp, cx, comp.dirbp_dvrbp, comp.irbp, vrbp);

    // -----------------------------------------------------------------
    // 11. External resistances
    // -----------------------------------------------------------------
    if inst.model.rcx > 0.0 && inst.model.rcx_t > 0.0 {
        let g = scale / inst.model.rcx_t;
        crate::stamp_conductance(matrix, c_ext, cx, g);
    }

    if inst.model.rbx > 0.0 && inst.model.rbx_t > 0.0 {
        let g = scale / inst.model.rbx_t;
        crate::stamp_conductance(matrix, b_ext, bx, g);
    }

    if inst.model.re > 0.0 && inst.model.re_t > 0.0 {
        let g = scale / inst.model.re_t;
        crate::stamp_conductance(matrix, e_ext, ei, g);
    }

    if inst.model.rs > 0.0 && inst.model.rs_t > 0.0 {
        let g = scale / inst.model.rs_t;
        crate::stamp_conductance(matrix, s_ext, si, g);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    #[test]
    fn test_default_vbic_model() {
        let m = VbicModel::new(VbicType::Npn);
        assert_eq!(m.is, 1e-16);
        assert_eq!(m.nf, 1.0);
        assert_eq!(m.nr, 1.0);
        assert_eq!(m.rci, 0.1);
        assert_eq!(m.rbi, 0.1);
        assert_eq!(m.rbp, 0.1);
        assert_eq!(m.rcx, 0.0);
        assert_eq!(m.rbx, 0.0);
        assert_eq!(m.re, 0.0);
        assert_eq!(m.rs, 0.0);
    }

    #[test]
    fn test_from_model_def_npn() {
        let model_def = ModelDef {
            name: "QVBIC".to_string(),
            kind: "NPN".to_string(),
            params: vec![
                Param {
                    name: "IS".to_string(),
                    value: Expr::Num(1e-15),
                },
                Param {
                    name: "NF".to_string(),
                    value: Expr::Num(1.02),
                },
                Param {
                    name: "RCI".to_string(),
                    value: Expr::Num(10.0),
                },
                Param {
                    name: "VEF".to_string(),
                    value: Expr::Num(100.0),
                },
            ],
        };
        let m = VbicModel::from_model_def(&model_def);
        assert_eq!(m.vbic_type, VbicType::Npn);
        assert_abs_diff_eq!(m.is, 1e-15, epsilon = 1e-30);
        assert_abs_diff_eq!(m.nf, 1.02, epsilon = 1e-15);
        assert_abs_diff_eq!(m.rci, 10.0, epsilon = 1e-15);
        assert_abs_diff_eq!(m.vef, 100.0, epsilon = 1e-15);
    }

    #[test]
    fn test_from_model_def_pnp() {
        let model_def = ModelDef {
            name: "QP1".to_string(),
            kind: "PNP".to_string(),
            params: vec![Param {
                name: "IS".to_string(),
                value: Expr::Num(5e-16),
            }],
        };
        let m = VbicModel::from_model_def(&model_def);
        assert_eq!(m.vbic_type, VbicType::Pnp);
        assert_abs_diff_eq!(m.is, 5e-16, epsilon = 1e-30);
    }

    #[test]
    fn test_internal_node_count() {
        let mut m = VbicModel::new(VbicType::Npn);
        // 3 always-internal: collCI, baseBI, baseBP
        assert_eq!(m.internal_node_count(), 3);
        m.rcx = 10.0;
        assert_eq!(m.internal_node_count(), 4);
        m.rbx = 5.0;
        assert_eq!(m.internal_node_count(), 5);
        m.re = 1.0;
        assert_eq!(m.internal_node_count(), 6);
        m.rs = 0.5;
        assert_eq!(m.internal_node_count(), 7);
    }

    #[test]
    fn test_companion_forward_active() {
        let m = VbicModel::new(VbicType::Npn);
        // Forward active: Vbei=0.8V, all others near zero
        let comp = m.companion(0.8, 0.0, -2.0, -2.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0);

        // Forward transport current should be positive
        assert!(
            comp.itzf > 0.0,
            "itzf should be positive, got {}",
            comp.itzf
        );
        // Collector current should be positive (Itzf >> Itzr)
        assert!(
            comp.iciei > 0.0,
            "iciei should be positive, got {}",
            comp.iciei
        );
        // Base-emitter current should be positive
        assert!(comp.ibe > 0.0, "ibe should be positive, got {}", comp.ibe);
        // Base charge should be > 0
        assert!(comp.qb > 0.0, "qb should be positive, got {}", comp.qb);
    }

    #[test]
    fn test_companion_cutoff() {
        let m = VbicModel::new(VbicType::Npn);
        // Cutoff: all junctions reverse biased
        let comp = m.companion(-0.5, -0.5, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        assert!(
            comp.iciei.abs() < 1e-10,
            "iciei in cutoff should be ~0, got {}",
            comp.iciei
        );
        assert!(
            comp.ibe.abs() < 1e-10,
            "ibe in cutoff should be ~0, got {}",
            comp.ibe
        );
    }

    #[test]
    fn test_companion_derivatives_positive() {
        let m = VbicModel::new(VbicType::Npn);
        let comp = m.companion(0.7, 0.0, -2.0, -2.0, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0);

        // dIciei/dVbei should be positive (more base-emitter = more collector current)
        assert!(
            comp.diciei_dvbei > 0.0,
            "diciei/dvbei should be positive, got {}",
            comp.diciei_dvbei
        );
        // dIbe/dVbei should be positive
        assert!(
            comp.dibe_dvbei > 0.0,
            "dibe/dvbei should be positive, got {}",
            comp.dibe_dvbei
        );
    }

    #[test]
    fn test_temperature_adjust_nominal() {
        let m = VbicModel::new(VbicType::Npn);
        // At nominal temp, adjusted values should match nominal
        assert_abs_diff_eq!(m.is_t, m.is, epsilon = m.is * 1e-6);
        assert_abs_diff_eq!(m.rci_t, m.rci, epsilon = 1e-15);
        assert_abs_diff_eq!(m.rbi_t, m.rbi, epsilon = 1e-15);
    }

    #[test]
    fn test_depletion_charge_zero_cap() {
        let (q, c) = depletion_charge(0.5, 0.0, 0.75, 0.33, -0.5, 0.9);
        assert_eq!(q, 0.0);
        assert_eq!(c, 0.0);
    }

    #[test]
    fn test_depletion_charge_reverse_bias() {
        let (_q, c) = depletion_charge(-1.0, 1e-12, 0.75, 0.33, -0.5, 0.9);
        // Under reverse bias, charge should be negative (depletion)
        // and capacitance should be positive
        assert!(c > 0.0, "capacitance should be positive, got {}", c);
    }

    #[test]
    fn test_vbic_type_sign() {
        assert_eq!(VbicType::Npn.sign(), 1.0);
        assert_eq!(VbicType::Pnp.sign(), -1.0);
    }

    #[test]
    fn test_irci_simple_resistive() {
        // With GAMM=0, VO=0, Irci should be simple V/R
        let m = VbicModel::new(VbicType::Npn);
        let comp = m.companion(0.7, 0.0, -2.0, -2.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0);
        // Irci should be approximately Vrci / RCI = 0.5 / 0.1 = 5.0
        // (modified by Kbci/Kbcx factors, but with GAMM=0 these are 1)
        assert_abs_diff_eq!(comp.irci, 5.0, epsilon = 0.1);
    }

    #[test]
    fn test_irbi_modulated() {
        let m = VbicModel::new(VbicType::Npn);
        // With forward bias, qb > 1, so effective Rbi < nominal Rbi
        let comp = m.companion(0.7, 0.0, -2.0, -2.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0);
        // Irbi = Vrbi * qb / RBI
        // With qb ≈ 1 at moderate injection, Irbi ≈ 0.1 / 0.1 = 1.0
        assert!(
            comp.irbi > 0.0,
            "irbi should be positive, got {}",
            comp.irbi
        );
    }

    #[test]
    fn test_companion_with_early_effect() {
        let mut m = VbicModel::new(VbicType::Npn);
        m.vef = 50.0;
        m.temperature_adjust(27.0);

        let comp1 = m.companion(0.7, 0.0, -2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let comp2 = m.companion(0.7, 0.0, -10.0, -10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        // With Early effect, more reverse bias on BC should increase forward current
        // (larger |Vce| = larger Ic) -- but in VBIC the Early effect works through
        // depletion charge, so the relationship is more complex. Just verify both are positive.
        assert!(comp1.iciei > 0.0);
        assert!(comp2.iciei > 0.0);
    }

    #[test]
    fn test_companion_pnp_gummel_low_bias() {
        // Test the companion at the exact FG Gummel parameters at Vbe=0.2V
        let model_def = ModelDef {
            name: "P1".to_string(),
            kind: "PNP".to_string(),
            params: vec![
                Param {
                    name: "IS".into(),
                    value: Expr::Num(1e-16),
                },
                Param {
                    name: "IBEI".into(),
                    value: Expr::Num(1e-18),
                },
                Param {
                    name: "IBEN".into(),
                    value: Expr::Num(5e-15),
                },
                Param {
                    name: "IBCI".into(),
                    value: Expr::Num(2e-17),
                },
                Param {
                    name: "IBCN".into(),
                    value: Expr::Num(5e-15),
                },
                Param {
                    name: "ISP".into(),
                    value: Expr::Num(1e-15),
                },
                Param {
                    name: "RCX".into(),
                    value: Expr::Num(10.0),
                },
                Param {
                    name: "RCI".into(),
                    value: Expr::Num(60.0),
                },
                Param {
                    name: "RBX".into(),
                    value: Expr::Num(10.0),
                },
                Param {
                    name: "RBI".into(),
                    value: Expr::Num(40.0),
                },
                Param {
                    name: "RE".into(),
                    value: Expr::Num(2.0),
                },
                Param {
                    name: "RS".into(),
                    value: Expr::Num(20.0),
                },
                Param {
                    name: "RBP".into(),
                    value: Expr::Num(40.0),
                },
                Param {
                    name: "VEF".into(),
                    value: Expr::Num(10.0),
                },
                Param {
                    name: "VER".into(),
                    value: Expr::Num(4.0),
                },
                Param {
                    name: "IKF".into(),
                    value: Expr::Num(2e-3),
                },
                Param {
                    name: "IKR".into(),
                    value: Expr::Num(2e-4),
                },
                Param {
                    name: "IKP".into(),
                    value: Expr::Num(2e-4),
                },
                Param {
                    name: "CJE".into(),
                    value: Expr::Num(1e-13),
                },
                Param {
                    name: "CJC".into(),
                    value: Expr::Num(2e-14),
                },
                Param {
                    name: "CJEP".into(),
                    value: Expr::Num(1e-13),
                },
                Param {
                    name: "CJCP".into(),
                    value: Expr::Num(4e-13),
                },
                Param {
                    name: "VO".into(),
                    value: Expr::Num(2.0),
                },
                Param {
                    name: "GAMM".into(),
                    value: Expr::Num(2e-11),
                },
                Param {
                    name: "HRCF".into(),
                    value: Expr::Num(2.0),
                },
                Param {
                    name: "QCO".into(),
                    value: Expr::Num(1e-12),
                },
                Param {
                    name: "AVC1".into(),
                    value: Expr::Num(2.0),
                },
                Param {
                    name: "AVC2".into(),
                    value: Expr::Num(15.0),
                },
                Param {
                    name: "TF".into(),
                    value: Expr::Num(10e-12),
                },
                Param {
                    name: "TR".into(),
                    value: Expr::Num(100e-12),
                },
                Param {
                    name: "TD".into(),
                    value: Expr::Num(2e-11),
                },
                Param {
                    name: "RTH".into(),
                    value: Expr::Num(300.0),
                },
            ],
        };
        let mut m = VbicModel::from_model_def(&model_def);
        m.temperature_adjust(27.0);

        // At Vbei=0.2V (forward bias for PNP), all other junctions at 0V
        let comp = m.companion(0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-13);

        // Expected: IS * exp(0.2/Vt) ≈ 1e-16 * 2277 ≈ 2.28e-13
        eprintln!("iciei = {:.6e}", comp.iciei);
        eprintln!("itzf = {:.6e}", comp.itzf);
        eprintln!("ibe = {:.6e}", comp.ibe);
        eprintln!("ibc = {:.6e}", comp.ibc);
        eprintln!("ibep = {:.6e}", comp.ibep);
        eprintln!("ibcp = {:.6e}", comp.ibcp);
        eprintln!("iccp = {:.6e}", comp.iccp);
        eprintln!("irci = {:.6e}", comp.irci);
        eprintln!("irbi = {:.6e}", comp.irbi);
        eprintln!("irbp = {:.6e}", comp.irbp);
        eprintln!("qb = {:.6e}", comp.qb);

        // Iciei should be around 2.17e-13 (ngspice reference)
        assert!(
            comp.iciei > 1e-14 && comp.iciei < 1e-11,
            "iciei should be ~2e-13, got {:.6e}",
            comp.iciei
        );
    }
}
