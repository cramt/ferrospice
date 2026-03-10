//! Shared Newton-Raphson device stamping logic.
//!
//! This module extracts the nonlinear device companion stamping loop that is
//! shared between DC operating point (`simulate.rs`) and transient analysis
//! (`transient.rs`). Each device type follows the same pattern: extract terminal
//! voltages → apply voltage limiting → compute companion model → stamp MNA.

use std::cell::RefCell;

use crate::LinearSystem;
use crate::bjt::stamp_bjt;
use crate::bsim3::{bsim3_companion, bsim3_limit, stamp_bsim3};
use crate::bsim3soi_dd::{bsim3soi_dd_companion, bsim3soi_dd_limit, stamp_bsim3soi_dd};
use crate::bsim3soi_fd::{bsim3soi_fd_companion, bsim3soi_fd_limit, stamp_bsim3soi_fd};
use crate::bsim3soi_pd::{bsim3soi_pd_companion, bsim3soi_pd_limit, stamp_bsim3soi_pd};
use crate::bsim4::{bsim4_companion, bsim4_limit, stamp_bsim4};
use crate::diode::{VT_NOM, pnjlim, vcrit};
use crate::jfet::{jfet_limit, stamp_jfet};
use crate::mna::{MnaSystem, stamp_conductance};
use crate::mosfet::{mos_limit, stamp_mosfet};
use crate::vbic::stamp_vbic_with_voltages;

/// Stamp a current source into the RHS vector.
/// Current flows from ni (pos) to nj (neg) externally:
/// subtract from ni, add to nj.
pub(crate) fn stamp_current_source(
    rhs: &mut [f64],
    ni: Option<usize>,
    nj: Option<usize>,
    i_val: f64,
) {
    if let Some(i) = ni {
        rhs[i] -= i_val;
    }
    if let Some(j) = nj {
        rhs[j] += i_val;
    }
}

/// FET voltage limiting (ngspice `DEVfetlim`).
///
/// Prevents large voltage jumps during NR iteration that could cause
/// divergence in FET device models. Used by BSIM3, BSIM4, and SOI variants.
pub(crate) fn fetlim(vnew: f64, vold: f64, vto: f64) -> f64 {
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
    } else if delv.abs() >= vtsthi {
        if vold < vnew {
            vold + vtsthi
        } else {
            vold - vtsthi
        }
    } else {
        vnew
    }
}

/// Simplified PN junction voltage limiting for BSIM models.
///
/// Unlike `diode::pnjlim` which takes explicit `vt` and `vcrit` parameters,
/// this version uses hardcoded defaults (IS=1e-14, Vt at 300.15K) matching
/// ngspice's BSIM3/BSIM4 implementations.
pub(crate) fn bsim_pnjlim(vnew: f64, vold: f64) -> f64 {
    const KBOQ: f64 = 8.617_333_262e-5;
    const TEMP_DEFAULT: f64 = 300.15;
    let vt = KBOQ * TEMP_DEFAULT;
    let vcrit_val = vt * (vt / (std::f64::consts::SQRT_2 * 1e-14)).ln();
    if vnew > vcrit_val && (vnew - vold).abs() > 2.0 * vt {
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

/// Tracks previous junction/terminal voltages for NR voltage limiting.
///
/// Each device type needs its previous voltages preserved between NR iterations
/// to compute limiting deltas. This struct bundles all device voltage state and
/// the precomputed critical voltages.
pub(crate) struct DeviceVoltageState {
    prev_jct: RefCell<Vec<f64>>,
    vcrits: Vec<f64>,
    prev_bjt: RefCell<Vec<(f64, f64)>>,
    prev_mos: RefCell<Vec<(f64, f64)>>,
    prev_mos6: RefCell<Vec<(f64, f64)>>,
    prev_jfet: RefCell<Vec<(f64, f64)>>,
    jfet_vcrits: Vec<f64>,
    prev_bsim3: RefCell<Vec<(f64, f64, f64)>>,
    prev_bsim3soi_pd: RefCell<Vec<(f64, f64, f64, f64)>>,
    prev_bsim3soi_fd: RefCell<Vec<(f64, f64, f64, f64)>>,
    prev_bsim3soi_dd: RefCell<Vec<(f64, f64, f64, f64)>>,
    prev_bsim4: RefCell<Vec<(f64, f64, f64)>>,
    prev_vbic: RefCell<Vec<(f64, f64)>>,
}

impl DeviceVoltageState {
    /// Create state with all previous voltages set to zero (for DC operating point).
    pub fn new_zero(mna: &MnaSystem) -> Self {
        let vcrits: Vec<f64> = mna
            .diodes
            .iter()
            .map(|d| vcrit(d.model.n * VT_NOM, d.model.is))
            .collect();
        let jfet_vcrits: Vec<f64> = mna
            .jfets
            .iter()
            .map(|j| vcrit(j.model.n * VT_NOM, j.model.is))
            .collect();

        Self {
            prev_jct: RefCell::new(vec![0.0; mna.diodes.len()]),
            vcrits,
            prev_bjt: RefCell::new(vec![(0.0, 0.0); mna.bjts.len()]),
            prev_mos: RefCell::new(vec![(0.0, 0.0); mna.mosfets.len()]),
            prev_mos6: RefCell::new(vec![(0.0, 0.0); mna.mos6s.len()]),
            prev_jfet: RefCell::new(vec![(0.0, 0.0); mna.jfets.len()]),
            jfet_vcrits,
            prev_bsim3: RefCell::new(vec![(0.0, 0.0, 0.0); mna.bsim3s.len()]),
            prev_bsim3soi_pd: RefCell::new(vec![(0.0, 0.0, 0.0, 0.0); mna.bsim3soi_pds.len()]),
            prev_bsim3soi_fd: RefCell::new(vec![(0.0, 0.0, 0.0, 0.0); mna.bsim3soi_fds.len()]),
            prev_bsim3soi_dd: RefCell::new(vec![(0.0, 0.0, 0.0, 0.0); mna.bsim3soi_dds.len()]),
            prev_bsim4: RefCell::new(vec![(0.0, 0.0, 0.0); mna.bsim4s.len()]),
            prev_vbic: RefCell::new(vec![(0.0, 0.0); mna.vbics.len()]),
        }
    }

    /// Create state initialized from a previous solution vector (for transient).
    pub fn from_solution(mna: &MnaSystem, prev_solution: &[f64]) -> Self {
        let vcrits: Vec<f64> = mna
            .diodes
            .iter()
            .map(|d| vcrit(d.model.n * VT_NOM, d.model.is))
            .collect();
        let jfet_vcrits: Vec<f64> = mna
            .jfets
            .iter()
            .map(|j| vcrit(j.model.n * VT_NOM, j.model.is))
            .collect();

        let prev_jct = mna
            .diodes
            .iter()
            .map(|d| {
                let (ja, jc) = if d.internal_idx.is_some() {
                    (d.internal_idx, d.cathode_idx)
                } else {
                    (d.anode_idx, d.cathode_idx)
                };
                let va = ja.map(|i| prev_solution[i]).unwrap_or(0.0);
                let vc = jc.map(|i| prev_solution[i]).unwrap_or(0.0);
                va - vc
            })
            .collect();

        Self {
            prev_jct: RefCell::new(prev_jct),
            vcrits,
            prev_bjt: RefCell::new(
                mna.bjts
                    .iter()
                    .map(|b| b.junction_voltages(prev_solution))
                    .collect(),
            ),
            prev_mos: RefCell::new(
                mna.mosfets
                    .iter()
                    .map(|m| {
                        let (vgs, vds, _vbs) = m.terminal_voltages(prev_solution);
                        (vgs, vds)
                    })
                    .collect(),
            ),
            prev_mos6: RefCell::new(
                mna.mos6s
                    .iter()
                    .map(|m| {
                        let (vgs, vds, _vbs) = m.terminal_voltages(prev_solution);
                        (vgs, vds)
                    })
                    .collect(),
            ),
            prev_jfet: RefCell::new(
                mna.jfets
                    .iter()
                    .map(|j| j.junction_voltages(prev_solution))
                    .collect(),
            ),
            jfet_vcrits,
            prev_bsim3: RefCell::new(
                mna.bsim3s
                    .iter()
                    .map(|b| b.terminal_voltages(prev_solution))
                    .collect(),
            ),
            prev_bsim3soi_pd: RefCell::new(
                mna.bsim3soi_pds
                    .iter()
                    .map(|b| b.terminal_voltages(prev_solution))
                    .collect(),
            ),
            prev_bsim3soi_fd: RefCell::new(
                mna.bsim3soi_fds
                    .iter()
                    .map(|b| b.terminal_voltages(prev_solution))
                    .collect(),
            ),
            prev_bsim3soi_dd: RefCell::new(
                mna.bsim3soi_dds
                    .iter()
                    .map(|b| b.terminal_voltages(prev_solution))
                    .collect(),
            ),
            prev_bsim4: RefCell::new(
                mna.bsim4s
                    .iter()
                    .map(|b| b.terminal_voltages(prev_solution))
                    .collect(),
            ),
            prev_vbic: RefCell::new(vec![(0.0, 0.0); mna.vbics.len()]),
        }
    }

    /// Stamp all nonlinear device companion models into the MNA system.
    ///
    /// This is the shared NR load logic used by both DC operating point and
    /// transient analysis. For each device type: extract terminal voltages from
    /// the current solution, apply voltage limiting, compute the companion model,
    /// and stamp it into the linear system.
    pub fn stamp_devices(&self, solution: &[f64], system: &mut LinearSystem, mna: &MnaSystem) {
        // Diodes
        {
            let mut prev = self.prev_jct.borrow_mut();
            for (di, diode) in mna.diodes.iter().enumerate() {
                let (jct_anode, jct_cathode) = if diode.internal_idx.is_some() {
                    (diode.internal_idx, diode.cathode_idx)
                } else {
                    (diode.anode_idx, diode.cathode_idx)
                };

                let v_anode = jct_anode.map(|i| solution[i]).unwrap_or(0.0);
                let v_cathode = jct_cathode.map(|i| solution[i]).unwrap_or(0.0);
                let mut v_jct = v_anode - v_cathode;

                v_jct = pnjlim(v_jct, prev[di], diode.model.n * VT_NOM, self.vcrits[di]);
                prev[di] = v_jct;

                let (g_d, i_eq) = diode.model.companion(v_jct);
                stamp_conductance(&mut system.matrix, jct_anode, jct_cathode, g_d);
                stamp_current_source(&mut system.rhs, jct_anode, jct_cathode, i_eq);

                if let Some(int_idx) = diode.internal_idx {
                    let g_rs = 1.0 / diode.model.rs;
                    stamp_conductance(&mut system.matrix, diode.anode_idx, Some(int_idx), g_rs);
                }
            }
        }

        // BJTs
        {
            let mut prev = self.prev_bjt.borrow_mut();
            for (bi, bjt) in mna.bjts.iter().enumerate() {
                let (raw_vbe, raw_vbc) = bjt.junction_voltages(solution);

                let vbe = bjt.model.limit_vbe(raw_vbe, prev[bi].0);
                let vbc = bjt.model.limit_vbc(raw_vbc, prev[bi].1);
                prev[bi] = (vbe, vbc);

                let comp = bjt.model.companion(vbe, vbc);
                stamp_bjt(&mut system.matrix, &mut system.rhs, bjt, &comp);
            }
        }

        // MOSFETs (level 1)
        {
            let mut prev = self.prev_mos.borrow_mut();
            for (mi, mos) in mna.mosfets.iter().enumerate() {
                let (raw_vgs, raw_vds, vbs) = mos.terminal_voltages(solution);

                let (vgs, vds) = mos_limit(raw_vgs, raw_vds, prev[mi].0, prev[mi].1, mos.model.vto);
                prev[mi] = (vgs, vds);

                let mut eff_model = mos.model.clone();
                eff_model.kp = mos.beta();

                let comp = eff_model.companion(vgs, vds, vbs);
                stamp_mosfet(&mut system.matrix, &mut system.rhs, mos, &comp);
            }
        }

        // MOS6 (level 6)
        {
            let mut prev = self.prev_mos6.borrow_mut();
            for (mi, mos) in mna.mos6s.iter().enumerate() {
                let (raw_vgs, raw_vds, vbs) = mos.terminal_voltages(solution);

                let (vgs, vds) = mos_limit(raw_vgs, raw_vds, prev[mi].0, prev[mi].1, mos.model.vto);
                prev[mi] = (vgs, vds);

                let betac = mos.betac();
                let comp = mos.model.companion(vgs, vds, vbs, betac);
                crate::mos6::stamp_mos6(&mut system.matrix, &mut system.rhs, mos, &comp);
            }
        }

        // JFETs
        {
            let mut prev = self.prev_jfet.borrow_mut();
            for (ji, jfet) in mna.jfets.iter().enumerate() {
                let (raw_vgs, raw_vgd) = jfet.junction_voltages(solution);

                let vt = jfet.model.n * VT_NOM;
                let (vgs, vgd) = jfet_limit(
                    raw_vgs,
                    raw_vgd,
                    prev[ji].0,
                    prev[ji].1,
                    vt,
                    self.jfet_vcrits[ji],
                );
                prev[ji] = (vgs, vgd);

                let comp = jfet.model.companion(vgs, vgd);
                stamp_jfet(&mut system.matrix, &mut system.rhs, jfet, &comp);
            }
        }

        // BSIM3
        {
            let mut prev = self.prev_bsim3.borrow_mut();
            for (bi, bsim) in mna.bsim3s.iter().enumerate() {
                let (raw_vgs, raw_vds, raw_vbs) = bsim.terminal_voltages(solution);

                let (vgs, vds, vbs) = bsim3_limit(
                    raw_vgs,
                    raw_vds,
                    raw_vbs,
                    prev[bi].0,
                    prev[bi].1,
                    prev[bi].2,
                    bsim.vth0_inst,
                );
                prev[bi] = (vgs, vds, vbs);

                let comp = bsim3_companion(vgs, vds, vbs, &bsim.size_params, &bsim.model);
                stamp_bsim3(&mut system.matrix, &mut system.rhs, bsim, &comp);
            }
        }

        // BSIM3SOI-PD
        {
            let mut prev = self.prev_bsim3soi_pd.borrow_mut();
            for (bi, bsim) in mna.bsim3soi_pds.iter().enumerate() {
                let (raw_vgs, raw_vds, raw_vbs, raw_ves) = bsim.terminal_voltages(solution);

                let (vgs, vds, vbs, ves) = bsim3soi_pd_limit(
                    raw_vgs,
                    raw_vds,
                    raw_vbs,
                    raw_ves,
                    prev[bi].0,
                    prev[bi].1,
                    prev[bi].2,
                    prev[bi].3,
                    bsim.vth0_inst,
                );
                prev[bi] = (vgs, vds, vbs, ves);

                let comp =
                    bsim3soi_pd_companion(vgs, vds, vbs, ves, &bsim.size_params, &bsim.model);
                stamp_bsim3soi_pd(&mut system.matrix, &mut system.rhs, bsim, &comp);
            }
        }

        // BSIM3SOI-FD
        {
            let mut prev = self.prev_bsim3soi_fd.borrow_mut();
            for (bi, bsim) in mna.bsim3soi_fds.iter().enumerate() {
                let (raw_vgs, raw_vds, raw_vbs, raw_ves) = bsim.terminal_voltages(solution);

                let (vgs, vds, vbs, ves) = bsim3soi_fd_limit(
                    raw_vgs,
                    raw_vds,
                    raw_vbs,
                    raw_ves,
                    prev[bi].0,
                    prev[bi].1,
                    prev[bi].2,
                    prev[bi].3,
                    bsim.vth0_inst,
                );
                prev[bi] = (vgs, vds, vbs, ves);

                let comp =
                    bsim3soi_fd_companion(vgs, vds, vbs, ves, &bsim.size_params, &bsim.model);
                stamp_bsim3soi_fd(&mut system.matrix, &mut system.rhs, bsim, &comp);
            }
        }

        // BSIM3SOI-DD
        {
            let mut prev = self.prev_bsim3soi_dd.borrow_mut();
            for (bi, bsim) in mna.bsim3soi_dds.iter().enumerate() {
                let (raw_vgs, raw_vds, raw_vbs, raw_ves) = bsim.terminal_voltages(solution);

                let (vgs, vds, vbs, ves) = bsim3soi_dd_limit(
                    raw_vgs,
                    raw_vds,
                    raw_vbs,
                    raw_ves,
                    prev[bi].0,
                    prev[bi].1,
                    prev[bi].2,
                    prev[bi].3,
                    bsim.vth0_inst,
                );
                prev[bi] = (vgs, vds, vbs, ves);

                let comp =
                    bsim3soi_dd_companion(vgs, vds, vbs, ves, &bsim.size_params, &bsim.model);
                stamp_bsim3soi_dd(&mut system.matrix, &mut system.rhs, bsim, &comp);
            }
        }

        // BSIM4
        {
            let mut prev = self.prev_bsim4.borrow_mut();
            for (bi, bsim) in mna.bsim4s.iter().enumerate() {
                let (raw_vgs, raw_vds, raw_vbs) = bsim.terminal_voltages(solution);

                let (vgs, vds, vbs) = bsim4_limit(
                    raw_vgs,
                    raw_vds,
                    raw_vbs,
                    prev[bi].0,
                    prev[bi].1,
                    prev[bi].2,
                    bsim.size_params.vth0,
                );
                prev[bi] = (vgs, vds, vbs);

                let comp = bsim4_companion(vgs, vds, vbs, &bsim.size_params, &bsim.model);
                stamp_bsim4(&mut system.matrix, &mut system.rhs, bsim, &comp);
            }
        }

        // VBIC
        {
            let mut prev = self.prev_vbic.borrow_mut();
            for (vi, vbic) in mna.vbics.iter().enumerate() {
                let (raw_vbei, vbex, raw_vbci, vbcx, vbep, vrci, vrbi, vrbp, vbcp) =
                    vbic.junction_voltages(solution);

                let vbei = vbic.model.limit_vbei(raw_vbei, prev[vi].0);
                let vbci = vbic.model.limit_vbci(raw_vbci, prev[vi].1);
                prev[vi] = (vbei, vbci);

                let comp = vbic
                    .model
                    .companion(vbei, vbex, vbci, vbcx, vbep, vrci, vrbi, vrbp, vbcp);
                stamp_vbic_with_voltages(
                    &mut system.matrix,
                    &mut system.rhs,
                    vbic,
                    &comp,
                    vbei,
                    vbex,
                    vbci,
                    vbcx,
                    vbep,
                    vrci,
                    vrbi,
                    vrbp,
                    vbcp,
                );
            }
        }
    }
}
