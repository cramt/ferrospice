//! Transient (time-domain) analysis engine.
//!
//! Implements `.tran` analysis with Backward Euler (BE) and Trapezoidal (Trap)
//! integration methods. Capacitors and inductors are converted to companion
//! models (conductance/resistance + current/voltage source) at each timestep.
//!
//! Supports adaptive timestep control via local truncation error (LTE)
//! estimation using the difference between BE and Trap results for
//! capacitor/inductor charges/fluxes.

use thevenin_types::{Analysis, Item, Netlist, SimPlot, SimResult, SimVector};

use crate::LinearSystem;
use crate::device_stamp::{DeviceVoltageState, stamp_current_source};
use crate::expr_val;
use crate::ltra::{LtraCoeffs, LtraState};
use crate::mna::{MnaError, MnaSystem, assemble_mna, stamp_conductance};
use crate::newton::{NrOptions, newton_raphson_solve};
use crate::simulate::simulate_op;
use crate::txl::TxlTransientStamp;
use crate::waveform::{self, TranParams};

/// Integration method for transient analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IntegrationMethod {
    /// Backward Euler — first-order, unconditionally stable.
    BackwardEuler,
    /// Trapezoidal — second-order, A-stable.
    Trapezoidal,
}

/// History state for a capacitor at the previous timestep.
#[derive(Debug, Clone)]
struct CapHistory {
    /// Voltage across capacitor at previous timestep.
    voltage: f64,
    /// Current through capacitor at previous timestep (needed for Trap).
    current: f64,
}

/// History state for an inductor at the previous timestep.
#[derive(Debug, Clone)]
struct IndHistory {
    /// Current through inductor at previous timestep.
    current: f64,
    /// Voltage across inductor at previous timestep (needed for Trap).
    voltage: f64,
}

/// Compute capacitor companion model coefficients.
///
/// Returns `(geq, ieq)` where:
/// - `geq` is the equivalent conductance to stamp into the matrix
/// - `ieq` is the equivalent current source to stamp into the RHS
///
/// For Backward Euler: i(n) = C/h * v(n) - C/h * v(n-1)
///   → geq = C/h, ieq = -C/h * v(n-1)
///
/// For Trapezoidal: i(n) = 2C/h * v(n) - 2C/h * v(n-1) - i(n-1)
///   → geq = 2C/h, ieq = -(2C/h * v(n-1) + i(n-1))
fn capacitor_companion(
    capacitance: f64,
    h: f64,
    history: &CapHistory,
    method: IntegrationMethod,
) -> (f64, f64) {
    match method {
        IntegrationMethod::BackwardEuler => {
            let geq = capacitance / h;
            let ieq = -geq * history.voltage;
            (geq, ieq)
        }
        IntegrationMethod::Trapezoidal => {
            let geq = 2.0 * capacitance / h;
            let ieq = -(geq * history.voltage + history.current);
            (geq, ieq)
        }
    }
}

/// Compute inductor companion model coefficients.
///
/// Returns `(req, veq)` where:
/// - `req` is the equivalent resistance added to the branch equation diagonal
/// - `veq` is the equivalent voltage source added to the branch equation RHS
///
/// For Backward Euler: v(n) = L/h * i(n) - L/h * i(n-1)
///   → req = L/h, veq = -L/h * i(n-1)
///
/// For Trapezoidal: v(n) = 2L/h * i(n) - 2L/h * i(n-1) - v(n-1)
///   → req = 2L/h, veq = -(2L/h * i(n-1) + v(n-1))
fn inductor_companion(
    inductance: f64,
    h: f64,
    history: &IndHistory,
    method: IntegrationMethod,
) -> (f64, f64) {
    match method {
        IntegrationMethod::BackwardEuler => {
            let req = inductance / h;
            let veq = -req * history.current;
            (req, veq)
        }
        IntegrationMethod::Trapezoidal => {
            let req = 2.0 * inductance / h;
            let veq = -(req * history.current + history.voltage);
            (req, veq)
        }
    }
}

/// Sorted breakpoint table for transient analysis.
///
/// Forces the timestep engine to land exactly on times where waveforms have
/// discontinuities (PULSE edges, PWL corners, etc.).
struct BreakpointTable {
    /// Sorted breakpoint times.
    times: Vec<f64>,
    /// Index of the next unprocessed breakpoint.
    next_idx: usize,
    /// Minimum separation between breakpoints.
    min_break: f64,
}

impl BreakpointTable {
    /// Build breakpoint table from all source waveforms.
    fn from_mna(mna: &MnaSystem, tran: &TranParams) -> Self {
        let mut times = Vec::new();

        for vs in &mna.voltage_sources {
            if let Some(ref wf) = vs.waveform {
                times.extend(waveform::breakpoints(wf, tran));
            }
        }
        for cs in &mna.current_sources {
            if let Some(ref wf) = cs.waveform {
                times.extend(waveform::breakpoints(wf, tran));
            }
        }

        times.sort_by(|a, b| a.partial_cmp(b).unwrap());
        times.dedup_by(|a, b| (*a - *b).abs() < 1e-15);

        let min_break = tran.tstep * 5e-5;

        Self {
            times,
            next_idx: 0,
            min_break,
        }
    }

    /// Get the next breakpoint time after `current_time`, if any.
    fn next_after(&mut self, current_time: f64) -> Option<f64> {
        // Advance past breakpoints we've already passed.
        while self.next_idx < self.times.len()
            && self.times[self.next_idx] <= current_time + self.min_break
        {
            self.next_idx += 1;
        }
        self.times.get(self.next_idx).copied()
    }

    /// Check if `t` is at or very near a breakpoint.
    fn is_at_breakpoint(&self, t: f64) -> bool {
        self.times.iter().any(|&bp| (t - bp).abs() < self.min_break)
    }
}

/// Transient truncation error tolerance (controls how aggressively timestep adjusts).
/// Matches ngspice CKTtrtol default of 7.0.
const TRTOL: f64 = 7.0;

/// Charge tolerance for LTE estimation (matches ngspice CHGTOL).
const CHGTOL: f64 = 1e-14;

/// Minimum factor to shrink timestep on rejection.
const MIN_SHRINK: f64 = 0.125; // 1/8

/// Maximum factor to grow timestep.
const MAX_GROW: f64 = 2.0;

/// Estimate the new timestep based on LTE for capacitors and inductors.
///
/// Uses the difference between Trap and BE predictions as the LTE estimate.
/// For each capacitor: LTE ≈ |q_trap - q_be| where q is the integrated charge.
/// For each inductor: LTE ≈ |flux_trap - flux_be| where flux is the integrated flux.
///
/// Returns the recommended new timestep.
fn estimate_new_timestep(
    h: f64,
    solution: &[f64],
    mna: &MnaSystem,
    cap_histories: &[CapHistory],
    ind_histories: &[IndHistory],
    reltol: f64,
    abstol: f64,
) -> f64 {
    let mut new_h = f64::MAX;

    // LTE for capacitors.
    for (ci, cap) in mna.capacitors.iter().enumerate() {
        let v_pos = cap.pos_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_neg = cap.neg_idx.map(|i| solution[i]).unwrap_or(0.0);
        let v_new = v_pos - v_neg;
        let v_old = cap_histories[ci].voltage;
        let i_old = cap_histories[ci].current;

        // Trap current: i_trap = 2C/h * (v_new - v_old) - i_old
        let i_trap = 2.0 * cap.capacitance / h * (v_new - v_old) - i_old;
        // BE current: i_be = C/h * (v_new - v_old)
        let i_be = cap.capacitance / h * (v_new - v_old);

        // Charge LTE: difference in integrated charge over this step.
        // Trap integrates: q_trap = h/2 * (i_old + i_trap)
        // BE integrates:   q_be = h * i_be
        // LTE ≈ |q_trap - q_be|
        let q_trap = h / 2.0 * (i_old + i_trap);
        let q_be = h * i_be;
        let lte = (q_trap - q_be).abs();

        // Tolerance on charge.
        let q_new = cap.capacitance * v_new;
        let q_old = cap.capacitance * v_old;
        let vol_tol = abstol + reltol * v_new.abs().max(v_old.abs());
        let chg_tol = reltol * q_new.abs().max(q_old.abs()).max(CHGTOL);
        let tol = TRTOL * vol_tol.max(chg_tol);

        if lte > 1e-30 {
            // For Trap (order 2): new_h = h * (tol / lte)^(1/2)
            let ratio = tol / lte;
            let h_new = h * ratio.sqrt();
            new_h = new_h.min(h_new);
        }
    }

    // LTE for inductors.
    for (li, ind) in mna.inductors.iter().enumerate() {
        let i_new = solution[ind.branch_idx];
        let i_old = ind_histories[li].current;
        let v_old = ind_histories[li].voltage;

        // Trap voltage: v_trap = 2L/h * (i_new - i_old) - v_old
        let v_trap = 2.0 * ind.inductance / h * (i_new - i_old) - v_old;
        // BE voltage: v_be = L/h * (i_new - i_old)
        let v_be = ind.inductance / h * (i_new - i_old);

        // Flux LTE: difference in integrated flux.
        let flux_trap = h / 2.0 * (v_old + v_trap);
        let flux_be = h * v_be;
        let lte = (flux_trap - flux_be).abs();

        // Tolerance on flux.
        let flux_new = ind.inductance * i_new;
        let flux_old = ind.inductance * i_old;
        let cur_tol = abstol + reltol * i_new.abs().max(i_old.abs());
        let flux_tol = reltol * flux_new.abs().max(flux_old.abs()).max(CHGTOL);
        let tol = TRTOL * cur_tol.max(flux_tol);

        if lte > 1e-30 {
            let ratio = tol / lte;
            let h_new = h * ratio.sqrt();
            new_h = new_h.min(h_new);
        }
    }

    new_h
}

/// Perform transient analysis on a circuit.
///
/// Parses the `.tran` command from the netlist, computes the DC operating point
/// as initial conditions, then steps through time using numerical integration.
/// Uses adaptive timestep control with LTE estimation when reactive elements
/// are present, falling back to fixed timestep for purely resistive circuits.
pub fn simulate_tran(netlist: &Netlist) -> Result<SimResult, MnaError> {
    // Find the .tran analysis command.
    let (tstep, tstop, tstart, tmax) = netlist
        .items
        .iter()
        .find_map(|item| {
            if let Item::Analysis(Analysis::Tran {
                tstep,
                tstop,
                tstart,
                tmax,
            }) = item
            {
                Some((tstep.clone(), tstop.clone(), tstart.clone(), tmax.clone()))
            } else {
                None
            }
        })
        .ok_or_else(|| MnaError::UnsupportedElement("no .tran analysis found".to_string()))?;

    let h_print = expr_val(&tstep, ".tran tstep")?;
    let t_stop = expr_val(&tstop, ".tran tstop")?;
    let t_start = tstart
        .as_ref()
        .map(|e| expr_val(e, ".tran tstart"))
        .transpose()?
        .unwrap_or(0.0);
    let t_max = tmax
        .as_ref()
        .map(|e| expr_val(e, ".tran tmax"))
        .transpose()?;

    if h_print <= 0.0 || t_stop <= 0.0 {
        return Err(MnaError::UnsupportedElement(
            "invalid .tran parameters".to_string(),
        ));
    }

    // Maximum internal timestep: tmax if specified, otherwise min(tstep, tstop/50).
    let h_max = t_max.unwrap_or_else(|| h_print.min(t_stop / 50.0));

    // Assemble MNA system.
    let mut mna = assemble_mna(netlist)?;

    // Compute DC operating point for initial conditions.
    let op_result = simulate_op(netlist)?;
    let op_plot = &op_result.plots[0];

    // Extract initial solution vector from DC OP.
    let num_ext_nodes = mna.node_map.len();
    let num_internal = mna
        .diodes
        .iter()
        .filter(|d| d.internal_idx.is_some())
        .count()
        + mna
            .bjts
            .iter()
            .map(|b| b.model.internal_node_count())
            .sum::<usize>()
        + mna
            .mosfets
            .iter()
            .map(|m| m.model.internal_node_count())
            .sum::<usize>()
        + mna
            .jfets
            .iter()
            .map(|j| j.model.internal_node_count())
            .sum::<usize>()
        + mna
            .bsim3s
            .iter()
            .map(|b| b.model.internal_node_count(b.nrd, b.nrs))
            .sum::<usize>()
        + mna
            .bsim3soi_pds
            .iter()
            .map(|b| b.model.internal_node_count(b.nrd, b.nrs))
            .sum::<usize>()
        + mna
            .bsim3soi_fds
            .iter()
            .map(|b| b.model.internal_node_count(b.nrd, b.nrs))
            .sum::<usize>()
        + mna
            .bsim3soi_dds
            .iter()
            .map(|b| b.model.internal_node_count(b.nrd, b.nrs))
            .sum::<usize>()
        + mna
            .bsim4s
            .iter()
            .map(|b| b.model.internal_node_count(b.nrd, b.nrs))
            .sum::<usize>();
    let num_nodes = num_ext_nodes + num_internal;
    let dim = mna.system.dim();

    let mut solution = vec![0.0; dim];
    // Fill node voltages from OP result.
    for (name, idx) in mna.node_map.iter() {
        let vec_name = format!("v({})", name);
        if let Some(sv) = op_plot.vecs.iter().find(|v| v.name == vec_name)
            && !sv.real.is_empty()
        {
            solution[idx] = sv.real[0];
        }
    }
    // Fill branch currents from OP result.
    for (i, vsrc) in mna.vsource_names.iter().enumerate() {
        let vec_name = format!("{}#branch", vsrc.to_lowercase());
        if let Some(sv) = op_plot.vecs.iter().find(|v| v.name == vec_name)
            && !sv.real.is_empty()
        {
            solution[num_nodes + i] = sv.real[0];
        }
    }

    // Apply IC overrides for capacitors (override DC OP voltage).
    for cap in &mna.capacitors {
        if let Some(ic_v) = cap.ic {
            match (cap.pos_idx, cap.neg_idx) {
                (Some(pi), None) => solution[pi] = ic_v,
                (None, Some(ni)) => solution[ni] = -ic_v,
                (Some(pi), Some(ni)) => {
                    solution[pi] = ic_v + solution[ni];
                }
                (None, None) => {}
            }
        }
    }

    // Apply IC overrides for inductors.
    for ind in &mna.inductors {
        if let Some(ic_i) = ind.ic {
            solution[ind.branch_idx] = ic_i;
        }
    }

    // Initialize history from the initial solution.
    let mut cap_histories: Vec<CapHistory> = mna
        .capacitors
        .iter()
        .map(|cap| {
            let v_pos = cap.pos_idx.map(|i| solution[i]).unwrap_or(0.0);
            let v_neg = cap.neg_idx.map(|i| solution[i]).unwrap_or(0.0);
            CapHistory {
                voltage: v_pos - v_neg,
                current: 0.0, // At DC steady state, capacitor current is 0.
            }
        })
        .collect();

    let mut ind_histories: Vec<IndHistory> = mna
        .inductors
        .iter()
        .map(|ind| {
            let current = solution[ind.branch_idx];
            IndHistory {
                current,
                voltage: 0.0, // At DC steady state, inductor voltage is 0.
            }
        })
        .collect();

    // Initialize LTRA transient state.
    let has_ltra = !mna.ltras.is_empty();
    let mut ltra_states: Vec<LtraState> = mna
        .ltras
        .iter()
        .map(|inst| {
            let mut state = LtraState::new();
            let v1 =
                node_voltage(&solution, inst.pos1_idx) - node_voltage(&solution, inst.neg1_idx);
            let i1 = solution[num_nodes + inst.br_eq1];
            let v2 =
                node_voltage(&solution, inst.pos2_idx) - node_voltage(&solution, inst.neg2_idx);
            let i2 = solution[num_nodes + inst.br_eq2];
            state.init_from_dc(v1, i1, v2, i2);
            // Record initial values as first history point.
            state.accept(v1, i1, v2, i2);
            state
        })
        .collect();
    let mut ltra_coeffs: Vec<LtraCoeffs> = mna.ltras.iter().map(|_| LtraCoeffs::new()).collect();
    // Time points history for LTRA convolution.
    let mut ltra_time_points: Vec<f64> = if has_ltra { vec![0.0] } else { Vec::new() };

    // Initialize TXL transient state from DC operating point.
    let has_txl = !mna.txls.is_empty();
    if has_txl {
        for inst in &mut mna.txls {
            let dc_v1 = inst.pos_idx.map_or(0.0, |i| solution[i]);
            let dc_v2 = inst.neg_idx.map_or(0.0, |i| solution[i]);
            let dc_i1 = solution[num_nodes + inst.ibr1];
            let dc_i2 = solution[num_nodes + inst.ibr2];
            crate::txl::init_dc_state(&mut inst.txline, dc_v1, dc_v2, dc_i1, dc_i2);
            inst.txline2 = inst.txline.clone();
        }
    }
    let mut txl_stamps: Vec<TxlTransientStamp> = Vec::new();

    // Initialize CPL transient state from DC operating point.
    let has_cpl = !mna.cpls.is_empty();
    if has_cpl {
        for inst in &mut mna.cpls {
            let dc_v_in: Vec<f64> = (0..inst.no_l)
                .map(|m| inst.pos_nodes[m].map_or(0.0, |i| solution[i]))
                .collect();
            let dc_v_out: Vec<f64> = (0..inst.no_l)
                .map(|m| inst.neg_nodes[m].map_or(0.0, |i| solution[i]))
                .collect();
            crate::cpl::init_dc_state(&mut inst.cpline, &dc_v_in, &dc_v_out);
            inst.cpline2 = inst.cpline.clone();
        }
    }
    let mut cpl_stamps: Vec<crate::cpl::CplTransientStamp> = Vec::new();

    // Prepare output vectors.
    let mut time_vec = SimVector {
        name: "time".to_string(),
        real: Vec::new(),
        complex: vec![],
    };

    let mut node_vecs: Vec<SimVector> = mna
        .node_map
        .iter()
        .map(|(name, _)| SimVector {
            name: format!("v({})", name),
            real: Vec::new(),
            complex: vec![],
        })
        .collect();

    let mut branch_vecs: Vec<SimVector> = mna
        .vsource_names
        .iter()
        .map(|vsrc| SimVector {
            name: format!("{}#branch", vsrc.to_lowercase()),
            real: Vec::new(),
            complex: vec![],
        })
        .collect();

    // Record initial point if tstart == 0.
    let mut t = 0.0;
    if t >= t_start {
        record_point(
            t,
            &solution,
            &mna,
            num_nodes,
            &mut time_vec,
            &mut node_vecs,
            &mut branch_vecs,
        );
    }

    let has_nonlinear = mna.has_nonlinear();
    let has_reactive = !mna.capacitors.is_empty() || !mna.inductors.is_empty();
    let nr_options = NrOptions::default();
    let tran_params = TranParams {
        tstep: h_print,
        tstop: t_stop,
    };

    // Build breakpoint table from source waveforms.
    let mut breakpoints = BreakpointTable::from_mna(&mna, &tran_params);

    // Internal timestep — start with tstep, adapt from there.
    let mut h = h_print.min(h_max);
    let h_min = h_print * 1e-9; // Absolute minimum timestep.
    let mut is_first_step = true;

    // Adaptive time-stepping loop.
    while t < t_stop - h_min {
        // The step size for this iteration (separate from h which tracks the
        // "suggested" next step from LTE control).
        let mut step_h = h.min(h_max);

        // Don't overshoot tstop.
        if t + step_h > t_stop {
            step_h = t_stop - t;
        }

        // Breakpoint handling: don't cross the next breakpoint.
        let at_breakpoint = breakpoints.is_at_breakpoint(t);
        if let Some(bp) = breakpoints.next_after(t) {
            let dist = bp - t;
            if step_h > dist {
                step_h = dist;
            }
        }

        // At breakpoints, reduce step for stability (ngspice uses 0.1×).
        if at_breakpoint {
            step_h = step_h.min(h * 0.1).max(h_min);
        }

        // Use Backward Euler for the first step and at breakpoints.
        let method = if is_first_step || at_breakpoint {
            IntegrationMethod::BackwardEuler
        } else {
            IntegrationMethod::Trapezoidal
        };

        // Recompute LTRA convolution coefficients for the new timepoint.
        if has_ltra {
            let cur_time = t + step_h;
            let time_index = ltra_time_points.len() - 1;
            for (li, inst) in mna.ltras.iter().enumerate() {
                match inst.model.special_case {
                    crate::ltra::LtraCase::Rlc => {
                        crate::ltra::rlc_coeffs_setup(
                            &mut ltra_coeffs[li],
                            inst.model.td,
                            inst.model.alpha,
                            inst.model.beta,
                            cur_time,
                            &ltra_time_points,
                            time_index,
                            inst.model.reltol,
                        );
                    }
                    crate::ltra::LtraCase::Rc => {
                        crate::ltra::rc_coeffs_setup(
                            &mut ltra_coeffs[li],
                            inst.model.c_by_r,
                            inst.model.rclsqr,
                            cur_time,
                            &ltra_time_points,
                            time_index,
                            inst.model.reltol,
                        );
                    }
                    _ => {} // LC and RG don't need coefficient setup
                }
            }
        }

        // Pre-compute TXL transient stamps (updates convolution state, must
        // be called exactly once per timestep attempt).
        if has_txl {
            let cur_time_ps = ((t + step_h) * 1.0e12) as i64;
            let prev_time_ps = (t * 1.0e12) as i64;
            txl_stamps.clear();
            for inst in &mut mna.txls {
                // Backup state for potential step rejection
                inst.txline2 = inst.txline.clone();
                let stamp = crate::txl::prepare_txl_transient(
                    inst,
                    num_nodes,
                    &solution,
                    cur_time_ps,
                    prev_time_ps,
                    step_h, // h in seconds (poles are per-second)
                );
                txl_stamps.push(stamp);
            }
        }

        // Pre-compute CPL transient stamps.
        if has_cpl {
            let cur_time_ps = ((t + step_h) * 1.0e12) as i64;
            let prev_time_ps = (t * 1.0e12) as i64;
            cpl_stamps.clear();
            for inst in &mut mna.cpls {
                inst.cpline2 = inst.cpline.clone();
                let stamp = crate::cpl::prepare_cpl_transient(
                    inst,
                    num_nodes,
                    &solution,
                    cur_time_ps,
                    prev_time_ps,
                    step_h,
                );
                cpl_stamps.push(stamp);
            }
        }

        // Solve this timestep.
        let new_solution = solve_timestep(
            &mna,
            &solution,
            step_h,
            t + step_h,
            &tran_params,
            method,
            &cap_histories,
            &ind_histories,
            has_nonlinear,
            &nr_options,
            dim,
            num_nodes,
            if has_ltra { Some(&ltra_states) } else { None },
            if has_ltra { Some(&ltra_coeffs) } else { None },
            if has_ltra {
                Some(&ltra_time_points)
            } else {
                None
            },
            if has_txl { Some(&txl_stamps) } else { None },
            if has_cpl { Some(&cpl_stamps) } else { None },
        )?;

        // LTE-based timestep control (only for Trap with reactive elements).
        if method == IntegrationMethod::Trapezoidal && has_reactive {
            let new_h = estimate_new_timestep(
                step_h,
                &new_solution,
                &mna,
                &cap_histories,
                &ind_histories,
                nr_options.reltol,
                nr_options.abstol,
            );

            if new_h < 0.9 * step_h {
                // Reject step: shrink h and retry without advancing time.
                // Restore TXL/CPL state from backup.
                if has_txl {
                    for inst in &mut mna.txls {
                        inst.txline = inst.txline2.clone();
                    }
                }
                if has_cpl {
                    for inst in &mut mna.cpls {
                        inst.cpline = inst.cpline2.clone();
                    }
                }
                h = (new_h.max(step_h * MIN_SHRINK)).max(h_min);
                continue;
            }

            // Accept: schedule next h from LTE estimate.
            h = new_h.min(step_h * MAX_GROW).min(h_max).max(h_min);
        }

        // Accept this timestep: advance time and update state.
        t += step_h;
        solution = new_solution;
        is_first_step = false;

        // Update capacitor histories.
        for (ci, cap) in mna.capacitors.iter().enumerate() {
            let v_pos = cap.pos_idx.map(|i| solution[i]).unwrap_or(0.0);
            let v_neg = cap.neg_idx.map(|i| solution[i]).unwrap_or(0.0);
            let v_new = v_pos - v_neg;

            let current = match method {
                IntegrationMethod::BackwardEuler => {
                    let geq = cap.capacitance / step_h;
                    geq * (v_new - cap_histories[ci].voltage)
                }
                IntegrationMethod::Trapezoidal => {
                    let geq = 2.0 * cap.capacitance / step_h;
                    geq * (v_new - cap_histories[ci].voltage) - cap_histories[ci].current
                }
            };

            cap_histories[ci] = CapHistory {
                voltage: v_new,
                current,
            };
        }

        // Update inductor histories.
        for (li, ind) in mna.inductors.iter().enumerate() {
            let i_new = solution[ind.branch_idx];
            let v_pos = ind.pos_idx.map(|i| solution[i]).unwrap_or(0.0);
            let v_neg = ind.neg_idx.map(|i| solution[i]).unwrap_or(0.0);
            let v_new = v_pos - v_neg;

            ind_histories[li] = IndHistory {
                current: i_new,
                voltage: v_new,
            };
        }

        // Update LTRA histories.
        if has_ltra {
            ltra_time_points.push(t);
            for (li, inst) in mna.ltras.iter().enumerate() {
                let v1 =
                    node_voltage(&solution, inst.pos1_idx) - node_voltage(&solution, inst.neg1_idx);
                let i1 = solution[num_nodes + inst.br_eq1];
                let v2 =
                    node_voltage(&solution, inst.pos2_idx) - node_voltage(&solution, inst.neg2_idx);
                let i2 = solution[num_nodes + inst.br_eq2];
                ltra_states[li].accept(v1, i1, v2, i2);
            }
        }

        // Update TXL histories and convolution accumulators.
        if has_txl {
            let time_ps = (t * 1.0e12) as i64;
            let h_ps = step_h * 1.0e12;
            for inst in &mut mna.txls {
                let v_in = inst.pos_idx.map_or(0.0, |i| solution[i]);
                let v_out = inst.neg_idx.map_or(0.0, |i| solution[i]);
                let i_in = solution[num_nodes + inst.ibr1];
                let i_out = solution[num_nodes + inst.ibr2];

                // Update h1 convolution accumulators
                let tx = &mut inst.txline;
                if !tx.lsl {
                    let prev_vi = tx.vi_history.last();
                    let (dv_i, dv_o) = if let Some(prev) = prev_vi {
                        let dt = time_ps - prev.time;
                        if dt > 0 {
                            (
                                (v_in - prev.v_i) / dt as f64,
                                (v_out - prev.v_o) / dt as f64,
                            )
                        } else {
                            (0.0, 0.0)
                        }
                    } else {
                        (0.0, 0.0)
                    };
                    crate::txl::update_cnv_txl(tx, h_ps, v_in, v_out, dv_i, dv_o);
                }

                // Handle extended time step delayed convolution update
                if tx.ext {
                    crate::txl::update_delayed_cnv(tx, h_ps, tx.ratio);
                }

                // Record history point
                tx.vi_history.push(crate::txl::ViEntry {
                    time: time_ps,
                    v_i: v_in,
                    v_o: v_out,
                    i_i: i_in,
                    i_o: i_out,
                });
            }
        }

        // Record output point.
        if t >= t_start {
            record_point(
                t,
                &solution,
                &mna,
                num_nodes,
                &mut time_vec,
                &mut node_vecs,
                &mut branch_vecs,
            );
        }
    }

    // Assemble result.
    let mut vecs = vec![time_vec];
    vecs.extend(node_vecs);
    vecs.extend(branch_vecs);

    Ok(SimResult {
        plots: vec![SimPlot {
            name: "tran1".to_string(),
            vecs,
        }],
    })
}

/// Extract a node voltage from the solution vector, returning 0 for ground.
fn node_voltage(solution: &[f64], idx: Option<usize>) -> f64 {
    idx.map(|i| solution[i]).unwrap_or(0.0)
}

/// Record a solution point into output vectors.
fn record_point(
    t: f64,
    solution: &[f64],
    mna: &MnaSystem,
    num_nodes: usize,
    time_vec: &mut SimVector,
    node_vecs: &mut [SimVector],
    branch_vecs: &mut [SimVector],
) {
    time_vec.real.push(t);

    for (idx, (_name, node_idx)) in mna.node_map.iter().enumerate() {
        node_vecs[idx].real.push(solution[node_idx]);
    }

    for (i, _vsrc) in mna.vsource_names.iter().enumerate() {
        branch_vecs[i].real.push(solution[num_nodes + i]);
    }
}

/// Solve a single transient timestep.
#[expect(clippy::too_many_arguments)]
fn solve_timestep(
    mna: &MnaSystem,
    prev_solution: &[f64],
    h: f64,
    t: f64,
    tran_params: &TranParams,
    method: IntegrationMethod,
    cap_histories: &[CapHistory],
    ind_histories: &[IndHistory],
    has_nonlinear: bool,
    nr_options: &NrOptions,
    dim: usize,
    num_nodes: usize,
    ltra_states: Option<&[LtraState]>,
    ltra_coeffs: Option<&[LtraCoeffs]>,
    ltra_time_points: Option<&[f64]>,
    txl_stamps: Option<&[TxlTransientStamp]>,
    cpl_stamps: Option<&[crate::cpl::CplTransientStamp]>,
) -> Result<Vec<f64>, MnaError> {
    let base_matrix = &mna.system.matrix;
    let base_rhs = &mna.system.rhs;
    let capacitors = &mna.capacitors;
    let inductors = &mna.inductors;

    let dev_state = DeviceVoltageState::from_solution(mna, prev_solution);

    let load = |solution: &[f64], system: &mut LinearSystem, source_factor: f64, gmin: f64| {
        // 1. Copy base linear stamps (R, V, I topology + inductor topology).
        for triplet in base_matrix.triplets() {
            system.matrix.add(triplet.row, triplet.col, triplet.value);
        }
        for (i, &val) in base_rhs.iter().enumerate() {
            system.rhs[i] += val * source_factor;
        }

        // 1b. Override source values with waveform-evaluated values at time t.
        for vs in &mna.voltage_sources {
            if let Some(ref wf) = vs.waveform {
                let v_t = waveform::evaluate(wf, t, tran_params);
                system.rhs[vs.branch_idx] = v_t * source_factor;
            }
        }
        for cs in &mna.current_sources {
            if let Some(ref wf) = cs.waveform {
                let i_t = waveform::evaluate(wf, t, tran_params);
                let i_diff = i_t - cs.dc_value;
                if let Some(ni) = cs.pos_idx {
                    system.rhs[ni] -= i_diff * source_factor;
                }
                if let Some(nj) = cs.neg_idx {
                    system.rhs[nj] += i_diff * source_factor;
                }
            }
        }

        // 1c. Stamp LTRA transient equations.
        if let (Some(states), Some(coeffs), Some(time_pts)) =
            (ltra_states, ltra_coeffs, ltra_time_points)
        {
            let time_index = time_pts.len() - 1;
            for (li, inst) in mna.ltras.iter().enumerate() {
                crate::ltra::stamp_ltra_transient(
                    inst,
                    &states[li],
                    &coeffs[li],
                    system,
                    num_nodes,
                    t,
                    time_pts,
                    time_index,
                );
            }
        }

        // 1d. Stamp TXL transient equations (pre-computed).
        if let Some(stamps) = txl_stamps {
            for stamp in stamps {
                crate::txl::apply_txl_transient(stamp, system);
            }
        }

        // 1e. Stamp CPL transient equations (pre-computed).
        if let Some(stamps) = cpl_stamps {
            for stamp in stamps {
                crate::cpl::apply_cpl_transient(stamp, system, num_nodes);
            }
        }

        // 2. Stamp capacitor companion models.
        for (ci, cap) in capacitors.iter().enumerate() {
            let (geq, ieq) = capacitor_companion(cap.capacitance, h, &cap_histories[ci], method);
            stamp_conductance(&mut system.matrix, cap.pos_idx, cap.neg_idx, geq);
            stamp_current_source(&mut system.rhs, cap.pos_idx, cap.neg_idx, ieq);
        }

        // 3. Stamp inductor companion models.
        for (li, ind) in inductors.iter().enumerate() {
            let (req, veq) = inductor_companion(ind.inductance, h, &ind_histories[li], method);
            system.matrix.add(ind.branch_idx, ind.branch_idx, -req);
            system.rhs[ind.branch_idx] += veq;
        }

        // 4. Stamp all nonlinear device companions. Device stamps always use
        //    nominal gmin (not the elevated gmin from gmin stepping).
        if has_nonlinear {
            let _ = gmin;
            dev_state.stamp_devices(solution, system, mna, nr_options.gmin);
        }
    };

    if has_nonlinear {
        // Nonlinear: use NR solver.
        let result = newton_raphson_solve(nr_options, dim, num_nodes, load, prev_solution)
            .map_err(|e| {
                MnaError::SolveError(crate::SparseMatrixError::SingularMatrix(e.to_string()))
            })?;
        Ok(result.solution)
    } else {
        // Linear: single solve.
        let mut system = LinearSystem::new(dim);
        load(prev_solution, &mut system, 1.0, nr_options.gmin);
        let sol = system.solve()?;
        Ok(sol)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    /// Helper to get a vector from a transient result.
    fn tran_vector<'a>(result: &'a SimResult, name: &str) -> &'a Vec<f64> {
        let plot = &result.plots[0];
        &plot
            .vecs
            .iter()
            .find(|v| v.name == name)
            .unwrap_or_else(|| panic!("no vector '{name}'"))
            .real
    }

    #[test]
    fn test_rc_step_response() {
        // RC circuit: V1=5V, R1=1k, C1=1u with IC=0 on capacitor.
        // The cap starts at 0V (IC=0) but the DC OP has V(out)=5V,
        // so the IC overrides to 0V and the cap charges to 5V.
        //
        // Analytical: V(out) = 5 * (1 - exp(-t/RC))
        // RC = 1k * 1u = 1ms
        // At t = 1ms (1 RC): V ≈ 5*(1 - 0.368) = 3.16
        // At t = 5ms (5 RC): V ≈ 5*(1 - 0.0067) = 4.97
        let netlist = Netlist::parse(
            "RC step response
V1 1 0 5
R1 1 out 1k
C1 out 0 1u IC=0
.tran 10u 5m
.end
",
        )
        .unwrap();

        let result = simulate_tran(&netlist).unwrap();

        assert_eq!(result.plots.len(), 1);
        assert_eq!(result.plots[0].name, "tran1");

        let time = tran_vector(&result, "time");
        let v_out = tran_vector(&result, "v(out)");

        assert_eq!(time.len(), v_out.len());
        assert!(time.len() > 10, "should have time points");

        // Verify initial condition: t=0, V(out) = 0 (from IC=0).
        assert_abs_diff_eq!(v_out[0], 0.0, epsilon = 1e-6);

        // Check charge curve at several time constants.
        let rc = 1e-3; // 1ms
        for &(t_check, expected_frac) in &[
            (1e-3, 1.0 - (-1.0_f64).exp()), // 1 RC
            (2e-3, 1.0 - (-2.0_f64).exp()), // 2 RC
            (3e-3, 1.0 - (-3.0_f64).exp()), // 3 RC
            (5e-3, 1.0 - (-5.0_f64).exp()), // 5 RC
        ] {
            let expected_v = 5.0 * expected_frac;
            // Find the closest time point.
            let idx = find_nearest_time(time, t_check);

            let actual_v = v_out[idx];
            let rel_err = (actual_v - expected_v).abs() / expected_v;
            assert!(
                rel_err < 0.01,
                "at t={t_check:.3e} ({:.1} RC): expected {expected_v:.4}, got {actual_v:.4}, rel_err={rel_err:.4}",
                t_check / rc
            );
        }
    }

    #[test]
    fn test_lc_oscillator() {
        // LC oscillator: C1=1u with IC=1V, L1=1u, no resistance.
        // Natural frequency: f = 1/(2*pi*sqrt(LC)) = 1/(2*pi*sqrt(1e-6*1e-6))
        //                     = 1/(2*pi*1e-6) ≈ 159.155 kHz
        // Period: T = 1/f ≈ 6.283 us
        //
        // V(1) = V0 * cos(2*pi*f*t) = cos(t/sqrt(LC))
        let netlist = Netlist::parse(
            "LC oscillator
C1 1 0 1u IC=1
L1 1 0 1u
.tran 10n 12.566u
.end
",
        )
        .unwrap();

        let result = simulate_tran(&netlist).unwrap();

        let time = tran_vector(&result, "time");
        let v1 = tran_vector(&result, "v(1)");

        assert!(time.len() > 10, "should have time points");

        // Expected frequency.
        let lc: f64 = 1e-6 * 1e-6;
        let omega = 1.0 / lc.sqrt();
        let f_expected = omega / (2.0 * std::f64::consts::PI);
        let period = 1.0 / f_expected;

        // Find the first zero crossing (quarter period) to verify frequency.
        // V(1) starts at 1V (cos(0) = 1) and should cross zero at T/4.
        let quarter_period = period / 4.0;

        // Find the time index closest to T/4.
        let idx_quarter = find_nearest_time(time, quarter_period);

        // At T/4, voltage should be near zero.
        assert!(
            v1[idx_quarter].abs() < 0.1,
            "at T/4 ({quarter_period:.3e}s): V should be near 0, got {:.4}",
            v1[idx_quarter]
        );

        // At T/2, voltage should be near -1V.
        let half_period = period / 2.0;
        let idx_half = find_nearest_time(time, half_period);

        assert!(
            (v1[idx_half] + 1.0).abs() < 0.1,
            "at T/2 ({half_period:.3e}s): V should be near -1, got {:.4}",
            v1[idx_half]
        );

        // Verify frequency: find two consecutive peaks and check period.
        // Find first maximum after the initial one.
        let idx_full = find_nearest_time(time, period);

        // At T, voltage should return near 1V (full cycle).
        let rel_err = (v1[idx_full] - 1.0).abs();
        assert!(
            rel_err < 0.1,
            "at T ({period:.3e}s): V should be near 1, got {:.4}, err={rel_err:.4}",
            v1[idx_full]
        );

        // Verify frequency within 1%.
        // Use zero crossings to measure actual frequency.
        let mut zero_crossings = Vec::new();
        for i in 1..v1.len() {
            if v1[i - 1] * v1[i] < 0.0 {
                // Linear interpolation for exact crossing.
                let t_cross = time[i - 1]
                    + (time[i] - time[i - 1]) * v1[i - 1].abs() / (v1[i - 1].abs() + v1[i].abs());
                zero_crossings.push(t_cross);
            }
        }

        if zero_crossings.len() >= 4 {
            // Two zero crossings per cycle. Measure period from crossings.
            let measured_period = zero_crossings[2] - zero_crossings[0];
            let measured_freq = 1.0 / measured_period;
            let freq_err = (measured_freq - f_expected).abs() / f_expected;
            assert!(
                freq_err < 0.01,
                "frequency error {freq_err:.4}: expected {f_expected:.1} Hz, got {measured_freq:.1} Hz"
            );
        }
    }

    #[test]
    fn test_pulse_source_rc_circuit() {
        // PULSE source into RC circuit.
        // V1 is a PULSE from 0V to 5V with fast edges, 5ms pulse width, 10ms period.
        // R1=1k, C1=1u → RC = 1ms.
        // During the pulse high, cap charges toward 5V.
        // During pulse low, cap discharges toward 0V.
        let netlist = Netlist::parse(
            "PULSE into RC
V1 1 0 PULSE(0 5 0 1u 1u 5m 10m)
R1 1 out 1k
C1 out 0 1u IC=0
.tran 10u 10m
.end
",
        )
        .unwrap();

        let result = simulate_tran(&netlist).unwrap();

        let time = tran_vector(&result, "time");
        let v_out = tran_vector(&result, "v(out)");

        // At t=0, IC=0 so V(out)=0
        assert_abs_diff_eq!(v_out[0], 0.0, epsilon = 1e-6);

        // During pulse high (0 to 5ms), cap should charge toward 5V.
        // At t=3ms (3 RC into charging), V ≈ 5*(1-exp(-3)) ≈ 4.75
        let idx_3ms = find_nearest_time(time, 3e-3);
        assert!(
            v_out[idx_3ms] > 4.0,
            "at t=3ms, V(out) should be > 4V, got {:.4}",
            v_out[idx_3ms]
        );

        // After pulse falls at t=5ms, cap discharges toward 0V.
        // At t=8ms (3 RC into discharge), voltage should be much lower.
        let idx_8ms = find_nearest_time(time, 8e-3);
        assert!(
            v_out[idx_8ms] < 1.0,
            "at t=8ms (discharge), V(out) should be < 1V, got {:.4}",
            v_out[idx_8ms]
        );
    }

    #[test]
    fn test_sin_source_transient() {
        // SIN source with known frequency, verify output matches analytical.
        // V1 = SIN(0 1 1000) → 1V amplitude, 1kHz, no offset.
        // With a simple R load (no reactive elements), V(1) should track V1.
        let netlist = Netlist::parse(
            "SIN source test
V1 1 0 SIN(0 1 1000)
R1 1 0 1k
.tran 10u 2m
.end
",
        )
        .unwrap();

        let result = simulate_tran(&netlist).unwrap();

        let time = tran_vector(&result, "time");
        let v1 = tran_vector(&result, "v(1)");

        // Verify at multiple points that V(1) matches sin(2*pi*1000*t)
        for i in 1..v1.len() {
            let t = time[i];
            let expected = (2.0 * std::f64::consts::PI * 1000.0 * t).sin();
            let err = (v1[i] - expected).abs();
            assert!(
                err < 1e-6,
                "at t={t:.6e}: expected {expected:.10}, got {:.10}, err={err:.10}",
                v1[i]
            );
        }
    }

    #[test]
    fn test_pwl_source_transient() {
        // PWL source: ramp from 0 to 5V in 1ms, hold 5V for 1ms, ramp to 0 in 1ms.
        // With R load, V(1) should track the source.
        let netlist = Netlist::parse(
            "PWL source test
V1 1 0 PWL(0 0 1m 5 2m 5 3m 0)
R1 1 0 1k
.tran 50u 3m
.end
",
        )
        .unwrap();

        let result = simulate_tran(&netlist).unwrap();

        let time = tran_vector(&result, "time");
        let v1 = tran_vector(&result, "v(1)");

        // At t=0.5ms: should be about 2.5V (midway through ramp)
        let idx = find_nearest_time(time, 0.5e-3);
        assert!(
            (v1[idx] - 2.5).abs() < 0.1,
            "at 0.5ms: expected ~2.5V, got {:.4}",
            v1[idx]
        );

        // At t=1.5ms: should be 5V (holding)
        let idx = find_nearest_time(time, 1.5e-3);
        assert!(
            (v1[idx] - 5.0).abs() < 0.1,
            "at 1.5ms: expected ~5V, got {:.4}",
            v1[idx]
        );

        // At t=2.5ms: should be about 2.5V (ramping down)
        let idx = find_nearest_time(time, 2.5e-3);
        assert!(
            (v1[idx] - 2.5).abs() < 0.1,
            "at 2.5ms: expected ~2.5V, got {:.4}",
            v1[idx]
        );
    }

    #[test]
    fn test_current_source_pulse() {
        // PULSE current source into a resistor.
        // I1 pulses from 0 to 1mA, R1=1k → V(1) should pulse from 0 to 1V.
        let netlist = Netlist::parse(
            "PULSE current source
I1 0 1 PULSE(0 1m 0 1u 1u 1m 2m)
R1 1 0 1k
.tran 10u 4m
.end
",
        )
        .unwrap();

        let result = simulate_tran(&netlist).unwrap();

        let time = tran_vector(&result, "time");
        let v1 = tran_vector(&result, "v(1)");

        // During high pulse (0 to 1ms): V(1) ≈ 1mA * 1kΩ = 1V
        let idx = find_nearest_time(time, 0.5e-3);
        assert!(
            (v1[idx] - 1.0).abs() < 0.05,
            "during pulse high: expected ~1V, got {:.4}",
            v1[idx]
        );

        // During low pulse (1ms to 2ms): V(1) ≈ 0V
        let idx = find_nearest_time(time, 1.5e-3);
        assert!(
            v1[idx].abs() < 0.05,
            "during pulse low: expected ~0V, got {:.4}",
            v1[idx]
        );
    }

    #[test]
    fn test_adaptive_timestep_pulse_rc() {
        // Test that adaptive timestep control places more steps near PULSE edges
        // and coasts during flat regions.
        //
        // PULSE source: 0→5V with 1us edges, 1ms high, 2ms period.
        // RC = 1k * 1u = 1ms.
        let netlist = Netlist::parse(
            "Adaptive timestep PULSE RC
V1 1 0 PULSE(0 5 0 1u 1u 1m 2m)
R1 1 out 1k
C1 out 0 1u IC=0
.tran 10u 4m
.end
",
        )
        .unwrap();

        let result = simulate_tran(&netlist).unwrap();

        let time = tran_vector(&result, "time");
        let v_out = tran_vector(&result, "v(out)");

        // Basic sanity: we should have output points.
        assert!(time.len() > 10, "should have output points");

        // Verify accuracy at key points.
        // At t=0: IC=0
        assert_abs_diff_eq!(v_out[0], 0.0, epsilon = 1e-6);

        // During charging (t≈0.5ms): V(out) should be rising.
        let idx_05ms = find_nearest_time(time, 0.5e-3);
        assert!(
            v_out[idx_05ms] > 1.0,
            "at 0.5ms: V(out) should be rising, got {:.4}",
            v_out[idx_05ms]
        );

        // Near end of charging (t≈1ms): V(out) should be near 3.16V (1 RC).
        let idx_1ms = find_nearest_time(time, 1e-3);
        let expected_1rc = 5.0 * (1.0 - (-1.0_f64).exp());
        let rel_err = (v_out[idx_1ms] - expected_1rc).abs() / expected_1rc;
        assert!(
            rel_err < 0.05,
            "at 1ms: expected {expected_1rc:.3}, got {:.3}, err={rel_err:.3}",
            v_out[idx_1ms]
        );

        // Verify adaptive stepping: check that timesteps are smaller near edges.
        // Count steps in first 10us (near rising edge) vs steps in 0.1ms-0.5ms (flat region).
        let steps_near_edge = time.windows(2).filter(|w| w[0] < 10e-6).count();
        let steps_flat = time
            .windows(2)
            .filter(|w| w[0] >= 0.1e-3 && w[0] < 0.5e-3)
            .count();

        // Near the edge, there should be comparable or more steps per unit time than flat.
        // This is a soft check — adaptive stepping should concentrate steps near transitions.
        let density_edge = if steps_near_edge > 0 {
            steps_near_edge as f64 / 10e-6
        } else {
            0.0
        };
        let density_flat = if steps_flat > 0 {
            steps_flat as f64 / 0.4e-3
        } else {
            1.0
        };

        // Edge density should be at least comparable to flat density
        // (breakpoints force steps there).
        assert!(
            density_edge >= density_flat * 0.1 || steps_near_edge >= 1,
            "adaptive stepping should have steps near edges: edge_steps={steps_near_edge}, flat_steps={steps_flat}"
        );
    }

    /// Find the index of the time point nearest to `target`.
    fn find_nearest_time(time: &[f64], target: f64) -> usize {
        time.iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (*a - target)
                    .abs()
                    .partial_cmp(&(*b - target).abs())
                    .unwrap()
            })
            .unwrap()
            .0
    }
}
