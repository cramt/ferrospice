//! Transient (time-domain) analysis engine.
//!
//! Implements `.tran` analysis with Backward Euler (BE) and Trapezoidal (Trap)
//! integration methods. Capacitors and inductors are converted to companion
//! models (conductance/resistance + current/voltage source) at each timestep.

use ferrospice_netlist::{Analysis, Expr, Item, Netlist, SimPlot, SimResult, SimVector};

use crate::LinearSystem;
use crate::diode::{VT_NOM, pnjlim, vcrit};
use crate::mna::{MnaError, MnaSystem, assemble_mna, stamp_conductance};
use crate::newton::{NrOptions, newton_raphson_solve};
use crate::simulate::simulate_op;
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

/// Stamp a current source into the RHS vector.
fn stamp_current_source(rhs: &mut [f64], ni: Option<usize>, nj: Option<usize>, i_val: f64) {
    if let Some(i) = ni {
        rhs[i] -= i_val;
    }
    if let Some(j) = nj {
        rhs[j] += i_val;
    }
}

/// Perform transient analysis on a circuit.
///
/// Parses the `.tran` command from the netlist, computes the DC operating point
/// as initial conditions, then steps through time using numerical integration.
pub fn simulate_tran(netlist: &Netlist) -> Result<SimResult, MnaError> {
    // Find the .tran analysis command.
    let (tstep, tstop, tstart, _tmax) = netlist
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

    let h = expr_val(&tstep, ".tran tstep")?;
    let t_stop = expr_val(&tstop, ".tran tstop")?;
    let t_start = tstart
        .as_ref()
        .map(|e| expr_val(e, ".tran tstart"))
        .transpose()?
        .unwrap_or(0.0);

    if h <= 0.0 || t_stop <= 0.0 {
        return Err(MnaError::UnsupportedElement(
            "invalid .tran parameters".to_string(),
        ));
    }

    // Assemble MNA system.
    let mna = assemble_mna(netlist)?;

    // Compute DC operating point for initial conditions.
    let op_result = simulate_op(netlist)?;
    let op_plot = &op_result.plots[0];

    // Extract initial solution vector from DC OP.
    let num_ext_nodes = mna.node_map.len();
    let num_internal = mna
        .diodes
        .iter()
        .filter(|d| d.internal_idx.is_some())
        .count();
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
                    // Both non-ground: set positive node relative to negative.
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

    // Precompute diode vcrits for NR.
    let vcrits: Vec<f64> = mna
        .diodes
        .iter()
        .map(|d| vcrit(d.model.n * VT_NOM, d.model.is))
        .collect();

    let has_nonlinear = !mna.diodes.is_empty();
    let nr_options = NrOptions::default();
    let tran_params = TranParams {
        tstep: h,
        tstop: t_stop,
    };

    // Time-stepping loop.
    while t < t_stop - h * 1e-9 {
        t += h;

        // Use Backward Euler for the first step, then Trapezoidal.
        let method = if time_vec.real.len() <= 1 {
            IntegrationMethod::BackwardEuler
        } else {
            IntegrationMethod::Trapezoidal
        };

        // Solve this timestep.
        solution = solve_timestep(
            &mna,
            &solution,
            h,
            t,
            &tran_params,
            method,
            &cap_histories,
            &ind_histories,
            &vcrits,
            has_nonlinear,
            &nr_options,
            dim,
            num_nodes,
        )?;

        // Update histories from new solution.
        for (ci, cap) in mna.capacitors.iter().enumerate() {
            let v_pos = cap.pos_idx.map(|i| solution[i]).unwrap_or(0.0);
            let v_neg = cap.neg_idx.map(|i| solution[i]).unwrap_or(0.0);
            let v_new = v_pos - v_neg;

            // Compute current through capacitor.
            let current = match method {
                IntegrationMethod::BackwardEuler => {
                    let geq = cap.capacitance / h;
                    geq * (v_new - cap_histories[ci].voltage)
                }
                IntegrationMethod::Trapezoidal => {
                    let geq = 2.0 * cap.capacitance / h;
                    geq * (v_new - cap_histories[ci].voltage) - cap_histories[ci].current
                }
            };

            cap_histories[ci] = CapHistory {
                voltage: v_new,
                current,
            };
        }

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
    vcrits: &[f64],
    has_nonlinear: bool,
    nr_options: &NrOptions,
    dim: usize,
    num_nodes: usize,
) -> Result<Vec<f64>, MnaError> {
    let base_matrix = &mna.system.matrix;
    let base_rhs = &mna.system.rhs;
    let diodes = &mna.diodes;
    let capacitors = &mna.capacitors;
    let inductors = &mna.inductors;

    // Track previous junction voltages for diode pnjlim.
    let prev_jct_voltages = std::cell::RefCell::new(
        diodes
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
            .collect::<Vec<f64>>(),
    );

    let load = |solution: &[f64], system: &mut LinearSystem, source_factor: f64| {
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
                // Replace the DC value in the branch equation RHS.
                system.rhs[vs.branch_idx] = v_t * source_factor;
            }
        }
        for cs in &mna.current_sources {
            if let Some(ref wf) = cs.waveform {
                let i_t = waveform::evaluate(wf, t, tran_params);
                // Undo DC stamp and apply transient value.
                let i_diff = i_t - cs.dc_value;
                if let Some(ni) = cs.pos_idx {
                    system.rhs[ni] -= i_diff * source_factor;
                }
                if let Some(nj) = cs.neg_idx {
                    system.rhs[nj] += i_diff * source_factor;
                }
            }
        }

        // 2. Stamp capacitor companion models.
        for (ci, cap) in capacitors.iter().enumerate() {
            let (geq, ieq) = capacitor_companion(cap.capacitance, h, &cap_histories[ci], method);
            stamp_conductance(&mut system.matrix, cap.pos_idx, cap.neg_idx, geq);
            // ieq is the history current source (flows from pos to neg).
            stamp_current_source(&mut system.rhs, cap.pos_idx, cap.neg_idx, ieq);
        }

        // 3. Stamp inductor companion models.
        for (li, ind) in inductors.iter().enumerate() {
            let (req, veq) = inductor_companion(ind.inductance, h, &ind_histories[li], method);
            // Add -req on the branch equation diagonal.
            system.matrix.add(ind.branch_idx, ind.branch_idx, -req);
            // Add veq to the branch equation RHS.
            system.rhs[ind.branch_idx] += veq;
        }

        // 4. Stamp diode companion models (same as DC NR).
        if has_nonlinear {
            let mut prev = prev_jct_voltages.borrow_mut();
            for (di, diode) in diodes.iter().enumerate() {
                let (jct_anode, jct_cathode) = if diode.internal_idx.is_some() {
                    (diode.internal_idx, diode.cathode_idx)
                } else {
                    (diode.anode_idx, diode.cathode_idx)
                };

                let v_anode = jct_anode.map(|i| solution[i]).unwrap_or(0.0);
                let v_cathode = jct_cathode.map(|i| solution[i]).unwrap_or(0.0);
                let mut v_jct = v_anode - v_cathode;

                v_jct = pnjlim(v_jct, prev[di], diode.model.n * VT_NOM, vcrits[di]);
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
        load(prev_solution, &mut system, 1.0);
        let sol = system.solve()?;
        Ok(sol)
    }
}

fn expr_val(expr: &Expr, context: &str) -> Result<f64, MnaError> {
    match expr {
        Expr::Num(v) => Ok(*v),
        _ => Err(MnaError::NonNumericValue {
            element: context.to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

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
        assert!(time.len() > 100, "should have many time points");

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
            let idx = time
                .iter()
                .position(|&t| (t - t_check).abs() < 1e-6)
                .unwrap_or_else(|| {
                    // Find nearest.
                    time.iter()
                        .enumerate()
                        .min_by(|(_, a), (_, b)| {
                            (*a - t_check)
                                .abs()
                                .partial_cmp(&(*b - t_check).abs())
                                .unwrap()
                        })
                        .unwrap()
                        .0
                });

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

        assert!(time.len() > 100, "should have many time points");

        // Expected frequency.
        let lc: f64 = 1e-6 * 1e-6;
        let omega = 1.0 / lc.sqrt();
        let f_expected = omega / (2.0 * std::f64::consts::PI);
        let period = 1.0 / f_expected;

        // Find the first zero crossing (quarter period) to verify frequency.
        // V(1) starts at 1V (cos(0) = 1) and should cross zero at T/4.
        let quarter_period = period / 4.0;

        // Find the time index closest to T/4.
        let idx_quarter = time
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (*a - quarter_period)
                    .abs()
                    .partial_cmp(&(*b - quarter_period).abs())
                    .unwrap()
            })
            .unwrap()
            .0;

        // At T/4, voltage should be near zero.
        assert!(
            v1[idx_quarter].abs() < 0.1,
            "at T/4 ({quarter_period:.3e}s): V should be near 0, got {:.4}",
            v1[idx_quarter]
        );

        // At T/2, voltage should be near -1V.
        let half_period = period / 2.0;
        let idx_half = time
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (*a - half_period)
                    .abs()
                    .partial_cmp(&(*b - half_period).abs())
                    .unwrap()
            })
            .unwrap()
            .0;

        assert!(
            (v1[idx_half] + 1.0).abs() < 0.1,
            "at T/2 ({half_period:.3e}s): V should be near -1, got {:.4}",
            v1[idx_half]
        );

        // Verify frequency: find two consecutive peaks and check period.
        // Find first maximum after the initial one.
        let idx_full = time
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (*a - period)
                    .abs()
                    .partial_cmp(&(*b - period).abs())
                    .unwrap()
            })
            .unwrap()
            .0;

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
