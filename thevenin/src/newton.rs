use crate::{LinearSystem, SparseMatrixError};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum NrError {
    #[error("Newton-Raphson failed to converge after {iterations} iterations")]
    NoConvergence { iterations: usize },
    #[error("failed to solve linear system: {0}")]
    SolveError(#[from] SparseMatrixError),
}

/// Newton-Raphson convergence options matching ngspice defaults.
///
/// See ngspice `src/spicelib/devices/cktinit.c` for default values.
#[derive(Debug, Clone)]
pub struct NrOptions {
    /// Absolute current tolerance (ngspice ABSTOL, default 1e-12).
    pub abstol: f64,
    /// Relative tolerance (ngspice RELTOL, default 1e-3).
    pub reltol: f64,
    /// Absolute voltage tolerance (ngspice VNTOL, default 1e-6).
    pub vntol: f64,
    /// Maximum iterations for DC operating point (ngspice ITL1, default 100).
    pub itl1: usize,
    /// Maximum iterations per step during Gmin/source stepping (ngspice ITL2, default 50).
    /// We use 200 (vs ngspice's 50) because slow-converging circuits (e.g. CMOS drivers
    /// with LAMBDA=0 MOSFETs feeding LTRA floating nodes) need more iterations to damp NR
    /// oscillations at each Gmin step.
    pub itl2: usize,
    /// Minimum conductance from each node to ground (ngspice GMIN, default 1e-12).
    /// Used by device models in junction conductance computations.
    pub gmin: f64,
    /// Diagonal Gmin added from every node to ground by the solver
    /// (ngspice `CKTdiagGmin`).  In ngspice this starts at 0 and is only
    /// elevated during Gmin stepping.  We default to `gmin` for backward
    /// compatibility; the DC sweep code sets it to 0 so that the only Gmin
    /// on device nodes comes from the device model equations (matching
    /// ngspice behaviour where `CKTdiagGmin` stays 0 when the initial
    /// NIiter converges).
    pub diag_gmin: f64,
}

impl Default for NrOptions {
    fn default() -> Self {
        Self {
            abstol: 1e-12,
            reltol: 1e-3,
            vntol: 1e-6,
            itl1: 100,
            itl2: 200,
            gmin: 1e-12,
            diag_gmin: 1e-12,
        }
    }
}

/// Result of Newton-Raphson iteration.
#[derive(Debug)]
pub struct NrResult {
    /// Final solution vector.
    pub solution: Vec<f64>,
    /// Number of iterations performed (total across all attempts).
    pub iterations: usize,
    /// Whether convergence was achieved.
    pub converged: bool,
}

/// Check convergence of NR iteration using ngspice-style criteria.
///
/// For node voltages (indices 0..num_nodes):
///   |v_new - v_old| <= reltol * max(|v_new|, |v_old|) + vntol
///
/// For branch currents (indices num_nodes..):
///   |i_new - i_old| <= reltol * max(|i_new|, |i_old|) + abstol
fn check_convergence(old: &[f64], new: &[f64], num_nodes: usize, options: &NrOptions) -> bool {
    for i in 0..old.len() {
        let diff = (new[i] - old[i]).abs();
        let tol = if i < num_nodes {
            options.reltol * new[i].abs().max(old[i].abs()) + options.vntol
        } else {
            options.reltol * new[i].abs().max(old[i].abs()) + options.abstol
        };
        if diff > tol {
            return false;
        }
    }
    true
}

/// Parameters for a single NR attempt.
struct NrAttempt {
    gmin: f64,
    source_factor: f64,
    max_iters: usize,
}

/// Run NR iteration with given attempt parameters.
///
/// Returns `Ok(NrResult)` if converged, `Err` otherwise.
fn try_nr<F>(
    options: &NrOptions,
    dim: usize,
    num_nodes: usize,
    load_system: &F,
    initial_guess: &[f64],
    attempt: &NrAttempt,
) -> Result<NrResult, NrError>
where
    F: Fn(&[f64], &mut LinearSystem, f64, f64),
{
    let mut solution = initial_guess.to_vec();
    let mut system = LinearSystem::new(dim);
    let mut total_iters = 0;
    for iter in 0..attempt.max_iters {
        system.matrix.clear();
        system.rhs.fill(0.0);
        load_system(&solution, &mut system, attempt.source_factor, attempt.gmin);

        // Add Gmin from each node to ground for numerical stability.
        for i in 0..num_nodes {
            system.matrix.add(i, i, attempt.gmin);
        }

        let new_solution = system.solve()?;
        total_iters = iter + 1;

        if iter > 0 && check_convergence(&solution, &new_solution, num_nodes, options) {
            return Ok(NrResult {
                solution: new_solution,
                iterations: total_iters,
                converged: true,
            });
        }

        solution = new_solution;
    }

    Err(NrError::NoConvergence {
        iterations: total_iters,
    })
}

/// Gmin stepping fallback.
///
/// Start with an elevated Gmin (1e-2), converge, then reduce Gmin by 10x
/// each step until reaching the target Gmin. Uses the previous step's
/// solution as initial guess for the next step.
///
/// Matches ngspice SPICE3-style Gmin stepping from `cktop.c`.
fn gmin_stepping<F>(
    options: &NrOptions,
    dim: usize,
    num_nodes: usize,
    load_system: &F,
    initial_guess: &[f64],
) -> Result<NrResult, NrError>
where
    F: Fn(&[f64], &mut LinearSystem, f64, f64),
{
    let mut gmin = 1e-2;
    let gmin_factor = 10.0; // matches ngspice SPICE3-style gmin stepping
    let mut solution = initial_guess.to_vec();
    let mut total_iters = 0;

    while gmin >= options.gmin * 0.9 {
        let attempt = NrAttempt {
            gmin,
            source_factor: 1.0,
            max_iters: options.itl2,
        };
        let result = try_nr(options, dim, num_nodes, load_system, &solution, &attempt)?;
        total_iters += result.iterations;
        solution = result.solution;
        gmin /= gmin_factor;
    }

    // Final solve with target Gmin.
    let attempt = NrAttempt {
        gmin: options.gmin,
        source_factor: 1.0,
        max_iters: options.itl2,
    };
    let result = try_nr(options, dim, num_nodes, load_system, &solution, &attempt)?;
    total_iters += result.iterations;

    Ok(NrResult {
        solution: result.solution,
        iterations: total_iters,
        converged: true,
    })
}

/// Source stepping fallback with elevated Gmin.
///
/// Ramps all independent sources from 0 to full value while keeping an
/// elevated diagonal Gmin (1e-2) to prevent near-singular matrices from
/// floating nodes (e.g. LTRA ports connected only through reactive elements).
/// After all source steps converge, reduces Gmin to the target value using
/// the same schedule as `gmin_stepping`.
///
/// This combined approach handles circuits where plain gmin_stepping diverges
/// (too many iterations needed) and plain source_stepping fails with a
/// singular matrix at Gmin = options.gmin.
fn source_stepping<F>(
    options: &NrOptions,
    dim: usize,
    num_nodes: usize,
    load_system: &F,
    initial_guess: &[f64],
) -> Result<NrResult, NrError>
where
    F: Fn(&[f64], &mut LinearSystem, f64, f64),
{
    // Use elevated Gmin throughout source ramping to prevent near-singular
    // matrices caused by floating nodes (e.g. LTRA-coupled node pairs where
    // the only path to ground is the tiny diagonal Gmin=1e-12).
    let gmin_elevated = 1e-2_f64;
    let num_steps = 10;
    let mut solution = initial_guess.to_vec();
    let mut total_iters = 0;

    // Phase 1: ramp sources with elevated Gmin.
    for step in 0..=num_steps {
        let factor = step as f64 / num_steps as f64;
        let attempt = NrAttempt {
            gmin: gmin_elevated,
            source_factor: factor,
            max_iters: options.itl2,
        };
        let result = try_nr(options, dim, num_nodes, load_system, &solution, &attempt)?;
        total_iters += result.iterations;
        solution = result.solution;
    }

    // Phase 2: reduce Gmin from elevated level to target (gmin_stepping schedule).
    let gmin_factor = 10.0;
    let mut gmin = gmin_elevated / gmin_factor;
    while gmin >= options.gmin * 0.9 {
        let attempt = NrAttempt {
            gmin,
            source_factor: 1.0,
            max_iters: options.itl2,
        };
        let result = try_nr(options, dim, num_nodes, load_system, &solution, &attempt)?;
        total_iters += result.iterations;
        solution = result.solution;
        gmin /= gmin_factor;
    }

    // Final solve at target Gmin.
    let attempt = NrAttempt {
        gmin: options.gmin,
        source_factor: 1.0,
        max_iters: options.itl2,
    };
    let result = try_nr(options, dim, num_nodes, load_system, &solution, &attempt)?;
    total_iters += result.iterations;

    Ok(NrResult {
        solution: result.solution,
        iterations: total_iters,
        converged: true,
    })
}

/// Solve a nonlinear system using source stepping directly, bypassing direct NR
/// and Gmin stepping.
///
/// Use this for circuits with transmission lines (LTRA/TXL) combined with
/// cascaded MOSFET stages.  Without voltage-step limiting, Gmin stepping can
/// converge to spurious negative-voltage fixed points in cascaded-inverter
/// circuits.  Source stepping avoids this by ramping all independent sources
/// from zero, ensuring NR always follows the physical solution trajectory.
pub fn source_stepping_solve<F>(
    options: &NrOptions,
    dim: usize,
    num_nodes: usize,
    load_system: F,
    initial_guess: &[f64],
) -> Result<NrResult, NrError>
where
    F: Fn(&[f64], &mut LinearSystem, f64, f64),
{
    source_stepping(options, dim, num_nodes, &load_system, initial_guess)
}

/// Solve a potentially nonlinear system using Newton-Raphson iteration.
///
/// # Arguments
///
/// * `options` - Convergence parameters (tolerances, iteration limits).
/// * `dim` - System dimension (number of unknowns).
/// * `num_nodes` - Number of node voltage entries (for convergence check;
///   indices 0..num_nodes are voltages, num_nodes.. are branch currents).
/// * `load_system` - Callback that fills the linear system for a given solution
///   estimate. Called as `load_system(solution, system, source_factor, gmin)` where
///   `source_factor` is 1.0 for normal operation and < 1.0 during source stepping,
///   and `gmin` is the current minimum conductance (elevated during Gmin stepping).
///   The system is zeroed before each call.
/// * `initial_guess` - Starting solution vector (length `dim`).
///
/// # Convergence strategy
///
/// 1. Try direct NR iteration (up to ITL1 iterations).
/// 2. If that fails, try Gmin stepping (elevated Gmin, gradually reduced).
/// 3. If that also fails, try source stepping (ramp sources from 0 to full).
pub fn newton_raphson_solve<F>(
    options: &NrOptions,
    dim: usize,
    num_nodes: usize,
    load_system: F,
    initial_guess: &[f64],
) -> Result<NrResult, NrError>
where
    F: Fn(&[f64], &mut LinearSystem, f64, f64),
{
    // Try direct NR first.
    // Use options.diag_gmin for the diagonal Gmin.  For DC OP this is 0
    // (matching ngspice CKTdiagGmin initial value); for transient timesteps
    // it may be elevated to options.gmin for regularization.
    let attempt = NrAttempt {
        gmin: options.diag_gmin,
        source_factor: 1.0,
        max_iters: options.itl1,
    };
    if let Ok(result) = try_nr(
        options,
        dim,
        num_nodes,
        &load_system,
        initial_guess,
        &attempt,
    ) {
        return Ok(result);
    }

    // Fallback 1: Gmin stepping.
    if let Ok(result) = gmin_stepping(options, dim, num_nodes, &load_system, initial_guess) {
        return Ok(result);
    }

    // Fallback 2: Source stepping.
    source_stepping(options, dim, num_nodes, &load_system, initial_guess)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    use approx::assert_abs_diff_eq;

    /// Thermal voltage at room temperature (300K).
    const VT: f64 = 0.02585;
    /// Diode saturation current.
    const IS: f64 = 1e-14;

    /// Diode current: I = IS * (exp(V/Vt) - 1)
    fn diode_current(v: f64) -> f64 {
        IS * (safe_exp(v / VT) - 1.0)
    }

    /// Diode conductance: dI/dV = IS/Vt * exp(V/Vt)
    fn diode_conductance(v: f64) -> f64 {
        IS / VT * safe_exp(v / VT)
    }

    /// Safe exponential that clamps the argument to avoid overflow.
    fn safe_exp(x: f64) -> f64 {
        if x > 500.0 {
            (500.0_f64).exp()
        } else {
            x.exp()
        }
    }

    /// Test NR on a simple diode circuit: V1=1V, R1=1k, D1.
    ///
    /// Circuit: V1(=1V) --- R1(=1k) --- anode(D1) --- cathode(ground)
    ///
    /// MNA system (2 unknowns: V(1)=node 1 voltage, I(V1)=branch current):
    /// Node 1: V1 is connected via branch equation.
    /// Node 2 (diode anode): KCL: (V2-V1)/R1 + I_diode(V2) = 0
    ///
    /// Actually, let's use a simpler setup:
    /// 3 unknowns: V(1), V(2), I(V1)
    /// Node 1: V1 branch current + G*(V1-V2) = 0 => I_branch + G*V1 - G*V2 = 0
    /// Node 2: G*(V2-V1) + I_diode(V2) = 0 => -G*V1 + G*V2 + I_diode(V2) = 0
    /// Branch: V1 - 0 = 1.0 => V1 = 1.0
    ///
    /// Where G = 1/R = 1/1000
    #[test]
    fn test_diode_circuit_convergence() {
        let g = 1.0 / 1000.0; // R = 1k
        let v_source = 1.0;
        let dim = 3; // V(1), V(2), I(V1)
        let num_nodes = 2;
        let options = NrOptions::default();

        let load = |solution: &[f64], system: &mut LinearSystem, source_factor: f64, _gmin: f64| {
            let v2 = solution[1]; // diode voltage

            // Resistor R1 between nodes 0 and 1 (matrix indices)
            system.matrix.add(0, 0, g);
            system.matrix.add(0, 1, -g);
            system.matrix.add(1, 0, -g);
            system.matrix.add(1, 1, g);

            // Voltage source V1: node 0 to ground, branch index 2
            system.matrix.add(0, 2, 1.0);
            system.matrix.add(2, 0, 1.0);
            system.rhs[2] = v_source * source_factor;

            // Diode D1 between node 1 (anode) and ground (cathode)
            // Companion model: I_eq + g_d * V2
            // where g_d = dI/dV at V2, I_eq = I(V2) - g_d * V2
            let g_d = diode_conductance(v2);
            let i_d = diode_current(v2);
            let i_eq = i_d - g_d * v2;

            // Stamp conductance from node 1 to ground
            system.matrix.add(1, 1, g_d);
            // Stamp equivalent current source at node 1
            system.rhs[1] -= i_eq;
        };

        let initial = vec![0.0; dim];
        let result = newton_raphson_solve(&options, dim, num_nodes, load, &initial).unwrap();

        assert!(result.converged);
        assert!(
            result.iterations < 20,
            "expected < 20 iterations, got {}",
            result.iterations
        );

        let v1 = result.solution[0];
        let v_diode = result.solution[1];

        // V1 should be the source voltage
        assert_abs_diff_eq!(v1, v_source, epsilon = 1e-9);

        // Diode forward voltage should be ~0.6-0.7V
        assert!(
            v_diode > 0.5 && v_diode < 0.8,
            "diode voltage {v_diode} not in expected range 0.5-0.8V"
        );

        // Verify KCL at node 2: current through R = current through diode
        let i_r = (v1 - v_diode) * g;
        let i_d = diode_current(v_diode);
        assert_abs_diff_eq!(i_r, i_d, epsilon = 1e-9);
    }

    /// Test that a purely linear system converges in 2 iterations
    /// (first iteration gets the answer, second confirms convergence).
    #[test]
    fn test_linear_system_converges_immediately() {
        let options = NrOptions::default();
        let dim = 2;
        let num_nodes = 1;

        // V1=5V, R1=1k to ground
        // V(1) = 5V, I(V1) = -5mA
        let load =
            |_solution: &[f64], system: &mut LinearSystem, source_factor: f64, _gmin: f64| {
                let g = 1.0 / 1000.0;
                // Resistor from node 0 to ground
                system.matrix.add(0, 0, g);
                // Voltage source: node 0, branch 1
                system.matrix.add(0, 1, 1.0);
                system.matrix.add(1, 0, 1.0);
                system.rhs[1] = 5.0 * source_factor;
            };

        let initial = vec![0.0; dim];
        let result = newton_raphson_solve(&options, dim, num_nodes, load, &initial).unwrap();

        assert!(result.converged);
        assert_eq!(result.iterations, 2); // solve once, confirm on second
        assert_abs_diff_eq!(result.solution[0], 5.0, epsilon = 1e-9);
        assert_abs_diff_eq!(result.solution[1], -5e-3, epsilon = 1e-9);
    }

    /// Test convergence check function directly.
    #[test]
    fn test_convergence_check() {
        let options = NrOptions::default();

        // Identical solutions should converge.
        let a = vec![1.0, 2.0, 0.001];
        assert!(check_convergence(&a, &a, 2, &options));

        // Small voltage change within vntol.
        let b = vec![1.0, 2.0 + 1e-7, 0.001];
        assert!(check_convergence(&a, &b, 2, &options));

        // Large voltage change — should not converge.
        let c = vec![1.0, 2.1, 0.001];
        assert!(!check_convergence(&a, &c, 2, &options));

        // Small current change within abstol.
        let d = vec![1.0, 2.0, 0.001 + 1e-13];
        assert!(check_convergence(&a, &d, 2, &options));

        // Large current change — should not converge.
        let e = vec![1.0, 2.0, 0.002];
        assert!(!check_convergence(&a, &e, 2, &options));
    }

    /// Test that Gmin stepping can solve a circuit that direct NR cannot
    /// easily handle (e.g., diode with high source voltage).
    #[test]
    fn test_diode_high_voltage_converges() {
        let g = 1.0 / 100.0; // R = 100 ohms
        let v_source = 10.0; // Higher voltage
        let dim = 3;
        let num_nodes = 2;
        let options = NrOptions::default();

        let load = |solution: &[f64], system: &mut LinearSystem, source_factor: f64, _gmin: f64| {
            let v2 = solution[1];

            // Resistor
            system.matrix.add(0, 0, g);
            system.matrix.add(0, 1, -g);
            system.matrix.add(1, 0, -g);
            system.matrix.add(1, 1, g);

            // Voltage source
            system.matrix.add(0, 2, 1.0);
            system.matrix.add(2, 0, 1.0);
            system.rhs[2] = v_source * source_factor;

            // Diode companion model
            let g_d = diode_conductance(v2);
            let i_d = diode_current(v2);
            let i_eq = i_d - g_d * v2;
            system.matrix.add(1, 1, g_d);
            system.rhs[1] -= i_eq;
        };

        let initial = vec![0.0; dim];
        let result = newton_raphson_solve(&options, dim, num_nodes, load, &initial).unwrap();

        assert!(result.converged);

        let v_diode = result.solution[1];
        // With 10V source and 100 ohm resistor, diode voltage should still be ~0.6-0.8V
        assert!(
            v_diode > 0.5 && v_diode < 0.9,
            "diode voltage {v_diode} not in expected range"
        );

        // Verify: I_R = I_D (within NR convergence tolerance)
        let i_r = (result.solution[0] - v_diode) * g;
        let i_d = diode_current(v_diode);
        assert_abs_diff_eq!(i_r, i_d, epsilon = 1e-4);
    }
}
