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
    pub itl2: usize,
    /// Minimum conductance from each node to ground (ngspice GMIN, default 1e-12).
    pub gmin: f64,
}

impl Default for NrOptions {
    fn default() -> Self {
        Self {
            abstol: 1e-12,
            reltol: 1e-3,
            vntol: 1e-6,
            itl1: 100,
            itl2: 50,
            gmin: 1e-12,
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
    F: Fn(&[f64], &mut LinearSystem, f64),
{
    let mut solution = initial_guess.to_vec();
    let mut system = LinearSystem::new(dim);
    let mut total_iters = 0;

    for iter in 0..attempt.max_iters {
        system.matrix.clear();
        system.rhs.fill(0.0);
        load_system(&solution, &mut system, attempt.source_factor);

        // Add Gmin from each node to ground for numerical stability.
        for i in 0..num_nodes {
            system.matrix.add(i, i, attempt.gmin);
        }

        let new_solution = system.solve()?;
        total_iters = iter + 1;

        // Check convergence (skip first iteration — need two solutions to compare).
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
    F: Fn(&[f64], &mut LinearSystem, f64),
{
    let mut gmin = 1e-2;
    let gmin_factor = 10.0;
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

/// Source stepping fallback.
///
/// Ramp all independent sources from 0 to full value in steps,
/// converging at each step. The `load_system` callback receives
/// a `source_factor` (0.0 to 1.0) that scales independent sources.
///
/// Matches ngspice SPICE3-style source stepping from `cktop.c`.
fn source_stepping<F>(
    options: &NrOptions,
    dim: usize,
    num_nodes: usize,
    load_system: &F,
    initial_guess: &[f64],
) -> Result<NrResult, NrError>
where
    F: Fn(&[f64], &mut LinearSystem, f64),
{
    let num_steps = 10;
    let mut solution = initial_guess.to_vec();
    let mut total_iters = 0;

    for step in 0..=num_steps {
        let factor = step as f64 / num_steps as f64;
        let attempt = NrAttempt {
            gmin: options.gmin,
            source_factor: factor,
            max_iters: options.itl2,
        };
        let result = try_nr(options, dim, num_nodes, load_system, &solution, &attempt)?;
        total_iters += result.iterations;
        solution = result.solution;
    }

    Ok(NrResult {
        solution,
        iterations: total_iters,
        converged: true,
    })
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
///   estimate. Called as `load_system(solution, system, source_factor)` where
///   `source_factor` is 1.0 for normal operation and < 1.0 during source stepping.
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
    F: Fn(&[f64], &mut LinearSystem, f64),
{
    // Try direct NR first.
    let attempt = NrAttempt {
        gmin: options.gmin,
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

        let load = |solution: &[f64], system: &mut LinearSystem, source_factor: f64| {
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
        let load = |_solution: &[f64], system: &mut LinearSystem, source_factor: f64| {
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

        let load = |solution: &[f64], system: &mut LinearSystem, source_factor: f64| {
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
