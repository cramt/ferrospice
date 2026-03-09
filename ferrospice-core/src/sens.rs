use faer::Mat;
use faer::linalg::solvers::FullPivLu;
use faer::prelude::Solve;

use ferrospice_netlist::{Analysis, Item, Netlist, SimPlot, SimResult, SimVector};

use crate::LinearSystem;
use crate::bjt::stamp_bjt;
use crate::jfet::stamp_jfet;
use crate::mna::{MnaError, MnaSystem, assemble_mna};
use crate::mosfet::stamp_mosfet;
use crate::simulate::solve_op_raw;
use crate::tf::build_jacobian;

/// Solve a transposed system Y^T * x = rhs.
fn solve_transposed(jacobian: &Mat<f64>, rhs: &[f64]) -> Result<Vec<f64>, MnaError> {
    let dim = rhs.len();
    let jac_t = jacobian.transpose().to_owned();
    let lu = FullPivLu::new(jac_t.as_ref());
    let mut b = Mat::zeros(dim, 1);
    for (i, &val) in rhs.iter().enumerate() {
        b[(i, 0)] = val;
    }
    let x = lu.solve(&b);
    let result: Vec<f64> = (0..dim).map(|i| x[(i, 0)]).collect();
    if result.iter().any(|v: &f64| v.is_nan() || v.is_infinite()) {
        return Err(MnaError::SolveError(
            crate::SparseMatrixError::SingularMatrix("adjoint solve singular".to_string()),
        ));
    }
    Ok(result)
}

/// Compute sensitivity for a conductance parameter change.
///
/// When conductance G between nodes (pos, neg) changes by dG:
/// - The matrix changes: dY has conductance stamp with value dG
/// - delta_residual = -dY * x_op
/// - S = lambda^T * delta_residual = -dG * V_across * (lambda[pos] - lambda[neg])
fn conductance_sensitivity(
    lambda: &[f64],
    pos: Option<usize>,
    neg: Option<usize>,
    dg: f64,
    solution: &[f64],
) -> f64 {
    let v_pos = pos.map(|i| solution[i]).unwrap_or(0.0);
    let v_neg = neg.map(|i| solution[i]).unwrap_or(0.0);
    let v_across = v_pos - v_neg;
    let lambda_pos = pos.map(|i| lambda[i]).unwrap_or(0.0);
    let lambda_neg = neg.map(|i| lambda[i]).unwrap_or(0.0);

    // S = lambda^T * (-dY/dp * x) = -dG * V_across * (lambda[pos] - lambda[neg])
    -dg * v_across * (lambda_pos - lambda_neg)
}

/// Compute sensitivity using numerical perturbation of device stamps.
///
/// For complex nonlinear devices, perturb the parameter and recompute stamps.
/// Uses the residual method: S = lambda^T * (delta_b - delta_Y * x) / delta_p
fn numerical_device_sensitivity(
    lambda: &[f64],
    solution: &[f64],
    _mna: &MnaSystem,
    dim: usize,
    stamp_original: impl Fn(&mut crate::SparseMatrix, &mut [f64]),
    stamp_perturbed: impl Fn(&mut crate::SparseMatrix, &mut [f64]),
    delta_p: f64,
) -> f64 {
    // Build original stamps
    let mut sys_orig = LinearSystem::new(dim);
    stamp_original(&mut sys_orig.matrix, &mut sys_orig.rhs);

    // Build perturbed stamps
    let mut sys_pert = LinearSystem::new(dim);
    stamp_perturbed(&mut sys_pert.matrix, &mut sys_pert.rhs);

    // Compute delta_Y * x_op for both
    let y_orig_x = mat_vec_product(&sys_orig.matrix, solution);
    let y_pert_x = mat_vec_product(&sys_pert.matrix, solution);

    // delta_residual = (b_pert - Y_pert*x) - (b_orig - Y_orig*x)
    //                = (b_pert - b_orig) - (Y_pert - Y_orig)*x
    let mut delta_residual = vec![0.0; dim];
    for i in 0..dim {
        let delta_b = sys_pert.rhs[i] - sys_orig.rhs[i];
        let delta_yx = y_pert_x[i] - y_orig_x[i];
        delta_residual[i] = delta_b - delta_yx;
    }

    // S = lambda^T * delta_residual / delta_p
    let dot: f64 = lambda
        .iter()
        .zip(delta_residual.iter())
        .map(|(l, d)| l * d)
        .sum();
    dot / delta_p
}

/// Compute matrix-vector product Y * x using sparse triplets.
fn mat_vec_product(matrix: &crate::SparseMatrix, x: &[f64]) -> Vec<f64> {
    let dim = matrix.dim();
    let mut result = vec![0.0; dim];
    for t in matrix.triplets() {
        result[t.row] += t.value * x[t.col];
    }
    result
}

/// Parse an output specification for sensitivity.
fn parse_sens_output(
    output: &str,
    mna: &MnaSystem,
) -> Result<(Option<usize>, Option<usize>, bool), MnaError> {
    let lower = output.to_lowercase();
    if lower.starts_with("v(") && lower.ends_with(')') {
        let inner = &lower[2..lower.len() - 1];
        let parts: Vec<&str> = inner.split(',').collect();
        let pos = mna.node_map.get(parts[0]);
        let neg = if parts.len() > 1 {
            mna.node_map.get(parts[1])
        } else {
            None
        };
        Ok((pos, neg, true))
    } else if lower.starts_with("i(") && lower.ends_with(')') {
        let src_name = &lower[2..lower.len() - 1];
        let num_nodes = mna.total_num_nodes();
        if let Some(pos) = mna
            .vsource_names
            .iter()
            .position(|n| n.eq_ignore_ascii_case(src_name))
        {
            Ok((Some(num_nodes + pos), None, false))
        } else {
            Err(MnaError::UnsupportedElement(format!(
                "source '{src_name}' not found"
            )))
        }
    } else {
        Err(MnaError::UnsupportedElement(format!(
            "unsupported sens output: {output}"
        )))
    }
}

/// Compute DC sensitivity analysis (.sens).
///
/// Computes d(output)/d(parameter) for each element parameter using the adjoint method:
/// 1. Solve Y*x = b (DC operating point)
/// 2. Build Jacobian Y
/// 3. Solve Y^T * lambda = e_out (adjoint)
/// 4. For each parameter p: S(p) = lambda^T * (db/dp - dY/dp * x)
pub fn simulate_sens(netlist: &Netlist) -> Result<SimResult, MnaError> {
    // Parse .sens analysis
    let sens_output = netlist
        .items
        .iter()
        .find_map(|item| {
            if let Item::Analysis(Analysis::Sens { output }) = item {
                Some(output.clone())
            } else {
                None
            }
        })
        .ok_or_else(|| MnaError::UnsupportedElement("no .sens analysis found".to_string()))?;

    let output_var = &sens_output[0];

    // Check DC vs AC (only DC supported currently)
    let is_ac = sens_output.len() > 1 && sens_output[1].eq_ignore_ascii_case("ac");
    if is_ac {
        return Err(MnaError::UnsupportedElement(
            "AC sensitivity not yet supported".to_string(),
        ));
    }

    // Assemble MNA and solve DC OP
    let mna = assemble_mna(netlist)?;
    let solution = solve_op_raw(&mna)?;
    let dim = mna.system.dim();
    let num_nodes = mna.total_num_nodes();

    // Build Jacobian
    let jacobian = build_jacobian(&mna, &solution);

    // Parse output spec and build selection vector
    let (out_pos, out_neg, _out_is_voltage) = parse_sens_output(output_var, &mna)?;
    let mut e_out = vec![0.0; dim];
    if let Some(i) = out_pos {
        e_out[i] = 1.0;
    }
    if let Some(i) = out_neg {
        e_out[i] = -1.0;
    }

    // Solve adjoint: Y^T * lambda = e_out
    let lambda = solve_transposed(&jacobian, &e_out)?;

    // Enumerate parameters and compute sensitivities
    let mut sens_names: Vec<String> = Vec::new();
    let mut sens_values: Vec<f64> = Vec::new();

    // === Resistors ===
    for r in &mna.resistors {
        let g = 1.0 / r.resistance;
        let name = r.name.to_lowercase();

        // Sensitivity to resistance R: dG/dR = -1/R^2
        let dg_dr = -g * g;
        sens_names.push(name.clone());
        sens_values.push(conductance_sensitivity(
            &lambda, r.pos_idx, r.neg_idx, dg_dr, &solution,
        ));

        // Sensitivity to multiplier m (default 1): dG/dm = G/m
        let m = 1.0;
        let dg_dm = g / m;
        sens_names.push(format!("{name}_m"));
        sens_values.push(conductance_sensitivity(
            &lambda, r.pos_idx, r.neg_idx, dg_dm, &solution,
        ));

        // Sensitivity to scale (default 1): dG/dscale = -G/scale
        let scale = 1.0;
        let dg_dscale = -g / scale;
        sens_names.push(format!("{name}_scale"));
        sens_values.push(conductance_sensitivity(
            &lambda, r.pos_idx, r.neg_idx, dg_dscale, &solution,
        ));
    }

    // === Voltage sources ===
    for (i, vsrc) in mna.vsource_names.iter().enumerate() {
        let branch_idx = num_nodes + i;
        // dS/dV: db/dV = 1 at branch equation row
        // S = lambda^T * e_branch = lambda[branch_idx]
        sens_names.push(vsrc.to_lowercase());
        sens_values.push(lambda[branch_idx]);
    }

    // === Current sources ===
    for cs in &mna.current_sources {
        let name = cs.name.to_lowercase();

        // dS/dI: current source stamp is rhs[pos] -= I, rhs[neg] += I
        // db/dI: rhs[pos] -= 1, rhs[neg] += 1
        // S = lambda^T * db/dI = -lambda[pos] + lambda[neg]
        let lp = cs.pos_idx.map(|i| lambda[i]).unwrap_or(0.0);
        let ln = cs.neg_idx.map(|i| lambda[i]).unwrap_or(0.0);
        sens_names.push(name.clone());
        sens_values.push(-lp + ln);

        // Sensitivity to multiplier m (default 1): I_eff = I_dc * m
        // dI_eff/dm = I_dc → same stamp pattern scaled by I_dc
        let m = 1.0;
        let di_dm = cs.dc_value / m;
        sens_names.push(format!("{name}_m"));
        sens_values.push(di_dm * (-lp + ln));
    }

    // === BJT model/instance parameters (numerical) ===
    for bjt in &mna.bjts {
        let name = bjt.name.to_lowercase();
        let (vbe, vbc) = bjt.junction_voltages(&solution);

        // Model parameters to compute sensitivity for
        type BjtSetter = Box<dyn Fn(&mut crate::bjt::BjtModel, f64)>;
        let model_params: &[(&str, f64, BjtSetter)] = &[
            ("bf", bjt.model.bf, Box::new(|m, v| m.bf = v)),
            ("br", bjt.model.br, Box::new(|m, v| m.br = v)),
            ("is", bjt.model.is, Box::new(|m, v| m.is = v)),
            ("vaf", bjt.model.vaf, Box::new(|m, v| m.vaf = v)),
            ("nf", bjt.model.nf, Box::new(|m, v| m.nf = v)),
            ("nr", bjt.model.nr, Box::new(|m, v| m.nr = v)),
            ("ne", bjt.model.ne, Box::new(|m, v| m.ne = v)),
            ("nc", bjt.model.nc, Box::new(|m, v| m.nc = v)),
            ("rb", bjt.model.rb, Box::new(|m, v| m.rb = v)),
            ("rc", bjt.model.rc, Box::new(|m, v| m.rc = v)),
            ("re", bjt.model.re, Box::new(|m, v| m.re = v)),
        ];

        for (param_name, param_value, setter) in model_params {
            let delta = if param_value.abs() > 1e-20 {
                param_value * 1e-6
            } else {
                1e-10
            };

            let s = numerical_device_sensitivity(
                &lambda,
                &solution,
                &mna,
                dim,
                |matrix, rhs| {
                    let comp = bjt.model.companion(vbe, vbc);
                    stamp_bjt(matrix, rhs, bjt, &comp);
                },
                |matrix, rhs| {
                    let mut model = bjt.model.clone();
                    setter(&mut model, *param_value + delta);
                    let comp = model.companion(vbe, vbc);
                    stamp_bjt(matrix, rhs, bjt, &comp);
                },
                delta,
            );

            sens_names.push(format!("{name}:{param_name}"));
            sens_values.push(s);
        }

        // Instance parameters: area, m
        let instance_params: &[(&str, f64)] = &[("area", bjt.area), ("m", bjt.m)];

        for (param_name, param_value) in instance_params {
            let delta = if param_value.abs() > 1e-20 {
                param_value * 1e-6
            } else {
                1e-10
            };

            let s = numerical_device_sensitivity(
                &lambda,
                &solution,
                &mna,
                dim,
                |matrix, rhs| {
                    let comp = bjt.model.companion(vbe, vbc);
                    stamp_bjt(matrix, rhs, bjt, &comp);
                },
                |matrix, rhs| {
                    let mut inst = bjt.clone();
                    match *param_name {
                        "area" => inst.area = *param_value + delta,
                        "m" => inst.m = *param_value + delta,
                        _ => {}
                    }
                    let comp = inst.model.companion(vbe, vbc);
                    stamp_bjt(matrix, rhs, &inst, &comp);
                },
                delta,
            );

            sens_names.push(format!("{name}_{param_name}"));
            sens_values.push(s);
        }
    }

    // === MOSFET parameters (numerical) ===
    for mos in &mna.mosfets {
        let name = mos.name.to_lowercase();
        let (vgs, vds, vbs) = mos.terminal_voltages(&solution);

        type MosSetter = Box<dyn Fn(&mut crate::mosfet::MosfetModel, f64)>;
        let mos_params: &[(&str, f64, MosSetter)] = &[
            ("vto", mos.model.vto, Box::new(|m, v| m.vto = v)),
            ("kp", mos.model.kp, Box::new(|m, v| m.kp = v)),
            ("lambda", mos.model.lambda, Box::new(|m, v| m.lambda = v)),
        ];

        for (param_name, param_value, setter) in mos_params {
            let delta = if param_value.abs() > 1e-20 {
                param_value * 1e-6
            } else {
                1e-10
            };

            let s = numerical_device_sensitivity(
                &lambda,
                &solution,
                &mna,
                dim,
                |matrix, rhs| {
                    let mut eff = mos.model.clone();
                    eff.kp = mos.beta();
                    let comp = eff.companion(vgs, vds, vbs);
                    stamp_mosfet(matrix, rhs, mos, &comp);
                },
                |matrix, rhs| {
                    let mut model = mos.model.clone();
                    setter(&mut model, *param_value + delta);
                    let l_eff = mos.l - 2.0 * model.ld;
                    if l_eff > 0.0 {
                        model.kp = model.kp * mos.w / l_eff;
                    }
                    let comp = model.companion(vgs, vds, vbs);
                    stamp_mosfet(matrix, rhs, mos, &comp);
                },
                delta,
            );

            sens_names.push(format!("{name}:{param_name}"));
            sens_values.push(s);
        }
    }

    // === JFET parameters (numerical) ===
    for jfet in &mna.jfets {
        let name = jfet.name.to_lowercase();
        let (vgs, vgd) = jfet.junction_voltages(&solution);

        type JfetSetter = Box<dyn Fn(&mut crate::jfet::JfetModel, f64)>;
        let jfet_params: &[(&str, f64, JfetSetter)] = &[
            ("vto", jfet.model.vto, Box::new(|m, v| m.vto = v)),
            ("beta", jfet.model.beta, Box::new(|m, v| m.beta = v)),
            ("lambda", jfet.model.lambda, Box::new(|m, v| m.lambda = v)),
        ];

        for (param_name, param_value, setter) in jfet_params {
            let delta = if param_value.abs() > 1e-20 {
                param_value * 1e-6
            } else {
                1e-10
            };

            let s = numerical_device_sensitivity(
                &lambda,
                &solution,
                &mna,
                dim,
                |matrix, rhs| {
                    let comp = jfet.model.companion(vgs, vgd);
                    stamp_jfet(matrix, rhs, jfet, &comp);
                },
                |matrix, rhs| {
                    let mut model = jfet.model.clone();
                    setter(&mut model, *param_value + delta);
                    let comp = model.companion(vgs, vgd);
                    stamp_jfet(matrix, rhs, jfet, &comp);
                },
                delta,
            );

            sens_names.push(format!("{name}:{param_name}"));
            sens_values.push(s);
        }
    }

    // Build result: one data row with all sensitivities
    let mut vecs: Vec<SimVector> = sens_names
        .into_iter()
        .zip(sens_values)
        .map(|(name, value)| SimVector {
            name,
            real: vec![value],
            complex: vec![],
        })
        .collect();

    // Add index vector at the beginning
    vecs.insert(
        0,
        SimVector {
            name: "Index".to_string(),
            real: vec![0.0],
            complex: vec![],
        },
    );

    Ok(SimResult {
        plots: vec![SimPlot {
            name: "sens1".to_string(),
            vecs,
        }],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    fn sens_value(result: &SimResult, name: &str) -> f64 {
        result.plots[0]
            .vecs
            .iter()
            .find(|v| v.name == name)
            .unwrap_or_else(|| {
                let names: Vec<_> = result.plots[0].vecs.iter().map(|v| &v.name).collect();
                panic!("no vector '{name}', available: {names:?}")
            })
            .real[0]
    }

    #[test]
    fn test_sens_dc_simple_ir() {
        // sens-dc-1: I1=42mA, R1=1k → V(1) = 42V
        // S[R1] = dV/dR = I = 42mA
        // S[I1] = dV/dI = R = 1k
        // S[R1_m] = -V/m = -42
        // S[I1_m] = I_dc * R / m = 42
        // S[R1_scale] = V/scale = 42
        let netlist = Netlist::parse(
            "sens dc test
i1 0 1 DC 42m
r1 1 0 1k
.sens v(1) dc
.end
",
        )
        .unwrap();

        let result = simulate_sens(&netlist).unwrap();

        assert_abs_diff_eq!(sens_value(&result, "i1"), 1000.0, epsilon = 1e-6);
        assert_abs_diff_eq!(sens_value(&result, "i1_m"), 42.0, epsilon = 1e-6);
        assert_abs_diff_eq!(sens_value(&result, "r1"), 0.042, epsilon = 1e-9);
        assert_abs_diff_eq!(sens_value(&result, "r1_m"), -42.0, epsilon = 1e-6);
        assert_abs_diff_eq!(sens_value(&result, "r1_scale"), 42.0, epsilon = 1e-6);
    }

    #[test]
    fn test_sens_dc_resistor_network() {
        // sens-dc-2: V1=42V, R1=1k, R2=1.5k, R3=2.2k, R4=3.3k, R5=1.8k, Rx=2.7k
        // Sensitivity of V(4) to each component
        let netlist = Netlist::parse(
            "sens dc resistor network
v1 1 0 DC 42
r1 1 2 1k
r2 2 3 1.5k
r3 3 0 2.2k
r4 3 4 3.3k
r5 4 0 1.8k
rx 2 4 2.7k
.sens v(4) dc
.end
",
        )
        .unwrap();

        let result = simulate_sens(&netlist).unwrap();

        // Golden values from ngspice
        assert_abs_diff_eq!(
            sens_value(&result, "v1"),
            0.2933065931311105,
            epsilon = 1e-6
        );
        assert_abs_diff_eq!(
            sens_value(&result, "r1"),
            -0.004104514242490192,
            epsilon = 1e-9
        );
        assert_abs_diff_eq!(
            sens_value(&result, "r2"),
            -2.954317312283796e-4,
            epsilon = 1e-9
        );
        assert_abs_diff_eq!(
            sens_value(&result, "r3"),
            0.001129248920024266,
            epsilon = 1e-9
        );
        assert_abs_diff_eq!(
            sens_value(&result, "r4"),
            -2.00582596465584e-4,
            epsilon = 1e-9
        );
        assert_abs_diff_eq!(
            sens_value(&result, "r5"),
            0.003755608696037442,
            epsilon = 1e-9
        );
        assert_abs_diff_eq!(
            sens_value(&result, "rx"),
            -0.001494392173796886,
            epsilon = 1e-9
        );
    }

    #[test]
    fn test_sens_dc_voltage_divider() {
        // V1=10V, R1=1k, R2=1k → V(mid) = 5V
        // S[V1] = R2/(R1+R2) = 0.5
        // S[R1] = dV/dR1 = -V1*R2/(R1+R2)^2 = -10*1000/4e6 = -2.5e-3
        // S[R2] = dV/dR2 = V1*R1/(R1+R2)^2 = 10*1000/4e6 = 2.5e-3
        let netlist = Netlist::parse(
            "sens voltage divider
v1 1 0 10
r1 1 mid 1k
r2 mid 0 1k
.sens v(mid) dc
.end
",
        )
        .unwrap();

        let result = simulate_sens(&netlist).unwrap();

        assert_abs_diff_eq!(sens_value(&result, "v1"), 0.5, epsilon = 1e-9);
        assert_abs_diff_eq!(sens_value(&result, "r1"), -2.5e-3, epsilon = 1e-9);
        assert_abs_diff_eq!(sens_value(&result, "r2"), 2.5e-3, epsilon = 1e-9);
    }
}
