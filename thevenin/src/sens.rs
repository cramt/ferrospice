use faer::Mat;
use faer::linalg::solvers::FullPivLu;
use faer::prelude::Solve;

use thevenin_types::{Analysis, Item, Netlist, SimPlot, SimResult, SimVector};

use crate::LinearSystem;
use crate::bjt::stamp_bjt;
use crate::jfet::stamp_jfet;
use crate::mna::{MnaError, MnaSystem, assemble_mna};
use crate::mosfet::stamp_mosfet;
use crate::newton::NrOptions;
use crate::simulate::solve_op_raw;
use crate::tf::build_jacobian;

/// Solve Y * delta_E = z and extract the output sensitivity.
///
/// The direct sensitivity method (as used by ngspice): for each parameter p, build
/// the perturbation z = delta_b - delta_Y*x, then solve Y * delta_E = z and read
/// off the output change delta_E[out_pos] - delta_E[out_neg].
///
/// This is numerically more stable than the adjoint method for circuits where
/// the output is a small difference of large node voltages, because the cancellation
/// happens naturally inside the sparse solve rather than in a dot product with a
/// possibly-inaccurate adjoint vector.
fn extract_output(delta_e: &Mat<f64>, out_pos: Option<usize>, out_neg: Option<usize>) -> f64 {
    let v_pos = out_pos.map(|i| delta_e[(i, 0)]).unwrap_or(0.0);
    let v_neg = out_neg.map(|i| delta_e[(i, 0)]).unwrap_or(0.0);
    v_pos - v_neg
}

/// Solve Y * x = rhs using a pre-factored LU.
fn solve_forward(lu: &FullPivLu<f64>, rhs: &[f64]) -> Mat<f64> {
    let dim = rhs.len();
    let mut b = Mat::zeros(dim, 1);
    for (i, &val) in rhs.iter().enumerate() {
        b[(i, 0)] = val;
    }
    lu.solve(&b)
}

/// Compute direct sensitivity for a conductance parameter change.
///
/// When conductance G between nodes (pos, neg) changes by dG:
///   z[pos] = -dG * V_across,  z[neg] = +dG * V_across
/// Solve Y * delta_E = z, return (delta_E[out_pos] - delta_E[out_neg]) / (dG/G * param)
fn conductance_sensitivity_direct(
    lu: &FullPivLu<f64>,
    pos: Option<usize>,
    neg: Option<usize>,
    dg: f64,
    solution: &[f64],
    out_pos: Option<usize>,
    out_neg: Option<usize>,
) -> f64 {
    let v_pos = pos.map(|i| solution[i]).unwrap_or(0.0);
    let v_neg = neg.map(|i| solution[i]).unwrap_or(0.0);
    let v_across = v_pos - v_neg;

    let dim = solution.len();
    let mut z = vec![0.0f64; dim];
    // z = -(dY * x): conductance change dG between pos-neg
    // dY * x: at pos gives +dG*(V_pos - V_neg), at neg gives -dG*(V_pos - V_neg)
    if let Some(p) = pos {
        z[p] -= dg * v_across;
    }
    if let Some(n) = neg {
        z[n] += dg * v_across;
    }

    let delta_e = solve_forward(lu, &z);
    extract_output(&delta_e, out_pos, out_neg)
}

/// Compute sensitivity using direct method: solve Y * delta_E = z.
///
/// For nonlinear devices, perturb the parameter and recompute stamps.
/// Uses: z = (b_pert - Y_pert*x) - (b_orig - Y_orig*x) = delta_b - delta_Y*x
/// Then solves Y * delta_E = z and returns (delta_E[out_pos] - delta_E[out_neg]) / delta_p.
#[allow(clippy::too_many_arguments)]
fn numerical_device_sensitivity_direct(
    lu: &FullPivLu<f64>,
    solution: &[f64],
    dim: usize,
    stamp_original: impl Fn(&mut crate::SparseMatrix, &mut [f64]),
    stamp_perturbed: impl Fn(&mut crate::SparseMatrix, &mut [f64]),
    delta_p: f64,
    out_pos: Option<usize>,
    out_neg: Option<usize>,
) -> f64 {
    // Build original stamps
    let mut sys_orig = LinearSystem::new(dim);
    stamp_original(&mut sys_orig.matrix, &mut sys_orig.rhs);

    // Build perturbed stamps
    let mut sys_pert = LinearSystem::new(dim);
    stamp_perturbed(&mut sys_pert.matrix, &mut sys_pert.rhs);

    // z = (b_pert - Y_pert*x) - (b_orig - Y_orig*x)
    let y_orig_x = mat_vec_product(&sys_orig.matrix, solution);
    let y_pert_x = mat_vec_product(&sys_pert.matrix, solution);
    let mut z = vec![0.0f64; dim];
    for i in 0..dim {
        let delta_b = sys_pert.rhs[i] - sys_orig.rhs[i];
        let delta_yx = y_pert_x[i] - y_orig_x[i];
        z[i] = delta_b - delta_yx;
    }

    let delta_e = solve_forward(lu, &z);
    extract_output(&delta_e, out_pos, out_neg) / delta_p
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
/// Computes d(output)/d(parameter) for each element parameter using the direct method
/// (matching ngspice's cktsens.c approach):
/// 1. Solve Y*x = b (DC operating point)
/// 2. Build Jacobian Y and factor it (LU decomposition)
/// 3. For each parameter p:
///    a. Compute z = delta_b - delta_Y*x (perturbation of companion stamps)
///    b. Solve Y * delta_E = z (forward solve with same LU)
///    c. S(p) = (delta_E[out_pos] - delta_E[out_neg]) / delta_p
///
/// The direct method is numerically more stable than the adjoint method for
/// circuits where the sensitivity is a small difference of large contributions,
/// because the near-cancellation happens inside the sparse linear solve rather
/// than in a dot product with a potentially-inaccurate adjoint vector.
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

    // Assemble MNA and solve DC OP.
    // Use diag_gmin = 0 to match ngspice's CKTdiagGmin = 0 for DC analysis.
    let mna = assemble_mna(netlist)?;
    let solution = if !mna.has_nonlinear() {
        solve_op_raw(&mna)?
    } else {
        let opts = NrOptions { diag_gmin: 0.0, ..NrOptions::default() };
        crate::simulate::solve_nonlinear_op(&mna, &opts)?
    };
    let dim = mna.system.dim();
    let num_nodes = mna.total_num_nodes();

    // Build the Jacobian and factor it using dense LU.
    let jacobian = build_jacobian(&mna, &solution);
    let lu = FullPivLu::new(jacobian.as_ref());

    // Parse output spec
    let (out_pos, out_neg, _out_is_voltage) = parse_sens_output(output_var, &mna)?;

    // Enumerate parameters and compute sensitivities.
    // ngspice outputs device types in alphabetical order by first letter:
    //   B (J)FET → M (MOSFET) → Q (BJT) → R (Resistor) → V (Vsrc) → I (Isrc)
    // Wait — actual ngspice order is Q→R→V for diffpair (no JFETs/MOSFETs/Isrcs there).
    // General rule: alphabetical by device-type letter, then alphabetical within each type.
    // Device type letters: J=JFET, M=MOSFET, Q=BJT, R=Resistor, V=Vsource, I=Isource
    // Alphabetical order: I < J < M < Q < R < V  → but ngspice puts Q before R before V.
    // Empirically from diffpair.out: Q then R then V. So order is Q→R→V (no I/J/M in that test).
    // For generality we use: Q → M → J → R → V → I (each group sorted by name).

    let mut sens_names: Vec<String> = Vec::new();
    let mut sens_values: Vec<f64> = Vec::new();

    // Helper to collect (name, sensitivity) pairs for a group, then sort and append.
    macro_rules! append_group {
        ($pairs:expr) => {{
            let mut pairs: Vec<(String, f64)> = $pairs;
            pairs.sort_by(|a, b| a.0.cmp(&b.0));
            for (n, v) in pairs {
                sens_names.push(n);
                sens_values.push(v);
            }
        }};
    }

    // === BJT model/instance parameters (numerical) ===
    // Circuit temperature in Kelvin (used for temperature-dependent sensitivities).
    let ckt_temp_k = crate::netlist_temp(netlist) + 273.15;

    {
        let mut bjt_pairs: Vec<(String, f64)> = Vec::new();
        for bjt in &mna.bjts {
            let name = bjt.name.to_lowercase();
            let (vbe, vbc) = bjt.junction_voltages(&solution);

            // Model parameters in ngspice alphabetical order.
            // Note: rc and re are NOT included (ngspice doesn't report them in .sens).
            // For tnom: use companion_at(ckt_temp_k) to capture IS temperature dependence.
            // For eg, fc, xti: sensitivity is theoretically 0 at T=TNOM; companion() used.
            type BjtSetter = Box<dyn Fn(&mut crate::bjt::BjtModel, f64)>;
            let model_params: &[(&str, f64, BjtSetter)] = &[
                ("bf",   bjt.model.bf,   Box::new(|m, v| m.bf = v)),
                ("br",   bjt.model.br,   Box::new(|m, v| m.br = v)),
                ("eg",   bjt.model.eg,   Box::new(|m, v| m.eg = v)),
                ("fc",   bjt.model.fc,   Box::new(|m, v| m.fc = v)),
                ("is",   bjt.model.is,   Box::new(|m, v| m.is = v)),
                ("nc",   bjt.model.nc,   Box::new(|m, v| m.nc = v)),
                ("ne",   bjt.model.ne,   Box::new(|m, v| m.ne = v)),
                ("nf",   bjt.model.nf,   Box::new(|m, v| m.nf = v)),
                ("nr",   bjt.model.nr,   Box::new(|m, v| m.nr = v)),
                ("rb",   bjt.model.rb,   Box::new(|m, v| m.rb = v)),
                ("rbm",  bjt.model.rbm,  Box::new(|m, v| m.rbm = v)),
                ("tnom", bjt.model.tnom, Box::new(|m, v| m.tnom = v)),
                ("vaf",  bjt.model.vaf,  Box::new(|m, v| m.vaf = v)),
                ("xti",  bjt.model.xti,  Box::new(|m, v| m.xti = v)),
            ];

            for (param_name, param_value, setter) in model_params {
                let delta = if param_value.abs() > 1e-20 {
                    param_value * 1e-6
                } else {
                    1e-10
                };

                // For tnom sensitivity, use companion_at to capture IS temperature dependence.
                // For other params, use companion() (which uses VT_NOM / self.is directly).
                let use_temp_companion = *param_name == "tnom";
                let s = numerical_device_sensitivity_direct(
                    &lu,
                    &solution,
                    dim,
                    |matrix, rhs| {
                        let comp = if use_temp_companion {
                            bjt.model.companion_at(vbe, vbc, ckt_temp_k)
                        } else {
                            bjt.model.companion(vbe, vbc)
                        };
                        stamp_bjt(matrix, rhs, bjt, &comp);
                    },
                    |matrix, rhs| {
                        let mut inst = bjt.clone();
                        setter(&mut inst.model, *param_value + delta);
                        let comp = if use_temp_companion {
                            inst.model.companion_at(vbe, vbc, ckt_temp_k)
                        } else {
                            inst.model.companion(vbe, vbc)
                        };
                        stamp_bjt(matrix, rhs, &inst, &comp);
                    },
                    delta,
                    out_pos,
                    out_neg,
                );

                bjt_pairs.push((format!("{name}:{param_name}"), s));
            }

            // Instance parameters in ngspice order: area, areab, areac, m, temp
            // areab and areac only affect junction capacitances (no DC effect → ~0 sensitivity).
            // temp uses companion_at to capture IS/VT temperature dependence.
            let instance_params: &[(&str, f64)] = &[
                ("area",  bjt.area),
                ("areab", bjt.areab),
                ("areac", bjt.areac),
                ("m",     bjt.m),
                ("temp",  bjt.temp),
            ];

            for (param_name, param_value) in instance_params {
                // For temp, use circuit temperature as base if not explicitly set.
                let base_temp = if *param_name == "temp" && param_value.is_nan() {
                    ckt_temp_k - 273.15 // °C
                } else {
                    *param_value
                };

                let delta = if base_temp.abs() > 1e-20 {
                    base_temp * 1e-6
                } else {
                    1e-6
                };

                let s = numerical_device_sensitivity_direct(
                    &lu,
                    &solution,
                    dim,
                    |matrix, rhs| {
                        let comp = if *param_name == "temp" {
                            bjt.model.companion_at(vbe, vbc, ckt_temp_k)
                        } else {
                            bjt.model.companion(vbe, vbc)
                        };
                        stamp_bjt(matrix, rhs, bjt, &comp);
                    },
                    |matrix, rhs| {
                        let mut inst = bjt.clone();
                        match *param_name {
                            "area"  => inst.area  = base_temp + delta,
                            "areab" => inst.areab = base_temp + delta,
                            "areac" => inst.areac = base_temp + delta,
                            "m"     => inst.m     = base_temp + delta,
                            "temp"  => {
                                // Perturb device temperature; recompute IS_t
                                let new_temp_k = (base_temp + delta) + 273.15;
                                let comp = inst.model.companion_at(vbe, vbc, new_temp_k);
                                stamp_bjt(matrix, rhs, &inst, &comp);
                                return;
                            }
                            _ => {}
                        }
                        let comp = inst.model.companion(vbe, vbc);
                        stamp_bjt(matrix, rhs, &inst, &comp);
                    },
                    delta,
                    out_pos,
                    out_neg,
                );

                bjt_pairs.push((format!("{name}_{param_name}"), s));
            }
        }
        append_group!(bjt_pairs);
    }

    // === MOSFET parameters (numerical) ===
    {
        let mut mos_pairs: Vec<(String, f64)> = Vec::new();
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

                let s = numerical_device_sensitivity_direct(
                    &lu,
                    &solution,
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
                    out_pos,
                    out_neg,
                );

                mos_pairs.push((format!("{name}:{param_name}"), s));
            }
        }
        append_group!(mos_pairs);
    }

    // === JFET parameters (numerical) ===
    {
        let mut jfet_pairs: Vec<(String, f64)> = Vec::new();
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

                let s = numerical_device_sensitivity_direct(
                    &lu,
                    &solution,
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
                    out_pos,
                    out_neg,
                );

                jfet_pairs.push((format!("{name}:{param_name}"), s));
            }
        }
        append_group!(jfet_pairs);
    }

    // === Resistors ===
    {
        let mut res_pairs: Vec<(String, f64)> = Vec::new();
        for r in &mna.resistors {
            let g = 1.0 / r.resistance;
            let name = r.name.to_lowercase();

            // Sensitivity to resistance R: dG/dR = -1/R^2
            let dg_dr = -g * g;
            res_pairs.push((name.clone(), conductance_sensitivity_direct(
                &lu, r.pos_idx, r.neg_idx, dg_dr, &solution, out_pos, out_neg,
            )));

            // Sensitivity to multiplier m (default 1): dG/dm = G/m
            let dg_dm = g; // m=1
            res_pairs.push((format!("{name}_m"), conductance_sensitivity_direct(
                &lu, r.pos_idx, r.neg_idx, dg_dm, &solution, out_pos, out_neg,
            )));

            // Sensitivity to scale (default 1): dG/dscale = -G/scale
            let dg_dscale = -g; // scale=1
            res_pairs.push((format!("{name}_scale"), conductance_sensitivity_direct(
                &lu, r.pos_idx, r.neg_idx, dg_dscale, &solution, out_pos, out_neg,
            )));
        }
        append_group!(res_pairs);
    }

    // === Voltage sources ===
    // dV(out)/dV_k: perturbing voltage source k by 1V → RHS row branch_idx changes by 1
    // z[branch_idx] = 1, all others 0 → solve Y * delta_E = z
    {
        let mut vsrc_pairs: Vec<(String, f64)> = Vec::new();
        for (i, vsrc) in mna.vsource_names.iter().enumerate() {
            let branch_idx = num_nodes + i;
            let mut z = vec![0.0f64; dim];
            z[branch_idx] = 1.0;
            let delta_e = solve_forward(&lu, &z);
            let s = extract_output(&delta_e, out_pos, out_neg);
            vsrc_pairs.push((vsrc.to_lowercase(), s));
        }
        append_group!(vsrc_pairs);
    }

    // === Current sources ===
    // dV(out)/dI_k: perturbing current source k by 1A → RHS[pos] += 1, RHS[neg] -= 1
    {
        let mut isrc_pairs: Vec<(String, f64)> = Vec::new();
        for cs in &mna.current_sources {
            let name = cs.name.to_lowercase();
            let mut z = vec![0.0f64; dim];
            // SPICE convention: current exits n+ (rhs[pos] -= I) and enters n-  (rhs[neg] += I)
            // So when I increases by 1: delta_b[pos] = -1, delta_b[neg] = +1
            if let Some(p) = cs.pos_idx { z[p] -= 1.0; }
            if let Some(n) = cs.neg_idx { z[n] += 1.0; }
            let delta_e = solve_forward(&lu, &z);
            let s = extract_output(&delta_e, out_pos, out_neg);
            isrc_pairs.push((name.clone(), s));
            let di_dm = cs.dc_value; // m=1
            isrc_pairs.push((format!("{name}_m"), s * di_dm));
        }
        append_group!(isrc_pairs);
    }

    // Build result: one data row with all sensitivities
    let vecs: Vec<SimVector> = sens_names
        .into_iter()
        .zip(sens_values)
        .map(|(name, value)| SimVector {
            name,
            real: vec![value],
            complex: vec![],
        })
        .collect();

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
