use faer::Mat;
use faer::linalg::solvers::FullPivLu;
use faer::prelude::Solve;

use ferrospice_netlist::{Analysis, ElementKind, Item, Netlist, SimPlot, SimResult, SimVector};

use crate::LinearSystem;
use crate::bjt::stamp_bjt;
use crate::jfet::stamp_jfet;
use crate::mna::{MnaError, MnaSystem, assemble_mna, stamp_conductance};
use crate::mosfet::stamp_mosfet;
use crate::simulate::{solve_op_raw, total_num_nodes};

/// Input source classification for .tf analysis.
enum InputSource {
    Voltage {
        branch_idx: usize,
    },
    Current {
        pos_idx: Option<usize>,
        neg_idx: Option<usize>,
    },
}

/// Parse an output specification like "v(5)", "v(5,4)", or "i(vsrc)".
/// Returns (positive_index, negative_index, is_voltage).
fn parse_output_spec(
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
        let num_nodes = total_num_nodes(mna);
        if let Some(pos) = mna
            .vsource_names
            .iter()
            .position(|n| n.eq_ignore_ascii_case(src_name))
        {
            Ok((Some(num_nodes + pos), None, false))
        } else {
            Err(MnaError::UnsupportedElement(format!(
                "source '{src_name}' not found for current output"
            )))
        }
    } else {
        Err(MnaError::UnsupportedElement(format!(
            "unsupported output spec: {output}"
        )))
    }
}

/// Find the input source type and indices.
fn find_input_source(
    input: &str,
    mna: &MnaSystem,
    netlist: &Netlist,
) -> Result<InputSource, MnaError> {
    let num_nodes = total_num_nodes(mna);

    // Check voltage sources
    if let Some(pos) = mna
        .vsource_names
        .iter()
        .position(|n| n.eq_ignore_ascii_case(input))
    {
        return Ok(InputSource::Voltage {
            branch_idx: num_nodes + pos,
        });
    }

    // Check current sources in the netlist
    for elem in netlist.elements() {
        if elem.name.eq_ignore_ascii_case(input)
            && let ElementKind::CurrentSource { pos, neg, .. } = &elem.kind
        {
            return Ok(InputSource::Current {
                pos_idx: mna.node_map.get(&pos.to_lowercase()),
                neg_idx: mna.node_map.get(&neg.to_lowercase()),
            });
        }
    }

    Err(MnaError::UnsupportedElement(format!(
        "input source '{input}' not found"
    )))
}

/// Build the linearized Jacobian matrix at the DC operating point.
/// For linear circuits, this is just the MNA matrix.
/// For nonlinear circuits, it includes linearized device conductances.
pub(crate) fn build_jacobian(mna: &MnaSystem, solution: &[f64]) -> Mat<f64> {
    let dim = mna.system.dim();
    let mut system = LinearSystem::new(dim);

    // Copy base linear stamps (R, V, I, L, controlled sources)
    for triplet in mna.system.matrix.triplets() {
        system.matrix.add(triplet.row, triplet.col, triplet.value);
    }

    // Add diode linearized conductances
    for diode in &mna.diodes {
        let (jct_a, jct_c) = if diode.internal_idx.is_some() {
            (diode.internal_idx, diode.cathode_idx)
        } else {
            (diode.anode_idx, diode.cathode_idx)
        };
        let v_a = jct_a.map(|i| solution[i]).unwrap_or(0.0);
        let v_c = jct_c.map(|i| solution[i]).unwrap_or(0.0);
        let v_jct = v_a - v_c;
        let (g_d, _) = diode.model.companion(v_jct);
        stamp_conductance(&mut system.matrix, jct_a, jct_c, g_d);

        if let Some(int_idx) = diode.internal_idx {
            let g_rs = 1.0 / diode.model.rs;
            stamp_conductance(&mut system.matrix, diode.anode_idx, Some(int_idx), g_rs);
        }
    }

    // Add BJT linearized conductances
    for bjt in &mna.bjts {
        let (vbe, vbc) = bjt.junction_voltages(solution);
        let comp = bjt.model.companion(vbe, vbc);
        stamp_bjt(&mut system.matrix, &mut system.rhs, bjt, &comp);
    }

    // Add MOSFET linearized conductances
    for mos in &mna.mosfets {
        let (vgs, vds, vbs) = mos.terminal_voltages(solution);
        let mut eff_model = mos.model.clone();
        eff_model.kp = mos.beta();
        let comp = eff_model.companion(vgs, vds, vbs);
        stamp_mosfet(&mut system.matrix, &mut system.rhs, mos, &comp);
    }

    // Add JFET linearized conductances
    for jfet in &mna.jfets {
        let (vgs, vgd) = jfet.junction_voltages(solution);
        let comp = jfet.model.companion(vgs, vgd);
        stamp_jfet(&mut system.matrix, &mut system.rhs, jfet, &comp);
    }

    system.matrix.to_dense()
}

/// Solve Y*x = rhs using LU factorization of the given dense matrix.
fn solve_dense(jacobian: &Mat<f64>, rhs: &[f64]) -> Vec<f64> {
    let dim = rhs.len();
    let lu = FullPivLu::new(jacobian.as_ref());
    let mut b = Mat::zeros(dim, 1);
    for (i, &val) in rhs.iter().enumerate() {
        b[(i, 0)] = val;
    }
    let x = lu.solve(&b);
    (0..dim).map(|i| x[(i, 0)]).collect()
}

/// Compute transfer function analysis (.tf).
///
/// For each `.tf output input` command, computes:
/// 1. Transfer function (output/input gain)
/// 2. Output impedance at the output port
/// 3. Input impedance at the input source
///
/// Algorithm follows ngspice tfanal.c:
/// - Compute DC operating point
/// - Build linearized Jacobian Y at the OP
/// - Inject unit excitation at input, solve Y*x = rhs, extract TF + input Z
/// - Inject unit excitation at output, solve Y*x = rhs, extract output Z
pub fn simulate_tf(netlist: &Netlist) -> Result<SimResult, MnaError> {
    // Collect all .tf analyses
    let tf_params: Vec<_> = netlist
        .items
        .iter()
        .filter_map(|item| {
            if let Item::Analysis(Analysis::Tf { output, input }) = item {
                Some((output.clone(), input.clone()))
            } else {
                None
            }
        })
        .collect();

    if tf_params.is_empty() {
        return Err(MnaError::UnsupportedElement(
            "no .tf analysis found".to_string(),
        ));
    }

    // Assemble MNA and solve DC OP
    let mna = assemble_mna(netlist)?;
    let solution = solve_op_raw(&mna)?;

    // Build linearized Jacobian at OP
    let jacobian = build_jacobian(&mna, &solution);

    let dim = mna.system.dim();
    let mut plots = Vec::new();

    for (output, input) in &tf_params {
        // Parse output spec
        let (out_pos, out_neg, out_is_voltage) = parse_output_spec(output, &mna)?;

        // Parse input spec
        let in_src = find_input_source(input, &mna, netlist)?;

        // === Transfer function + input impedance ===
        // Inject unit excitation at input
        let mut rhs1 = vec![0.0; dim];
        match &in_src {
            InputSource::Current { pos_idx, neg_idx } => {
                if let Some(i) = pos_idx {
                    rhs1[*i] -= 1.0;
                }
                if let Some(i) = neg_idx {
                    rhs1[*i] += 1.0;
                }
            }
            InputSource::Voltage { branch_idx } => {
                rhs1[*branch_idx] += 1.0;
            }
        }

        let x1 = solve_dense(&jacobian, &rhs1);

        // Extract transfer function
        let tf_value = if out_is_voltage {
            let v_pos = out_pos.map(|i| x1[i]).unwrap_or(0.0);
            let v_neg = out_neg.map(|i| x1[i]).unwrap_or(0.0);
            v_pos - v_neg
        } else {
            x1[out_pos.unwrap()]
        };

        // Extract input impedance
        let input_impedance = match &in_src {
            InputSource::Current { pos_idx, neg_idx } => {
                let v_pos = pos_idx.map(|i| x1[i]).unwrap_or(0.0);
                let v_neg = neg_idx.map(|i| x1[i]).unwrap_or(0.0);
                v_neg - v_pos
            }
            InputSource::Voltage { branch_idx } => {
                let branch_current = x1[*branch_idx];
                if branch_current.abs() < 1e-20 {
                    1e20
                } else {
                    -1.0 / branch_current
                }
            }
        };

        // === Output impedance ===
        let mut rhs2 = vec![0.0; dim];
        if out_is_voltage {
            if let Some(i) = out_pos {
                rhs2[i] -= 1.0;
            }
            if let Some(i) = out_neg {
                rhs2[i] += 1.0;
            }
        } else {
            rhs2[out_pos.unwrap()] += 1.0;
        }

        let x2 = solve_dense(&jacobian, &rhs2);

        let output_impedance = if out_is_voltage {
            let v_pos = out_pos.map(|i| x2[i]).unwrap_or(0.0);
            let v_neg = out_neg.map(|i| x2[i]).unwrap_or(0.0);
            v_neg - v_pos
        } else {
            let branch_val = x2[out_pos.unwrap()];
            if branch_val.abs() < 1e-20 {
                1e20
            } else {
                1.0 / branch_val
            }
        };

        // Build result vectors
        let vecs = vec![
            SimVector {
                name: "transfer_function".to_string(),
                real: vec![tf_value],
                complex: vec![],
            },
            SimVector {
                name: format!("output_impedance_at_{}", output.to_lowercase()),
                real: vec![output_impedance],
                complex: vec![],
            },
            SimVector {
                name: format!("{}#input_impedance", input.to_lowercase()),
                real: vec![input_impedance],
                complex: vec![],
            },
        ];

        plots.push(SimPlot {
            name: "tf1".to_string(),
            vecs,
        });
    }

    Ok(SimResult { plots })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    fn tf_value(result: &SimResult, plot_idx: usize, name: &str) -> f64 {
        result.plots[plot_idx]
            .vecs
            .iter()
            .find(|v| v.name == name)
            .unwrap_or_else(|| panic!("no vector '{name}'"))
            .real[0]
    }

    #[test]
    fn test_tf_voltage_divider() {
        // V1=10V, R1=1k, R2=1k → TF from V1 to V(mid) = 0.5
        // Input impedance = R1 + R2 = 2k (seen by V1)
        // Output impedance at V(mid) = R1||R2 = 500
        let netlist = Netlist::parse(
            "TF test - voltage divider
V1 1 0 10
R1 1 mid 1k
R2 mid 0 1k
.tf v(mid) V1
.end
",
        )
        .unwrap();

        let result = simulate_tf(&netlist).unwrap();
        assert_eq!(result.plots.len(), 1);

        let tf = tf_value(&result, 0, "transfer_function");
        assert_abs_diff_eq!(tf, 0.5, epsilon = 1e-9);

        let z_out = tf_value(&result, 0, "output_impedance_at_v(mid)");
        assert_abs_diff_eq!(z_out, 500.0, epsilon = 1e-6);

        let z_in = tf_value(&result, 0, "v1#input_impedance");
        assert_abs_diff_eq!(z_in, 2000.0, epsilon = 1e-6);
    }

    #[test]
    fn test_tf_current_source_input() {
        // I1=1mA, R1=2k → V(1) = 2V
        // TF = V(1)/I1 = R1 = 2k (transresistance)
        // Input impedance = R1 = 2k (looking into the current source)
        // Output impedance at V(1) = R1 = 2k
        let netlist = Netlist::parse(
            "TF test - current source
I1 0 1 1m
R1 1 0 2k
.tf v(1) I1
.end
",
        )
        .unwrap();

        let result = simulate_tf(&netlist).unwrap();

        let tf = tf_value(&result, 0, "transfer_function");
        assert_abs_diff_eq!(tf, 2000.0, epsilon = 1e-6);

        let z_in = tf_value(&result, 0, "i1#input_impedance");
        assert_abs_diff_eq!(z_in, 2000.0, epsilon = 1e-6);

        let z_out = tf_value(&result, 0, "output_impedance_at_v(1)");
        assert_abs_diff_eq!(z_out, 2000.0, epsilon = 1e-6);
    }

    #[test]
    fn test_tf_differential_output() {
        // V1=10V, R1=1k from 1 to a, R2=2k from 1 to b
        // V(a) = 10V (same node as V1+ since R1 from 1 to a)
        // Actually: V1 at node 1 = 10V, R1 from 1→a, R2 from a→0
        // V(a) = V1 * R2/(R1+R2) = 10 * 2/3 = 6.667V
        // TF of V(1,a) = V(1) - V(a) = 10 - 6.667 = 3.333V per 10V = 1/3
        let netlist = Netlist::parse(
            "TF test - differential output
V1 1 0 10
R1 1 a 1k
R2 a 0 2k
.tf v(1,a) V1
.end
",
        )
        .unwrap();

        let result = simulate_tf(&netlist).unwrap();

        let tf = tf_value(&result, 0, "transfer_function");
        // V(1) - V(a): V(1)=10V (forced by V1), V(a) = 10*2k/(1k+2k) = 20/3
        // V(1,a) = 10 - 20/3 = 10/3
        // TF = (10/3)/10 = 1/3 per unit voltage change at V1
        // Actually: injecting 1V at V1 → V(1) increases by 1, V(a) increases by 2/3
        // TF = 1 - 2/3 = 1/3
        assert_abs_diff_eq!(tf, 1.0 / 3.0, epsilon = 1e-9);
    }

    #[test]
    fn test_tf_multiple_analyses() {
        // Two .tf commands in same netlist
        let netlist = Netlist::parse(
            "Multiple TF
V1 1 0 10
R1 1 mid 1k
R2 mid 0 1k
.tf v(mid) V1
.tf v(1) V1
.end
",
        )
        .unwrap();

        let result = simulate_tf(&netlist).unwrap();
        assert_eq!(result.plots.len(), 2);

        let tf1 = tf_value(&result, 0, "transfer_function");
        assert_abs_diff_eq!(tf1, 0.5, epsilon = 1e-9);

        let tf2 = tf_value(&result, 1, "transfer_function");
        assert_abs_diff_eq!(tf2, 1.0, epsilon = 1e-9);
    }
}
