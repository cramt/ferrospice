use ferrospice_netlist::{Analysis, Expr, Item, Netlist, SimPlot, SimResult, SimVector};

use crate::mna::{MnaError, MnaSystem, assemble_mna};

/// Compute the DC operating point of a circuit.
///
/// Assembles the MNA system from the netlist, solves it, and returns
/// a `SimResult` containing node voltages and voltage source branch currents.
pub fn simulate_op(netlist: &Netlist) -> Result<SimResult, MnaError> {
    let mna = assemble_mna(netlist)?;
    let solution = mna.solve()?;

    let mut vecs = Vec::new();

    // Node voltages
    for (name, _idx) in mna.node_map.iter() {
        let v = solution.voltage(name).unwrap_or(0.0);
        vecs.push(SimVector {
            name: format!("v({})", name),
            real: vec![v],
            complex: vec![],
        });
    }

    // Voltage source branch currents
    for vsrc in &mna.vsource_names {
        let i = solution.branch_current(vsrc).unwrap_or(0.0);
        vecs.push(SimVector {
            name: format!("{}#branch", vsrc.to_lowercase()),
            real: vec![i],
            complex: vec![],
        });
    }

    Ok(SimResult {
        plots: vec![SimPlot {
            name: "op1".to_string(),
            vecs,
        }],
    })
}

/// Information about a source being swept in a DC analysis.
enum SweepSource {
    /// Voltage source: the branch equation RHS index.
    Voltage { rhs_index: usize },
    /// Current source: node indices (pos, neg) in the RHS vector.
    Current {
        ni: Option<usize>,
        nj: Option<usize>,
    },
}

/// Find the sweep source in the MNA system and return info for modifying the RHS.
fn find_sweep_source(mna: &MnaSystem, src_name: &str) -> Result<SweepSource, MnaError> {
    let src_lower = src_name.to_lowercase();

    // Check voltage sources (including inductors modeled as 0V sources)
    if let Some(pos) = mna
        .vsource_names
        .iter()
        .position(|n| n.to_lowercase() == src_lower)
    {
        return Ok(SweepSource::Voltage {
            rhs_index: mna.node_map.len() + pos,
        });
    }

    // Check current sources — need to find node connections from the netlist
    // This is handled by the caller providing the node info
    Err(MnaError::UnsupportedElement(format!(
        "sweep source '{src_name}' not found in circuit"
    )))
}

/// Find a current source's node indices in the MNA system.
fn find_current_source_sweep(
    mna: &MnaSystem,
    netlist: &Netlist,
    src_name: &str,
) -> Result<SweepSource, MnaError> {
    let src_lower = src_name.to_lowercase();
    for element in netlist.elements() {
        if element.name.to_lowercase() == src_lower
            && let ferrospice_netlist::ElementKind::CurrentSource { pos, neg, .. } = &element.kind
        {
            return Ok(SweepSource::Current {
                ni: mna.node_map.get(pos),
                nj: mna.node_map.get(neg),
            });
        }
    }
    Err(MnaError::UnsupportedElement(format!(
        "sweep source '{src_name}' not found in circuit"
    )))
}

/// Resolve a sweep source — tries voltage sources first, then current sources.
fn resolve_sweep_source(
    mna: &MnaSystem,
    netlist: &Netlist,
    src_name: &str,
) -> Result<SweepSource, MnaError> {
    match find_sweep_source(mna, src_name) {
        Ok(info) => Ok(info),
        Err(_) => find_current_source_sweep(mna, netlist, src_name),
    }
}

/// Generate sweep points from start to stop with given step.
fn generate_sweep_points(start: f64, stop: f64, step: f64) -> Vec<f64> {
    if step == 0.0 || (stop - start).signum() != step.signum() {
        return vec![start];
    }
    let mut points = Vec::new();
    let mut v = start;
    if step > 0.0 {
        while v <= stop + step * 1e-9 {
            points.push(v);
            v += step;
        }
    } else {
        while v >= stop + step * 1e-9 {
            points.push(v);
            v += step;
        }
    }
    points
}

/// Set the swept source value in the MNA system's RHS.
///
/// `original_rhs` is the RHS after initial assembly (before any sweep modifications).
/// For voltage sources, we simply set the branch equation RHS.
/// For current sources, we undo the original stamp and apply the new value.
fn set_source_value(
    mna: &mut MnaSystem,
    sweep: &SweepSource,
    new_value: f64,
    original_rhs: &[f64],
    original_dc: f64,
) {
    match sweep {
        SweepSource::Voltage { rhs_index } => {
            mna.system.rhs[*rhs_index] = new_value;
        }
        SweepSource::Current { ni, nj } => {
            // Undo original stamp and apply new value.
            // Original: rhs[ni] -= old_val, rhs[nj] += old_val
            // Restore to base then apply new:
            if let Some(i) = ni {
                mna.system.rhs[*i] = original_rhs[*i] + original_dc - new_value;
            }
            if let Some(j) = nj {
                mna.system.rhs[*j] = original_rhs[*j] - original_dc + new_value;
            }
        }
    }
}

/// Get the original DC value of a source from the netlist.
fn get_source_dc_value(netlist: &Netlist, src_name: &str) -> f64 {
    let src_lower = src_name.to_lowercase();
    for element in netlist.elements() {
        if element.name.to_lowercase() == src_lower {
            match &element.kind {
                ferrospice_netlist::ElementKind::VoltageSource { source, .. }
                | ferrospice_netlist::ElementKind::CurrentSource { source, .. } => {
                    return source
                        .dc
                        .as_ref()
                        .and_then(|e| if let Expr::Num(v) = e { Some(*v) } else { None })
                        .unwrap_or(0.0);
                }
                _ => {}
            }
        }
    }
    0.0
}

/// Perform a DC sweep analysis.
///
/// Sweeps a source across a range of values, computing the DC operating point
/// at each step. Supports single and double (nested) sweeps.
pub fn simulate_dc(netlist: &Netlist) -> Result<SimResult, MnaError> {
    // Find the .dc analysis command in the netlist.
    let (src, start, stop, step, src2) = netlist
        .items
        .iter()
        .find_map(|item| {
            if let Item::Analysis(Analysis::Dc {
                src,
                start,
                stop,
                step,
                src2,
            }) = item
            {
                Some((
                    src.clone(),
                    start.clone(),
                    stop.clone(),
                    step.clone(),
                    src2.clone(),
                ))
            } else {
                None
            }
        })
        .ok_or_else(|| MnaError::UnsupportedElement("no .dc analysis found".to_string()))?;

    let start_val = expr_val(&start, ".dc")?;
    let stop_val = expr_val(&stop, ".dc")?;
    let step_val = expr_val(&step, ".dc")?;

    // Assemble the MNA system.
    let mut mna = assemble_mna(netlist)?;
    let original_rhs = mna.system.rhs.clone();

    // Resolve the primary sweep source.
    let sweep1 = resolve_sweep_source(&mna, netlist, &src)?;
    let original_dc1 = get_source_dc_value(netlist, &src);
    let points1 = generate_sweep_points(start_val, stop_val, step_val);

    // Prepare sweep source vector (the independent variable).
    let sweep_var_name = src.to_lowercase();

    // Initialize result vectors: sweep variable + node voltages + branch currents.
    let mut vecs = Vec::new();
    vecs.push(SimVector {
        name: sweep_var_name,
        real: Vec::new(),
        complex: vec![],
    });
    for (name, _) in mna.node_map.iter() {
        vecs.push(SimVector {
            name: format!("v({})", name),
            real: Vec::new(),
            complex: vec![],
        });
    }
    for vsrc in &mna.vsource_names {
        vecs.push(SimVector {
            name: format!("{}#branch", vsrc.to_lowercase()),
            real: Vec::new(),
            complex: vec![],
        });
    }

    // Resolve optional second sweep source.
    let sweep2_info = if let Some(ref s2) = src2 {
        let sweep2 = resolve_sweep_source(&mna, netlist, &s2.src)?;
        let original_dc2 = get_source_dc_value(netlist, &s2.src);
        let start2 = expr_val(&s2.start, ".dc")?;
        let stop2 = expr_val(&s2.stop, ".dc")?;
        let step2 = expr_val(&s2.step, ".dc")?;
        let points2 = generate_sweep_points(start2, stop2, step2);
        Some((sweep2, original_dc2, points2))
    } else {
        None
    };

    if let Some((ref sweep2, original_dc2, ref points2)) = sweep2_info {
        // Double sweep: outer loop is src2, inner loop is src1
        for &v2 in points2 {
            set_source_value(&mut mna, sweep2, v2, &original_rhs, original_dc2);
            // Save RHS after setting src2 (base for src1 sweep)
            let rhs_after_src2 = mna.system.rhs.clone();
            for &v1 in &points1 {
                set_source_value(&mut mna, &sweep1, v1, &rhs_after_src2, original_dc1);
                vecs[0].real.push(v1);
                collect_solution_into(&mna, &mut vecs[1..])?;
            }
        }
    } else {
        // Single sweep
        for &v1 in &points1 {
            set_source_value(&mut mna, &sweep1, v1, &original_rhs, original_dc1);
            vecs[0].real.push(v1);
            collect_solution_into(&mna, &mut vecs[1..])?;
        }
    }

    Ok(SimResult {
        plots: vec![SimPlot {
            name: "dc1".to_string(),
            vecs,
        }],
    })
}

/// Solve the MNA system and append the solution to result vectors.
fn collect_solution_into(mna: &MnaSystem, vecs: &mut [SimVector]) -> Result<(), MnaError> {
    let solution = mna.solve()?;
    let mut idx = 0;

    for (name, _) in mna.node_map.iter() {
        let v = solution.voltage(name).unwrap_or(0.0);
        vecs[idx].real.push(v);
        idx += 1;
    }

    for vsrc in &mna.vsource_names {
        let i = solution.branch_current(vsrc).unwrap_or(0.0);
        vecs[idx].real.push(i);
        idx += 1;
    }

    Ok(())
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
    use ferrospice_netlist::Netlist;

    /// Helper to get a node voltage from the OP result.
    fn op_voltage(result: &SimResult, node: &str) -> f64 {
        let plot = &result.plots[0];
        let name = format!("v({})", node);
        plot.vecs
            .iter()
            .find(|v| v.name == name)
            .unwrap_or_else(|| panic!("no vector {name}"))
            .real[0]
    }

    /// Helper to get a branch current from the OP result.
    fn op_branch_current(result: &SimResult, vsrc: &str) -> f64 {
        let plot = &result.plots[0];
        let name = format!("{}#branch", vsrc.to_lowercase());
        plot.vecs
            .iter()
            .find(|v| v.name == name)
            .unwrap_or_else(|| panic!("no vector {name}"))
            .real[0]
    }

    #[test]
    fn test_voltage_divider_op() {
        // V1=5V, R1=1k, R2=1k → V(mid) = 2.5V
        let netlist = Netlist::parse(
            "Voltage divider
V1 1 0 5
R1 1 mid 1k
R2 mid 0 1k
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        assert_eq!(result.plots.len(), 1);
        assert_eq!(result.plots[0].name, "op1");

        assert_abs_diff_eq!(op_voltage(&result, "1"), 5.0, epsilon = 1e-9);
        assert_abs_diff_eq!(op_voltage(&result, "mid"), 2.5, epsilon = 1e-9);

        // I(V1) = V/R_total = 5/2000 = 2.5mA (negative by convention)
        let i_v1 = op_branch_current(&result, "V1");
        assert_abs_diff_eq!(i_v1.abs(), 2.5e-3, epsilon = 1e-9);
    }

    #[test]
    fn test_current_source_with_series_resistors() {
        // I1=1mA from 0 to 1, R1=1k (1 to mid), R2=2k (mid to 0)
        // Total R = 3k, V(1) = 1e-3 * 3000 = 3.0V, V(mid) = 1e-3 * 2000 = 2.0V
        let netlist = Netlist::parse(
            "Series resistors with current source
I1 0 1 1m
R1 1 mid 1k
R2 mid 0 2k
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        assert_abs_diff_eq!(op_voltage(&result, "1"), 3.0, epsilon = 1e-9);
        assert_abs_diff_eq!(op_voltage(&result, "mid"), 2.0, epsilon = 1e-9);
    }

    /// Helper to get a vector from a DC sweep result.
    fn dc_vector<'a>(result: &'a SimResult, name: &str) -> &'a Vec<f64> {
        let plot = &result.plots[0];
        &plot
            .vecs
            .iter()
            .find(|v| v.name == name)
            .unwrap_or_else(|| panic!("no vector {name}"))
            .real
    }

    #[test]
    fn test_dc_sweep_voltage_source() {
        // Sweep V1 from 0 to 5V in 1V steps across a 1k resistor
        // Expect 6 data points, I(V1) = -V/R at each point
        let netlist = Netlist::parse(
            "DC sweep test
V1 1 0 0
R1 1 0 1k
.dc V1 0 5 1
.end
",
        )
        .unwrap();

        let result = simulate_dc(&netlist).unwrap();
        assert_eq!(result.plots.len(), 1);
        assert_eq!(result.plots[0].name, "dc1");

        let v1_sweep = dc_vector(&result, "v1");
        assert_eq!(v1_sweep.len(), 6);
        assert_abs_diff_eq!(v1_sweep[0], 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(v1_sweep[5], 5.0, epsilon = 1e-12);

        let v_node = dc_vector(&result, "v(1)");
        assert_eq!(v_node.len(), 6);
        for i in 0..6 {
            let expected_v = i as f64;
            assert_abs_diff_eq!(v_node[i], expected_v, epsilon = 1e-9);
        }

        let i_v1 = dc_vector(&result, "v1#branch");
        assert_eq!(i_v1.len(), 6);
        for i in 0..6 {
            let expected_i = -(i as f64) / 1000.0;
            assert_abs_diff_eq!(i_v1[i], expected_i, epsilon = 1e-9);
        }
    }

    #[test]
    fn test_dc_sweep_double() {
        // Double sweep: V1 from 0 to 2V (step 1), V2 from 0 to 1V (step 1)
        // R1 between node 1 and node 2
        // V1 at node 1, V2 at node 2
        // I = (V1 - V2) / R
        let netlist = Netlist::parse(
            "Double DC sweep
V1 1 0 0
V2 2 0 0
R1 1 2 1k
.dc V1 0 2 1 V2 0 1 1
.end
",
        )
        .unwrap();

        let result = simulate_dc(&netlist).unwrap();

        let v1_sweep = dc_vector(&result, "v1");
        // 3 points for V1 * 2 points for V2 = 6 total
        assert_eq!(v1_sweep.len(), 6);

        // Outer loop is V2, inner is V1
        // V2=0: V1=0,1,2 → I = 0, -1m, -2m
        // V2=1: V1=0,1,2 → I = 1m, 0, -1m
        let i_v1 = dc_vector(&result, "v1#branch");
        assert_eq!(i_v1.len(), 6);

        // V2=0, V1=0: I(V1) = -(0-0)/1k = 0
        assert_abs_diff_eq!(i_v1[0], 0.0, epsilon = 1e-9);
        // V2=0, V1=1: I(V1) = -(1-0)/1k = -1mA
        assert_abs_diff_eq!(i_v1[1], -1e-3, epsilon = 1e-9);
        // V2=0, V1=2: I(V1) = -(2-0)/1k = -2mA
        assert_abs_diff_eq!(i_v1[2], -2e-3, epsilon = 1e-9);
        // V2=1, V1=0: I(V1) = -(0-1)/1k = 1mA
        assert_abs_diff_eq!(i_v1[3], 1e-3, epsilon = 1e-9);
        // V2=1, V1=1: I(V1) = -(1-1)/1k = 0
        assert_abs_diff_eq!(i_v1[4], 0.0, epsilon = 1e-9);
        // V2=1, V1=2: I(V1) = -(2-1)/1k = -1mA
        assert_abs_diff_eq!(i_v1[5], -1e-3, epsilon = 1e-9);
    }

    #[test]
    fn test_res_simple_op() {
        // Port of ngspice-upstream/tests/resistance/res_simple.cir
        // R1=10k, V1=1V DC → I(V1) = -1/10000 = -0.0001A
        let netlist = Netlist::parse(
            "A simple resistor with a voltage source
R1 1 0 10k
V1 1 0 DC 1
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        assert_abs_diff_eq!(op_voltage(&result, "1"), 1.0, epsilon = 1e-9);

        let i_v1 = op_branch_current(&result, "V1");
        assert_abs_diff_eq!(i_v1, -1e-4, epsilon = 1e-9);
    }
}
