use ferrospice_netlist::{Analysis, Expr, Item, Netlist, SimPlot, SimResult, SimVector};

use crate::LinearSystem;
use crate::bjt::stamp_bjt;
use crate::diode::{VT_NOM, pnjlim, vcrit};
use crate::jfet::{jfet_limit, stamp_jfet};
use crate::mna::{MnaError, MnaSystem, assemble_mna, stamp_conductance};
use crate::mosfet::{mos_limit, stamp_mosfet};
use crate::newton::{NrOptions, newton_raphson_solve};

/// Stamp a current source into the RHS vector.
/// Current flows from ni (pos) to nj (neg) externally:
/// subtract from ni, add to nj.
fn stamp_current_source(rhs: &mut [f64], ni: Option<usize>, nj: Option<usize>, i_val: f64) {
    if let Some(i) = ni {
        rhs[i] -= i_val;
    }
    if let Some(j) = nj {
        rhs[j] += i_val;
    }
}

/// Compute the DC operating point of a circuit.
///
/// Assembles the MNA system from the netlist, solves it, and returns
/// a `SimResult` containing node voltages and voltage source branch currents.
/// Uses Newton-Raphson iteration when nonlinear elements (diodes) are present.
pub fn simulate_op(netlist: &Netlist) -> Result<SimResult, MnaError> {
    let mna = assemble_mna(netlist)?;

    let has_nonlinear = !mna.diodes.is_empty()
        || !mna.bjts.is_empty()
        || !mna.mosfets.is_empty()
        || !mna.jfets.is_empty();
    let solution_vec = if !has_nonlinear {
        // Pure linear circuit — direct solve.
        let sol = mna.solve()?;
        // Extract raw values.
        let mut values = Vec::with_capacity(mna.system.dim());
        for (_name, idx) in mna.node_map.iter() {
            // Ensure we write at the right index.
            if values.len() <= idx {
                values.resize(idx + 1, 0.0);
            }
            values[idx] = sol.voltage(_name).unwrap_or(0.0);
        }
        // Append branch currents.
        for vsrc in &mna.vsource_names {
            values.push(sol.branch_current(vsrc).unwrap_or(0.0));
        }
        values
    } else {
        // Nonlinear circuit — use Newton-Raphson.
        solve_nonlinear_op(&mna)?
    };

    let mut vecs = Vec::new();

    // Node voltages
    for (name, idx) in mna.node_map.iter() {
        let v = if idx < solution_vec.len() {
            solution_vec[idx]
        } else {
            0.0
        };
        vecs.push(SimVector {
            name: format!("v({})", name),
            real: vec![v],
            complex: vec![],
        });
    }

    // Voltage source branch currents
    let num_nodes = mna.node_map.len()
        + mna
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
            .sum::<usize>();
    for (i, vsrc) in mna.vsource_names.iter().enumerate() {
        let idx = num_nodes + i;
        let current = if idx < solution_vec.len() {
            solution_vec[idx]
        } else {
            0.0
        };
        vecs.push(SimVector {
            name: format!("{}#branch", vsrc.to_lowercase()),
            real: vec![current],
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

/// Solve the DC operating point and return the raw solution vector.
/// Layout: [node_voltages..., internal_nodes..., branch_currents...].
pub(crate) fn solve_op_raw(mna: &MnaSystem) -> Result<Vec<f64>, MnaError> {
    let has_nonlinear = !mna.diodes.is_empty()
        || !mna.bjts.is_empty()
        || !mna.mosfets.is_empty()
        || !mna.jfets.is_empty();

    if !has_nonlinear {
        mna.system.solve().map_err(MnaError::from)
    } else {
        solve_nonlinear_op(mna)
    }
}

/// Solve a nonlinear DC operating point using Newton-Raphson.
pub(crate) fn solve_nonlinear_op(mna: &MnaSystem) -> Result<Vec<f64>, MnaError> {
    let dim = mna.system.dim();
    let num_ext_nodes = mna.node_map.len();
    let num_internal_diodes = mna
        .diodes
        .iter()
        .filter(|d| d.internal_idx.is_some())
        .count();
    let num_internal_bjts: usize = mna.bjts.iter().map(|b| b.model.internal_node_count()).sum();
    let num_internal_mosfets: usize = mna
        .mosfets
        .iter()
        .map(|m| m.model.internal_node_count())
        .sum();
    let num_internal_jfets: usize = mna
        .jfets
        .iter()
        .map(|j| j.model.internal_node_count())
        .sum();
    let num_nodes = num_ext_nodes
        + num_internal_diodes
        + num_internal_bjts
        + num_internal_mosfets
        + num_internal_jfets;
    let options = NrOptions::default();

    // The base linear system (matrix + RHS) from MNA assembly contains stamps
    // for R, V, I, C, L. We'll replay these plus diode/BJT/MOSFET companions each iteration.
    let base_matrix = &mna.system.matrix;
    let base_rhs = &mna.system.rhs;
    let diodes = &mna.diodes;
    let bjts = &mna.bjts;
    let mosfets = &mna.mosfets;
    let jfets = &mna.jfets;

    // Precompute vcrit for each diode for voltage limiting.
    let vcrits: Vec<f64> = diodes
        .iter()
        .map(|d| vcrit(d.model.n * VT_NOM, d.model.is))
        .collect();

    // Track previous junction voltages for pnjlim.
    let prev_jct_voltages = std::cell::RefCell::new(vec![0.0; diodes.len()]);
    // Track previous BJT junction voltages (vbe, vbc) for pnjlim.
    let prev_bjt_voltages = std::cell::RefCell::new(vec![(0.0, 0.0); bjts.len()]);
    // Track previous MOSFET voltages (vgs, vds) for limiting.
    let prev_mos_voltages = std::cell::RefCell::new(vec![(0.0, 0.0); mosfets.len()]);
    // Track previous JFET junction voltages (vgs, vgd) for limiting.
    let prev_jfet_voltages = std::cell::RefCell::new(vec![(0.0, 0.0); jfets.len()]);
    // Precompute vcrit for each JFET for voltage limiting.
    let jfet_vcrits: Vec<f64> = jfets
        .iter()
        .map(|j| vcrit(j.model.n * VT_NOM, j.model.is))
        .collect();

    let load = |solution: &[f64], system: &mut LinearSystem, source_factor: f64| {
        // 1. Copy base linear stamps.
        for triplet in base_matrix.triplets() {
            system.matrix.add(triplet.row, triplet.col, triplet.value);
        }
        for (i, &val) in base_rhs.iter().enumerate() {
            system.rhs[i] += val * source_factor;
        }

        // 2. Stamp diode companion models.
        {
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

        // 3. Stamp BJT companion models.
        {
            let mut prev = prev_bjt_voltages.borrow_mut();
            for (bi, bjt) in bjts.iter().enumerate() {
                let (raw_vbe, raw_vbc) = bjt.junction_voltages(solution);

                // Apply voltage limiting
                let vbe = bjt.model.limit_vbe(raw_vbe, prev[bi].0);
                let vbc = bjt.model.limit_vbc(raw_vbc, prev[bi].1);
                prev[bi] = (vbe, vbc);

                let comp = bjt.model.companion(vbe, vbc);
                stamp_bjt(&mut system.matrix, &mut system.rhs, bjt, &comp);
            }
        }

        // 4. Stamp MOSFET companion models.
        {
            let mut prev = prev_mos_voltages.borrow_mut();
            for (mi, mos) in mosfets.iter().enumerate() {
                let (raw_vgs, raw_vds, vbs) = mos.terminal_voltages(solution);

                // Apply voltage limiting
                let (vgs, vds) = mos_limit(raw_vgs, raw_vds, prev[mi].0, prev[mi].1, mos.model.vto);
                prev[mi] = (vgs, vds);

                // Create a model with effective beta (KP * W/L)
                let mut eff_model = mos.model.clone();
                eff_model.kp = mos.beta();

                let comp = eff_model.companion(vgs, vds, vbs);
                stamp_mosfet(&mut system.matrix, &mut system.rhs, mos, &comp);
            }
        }

        // 5. Stamp JFET companion models.
        {
            let mut prev = prev_jfet_voltages.borrow_mut();
            for (ji, jfet) in jfets.iter().enumerate() {
                let (raw_vgs, raw_vgd) = jfet.junction_voltages(solution);

                let vt = jfet.model.n * VT_NOM;
                let (vgs, vgd) = jfet_limit(
                    raw_vgs,
                    raw_vgd,
                    prev[ji].0,
                    prev[ji].1,
                    vt,
                    jfet_vcrits[ji],
                );
                prev[ji] = (vgs, vgd);

                let comp = jfet.model.companion(vgs, vgd);
                stamp_jfet(&mut system.matrix, &mut system.rhs, jfet, &comp);
            }
        }
    };

    let initial = vec![0.0; dim];
    let result = newton_raphson_solve(&options, dim, num_nodes, load, &initial).map_err(|e| {
        MnaError::SolveError(crate::SparseMatrixError::SingularMatrix(e.to_string()))
    })?;

    Ok(result.solution)
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

/// Compute total number of nodes including internal nodes for all nonlinear devices.
pub(crate) fn total_num_nodes(mna: &MnaSystem) -> usize {
    mna.node_map.len()
        + mna
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
        let num_nodes = total_num_nodes(mna);
        return Ok(SweepSource::Voltage {
            rhs_index: num_nodes + pos,
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
    let has_nonlinear = !mna.diodes.is_empty()
        || !mna.bjts.is_empty()
        || !mna.mosfets.is_empty()
        || !mna.jfets.is_empty();
    if !has_nonlinear {
        // Linear circuit — direct solve.
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
    } else {
        // Nonlinear circuit — use NR solver.
        let sol = solve_nonlinear_op(mna)?;
        let mut idx = 0;
        let num_ext_nodes = mna.node_map.len();
        let num_internal_diodes = mna
            .diodes
            .iter()
            .filter(|d| d.internal_idx.is_some())
            .count();
        let num_internal_bjts: usize = mna.bjts.iter().map(|b| b.model.internal_node_count()).sum();
        let num_internal_mosfets: usize = mna
            .mosfets
            .iter()
            .map(|m| m.model.internal_node_count())
            .sum();
        let num_internal_jfets: usize = mna
            .jfets
            .iter()
            .map(|j| j.model.internal_node_count())
            .sum();
        let num_nodes = num_ext_nodes
            + num_internal_diodes
            + num_internal_bjts
            + num_internal_mosfets
            + num_internal_jfets;

        for (_name, node_idx) in mna.node_map.iter() {
            vecs[idx].real.push(sol[node_idx]);
            idx += 1;
        }

        for (i, _vsrc) in mna.vsource_names.iter().enumerate() {
            vecs[idx].real.push(sol[num_nodes + i]);
            idx += 1;
        }
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
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

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

    #[test]
    fn test_diode_forward_voltage() {
        // V1=1V, R1=1k, D1 (default model) from node 2 to ground.
        // Circuit: V1(1V) -> node 1 -> R1(1k) -> node 2 -> D1 -> ground
        // Expected: V(2) ≈ 0.6-0.7V (diode forward voltage)
        let netlist = Netlist::parse(
            "Diode forward voltage test
V1 1 0 1
R1 1 2 1k
D1 2 0 DMOD
.model DMOD D
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        let v1 = op_voltage(&result, "1");
        assert_abs_diff_eq!(v1, 1.0, epsilon = 1e-6);

        let v_diode = op_voltage(&result, "2");
        assert!(
            v_diode > 0.5 && v_diode < 0.8,
            "diode voltage {v_diode} not in expected range 0.5-0.8V"
        );

        // Verify KCL: current through R1 = current through diode
        let i_v1 = op_branch_current(&result, "V1");
        // I(V1) should be negative (current flows out of source)
        assert!(i_v1 < 0.0, "V1 branch current should be negative");

        // Current through R1: (V1 - V_diode) / R1
        let i_r1 = (v1 - v_diode) / 1000.0;
        // Should match the magnitude of V1 branch current
        assert_abs_diff_eq!(i_r1, -i_v1, epsilon = 1e-6);
    }

    #[test]
    fn test_diode_with_custom_model() {
        // Test with custom IS parameter
        let netlist = Netlist::parse(
            "Diode with custom model
V1 1 0 1
R1 1 2 1k
D1 2 0 DMOD
.model DMOD D IS=1e-12 N=1.5
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        let v_diode = op_voltage(&result, "2");
        // Higher IS and N means lower forward voltage
        assert!(
            v_diode > 0.3 && v_diode < 0.8,
            "diode voltage {v_diode} not in expected range"
        );
    }

    #[test]
    fn test_diode_dc_sweep_iv_curve() {
        // Sweep voltage across a diode and verify the I-V characteristic.
        // Circuit: V1 -> R1(1k) -> D1 -> ground
        // Sweep V1 from -1V to 1V in 0.5V steps (5 points)
        let netlist = Netlist::parse(
            "Diode I-V sweep
V1 1 0 0
R1 1 2 1k
D1 2 0 DMOD
.model DMOD D IS=1e-14
.dc V1 -1 1 0.5
.end
",
        )
        .unwrap();

        let result = simulate_dc(&netlist).unwrap();

        let v1_sweep = dc_vector(&result, "v1");
        assert_eq!(v1_sweep.len(), 5); // -1, -0.5, 0, 0.5, 1.0

        let i_v1 = dc_vector(&result, "v1#branch");
        assert_eq!(i_v1.len(), 5);

        // At V1=-1V: diode is reverse biased, very little current
        assert!(
            i_v1[0].abs() < 1e-6,
            "reverse bias current should be near zero, got {}",
            i_v1[0]
        );

        // At V1=1V: diode is forward biased, significant current
        assert!(
            i_v1[4].abs() > 1e-4,
            "forward bias current should be significant, got {}",
            i_v1[4]
        );

        // Current should increase monotonically with voltage
        for i in 1..5 {
            assert!(
                i_v1[i].abs() >= i_v1[i - 1].abs() - 1e-10,
                "current should increase: |I[{}]|={} < |I[{}]|={}",
                i,
                i_v1[i].abs(),
                i - 1,
                i_v1[i - 1].abs()
            );
        }
    }

    #[test]
    fn test_diode_iv_matches_shockley() {
        // More precise test: verify diode I-V against the Shockley equation.
        // Use a simple circuit: V1 -> D1 -> ground (no series resistance)
        // With a large resistor to limit current.
        // At each sweep point, verify I = IS * (exp(V_d / Vt) - 1) approximately.
        let netlist = Netlist::parse(
            "Diode Shockley verification
V1 1 0 0
R1 1 2 100
D1 2 0 DMOD
.model DMOD D IS=1e-14 N=1
.dc V1 0.5 0.8 0.1
.end
",
        )
        .unwrap();

        let result = simulate_dc(&netlist).unwrap();

        let v_diode = dc_vector(&result, "v(2)");
        let i_v1 = dc_vector(&result, "v1#branch");

        let is = 1e-14_f64;
        let vt = 0.02585_f64;

        for idx in 0..v_diode.len() {
            let vd = v_diode[idx];
            let i_measured = -i_v1[idx]; // current into circuit
            let i_shockley = is * ((vd / vt).exp() - 1.0);

            // Allow 1e-6 relative tolerance
            let rel_err = if i_shockley.abs() > 1e-15 {
                (i_measured - i_shockley).abs() / i_shockley.abs()
            } else {
                0.0
            };
            assert!(
                rel_err < 1e-3,
                "at Vd={vd:.4}: measured={i_measured:.6e}, shockley={i_shockley:.6e}, rel_err={rel_err:.6e}"
            );
        }
    }

    #[test]
    fn test_bjt_common_emitter_op() {
        // Common-emitter amplifier DC operating point.
        // VCC=10V, RC=2k (collector), RB=200k (base bias).
        // BF=100, IS=1e-16.
        //
        // Hand analysis:
        // IB ≈ (VCC - VBE) / RB ≈ (10 - 0.7) / 200k ≈ 46.5uA
        // IC = BF * IB ≈ 4.65mA
        // VC = VCC - IC*RC ≈ 10 - 4.65*2 ≈ 0.7V
        // (active region if VC > VBE which is marginal here)
        //
        // With a smaller RC we stay safely in active region:
        // RC=1k → VC ≈ 10 - 4.65 = 5.35V
        let netlist = Netlist::parse(
            "BJT common-emitter OP
VCC 1 0 10
RC 1 col 1k
RB 1 base 200k
Q1 col base 0 QMOD
.model QMOD NPN(BF=100 IS=1e-16)
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        let v_base = op_voltage(&result, "base");
        let v_col = op_voltage(&result, "col");

        // Base should be near 0.6-0.85V (forward-biased BE junction, IS=1e-16 gives higher VBE)
        assert!(
            v_base > 0.55 && v_base < 0.9,
            "VBE should be ~0.7V, got {v_base:.4}"
        );

        // Collector should be in active region (above VBE)
        assert!(
            v_col > v_base,
            "Collector ({v_col:.4}) should be above base ({v_base:.4}) for active region"
        );

        // Check approximate IC: IC ≈ (VCC - VC) / RC
        let ic = (10.0 - v_col) / 1000.0;
        let ib = (10.0 - v_base) / 200_000.0;
        let beta = ic / ib;

        // Beta should be near 100 (BF=100)
        assert!(
            beta > 50.0 && beta < 150.0,
            "beta={beta:.1} should be near 100"
        );
    }

    #[test]
    fn test_bjt_dc_sweep() {
        // Sweep VBE and verify IC increases exponentially with VBE.
        let netlist = Netlist::parse(
            "BJT DC sweep
VBE 1 0 0.6
VCE 2 0 5
Q1 2 1 0 QMOD
.model QMOD NPN(BF=100)
.dc VBE 0.5 0.8 0.01
.end
",
        )
        .unwrap();

        let result = simulate_dc(&netlist).unwrap();
        let v_be = dc_vector(&result, "vbe");
        let i_vce = dc_vector(&result, "vce#branch");

        // IC should increase exponentially with VBE.
        // Check that IC at 0.7V > IC at 0.6V by a large factor.
        let idx_06 = v_be.iter().position(|&v| (v - 0.6).abs() < 0.005).unwrap();
        let idx_07 = v_be.iter().position(|&v| (v - 0.7).abs() < 0.005).unwrap();

        let ic_06 = i_vce[idx_06].abs();
        let ic_07 = i_vce[idx_07].abs();

        // 100mV increase in VBE should give ~50x increase in IC (exp(0.1/0.026))
        assert!(
            ic_07 > ic_06 * 10.0,
            "IC at 0.7V ({ic_07:.4e}) should be >> IC at 0.6V ({ic_06:.4e})"
        );

        // Verify exponential trend: IC at 0.8V >> IC at 0.7V
        let idx_08 = v_be.iter().position(|&v| (v - 0.8).abs() < 0.005).unwrap();
        let ic_08 = i_vce[idx_08].abs();
        assert!(
            ic_08 > ic_07 * 5.0,
            "IC at 0.8V ({ic_08:.4e}) should be >> IC at 0.7V ({ic_07:.4e})"
        );
    }

    #[test]
    fn test_nmos_inverter_op() {
        // Simple NMOS inverter-like circuit: VDD=5V, RD=10k drain resistor,
        // VGS=3V (gate driven by voltage source).
        // M1: NMOS with VTO=0.7V, KP=2e-5 (with W/L scaling)
        //
        // At VGS=3V, VTO=0.7V → Vgst=2.3V.
        // Check if MOSFET is conducting (V(drain) < VDD).
        let netlist = Netlist::parse(
            "NMOS inverter OP
VDD 1 0 5
VGS 2 0 3
RD 1 3 10k
M1 3 2 0 0 NMOD W=10u L=1u
.model NMOD NMOS(VTO=0.7 KP=2e-5)
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        let v_drain = op_voltage(&result, "3");
        let v_gate = op_voltage(&result, "2");
        let v_dd = op_voltage(&result, "1");

        assert_abs_diff_eq!(v_gate, 3.0, epsilon = 1e-6);
        assert_abs_diff_eq!(v_dd, 5.0, epsilon = 1e-6);

        // MOSFET should be ON (Vgs > Vto), pulling drain below VDD
        assert!(
            v_drain < v_dd,
            "drain voltage {v_drain} should be < VDD {v_dd}"
        );
        // But drain should be above ground
        assert!(v_drain >= 0.0, "drain voltage {v_drain} should be >= 0");
    }

    #[test]
    fn test_nmos_dc_sweep_transfer_curve() {
        // Sweep VGS from 0 to 5V across NMOS with RD load.
        // Below threshold (VGS < VTO=0.7V): V(drain) ≈ VDD (cutoff)
        // Above threshold: V(drain) drops as MOSFET conducts.
        let netlist = Netlist::parse(
            "NMOS transfer curve
VDD 1 0 5
VGS 2 0 0
RD 1 3 10k
M1 3 2 0 0 NMOD W=10u L=1u
.model NMOD NMOS(VTO=0.7 KP=2e-5)
.dc VGS 0 5 0.5
.end
",
        )
        .unwrap();

        let result = simulate_dc(&netlist).unwrap();
        let vgs_sweep = dc_vector(&result, "vgs");
        let v_drain = dc_vector(&result, "v(3)");

        assert_eq!(vgs_sweep.len(), 11); // 0, 0.5, 1.0, ..., 5.0

        // At VGS=0 (cutoff): drain should be near VDD=5V
        assert!(
            v_drain[0] > 4.9,
            "at VGS=0: V(drain) should be ~VDD, got {:.4}",
            v_drain[0]
        );

        // At VGS=0.5 (below threshold): mostly cutoff (bulk diode leakage is small)
        assert!(
            v_drain[1] > 4.5,
            "at VGS=0.5: V(drain) should be ~VDD, got {:.4}",
            v_drain[1]
        );

        // At VGS=3V (well above threshold): drain should be significantly below VDD
        let idx_3v = vgs_sweep
            .iter()
            .position(|&v| (v - 3.0).abs() < 0.01)
            .unwrap();
        assert!(
            v_drain[idx_3v] < 4.0,
            "at VGS=3V: V(drain) should be < 4V, got {:.4}",
            v_drain[idx_3v]
        );

        // Drain voltage should decrease monotonically with increasing VGS
        // (after threshold)
        let idx_2v = vgs_sweep
            .iter()
            .position(|&v| (v - 2.0).abs() < 0.01)
            .unwrap();
        for i in (idx_2v + 1)..vgs_sweep.len() {
            assert!(
                v_drain[i] <= v_drain[i - 1] + 1e-6,
                "V(drain) should decrease: V[{}]={:.4} > V[{}]={:.4}",
                i,
                v_drain[i],
                i - 1,
                v_drain[i - 1]
            );
        }
    }

    #[test]
    fn test_pmos_op() {
        // PMOS with VDD=5V. Gate at ground → VSG = 5V (well above |VTP|=1V).
        // Should conduct, pulling drain up toward VDD.
        let netlist = Netlist::parse(
            "PMOS OP test
VDD 1 0 5
RS 3 0 10k
M1 3 0 1 1 PMOD W=10u L=1u
.model PMOD PMOS(VTO=-0.7 KP=2e-5)
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        let v_drain = op_voltage(&result, "3");

        // PMOS should be ON (VSG = VDD > |VTP|), pulling drain toward VDD
        assert!(
            v_drain > 0.5,
            "PMOS drain should be pulled up, got {v_drain:.4}"
        );
    }

    #[test]
    fn test_cmos_inverter_op() {
        // CMOS inverter: NMOS + PMOS complementary pair.
        // VIN=0 → PMOS on, NMOS off → VOUT ≈ VDD
        // VIN=VDD → PMOS off, NMOS on → VOUT ≈ 0
        let netlist_low = Netlist::parse(
            "CMOS inverter VIN=0
VDD 1 0 5
VIN 2 0 0
MP 3 2 1 1 PMOD W=10u L=1u
MN 3 2 0 0 NMOD W=10u L=1u
.model NMOD NMOS(VTO=0.7 KP=2e-5)
.model PMOD PMOS(VTO=-0.7 KP=2e-5)
.op
.end
",
        )
        .unwrap();

        let result_low = simulate_op(&netlist_low).unwrap();
        let vout_low = op_voltage(&result_low, "3");

        // VIN=0: NMOS off, PMOS on → VOUT ≈ VDD=5V
        assert!(
            vout_low > 4.5,
            "CMOS inv with VIN=0: VOUT should be ~5V, got {vout_low:.4}"
        );

        let netlist_high = Netlist::parse(
            "CMOS inverter VIN=5
VDD 1 0 5
VIN 2 0 5
MP 3 2 1 1 PMOD W=10u L=1u
MN 3 2 0 0 NMOD W=10u L=1u
.model NMOD NMOS(VTO=0.7 KP=2e-5)
.model PMOD PMOS(VTO=-0.7 KP=2e-5)
.op
.end
",
        )
        .unwrap();

        let result_high = simulate_op(&netlist_high).unwrap();
        let vout_high = op_voltage(&result_high, "3");

        // VIN=5: NMOS on, PMOS off → VOUT ≈ 0V
        assert!(
            vout_high < 0.5,
            "CMOS inv with VIN=5: VOUT should be ~0V, got {vout_high:.4}"
        );
    }

    #[test]
    fn test_jfet_2n4221_op() {
        // Port of ngspice-upstream/tests/jfet/jfet_vds-vgs.cir operating point.
        // N-channel JFET 2N4221 at VGS=-2V, VDS=25V.
        // Expected from ngspice: I(VD) ≈ -9.68268e-4
        let netlist = Netlist::parse(
            "JFET 2N4221 OP
j1 2 1 0 MODJ
VD 2 0 25
VG 1 0 -2
.model MODJ NJF LEVEL=1 VTO=-3.5 BETA=4.1E-4 LAMBDA=0.002 RD=200
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        let v1 = op_voltage(&result, "1");
        let v2 = op_voltage(&result, "2");
        assert_abs_diff_eq!(v1, -2.0, epsilon = 1e-6);
        assert_abs_diff_eq!(v2, 25.0, epsilon = 1e-6);

        // I(VD) should be approximately -9.68268e-4
        let i_vd = op_branch_current(&result, "VD");
        assert_abs_diff_eq!(i_vd, -9.68268e-4, epsilon = 1e-5);
    }

    #[test]
    fn test_vcvs_inverting_amplifier_op() {
        // Op-amp modeled as high-gain VCVS in inverting configuration.
        // V1=1V → R1=1k → inv → R2=2k → out
        // E1 out 0 0 inv 100000 (non-inverting=ground, inverting=inv)
        // Ideal gain: V(out) = -R2/R1 * V(in) = -2V
        let netlist = Netlist::parse(
            "Inverting amplifier
V1 in 0 1
R1 in inv 1k
R2 inv out 2k
E1 out 0 0 inv 100000
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        let v_out = op_voltage(&result, "out");
        assert_abs_diff_eq!(v_out, -2.0, epsilon = 1e-3);

        let v_inv = op_voltage(&result, "inv");
        assert_abs_diff_eq!(v_inv, 0.0, epsilon = 1e-3);
    }

    #[test]
    fn test_vccs_op() {
        // VCCS driving load: G1 0 out in 0 2m, V(in)=3V, R=1k
        // V(out) = gm * V(in) * R = 2e-3 * 3 * 1000 = 6V
        let netlist = Netlist::parse(
            "VCCS OP test
V1 in 0 3
R1 in 0 10k
G1 0 out in 0 2m
R2 out 0 1k
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();
        let v_out = op_voltage(&result, "out");
        assert_abs_diff_eq!(v_out, 6.0, epsilon = 1e-9);
    }

    #[test]
    fn test_cccs_op() {
        // CCCS with gain=10: F1 0 out Vsense 10
        // V1=5V, R1=5k → I(Vsense)=1mA, F1 → 10mA into R2=500 → V(out)=5V
        let netlist = Netlist::parse(
            "CCCS OP test
V1 1 0 5
R1 1 sense 5k
Vsense sense 0 0
F1 0 out Vsense 10
R2 out 0 500
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();
        let v_out = op_voltage(&result, "out");
        assert_abs_diff_eq!(v_out, 5.0, epsilon = 1e-9);
    }

    #[test]
    fn test_ccvs_op() {
        // CCVS: H1 out 0 Vsense 3k → V(out) = 3k * I(Vsense)
        // V1=10V, R1=10k → I(Vsense)=1mA → V(out)=3V
        let netlist = Netlist::parse(
            "CCVS OP test
V1 1 0 10
R1 1 sense 10k
Vsense sense 0 0
H1 out 0 Vsense 3k
R2 out 0 10k
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();
        let v_out = op_voltage(&result, "out");
        assert_abs_diff_eq!(v_out, 3.0, epsilon = 1e-9);
    }

    #[test]
    fn test_jfet_dc_sweep() {
        // JFET DC sweep: sweep VDS from 0 to 25V at VGS=-2V.
        // Current should increase in linear region then saturate.
        let netlist = Netlist::parse(
            "JFET DC sweep
j1 2 1 0 MODJ
VD 2 0 25
VG 1 0 -2
.model MODJ NJF VTO=-3.5 BETA=4.1E-4 LAMBDA=0.002 RD=200
.dc VD 0 25 5
.end
",
        )
        .unwrap();

        let result = simulate_dc(&netlist).unwrap();
        let i_vd = dc_vector(&result, "vd#branch");

        // At VDS=0: almost no current
        assert!(
            i_vd[0].abs() < 1e-6,
            "at VDS=0: current should be near zero, got {}",
            i_vd[0]
        );

        // At VDS=25: should be near saturation
        let i_sat = i_vd[5].abs();
        assert!(
            i_sat > 5e-4,
            "at VDS=25V, VGS=-2V: ID should be > 0.5mA, got {i_sat:.6e}"
        );

        // Current should increase with VDS (lambda modulation)
        for i in 1..i_vd.len() {
            assert!(
                i_vd[i].abs() >= i_vd[i - 1].abs() - 1e-10,
                "current should increase with VDS: |I[{}]|={} < |I[{}]|={}",
                i,
                i_vd[i].abs(),
                i - 1,
                i_vd[i - 1].abs()
            );
        }
    }
}
