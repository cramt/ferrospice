use faer::Mat;

use ferrospice_netlist::{
    Analysis, Complex, ElementKind, Item, Netlist, PzAnalysisType, PzInputType, SimPlot, SimResult,
    SimVector,
};

use crate::mna::{MnaError, MnaSystem, assemble_mna};
use crate::simulate::solve_op_raw;
use crate::tf::build_jacobian;

/// Perform pole-zero analysis (.pz).
///
/// Computes poles and/or zeros of the transfer function between specified
/// input and output ports. Uses the eigenvalue approach:
/// - Build G (conductance Jacobian) and C (capacitance) matrices
/// - Poles = -1/eigenvalues of G^{-1}*C (for non-zero eigenvalues)
/// - Zeros via cofactor submatrix eigenvalues
pub fn simulate_pz(netlist: &Netlist) -> Result<SimResult, MnaError> {
    // Find .pz analysis command.
    let pz_params = netlist
        .items
        .iter()
        .find_map(|item| {
            if let Item::Analysis(Analysis::Pz {
                node_i,
                node_g,
                node_j,
                node_k,
                input_type,
                analysis_type,
            }) = item
            {
                Some((
                    node_i.clone(),
                    node_g.clone(),
                    node_j.clone(),
                    node_k.clone(),
                    *input_type,
                    *analysis_type,
                ))
            } else {
                None
            }
        })
        .ok_or_else(|| MnaError::UnsupportedElement("no .pz analysis found".to_string()))?;

    let (node_i, node_g, node_j, node_k, input_type, analysis_type) = pz_params;

    // Assemble MNA and solve DC operating point.
    let mna = assemble_mna(netlist)?;
    let solution = solve_op_raw(&mna)?;
    let dim = mna.system.dim();

    // Build G matrix (Jacobian at DC OP) and adjust for AC resistance.
    let mut g = build_jacobian(&mna, &solution);

    // Adjust G for resistors with ac= parameter.
    for r in &mna.resistors {
        if let Some(ac_r) = r.ac_resistance
            && (ac_r - r.resistance).abs() > 1e-20
        {
            let g_dc = 1.0 / r.resistance;
            let g_ac = 1.0 / ac_r;
            let delta = g_ac - g_dc;
            stamp_dense_conductance(&mut g, r.pos_idx, r.neg_idx, delta);
        }
    }

    // Build C matrix (capacitance/inductance susceptances without omega).
    let c = build_capacitance_matrix(&mna, &solution);

    // Find poles and/or zeros.
    let mut plots = Vec::new();

    let do_poles = analysis_type == PzAnalysisType::Pol || analysis_type == PzAnalysisType::Pz;
    let do_zeros = analysis_type == PzAnalysisType::Zer || analysis_type == PzAnalysisType::Pz;

    if do_poles {
        let mut poles = find_eigenvalue_roots(&g, &c, dim)?;
        // Sort poles by real part (most negative first).
        poles.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        if !poles.is_empty() {
            let mut pole_vecs = Vec::new();
            for (i, &(re, im)) in poles.iter().enumerate() {
                pole_vecs.push(SimVector {
                    name: format!("pole({})", i + 1),
                    real: vec![],
                    complex: vec![Complex { re, im }],
                });
            }
            plots.push(SimPlot {
                name: "pz1".to_string(),
                vecs: pole_vecs,
            });
        }
    }

    if do_zeros {
        let zeros = find_zeros(
            &g, &c, dim, &mna, netlist, &node_i, &node_g, &node_j, &node_k, input_type,
        )?;

        if !zeros.is_empty() {
            let mut zero_vecs = Vec::new();
            for (i, &(re, im)) in zeros.iter().enumerate() {
                zero_vecs.push(SimVector {
                    name: format!("zero({})", i + 1),
                    real: vec![],
                    complex: vec![Complex { re, im }],
                });
            }
            plots.push(SimPlot {
                name: if do_poles {
                    "pz2".to_string()
                } else {
                    "pz1".to_string()
                },
                vecs: zero_vecs,
            });
        }
    }

    Ok(SimResult { plots })
}

/// Stamp a conductance value into a dense matrix (like stamp_conductance but for Mat<f64>).
fn stamp_dense_conductance(matrix: &mut Mat<f64>, ni: Option<usize>, nj: Option<usize>, g: f64) {
    if let Some(i) = ni {
        matrix[(i, i)] += g;
    }
    if let Some(j) = nj {
        matrix[(j, j)] += g;
    }
    if let (Some(i), Some(j)) = (ni, nj) {
        matrix[(i, j)] -= g;
        matrix[(j, i)] -= g;
    }
}

/// Build the capacitance matrix (C) from MNA system components.
///
/// This is the imaginary part of the AC MNA matrix divided by omega,
/// containing capacitor susceptances, inductor reactances, and
/// nonlinear device junction capacitances.
fn build_capacitance_matrix(mna: &MnaSystem, op_solution: &[f64]) -> Mat<f64> {
    let dim = mna.system.dim();
    let mut c = Mat::<f64>::zeros(dim, dim);

    // Capacitors: C between nodes.
    for cap in &mna.capacitors {
        stamp_dense_conductance(&mut c, cap.pos_idx, cap.neg_idx, cap.capacitance);
    }

    // Inductors: -L on branch diagonal.
    for ind in &mna.inductors {
        c[(ind.branch_idx, ind.branch_idx)] -= ind.inductance;
    }

    // Diode junction capacitances.
    for diode in &mna.diodes {
        let (jct_a, jct_c) = if diode.internal_idx.is_some() {
            (diode.internal_idx, diode.cathode_idx)
        } else {
            (diode.anode_idx, diode.cathode_idx)
        };
        let v_a = jct_a.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_c = jct_c.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_jct = v_a - v_c;

        if diode.model.cjo > 0.0 {
            let cj = diode.model.junction_capacitance(v_jct);
            stamp_dense_conductance(&mut c, jct_a, jct_c, cj);
        }
        if diode.model.tt > 0.0 {
            let (g_d, _) = diode.model.companion(v_jct);
            let cd = diode.model.tt * g_d;
            stamp_dense_conductance(&mut c, jct_a, jct_c, cd);
        }
    }

    // BJT junction capacitances.
    for bjt in &mna.bjts {
        let (vbe, vbc) = bjt.junction_voltages(op_solution);
        let comp = bjt.model.companion(vbe, vbc);
        let m = bjt.m * bjt.area;
        let bp = bjt.base_prime_idx;
        let cp = bjt.col_prime_idx;
        let ep = bjt.emit_prime_idx;

        // B-E junction cap: Cje + Tf*gbe (recover gbe from gpi = gbe/BF)
        let cap_be = bjt.model.cap_be(vbe);
        let gbe = comp.gpi * bjt.model.bf;
        let total_cap_be = cap_be + bjt.model.tf * gbe;
        if total_cap_be > 0.0 {
            stamp_dense_conductance(&mut c, bp, ep, m * total_cap_be);
        }
        // B-C junction cap: Cjc + Tr*gbc (recover gbc from gmu = gmu/BR)
        let cap_bc = bjt.model.cap_bc(vbc);
        let gbc = comp.gmu * bjt.model.br;
        let total_cap_bc = cap_bc + bjt.model.tr * gbc;
        if total_cap_bc > 0.0 {
            stamp_dense_conductance(&mut c, bp, cp, m * total_cap_bc);
        }
        // Substrate cap CJS
        if bjt.model.cjs > 0.0 {
            stamp_dense_conductance(&mut c, bjt.col_idx, None, m * bjt.model.cjs);
        }
    }

    // MOSFET gate overlap and junction capacitances.
    for mos in &mna.mosfets {
        let (vgs, vds, vbs) = mos.terminal_voltages(op_solution);
        let mut eff_model = mos.model.clone();
        eff_model.kp = mos.beta();
        let comp = eff_model.companion(vgs, vds, vbs);
        let _ = comp; // we use direct cap formulas

        let gp = mos.gate_idx;
        let dp = mos.drain_prime_idx;
        let sp = mos.source_prime_idx;
        let bp_idx = mos.bulk_idx;

        // Gate overlap caps
        let cgso = mos.model.cgso * mos.w;
        let cgdo = mos.model.cgdo * mos.w;
        let cgbo = mos.model.cgbo * mos.l;
        if cgso > 0.0 {
            stamp_dense_conductance(&mut c, gp, sp, cgso);
        }
        if cgdo > 0.0 {
            stamp_dense_conductance(&mut c, gp, dp, cgdo);
        }
        if cgbo > 0.0 {
            stamp_dense_conductance(&mut c, gp, bp_idx, cgbo);
        }

        // Junction caps
        if mos.model.cbd > 0.0 {
            stamp_dense_conductance(&mut c, dp, bp_idx, mos.model.cbd);
        }
        if mos.model.cbs > 0.0 {
            stamp_dense_conductance(&mut c, sp, bp_idx, mos.model.cbs);
        }
    }

    // JFET gate capacitances.
    for jfet in &mna.jfets {
        let gp = jfet.gate_idx;
        let dp = jfet.drain_prime_idx;
        let sp = jfet.source_prime_idx;

        if jfet.model.cgs > 0.0 {
            stamp_dense_conductance(&mut c, gp, sp, jfet.model.cgs);
        }
        if jfet.model.cgd > 0.0 {
            stamp_dense_conductance(&mut c, gp, dp, jfet.model.cgd);
        }
    }

    c
}

/// Find roots of det(G + sC) = 0 using eigenvalue decomposition.
///
/// Two approaches based on whether G is invertible:
/// 1. If G is invertible: eigenvalues of G^{-1}*C, poles = -1/λ for non-zero λ
/// 2. If G is singular: use generalized eigenvalue decomposition
fn find_eigenvalue_roots(
    g: &Mat<f64>,
    c: &Mat<f64>,
    dim: usize,
) -> Result<Vec<(f64, f64)>, MnaError> {
    if dim == 0 {
        return Ok(vec![]);
    }

    // Special case: 1×1 matrix.
    if dim == 1 {
        let gv = g[(0, 0)];
        let cv = c[(0, 0)];
        if cv.abs() < 1e-30 {
            return Ok(vec![]); // No frequency-dependent part.
        }
        // g + s*c = 0 → s = -g/c
        return Ok(vec![(-gv / cv, 0.0)]);
    }

    // Try the standard approach: A = G^{-1} * C, then eigenvalues of A.
    // Check if G is well-conditioned enough by trying to solve.
    use faer::linalg::solvers::FullPivLu;
    use faer::prelude::Solve;

    let lu = FullPivLu::new(g.as_ref());

    // Check if G is approximately singular by examining the LU factor.
    // We test by solving G*x = e_1 and checking for NaN/Inf.
    let mut test_rhs = Mat::<f64>::zeros(dim, 1);
    test_rhs[(0, 0)] = 1.0;
    let test_sol = lu.solve(&test_rhs);
    let g_is_invertible = (0..dim).all(|i| test_sol[(i, 0)].is_finite());

    if g_is_invertible {
        let a = lu.solve(c);

        let eigenvalues = a.eigenvalues().map_err(|e| {
            MnaError::UnsupportedElement(format!("eigenvalue computation failed: {e:?}"))
        })?;

        let mut roots = Vec::new();

        // Use a relative threshold: eigenvalues smaller than eps * max_eigenvalue
        // are considered zero (numerical noise from the matrix solve).
        let max_mag = eigenvalues
            .iter()
            .map(|ev| (ev.re * ev.re + ev.im * ev.im).sqrt())
            .fold(0.0_f64, f64::max);
        let threshold = max_mag * 1e-10;

        for ev in &eigenvalues {
            let mag = (ev.re * ev.re + ev.im * ev.im).sqrt();
            if mag > threshold {
                // pole = -1/lambda
                let mag_sq = mag * mag;
                let re = -ev.re / mag_sq;
                let im = ev.im / mag_sq;
                roots.push((re, im));
            }
        }

        // Remove conjugate pairs (keep one representative).
        remove_conjugate_duplicates(&mut roots);

        return Ok(roots);
    }

    // G is singular — use generalized eigenvalue decomposition.
    let neg_c = c * faer::Scale(-1.0_f64);

    let gevd = g
        .generalized_eigen(&neg_c)
        .map_err(|e| MnaError::UnsupportedElement(format!("GEVD failed: {e:?}")))?;

    let s_a = gevd.S_a();
    let s_b = gevd.S_b();

    let mut roots = Vec::new();
    let threshold = 1e-15;

    for i in 0..dim {
        let alpha = s_a.column_vector()[i];
        let beta = s_b.column_vector()[i];
        let beta_mag = (beta.re * beta.re + beta.im * beta.im).sqrt();

        if beta_mag < threshold {
            continue; // Infinite eigenvalue — skip.
        }

        let denom = beta.re * beta.re + beta.im * beta.im;
        let re = (alpha.re * beta.re + alpha.im * beta.im) / denom;
        let im = (alpha.im * beta.re - alpha.re * beta.im) / denom;

        let alpha_mag = (alpha.re * alpha.re + alpha.im * alpha.im).sqrt();
        if alpha_mag < threshold * beta_mag {
            roots.push((0.0, 0.0));
            continue;
        }

        roots.push((re, im));
    }

    remove_conjugate_duplicates(&mut roots);

    Ok(roots)
}

/// Remove conjugate duplicates from root list.
/// When faer returns both members of a conjugate pair (a+bi) and (a-bi),
/// keep only one (the one with positive imaginary part).
/// Does NOT remove repeated real roots (which are physically meaningful).
fn remove_conjugate_duplicates(roots: &mut Vec<(f64, f64)>) {
    let tol = 1e-8;

    // First snap near-real roots to real.
    for root in roots.iter_mut() {
        if root.1.abs() < tol * root.0.abs().max(1.0) {
            root.1 = 0.0;
        }
    }

    // Remove conjugate pairs: for each complex root (a, b) with b > 0,
    // remove any matching (a, -b).
    let mut i = 0;
    while i < roots.len() {
        if roots[i].1.abs() > tol {
            // Complex root — look for its conjugate.
            let mut j = i + 1;
            while j < roots.len() {
                let (ri, ii) = roots[i];
                let (rj, ij) = roots[j];
                let scale = ri.abs().max(ii.abs()).max(1.0);
                if (ri - rj).abs() < tol * scale && (ii + ij).abs() < tol * scale {
                    // j is the conjugate of i — remove j, keep the one with im > 0.
                    if ii < 0.0 {
                        roots[i].1 = -roots[i].1; // keep positive im
                    }
                    roots.remove(j);
                    break;
                }
                j += 1;
            }
        }
        i += 1;
    }
}

/// Find zeros of the transfer function using cofactor submatrix eigenvalues.
#[expect(clippy::too_many_arguments)]
fn find_zeros(
    g: &Mat<f64>,
    c: &Mat<f64>,
    dim: usize,
    mna: &MnaSystem,
    netlist: &Netlist,
    node_i: &str,
    node_g: &str,
    node_j: &str,
    node_k: &str,
    input_type: PzInputType,
) -> Result<Vec<(f64, f64)>, MnaError> {
    let num_nodes = mna.total_num_nodes();

    // Determine which row to remove (input side).
    let remove_row = match input_type {
        PzInputType::Cur => {
            // Current input: remove row of input node.
            resolve_node_index(node_i, node_g, mna)
        }
        PzInputType::Vol => {
            // Voltage input: remove row of the voltage source branch
            // connecting the input nodes.
            find_input_vsource_branch(node_i, node_g, mna, netlist, num_nodes)
        }
    }?;

    // Determine which column to remove (output side).
    let remove_col = resolve_node_index(node_j, node_k, mna)?;

    if dim < 2 {
        return Ok(vec![]);
    }

    // Build submatrices by removing the specified row and column.
    let sub_dim = dim - 1;
    let mut g_sub = Mat::<f64>::zeros(sub_dim, sub_dim);
    let mut c_sub = Mat::<f64>::zeros(sub_dim, sub_dim);

    let mut si = 0;
    for i in 0..dim {
        if i == remove_row {
            continue;
        }
        let mut sj = 0;
        for j in 0..dim {
            if j == remove_col {
                continue;
            }
            g_sub[(si, sj)] = g[(i, j)];
            c_sub[(si, sj)] = c[(i, j)];
            sj += 1;
        }
        si += 1;
    }

    let mut zeros = find_eigenvalue_roots(&g_sub, &c_sub, sub_dim)?;
    // Sort zeros by real part.
    zeros.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    Ok(zeros)
}

/// Resolve a node name to its matrix index.
/// For ground reference (node "0"), returns None in the node_map,
/// so we use the other node directly.
fn resolve_node_index(node_pos: &str, node_neg: &str, mna: &MnaSystem) -> Result<usize, MnaError> {
    let pos_is_ground = node_pos == "0";
    let neg_is_ground = node_neg == "0";

    if pos_is_ground && neg_is_ground {
        return Err(MnaError::UnsupportedElement(
            "both nodes are ground in PZ analysis".to_string(),
        ));
    }

    if neg_is_ground {
        // Single-ended: use positive node index.
        mna.node_map.get(&node_pos.to_lowercase()).ok_or_else(|| {
            MnaError::UnsupportedElement(format!("node '{node_pos}' not found in circuit"))
        })
    } else if pos_is_ground {
        mna.node_map.get(&node_neg.to_lowercase()).ok_or_else(|| {
            MnaError::UnsupportedElement(format!("node '{node_neg}' not found in circuit"))
        })
    } else {
        // Differential: use the positive node.
        // TODO: proper differential handling would need both nodes.
        mna.node_map.get(&node_pos.to_lowercase()).ok_or_else(|| {
            MnaError::UnsupportedElement(format!("node '{node_pos}' not found in circuit"))
        })
    }
}

/// Find the branch index of a voltage source connecting the input nodes.
fn find_input_vsource_branch(
    node_i: &str,
    node_g: &str,
    mna: &MnaSystem,
    netlist: &Netlist,
    num_nodes: usize,
) -> Result<usize, MnaError> {
    let node_i_lower = node_i.to_lowercase();
    let node_g_lower = node_g.to_lowercase();

    // Search through netlist elements for a voltage source at the input nodes.
    for elem in netlist.elements() {
        if let ElementKind::VoltageSource { pos, neg, .. } = &elem.kind {
            let pos_lower = pos.to_lowercase();
            let neg_lower = neg.to_lowercase();
            if (pos_lower == node_i_lower && neg_lower == node_g_lower)
                || (pos_lower == node_g_lower && neg_lower == node_i_lower)
            {
                // Found the voltage source. Get its branch index.
                if let Some(vs_idx) = mna
                    .vsource_names
                    .iter()
                    .position(|n| n.eq_ignore_ascii_case(&elem.name))
                {
                    return Ok(num_nodes + vs_idx);
                }
            }
        }
    }

    Err(MnaError::UnsupportedElement(format!(
        "no voltage source found between nodes '{node_i}' and '{node_g}' for PZ voltage input"
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use approx::assert_relative_eq;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    fn get_poles(result: &SimResult) -> Vec<(f64, f64)> {
        let mut poles = Vec::new();
        for plot in &result.plots {
            for vec in &plot.vecs {
                if vec.name.starts_with("pole(") {
                    if let Some(c) = vec.complex.first() {
                        poles.push((c.re, c.im));
                    }
                }
            }
        }
        poles
    }

    fn get_zeros(result: &SimResult) -> Vec<(f64, f64)> {
        let mut zeros = Vec::new();
        for plot in &result.plots {
            for vec in &plot.vecs {
                if vec.name.starts_with("zero(") {
                    if let Some(c) = vec.complex.first() {
                        zeros.push((c.re, c.im));
                    }
                }
            }
        }
        zeros
    }

    #[test]
    fn test_pz_simple_rc_poles() {
        // Simple RC filter: R1=1k to ground, R2=1k to ground, C1=1pF between nodes.
        // Expected: 1 pole at -5e8.
        let netlist = Netlist::parse(
            "simple pz test
r1 1 0 1k
r2 2 0 1k
c1 1 2 1.0e-12
.pz 1 0 2 0 cur pz
.end
",
        )
        .unwrap();

        let result = simulate_pz(&netlist).unwrap();
        let poles = get_poles(&result);
        assert_eq!(poles.len(), 1);
        assert_relative_eq!(poles[0].0, -5.0e8, max_relative = 1e-6);
        assert_abs_diff_eq!(poles[0].1, 0.0, epsilon = 1.0);

        let zeros = get_zeros(&result);
        assert_eq!(zeros.len(), 1);
        assert_abs_diff_eq!(zeros[0].0, 0.0, epsilon = 1.0);
        assert_abs_diff_eq!(zeros[0].1, 0.0, epsilon = 1.0);
    }

    #[test]
    fn test_pz_rc_lowpass_voltage() {
        // RC lowpass: V1 at input, R1=1k, C1=10p to ground.
        // Expected: 1 pole at -1e8.
        let netlist = Netlist::parse(
            "RC filter
v1 1 0 0 ac 1.0
r1 1 2 1k
c1 2 0 10p
.pz 1 0 2 0 vol pz
.end
",
        )
        .unwrap();

        let result = simulate_pz(&netlist).unwrap();
        let poles = get_poles(&result);
        assert_eq!(poles.len(), 1);
        assert_relative_eq!(poles[0].0, -1.0e8, max_relative = 1e-6);
        assert_abs_diff_eq!(poles[0].1, 0.0, epsilon = 1.0);

        // No zeros expected for this simple lowpass.
        let zeros = get_zeros(&result);
        assert_eq!(zeros.len(), 0);
    }

    #[test]
    fn test_pz_bridge_t() {
        // Bridge-T filter with 2 poles and 2 zeros.
        let netlist = Netlist::parse(
            "BRIDGE-T FILTER
V1 1 0 12 AC 1
C1 1 2 1U
C2 2 3 1U
R3 2 0 1K
R4 1 3 1K
.PZ 1 0 3 0 VOL PZ
.end
",
        )
        .unwrap();

        let result = simulate_pz(&netlist).unwrap();
        let poles = get_poles(&result);
        assert_eq!(poles.len(), 2);

        // Sort by magnitude of real part.
        let mut sorted_poles = poles.clone();
        sorted_poles.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        assert_relative_eq!(sorted_poles[0].0, -2618.034, max_relative = 1e-4);
        assert_relative_eq!(sorted_poles[1].0, -381.966, max_relative = 1e-4);

        let zeros = get_zeros(&result);
        assert_eq!(zeros.len(), 2);
        for z in &zeros {
            assert_relative_eq!(z.0, -1000.0, max_relative = 1e-4);
            assert_abs_diff_eq!(z.1, 0.0, epsilon = 1.0);
        }
    }

    #[test]
    fn test_pz_multistage() {
        // Three-stage VCVS cascade with 3 poles.
        let netlist = Netlist::parse(
            "Multistage filter
v1 1 0 0 ac 1.0
r1 1 2 1k
c1 2 0 10p
e2 3 0 2 0 10
r2 3 4 1k
c2 4 0 1.25p
e3 5 0 4 0 10
r3 5 6 1k
c3 6 0 .02p
.pz 1 0 6 0 vol pz
.end
",
        )
        .unwrap();

        let result = simulate_pz(&netlist).unwrap();
        let poles = get_poles(&result);
        assert_eq!(poles.len(), 3);

        let mut sorted_poles: Vec<f64> = poles.iter().map(|p| p.0).collect();
        sorted_poles.sort_by(|a, b| a.partial_cmp(b).unwrap());

        assert_relative_eq!(sorted_poles[0], -5.0e10, max_relative = 1e-4);
        assert_relative_eq!(sorted_poles[1], -8.0e8, max_relative = 1e-4);
        assert_relative_eq!(sorted_poles[2], -1.0e8, max_relative = 1e-4);
    }

    #[test]
    fn test_pz_poles_only() {
        // Test POL (poles only) mode.
        let netlist = Netlist::parse(
            "test pz
iin 1 0 ac
r1 1 0 1
l1 1 0 0.05
gm2 2 0 1 0 1
r2 2 0 1
l2 2 0 0.05
gm3 3 0 2 0 1
r3 3 0 1
l3 3 0 0.05
.pz 1 0 3 0 cur pol
.end
",
        )
        .unwrap();

        let result = simulate_pz(&netlist).unwrap();
        let poles = get_poles(&result);
        assert_eq!(poles.len(), 3);

        for p in &poles {
            assert_relative_eq!(p.0, -20.0, max_relative = 1e-4);
            assert_abs_diff_eq!(p.1, 0.0, epsilon = 1.0);
        }

        // Should be no zeros (POL mode).
        let zeros = get_zeros(&result);
        assert_eq!(zeros.len(), 0);
    }

    #[test]
    fn test_pz_ac_resistance() {
        // AC resistance: R1 has dc=1k, ac=4k. PZ should use ac value.
        // tau = (4k || 4k) * 1n = 2k * 1n = 2e-6
        // pole = -1/tau = -5e5
        let netlist = Netlist::parse(
            "ac resistance test
Vin 1 0 dc 5.0 ac 3.0
R1 1 2 1k ac=4k
C2 2 0 1n
R2 2 0 4k
.pz 1 0 2 0 vol pz
.end
",
        )
        .unwrap();

        let result = simulate_pz(&netlist).unwrap();
        let poles = get_poles(&result);
        assert_eq!(poles.len(), 1);
        assert_relative_eq!(poles[0].0, -5.0e5, max_relative = 1e-6);
    }
}
