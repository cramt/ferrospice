use faer::linalg::solvers::Solve as _;
use faer::Mat;

use thevenin_types::{
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

/// Find roots of det(G + sC) = 0 using Schur complement + generalized eigenvalue decomposition.
///
/// When the C matrix has zero rows (purely algebraic variables such as voltage nodes in
/// inductor-only circuits), the full GEVD is ill-conditioned because the QZ algorithm
/// cannot accurately compute large-magnitude eigenvalues when G and C entries span many
/// orders of magnitude.
///
/// Strategy: identify "reactive" rows (those where C has at least one non-zero entry)
/// and "algebraic" rows (where the entire C row is zero).  Eliminate the algebraic
/// variables via Schur complement, reducing the problem to a smaller, well-conditioned
/// GEVD on the reactive subspace only.
///
/// For circuits with no purely-algebraic rows (e.g. RC only), the full GEVD is used.
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

    // Classify rows into reactive (any C entry non-zero) and algebraic (C row all zero).
    let c_abs_max = {
        let mut m = 0.0_f64;
        for i in 0..dim {
            for j in 0..dim {
                m = m.max(c[(i, j)].abs());
            }
        }
        m
    };
    let c_zero_tol = c_abs_max * 1e-12;

    let is_reactive = |row: usize| -> bool {
        (0..dim).any(|j| c[(row, j)].abs() > c_zero_tol || c[(j, row)].abs() > c_zero_tol)
    };

    let reactive: Vec<usize> = (0..dim).filter(|&i| is_reactive(i)).collect();
    let algebraic: Vec<usize> = (0..dim).filter(|&i| !is_reactive(i)).collect();

    // If there are algebraic variables, attempt Schur complement elimination.
    // This reduces the GEVD to the reactive subspace only, which is much better
    // conditioned for circuits with widely-varying component values (e.g. pz2
    // with R from 1e7 to 1e9 Ω alongside inductors).
    //
    // Prerequisite: Gaa (the algebraic-to-algebraic block of G) must be non-singular.
    // It is singular when algebraic variables are voltage-source branch currents,
    // because G[iv, iv] = 0 for all voltage sources.  In that case we fall back
    // to the full GEVD on the original system.
    if !algebraic.is_empty() && !reactive.is_empty() {
        let n_a = algebraic.len();
        let n_q = reactive.len();

        // Build partitioned blocks.
        //   [Gaa Gaq] [xa]       [0   0  ] [xa]
        //   [Gqa Gqq] [xq] = s * [0   Cqq] [xq]
        let mut gaa = Mat::<f64>::zeros(n_a, n_a);
        let mut gaq = Mat::<f64>::zeros(n_a, n_q);
        let mut gqa = Mat::<f64>::zeros(n_q, n_a);
        let mut gqq = Mat::<f64>::zeros(n_q, n_q);
        let mut cqq = Mat::<f64>::zeros(n_q, n_q);

        for (ia, &i) in algebraic.iter().enumerate() {
            for (ja, &j) in algebraic.iter().enumerate() {
                gaa[(ia, ja)] = g[(i, j)];
            }
            for (jq, &j) in reactive.iter().enumerate() {
                gaq[(ia, jq)] = g[(i, j)];
            }
        }
        for (iq, &i) in reactive.iter().enumerate() {
            for (ja, &j) in algebraic.iter().enumerate() {
                gqa[(iq, ja)] = g[(i, j)];
            }
            for (jq, &j) in reactive.iter().enumerate() {
                gqq[(iq, jq)] = g[(i, j)];
                cqq[(iq, jq)] = c[(i, j)];
            }
        }

        // Guard: if any row of Gaa is all-zero, Gaa is singular (typical for V-source
        // branch current rows where G[iv, iv] = 0). Fall back to full GEVD.
        let g_max = {
            let mut m = 1e-30_f64;
            for i in 0..dim {
                for j in 0..dim {
                    m = m.max(g[(i, j)].abs());
                }
            }
            m
        };
        let zero_row_tol = g_max * 1e-12;
        let gaa_ok = (0..n_a)
            .all(|ia| (0..n_a).any(|ja| gaa[(ia, ja)].abs() > zero_row_tol));

        if gaa_ok {
            // Compute Schur complement: S = Gqq - Gqa * Gaa^{-1} * Gaq.
            // Solve Gaa * X = Gaq for X, then S = Gqq - Gqa * X.
            let gaa_lu = gaa.partial_piv_lu();
            let x = gaa_lu.solve(&gaq); // X = Gaa^{-1} * Gaq, shape n_a × n_q

            // If Gaa was effectively singular (e.g., a node's self-conductance was
            // in the removed column for zeros), the LU solution contains NaN/inf.
            // In that case fall through to full GEVD.
            let x_ok = (0..n_a).all(|ia| (0..n_q).all(|jq| x[(ia, jq)].is_finite()));
            if !x_ok {
                return find_eigenvalue_roots_gevd(g, c, dim);
            }

            let s = &gqq - &gqa * &x; // Schur complement, shape n_q × n_q

            // Reduced eigenvalue problem: det(S + λ·Cqq) = 0.
            // Since Cqq is always invertible (reactive elements have non-zero C),
            // convert to standard form: A = -Cqq^{-1} · S, poles = eigenvalues(A).
            let cqq_lu = cqq.partial_piv_lu();
            let neg_s = &s * faer::Scale(-1.0_f64);
            let a = cqq_lu.solve(&neg_s); // A = Cqq^{-1} · (-S), poles = eig(A)

            // Apply diagonal similarity balancing D^{-1}·A·D so sub-diagonal
            // magnitudes ≤ diagonal magnitudes.  For VCVS-cascaded circuits the
            // raw A has sub-diagonal entries >> diagonal, which prevents QR
            // convergence.  The similarity preserves eigenvalues exactly.
            let a_bal = balance_for_eigenvalues(&a, n_q);

            // Normalise by spectral scale to keep entries O(1).
            let scale = {
                let mut m = 1e-30_f64;
                for i in 0..n_q {
                    m = m.max(a_bal[(i, i)].abs());
                }
                m
            };
            let a_norm = &a_bal * faer::Scale(1.0 / scale);

            let evs = a_norm
                .eigenvalues()
                .map_err(|e| MnaError::UnsupportedElement(format!("EVD failed: {e:?}")))?;

            let mut roots: Vec<(f64, f64)> =
                evs.iter().map(|c| (c.re * scale, c.im * scale)).collect();
            remove_conjugate_duplicates(&mut roots);
            return Ok(roots);
        }
    }

    // Fall back: full GEVD on the original system.
    find_eigenvalue_roots_gevd(g, c, dim)
}

/// Apply a diagonal similarity transformation D^{-1}·A·D that reduces sub-diagonal
/// magnitudes to be ≤ diagonal magnitudes (iterative Parlett-Reinsch-style balancing).
///
/// This does not change eigenvalues but improves QR/QZ convergence for matrices where
/// off-diagonal entries (e.g. from high-gain VCVS stages) dominate the diagonal.
fn balance_for_eigenvalues(a: &Mat<f64>, n: usize) -> Mat<f64> {
    // Start with d[i] = 1 (identity scaling).
    let mut d = vec![1.0_f64; n];

    // Iterate: for each row, find the max sub-diagonal magnitude and scale d[i]
    // so that the sub-diagonal → diagonal ratio becomes ≤ 1.
    // 5 passes is sufficient for practical circuit matrices.
    for _ in 0..10 {
        for i in 0..n {
            let diag = a[(i, i)].abs().max(1e-100);
            for j in 0..n {
                if j == i {
                    continue;
                }
                let off_scaled = a[(i, j)].abs() * d[j] / d[i];
                if off_scaled > diag {
                    // Scale d[i] up to bring ratio to 1.
                    d[i] = a[(i, j)].abs() * d[j] / diag;
                }
            }
        }
    }

    // Build balanced matrix A' = D^{-1} · A · D.
    let mut out = Mat::<f64>::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            out[(i, j)] = a[(i, j)] * d[j] / d[i];
        }
    }
    out
}

/// Core GEVD-based root finder: solves det(G_r + s*C_r) = 0 where both G_r and C_r
/// are the (possibly reduced) matrices with no purely-algebraic zero rows in C_r.
fn find_eigenvalue_roots_gevd(
    g: &Mat<f64>,
    c: &Mat<f64>,
    dim: usize,
) -> Result<Vec<(f64, f64)>, MnaError> {
    if dim == 0 {
        return Ok(vec![]);
    }

    // Special case: faer's generalized_eigen panics on 1×1 matrices.
    if dim == 1 {
        let gv = g[(0, 0)];
        let cv = c[(0, 0)];
        if cv.abs() < 1e-30 {
            return Ok(vec![]);
        }
        return Ok(vec![(-gv / cv, 0.0)]);
    }

    // Generalized eigenvalue decomposition: G*v = λ*(-C)*v.
    // Eigenvalues λ = α/β satisfy det(G + λC) = 0, i.e. λ = pole directly.
    let neg_c = c * faer::Scale(-1.0_f64);

    let gevd = g
        .generalized_eigen(&neg_c)
        .map_err(|e| MnaError::UnsupportedElement(format!("GEVD failed: {e:?}")))?;

    let s_a = gevd.S_a();
    let s_b = gevd.S_b();

    // Threshold for identifying "infinite" eigenvalues (β → 0).
    let mut c_norm = 1e-30_f64;
    for i in 0..dim {
        for j in 0..dim {
            c_norm = c_norm.max(neg_c[(i, j)].abs());
        }
    }
    let threshold = c_norm * f64::EPSILON * 1e6;

    let mut roots = Vec::new();

    for i in 0..dim {
        let alpha = s_a.column_vector()[i];
        let beta = s_b.column_vector()[i];
        let beta_mag = (beta.re * beta.re + beta.im * beta.im).sqrt();

        if beta_mag < threshold {
            continue; // Infinite eigenvalue — skip.
        }

        let alpha_mag = (alpha.re * alpha.re + alpha.im * alpha.im).sqrt();
        if alpha_mag < threshold {
            roots.push((0.0, 0.0));
            continue;
        }

        let denom = beta.re * beta.re + beta.im * beta.im;
        let re = (alpha.re * beta.re + alpha.im * beta.im) / denom;
        let im = (alpha.im * beta.re - alpha.re * beta.im) / denom;

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
    // Snap to real: if |im| is less than 0.1% of the pole magnitude, treat as real.
    // This threshold is deliberately generous: the GEVD (QZ algorithm) may return
    // repeated real poles as complex conjugate pairs with small but non-negligible
    // imaginary noise. Snapping those to real prevents them from being incorrectly
    // deduplicated as "one complex pole" instead of "two real poles".
    let snap_frac = 1e-3;
    for root in roots.iter_mut() {
        let mag = (root.0 * root.0 + root.1 * root.1).sqrt().max(1.0);
        if root.1.abs() < snap_frac * mag {
            root.1 = 0.0;
        }
    }

    // Remove conjugate pairs: for each truly complex root (a, b) with b > 0,
    // remove any matching (a, -b).
    let tol = 1e-6;
    let mut i = 0;
    while i < roots.len() {
        if roots[i].1 > 0.0 {
            // Complex root with positive imaginary part — look for its conjugate.
            let mut j = i + 1;
            while j < roots.len() {
                let (ri, ii) = roots[i];
                let (rj, ij) = roots[j];
                let scale = (ri.abs() + ii.abs()).max(1.0);
                if (ri - rj).abs() < tol * scale && (ii + ij).abs() < tol * scale {
                    // j is the conjugate of i — remove j, keep i (positive im).
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
