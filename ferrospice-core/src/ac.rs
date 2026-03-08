use std::f64::consts::PI;

use ferrospice_netlist::{
    AcVariation, Analysis, Complex, Expr, Item, Netlist, SimPlot, SimResult, SimVector,
};

use crate::mna::{MnaError, MnaSystem, assemble_mna};
use crate::simulate::simulate_op;
use crate::sparse::ComplexLinearSystem;

/// Perform AC small-signal analysis.
///
/// 1. Compute DC operating point to linearize nonlinear devices.
/// 2. Build the complex MNA matrix (G + jωC) at each frequency.
/// 3. Apply AC source excitation.
/// 4. Solve for complex node voltages.
pub fn simulate_ac(netlist: &Netlist) -> Result<SimResult, MnaError> {
    // Find the .ac analysis command.
    let (variation, n, fstart, fstop) = netlist
        .items
        .iter()
        .find_map(|item| {
            if let Item::Analysis(Analysis::Ac {
                variation,
                n,
                fstart,
                fstop,
            }) = item
            {
                Some((*variation, *n, fstart.clone(), fstop.clone()))
            } else {
                None
            }
        })
        .ok_or_else(|| MnaError::UnsupportedElement("no .ac analysis found".to_string()))?;

    let fstart_val = expr_val(&fstart, ".ac")?;
    let fstop_val = expr_val(&fstop, ".ac")?;

    // Compute DC operating point (linearizes nonlinear devices).
    let op_result = simulate_op(netlist)?;

    // Assemble the DC MNA system for its structure.
    let mna = assemble_mna(netlist)?;

    // Extract DC operating point solution for linearization.
    let op_solution = extract_op_solution(&op_result, &mna);

    // Generate frequency sweep points.
    let frequencies = generate_ac_sweep(variation, n, fstart_val, fstop_val);

    // Build result vectors.
    let mut freq_vec = SimVector {
        name: "frequency".to_string(),
        real: Vec::with_capacity(frequencies.len()),
        complex: vec![],
    };

    let mut node_vecs: Vec<SimVector> = mna
        .node_map
        .iter()
        .map(|(name, _)| SimVector {
            name: format!("v({})", name),
            real: vec![],
            complex: Vec::with_capacity(frequencies.len()),
        })
        .collect();

    let mut branch_vecs: Vec<SimVector> = mna
        .vsource_names
        .iter()
        .map(|name| SimVector {
            name: format!("{}#branch", name.to_lowercase()),
            real: vec![],
            complex: Vec::with_capacity(frequencies.len()),
        })
        .collect();

    // Sweep each frequency point.
    for &freq in &frequencies {
        let omega = 2.0 * PI * freq;

        // Build and solve the complex MNA system at this frequency.
        let solution = solve_ac_point(&mna, &op_solution, omega, netlist)?;

        freq_vec.real.push(freq);

        // Collect node voltages (complex).
        for (i, (_name, node_idx)) in mna.node_map.iter().enumerate() {
            let (re, im) = solution[node_idx];
            node_vecs[i].complex.push(Complex { re, im });
        }

        // Collect branch currents (complex).
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
                .sum::<usize>();
        for (i, _vsrc) in mna.vsource_names.iter().enumerate() {
            let idx = num_nodes + i;
            let (re, im) = solution[idx];
            branch_vecs[i].complex.push(Complex { re, im });
        }
    }

    let mut vecs = vec![freq_vec];
    vecs.extend(node_vecs);
    vecs.extend(branch_vecs);

    Ok(SimResult {
        plots: vec![SimPlot {
            name: "ac1".to_string(),
            vecs,
        }],
    })
}

/// Solve the AC MNA system at a single frequency point.
///
/// Builds the complex matrix: real part = G (conductances + voltage source stamps),
/// imaginary part = ωC (susceptances from capacitors/inductors).
fn solve_ac_point(
    mna: &MnaSystem,
    op_solution: &[f64],
    omega: f64,
    netlist: &Netlist,
) -> Result<Vec<(f64, f64)>, MnaError> {
    let num_ext_nodes = mna.node_map.len();
    let num_internal = mna
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
            .sum::<usize>();
    let num_nodes = num_ext_nodes + num_internal;
    let dim = num_nodes + mna.vsource_names.len();

    let mut sys = ComplexLinearSystem::new(dim);

    // 1. Stamp DC conductance matrix (real part) from base MNA stamps.
    // This includes resistors and voltage source topology.
    for triplet in mna.system.matrix.triplets() {
        sys.real.add(triplet.row, triplet.col, triplet.value);
    }

    // 2. Stamp capacitor AC admittance: jωC between nodes.
    for cap in &mna.capacitors {
        let bc = omega * cap.capacitance; // susceptance = ωC
        stamp_imag_conductance(&mut sys.imag, cap.pos_idx, cap.neg_idx, bc);
    }

    // 3. Stamp inductor AC impedance: jωL on the branch equation.
    // Inductor topology stamps (±1 entries) are already in the real matrix from DC.
    // For AC, the branch equation becomes: V(pos) - V(neg) = jωL * I_branch
    // → V(pos) - V(neg) - jωL * I_branch = 0
    // The -jωL goes on the diagonal of the branch equation (imaginary part).
    for ind in &mna.inductors {
        sys.imag
            .add(ind.branch_idx, ind.branch_idx, -omega * ind.inductance);
    }

    // 4. Stamp diode small-signal conductance at DC operating point.
    for diode in &mna.diodes {
        let (jct_anode, jct_cathode) = if diode.internal_idx.is_some() {
            (diode.internal_idx, diode.cathode_idx)
        } else {
            (diode.anode_idx, diode.cathode_idx)
        };

        // Get junction voltage from DC solution.
        let v_anode = jct_anode.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_cathode = jct_cathode.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_jct = v_anode - v_cathode;

        // Small-signal conductance: dI/dV at operating point.
        let (g_d, _i_eq) = diode.model.companion(v_jct);
        crate::stamp_conductance(&mut sys.real, jct_anode, jct_cathode, g_d);

        // Series resistance (already in DC stamps if RS > 0).
        if let Some(int_idx) = diode.internal_idx {
            let g_rs = 1.0 / diode.model.rs;
            crate::stamp_conductance(&mut sys.real, diode.anode_idx, Some(int_idx), g_rs);
        }

        // Junction capacitance: jω * Cj
        if diode.model.cjo > 0.0 {
            let cj = diode.model.junction_capacitance(v_jct);
            let bc = omega * cj;
            stamp_imag_conductance(&mut sys.imag, jct_anode, jct_cathode, bc);
        }

        // Transit time capacitance: jω * τ * gd
        if diode.model.tt > 0.0 {
            let cd = diode.model.tt * g_d;
            let bc = omega * cd;
            stamp_imag_conductance(&mut sys.imag, jct_anode, jct_cathode, bc);
        }
    }

    // 5. Stamp BJT small-signal model at DC operating point.
    for bjt in &mna.bjts {
        let (vbe, vbc) = bjt.junction_voltages(op_solution);
        let comp = bjt.model.companion(vbe, vbc);
        let m = bjt.m * bjt.area;

        let bp = bjt.base_prime_idx;
        let cp = bjt.col_prime_idx;
        let ep = bjt.emit_prime_idx;

        // gpi conductance b'-e' (real)
        crate::stamp_conductance(&mut sys.real, bp, ep, m * comp.gpi);
        // gmu conductance b'-c' (real)
        crate::stamp_conductance(&mut sys.real, bp, cp, m * comp.gmu);
        // go conductance c'-e' (real)
        crate::stamp_conductance(&mut sys.real, cp, ep, m * comp.go);

        // gm VCCS: V_be controls current from e' to c' (real)
        let gm_scaled = m * comp.gm;
        if let Some(c) = cp {
            if let Some(b) = bp {
                sys.real.add(c, b, gm_scaled);
            }
            if let Some(e) = ep {
                sys.real.add(c, e, -gm_scaled);
            }
        }
        if let Some(e) = ep {
            if let Some(b) = bp {
                sys.real.add(e, b, -gm_scaled);
            }
            sys.real.add(e, e, gm_scaled);
        }

        // Series resistances (real)
        if bjt.model.rb > 0.0 {
            let gx = 1.0 / bjt.model.rb;
            crate::stamp_conductance(&mut sys.real, bjt.base_idx, bp, m * gx);
        }
        if bjt.model.rc > 0.0 {
            let gcpr = 1.0 / bjt.model.rc;
            crate::stamp_conductance(&mut sys.real, bjt.col_idx, cp, m * gcpr);
        }
        if bjt.model.re > 0.0 {
            let gepr = 1.0 / bjt.model.re;
            crate::stamp_conductance(&mut sys.real, bjt.emit_idx, ep, m * gepr);
        }

        // B-E junction capacitance (imaginary): jω * (Cje + Tf*gbe)
        let cap_be = bjt.model.cap_be(vbe);
        let gbe = comp.gpi * bjt.model.bf; // recover gbe from gpi = gbe/BF
        let total_cap_be = cap_be + bjt.model.tf * gbe;
        if total_cap_be > 0.0 {
            stamp_imag_conductance(&mut sys.imag, bp, ep, omega * m * total_cap_be);
        }

        // B-C junction capacitance (imaginary): jω * (Cjc + Tr*gbc)
        let cap_bc = bjt.model.cap_bc(vbc);
        let gbc = comp.gmu * bjt.model.br; // recover gbc from gmu = gbc/BR
        let total_cap_bc = cap_bc + bjt.model.tr * gbc;
        if total_cap_bc > 0.0 {
            stamp_imag_conductance(&mut sys.imag, bp, cp, omega * m * total_cap_bc);
        }
    }

    // 5b. Stamp MOSFET small-signal model at DC operating point.
    for mos in &mna.mosfets {
        let (vgs, vds, vbs) = mos.terminal_voltages(op_solution);

        // Compute companion with effective beta
        let mut eff_model = mos.model.clone();
        eff_model.kp = mos.beta();
        let comp = eff_model.companion(vgs, vds, vbs);
        let m = mos.m;

        let dp = mos.drain_prime_idx;
        let g = mos.gate_idx;
        let sp = mos.source_prime_idx;
        let b = mos.bulk_idx;

        let (xnrm, xrev) = if comp.mode > 0 {
            (1.0, 0.0)
        } else {
            (0.0, 1.0)
        };

        // gds conductance d'-s' (real)
        crate::stamp_conductance(&mut sys.real, dp, sp, m * comp.gds);

        // gm VCCS: Vgs controls current s'→d' (real)
        let gm_scaled = m * comp.gm;
        if let Some(d) = dp {
            if let Some(gate) = g {
                sys.real.add(d, gate, (xnrm - xrev) * gm_scaled);
            }
            if let Some(s) = sp {
                sys.real.add(d, s, -(xnrm - xrev) * gm_scaled);
            }
        }
        if let Some(s) = sp {
            if let Some(gate) = g {
                sys.real.add(s, gate, -(xnrm - xrev) * gm_scaled);
            }
            sys.real.add(s, s, (xnrm - xrev) * gm_scaled);
        }

        // gmbs: Vbs controls current s'→d' (real)
        let gmbs_scaled = m * comp.gmbs;
        if let Some(d) = dp {
            if let Some(bulk) = b {
                sys.real.add(d, bulk, (xnrm - xrev) * gmbs_scaled);
            }
            if let Some(s) = sp {
                sys.real.add(d, s, -(xnrm - xrev) * gmbs_scaled);
            }
        }
        if let Some(s) = sp {
            if let Some(bulk) = b {
                sys.real.add(s, bulk, -(xnrm - xrev) * gmbs_scaled);
            }
            sys.real.add(s, s, (xnrm - xrev) * gmbs_scaled);
        }

        // gbd conductance b-d' (real)
        crate::stamp_conductance(&mut sys.real, b, dp, m * comp.gbd);

        // gbs conductance b-s' (real)
        crate::stamp_conductance(&mut sys.real, b, sp, m * comp.gbs);

        // Series resistances (real)
        if mos.model.rd > 0.0 {
            let grd = 1.0 / mos.model.rd;
            crate::stamp_conductance(&mut sys.real, mos.drain_idx, dp, m * grd);
        }
        if mos.model.rs > 0.0 {
            let grs = 1.0 / mos.model.rs;
            crate::stamp_conductance(&mut sys.real, mos.source_idx, sp, m * grs);
        }

        // Gate overlap capacitances (imaginary)
        let cgso = mos.model.cgso * mos.w;
        let cgdo = mos.model.cgdo * mos.w;
        let l_eff = (mos.l - 2.0 * mos.model.ld).max(1e-12);
        let cgbo = mos.model.cgbo * l_eff;

        if cgso > 0.0 {
            stamp_imag_conductance(&mut sys.imag, g, sp, omega * m * cgso);
        }
        if cgdo > 0.0 {
            stamp_imag_conductance(&mut sys.imag, g, dp, omega * m * cgdo);
        }
        if cgbo > 0.0 {
            stamp_imag_conductance(&mut sys.imag, g, b, omega * m * cgbo);
        }

        // Bulk junction capacitances (imaginary)
        if mos.model.cbd > 0.0 {
            stamp_imag_conductance(&mut sys.imag, b, dp, omega * m * mos.model.cbd);
        }
        if mos.model.cbs > 0.0 {
            stamp_imag_conductance(&mut sys.imag, b, sp, omega * m * mos.model.cbs);
        }
    }

    // 6. Apply AC source excitation to RHS.
    for element in netlist.elements() {
        match &element.kind {
            ferrospice_netlist::ElementKind::VoltageSource { source, .. } => {
                if let Some(ac_spec) = &source.ac {
                    let mag = expr_val_or(&ac_spec.mag, 0.0);
                    let phase_deg = ac_spec
                        .phase
                        .as_ref()
                        .map(|e| expr_val_or(e, 0.0))
                        .unwrap_or(0.0);
                    let phase_rad = phase_deg * PI / 180.0;
                    let ac_real = mag * phase_rad.cos();
                    let ac_imag = mag * phase_rad.sin();

                    // Find the branch index for this voltage source.
                    if let Some(branch_pos) = mna
                        .vsource_names
                        .iter()
                        .position(|n| n.to_lowercase() == element.name.to_lowercase())
                    {
                        let branch_idx = num_nodes + branch_pos;
                        sys.rhs_real[branch_idx] += ac_real;
                        sys.rhs_imag[branch_idx] += ac_imag;
                    }
                }
            }
            ferrospice_netlist::ElementKind::CurrentSource { pos, neg, source } => {
                if let Some(ac_spec) = &source.ac {
                    let mag = expr_val_or(&ac_spec.mag, 0.0);
                    let phase_deg = ac_spec
                        .phase
                        .as_ref()
                        .map(|e| expr_val_or(e, 0.0))
                        .unwrap_or(0.0);
                    let phase_rad = phase_deg * PI / 180.0;
                    let ac_real = mag * phase_rad.cos();
                    let ac_imag = mag * phase_rad.sin();

                    let ni = mna.node_map.get(pos);
                    let nj = mna.node_map.get(neg);

                    // Current source stamps on RHS (same convention as DC).
                    if let Some(i) = ni {
                        sys.rhs_real[i] -= ac_real;
                        sys.rhs_imag[i] -= ac_imag;
                    }
                    if let Some(j) = nj {
                        sys.rhs_real[j] += ac_real;
                        sys.rhs_imag[j] += ac_imag;
                    }
                }
            }
            _ => {}
        }
    }

    sys.solve().map_err(MnaError::SolveError)
}

/// Stamp an imaginary-part conductance (susceptance) between two nodes.
fn stamp_imag_conductance(
    matrix: &mut crate::SparseMatrix,
    ni: Option<usize>,
    nj: Option<usize>,
    b: f64,
) {
    if let Some(i) = ni {
        matrix.add(i, i, b);
    }
    if let Some(j) = nj {
        matrix.add(j, j, b);
    }
    if let (Some(i), Some(j)) = (ni, nj) {
        matrix.add(i, j, -b);
        matrix.add(j, i, -b);
    }
}

/// Generate frequency sweep points for AC analysis.
fn generate_ac_sweep(variation: AcVariation, n: u32, fstart: f64, fstop: f64) -> Vec<f64> {
    let mut freqs = Vec::new();

    match variation {
        AcVariation::Dec => {
            // n points per decade, logarithmic spacing.
            let log_ratio = (fstop / fstart).ln();
            let total_points = (n as f64 * log_ratio / 10.0_f64.ln()).round() as usize + 1;
            let delta = (10.0_f64.ln()) / n as f64;
            let mut f = fstart;
            for _ in 0..total_points {
                freqs.push(f);
                f *= delta.exp();
            }
        }
        AcVariation::Oct => {
            // n points per octave, logarithmic spacing.
            let log_ratio = (fstop / fstart).ln();
            let total_points = (n as f64 * log_ratio / 2.0_f64.ln()).round() as usize + 1;
            let delta = (2.0_f64.ln()) / n as f64;
            let mut f = fstart;
            for _ in 0..total_points {
                freqs.push(f);
                f *= delta.exp();
            }
        }
        AcVariation::Lin => {
            // n total points, linearly spaced.
            if n <= 1 {
                freqs.push(fstart);
            } else {
                let step = (fstop - fstart) / (n - 1) as f64;
                for i in 0..n {
                    freqs.push(fstart + i as f64 * step);
                }
            }
        }
    }

    freqs
}

/// Extract DC operating point solution values in matrix-index order.
fn extract_op_solution(op_result: &SimResult, mna: &MnaSystem) -> Vec<f64> {
    let num_ext_nodes = mna.node_map.len();
    let num_internal = mna
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
            .sum::<usize>();
    let num_nodes = num_ext_nodes + num_internal;
    let dim = num_nodes + mna.vsource_names.len();
    let mut solution = vec![0.0; dim];

    if let Some(plot) = op_result.plots.first() {
        for (name, idx) in mna.node_map.iter() {
            let vec_name = format!("v({})", name);
            if let Some(v) = plot.vecs.iter().find(|v| v.name == vec_name)
                && let Some(&val) = v.real.first()
            {
                solution[idx] = val;
            }
        }

        for (i, vsrc) in mna.vsource_names.iter().enumerate() {
            let vec_name = format!("{}#branch", vsrc.to_lowercase());
            if let Some(v) = plot.vecs.iter().find(|v| v.name == vec_name)
                && let Some(&val) = v.real.first()
            {
                solution[num_nodes + i] = val;
            }
        }
    }

    solution
}

fn expr_val(expr: &Expr, context: &str) -> Result<f64, MnaError> {
    match expr {
        Expr::Num(v) => Ok(*v),
        _ => Err(MnaError::NonNumericValue {
            element: context.to_string(),
        }),
    }
}

fn expr_val_or(expr: &Expr, default: f64) -> f64 {
    match expr {
        Expr::Num(v) => *v,
        _ => default,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_ac_sweep_dec() {
        let freqs = generate_ac_sweep(AcVariation::Dec, 10, 1e3, 1e6);
        // 3 decades × 10 pts/decade + 1 = 31 points
        assert_eq!(freqs.len(), 31);
        assert_abs_diff_eq!(freqs[0], 1e3, epsilon = 1.0);
        assert_abs_diff_eq!(freqs[10], 1e4, epsilon = 1.0);
        assert_abs_diff_eq!(freqs[20], 1e5, epsilon = 10.0);
        assert_abs_diff_eq!(freqs[30], 1e6, epsilon = 100.0);
    }

    #[test]
    fn test_ac_sweep_oct() {
        let freqs = generate_ac_sweep(AcVariation::Oct, 1, 100.0, 800.0);
        // 3 octaves (100→200→400→800) × 1 pt/oct + 1 = 4
        assert_eq!(freqs.len(), 4);
        assert_abs_diff_eq!(freqs[0], 100.0, epsilon = 0.1);
        assert_abs_diff_eq!(freqs[1], 200.0, epsilon = 0.1);
        assert_abs_diff_eq!(freqs[2], 400.0, epsilon = 0.1);
        assert_abs_diff_eq!(freqs[3], 800.0, epsilon = 0.1);
    }

    #[test]
    fn test_ac_sweep_lin() {
        let freqs = generate_ac_sweep(AcVariation::Lin, 5, 0.0, 100.0);
        assert_eq!(freqs.len(), 5);
        assert_abs_diff_eq!(freqs[0], 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(freqs[1], 25.0, epsilon = 1e-10);
        assert_abs_diff_eq!(freqs[4], 100.0, epsilon = 1e-10);
    }

    #[test]
    fn test_rc_lowpass_3db() {
        // RC lowpass: R=1k, C=1uF, f_3dB = 1/(2πRC) ≈ 159.15 Hz
        // At f_3dB, magnitude should be 1/√2 ≈ 0.7071
        let netlist = Netlist::parse(
            "RC lowpass filter
V1 1 0 DC 0 AC 1
R1 1 2 1k
C1 2 0 1u
.ac DEC 100 1 100k
.end
",
        )
        .unwrap();

        let result = simulate_ac(&netlist).unwrap();
        assert_eq!(result.plots.len(), 1);
        assert_eq!(result.plots[0].name, "ac1");

        let plot = &result.plots[0];
        let freq_vec = plot.vecs.iter().find(|v| v.name == "frequency").unwrap();
        let v2_vec = plot.vecs.iter().find(|v| v.name == "v(2)").unwrap();

        assert!(!freq_vec.real.is_empty());
        assert!(!v2_vec.complex.is_empty());

        // Find the frequency closest to f_3dB = 1/(2π × 1000 × 1e-6) ≈ 159.15 Hz
        let r = 1000.0;
        let c = 1e-6;
        let f_3db = 1.0 / (2.0 * PI * r * c);

        // Find the -3dB point by checking magnitude
        let mut found_3db = false;
        for (i, &freq) in freq_vec.real.iter().enumerate() {
            let re = v2_vec.complex[i].re;
            let im = v2_vec.complex[i].im;
            let mag = (re * re + im * im).sqrt();

            if (freq - f_3db).abs() / f_3db < 0.05 {
                // Within 5% of f_3dB
                let expected_mag = 1.0 / 2.0_f64.sqrt(); // 0.7071
                assert!(
                    (mag - expected_mag).abs() < 0.02,
                    "at f={freq:.1}Hz: mag={mag:.4}, expected ~{expected_mag:.4}"
                );
                found_3db = true;
            }
        }
        assert!(found_3db, "did not find frequency near f_3dB={f_3db:.1}Hz");

        // At low frequencies, magnitude should be close to 1.0
        let low_freq_mag = {
            let re = v2_vec.complex[0].re;
            let im = v2_vec.complex[0].im;
            (re * re + im * im).sqrt()
        };
        assert!(
            low_freq_mag > 0.99,
            "low freq magnitude should be ~1.0, got {low_freq_mag}"
        );

        // At high frequencies, magnitude should be much less than 1.0
        let last = v2_vec.complex.len() - 1;
        let high_freq_mag = {
            let re = v2_vec.complex[last].re;
            let im = v2_vec.complex[last].im;
            (re * re + im * im).sqrt()
        };
        assert!(
            high_freq_mag < 0.01,
            "high freq magnitude should be <<1.0, got {high_freq_mag}"
        );
    }

    #[test]
    fn test_ac_rl_highpass() {
        // RL highpass: R=1k in series, L=1H to ground → V(out) across R
        // Actually simpler: V1 -> R1 -> node 2, L1 from node 2 to ground
        // V(2)/V(1) = jωL / (R + jωL) = jωL/(R+jωL)
        // |H(f)| = ωL / √(R² + (ωL)²)
        // f_3dB = R/(2πL) = 1000/(2π×0.1) ≈ 1591.5 Hz
        let netlist = Netlist::parse(
            "RL highpass
V1 1 0 DC 0 AC 1
R1 1 2 1k
L1 2 0 0.1
.ac DEC 10 100 100k
.end
",
        )
        .unwrap();

        let result = simulate_ac(&netlist).unwrap();
        let plot = &result.plots[0];
        let v2_vec = plot.vecs.iter().find(|v| v.name == "v(2)").unwrap();

        // V(2) is across the inductor:
        // V(2) = V1 * jωL / (R + jωL) → as ω→∞, V(2) → V1
        // As ω→0, V(2) → 0 (inductor is short)

        // At low frequencies
        let low_mag = {
            let re = v2_vec.complex[0].re;
            let im = v2_vec.complex[0].im;
            (re * re + im * im).sqrt()
        };
        assert!(
            low_mag < 0.1,
            "low freq magnitude should be small, got {low_mag}"
        );

        // At high frequencies
        let last = v2_vec.complex.len() - 1;
        let high_mag = {
            let re = v2_vec.complex[last].re;
            let im = v2_vec.complex[last].im;
            (re * re + im * im).sqrt()
        };
        assert!(
            high_mag > 0.9,
            "high freq magnitude should be ~1.0, got {high_mag}"
        );
    }
}
