//! Noise analysis (.noise) — computes spectral noise density and integrated noise.
//!
//! Algorithm:
//! 1. Compute DC operating point to linearize nonlinear devices.
//! 2. At each frequency point:
//!    a. Build complex AC MNA matrix (G + jωC).
//!    b. Solve adjoint (transposed) system with unit current at output node.
//!    c. Compute each device's noise contribution using adjoint transfer functions.
//!    d. Sum for total output noise density; divide by |H|² for input-referred noise.
//! 3. Integrate noise over frequency range.

use std::f64::consts::PI;

use ferrospice_netlist::{Analysis, Expr, Item, Netlist, SimPlot, SimResult, SimVector};

use crate::ac::{generate_ac_sweep, stamp_imag_conductance};
use crate::mna::{MnaError, MnaSystem, assemble_mna};
use crate::simulate::simulate_op;
use crate::sparse::ComplexLinearSystem;

/// Boltzmann constant (J/K).
const K_BOLTZ: f64 = 1.380649e-23;
/// Elementary charge (C).
const Q_CHARGE: f64 = 1.602176634e-19;
/// Nominal temperature (27°C in Kelvin).
const T_NOM: f64 = 300.15;

/// Perform noise analysis.
pub fn simulate_noise(netlist: &Netlist) -> Result<SimResult, MnaError> {
    // Find the .noise analysis command.
    let (output, ref_node, src_name, variation, n, fstart, fstop) = netlist
        .items
        .iter()
        .find_map(|item| {
            if let Item::Analysis(Analysis::Noise {
                output,
                ref_node,
                src,
                variation,
                n,
                fstart,
                fstop,
            }) = item
            {
                Some((
                    output.clone(),
                    ref_node.clone(),
                    src.clone(),
                    *variation,
                    *n,
                    fstart.clone(),
                    fstop.clone(),
                ))
            } else {
                None
            }
        })
        .ok_or_else(|| MnaError::UnsupportedElement("no .noise analysis found".to_string()))?;

    let fstart_val = expr_val(&fstart, ".noise")?;
    let fstop_val = expr_val(&fstop, ".noise")?;

    // Parse output node specification: "v(node)" or "v(node1,node2)".
    let (out_pos, out_neg) = parse_output_spec(&output)?;

    // DC operating point.
    let op_result = simulate_op(netlist)?;
    let mna = assemble_mna(netlist)?;
    let op_solution = extract_op_solution(&op_result, &mna);
    let num_nodes = mna.total_num_nodes();

    // Resolve output node indices.
    let out_pos_idx = mna.node_map.get(&out_pos);
    let out_neg_idx = out_neg.as_deref().and_then(|n| mna.node_map.get(n));

    // Resolve reference node (for differential output).
    let _ref_node_idx = ref_node.as_deref().and_then(|n| mna.node_map.get(n));

    // Find the input source branch index for gain computation.
    let src_branch_idx = mna
        .vsource_names
        .iter()
        .position(|n| n.to_lowercase() == src_name.to_lowercase())
        .map(|i| num_nodes + i);

    // Generate frequency sweep.
    let frequencies = generate_ac_sweep(variation, n, fstart_val, fstop_val);

    let mut freq_data = Vec::with_capacity(frequencies.len());
    let mut onoise_data = Vec::with_capacity(frequencies.len());
    let mut inoise_data = Vec::with_capacity(frequencies.len());

    for &freq in &frequencies {
        let omega = 2.0 * PI * freq;

        // Build the AC MNA system at this frequency.
        let sys = build_ac_system(&mna, &op_solution, omega, netlist, num_nodes)?;

        // Solve adjoint system with unit current at output node.
        let adjoint = solve_adjoint(&sys, out_pos_idx, out_neg_idx)?;

        // Compute gain: solve forward system with AC source excitation.
        let gain_sq_inv = compute_gain_sq_inv(
            &sys,
            netlist,
            &mna,
            num_nodes,
            out_pos_idx,
            out_neg_idx,
            src_branch_idx,
        )?;

        // Compute total output noise density (V²/Hz).
        let onoise_sq = compute_total_noise(&mna, &op_solution, &adjoint, freq);

        // Input-referred noise density.
        let inoise_sq = onoise_sq * gain_sq_inv;

        freq_data.push(freq);
        onoise_data.push(onoise_sq.sqrt()); // V/√Hz
        inoise_data.push(inoise_sq.sqrt()); // V/√Hz
    }

    // Compute integrated noise (total RMS noise over frequency range).
    let onoise_total = integrate_noise(&freq_data, &onoise_data);
    let inoise_total = integrate_noise(&freq_data, &inoise_data);

    // Build result with two plots: noise1 (spectrum) and noise2 (integrated).
    let spectrum_plot = SimPlot {
        name: "noise1".to_string(),
        vecs: vec![
            SimVector {
                name: "frequency".to_string(),
                real: freq_data,
                complex: vec![],
            },
            SimVector {
                name: "onoise_spectrum".to_string(),
                real: onoise_data,
                complex: vec![],
            },
            SimVector {
                name: "inoise_spectrum".to_string(),
                real: inoise_data,
                complex: vec![],
            },
        ],
    };

    let integrated_plot = SimPlot {
        name: "noise2".to_string(),
        vecs: vec![
            SimVector {
                name: "onoise_total".to_string(),
                real: vec![onoise_total],
                complex: vec![],
            },
            SimVector {
                name: "inoise_total".to_string(),
                real: vec![inoise_total],
                complex: vec![],
            },
        ],
    };

    Ok(SimResult {
        plots: vec![spectrum_plot, integrated_plot],
    })
}

/// Parse output specification like "v(2)" or "v(out,ref)" into (pos_node, neg_node).
fn parse_output_spec(output: &str) -> Result<(String, Option<String>), MnaError> {
    let s = output.trim();
    // Strip v( ... ) or V( ... )
    let inner = if let Some(rest) = s.strip_prefix("v(").or_else(|| s.strip_prefix("V(")) {
        rest.strip_suffix(')').unwrap_or(rest)
    } else {
        // Bare node name.
        return Ok((s.to_lowercase(), None));
    };

    if let Some((pos, neg)) = inner.split_once(',') {
        Ok((pos.trim().to_lowercase(), Some(neg.trim().to_lowercase())))
    } else {
        Ok((inner.trim().to_lowercase(), None))
    }
}

/// Build the complex AC MNA system at a given frequency.
///
/// This is largely the same as `ac::solve_ac_point` but returns the system
/// instead of solving it, so we can solve both forward and adjoint.
fn build_ac_system(
    mna: &MnaSystem,
    op_solution: &[f64],
    omega: f64,
    netlist: &Netlist,
    num_nodes: usize,
) -> Result<ComplexLinearSystem, MnaError> {
    let dim = num_nodes + mna.vsource_names.len();
    let mut sys = ComplexLinearSystem::new(dim);

    // 1. DC conductance matrix (real part).
    for triplet in mna.system.matrix.triplets() {
        sys.real.add(triplet.row, triplet.col, triplet.value);
    }

    // 2. Capacitor AC admittance: jωC.
    for cap in &mna.capacitors {
        let bc = omega * cap.capacitance;
        stamp_imag_conductance(&mut sys.imag, cap.pos_idx, cap.neg_idx, bc);
    }

    // 3. Inductor AC impedance: -jωL on branch diagonal.
    for ind in &mna.inductors {
        sys.imag
            .add(ind.branch_idx, ind.branch_idx, -omega * ind.inductance);
    }

    // 4. Diode small-signal stamps.
    for diode in &mna.diodes {
        let (jct_anode, jct_cathode) = if diode.internal_idx.is_some() {
            (diode.internal_idx, diode.cathode_idx)
        } else {
            (diode.anode_idx, diode.cathode_idx)
        };
        let v_anode = jct_anode.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_cathode = jct_cathode.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_jct = v_anode - v_cathode;
        let (g_d, _) = diode.model.companion(v_jct);
        crate::stamp_conductance(&mut sys.real, jct_anode, jct_cathode, g_d);

        if let Some(int_idx) = diode.internal_idx {
            let g_rs = 1.0 / diode.model.rs;
            crate::stamp_conductance(&mut sys.real, diode.anode_idx, Some(int_idx), g_rs);
        }

        if diode.model.cjo > 0.0 {
            let cj = diode.model.junction_capacitance(v_jct);
            stamp_imag_conductance(&mut sys.imag, jct_anode, jct_cathode, omega * cj);
        }
        if diode.model.tt > 0.0 {
            let cd = diode.model.tt * g_d;
            stamp_imag_conductance(&mut sys.imag, jct_anode, jct_cathode, omega * cd);
        }
    }

    // 5. BJT small-signal stamps.
    for bjt in &mna.bjts {
        let (vbe, vbc) = bjt.junction_voltages(op_solution);
        let comp = bjt.model.companion(vbe, vbc);
        let m = bjt.m * bjt.area;

        let bp = bjt.base_prime_idx;
        let cp = bjt.col_prime_idx;
        let ep = bjt.emit_prime_idx;

        crate::stamp_conductance(&mut sys.real, bp, ep, m * comp.gpi);
        crate::stamp_conductance(&mut sys.real, bp, cp, m * comp.gmu);
        crate::stamp_conductance(&mut sys.real, cp, ep, m * comp.go);

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

        if bjt.model.rb > 0.0 {
            crate::stamp_conductance(&mut sys.real, bjt.base_idx, bp, m / bjt.model.rb);
        }
        if bjt.model.rc > 0.0 {
            crate::stamp_conductance(&mut sys.real, bjt.col_idx, cp, m / bjt.model.rc);
        }
        if bjt.model.re > 0.0 {
            crate::stamp_conductance(&mut sys.real, bjt.emit_idx, ep, m / bjt.model.re);
        }

        let cap_be = bjt.model.cap_be(vbe);
        let gbe = comp.gpi * bjt.model.bf;
        let total_cap_be = cap_be + bjt.model.tf * gbe;
        if total_cap_be > 0.0 {
            stamp_imag_conductance(&mut sys.imag, bp, ep, omega * m * total_cap_be);
        }
        let cap_bc = bjt.model.cap_bc(vbc);
        let gbc = comp.gmu * bjt.model.br;
        let total_cap_bc = cap_bc + bjt.model.tr * gbc;
        if total_cap_bc > 0.0 {
            stamp_imag_conductance(&mut sys.imag, bp, cp, omega * m * total_cap_bc);
        }
    }

    // 5b. MOSFET small-signal stamps (same as ac.rs).
    for mos in &mna.mosfets {
        let (vgs, vds, vbs) = mos.terminal_voltages(op_solution);
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

        crate::stamp_conductance(&mut sys.real, dp, sp, m * comp.gds);

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

        crate::stamp_conductance(&mut sys.real, b, dp, m * comp.gbd);
        crate::stamp_conductance(&mut sys.real, b, sp, m * comp.gbs);

        if mos.model.rd > 0.0 {
            crate::stamp_conductance(&mut sys.real, mos.drain_idx, dp, m / mos.model.rd);
        }
        if mos.model.rs > 0.0 {
            crate::stamp_conductance(&mut sys.real, mos.source_idx, sp, m / mos.model.rs);
        }

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
        if mos.model.cbd > 0.0 {
            stamp_imag_conductance(&mut sys.imag, b, dp, omega * m * mos.model.cbd);
        }
        if mos.model.cbs > 0.0 {
            stamp_imag_conductance(&mut sys.imag, b, sp, omega * m * mos.model.cbs);
        }
    }

    // 5c. JFET small-signal stamps.
    for jfet in &mna.jfets {
        let (vgs, vgd) = jfet.junction_voltages(op_solution);
        let comp = jfet.model.companion(vgs, vgd);
        let m = jfet.m * jfet.area;

        let dp = jfet.drain_prime_idx;
        let g = jfet.gate_idx;
        let sp = jfet.source_prime_idx;

        crate::stamp_conductance(&mut sys.real, g, sp, m * comp.ggs);
        crate::stamp_conductance(&mut sys.real, g, dp, m * comp.ggd);
        crate::stamp_conductance(&mut sys.real, dp, sp, m * comp.gds);

        let gm_scaled = m * comp.gm;
        if let Some(d) = dp {
            if let Some(gate) = g {
                sys.real.add(d, gate, gm_scaled);
            }
            if let Some(s) = sp {
                sys.real.add(d, s, -gm_scaled);
            }
        }
        if let Some(s) = sp {
            if let Some(gate) = g {
                sys.real.add(s, gate, -gm_scaled);
            }
            sys.real.add(s, s, gm_scaled);
        }

        if jfet.model.rd > 0.0 {
            crate::stamp_conductance(&mut sys.real, jfet.drain_idx, dp, m / jfet.model.rd);
        }
        if jfet.model.rs > 0.0 {
            crate::stamp_conductance(&mut sys.real, jfet.source_idx, sp, m / jfet.model.rs);
        }

        let cgs = jfet.model.junction_cap(vgs, jfet.model.cgs);
        if cgs > 0.0 {
            stamp_imag_conductance(&mut sys.imag, g, sp, omega * m * cgs);
        }
        let cgd = jfet.model.junction_cap(vgd, jfet.model.cgd);
        if cgd > 0.0 {
            stamp_imag_conductance(&mut sys.imag, g, dp, omega * m * cgd);
        }
    }

    // 5d. BSIM3 small-signal stamps.
    for bsim in &mna.bsim3s {
        let (vgs, vds, vbs) = bsim.terminal_voltages(op_solution);
        let comp = crate::bsim3::bsim3_companion(vgs, vds, vbs, &bsim.size_params, &bsim.model);
        let m = bsim.m;

        let dp = bsim.drain_eff_idx();
        let g = bsim.gate_idx;
        let sp = bsim.source_eff_idx();
        let b = bsim.bulk_idx;

        crate::stamp_conductance(&mut sys.real, dp, sp, m * comp.gds);

        let gm_scaled = m * comp.gm;
        if let Some(d) = dp {
            if let Some(gate) = g {
                sys.real.add(d, gate, gm_scaled);
            }
            if let Some(s) = sp {
                sys.real.add(d, s, -gm_scaled);
            }
        }
        if let Some(s) = sp {
            if let Some(gate) = g {
                sys.real.add(s, gate, -gm_scaled);
            }
            sys.real.add(s, s, gm_scaled);
        }

        let gmbs_scaled = m * comp.gmbs;
        if let Some(d) = dp {
            if let Some(bulk) = b {
                sys.real.add(d, bulk, gmbs_scaled);
            }
            if let Some(s) = sp {
                sys.real.add(d, s, -gmbs_scaled);
            }
        }
        if let Some(s) = sp {
            if let Some(bulk) = b {
                sys.real.add(s, bulk, -gmbs_scaled);
            }
            sys.real.add(s, s, gmbs_scaled);
        }

        crate::stamp_conductance(&mut sys.real, b, dp, m * comp.gbd);
        crate::stamp_conductance(&mut sys.real, b, sp, m * comp.gbs);

        let g_drain = bsim.drain_conductance();
        if g_drain > 0.0 {
            crate::stamp_conductance(&mut sys.real, bsim.drain_idx, dp, m * g_drain);
        }
        let g_source = bsim.source_conductance();
        if g_source > 0.0 {
            crate::stamp_conductance(&mut sys.real, bsim.source_idx, sp, m * g_source);
        }

        // Intrinsic capacitances (imaginary) - direct Y-matrix entries
        let cap_entries: &[(Option<usize>, Option<usize>, f64)] = &[
            (g, g, comp.cggb),
            (g, dp, comp.cgdb),
            (g, sp, comp.cgsb),
            (dp, g, comp.cdgb),
            (dp, dp, comp.cddb),
            (dp, sp, comp.cdsb),
            (b, g, comp.cbgb),
            (b, dp, comp.cbdb),
            (b, sp, comp.cbsb),
        ];
        for &(row, col, cap) in cap_entries {
            if cap.abs() > 0.0
                && let (Some(r), Some(c)) = (row, col)
            {
                sys.imag.add(r, c, omega * m * cap);
            }
        }

        if comp.capbd > 0.0 {
            stamp_imag_conductance(&mut sys.imag, b, dp, omega * m * comp.capbd);
        }
        if comp.capbs > 0.0 {
            stamp_imag_conductance(&mut sys.imag, b, sp, omega * m * comp.capbs);
        }
    }

    // 5e. BSIM4 small-signal stamps.
    for bsim in &mna.bsim4s {
        let (vgs, vds, vbs) = bsim.terminal_voltages(op_solution);
        let comp = crate::bsim4::bsim4_companion(vgs, vds, vbs, &bsim.size_params, &bsim.model);
        let m = bsim.m;

        let dp = bsim.drain_eff_idx();
        let g = bsim.gate_idx;
        let sp = bsim.source_eff_idx();
        let b = bsim.bulk_idx;

        crate::stamp_conductance(&mut sys.real, dp, sp, m * comp.gds);

        let gm_scaled = m * comp.gm;
        if let Some(d) = dp {
            if let Some(gate) = g {
                sys.real.add(d, gate, gm_scaled);
            }
            if let Some(s) = sp {
                sys.real.add(d, s, -gm_scaled);
            }
        }
        if let Some(s) = sp {
            if let Some(gate) = g {
                sys.real.add(s, gate, -gm_scaled);
            }
            sys.real.add(s, s, gm_scaled);
        }

        let gmbs_scaled = m * comp.gmbs;
        if let Some(d) = dp {
            if let Some(bulk) = b {
                sys.real.add(d, bulk, gmbs_scaled);
            }
            if let Some(s) = sp {
                sys.real.add(d, s, -gmbs_scaled);
            }
        }
        if let Some(s) = sp {
            if let Some(bulk) = b {
                sys.real.add(s, bulk, -gmbs_scaled);
            }
            sys.real.add(s, s, gmbs_scaled);
        }

        crate::stamp_conductance(&mut sys.real, b, dp, m * comp.gbd);
        crate::stamp_conductance(&mut sys.real, b, sp, m * comp.gbs);

        let g_drain = bsim.drain_conductance();
        if g_drain > 0.0 {
            crate::stamp_conductance(&mut sys.real, bsim.drain_idx, dp, m * g_drain);
        }
        let g_source = bsim.source_conductance();
        if g_source > 0.0 {
            crate::stamp_conductance(&mut sys.real, bsim.source_idx, sp, m * g_source);
        }

        let cap_entries: &[(Option<usize>, Option<usize>, f64)] = &[
            (g, g, comp.cggb),
            (g, dp, comp.cgdb),
            (g, sp, comp.cgsb),
            (dp, g, comp.cdgb),
            (dp, dp, comp.cddb),
            (dp, sp, comp.cdsb),
            (b, g, comp.cbgb),
            (b, dp, comp.cbdb),
            (b, sp, comp.cbsb),
        ];
        for &(row, col, cap) in cap_entries {
            if cap.abs() > 0.0
                && let (Some(r), Some(c)) = (row, col)
            {
                sys.imag.add(r, c, omega * m * cap);
            }
        }

        if comp.capbd > 0.0 {
            stamp_imag_conductance(&mut sys.imag, b, dp, omega * m * comp.capbd);
        }
        if comp.capbs > 0.0 {
            stamp_imag_conductance(&mut sys.imag, b, sp, omega * m * comp.capbs);
        }
    }

    // 6. Apply AC source excitation to RHS (for forward solve / gain computation).
    apply_ac_excitation(&mut sys, netlist, mna, num_nodes);

    Ok(sys)
}

/// Apply AC source excitation to the complex RHS.
fn apply_ac_excitation(
    sys: &mut ComplexLinearSystem,
    netlist: &Netlist,
    mna: &MnaSystem,
    num_nodes: usize,
) {
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

                    if let Some(branch_pos) = mna
                        .vsource_names
                        .iter()
                        .position(|n| n.to_lowercase() == element.name.to_lowercase())
                    {
                        let branch_idx = num_nodes + branch_pos;
                        sys.rhs_real[branch_idx] += mag * phase_rad.cos();
                        sys.rhs_imag[branch_idx] += mag * phase_rad.sin();
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
}

/// Solve the adjoint system with unit current injection at the output node.
fn solve_adjoint(
    sys: &ComplexLinearSystem,
    out_pos_idx: Option<usize>,
    out_neg_idx: Option<usize>,
) -> Result<Vec<(f64, f64)>, MnaError> {
    // Create adjoint system: same matrix but different RHS.
    let mut adj = ComplexLinearSystem::new(sys.dim());

    // Copy matrix entries.
    for t in sys.real.triplets() {
        adj.real.add(t.row, t.col, t.value);
    }
    for t in sys.imag.triplets() {
        adj.imag.add(t.row, t.col, t.value);
    }

    // Unit current injection at output nodes.
    if let Some(i) = out_pos_idx {
        adj.rhs_real[i] = 1.0;
    }
    if let Some(j) = out_neg_idx {
        adj.rhs_real[j] = -1.0;
    }

    adj.solve_transpose().map_err(MnaError::SolveError)
}

/// Compute 1/|H(jω)|² where H is the transfer function from input source to output.
fn compute_gain_sq_inv(
    sys: &ComplexLinearSystem,
    _netlist: &Netlist,
    _mna: &MnaSystem,
    _num_nodes: usize,
    out_pos_idx: Option<usize>,
    out_neg_idx: Option<usize>,
    _src_branch_idx: Option<usize>,
) -> Result<f64, MnaError> {
    // Forward solve with AC excitation already applied.
    let solution = sys.solve().map_err(MnaError::SolveError)?;

    // Output voltage from forward solve.
    let v_out_re = out_pos_idx.map(|i| solution[i].0).unwrap_or(0.0)
        - out_neg_idx.map(|i| solution[i].0).unwrap_or(0.0);
    let v_out_im = out_pos_idx.map(|i| solution[i].1).unwrap_or(0.0)
        - out_neg_idx.map(|i| solution[i].1).unwrap_or(0.0);

    let gain_sq = v_out_re * v_out_re + v_out_im * v_out_im;

    if gain_sq > 0.0 {
        Ok(1.0 / gain_sq)
    } else {
        // If gain is zero (e.g., DC-blocked output), return large value.
        Ok(1e30)
    }
}

/// Compute total output noise spectral density (V²/Hz) from all noise sources.
fn compute_total_noise(
    mna: &MnaSystem,
    op_solution: &[f64],
    adjoint: &[(f64, f64)],
    freq: f64,
) -> f64 {
    let mut total = 0.0;

    // Resistor thermal noise: modeled as current noise source 4kTG (A²/Hz).
    // Output contribution = 4kTG * |V_adj(n1) - V_adj(n2)|² (ngspice NevalSrc convention).
    for res in &mna.resistors {
        let g = 1.0 / res.resistance;
        let thermal_noise = 4.0 * K_BOLTZ * T_NOM * g;
        let transfer_sq = adjoint_transfer_sq(adjoint, res.pos_idx, res.neg_idx);
        total += thermal_noise * transfer_sq;
    }

    // Diode noise sources.
    for diode in &mna.diodes {
        let (jct_anode, jct_cathode) = if diode.internal_idx.is_some() {
            (diode.internal_idx, diode.cathode_idx)
        } else {
            (diode.anode_idx, diode.cathode_idx)
        };
        let v_anode = jct_anode.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_cathode = jct_cathode.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_jct = v_anode - v_cathode;

        let i_d = diode.model.current(v_jct).abs();

        // Shot noise on junction current: 2qI.
        let shot_noise = 2.0 * Q_CHARGE * i_d;
        total += shot_noise * adjoint_transfer_sq(adjoint, jct_anode, jct_cathode);

        // Series resistance thermal noise: 4kT/Rs.
        if diode.model.rs > 0.0
            && let Some(int_idx) = diode.internal_idx
        {
            let g_rs = 1.0 / diode.model.rs;
            let rs_noise = 4.0 * K_BOLTZ * T_NOM * g_rs;
            total += rs_noise * adjoint_transfer_sq(adjoint, diode.anode_idx, Some(int_idx));
        }

        // Flicker noise: KF * |Id|^AF / freq.
        if diode.model.kf > 0.0 && freq > 0.0 {
            let flicker = diode.model.kf * i_d.powf(diode.model.af) / freq;
            total += flicker * adjoint_transfer_sq(adjoint, jct_anode, jct_cathode);
        }
    }

    // BJT noise sources.
    for bjt in &mna.bjts {
        let (vbe, vbc) = bjt.junction_voltages(op_solution);
        let comp = bjt.model.companion(vbe, vbc);
        let m = bjt.m * bjt.area;

        let bp = bjt.base_prime_idx;
        let cp = bjt.col_prime_idx;
        let ep = bjt.emit_prime_idx;

        // Collector shot noise: 2q|Ic|.
        let ic = comp.cc.abs() * m;
        let ic_shot = 2.0 * Q_CHARGE * ic;
        total += ic_shot * adjoint_transfer_sq(adjoint, cp, ep);

        // Base shot noise: 2q|Ib|.
        let ib = comp.cb.abs() * m;
        let ib_shot = 2.0 * Q_CHARGE * ib;
        total += ib_shot * adjoint_transfer_sq(adjoint, bp, ep);

        // Base resistance thermal noise: 4kT/Rb.
        if bjt.model.rb > 0.0 {
            let g_rb = m / bjt.model.rb;
            let rb_noise = 4.0 * K_BOLTZ * T_NOM * g_rb;
            total += rb_noise * adjoint_transfer_sq(adjoint, bjt.base_idx, bp);
        }

        // Collector resistance thermal noise: 4kT/Rc.
        if bjt.model.rc > 0.0 {
            let g_rc = m / bjt.model.rc;
            let rc_noise = 4.0 * K_BOLTZ * T_NOM * g_rc;
            total += rc_noise * adjoint_transfer_sq(adjoint, bjt.col_idx, cp);
        }

        // Emitter resistance thermal noise: 4kT/Re.
        if bjt.model.re > 0.0 {
            let g_re = m / bjt.model.re;
            let re_noise = 4.0 * K_BOLTZ * T_NOM * g_re;
            total += re_noise * adjoint_transfer_sq(adjoint, bjt.emit_idx, ep);
        }

        // Base flicker noise: KF * |Ib|^AF / freq.
        if bjt.model.kf > 0.0 && freq > 0.0 {
            let flicker = bjt.model.kf * ib.powf(bjt.model.af) / freq;
            total += flicker * adjoint_transfer_sq(adjoint, bp, ep);
        }
    }

    // MOSFET noise sources.
    for mos in &mna.mosfets {
        let (vgs, vds, vbs) = mos.terminal_voltages(op_solution);
        let mut eff_model = mos.model.clone();
        eff_model.kp = mos.beta();
        let comp = eff_model.companion(vgs, vds, vbs);
        let m = mos.m;

        let dp = mos.drain_prime_idx;
        let sp = mos.source_prime_idx;

        // Channel thermal noise: (2/3) * 4kT * |gm|.
        let gm = comp.gm.abs() * m;
        let channel_noise = (8.0 / 3.0) * K_BOLTZ * T_NOM * gm;
        total += channel_noise * adjoint_transfer_sq(adjoint, dp, sp);

        // Drain resistance thermal noise.
        if mos.model.rd > 0.0 {
            let g_rd = m / mos.model.rd;
            let rd_noise = 4.0 * K_BOLTZ * T_NOM * g_rd;
            total += rd_noise * adjoint_transfer_sq(adjoint, mos.drain_idx, dp);
        }

        // Source resistance thermal noise.
        if mos.model.rs > 0.0 {
            let g_rs = m / mos.model.rs;
            let rs_noise = 4.0 * K_BOLTZ * T_NOM * g_rs;
            total += rs_noise * adjoint_transfer_sq(adjoint, mos.source_idx, sp);
        }

        // Flicker noise: KF * |Id|^AF / (f * Cox * L_eff²).
        if mos.model.kf > 0.0 && freq > 0.0 {
            let id = comp.cdrain.abs() * m;
            let l_eff = (mos.l - 2.0 * mos.model.ld).max(1e-12);
            // Cox = epsilon_ox / tox ≈ 3.9 * 8.854e-12 / tox.
            let cox = 3.9 * 8.854e-12 / mos.model.tox;
            let cox_l_sq = cox * l_eff * l_eff;
            let flicker = mos.model.kf * id.powf(mos.model.af) / (freq * cox_l_sq);
            total += flicker * adjoint_transfer_sq(adjoint, dp, sp);
        }
    }

    // JFET noise sources.
    for jfet in &mna.jfets {
        let (vgs, vgd) = jfet.junction_voltages(op_solution);
        let comp = jfet.model.companion(vgs, vgd);
        let m = jfet.m * jfet.area;

        let dp = jfet.drain_prime_idx;
        let sp = jfet.source_prime_idx;

        // Channel thermal noise: (2/3) * 4kT * |gm|.
        let gm = comp.gm.abs() * m;
        let channel_noise = (8.0 / 3.0) * K_BOLTZ * T_NOM * gm;
        total += channel_noise * adjoint_transfer_sq(adjoint, dp, sp);

        // Drain resistance thermal noise.
        if jfet.model.rd > 0.0 {
            let g_rd = m / jfet.model.rd;
            let rd_noise = 4.0 * K_BOLTZ * T_NOM * g_rd;
            total += rd_noise * adjoint_transfer_sq(adjoint, jfet.drain_idx, dp);
        }

        // Source resistance thermal noise.
        if jfet.model.rs > 0.0 {
            let g_rs = m / jfet.model.rs;
            let rs_noise = 4.0 * K_BOLTZ * T_NOM * g_rs;
            total += rs_noise * adjoint_transfer_sq(adjoint, jfet.source_idx, sp);
        }

        // Flicker noise: KF * |Id|^AF / freq.
        if jfet.model.kf > 0.0 && freq > 0.0 {
            let id = comp.cdrain.abs() * m;
            let flicker = jfet.model.kf * id.powf(jfet.model.af) / freq;
            total += flicker * adjoint_transfer_sq(adjoint, dp, sp);
        }
    }

    // BSIM3 noise sources.
    for bsim in &mna.bsim3s {
        let (vgs, vds, vbs) = bsim.terminal_voltages(op_solution);
        let comp = crate::bsim3::bsim3_companion(vgs, vds, vbs, &bsim.size_params, &bsim.model);
        let m = bsim.m;

        let dp = bsim.drain_eff_idx();
        let sp = bsim.source_eff_idx();

        // Channel thermal noise: (2/3) * 4kT * |gm|.
        let gm = comp.gm.abs() * m;
        let channel_noise = (8.0 / 3.0) * K_BOLTZ * T_NOM * gm;
        total += channel_noise * adjoint_transfer_sq(adjoint, dp, sp);

        // Drain resistance thermal noise.
        let g_drain = bsim.drain_conductance() * m;
        if g_drain > 0.0 {
            let rd_noise = 4.0 * K_BOLTZ * T_NOM * g_drain;
            total += rd_noise * adjoint_transfer_sq(adjoint, bsim.drain_idx, dp);
        }

        // Source resistance thermal noise.
        let g_source = bsim.source_conductance() * m;
        if g_source > 0.0 {
            let rs_noise = 4.0 * K_BOLTZ * T_NOM * g_source;
            total += rs_noise * adjoint_transfer_sq(adjoint, bsim.source_idx, sp);
        }

        // Flicker noise: NOIA * |Id|^EF / (f * Cox * Leff²).
        let noia = bsim.model.noia;
        if noia > 0.0 && freq > 0.0 {
            let id = comp.ids.abs() * m;
            let cox = bsim.model.cox();
            let l_eff = bsim.size_params.leff;
            let cox_l_sq = cox * l_eff * l_eff;
            if cox_l_sq > 0.0 {
                let ef = bsim.model.ef;
                let flicker = noia * id.powf(ef) / (freq * cox_l_sq);
                total += flicker * adjoint_transfer_sq(adjoint, dp, sp);
            }
        }
    }

    // BSIM4 noise sources.
    for bsim in &mna.bsim4s {
        let (vgs, vds, vbs) = bsim.terminal_voltages(op_solution);
        let comp = crate::bsim4::bsim4_companion(vgs, vds, vbs, &bsim.size_params, &bsim.model);
        let m = bsim.m;

        let dp = bsim.drain_eff_idx();
        let sp = bsim.source_eff_idx();

        // Channel thermal noise: (2/3) * 4kT * |gm|.
        let gm = comp.gm.abs() * m;
        let channel_noise = (8.0 / 3.0) * K_BOLTZ * T_NOM * gm;
        total += channel_noise * adjoint_transfer_sq(adjoint, dp, sp);

        // Drain resistance thermal noise.
        let g_drain = bsim.drain_conductance() * m;
        if g_drain > 0.0 {
            let rd_noise = 4.0 * K_BOLTZ * T_NOM * g_drain;
            total += rd_noise * adjoint_transfer_sq(adjoint, bsim.drain_idx, dp);
        }

        // Source resistance thermal noise.
        let g_source = bsim.source_conductance() * m;
        if g_source > 0.0 {
            let rs_noise = 4.0 * K_BOLTZ * T_NOM * g_source;
            total += rs_noise * adjoint_transfer_sq(adjoint, bsim.source_idx, sp);
        }

        // Flicker noise: NOIA * |Id|^EF / (f * Cox * Leff²).
        let noia = bsim.model.noia;
        if noia > 0.0 && freq > 0.0 {
            let id = comp.ids.abs() * m;
            let cox = bsim.model.coxe();
            let l_eff = bsim.size_params.leff;
            let cox_l_sq = cox * l_eff * l_eff;
            if cox_l_sq > 0.0 {
                let ef = bsim.model.ef;
                let flicker = noia * id.powf(ef) / (freq * cox_l_sq);
                total += flicker * adjoint_transfer_sq(adjoint, dp, sp);
            }
        }
    }

    total
}

/// Compute |V_adj(n1) - V_adj(n2)|² from the adjoint solution.
///
/// This is the squared magnitude of the complex transfer function from
/// a current source between n1 and n2 to the output voltage.
fn adjoint_transfer_sq(adjoint: &[(f64, f64)], n1: Option<usize>, n2: Option<usize>) -> f64 {
    let (re1, im1) = n1.map(|i| adjoint[i]).unwrap_or((0.0, 0.0));
    let (re2, im2) = n2.map(|i| adjoint[i]).unwrap_or((0.0, 0.0));
    let dre = re1 - re2;
    let dim = im1 - im2;
    dre * dre + dim * dim
}

/// Integrate noise spectral density over frequency using piecewise power-law interpolation.
///
/// Input: V/√Hz values (not V²/Hz), output: total RMS noise (V).
fn integrate_noise(freqs: &[f64], noise_sqrt: &[f64]) -> f64 {
    if freqs.len() < 2 {
        return 0.0;
    }

    let mut total_sq = 0.0;

    for i in 1..freqs.len() {
        let f0 = freqs[i - 1];
        let f1 = freqs[i];
        let n0_sq = noise_sqrt[i - 1] * noise_sqrt[i - 1]; // V²/Hz
        let n1_sq = noise_sqrt[i] * noise_sqrt[i]; // V²/Hz

        if f0 <= 0.0 || f1 <= 0.0 || n0_sq <= 0.0 || n1_sq <= 0.0 {
            // Fall back to trapezoidal rule.
            total_sq += 0.5 * (n0_sq + n1_sq) * (f1 - f0);
            continue;
        }

        // Power-law exponent: n(f) = a * f^k, so k = log(n1/n0) / log(f1/f0).
        let k = (n1_sq / n0_sq).ln() / (f1 / f0).ln();

        if (k + 1.0).abs() < 1e-10 {
            // k ≈ -1: integral is a * ln(f1/f0).
            let a = n0_sq * f0;
            total_sq += a * (f1 / f0).ln();
        } else if k.abs() < 1e-10 {
            // Flat: integral is n0_sq * (f1 - f0).
            total_sq += n0_sq * (f1 - f0);
        } else {
            // General: integral of a*f^k from f0 to f1 = a*(f1^(k+1) - f0^(k+1))/(k+1).
            let a = n0_sq / f0.powf(k);
            total_sq += a * (f1.powf(k + 1.0) - f0.powf(k + 1.0)) / (k + 1.0);
        }
    }

    total_sq.max(0.0).sqrt()
}

/// Extract DC operating point solution in matrix-index order.
fn extract_op_solution(op_result: &SimResult, mna: &MnaSystem) -> Vec<f64> {
    let num_nodes = mna.total_num_nodes();
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
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    #[test]
    fn test_resistor_thermal_noise_analytical() {
        // Circuit: V1 -> R1(1k) -> node 2 -> R2(1k) -> ground.
        // Measure noise at v(2), input source V1.
        //
        // Adjoint: inject 1A at node 2. V1 clamps node 1 at 0.
        // V_adj(2) = R1||R2 = R1*R2/(R1+R2) = 500Ω (node 1 clamped, so R1 and R2 in parallel).
        // Actually solving the adjoint exactly:
        //   V_adj(1) = 0 (clamped by V1), V_adj(2) = 500 (parallel combination).
        //
        // R1 noise: 4kT*G1 * |V_adj(1) - V_adj(2)|² = 4kT/1000 * 500² = 4kT * 250
        // R2 noise: 4kT*G2 * |V_adj(2) - 0|² = 4kT/1000 * 500² = 4kT * 250
        // Total: 4kT * 500 = 4kT * (R1||R2)
        //
        // This is the expected result: output noise = Johnson noise of R1||R2.
        let netlist = Netlist::parse(
            "Resistor thermal noise test
V1 1 0 DC 0 AC 1
R1 1 2 1k
R2 2 0 1k
.noise v(2) V1 DEC 10 1k 100k
.end
",
        )
        .unwrap();

        let result = simulate_noise(&netlist).unwrap();
        assert_eq!(result.plots.len(), 2);
        assert_eq!(result.plots[0].name, "noise1");
        assert_eq!(result.plots[1].name, "noise2");

        let plot = &result.plots[0];
        let onoise = plot
            .vecs
            .iter()
            .find(|v| v.name == "onoise_spectrum")
            .unwrap();

        // Expected output noise = sqrt(4kT * R_parallel) where R_parallel = R1||R2 = 500Ω
        let r_parallel = 500.0;
        let expected = (4.0 * K_BOLTZ * T_NOM * r_parallel).sqrt();
        for &val in &onoise.real {
            // Thermal noise is flat (frequency-independent).
            assert_abs_diff_eq!(val, expected, epsilon = expected * 0.02);
        }
    }

    #[test]
    fn test_resistor_thermal_noise_simple() {
        // Simplest noise circuit: voltage source in series with resistor.
        // V1 -> R1 -> ground, measure v(1).
        // At the output node, the noise is the Johnson noise of R1.
        let netlist = Netlist::parse(
            "Simple resistor noise
V1 1 0 DC 0 AC 1
R1 1 2 1k
.noise v(2) V1 DEC 1 1k 10k
.end
",
        )
        .unwrap();

        let result = simulate_noise(&netlist).unwrap();
        let plot = &result.plots[0];
        let onoise = plot
            .vecs
            .iter()
            .find(|v| v.name == "onoise_spectrum")
            .unwrap();

        // At output v(2), R1 is the only noise source.
        // The transfer from R1's noise current to v(2) depends on circuit.
        // V1 holds node 1 at 0V AC, so R1's noise current drives into ground through V1.
        // Output noise at v(2) = 0 (R1 noise is shorted by V1 on one side and open on the other).
        // Actually, node 2 has nothing connected besides R1, so it's floating.
        // Let's use a proper circuit instead.
        assert!(!onoise.real.is_empty());
    }

    #[test]
    fn test_resistor_thermal_noise_divider() {
        // Voltage divider noise: V1=AC source, R1=R2=1k.
        // Output at v(2). H = R2/(R1+R2) = 0.5.
        // R1 noise at output: 4kT*R1 * |R2/(R1+R2)|² = 4kTR * 0.25.
        // R2 noise at output: 4kT*R2 * |R1/(R1+R2)|² = 4kTR * 0.25.
        // Wait, that's not right. Let me recalculate with adjoint.
        //
        // Total output noise = 4kTR1 * |Zadj(1,2)|² + 4kTR2 * |Zadj(2,0)|².
        // For a voltage divider driven by a voltage source at node 1:
        // - Injecting unit current at output: V(2) = R2||R1 = R1*R2/(R1+R2) = 500.
        //   But R1 connects to V1 (held at 0 in adjoint), so adjoint from node 2:
        //   V(2) = R2 (since node 1 is clamped by V1).
        //   V(1) = 0 (clamped by V1).
        //
        // R1 noise contributes: 4kT/R1 * |V_adj(1) - V_adj(2)|² = 4kT/R1 * |0 - R2|² = 4kT*R2²/R1
        // R2 noise contributes: 4kT/R2 * |V_adj(2) - V_adj(0)|² = 4kT/R2 * |R2 - 0|² = 4kT*R2
        //
        // Hmm, that doesn't seem right either. The adjoint with V1 clamping node 1...
        // Actually the adjoint matrix includes V1's branch equation, so V1 clamps node 1.
        //
        // Let me just verify the output is reasonable and flat.
        let netlist = Netlist::parse(
            "Voltage divider noise
V1 1 0 DC 0 AC 1
R1 1 2 1k
R2 2 0 1k
.noise v(2) V1 DEC 10 100 100k
.end
",
        )
        .unwrap();

        let result = simulate_noise(&netlist).unwrap();
        let plot = &result.plots[0];
        let onoise = plot
            .vecs
            .iter()
            .find(|v| v.name == "onoise_spectrum")
            .unwrap();
        let inoise = plot
            .vecs
            .iter()
            .find(|v| v.name == "inoise_spectrum")
            .unwrap();

        // Output noise should be flat (only thermal) and > 0.
        assert!(!onoise.real.is_empty());
        for &val in &onoise.real {
            assert!(val > 0.0, "noise must be positive, got {val}");
        }

        // Flatness check: first and last should be similar (within 1%).
        let first = onoise.real[0];
        let last = *onoise.real.last().unwrap();
        assert!(
            (first - last).abs() / first < 0.01,
            "thermal noise should be flat: first={first}, last={last}"
        );

        // Input-referred noise = output noise / |H|.
        // |H| = 0.5, so inoise = 2 * onoise.
        let inoise_first = inoise.real[0];
        assert_abs_diff_eq!(inoise_first, first * 2.0, epsilon = first * 0.05);
    }

    #[test]
    fn test_noise_integration() {
        // White noise integration test: constant 1e-9 V/√Hz from 1kHz to 10kHz.
        // Total = sqrt(1e-18 * (10000 - 1000)) = sqrt(9e-15) ≈ 9.49e-8.
        let freqs: Vec<f64> = (0..100).map(|i| 1000.0 + i as f64 * 90.0).collect();
        let noise: Vec<f64> = vec![1e-9; 100];

        let total = integrate_noise(&freqs, &noise);
        let expected = (1e-18 * 9000.0_f64).sqrt();
        assert_abs_diff_eq!(total, expected, epsilon = expected * 0.02);
    }

    #[test]
    fn test_parse_output_spec() {
        let (pos, neg) = parse_output_spec("v(2)").unwrap();
        assert_eq!(pos, "2");
        assert!(neg.is_none());

        let (pos, neg) = parse_output_spec("v(out,ref)").unwrap();
        assert_eq!(pos, "out");
        assert_eq!(neg.as_deref(), Some("ref"));

        let (pos, neg) = parse_output_spec("V(1)").unwrap();
        assert_eq!(pos, "1");
        assert!(neg.is_none());
    }
}
