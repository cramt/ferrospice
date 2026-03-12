use std::f64::consts::PI;

use thevenin_types::{
    AcVariation, Analysis, Complex, Item, Netlist, SimPlot, SimResult, SimVector,
};

use crate::expr_val;
use crate::expr_val_or;
use crate::mna::{MnaError, MnaSystem, assemble_mna};
use crate::simulate::{nr_options_from_netlist, solve_nonlinear_op};
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

    // Assemble the MNA system and solve DC operating point directly
    // to get the full solution vector including internal device nodes.
    let mna = assemble_mna(netlist)?;
    let nr_opts = nr_options_from_netlist(netlist);
    let op_solution = if mna.has_nonlinear() {
        solve_nonlinear_op(&mna, &nr_opts)?
    } else {
        mna.system.solve()?
    };

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
        let num_nodes = mna.total_num_nodes();
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
    let num_nodes = mna.total_num_nodes();
    let dim = num_nodes + mna.vsource_names.len();

    let mut sys = ComplexLinearSystem::new(dim);

    // Stamp base matrix, reactive elements, and all device small-signal models.
    stamp_ac_devices(mna, op_solution, omega, &mut sys);

    // Apply AC source excitation to RHS.
    apply_ac_excitation(&mut sys, netlist, mna, num_nodes);

    sys.solve().map_err(MnaError::SolveError)
}

/// Build a complex AC MNA system at a given frequency without solving it.
///
/// Used by noise analysis which needs the unsolved system to solve both
/// forward and adjoint (transposed) systems.
pub fn build_ac_system(
    mna: &MnaSystem,
    op_solution: &[f64],
    omega: f64,
    netlist: &Netlist,
    num_nodes: usize,
) -> ComplexLinearSystem {
    let dim = num_nodes + mna.vsource_names.len();
    let mut sys = ComplexLinearSystem::new(dim);

    stamp_ac_devices(mna, op_solution, omega, &mut sys);
    apply_ac_excitation(&mut sys, netlist, mna, num_nodes);

    sys
}

/// Stamp all AC small-signal device models into a complex MNA system.
///
/// This includes:
/// - DC conductance matrix (real part from base MNA stamps)
/// - Capacitor and inductor reactive stamps (imaginary part)
/// - Diode, BJT, MOSFET, JFET, BSIM3, BSIM4 small-signal models at the DC operating point
pub fn stamp_ac_devices(
    mna: &MnaSystem,
    op_solution: &[f64],
    omega: f64,
    sys: &mut ComplexLinearSystem,
) {
    // 1. Stamp DC conductance matrix (real part) from base MNA stamps.
    for triplet in mna.system.matrix.triplets() {
        sys.real.add(triplet.row, triplet.col, triplet.value);
    }

    // 2. Stamp capacitor AC admittance: jωC between nodes.
    for cap in &mna.capacitors {
        let bc = omega * cap.capacitance;
        stamp_imag_conductance(&mut sys.imag, cap.pos_idx, cap.neg_idx, bc);
    }

    // 3. Stamp inductor AC impedance: -jωL on branch diagonal.
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

        let v_anode = jct_anode.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_cathode = jct_cathode.map(|i| op_solution[i]).unwrap_or(0.0);
        let v_jct = v_anode - v_cathode;

        let (g_d, _i_eq) = diode.model.companion(v_jct);
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

    // 5. Stamp BJT small-signal model at DC operating point.
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

    // 5b. Stamp MOSFET small-signal model at DC operating point.
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

    // 5b2. Stamp MOS6 small-signal model at DC operating point.
    for mos in &mna.mos6s {
        let (vgs, vds, vbs) = mos.terminal_voltages(op_solution);
        let betac = mos.betac();
        let comp = mos.model.companion(vgs, vds, vbs, betac);
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

    // 5c. Stamp JFET small-signal model at DC operating point.
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

    // 5d. Stamp BSIM3 small-signal model at DC operating point.
    for bsim in &mna.bsim3s {
        let (vgs, vds, vbs) = bsim.terminal_voltages(op_solution);
        let comp = crate::bsim3::bsim3_companion(vgs, vds, vbs, &bsim.size_params, &bsim.model);
        stamp_bsim_ac(&bsim.ac_stamp(&comp), omega, sys);
    }

    // 5e. Stamp BSIM3SOI-PD small-signal model at DC operating point.
    for bsim in &mna.bsim3soi_pds {
        let (vgs, vds, vbs, ves) = bsim.terminal_voltages(op_solution);
        let comp = crate::bsim3soi_pd::bsim3soi_pd_companion(
            vgs,
            vds,
            vbs,
            ves,
            &bsim.size_params,
            &bsim.model,
        );
        stamp_bsim_ac(&bsim.ac_stamp(&comp), omega, sys);
    }

    // 5e2. Stamp BSIM3SOI-FD small-signal model at DC operating point.
    for bsim in &mna.bsim3soi_fds {
        let (vgs, vds, vbs, ves) = bsim.terminal_voltages(op_solution);
        let comp = crate::bsim3soi_fd::bsim3soi_fd_companion(
            vgs,
            vds,
            vbs,
            ves,
            &bsim.size_params,
            &bsim.model,
            bsim.body_idx.is_none(),
        );
        stamp_bsim_ac(&bsim.ac_stamp(&comp), omega, sys);
    }

    // 5e3. Stamp BSIM3SOI-DD small-signal model at DC operating point.
    for bsim in &mna.bsim3soi_dds {
        let (vgs, vds, vbs, ves) = bsim.terminal_voltages(op_solution);
        let comp = crate::bsim3soi_dd::bsim3soi_dd_companion(
            vgs,
            vds,
            vbs,
            ves,
            &bsim.size_params,
            &bsim.model,
        );
        stamp_bsim_ac(&bsim.ac_stamp(&comp), omega, sys);
    }

    // 5f. Stamp BSIM4 small-signal model at DC operating point.
    for bsim in &mna.bsim4s {
        let (vgs, vds, vbs) = bsim.terminal_voltages(op_solution);
        let comp = crate::bsim4::bsim4_companion(vgs, vds, vbs, &bsim.size_params, &bsim.model);
        stamp_bsim_ac(&bsim.ac_stamp(&comp), omega, sys);
    }

    // 5g. Stamp VBIC small-signal model at DC operating point.
    for vbic in &mna.vbics {
        let (vbei, vbex, vbci, vbcx, vbep, vrci, vrbi, vrbp, vbcp) =
            vbic.junction_voltages(op_solution);
        let comp = vbic
            .model
            .companion(vbei, vbex, vbci, vbcx, vbep, vrci, vrbi, vrbp, vbcp);
        let s = vbic.m * vbic.area;

        let bi = vbic.base_bi_idx;
        let bx = vbic.base_bx_idx;
        let bp = vbic.base_bp_idx;
        let ci = vbic.coll_ci_idx;
        let cx = vbic.coll_cx_idx;
        let ei = vbic.emit_ei_idx;
        let si = vbic.subs_si_idx;

        // Stamp all conductances into real part using stamp_conductance
        // Ibe: BI->EI
        crate::stamp_conductance(&mut sys.real, bi, ei, s * comp.dibe_dvbei);
        // Ibex: BX->EI
        crate::stamp_conductance(&mut sys.real, bx, ei, s * comp.dibex_dvbex);
        // Ibc: BI->CI
        crate::stamp_conductance(&mut sys.real, bi, ci, s * comp.dibc_dvbci);
        // Ibep: BX->BP (controlled by Vbep = V(BX)-V(BP))
        crate::stamp_conductance(&mut sys.real, bx, bp, s * comp.dibep_dvbep);
        // Ibcp: SI->BP (controlled by Vbcp = V(SI)-V(BP))
        crate::stamp_conductance(&mut sys.real, si, bp, s * comp.dibcp_dvbcp);
        // Irci: CX->CI
        crate::stamp_conductance(&mut sys.real, cx, ci, s * comp.dirci_dvrci);
        // Irbi: BX->BI
        crate::stamp_conductance(&mut sys.real, bx, bi, s * comp.dirbi_dvrbi);
        // Irbp: BP->CX
        crate::stamp_conductance(&mut sys.real, bp, cx, s * comp.dirbp_dvrbp);

        // VCCS terms (transconductances)
        // Iciei controlled by Vbei (bi-ei) -> current ci->ei
        stamp_vccs(&mut sys.real, ci, ei, bi, ei, s * comp.diciei_dvbei);
        // Iciei controlled by Vbci (bi-ci) -> current ci->ei
        stamp_vccs(&mut sys.real, ci, ei, bi, ci, s * comp.diciei_dvbci);
        // Ibc avalanche cross-term: controlled by Vbei -> current bi->ci
        if comp.dibc_dvbei.abs() > 0.0 {
            stamp_vccs(&mut sys.real, bi, ci, bi, ei, s * comp.dibc_dvbei);
        }
        // Iccp controlled by Vbep -> current bx->si (Vbep = V(BX)-V(BP))
        stamp_vccs(&mut sys.real, bx, si, bx, bp, s * comp.diccp_dvbep);
        // Iccp controlled by Vbcp -> current bx->si (Vbcp = V(SI)-V(BP))
        stamp_vccs(&mut sys.real, bx, si, si, bp, s * comp.diccp_dvbcp);
        // Irci cross-terms
        if comp.dirci_dvbci.abs() > 0.0 {
            stamp_vccs(&mut sys.real, cx, ci, bi, ci, s * comp.dirci_dvbci);
        }
        if comp.dirci_dvbcx.abs() > 0.0 {
            stamp_vccs(&mut sys.real, cx, ci, bi, cx, s * comp.dirci_dvbcx);
        }
        // Irbi cross-terms
        if comp.dirbi_dvbei.abs() > 0.0 {
            stamp_vccs(&mut sys.real, bx, bi, bi, ei, s * comp.dirbi_dvbei);
        }
        if comp.dirbi_dvbci.abs() > 0.0 {
            stamp_vccs(&mut sys.real, bx, bi, bi, ci, s * comp.dirbi_dvbci);
        }

        // External resistances
        if vbic.model.rcx > 0.0 && vbic.model.rcx_t > 0.0 {
            crate::stamp_conductance(&mut sys.real, vbic.coll_idx, cx, s / vbic.model.rcx_t);
        }
        if vbic.model.rbx > 0.0 && vbic.model.rbx_t > 0.0 {
            crate::stamp_conductance(&mut sys.real, vbic.base_idx, bx, s / vbic.model.rbx_t);
        }
        if vbic.model.re > 0.0 && vbic.model.re_t > 0.0 {
            crate::stamp_conductance(&mut sys.real, vbic.emit_idx, ei, s / vbic.model.re_t);
        }
        if vbic.model.rs > 0.0 && vbic.model.rs_t > 0.0 {
            crate::stamp_conductance(&mut sys.real, vbic.subs_idx, si, s / vbic.model.rs_t);
        }

        // Capacitances (imaginary part) — total charge derivatives including transit time
        // Qbe: BI-EI (main diagonal) + BI-CI (cross-term from transit time)
        let xqbe_vbei = omega * s * comp.cqbe_vbei;
        if xqbe_vbei.abs() > 0.0 {
            stamp_imag_conductance(&mut sys.imag, bi, ei, xqbe_vbei);
        }
        let xqbe_vbci = omega * s * comp.cqbe_vbci;
        if xqbe_vbci.abs() > 0.0 {
            // Cross-term: charge at BI-EI controlled by Vbci (BI-CI)
            stamp_vccs(&mut sys.imag, bi, ei, bi, ci, xqbe_vbci);
        }

        // Qbex: BX-EI (external BE depletion)
        let xqbex_vbex = omega * s * comp.cqbex_vbex;
        if xqbex_vbex.abs() > 0.0 {
            stamp_imag_conductance(&mut sys.imag, bx, ei, xqbex_vbex);
        }

        // Qbc: BI-CI
        let xqbc_vbci = omega * s * comp.cqbc_vbci;
        if xqbc_vbci.abs() > 0.0 {
            stamp_imag_conductance(&mut sys.imag, bi, ci, xqbc_vbci);
        }

        // Qbcx: BI-CX (Vbcx = V(BI)-V(CX))
        let xqbcx_vbcx = omega * s * comp.cqbcx_vbcx;
        if xqbcx_vbcx.abs() > 0.0 {
            stamp_imag_conductance(&mut sys.imag, bi, cx, xqbcx_vbcx);
        }

        // Qbep: BX-BP (main, Vbep = V(BX)-V(BP)) + cross from Vbci
        let xqbep_vbep = omega * s * comp.cqbep_vbep;
        if xqbep_vbep.abs() > 0.0 {
            stamp_imag_conductance(&mut sys.imag, bx, bp, xqbep_vbep);
        }
        let xqbep_vbci = omega * s * comp.cqbep_vbci;
        if xqbep_vbci.abs() > 0.0 {
            stamp_vccs(&mut sys.imag, bx, bp, bi, ci, xqbep_vbci);
        }

        // Qbcp: SI-BP (Vbcp = V(SI)-V(BP))
        let xqbcp_vbcp = omega * s * comp.cqbcp_vbcp;
        if xqbcp_vbcp.abs() > 0.0 {
            stamp_imag_conductance(&mut sys.imag, si, bp, xqbcp_vbcp);
        }

        // Overlap capacitances (external nodes)
        // Qbeo = CBEO * Vbe (B-E external)
        if vbic.model.cbeo > 0.0 {
            let xqbeo = omega * s * vbic.model.cbeo;
            stamp_imag_conductance(&mut sys.imag, vbic.base_idx, vbic.emit_idx, xqbeo);
        }
        // Qbco = CBCO * Vbc (B-C external)
        if vbic.model.cbco > 0.0 {
            let xqbco = omega * s * vbic.model.cbco;
            stamp_imag_conductance(&mut sys.imag, vbic.base_idx, vbic.coll_idx, xqbco);
        }
    }

    // 5h. Stamp MESA FET small-signal model at DC operating point.
    for mesa in &mna.mesas {
        let (vgs, vgd) = mesa.junction_voltages(op_solution);
        let comp = crate::mesa::mesa_companion(mesa, vgs, vgd, 1e-12);

        let f = omega / (2.0 * PI);
        let pre = &mesa.precomp;

        // Frequency-dependent lambda
        let lambda = if pre.delf == 0.0 {
            pre.t_lambda
        } else {
            pre.t_lambda
                + 0.5 * (pre.t_lambdahf - pre.t_lambda) * (1.0 + ((f - pre.fl) / pre.delf).tanh())
        };

        // Recompute gm/gds with frequency-dependent lambda
        let vds = vgs - vgd;
        let delidgch = comp.delidgch0 * (1.0 + lambda * vds);
        let delidvds = comp.delidvds0 * (1.0 + 2.0 * lambda * vds) - comp.delidvds1;
        let gm = (delidgch * comp.gm0 + comp.gm1) * comp.gm2;
        let gds = delidvds + comp.gds0;

        let dp = mesa.drain_prime_idx;
        let gp = mesa.gate_prime_idx;
        let sp = mesa.source_prime_idx;
        let spp = mesa.source_prm_prm_idx;
        let dpp = mesa.drain_prm_prm_idx;

        // Real stamps - series resistances (individual entries)
        if let Some(d) = mesa.drain_idx {
            sys.real.add(d, d, pre.drain_conduct);
        }
        if let Some(s) = mesa.source_idx {
            sys.real.add(s, s, pre.source_conduct);
        }
        if let Some(g) = mesa.gate_idx {
            sys.real.add(g, g, pre.gate_conduct);
        }
        if let (Some(d), Some(di)) = (mesa.drain_idx, dp) {
            sys.real.add(d, di, -pre.drain_conduct);
            sys.real.add(di, d, -pre.drain_conduct);
        }
        if let (Some(s), Some(si)) = (mesa.source_idx, sp) {
            sys.real.add(s, si, -pre.source_conduct);
            sys.real.add(si, s, -pre.source_conduct);
        }
        if let (Some(g), Some(gi)) = (mesa.gate_idx, gp) {
            sys.real.add(g, gi, -pre.gate_conduct);
            sys.real.add(gi, g, -pre.gate_conduct);
        }

        // spp, dpp diagonal (Gi + ggspp, Gf + ggdpp; spp/dpp terms 0 in DC)
        if let Some(spp_i) = spp {
            sys.real.add(spp_i, spp_i, pre.t_gi);
        }
        if let Some(dpp_i) = dpp {
            sys.real.add(dpp_i, dpp_i, pre.t_gf);
        }

        // Gate' diagonal
        if let Some(gp_i) = gp {
            sys.real
                .add(gp_i, gp_i, comp.ggd + comp.ggs + pre.gate_conduct);
        }
        // Drain' diagonal
        if let Some(dp_i) = dp {
            sys.real
                .add(dp_i, dp_i, gds + comp.ggd + pre.drain_conduct + pre.t_gf);
        }
        // Source' diagonal
        if let Some(sp_i) = sp {
            sys.real.add(
                sp_i,
                sp_i,
                gds + gm + comp.ggs + pre.source_conduct + pre.t_gi,
            );
        }

        // Off-diagonal real (individual matrix entries)
        if let (Some(gi), Some(di)) = (gp, dp) {
            sys.real.add(gi, di, -comp.ggd);
        }
        if let (Some(gi), Some(si)) = (gp, sp) {
            sys.real.add(gi, si, -comp.ggs);
        }
        if let (Some(di), Some(gi)) = (dp, gp) {
            sys.real.add(di, gi, gm - comp.ggd);
        }
        if let (Some(di), Some(si)) = (dp, sp) {
            sys.real.add(di, si, -gds - gm);
        }
        if let (Some(si), Some(gi)) = (sp, gp) {
            sys.real.add(si, gi, -comp.ggs - gm);
        }
        if let (Some(si), Some(di)) = (sp, dp) {
            sys.real.add(si, di, -gds);
        }

        // sp <-> spp, dp <-> dpp coupling
        if let (Some(si), Some(spi)) = (sp, spp) {
            sys.real.add(si, spi, -pre.t_gi);
            sys.real.add(spi, si, -pre.t_gi);
        }
        if let (Some(di), Some(dpi)) = (dp, dpp) {
            sys.real.add(di, dpi, -pre.t_gf);
            sys.real.add(dpi, di, -pre.t_gf);
        }

        // Imaginary stamps - capacitances on spp/dpp nodes
        let xgs = comp.capgs * omega;
        let xgd = comp.capgd * omega;
        if xgs != 0.0 {
            if let Some(spi) = spp {
                sys.imag.add(spi, spi, xgs);
            }
            if let Some(gi) = gp {
                sys.imag.add(gi, gi, xgs);
            }
            if let (Some(gi), Some(spi)) = (gp, spp) {
                sys.imag.add(gi, spi, -xgs);
                sys.imag.add(spi, gi, -xgs);
            }
        }
        if xgd != 0.0 {
            if let Some(dpi) = dpp {
                sys.imag.add(dpi, dpi, xgd);
            }
            if let Some(gi) = gp {
                sys.imag.add(gi, gi, xgd);
            }
            if let (Some(gi), Some(dpi)) = (gp, dpp) {
                sys.imag.add(gi, dpi, -xgd);
                sys.imag.add(dpi, gi, -xgd);
            }
        }
    }

    // 5i. Stamp MESFET small-signal model at DC operating point.
    for mes in &mna.mesfets {
        let (vgs, vgd) = mes.junction_voltages(op_solution);
        let comp = crate::mesfet::mesfet_companion(mes, vgs, vgd, 1e-12);

        let gate = mes.gate_idx;
        let dp = mes.drain_prime_idx.or(mes.drain_idx);
        let sp = mes.source_prime_idx.or(mes.source_idx);

        // Series resistance stamps (use stamp_conductance to handle ground nodes)
        if mes.drain_prime_idx.is_some() {
            let gdpr = mes.model.drain_conduct * mes.area;
            crate::mna::stamp_conductance(&mut sys.real, mes.drain_idx, mes.drain_prime_idx, gdpr);
        }
        if mes.source_prime_idx.is_some() {
            let gspr = mes.model.source_conduct * mes.area;
            crate::mna::stamp_conductance(
                &mut sys.real,
                mes.source_idx,
                mes.source_prime_idx,
                gspr,
            );
        }

        // Channel conductance stamps
        if let Some(g) = gate {
            sys.real.add(g, g, comp.ggs + comp.ggd);
        }
        if let Some(dpi) = dp {
            sys.real.add(dpi, dpi, comp.gds + comp.ggd);
        }
        if let Some(spi) = sp {
            sys.real.add(spi, spi, comp.gds + comp.gm + comp.ggs);
        }
        if let (Some(g), Some(dpi)) = (gate, dp) {
            sys.real.add(g, dpi, -comp.ggd);
            sys.real.add(dpi, g, comp.gm - comp.ggd);
        }
        if let (Some(g), Some(spi)) = (gate, sp) {
            sys.real.add(g, spi, -comp.ggs);
            sys.real.add(spi, g, -comp.ggs - comp.gm);
        }
        if let (Some(dpi), Some(spi)) = (dp, sp) {
            sys.real.add(dpi, spi, -comp.gds - comp.gm);
            sys.real.add(spi, dpi, -comp.gds);
        }

        // Capacitance stamps
        let xgs = comp.capgs * omega;
        let xgd = comp.capgd * omega;
        if xgs != 0.0 {
            stamp_imag_conductance(&mut sys.imag, gate, sp, xgs);
        }
        if xgd != 0.0 {
            stamp_imag_conductance(&mut sys.imag, gate, dp, xgd);
        }
    }

    // 5j. Stamp HFET small-signal model at DC operating point.
    for hfet_inst in &mna.hfets {
        let (vgs, vgd) = hfet_inst.junction_voltages(op_solution);
        let comp = crate::hfet::hfet_companion_full(hfet_inst, vgs, vgd, 1e-12);

        let gp = hfet_inst.gate_prime_idx.or(hfet_inst.gate_idx);
        let dp = hfet_inst.drain_prime_idx.or(hfet_inst.drain_idx);
        let sp = hfet_inst.source_prime_idx.or(hfet_inst.source_idx);

        // Series resistance stamps (use stamp_conductance to handle ground nodes)
        if hfet_inst.drain_prime_idx.is_some() {
            crate::mna::stamp_conductance(
                &mut sys.real,
                hfet_inst.drain_idx,
                hfet_inst.drain_prime_idx,
                hfet_inst.model.drain_conduct,
            );
        }
        if hfet_inst.source_prime_idx.is_some() {
            crate::mna::stamp_conductance(
                &mut sys.real,
                hfet_inst.source_idx,
                hfet_inst.source_prime_idx,
                hfet_inst.model.source_conduct,
            );
        }
        if hfet_inst.gate_prime_idx.is_some() {
            crate::mna::stamp_conductance(
                &mut sys.real,
                hfet_inst.gate_idx,
                hfet_inst.gate_prime_idx,
                hfet_inst.model.gate_conduct,
            );
        }

        // Channel conductance stamps
        if let Some(gpi) = gp {
            sys.real.add(gpi, gpi, comp.ggs + comp.ggd);
        }
        if let Some(dpi) = dp {
            sys.real.add(dpi, dpi, comp.gds + comp.ggd);
        }
        if let Some(spi) = sp {
            sys.real.add(spi, spi, comp.gds + comp.gm + comp.ggs);
        }
        if let (Some(gpi), Some(dpi)) = (gp, dp) {
            sys.real.add(gpi, dpi, -comp.ggd);
            sys.real.add(dpi, gpi, comp.gm - comp.ggd);
        }
        if let (Some(gpi), Some(spi)) = (gp, sp) {
            sys.real.add(gpi, spi, -comp.ggs);
            sys.real.add(spi, gpi, -comp.ggs - comp.gm);
        }
        if let (Some(dpi), Some(spi)) = (dp, sp) {
            sys.real.add(dpi, spi, -comp.gds - comp.gm);
            sys.real.add(spi, dpi, -comp.gds);
        }

        // Feedback resistances (ri, rf)
        if let (Some(dpi), Some(dppi)) = (hfet_inst.drain_prime_idx, hfet_inst.drain_prm_prm_idx) {
            let g = hfet_inst.model.gf;
            sys.real.add(dpi, dpi, g);
            sys.real.add(dppi, dppi, g);
            sys.real.add(dpi, dppi, -g);
            sys.real.add(dppi, dpi, -g);
        }
        if let (Some(spi), Some(sppi)) = (hfet_inst.source_prime_idx, hfet_inst.source_prm_prm_idx)
        {
            let g = hfet_inst.model.gi;
            sys.real.add(spi, spi, g);
            sys.real.add(sppi, sppi, g);
            sys.real.add(spi, sppi, -g);
            sys.real.add(sppi, spi, -g);
        }

        // Capacitance stamps
        let spp = hfet_inst.source_prm_prm_idx.or(sp);
        let dpp = hfet_inst.drain_prm_prm_idx.or(dp);
        let xgs = comp.capgs * omega;
        let xgd = comp.capgd * omega;
        if xgs != 0.0 {
            stamp_imag_conductance(&mut sys.imag, gp, spp, xgs);
        }
        if xgd != 0.0 {
            stamp_imag_conductance(&mut sys.imag, gp, dpp, xgd);
        }
    }
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
            thevenin_types::ElementKind::VoltageSource { source, .. } => {
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
            thevenin_types::ElementKind::CurrentSource { pos, neg, source } => {
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

/// AC small-signal data for a 4-terminal MOSFET-family device (BSIM3/BSIM4/SOI).
///
/// Captures the node indices, scale factor, conductances, and capacitances needed
/// to stamp a BSIM-family device into the complex AC MNA matrix.
pub struct BsimAcStamp {
    pub dp: Option<usize>,
    pub g: Option<usize>,
    pub sp: Option<usize>,
    pub b: Option<usize>,
    pub drain_idx: Option<usize>,
    pub source_idx: Option<usize>,
    pub m: f64,
    pub gds: f64,
    pub gm: f64,
    pub gmbs: f64,
    pub gbd: f64,
    pub gbs: f64,
    pub g_drain: f64,
    pub g_source: f64,
    pub cggb: f64,
    pub cgdb: f64,
    pub cgsb: f64,
    pub cbgb: f64,
    pub cbdb: f64,
    pub cbsb: f64,
    pub cdgb: f64,
    pub cddb: f64,
    pub cdsb: f64,
    pub capbd: f64,
    pub capbs: f64,
}

/// Stamp a BSIM-family device's small-signal model into the complex AC system.
///
/// Used by BSIM3, BSIM4, and future BSIM-SOI variants to avoid duplicating
/// the identical gm/gds/gmbs/cap stamping logic.
pub fn stamp_bsim_ac(stamp: &BsimAcStamp, omega: f64, sys: &mut ComplexLinearSystem) {
    let BsimAcStamp {
        dp,
        g,
        sp,
        b,
        drain_idx,
        source_idx,
        m,
        gds,
        gm,
        gmbs,
        gbd,
        gbs,
        g_drain,
        g_source,
        cggb,
        cgdb,
        cgsb,
        cbgb,
        cbdb,
        cbsb,
        cdgb,
        cddb,
        cdsb,
        capbd,
        capbs,
    } = *stamp;

    // gds conductance d'-s' (real)
    crate::stamp_conductance(&mut sys.real, dp, sp, m * gds);

    // gm VCCS: Vgs controls current s'→d' (real)
    let gm_scaled = m * gm;
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

    // gmbs: Vbs controls current s'→d' (real)
    let gmbs_scaled = m * gmbs;
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

    // gbd conductance b-d' (real)
    crate::stamp_conductance(&mut sys.real, b, dp, m * gbd);
    // gbs conductance b-s' (real)
    crate::stamp_conductance(&mut sys.real, b, sp, m * gbs);

    // Series resistances (real)
    if g_drain > 0.0 {
        crate::stamp_conductance(&mut sys.real, drain_idx, dp, m * g_drain);
    }
    if g_source > 0.0 {
        crate::stamp_conductance(&mut sys.real, source_idx, sp, m * g_source);
    }

    // Intrinsic capacitances (imaginary) - direct Y-matrix entries [G, D', S', B]
    let cap_entries: &[(Option<usize>, Option<usize>, f64)] = &[
        (g, g, cggb),
        (g, dp, cgdb),
        (g, sp, cgsb),
        (dp, g, cdgb),
        (dp, dp, cddb),
        (dp, sp, cdsb),
        (b, g, cbgb),
        (b, dp, cbdb),
        (b, sp, cbsb),
    ];
    for &(row, col, cap) in cap_entries {
        if cap.abs() > 0.0
            && let (Some(r), Some(c)) = (row, col)
        {
            sys.imag.add(r, c, omega * m * cap);
        }
    }

    // Junction capacitances (imaginary)
    if capbd > 0.0 {
        stamp_imag_conductance(&mut sys.imag, b, dp, omega * m * capbd);
    }
    if capbs > 0.0 {
        stamp_imag_conductance(&mut sys.imag, b, sp, omega * m * capbs);
    }
}

/// Stamp a VCCS (gm * V(ctrl+ - ctrl-) injected at out+ - out-) into a sparse matrix.
fn stamp_vccs(
    matrix: &mut crate::SparseMatrix,
    out_pos: Option<usize>,
    out_neg: Option<usize>,
    ctrl_pos: Option<usize>,
    ctrl_neg: Option<usize>,
    gm: f64,
) {
    if let Some(i) = out_pos {
        if let Some(cp) = ctrl_pos {
            matrix.add(i, cp, gm);
        }
        if let Some(cm) = ctrl_neg {
            matrix.add(i, cm, -gm);
        }
    }
    if let Some(j) = out_neg {
        if let Some(cp) = ctrl_pos {
            matrix.add(j, cp, -gm);
        }
        if let Some(cm) = ctrl_neg {
            matrix.add(j, cm, gm);
        }
    }
}

/// Stamp an imaginary-part conductance (susceptance) between two nodes.
pub fn stamp_imag_conductance(
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
pub fn generate_ac_sweep(variation: AcVariation, n: u32, fstart: f64, fstop: f64) -> Vec<f64> {
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
pub fn extract_op_solution(op_result: &SimResult, mna: &MnaSystem) -> Vec<f64> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

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
