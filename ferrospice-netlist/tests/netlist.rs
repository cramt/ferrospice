//! Integration tests for `spice-netlist`.
//!
//! Covers parsing, display (generation), and round-trip fidelity for every
//! element type, every analysis command, every waveform, and all common dot
//! commands.

#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use approx::assert_abs_diff_eq;
use ferrospice_netlist::{
    AcSpec, AcVariation, Analysis, DcSweep, ElementKind, Expr, Item, Netlist, Source, Waveform,
    format_si,
};
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse(src: &str) -> Netlist {
    Netlist::parse(src).unwrap_or_else(|e| panic!("parse failed: {e}"))
}

fn num(src: &str) -> f64 {
    ferrospice_netlist::parse::parse_spice_number(src)
        .unwrap_or_else(|| panic!("not a number: {src}"))
}

fn assert_num(expr: &Expr, expected: f64) {
    if let Expr::Num(v) = expr {
        assert_abs_diff_eq!(*v, expected, epsilon = expected.abs() * 1e-9 + 1e-30);
    } else {
        panic!("expected Num, got {expr:?}");
    }
}

/// Parse → Display → re-parse; check that the title and element count survive.
fn roundtrip(src: &str) -> (Netlist, Netlist) {
    let n1 = parse(src);
    let rendered = n1.to_string();
    let n2 = parse(&rendered);
    assert_eq!(n1.title, n2.title, "title changed after round-trip");
    (n1, n2)
}

// ---------------------------------------------------------------------------
// SI suffix / number parsing
// ---------------------------------------------------------------------------

#[test]
fn si_all_suffixes() {
    assert_abs_diff_eq!(num("1T"), 1e12, epsilon = 1e3);
    assert_abs_diff_eq!(num("1G"), 1e9, epsilon = 1e0);
    assert_abs_diff_eq!(num("1Meg"), 1e6, epsilon = 1e-3);
    assert_abs_diff_eq!(num("1MEG"), 1e6, epsilon = 1e-3);
    assert_abs_diff_eq!(num("1k"), 1e3, epsilon = 1e-6);
    assert_abs_diff_eq!(num("1K"), 1e3, epsilon = 1e-6);
    assert_abs_diff_eq!(num("1"), 1.0, epsilon = 1e-9);
    assert_abs_diff_eq!(num("1m"), 1e-3, epsilon = 1e-12);
    assert_abs_diff_eq!(num("1M"), 1e-3, epsilon = 1e-12); // M = milli!
    assert_abs_diff_eq!(num("1u"), 1e-6, epsilon = 1e-15);
    assert_abs_diff_eq!(num("1U"), 1e-6, epsilon = 1e-15);
    assert_abs_diff_eq!(num("1n"), 1e-9, epsilon = 1e-18);
    assert_abs_diff_eq!(num("1p"), 1e-12, epsilon = 1e-21);
    assert_abs_diff_eq!(num("1f"), 1e-15, epsilon = 1e-24);
    assert_abs_diff_eq!(num("1a"), 1e-18, epsilon = 1e-27);
    assert_abs_diff_eq!(num("1mil"), 25.4e-6, epsilon = 1e-15);
}

#[test]
fn si_format_round_trip_values() {
    for val in [0.0, 1.0, -1.5, 1e3, 2.5e3, 1e6, 1e-3, 47e-9, 100e-12] {
        let s = format_si(val);
        if val != 0.0 {
            assert_abs_diff_eq!(num(&s), val, epsilon = val.abs() * 1e-6);
        }
    }
}

#[test]
fn number_with_trailing_unit_letters() {
    // "1kohm" → K suffix → 1000.0; trailing "ohm" ignored
    assert_abs_diff_eq!(num("1kohm"), 1e3, epsilon = 1e-9);
    assert_abs_diff_eq!(num("10nF"), 10e-9, epsilon = 1e-18);
    assert_abs_diff_eq!(num("100uA"), 100e-6, epsilon = 1e-15);
}

#[test]
fn scientific_notation() {
    assert_abs_diff_eq!(num("1.5e-9"), 1.5e-9, epsilon = 1e-18);
    assert_abs_diff_eq!(num("2E3"), 2000.0, epsilon = 1e-9);
    assert_abs_diff_eq!(num("-3.3e-3"), -3.3e-3, epsilon = 1e-12);
}

// ---------------------------------------------------------------------------
// Pre-processing
// ---------------------------------------------------------------------------

#[test]
fn continuation_folds_into_previous_line() {
    let n = parse("Title\n.model DIODE D\n+ IS=1e-14 N=1.05\n.end");
    if let Item::Model(m) = &n.items[0] {
        assert_eq!(m.params.len(), 2);
        assert_eq!(m.params[0].name, "IS");
        assert_eq!(m.params[1].name, "N");
    } else {
        panic!("expected Model");
    }
}

#[test]
fn inline_dollar_comment_stripped() {
    let n = parse("T\nR1 a b 4.7k $ bypass\n.end");
    if let ElementKind::Resistor { value, .. } = &n.elements().next().unwrap().kind {
        assert_num(value, 4700.0);
    }
}

#[test]
fn inline_semicolon_comment_stripped() {
    let n = parse("T\nC1 a b 100n ; decoupling\n.end");
    if let ElementKind::Capacitor { value, .. } = &n.elements().next().unwrap().kind {
        assert_num(value, 100e-9);
    }
}

#[test]
fn dollar_inside_pulse_parens_not_stripped() {
    // Ensure inline comment detection doesn't fire inside parentheses
    let n = parse("T\nV1 a 0 PULSE(0 1 1n 1n 1n 5n 10n)\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        assert!(matches!(source.waveform, Some(Waveform::Pulse { .. })));
    }
}

#[test]
fn full_line_star_comment() {
    let n = parse("Title\n* this is a comment\nR1 a b 1k\n.end");
    assert_eq!(n.items.len(), 2);
    assert!(matches!(&n.items[0], Item::Comment(_)));
    assert!(matches!(&n.items[1], Item::Element(_)));
}

#[test]
fn title_is_first_non_empty_line() {
    let n = parse("\n\nMy Title\nR1 a b 1\n.end");
    assert_eq!(n.title, "My Title");
}

// ---------------------------------------------------------------------------
// Resistor / Capacitor / Inductor
// ---------------------------------------------------------------------------

#[test]
fn r_c_l_parsed() {
    let n = parse("T\nR1 a b 1k\nC2 c d 100n\nL3 e f 1u\n.end");
    let elems: Vec<_> = n.elements().collect();
    assert_eq!(elems.len(), 3);

    if let ElementKind::Resistor {
        pos,
        neg,
        value,
        params,
    } = &elems[0].kind
    {
        assert_eq!(pos, "a");
        assert_eq!(neg, "b");
        assert_num(value, 1e3);
        assert!(params.is_empty());
    } else {
        panic!()
    }

    if let ElementKind::Capacitor { value, .. } = &elems[1].kind {
        assert_num(value, 100e-9);
    } else {
        panic!()
    }

    if let ElementKind::Inductor { value, .. } = &elems[2].kind {
        assert_num(value, 1e-6);
    } else {
        panic!()
    }
}

#[test]
fn r_with_params() {
    let n = parse("T\nR1 a b 1k TC1=0.001 TC2=5e-6\n.end");
    if let ElementKind::Resistor { params, .. } = &n.elements().next().unwrap().kind {
        assert_eq!(params.len(), 2);
        assert_eq!(params[0].name, "TC1");
        assert_eq!(params[1].name, "TC2");
    }
}

#[test]
fn c_with_ic() {
    let n = parse("T\nC1 a b 1u IC=2.5\n.end");
    if let ElementKind::Capacitor { params, .. } = &n.elements().next().unwrap().kind {
        assert_eq!(params.len(), 1);
        assert_eq!(params[0].name, "IC");
        assert_num(&params[0].value, 2.5);
    }
}

// ---------------------------------------------------------------------------
// Voltage / current sources
// ---------------------------------------------------------------------------

#[test]
fn v_dc_bare_value() {
    // No "DC" keyword — bare numeric is DC
    let n = parse("T\nV1 vcc 0 5\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        assert_num(source.dc.as_ref().unwrap(), 5.0);
        assert!(source.ac.is_none());
        assert!(source.waveform.is_none());
    }
}

#[test]
fn v_dc_keyword() {
    let n = parse("T\nV1 vcc 0 DC 3.3\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        assert_num(source.dc.as_ref().unwrap(), 3.3);
    }
}

#[test]
fn v_ac_only() {
    let n = parse("T\nV1 a 0 AC 1\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        assert!(source.dc.is_none());
        let ac = source.ac.as_ref().unwrap();
        assert_num(&ac.mag, 1.0);
        assert!(ac.phase.is_none());
    }
}

#[test]
fn v_ac_with_phase() {
    let n = parse("T\nV1 a 0 AC 1 90\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        let ac = source.ac.as_ref().unwrap();
        assert_num(&ac.mag, 1.0);
        assert_num(ac.phase.as_ref().unwrap(), 90.0);
    }
}

#[test]
fn v_dc_and_ac() {
    let n = parse("T\nV1 a 0 DC 0 AC 1 0\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        assert!(source.dc.is_some());
        assert!(source.ac.is_some());
        assert!(source.waveform.is_none());
    }
}

#[test]
fn i_source() {
    let n = parse("T\nI1 a b DC 2m\n.end");
    if let ElementKind::CurrentSource { source, .. } = &n.elements().next().unwrap().kind {
        assert_num(source.dc.as_ref().unwrap(), 2e-3);
    }
}

// ---------------------------------------------------------------------------
// Waveforms
// ---------------------------------------------------------------------------

#[test]
fn waveform_pulse_all_args() {
    let n = parse("T\nV1 a 0 PULSE(0 5 1n 2n 3n 10n 20n)\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        if let Some(Waveform::Pulse {
            v1,
            v2,
            td,
            tr,
            tf,
            pw,
            per,
        }) = &source.waveform
        {
            assert_num(v1, 0.0);
            assert_num(v2, 5.0);
            assert_num(td.as_ref().unwrap(), 1e-9);
            assert_num(tr.as_ref().unwrap(), 2e-9);
            assert_num(tf.as_ref().unwrap(), 3e-9);
            assert_num(pw.as_ref().unwrap(), 10e-9);
            assert_num(per.as_ref().unwrap(), 20e-9);
        } else {
            panic!("expected Pulse")
        }
    }
}

#[test]
fn waveform_pulse_minimal() {
    let n = parse("T\nV1 a 0 PULSE(0 1)\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        if let Some(Waveform::Pulse { v1, v2, td, .. }) = &source.waveform {
            assert_num(v1, 0.0);
            assert_num(v2, 1.0);
            assert!(td.is_none());
        } else {
            panic!()
        }
    }
}

#[test]
fn waveform_sin() {
    let n = parse("T\nV1 a 0 SIN(0 5 1k 0 0 0)\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        if let Some(Waveform::Sin { v0, va, freq, .. }) = &source.waveform {
            assert_num(v0, 0.0);
            assert_num(va, 5.0);
            assert_num(freq.as_ref().unwrap(), 1e3);
        } else {
            panic!("expected Sin")
        }
    }
}

#[test]
fn waveform_exp() {
    let n = parse("T\nV1 a 0 EXP(0 5 1n 10n 50n 100n)\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        assert!(matches!(source.waveform, Some(Waveform::Exp { .. })));
    }
}

#[test]
fn waveform_pwl() {
    let n = parse("T\nV1 a 0 PWL(0 0 1n 5 2n 5 3n 0)\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        if let Some(Waveform::Pwl(pts)) = &source.waveform {
            assert_eq!(pts.len(), 4);
            assert_num(&pts[0].time, 0.0);
            assert_num(&pts[0].value, 0.0);
            assert_num(&pts[1].time, 1e-9);
            assert_num(&pts[1].value, 5.0);
        } else {
            panic!("expected Pwl")
        }
    }
}

#[test]
fn waveform_sffm() {
    let n = parse("T\nV1 a 0 SFFM(0 1 1k 100)\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        if let Some(Waveform::Sffm { v0, va, fc, fs, md }) = &source.waveform {
            assert_num(v0, 0.0);
            assert_num(va, 1.0);
            assert_num(fc.as_ref().unwrap(), 1e3);
            assert_num(fs.as_ref().unwrap(), 100.0);
            assert!(md.is_none());
        } else {
            panic!("expected Sffm")
        }
    }
}

#[test]
fn waveform_am() {
    let n = parse("T\nV1 a 0 AM(1 0 10k 1k)\n.end");
    if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
        if let Some(Waveform::Am { va, vo, fc, fs, td }) = &source.waveform {
            assert_num(va, 1.0);
            assert_num(vo, 0.0);
            assert_num(fc, 10e3);
            assert_num(fs, 1e3);
            assert!(td.is_none());
        } else {
            panic!("expected Am")
        }
    }
}

// ---------------------------------------------------------------------------
// Semiconductors
// ---------------------------------------------------------------------------

#[test]
fn diode() {
    let n = parse("T\nD1 anode cathode 1N4148\n.end");
    if let ElementKind::Diode {
        anode,
        cathode,
        model,
        params,
    } = &n.elements().next().unwrap().kind
    {
        assert_eq!(anode, "anode");
        assert_eq!(cathode, "cathode");
        assert_eq!(model, "1N4148");
        assert!(params.is_empty());
    } else {
        panic!()
    }
}

#[test]
fn bjt_npn_no_substrate() {
    let n = parse("T\nQ1 c b e BC547\n.end");
    if let ElementKind::Bjt {
        c,
        b,
        e,
        substrate,
        model,
        ..
    } = &n.elements().next().unwrap().kind
    {
        assert_eq!(c, "c");
        assert_eq!(b, "b");
        assert_eq!(e, "e");
        assert!(substrate.is_none());
        assert_eq!(model, "BC547");
    } else {
        panic!()
    }
}

#[test]
fn bjt_with_substrate() {
    let n = parse("T\nQ1 c b e sub BC547\n.end");
    if let ElementKind::Bjt {
        substrate, model, ..
    } = &n.elements().next().unwrap().kind
    {
        assert_eq!(substrate.as_deref(), Some("sub"));
        assert_eq!(model, "BC547");
    } else {
        panic!()
    }
}

#[test]
fn mosfet_with_params() {
    let n = parse("T\nM1 drain gate source bulk NMOD W=10u L=180n\n.end");
    if let ElementKind::Mosfet {
        d,
        g,
        s,
        bulk,
        model,
        params,
    } = &n.elements().next().unwrap().kind
    {
        assert_eq!(d, "drain");
        assert_eq!(g, "gate");
        assert_eq!(s, "source");
        assert_eq!(bulk, "bulk");
        assert_eq!(model, "NMOD");
        assert_eq!(params.len(), 2);
        assert_eq!(params[0].name, "W");
        assert_num(&params[0].value, 10e-6);
        assert_eq!(params[1].name, "L");
        assert_num(&params[1].value, 180e-9);
    } else {
        panic!()
    }
}

#[test]
fn jfet() {
    let n = parse("T\nJ1 d g s J2N3819\n.end");
    if let ElementKind::Jfet {
        d,
        g,
        s,
        model,
        params,
    } = &n.elements().next().unwrap().kind
    {
        assert_eq!(d, "d");
        assert_eq!(g, "g");
        assert_eq!(s, "s");
        assert_eq!(model, "J2N3819");
        assert!(params.is_empty());
    } else {
        panic!()
    }
}

// ---------------------------------------------------------------------------
// Controlled sources and coupling
// ---------------------------------------------------------------------------

#[test]
fn vcvs_e() {
    let n = parse("T\nE1 out+ out- in+ in- 100\n.end");
    if let ElementKind::Vcvs {
        out_pos,
        out_neg,
        in_pos,
        in_neg,
        gain,
    } = &n.elements().next().unwrap().kind
    {
        assert_eq!(out_pos, "out+");
        assert_eq!(out_neg, "out-");
        assert_eq!(in_pos, "in+");
        assert_eq!(in_neg, "in-");
        assert_num(gain, 100.0);
    } else {
        panic!()
    }
}

#[test]
fn cccs_f() {
    let n = parse("T\nF1 out+ out- Vsense 10\n.end");
    if let ElementKind::Cccs { vsrc, gain, .. } = &n.elements().next().unwrap().kind {
        assert_eq!(vsrc, "Vsense");
        assert_num(gain, 10.0);
    } else {
        panic!()
    }
}

#[test]
fn vccs_g() {
    let n = parse("T\nG1 o+ o- i+ i- 0.02\n.end");
    if let ElementKind::Vccs { gm, .. } = &n.elements().next().unwrap().kind {
        assert_num(gm, 0.02);
    } else {
        panic!()
    }
}

#[test]
fn ccvs_h() {
    let n = parse("T\nH1 o+ o- Vcs 1k\n.end");
    if let ElementKind::Ccvs { vsrc, rm, .. } = &n.elements().next().unwrap().kind {
        assert_eq!(vsrc, "Vcs");
        assert_num(rm, 1e3);
    } else {
        panic!()
    }
}

#[test]
fn mutual_inductance_k() {
    let n = parse("T\nK1 L1 L2 0.99\n.end");
    if let ElementKind::MutualCoupling { l1, l2, coupling } = &n.elements().next().unwrap().kind {
        assert_eq!(l1, "L1");
        assert_eq!(l2, "L2");
        assert_num(coupling, 0.99);
    } else {
        panic!()
    }
}

#[test]
fn behavioural_source_b() {
    let n = parse("T\nB1 out 0 V={sin(2*pi*1k*time)}\n.end");
    if let ElementKind::BehavioralSource { pos, neg, spec } = &n.elements().next().unwrap().kind {
        assert_eq!(pos, "out");
        assert_eq!(neg, "0");
        assert!(spec.starts_with("V="));
    } else {
        panic!()
    }
}

// ---------------------------------------------------------------------------
// Subcircuit call (X element)
// ---------------------------------------------------------------------------

#[test]
fn x_basic() {
    let n = parse("T\nX1 in out gnd OPAMP\n.end");
    if let ElementKind::SubcktCall {
        ports,
        subckt,
        params,
    } = &n.elements().next().unwrap().kind
    {
        assert_eq!(subckt, "OPAMP");
        assert_eq!(ports, &["in", "out", "gnd"]);
        assert!(params.is_empty());
    } else {
        panic!()
    }
}

#[test]
fn x_with_params_keyword() {
    let n = parse("T\nX1 a b FILTER PARAMS: fc=1k Q=0.707\n.end");
    if let ElementKind::SubcktCall { subckt, params, .. } = &n.elements().next().unwrap().kind {
        assert_eq!(subckt, "FILTER");
        assert_eq!(params.len(), 2);
        assert_eq!(params[0].name, "fc");
        assert_num(&params[0].value, 1e3);
        assert_eq!(params[1].name, "Q");
    } else {
        panic!()
    }
}

#[test]
fn x_with_inline_params() {
    // No PARAMS: keyword — key=val tokens recognised directly
    let n = parse("T\nX1 in out MYMOD gain=10 offset=0.5\n.end");
    if let ElementKind::SubcktCall { subckt, params, .. } = &n.elements().next().unwrap().kind {
        assert_eq!(subckt, "MYMOD");
        assert_eq!(params.len(), 2);
    } else {
        panic!()
    }
}

// ---------------------------------------------------------------------------
// Analysis commands
// ---------------------------------------------------------------------------

#[test]
fn analysis_op() {
    let n = parse("T\n.op\n.end");
    assert!(matches!(n.items[0], Item::Analysis(Analysis::Op)));
}

#[test]
fn analysis_dc() {
    let n = parse("T\n.dc Vin -5 5 0.1\n.end");
    if let Item::Analysis(Analysis::Dc {
        src,
        start,
        stop,
        step,
        src2,
    }) = &n.items[0]
    {
        assert_eq!(src, "Vin");
        assert_num(start, -5.0);
        assert_num(stop, 5.0);
        assert_num(step, 0.1);
        assert!(src2.is_none());
    } else {
        panic!()
    }
}

#[test]
fn analysis_dc_double_sweep() {
    let n = parse("T\n.dc V1 0 5 1 V2 0 3 0.5\n.end");
    if let Item::Analysis(Analysis::Dc { src, src2, .. }) = &n.items[0] {
        assert_eq!(src, "V1");
        let s2 = src2.as_ref().expect("should have src2");
        assert_eq!(s2.src, "V2");
        assert_num(&s2.start, 0.0);
        assert_num(&s2.stop, 3.0);
        assert_num(&s2.step, 0.5);
    } else {
        panic!()
    }
}

#[test]
fn analysis_tran() {
    let n = parse("T\n.tran 1n 100n\n.end");
    if let Item::Analysis(Analysis::Tran {
        tstep,
        tstop,
        tstart,
        tmax,
    }) = &n.items[0]
    {
        assert_num(tstep, 1e-9);
        assert_num(tstop, 100e-9);
        assert!(tstart.is_none());
        assert!(tmax.is_none());
    } else {
        panic!()
    }
}

#[test]
fn analysis_tran_with_tstart() {
    let n = parse("T\n.tran 1n 100n 10n\n.end");
    if let Item::Analysis(Analysis::Tran { tstart, tmax, .. }) = &n.items[0] {
        assert_num(tstart.as_ref().unwrap(), 10e-9);
        assert!(tmax.is_none());
    } else {
        panic!()
    }
}

#[test]
fn analysis_ac_dec() {
    let n = parse("T\n.ac DEC 100 1 100Meg\n.end");
    if let Item::Analysis(Analysis::Ac {
        variation,
        n: pts,
        fstart,
        fstop,
    }) = &n.items[0]
    {
        assert_eq!(*variation, AcVariation::Dec);
        assert_eq!(*pts, 100);
        assert_num(fstart, 1.0);
        assert_num(fstop, 100e6);
    } else {
        panic!()
    }
}

#[test]
fn analysis_ac_oct() {
    let n = parse("T\n.ac OCT 8 100 10k\n.end");
    if let Item::Analysis(Analysis::Ac { variation, .. }) = &n.items[0] {
        assert_eq!(*variation, AcVariation::Oct);
    } else {
        panic!()
    }
}

#[test]
fn analysis_ac_lin() {
    let n = parse("T\n.ac LIN 1000 0 1Meg\n.end");
    if let Item::Analysis(Analysis::Ac { variation, .. }) = &n.items[0] {
        assert_eq!(*variation, AcVariation::Lin);
    } else {
        panic!()
    }
}

#[test]
fn analysis_noise() {
    let n = parse("T\n.noise V(out) Vin DEC 10 1 1Meg\n.end");
    if let Item::Analysis(Analysis::Noise {
        output,
        ref_node,
        src,
        variation,
        n: pts,
        ..
    }) = &n.items[0]
    {
        assert_eq!(output, "V(out)");
        assert!(ref_node.is_none());
        assert_eq!(src, "Vin");
        assert_eq!(*variation, AcVariation::Dec);
        assert_eq!(*pts, 10);
    } else {
        panic!()
    }
}

#[test]
fn analysis_tf() {
    let n = parse("T\n.tf V(out) Vin\n.end");
    if let Item::Analysis(Analysis::Tf { output, input }) = &n.items[0] {
        assert_eq!(output, "V(out)");
        assert_eq!(input, "Vin");
    } else {
        panic!()
    }
}

#[test]
fn analysis_sens() {
    let n = parse("T\n.sens V(out)\n.end");
    if let Item::Analysis(Analysis::Sens { output }) = &n.items[0] {
        assert_eq!(output, &["V(out)"]);
    } else {
        panic!()
    }
}

// ---------------------------------------------------------------------------
// Model definitions
// ---------------------------------------------------------------------------

#[test]
fn model_diode() {
    let n = parse("T\n.model 1N4148 D IS=2.52e-9 RS=0.568 N=1.752\n.end");
    if let Item::Model(m) = &n.items[0] {
        assert_eq!(m.name, "1N4148");
        assert_eq!(m.kind, "D");
        assert_eq!(m.params.len(), 3);
        assert_eq!(m.params[0].name, "IS");
    } else {
        panic!()
    }
}

#[test]
fn model_nmos() {
    let n = parse("T\n.model NMOD NMOS VTO=1.5 KP=2e-5\n.end");
    if let Item::Model(m) = &n.items[0] {
        assert_eq!(m.kind, "NMOS");
        assert_eq!(m.params.len(), 2);
    } else {
        panic!()
    }
}

#[test]
fn model_npn_bjt() {
    let n = parse("T\n.model BC547 NPN BF=200 VAF=100\n.end");
    if let Item::Model(m) = &n.items[0] {
        assert_eq!(m.kind, "NPN");
    } else {
        panic!()
    }
}

// ---------------------------------------------------------------------------
// Subcircuit definitions
// ---------------------------------------------------------------------------

#[test]
fn subckt_simple() {
    let n = parse("T\n.subckt BUF in out\nR1 in out 0\n.ends BUF\n.end");
    if let Item::Subckt(s) = &n.items[0] {
        assert_eq!(s.name, "BUF");
        assert_eq!(s.ports, ["in", "out"]);
        assert_eq!(s.items.len(), 1);
        assert!(s.params.is_empty());
    } else {
        panic!()
    }
}

#[test]
fn subckt_with_params() {
    let n =
        parse("T\n.subckt FILTER in out PARAMS: fc=1k Q=0.707\nR1 in out 1k\n.ends FILTER\n.end");
    if let Item::Subckt(s) = &n.items[0] {
        assert_eq!(s.name, "FILTER");
        assert_eq!(s.params.len(), 2);
        assert_eq!(s.params[0].name, "fc");
        assert_num(&s.params[0].value, 1e3);
    } else {
        panic!()
    }
}

#[test]
fn subckt_nested() {
    let src = "T\n.subckt OUTER a b\n.subckt INNER x y\nR1 x y 1k\n.ends INNER\nX1 a b INNER\n.ends OUTER\n.end";
    let n = parse(src);
    if let Item::Subckt(outer) = &n.items[0] {
        assert_eq!(outer.name, "OUTER");
        // First item of OUTER should be a nested subckt
        assert!(outer.items.iter().any(|i| matches!(i, Item::Subckt(_))));
    } else {
        panic!()
    }
}

// ---------------------------------------------------------------------------
// Dot commands
// ---------------------------------------------------------------------------

#[test]
fn dot_param() {
    let n = parse("T\n.param Rval=4.7k Cval=100n\n.end");
    if let Item::Param(ps) = &n.items[0] {
        assert_eq!(ps.len(), 2);
        assert_eq!(ps[0].name, "Rval");
        assert_num(&ps[0].value, 4.7e3);
    } else {
        panic!()
    }
}

#[test]
fn dot_param_brace_expr() {
    let n = parse("T\n.param f0=1k Cval={1/(2*pi*f0*50)}\n.end");
    if let Item::Param(ps) = &n.items[0] {
        assert_eq!(ps.len(), 2);
        assert!(matches!(&ps[1].value, Expr::Brace(_)));
    } else {
        panic!()
    }
}

#[test]
fn dot_include() {
    let n = parse("T\n.include \"models/diodes.lib\"\n.end");
    if let Item::Include(path) = &n.items[0] {
        assert_eq!(path, "models/diodes.lib");
    } else {
        panic!()
    }
}

#[test]
fn dot_lib() {
    let n = parse("T\n.lib \"spice.lib\" nom\n.end");
    if let Item::Lib { file, entry } = &n.items[0] {
        assert_eq!(file, "spice.lib");
        assert_eq!(entry.as_deref(), Some("nom"));
    } else {
        panic!()
    }
}

#[test]
fn dot_global() {
    let n = parse("T\n.global vcc vdd gnd\n.end");
    if let Item::Global(nodes) = &n.items[0] {
        assert_eq!(nodes, &["vcc", "vdd", "gnd"]);
    } else {
        panic!()
    }
}

#[test]
fn dot_options() {
    let n = parse("T\n.options ABSTOL=1e-12 RELTOL=1e-6\n.end");
    if let Item::Options(opts) = &n.items[0] {
        assert_eq!(opts.len(), 2);
        assert_eq!(opts[0].name, "ABSTOL");
    } else {
        panic!()
    }
}

#[test]
fn dot_save() {
    let n = parse("T\n.save V(out) I(R1)\n.end");
    if let Item::Save(vecs) = &n.items[0] {
        assert_eq!(vecs, &["V(out)", "I(R1)"]);
    } else {
        panic!()
    }
}

// ---------------------------------------------------------------------------
// Expr variants
// ---------------------------------------------------------------------------

#[test]
fn param_expr_in_element() {
    let n = parse("T\nR1 a b Rval\n.end");
    if let ElementKind::Resistor { value, .. } = &n.elements().next().unwrap().kind {
        assert!(matches!(value, Expr::Param(s) if s == "Rval"));
    }
}

#[test]
fn brace_expr_in_element() {
    let n = parse("T\nR1 a b {2*Rval}\n.end");
    if let ElementKind::Resistor { value, .. } = &n.elements().next().unwrap().kind {
        assert!(matches!(value, Expr::Brace(s) if s == "2*Rval"));
    }
}

// ---------------------------------------------------------------------------
// Round-trip tests
// ---------------------------------------------------------------------------

#[test]
fn roundtrip_voltage_divider() {
    let src = "Voltage divider\nV1 in 0 DC 5\nR1 in mid 1k\nR2 mid 0 1k\n.op\n.end";
    let (n1, n2) = roundtrip(src);
    assert_eq!(n1.elements().count(), n2.elements().count());
    assert!(matches!(
        &n2.items.last().unwrap(),
        Item::Analysis(Analysis::Op)
    ));
}

#[test]
fn roundtrip_rc_transient() {
    let src = "RC step\nV1 in 0 PULSE(0 5 0 1n 1n 5n 10n)\nR1 in out 1k\nC1 out 0 100n\n.tran 0.1n 50n\n.end";
    let (_n1, n2) = roundtrip(src);
    // Pulse waveform survives round-trip
    if let ElementKind::VoltageSource { source, .. } = &n2.elements().next().unwrap().kind {
        assert!(matches!(source.waveform, Some(Waveform::Pulse { .. })));
    }
}

#[test]
fn roundtrip_mosfet_amp() {
    let src = "NMOS amp\n.model NMOD NMOS VTO=1 KP=2e-5\nM1 drain gate source 0 NMOD W=20u L=2u\nV1 drain 0 DC 5\nVgs gate 0 DC 2\n.op\n.end";
    let (n1, n2) = roundtrip(src);
    assert_eq!(n1.elements().count(), n2.elements().count());
}

#[test]
fn roundtrip_subckt_definition() {
    let src = "Subckt test\n.subckt VOLTDIV in out gnd\nR1 in out 1k\nR2 out gnd 1k\n.ends VOLTDIV\n.op\n.end";
    let (_n1, n2) = roundtrip(src);
    assert!(matches!(&n2.items[0], Item::Subckt(_)));
}

#[test]
fn roundtrip_pwl_source() {
    let src = "PWL test\nV1 a 0 PWL(0 0 1n 5 2n 5 3n 0)\n.tran 0.1n 4n\n.end";
    let (_n1, n2) = roundtrip(src);
    if let ElementKind::VoltageSource { source, .. } = &n2.elements().next().unwrap().kind {
        if let Some(Waveform::Pwl(pts)) = &source.waveform {
            assert_eq!(pts.len(), 4);
        } else {
            panic!()
        }
    }
}

#[test]
fn roundtrip_ac_analysis() {
    let src = "AC sweep\nV1 in 0 AC 1\nR1 in out 1k\nC1 out 0 100n\n.ac DEC 100 1 10Meg\n.end";
    let (_n1, n2) = roundtrip(src);
    assert!(matches!(
        &n2.items.last().unwrap(),
        Item::Analysis(Analysis::Ac {
            variation: AcVariation::Dec,
            ..
        })
    ));
}

// ---------------------------------------------------------------------------
// Error cases
// ---------------------------------------------------------------------------

#[test]
fn empty_input_is_error() {
    assert!(Netlist::parse("").is_err());
    assert!(Netlist::parse("   \n  ").is_err());
}

#[test]
fn missing_resistor_value_is_error() {
    assert!(Netlist::parse("T\nR1 a b\n.end").is_err());
}

#[test]
fn missing_mosfet_bulk_is_error() {
    // M needs d g s bulk model — only 3 nodes provided
    assert!(Netlist::parse("T\nM1 d g s NMOD\n.end").is_err());
}

// ---------------------------------------------------------------------------
// Netlist::elements / elements_mut
// ---------------------------------------------------------------------------

#[test]
fn elements_iter_skips_non_elements() {
    let n = parse("T\n* comment\n.param R=1k\nR1 a b 1k\nC1 a b 100n\n.op\n.end");
    assert_eq!(n.elements().count(), 2);
}

#[test]
fn elements_mut_rename() {
    let mut n = parse("T\nR1 a b 1k\n.end");
    for e in n.elements_mut() {
        e.name = "R99".into();
    }
    let rendered = n.to_string();
    assert!(rendered.contains("R99"));
}

// ---------------------------------------------------------------------------
// to_lines
// ---------------------------------------------------------------------------

#[test]
fn to_lines_ends_with_dot_end() {
    let n = parse("T\nR1 a b 1k\n.op\n.end");
    let lines = n.to_lines();
    assert_eq!(lines.last().unwrap(), ".end");
    assert!(lines[0].contains("T"));
}

// ---------------------------------------------------------------------------
// Display / generation
// ---------------------------------------------------------------------------

#[test]
fn display_dc_sweep_with_src2() {
    let a = Analysis::Dc {
        src: "V1".into(),
        start: Expr::Num(0.0),
        stop: Expr::Num(5.0),
        step: Expr::Num(0.1),
        src2: Some(DcSweep {
            src: "V2".into(),
            start: Expr::Num(0.0),
            stop: Expr::Num(3.0),
            step: Expr::Num(0.5),
        }),
    };
    let s = format!("{a}");
    assert!(s.starts_with(".dc V1"), "got: {s}");
    assert!(s.contains("V2"), "got: {s}");
}

#[test]
fn display_source_dc_and_ac() {
    let s = Source {
        dc: Some(Expr::Num(0.0)),
        ac: Some(AcSpec {
            mag: Expr::Num(1.0),
            phase: Some(Expr::Num(45.0)),
        }),
        waveform: None,
    };
    let out = s.to_string();
    assert!(out.contains("DC 0"), "got: {out}");
    assert!(out.contains("AC 1 45"), "got: {out}");
}
