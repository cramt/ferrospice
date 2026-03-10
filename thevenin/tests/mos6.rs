//! Integration tests for MOS Level 6 (Sakurai-Newton) model.
//!
//! Ported from ngspice-upstream/tests/mos6/simpleinv.cir.

#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use thevenin_types::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

const SIMPLEINV_CIR: &str = include_str!("fixtures/mos6/simpleinv.cir");

/// Test that the MOS6 CMOS inverter netlist parses and the DC operating point
/// converges. At VIN=0 (DC value), PMOS is on and NMOS is off, so V(11) should
/// be near VDD (5V) or near 0V depending on the model bias conditions.
#[test]
fn test_mos6_simpleinv_parses_and_op_converges() {
    let netlist = Netlist::parse(SIMPLEINV_CIR).unwrap();

    // Run the operating point (DC) analysis to verify convergence.
    let result = thevenin::simulate_op(&netlist).unwrap();
    assert!(!result.plots.is_empty(), "should produce at least one plot");

    let plot = &result.plots[0];
    // Check that VDD node (100) is at 5V.
    let v100 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(100)")
        .expect("should have v(100)");
    assert!(
        (v100.real[0] - 5.0).abs() < 0.01,
        "VDD should be 5V, got {}",
        v100.real[0]
    );
}

/// Test transient analysis of MOS6 CMOS inverter with PWL input ramp.
#[test]
fn test_mos6_simpleinv_transient() {
    let netlist = Netlist::parse(SIMPLEINV_CIR).unwrap();
    let result = thevenin::simulate_tran(&netlist).unwrap();

    assert!(!result.plots.is_empty(), "should produce at least one plot");
    let plot = &result.plots[0];
    assert_eq!(plot.name, "tran1");

    // Find v(1) (input) and v(11) (output).
    let v1 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(1)")
        .expect("should have v(1)");
    let v11 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(11)")
        .expect("should have v(11)");

    assert!(!v1.real.is_empty(), "v(1) should have data points");
    assert!(!v11.real.is_empty(), "v(11) should have data points");

    // Input starts at 0 and ramps to 5V.
    assert!(
        v1.real[0].abs() < 0.01,
        "v(1) should start near 0, got {}",
        v1.real[0]
    );
    let last_v1 = *v1.real.last().unwrap();
    assert!(last_v1 > 4.5, "v(1) should end near 5V, got {last_v1}");
}
