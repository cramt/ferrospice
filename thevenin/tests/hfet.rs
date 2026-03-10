#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use approx::assert_abs_diff_eq;
use thevenin::simulate_dc;
use thevenin_types::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// Port of ngspice-upstream/tests/hfet/id_vgs.cir
///
/// HFET Id versus Vgs characteristic.
/// DC sweep VDS from 0 to 1V, NHFET level=5.
#[test]
fn test_hfet_id_vgs() {
    let cir = include_str!("fixtures/hfet/id_vgs.cir");
    let netlist = Netlist::parse(cir).unwrap();
    let result = simulate_dc(&netlist).unwrap();

    let plot = &result.plots[0];
    let ids = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");

    // 101 sweep points (0 to 1.0 in steps of 0.01)
    assert_eq!(ids.real.len(), 101, "expected 101 sweep points");

    // At vds=0 (index 0): current ≈ 0 (3.77e-11)
    assert!(
        ids.real[0].abs() < 1e-6,
        "zero-bias current should be near zero"
    );

    // At vds=0.1 (index 10): ~-1.075e-4
    assert_abs_diff_eq!(ids.real[10], -1.075027e-4, epsilon = 2e-5);

    // At vds=0.5 (index 50): ~-1.687e-4
    assert_abs_diff_eq!(ids.real[50], -1.687463e-4, epsilon = 2e-5);

    // At vds=1.0 (index 100): ~-2.108e-4
    assert_abs_diff_eq!(ids.real[100], -2.108001e-4, epsilon = 2e-5);
}

/// Port of ngspice-upstream/tests/hfet/inverter.cir
///
/// DCFL inverter with NHFET subcircuits and transient analysis.
/// This test is ignored until transient analysis for HFET is verified.
#[test]
#[ignore = "HFET inverter requires transient analysis support"]
fn test_hfet_inverter() {
    let cir = include_str!("fixtures/hfet/inverter.cir");
    let _netlist = Netlist::parse(cir).unwrap();
    // TODO: implement transient test once HFET transient support is verified
}
