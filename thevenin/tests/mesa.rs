#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use approx::assert_abs_diff_eq;
use thevenin::simulate_dc;
use thevenin_types::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// Port of ngspice-upstream/tests/mesa/mesa.cir
///
/// DCFL GaAs MESFET inverter gate. DC sweep vin 0→3V.
/// Two MESA FETs: enhancement (vto=0.1) and depletion (vto=-1.0).
#[test]
fn test_mesa_dcfl_inverter() {
    let cir = include_str!("fixtures/mesa/mesa.cir");
    let netlist = Netlist::parse(cir).unwrap();
    let result = simulate_dc(&netlist).unwrap();

    let plot = &result.plots[0];
    let v2 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(2)")
        .expect("no v(2)");

    // 61 sweep points (0 to 3.0 in steps of 0.05)
    assert_eq!(
        v2.real.len(),
        61,
        "expected 61 sweep points, got {}",
        v2.real.len()
    );

    // Reference from ngspice: at vin=0 (index 0), v(2) ≈ 2.992148
    assert_abs_diff_eq!(v2.real[0], 2.992148, epsilon = 0.05);

    // At vin=1.0 (index 20), v(2) ≈ 1.945616
    assert_abs_diff_eq!(v2.real[20], 1.945616, epsilon = 0.05);

    // At vin=2.0 (index 40), v(2) ≈ 0.219048
    assert_abs_diff_eq!(v2.real[40], 0.219048, epsilon = 0.05);

    // At vin=3.0 (index 60), v(2) ≈ 0.160660
    assert_abs_diff_eq!(v2.real[60], 0.160660, epsilon = 0.05);
}

/// Port of ngspice-upstream/tests/mesa/mesa13.cir
///
/// Gate leakage test with level=2 MESA FET.
#[test]
fn test_mesa_gate_leakage() {
    let cir = include_str!("fixtures/mesa/mesa13.cir");
    let netlist = Netlist::parse(cir).unwrap();
    let result = simulate_dc(&netlist).unwrap();
    assert!(!result.plots.is_empty());
}

/// Port of ngspice-upstream/tests/mesa/mesa11.cir
///
/// Level 2 DC double sweep (vds 0→2, vgs -1.2→0).
#[test]
fn test_mesa_level2_double_sweep() {
    let cir = include_str!("fixtures/mesa/mesa11.cir");
    let netlist = Netlist::parse(cir).unwrap();
    let result = simulate_dc(&netlist).unwrap();
    assert!(!result.plots.is_empty());
}
