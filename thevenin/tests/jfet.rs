#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use approx::assert_abs_diff_eq;
use thevenin::simulate_dc;
use thevenin_types::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// Port of ngspice-upstream/tests/jfet/jfet_vds-vgs.cir.
///
/// Tests N-channel JFET 2N4221 DC sweep: VDS 0→25V, VGS -3→0V.
/// Validates drain current against ngspice reference output.
#[test]
fn test_jfet_2n4221_dc_sweep() {
    let cir = include_str!("fixtures/jfet/jfet_vds-vgs.cir");
    let netlist = Netlist::parse(cir).unwrap();
    let result = simulate_dc(&netlist).unwrap();

    let plot = &result.plots[0];
    let sweep_var = plot
        .vecs
        .iter()
        .find(|v| v.name == "v-sweep")
        .expect("no v-sweep variable");
    let i_vd = plot
        .vecs
        .iter()
        .find(|v| v.name == "vd#branch")
        .expect("no vd#branch");

    // Total expected points: 26 per VGS step × 4 VGS steps = 104
    assert_eq!(
        sweep_var.real.len(),
        104,
        "expected 104 sweep points, got {}",
        sweep_var.real.len()
    );

    // Reference data from ngspice-upstream/tests/jfet/jfet_vds-vgs.out
    // VGS=-3V (cutoff region: VTO=-3.5, so Vgst=0.5)
    // At VDS=0: ~0
    // At VDS=1: ~-1.027e-4 (saturation since vgst=0.5 < vds=1)
    // VGS=-3, VDS=25: ~-1.076e-4

    // First 26 points are VGS=-3
    // Check index 0 (VDS=0)
    assert!(
        i_vd.real[0].abs() < 1e-6,
        "VGS=-3, VDS=0: expected ~0, got {:.6e}",
        i_vd.real[0]
    );

    // Check index 1 (VDS=1, VGS=-3): should be in saturation
    let expected_1 = -1.027008e-4;
    assert_abs_diff_eq!(i_vd.real[1], expected_1, epsilon = 1e-6);

    // Check index 25 (VDS=25, VGS=-3)
    let expected_25 = -1.076206e-4;
    assert_abs_diff_eq!(i_vd.real[25], expected_25, epsilon = 1e-6);

    // Second 26 points are VGS=-2
    // Check index 27 (VDS=1, VGS=-2): should have larger current
    let expected_27 = -7.504957e-4;
    assert_abs_diff_eq!(i_vd.real[27], expected_27, epsilon = 1e-5);

    // VGS=-2, VDS=25: ~-9.682677e-4
    let expected_51 = -9.682677e-4;
    assert_abs_diff_eq!(i_vd.real[51], expected_51, epsilon = 1e-5);

    // Fourth 26 points are VGS=0 (fully on)
    // VGS=0, VDS=5: ~-5.062554e-3
    let expected_83 = -5.062554e-3;
    assert_abs_diff_eq!(i_vd.real[83], expected_83, epsilon = 1e-4);
}
