#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use approx::assert_abs_diff_eq;
use thevenin::simulate_dc;
use thevenin_types::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// Port of ngspice-upstream/tests/mes/subth.cir
///
/// MESFET level=1 subthreshold characteristics.
/// DC sweep VGS from -3V to 0V, measures drain current through vids.
#[test]
fn test_mesfet_subthreshold() {
    let cir = include_str!("fixtures/mes/subth.cir");
    let netlist = Netlist::parse(cir).unwrap();
    let result = simulate_dc(&netlist).unwrap();

    let plot = &result.plots[0];
    let ibranch = plot
        .vecs
        .iter()
        .find(|v| v.name == "vids#branch")
        .expect("no vids#branch");

    for (i, val) in ibranch.real.iter().enumerate() {
        if i % 10 == 0 || (i >= 33 && i <= 37) || i == 60 {
            eprintln!("  [{:3}] I = {:12.6e}", i, val);
        }
    }
    // 61 sweep points (-3 to 0 in steps of 0.05)
    assert_eq!(ibranch.real.len(), 61, "expected 61 sweep points");

    // In subthreshold (vgs << vt0=-1.3): current is very small (~pA)
    // At index 0 (vgs=-3.0): ~3.114e-12
    assert!(
        ibranch.real[0].abs() < 1e-9,
        "subthreshold current should be near zero"
    );

    // Near threshold (vgs=-1.25, index 35): ~1.31e-6
    assert_abs_diff_eq!(ibranch.real[35], 1.308946e-6, epsilon = 1e-6);

    // Above threshold (vgs=0, index 60): ~4.575e-4
    assert_abs_diff_eq!(ibranch.real[60], 4.575371e-4, epsilon = 5e-5);
}
