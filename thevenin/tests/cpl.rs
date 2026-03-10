use thevenin_types::Netlist;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// Helper to get a vector from a transient result.
fn tran_vector<'a>(result: &'a thevenin_types::SimResult, name: &str) -> &'a Vec<f64> {
    let plot = &result.plots[0];
    &plot
        .vecs
        .iter()
        .find(|v| v.name == name)
        .unwrap_or_else(|| panic!("no vector '{name}'"))
        .real
}

/// Test CPL DC operating point: 2-line CPL with R*length series resistance.
#[test]
fn test_cpl_dc_op() {
    let netlist = Netlist::parse(
        "CPL DC test
V1 1 0 DC 5
R1 1 2 100
P1 2 3 0 4 5 0 cmod
R2 4 0 100
R3 5 0 100
.model cmod cpl R=0.5 0 0.5 L=1e-9 0 1e-9 G=0 0 0 C=1e-12 0 1e-12 length=1
.op
.end
",
    )
    .expect("parse failed");

    let result = thevenin::simulate_op(&netlist).expect("op failed");
    let plot = &result.plots[0];

    // DC: CPL line 1 acts as R[0][0]*length = 0.5*1 = 0.5 ohm resistor
    // V1(5V) -> R1(100) -> CPL_line1(0.5) -> R2(100) -> GND
    // I = 5 / (100 + 0.5 + 100) = 5/200.5 ≈ 0.02494A
    // V(4) = I * R2 ≈ 2.494V
    let v4_vec = plot.vecs.iter().find(|v| v.name == "v(4)").unwrap();
    let v4 = v4_vec.real[0];
    let expected = 5.0 * 100.0 / 200.5;
    assert!(
        (v4 - expected).abs() < 0.05,
        "v(4) DC OP should be ~{expected:.3}, got {v4:.3}"
    );
}

/// IBM2 coupled transmission line transient test.
#[test]
fn test_cpl_ibm2_transient() {
    let cir = include_str!("fixtures/transmission/cpl_ibm2.cir");
    let netlist = Netlist::parse(cir).expect("parse failed");
    let result = thevenin::simulate_tran(&netlist).expect("tran failed");
    let time = tran_vector(&result, "time");

    // The simulation should complete without errors.
    assert!(
        time.len() > 50,
        "should have >50 timepoints, got {}",
        time.len()
    );
}

/// 4-line CPL transient test with R load.
#[test]
fn test_cpl3_4_line_transient() {
    let cir = include_str!("fixtures/transmission/cpl3_4_line.cir");
    let netlist = Netlist::parse(cir).expect("parse failed");
    let result = thevenin::simulate_tran(&netlist).expect("tran failed");
    let time = tran_vector(&result, "time");

    assert!(
        time.len() > 50,
        "should have >50 timepoints, got {}",
        time.len()
    );
}
