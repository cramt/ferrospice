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

/// Test TXL DC operating point.
/// A TXL line with R > 0 should act as R*length resistance in DC.
#[test]
fn test_txl_dc_op() {
    let netlist = Netlist::parse(
        "TXL DC test
V1 1 0 DC 5
R1 1 2 100
Y1 2 0 3 0 ymod
R2 3 0 100
.model ymod txl R=10 L=1e-9 G=0 C=1e-12 length=1
.op
.end
",
    )
    .expect("parse failed");

    let result = thevenin::simulate_op(&netlist).expect("op failed");
    let plot = &result.plots[0];

    // DC: TXL acts as R*length = 10*1 = 10 ohm resistor
    // V1(5V) -> R1(100) -> TXL(10) -> R2(100) -> GND
    // I = 5 / (100 + 10 + 100) = 5/210 ≈ 0.0238A
    // V(3) = I * R2 = 0.0238 * 100 ≈ 2.381V
    let v3_vec = plot.vecs.iter().find(|v| v.name == "v(3)").unwrap();
    let v3 = v3_vec.real[0];
    let expected = 5.0 * 100.0 / 210.0;
    assert!(
        (v3 - expected).abs() < 0.01,
        "v(3) DC OP should be ~{expected:.3}, got {v3:.3}"
    );
}

/// Simple TXL transient test with a step source.
#[test]
fn test_txl_transient_simple() {
    let netlist = Netlist::parse(
        "TXL transient test
V1 1 0 PULSE(0 5 0 0.1n 0.1n 50n 100n)
R1 1 2 50
Y1 2 0 3 0 ymod
R2 3 0 200
.model ymod txl R=12.45 L=8.972e-9 G=0 C=0.468e-12 length=16
.tran 0.2n 40n
.end
",
    )
    .expect("parse failed");

    let result = thevenin::simulate_tran(&netlist).expect("tran failed");
    let time = tran_vector(&result, "time");
    let v3 = tran_vector(&result, "v(3)");

    // td = sqrt(LC)*length = sqrt(8.972e-9 * 0.468e-12) * 16 ≈ 1.04ns
    // Signal should appear at port 2 after ~1ns delay

    // Before delay: v(3) should be ~0
    let idx_early = time.iter().position(|&t| t >= 0.5e-9).unwrap();
    assert!(
        v3[idx_early].abs() < 0.1,
        "v(3) should be small before td, got {}",
        v3[idx_early]
    );

    // After signal settles: v(3) should show significant voltage
    let idx_late = time.iter().position(|&t| t >= 10e-9).unwrap();
    assert!(
        v3[idx_late] > 0.5,
        "v(3) should show signal after settling, got {} at t={}",
        v3[idx_late],
        time[idx_late]
    );
}

/// Full ngspice txl1_1_line test (MOS driver + TXL line + C load).
#[test]
fn test_txl1_1_line_ngspice() {
    let cir = include_str!("fixtures/transmission/txl1_1_line.cir");
    let netlist = Netlist::parse(cir).expect("parse failed");
    let result = thevenin::simulate_tran(&netlist).expect("tran failed");
    let time = tran_vector(&result, "time");
    let _v2 = tran_vector(&result, "v(2)");
    let _v3 = tran_vector(&result, "v(3)");

    // The simulation should complete without errors.
    assert!(
        time.len() > 50,
        "should have >50 timepoints, got {}",
        time.len()
    );
}

/// Full ngspice txl2_3_line test (3 TXL lines in cascade with MOSFETs).
#[test]
fn test_txl2_3_line_ngspice() {
    let cir = include_str!("fixtures/transmission/txl2_3_line.cir");
    let netlist = Netlist::parse(cir).expect("parse failed");
    let result = thevenin::simulate_tran(&netlist).expect("tran failed");
    let time = tran_vector(&result, "time");

    assert!(
        time.len() > 50,
        "should have >50 timepoints, got {}",
        time.len()
    );
}
