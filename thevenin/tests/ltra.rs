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

/// Simple lossless transmission line (LC) test.
/// A step source drives one end, the other end is open-circuit.
/// The signal should appear at the far end after td.
#[test]
fn test_ltra_lossless_step() {
    let netlist = Netlist::parse(
        "LTRA lossless test
V1 1 0 PULSE(0 1 0 0.1n 0.1n 100n 200n)
R1 1 2 50
O1 2 0 3 0 tline
R2 3 0 200
.model tline ltra r=0 l=250e-9 c=100e-12 len=1
.tran 0.1n 60n
.end
",
    )
    .expect("parse failed");

    let result = thevenin::simulate_tran(&netlist).expect("tran failed");
    let time = tran_vector(&result, "time");
    let _v2 = tran_vector(&result, "v(2)");
    let v3 = tran_vector(&result, "v(3)");

    // td = sqrt(LC)*len = sqrt(250e-9 * 100e-12) * 1 = 5ns
    // At t < 5ns, v(3) should be ~0
    // At t > 5ns, v(3) should show the delayed signal

    // Find time index near t=2ns (before signal arrives at port 2)
    let idx_2ns = time.iter().position(|&t| t >= 2e-9).unwrap();
    assert!(
        v3[idx_2ns].abs() < 0.01,
        "v(3) should be ~0 before td, got {} at t={}",
        v3[idx_2ns],
        time[idx_2ns]
    );

    // After the delay + some settling, v(3) should show signal
    let idx_20ns = time.iter().position(|&t| t >= 20e-9).unwrap();
    assert!(
        v3[idx_20ns] > 0.1,
        "v(3) should show signal after td, got {} at t={}",
        v3[idx_20ns],
        time[idx_20ns]
    );
}

/// RLC lossy transmission line test.
/// Uses the same model parameters as the ngspice ltra1 test but with a simpler driver.
#[test]
fn test_ltra_rlc_simple() {
    let netlist = Netlist::parse(
        "LTRA RLC simple test
V1 1 0 PULSE(0 5 0 0.2n 0.2n 15n 32n)
R1 1 2 50
O1 2 0 3 0 lline
C1 3 0 0.025e-12
.model lline ltra r=12.45 l=8.972e-9 c=0.468e-12 len=16 rel=1
.tran 0.2n 40n
.end
",
    )
    .expect("parse failed");

    let result = thevenin::simulate_tran(&netlist).expect("tran failed");
    let time = tran_vector(&result, "time");
    let v3 = tran_vector(&result, "v(3)");

    // td = sqrt(LC)*len = sqrt(8.972e-9 * 0.468e-12) * 16 ≈ 1.04ns
    // Signal should appear at port 2 after ~1ns delay

    // Before delay: v(3) should be ~0
    let idx_early = time.iter().position(|&t| t >= 0.5e-9).unwrap();
    assert!(
        v3[idx_early].abs() < 0.1,
        "v(3) should be small before td, got {}",
        v3[idx_early]
    );

    // After signal arrives and settles: v(3) should show significant voltage
    let idx_late = time.iter().position(|&t| t >= 10e-9).unwrap();
    assert!(
        v3[idx_late] > 0.5,
        "v(3) should show signal after settling, got {} at t={}",
        v3[idx_late],
        time[idx_late]
    );
}

/// Test that LTRA DC operating point works correctly.
/// A simple LTRA with R > 0 should act like a resistor in DC.
#[test]
fn test_ltra_dc_op() {
    let netlist = Netlist::parse(
        "LTRA DC test
V1 1 0 DC 5
R1 1 2 100
O1 2 0 3 0 rline
R2 3 0 100
.model rline ltra r=10 l=1e-9 c=1e-12 len=1
.op
.end
",
    )
    .expect("parse failed");

    let result = thevenin::simulate_op(&netlist).expect("op failed");
    let plot = &result.plots[0];

    // DC: LTRA acts as R*len = 10*1 = 10 ohm resistor
    // V1(5V) -> R1(100) -> LTRA(10) -> R2(100) -> GND
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

/// Full ngspice ltra1_1_line test (MOS driver + lossy line + cap load).
/// Compare against reference output at key timepoints.
#[test]
fn test_ltra1_1_line_ngspice() {
    let cir = include_str!("fixtures/transmission/ltra1_1_line.cir");
    let netlist = Netlist::parse(cir).expect("parse failed");
    let result = thevenin::simulate_tran(&netlist).expect("tran failed");
    let time = tran_vector(&result, "time");
    let v2 = tran_vector(&result, "v(2)");
    let v3 = tran_vector(&result, "v(3)");

    // Reference: at t=0, v(2)~5.0, v(3)~5.0 (MOSFET pull-up)
    assert!(
        (v2[0] - 5.0).abs() < 0.5,
        "v(2) at t=0 should be ~5V, got {}",
        v2[0]
    );
    assert!(
        (v3[0] - 5.0).abs() < 0.5,
        "v(3) at t=0 should be ~5V, got {}",
        v3[0]
    );

    // The simulation should complete without errors.
    // Verify we have a reasonable number of output points.
    assert!(
        time.len() > 50,
        "should have >50 timepoints, got {}",
        time.len()
    );
}

/// Full ngspice ltra2_2_line test (MOS driver + 2 lossy lines in cascade).
#[test]
fn test_ltra2_2_line_ngspice() {
    let cir = include_str!("fixtures/transmission/ltra2_2_line.cir");
    let netlist = Netlist::parse(cir).expect("parse failed");
    let result = thevenin::simulate_tran(&netlist).expect("tran failed");
    let time = tran_vector(&result, "time");
    let _v2 = tran_vector(&result, "v(2)");
    let _v3 = tran_vector(&result, "v(3)");
    let _v4 = tran_vector(&result, "v(4)");
    let _v5 = tran_vector(&result, "v(5)");

    // Should complete with reasonable number of points.
    assert!(
        time.len() > 50,
        "should have >50 timepoints, got {}",
        time.len()
    );
}
