#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use ferrospice_core::simulate_op;
use ferrospice_netlist::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// Basic BSIM3 NMOS operating point test.
///
/// Simple NMOS with LEVEL=8 in saturation.
/// Validates that the simulator produces reasonable current.
#[test]
fn test_bsim3_nmos_op() {
    let cir = "\
BSIM3 NMOS OP Test
.model nmod nmos level=8 version=3.3 tox=1.5e-8 vth0=0.7 u0=670 k1=0.5 k2=0.0 vsat=8e4
M1 d g 0 0 nmod W=10e-6 L=1e-6
Vgs g 0 1.8
Vds d 0 1.8
.op
.end
";
    let netlist = Netlist::parse(cir).unwrap();
    let result = simulate_op(&netlist).unwrap();

    let plot = &result.plots[0];

    // Get drain current from Vds branch
    let i_vds = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");
    let ids = i_vds.real[0];

    // BSIM3 NMOS in saturation: Vgs=1.8V, Vds=1.8V, Vth~0.7V
    // Should have negative branch current (current flows into drain)
    // Expect roughly 100μA to 10mA range for W=10μ, L=1μ
    assert!(
        ids < 0.0,
        "Drain current should be negative (conventional): {}",
        ids
    );
    assert!(
        ids.abs() > 1e-6,
        "Should have significant drain current: {}",
        ids
    );
    assert!(
        ids.abs() < 0.1,
        "Current should be reasonable (< 100mA): {}",
        ids
    );
}

/// Basic BSIM3 NMOS DC sweep test.
///
/// Sweeps VDS from 0 to 2V at fixed VGS=1.8V.
/// Validates monotonic increasing drain current.
#[test]
fn test_bsim3_nmos_dc_sweep() {
    let cir = "\
BSIM3 NMOS DC Sweep Test
.model nmod nmos level=8 version=3.3 tox=1.5e-8 vth0=0.7 u0=670 k1=0.5 k2=0.0 vsat=8e4
M1 d g 0 0 nmod W=10e-6 L=1e-6
Vgs g 0 1.8
Vds d 0 0.0
.dc Vds 0.0 2.0 0.2
.end
";
    let netlist = Netlist::parse(cir).unwrap();
    let result = ferrospice_core::simulate_dc(&netlist).unwrap();

    let plot = &result.plots[0];

    let sweep_var = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds")
        .expect("no vds sweep variable");

    let i_vds = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");

    // Should have multiple sweep points
    assert!(
        sweep_var.real.len() >= 10,
        "Expected at least 10 sweep points"
    );

    // Current should increase (become more negative) as VDS increases from 0
    let first_current = i_vds.real[0].abs();
    let last_current = i_vds.real[i_vds.real.len() - 1].abs();

    // At VDS=0, current should be ~0; at VDS=2V (saturation), current should be significant
    assert!(
        last_current > first_current,
        "Current should increase with VDS: first={}, last={}",
        first_current,
        last_current
    );
}

/// BSIM3 PMOS operating point test.
#[test]
fn test_bsim3_pmos_op() {
    let cir = "\
BSIM3 PMOS OP Test
.model pmod pmos level=8 version=3.3 tox=1.5e-8 vth0=-0.7 u0=250 k1=0.5 k2=0.0 vsat=8e4
M1 d g vdd vdd pmod W=10e-6 L=1e-6
Vdd vdd 0 1.8
Vgs g 0 0.0
Vds d 0 0.0
.op
.end
";
    let netlist = Netlist::parse(cir).unwrap();
    let result = simulate_op(&netlist).unwrap();

    let plot = &result.plots[0];

    // Get drain current from Vds branch
    let i_vds = plot
        .vecs
        .iter()
        .find(|v| v.name == "vds#branch")
        .expect("no vds#branch");
    let ids = i_vds.real[0];

    // PMOS: VGS = 0 - 1.8 = -1.8V, VDS = 0 - 1.8 = -1.8V
    // Current should flow out of drain (positive branch current)
    assert!(
        ids > 0.0,
        "PMOS drain current should be positive (conventional): {}",
        ids
    );
    assert!(
        ids.abs() > 1e-6,
        "Should have significant drain current: {}",
        ids
    );
}
