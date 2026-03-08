#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use approx::assert_abs_diff_eq;
use ferrospice_core::{simulate_sens, simulate_tf};
use ferrospice_netlist::Netlist;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

fn tf_value(result: &ferrospice_netlist::SimResult, plot_idx: usize, name: &str) -> f64 {
    result.plots[plot_idx]
        .vecs
        .iter()
        .find(|v| v.name == name)
        .unwrap_or_else(|| panic!("no vector '{name}'"))
        .real[0]
}

fn sens_value(result: &ferrospice_netlist::SimResult, name: &str) -> f64 {
    result.plots[0]
        .vecs
        .iter()
        .find(|v| v.name == name)
        .unwrap_or_else(|| {
            let names: Vec<_> = result.plots[0].vecs.iter().map(|v| &v.name).collect();
            panic!("no vector '{name}', available: {names:?}")
        })
        .real[0]
}

/// Port of ngspice-upstream/tests/sensitivity/diffpair.cir — .tf analysis
///
/// Simple differential pair with two .tf commands:
/// .tf v(5) vcm — common-mode transfer function
/// .tf v(5) vdm — differential-mode transfer function
#[test]
fn test_diffpair_tf() {
    let netlist = Netlist::parse(
        "simple differential pair - transfer functions
.model qnl npn(bf=80 rb=100 tf=0.3n tr=6n cje=3p cjc=2p vaf=50)
.model qnr npn(bf=80 rb=100 tf=0.3n tr=6n cje=3p cjc=2p vaf=50)

q1 4 2 6 qnr
q2 5 3 6 qnl
rs1 11 2 1k
rs2 3 1 1k
rc1 4 8 10k
rc2 5 8 10k
q3 7 7 9 qnl
q4 6 7 9 qnr
rbias 7 8 20k

vcm 1 0 DC 0
vdm 1 11 DC 0
vcc 8 0 12
vee 9 0 -12

.tf v(5) vcm
.tf v(5) vdm
.end
",
    )
    .unwrap();

    let result = simulate_tf(&netlist).unwrap();
    assert_eq!(result.plots.len(), 2);

    // .tf v(5) vcm — common-mode TF
    // ngspice reference: transfer_function = -1.10341e-01
    let tf_cm = tf_value(&result, 0, "transfer_function");
    assert_abs_diff_eq!(tf_cm, -1.10341e-01, epsilon = 1e-2);

    // Output impedance at v(5) ≈ 9447 Ω
    let z_out_cm = tf_value(&result, 0, "output_impedance_at_v(5)");
    assert_abs_diff_eq!(z_out_cm, 9.447e3, epsilon = 1e2);

    // Input impedance of vcm ≈ 1.793 MΩ
    let z_in_cm = tf_value(&result, 0, "vcm#input_impedance");
    assert_abs_diff_eq!(z_in_cm, 1.793e6, epsilon = 1e4);

    // .tf v(5) vdm — differential-mode TF
    // ngspice reference: transfer_function = -8.78493e+01
    let tf_dm = tf_value(&result, 1, "transfer_function");
    assert_abs_diff_eq!(tf_dm, -8.78493e1, epsilon = 1.0);

    // Input impedance of vdm ≈ 8941 Ω
    let z_in_dm = tf_value(&result, 1, "vdm#input_impedance");
    assert_abs_diff_eq!(z_in_dm, 8.941e3, epsilon = 1e2);
}

/// Port of ngspice-upstream/tests/sensitivity/diffpair.cir — .sens analysis
///
/// Tests DC sensitivity of v(5,4) to selected component parameters.
#[test]
fn test_diffpair_sens() {
    let netlist = Netlist::parse(
        "simple differential pair - sensitivity
.model qnl npn(bf=80 rb=100 tf=0.3n tr=6n cje=3p cjc=2p vaf=50)
.model qnr npn(bf=80 rb=100 tf=0.3n tr=6n cje=3p cjc=2p vaf=50)

q1 4 2 6 qnr
q2 5 3 6 qnl
rs1 11 2 1k
rs2 3 1 1k
rc1 4 8 10k
rc2 5 8 10k
q3 7 7 9 qnl
q4 6 7 9 qnr
rbias 7 8 20k

vcm 1 0 DC 0
vdm 1 11 DC 0
vcc 8 0 12
vee 9 0 -12

.sens v(5,4)
.end
",
    )
    .unwrap();

    let result = simulate_sens(&netlist).unwrap();

    // Resistor sensitivities (from diffpair.out reference)
    // rc1 sensitivity: 6.031558e-04
    let s_rc1 = sens_value(&result, "rc1");
    assert_abs_diff_eq!(s_rc1, 6.031558e-4, epsilon = 1e-6);

    // rc2 sensitivity: -6.031558e-04
    let s_rc2 = sens_value(&result, "rc2");
    assert_abs_diff_eq!(s_rc2, -6.031558e-4, epsilon = 1e-6);

    // rs1 sensitivity: -1.346892e-03
    let s_rs1 = sens_value(&result, "rs1");
    assert_abs_diff_eq!(s_rs1, -1.346892e-3, epsilon = 1e-6);

    // rs2 sensitivity: 1.346892e-03
    let s_rs2 = sens_value(&result, "rs2");
    assert_abs_diff_eq!(s_rs2, 1.346892e-3, epsilon = 1e-6);

    // Voltage source sensitivity: vdm = -1.758090e+02
    let s_vdm = sens_value(&result, "vdm");
    assert_abs_diff_eq!(s_vdm, -1.758090e2, epsilon = 1.0);
}

/// Port of ngspice-upstream/tests/regression/sens/sens-dc-1
///
/// Simple I*R circuit: I1=42mA, R1=1k → V(1) = 42V
#[test]
fn test_sens_dc_1() {
    let netlist = Netlist::parse(
        "test sens dc
i1 0 1 DC 42m
r1 1 0 1k
.sens v(1) dc
.end
",
    )
    .unwrap();

    let result = simulate_sens(&netlist).unwrap();

    // Golden values from ngspice
    assert_abs_diff_eq!(sens_value(&result, "i1"), 1000.0, epsilon = 1e-3);
    assert_abs_diff_eq!(sens_value(&result, "r1"), 0.042, epsilon = 1e-6);
}

/// Port of ngspice-upstream/tests/regression/sens/sens-dc-2
///
/// Resistor network with voltage source
#[test]
fn test_sens_dc_2() {
    let netlist = Netlist::parse(
        "test sens dc
v1 1 0 DC 42
r1 1 2 1k
r2 2 3 1.5k
r3 3 0 2.2k
r4 3 4 3.3k
r5 4 0 1.8k
rx 2 4 2.7k
.sens v(4) dc
.end
",
    )
    .unwrap();

    let result = simulate_sens(&netlist).unwrap();

    // Golden values from ngspice
    assert_abs_diff_eq!(
        sens_value(&result, "v1"),
        0.2933065931311105,
        epsilon = 1e-6
    );
    assert_abs_diff_eq!(
        sens_value(&result, "r1"),
        -0.004104514242490192,
        epsilon = 1e-9
    );
    assert_abs_diff_eq!(
        sens_value(&result, "r2"),
        -2.954317312283796e-4,
        epsilon = 1e-9
    );
    assert_abs_diff_eq!(
        sens_value(&result, "r3"),
        0.001129248920024266,
        epsilon = 1e-9
    );
    assert_abs_diff_eq!(
        sens_value(&result, "r4"),
        -2.00582596465584e-4,
        epsilon = 1e-9
    );
    assert_abs_diff_eq!(
        sens_value(&result, "r5"),
        0.003755608696037442,
        epsilon = 1e-9
    );
    assert_abs_diff_eq!(
        sens_value(&result, "rx"),
        -0.001494392173796886,
        epsilon = 1e-9
    );
}
