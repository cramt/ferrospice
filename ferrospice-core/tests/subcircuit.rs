//! Integration tests for subcircuit expansion (.subckt / X).

use approx::assert_abs_diff_eq;
use ferrospice_netlist::Netlist;

/// Port of ngspice-upstream/tests/regression/subckt-processing/global-1.cir
/// Tests that multiple .global cards are accumulative and global nodes are
/// properly resolved inside subcircuit instances.
#[test]
fn test_global_nodes() {
    let netlist = Netlist::parse(
        "check treatment of multiple .global cards

.global n1001_t n1002_t
.global n1003_t

.subckt su1 1
R2        1 0  10k
R3  n1001_t 0  2000.0
.ends

.subckt su2 1
R2        1 0  10k
R3  n1002_t 0  4000.0
.ends

.subckt su3 1
R2        1 0  10k
R3  n1003_t 0  8000.0
.ends

i1  n1001_t 0  -3mA
R1  n1001_t 0  1000.0

i2  n1002_t 0  -5mA
R2  n1002_t 0  1000.0

i3  n1003_t 0  -9mA
R3  n1003_t 0  1000.0

X1  n1    su1
X2  n1    su2
X3  n1    su3
Rn  n1 0  100k

.op
.end
",
    )
    .unwrap();

    let result = ferrospice_core::simulate_op(&netlist).unwrap();
    let plot = &result.plots[0];

    // Global nodes should have expected voltages:
    // n1001_t: 3mA through (1k || 2k) = 2.0V
    // n1002_t: 5mA through (1k || 4k) = 4.0V
    // n1003_t: 9mA through (1k || 8k) = 8.0V
    let v_n1001 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(n1001_t)")
        .unwrap()
        .real[0];
    let v_n1002 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(n1002_t)")
        .unwrap()
        .real[0];
    let v_n1003 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(n1003_t)")
        .unwrap()
        .real[0];

    assert_abs_diff_eq!(v_n1001, 2.0, epsilon = 1e-9);
    assert_abs_diff_eq!(v_n1002, 4.0, epsilon = 1e-9);
    assert_abs_diff_eq!(v_n1003, 8.0, epsilon = 1e-9);
}

#[test]
fn test_voltage_divider_as_subcircuit() {
    // Voltage divider implemented as a subcircuit, output to ground.
    let netlist = Netlist::parse(
        "Voltage divider subcircuit test
.subckt VDIV in out
R1 in mid 1k
R2 mid out 1k
.ends VDIV
V1 1 0 10
X1 1 0 VDIV
.op
.end
",
    )
    .unwrap();

    let result = ferrospice_core::simulate_op(&netlist).unwrap();
    let plot = &result.plots[0];

    // Find v(1) and v(x1.mid)
    let v1 = plot.vecs.iter().find(|v| v.name == "v(1)").unwrap().real[0];
    let v_mid = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(x1.mid)")
        .unwrap()
        .real[0];

    assert_abs_diff_eq!(v1, 10.0, epsilon = 1e-6);
    assert_abs_diff_eq!(v_mid, 5.0, epsilon = 1e-9);
}

#[test]
fn test_voltage_divider_subcircuit_matches_direct() {
    // Direct circuit
    let direct = Netlist::parse(
        "Direct voltage divider
V1 1 0 5
R1 1 mid 1k
R2 mid 0 1k
.op
.end
",
    )
    .unwrap();
    let direct_result = ferrospice_core::simulate_op(&direct).unwrap();

    // Same circuit as subcircuit
    let subckt = Netlist::parse(
        "Subcircuit voltage divider
.subckt VDIV in out
R1 in mid 1k
R2 mid out 1k
.ends VDIV
V1 1 0 5
X1 1 0 VDIV
.op
.end
",
    )
    .unwrap();
    let subckt_result = ferrospice_core::simulate_op(&subckt).unwrap();

    // Both should give same V1 branch current
    let direct_iv1 = direct_result.plots[0]
        .vecs
        .iter()
        .find(|v| v.name == "v1#branch")
        .unwrap()
        .real[0];
    let subckt_iv1 = subckt_result.plots[0]
        .vecs
        .iter()
        .find(|v| v.name == "v1#branch")
        .unwrap()
        .real[0];

    assert_abs_diff_eq!(direct_iv1, subckt_iv1, epsilon = 1e-12);
}

#[test]
fn test_nested_subcircuit_simulation() {
    // Two-stage resistor chain using nested subcircuits
    let netlist = Netlist::parse(
        "Nested subcircuit simulation
.subckt RUNIT a b
R1 a b 1k
.ends RUNIT
.subckt CHAIN a b
X1 a mid RUNIT
X2 mid b RUNIT
.ends CHAIN
V1 in 0 10
X1 in 0 CHAIN
.op
.end
",
    )
    .unwrap();

    let result = ferrospice_core::simulate_op(&netlist).unwrap();
    let plot = &result.plots[0];

    let v_in = plot.vecs.iter().find(|v| v.name == "v(in)").unwrap().real[0];
    let v_mid = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(x1.mid)")
        .unwrap()
        .real[0];

    assert_abs_diff_eq!(v_in, 10.0, epsilon = 1e-9);
    assert_abs_diff_eq!(v_mid, 5.0, epsilon = 1e-9); // voltage divider midpoint
}

#[test]
fn test_subcircuit_with_parameter_substitution() {
    // Subcircuit with parameterized resistance
    let netlist = Netlist::parse(
        "Parameterized subcircuit
.subckt RLOAD p n PARAMS: rval=1k
R1 p n rval
.ends RLOAD
V1 1 0 10
X1 1 0 RLOAD PARAMS: rval=5k
.op
.end
",
    )
    .unwrap();

    let result = ferrospice_core::simulate_op(&netlist).unwrap();
    let plot = &result.plots[0];

    // I = V/R = 10/5000 = 2mA
    let i_v1 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v1#branch")
        .unwrap()
        .real[0];

    assert_abs_diff_eq!(i_v1.abs(), 2e-3, epsilon = 1e-9);
}

#[test]
fn test_multiple_subcircuit_instances() {
    // Three instances of the same subcircuit
    let netlist = Netlist::parse(
        "Multiple instances
.subckt RBUF in out
R1 in out 1k
.ends RBUF
V1 1 0 4
X1 1 2 RBUF
X2 2 3 RBUF
X3 3 4 RBUF
R_load 4 0 1k
.op
.end
",
    )
    .unwrap();

    let result = ferrospice_core::simulate_op(&netlist).unwrap();
    let plot = &result.plots[0];

    // 4 x 1k resistors in series, V=4V, I=4/4000=1mA
    // V(2) = 4 * 3/4 = 3V, V(3) = 4 * 2/4 = 2V, V(4) = 4 * 1/4 = 1V
    let v2 = plot.vecs.iter().find(|v| v.name == "v(2)").unwrap().real[0];
    let v3 = plot.vecs.iter().find(|v| v.name == "v(3)").unwrap().real[0];
    let v4 = plot.vecs.iter().find(|v| v.name == "v(4)").unwrap().real[0];

    assert_abs_diff_eq!(v2, 3.0, epsilon = 1e-9);
    assert_abs_diff_eq!(v3, 2.0, epsilon = 1e-9);
    assert_abs_diff_eq!(v4, 1.0, epsilon = 1e-9);
}

#[test]
fn test_subcircuit_with_diode() {
    // Subcircuit containing a diode (nonlinear)
    let netlist = Netlist::parse(
        "Subcircuit with diode
.subckt DCLAMP a k
.model DMOD D IS=1e-14
D1 a k DMOD
.ends DCLAMP
V1 in 0 5
R1 in out 1k
X1 out 0 DCLAMP
.op
.end
",
    )
    .unwrap();

    let result = ferrospice_core::simulate_op(&netlist).unwrap();
    let plot = &result.plots[0];

    let v_out = plot.vecs.iter().find(|v| v.name == "v(out)").unwrap().real[0];

    // Diode forward voltage should be ~0.6-0.8V
    assert!(
        v_out > 0.5 && v_out < 0.9,
        "v_out = {} expected ~0.6-0.8V",
        v_out
    );
}

#[test]
fn test_subcircuit_dc_sweep() {
    // DC sweep through subcircuit
    let netlist = Netlist::parse(
        "DC sweep with subcircuit
.subckt RLOAD p n
R1 p n 2k
.ends RLOAD
V1 in 0 DC 0
X1 in 0 RLOAD
.dc V1 0 10 5
.end
",
    )
    .unwrap();

    let result = ferrospice_core::simulate_dc(&netlist).unwrap();
    let plot = &result.plots[0];

    // 3 data points: V=0, V=5, V=10
    let v1_vec = plot.vecs.iter().find(|v| v.name == "v1").unwrap();
    assert_eq!(v1_vec.real.len(), 3);

    let i_vec = plot.vecs.iter().find(|v| v.name == "v1#branch").unwrap();
    // At V=10: I = 10/2000 = 5mA
    assert_abs_diff_eq!(i_vec.real[2].abs(), 5e-3, epsilon = 1e-9);
}
