//! Regression tests for .func user-defined function processing.
//!
//! Ported from ngspice-upstream/tests/regression/func/func-1.cir

use approx::assert_abs_diff_eq;
use thevenin_types::Netlist;

fn op_voltage(result: &thevenin_types::SimResult, node: &str) -> f64 {
    let plot = &result.plots[0];
    let key = format!("v({})", node.to_lowercase());
    plot.vecs
        .iter()
        .find(|v| v.name == key)
        .unwrap_or_else(|| panic!("node {key} not found in results"))
        .real[0]
}

/// func-1: Zero-parameter function.
#[test]
fn func_zero_params() {
    let netlist = Netlist::parse(
        "func zero params
.func foo0() 1013.0
V1 1 0 'foo0()'
R1 1 0 1k
.op
.end
",
    )
    .unwrap();
    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "1"), 1013.0, epsilon = 1e-9);
}

/// func-1: Single-parameter function.
#[test]
fn func_single_param() {
    let netlist = Netlist::parse(
        "func single param
.func bar1(p) p
V1 1 0 'bar1(42)'
R1 1 0 1k
.op
.end
",
    )
    .unwrap();
    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "1"), 42.0, epsilon = 1e-9);
}

/// func-1: Function with arithmetic.
#[test]
fn func_arithmetic() {
    let netlist = Netlist::parse(
        "func arithmetic
.func double(x) x*2
.func add(a,b) a+b
V1 1 0 'double(7)'
V2 2 0 'add(3,4)'
R1 1 0 1k
R2 2 0 1k
.op
.end
",
    )
    .unwrap();
    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "1"), 14.0, epsilon = 1e-9);
    assert_abs_diff_eq!(op_voltage(&result, "2"), 7.0, epsilon = 1e-9);
}

/// func-1: Function with ternary operator.
#[test]
fn func_ternary() {
    let netlist = Netlist::parse(
        "func ternary
.func clamp(x,lo,hi) x<lo ? lo : x>hi ? hi : x
V1 1 0 'clamp(5, 0, 10)'
V2 2 0 'clamp(-5, 0, 10)'
V3 3 0 'clamp(15, 0, 10)'
R1 1 0 1k
R2 2 0 1k
R3 3 0 1k
.op
.end
",
    )
    .unwrap();
    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "1"), 5.0, epsilon = 1e-9);
    assert_abs_diff_eq!(op_voltage(&result, "2"), 0.0, epsilon = 1e-9);
    assert_abs_diff_eq!(op_voltage(&result, "3"), 10.0, epsilon = 1e-9);
}

/// func-1: Parameter name shadowing — function param shadows .param.
#[test]
fn func_param_shadowing() {
    let netlist = Netlist::parse(
        "func shadowing
.param xoo = 100
.func fun1(a, xoo) a*xoo
V1 1 0 'fun1(3, 5)'
R1 1 0 1k
.op
.end
",
    )
    .unwrap();
    let result = thevenin::simulate_op(&netlist).unwrap();
    // fun1(3, 5) = 3*5 = 15 (xoo=5 shadows .param xoo=100)
    assert_abs_diff_eq!(op_voltage(&result, "1"), 15.0, epsilon = 1e-9);
}

/// func-1: Function calling built-in math.
#[test]
fn func_calling_builtins() {
    let netlist = Netlist::parse(
        "func builtins
.func hypot(a,b) sqrt(a*a + b*b)
V1 1 0 'hypot(3,4)'
R1 1 0 1k
.op
.end
",
    )
    .unwrap();
    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "1"), 5.0, epsilon = 1e-9);
}
