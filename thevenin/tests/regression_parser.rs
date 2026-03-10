//! Regression tests for expression parsing and evaluation.
//!
//! Ported from ngspice-upstream/tests/regression/parser/

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

/// minus-minus: Test that `2--3` is parsed and evaluated correctly as 5.
/// Ported from ngspice-upstream/tests/regression/parser/minus-minus.cir
#[test]
fn minus_minus_voltage_source() {
    let netlist = Netlist::parse(
        "minus-minus test
V1 1 0 '2--3'
R1 1 0 1k
.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "1"), 5.0, epsilon = 1e-9);
}

/// minus-minus: Test double-minus in B-source expression.
#[test]
fn minus_minus_bsource() {
    let netlist = Netlist::parse(
        "minus-minus bsource test
B1 1 0 V = 2--3
R1 1 0 1k
.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "1"), 5.0, epsilon = 1e-9);
}

/// xpressn-1: Basic arithmetic expressions in voltage sources.
/// Ported from ngspice-upstream/tests/regression/parser/xpressn-1.cir
#[test]
fn xpressn1_arithmetic() {
    let cases: &[(&str, f64)] = &[
        ("1+2", 3.0),
        ("1+2*3", 7.0),
        ("1-2.1", -1.1),
        ("1--1", 2.0),
        ("5+2/-4", 4.5),
        ("2**3", 8.0),
        ("2^3", 8.0),
        ("(1+2)*3", 9.0),
    ];

    for (i, (expr, expected)) in cases.iter().enumerate() {
        let netlist = Netlist::parse(&format!(
            "xpressn1 test {i}
V1 1 0 '{expr}'
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), *expected, epsilon = 1e-9);
    }
}

/// xpressn-1: Boolean and comparison operators.
#[test]
fn xpressn1_boolean() {
    let cases: &[(&str, f64)] = &[
        ("1&&1", 1.0),
        ("1&&0", 0.0),
        ("0&&1", 0.0),
        ("0&&0", 0.0),
        ("1||1", 1.0),
        ("1||0", 1.0),
        ("0||1", 1.0),
        ("0||0", 0.0),
        ("!0", 1.0),
        ("!1", 0.0),
        ("1>0", 1.0),
        ("0>1", 0.0),
        ("1>=1", 1.0),
        ("1<2", 1.0),
        ("2<1", 0.0),
        ("1<=1", 1.0),
        ("1==1", 1.0),
        ("1==0", 0.0),
        ("1!=0", 1.0),
        ("1!=1", 0.0),
    ];

    for (i, (expr, expected)) in cases.iter().enumerate() {
        let netlist = Netlist::parse(&format!(
            "xpressn1 bool test {i}
V1 1 0 '{expr}'
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), *expected, epsilon = 1e-9);
    }
}

/// xpressn-1: Ternary operator.
#[test]
fn xpressn1_ternary() {
    let cases: &[(&str, f64)] = &[
        ("1?2:3", 2.0),
        ("0?2:3", 3.0),
        ("1?2:0?4:5", 2.0),
        ("0?2:0?4:5", 5.0),
        ("0?2:1?4:5", 4.0),
    ];

    for (i, (expr, expected)) in cases.iter().enumerate() {
        let netlist = Netlist::parse(&format!(
            "xpressn1 ternary test {i}
V1 1 0 '{expr}'
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), *expected, epsilon = 1e-9);
    }
}

/// xpressn-1: Math functions.
#[test]
fn xpressn1_math_functions() {
    let cases: &[(&str, f64)] = &[
        ("sin(0)", 0.0),
        ("cos(0)", 1.0),
        ("tan(0)", 0.0),
        ("asin(0)", 0.0),
        ("acos(1)", 0.0),
        ("atan(0)", 0.0),
        ("exp(0)", 1.0),
        ("log(1)", 0.0),
        ("ln(1)", 0.0),
        ("log10(10)", 1.0),
        ("sqrt(4)", 2.0),
        ("sqr(3)", 9.0),
        ("abs(-5)", 5.0),
        ("sgn(3)", 1.0),
        ("sgn(-3)", -1.0),
        ("sgn(0)", 0.0),
        ("int(3.7)", 3.0),
        ("int(-3.7)", -3.0),
        ("floor(3.7)", 3.0),
        ("floor(-3.2)", -4.0),
        ("ceil(3.2)", 4.0),
        ("ceil(-3.7)", -3.0),
        ("min(2,3)", 2.0),
        ("max(2,3)", 3.0),
        ("pow(2,3)", 8.0),
        ("pwr(2,3)", 8.0),
    ];

    for (i, (expr, expected)) in cases.iter().enumerate() {
        let netlist = Netlist::parse(&format!(
            "xpressn1 func test {i}
V1 1 0 '{expr}'
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), *expected, epsilon = 1e-6);
    }
}

/// xpressn-2: Precision testing of trig/hyperbolic functions.
#[test]
fn xpressn2_precision() {
    let cases: &[(&str, f64)] = &[
        ("sin(0.1)", 0.1_f64.sin()),
        ("cos(0.1)", 0.1_f64.cos()),
        ("tan(0.1)", 0.1_f64.tan()),
        ("asin(0.1)", 0.1_f64.asin()),
        ("acos(0.1)", 0.1_f64.acos()),
        ("atan(0.1)", 0.1_f64.atan()),
        ("sinh(0.1)", 0.1_f64.sinh()),
        ("cosh(0.1)", 0.1_f64.cosh()),
        ("tanh(0.1)", 0.1_f64.tanh()),
        ("exp(0.7)", 0.7_f64.exp()),
        ("log(0.7)", 0.7_f64.ln()),
        ("sqrt(0.7)", 0.7_f64.sqrt()),
    ];

    for (i, (expr, expected)) in cases.iter().enumerate() {
        let netlist = Netlist::parse(&format!(
            "xpressn2 test {i}
V1 1 0 '{expr}'
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        let v = op_voltage(&result, "1");
        let rel_err = if expected.abs() > 1e-20 {
            (v - expected).abs() / expected.abs()
        } else {
            (v - expected).abs()
        };
        assert!(
            rel_err < 1e-12,
            "expr '{expr}': got {v}, expected {expected}, rel_err {rel_err}"
        );
    }
}

/// xpressn-3: nint/floor/ceil rounding functions.
#[test]
fn xpressn3_rounding() {
    let test_vals: [f64; 23] = [
        -2.9, -2.5, -2.1, -2.0, -1.9, -1.5, -1.1, -1.0, -0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9, 1.0,
        1.1, 1.5, 1.9, 2.0, 2.1, 2.5, 2.9,
    ];

    for &val in &test_vals {
        // nint
        let expected_nint = val.round();
        let netlist = Netlist::parse(&format!(
            "xpressn3 nint test
V1 1 0 'nint({val})'
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), expected_nint, epsilon = 1e-9);

        // floor
        let expected_floor = val.floor();
        let netlist = Netlist::parse(&format!(
            "xpressn3 floor test
V1 1 0 'floor({val})'
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), expected_floor, epsilon = 1e-9);

        // ceil
        let expected_ceil = val.ceil();
        let netlist = Netlist::parse(&format!(
            "xpressn3 ceil test
V1 1 0 'ceil({val})'
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), expected_ceil, epsilon = 1e-9);
    }
}

/// bxpressn-1: B-source expression evaluation.
/// Ported from ngspice-upstream/tests/regression/parser/bxpressn-1.cir
#[test]
fn bxpressn1_arithmetic() {
    let cases: &[(&str, f64)] = &[
        ("1+2", 3.0),
        ("1+2*3", 7.0),
        ("1-2.1", -1.1),
        ("1--1", 2.0),
        ("5+2/-4", 4.5),
        ("2**3", 8.0),
        ("2^3", 8.0),
        ("(1+2)*3", 9.0),
    ];

    for (i, (expr, expected)) in cases.iter().enumerate() {
        let netlist = Netlist::parse(&format!(
            "bxpressn1 test {i}
B1 1 0 V = {expr}
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), *expected, epsilon = 1e-9);
    }
}

/// bxpressn-1: Step and predicate functions in B-sources.
#[test]
fn bxpressn1_step_predicates() {
    let cases: &[(&str, f64)] = &[
        ("u(1)", 1.0),
        ("u(-1)", 0.0),
        ("u(0)", 0.0),
        ("u2(0.5)", 0.5),
        ("u2(-1)", 0.0),
        ("u2(2)", 1.0),
        ("uramp(2)", 2.0),
        ("uramp(-1)", 0.0),
        ("eq0(0)", 1.0),
        ("eq0(1)", 0.0),
        ("ne0(0)", 0.0),
        ("ne0(1)", 1.0),
        ("gt0(1)", 1.0),
        ("gt0(-1)", 0.0),
        ("lt0(-1)", 1.0),
        ("lt0(1)", 0.0),
        ("ge0(0)", 1.0),
        ("ge0(-1)", 0.0),
        ("le0(0)", 1.0),
        ("le0(1)", 0.0),
    ];

    for (i, (expr, expected)) in cases.iter().enumerate() {
        let netlist = Netlist::parse(&format!(
            "bxpressn1 step test {i}
B1 1 0 V = {expr}
R1 1 0 1k
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), *expected, epsilon = 1e-9);
    }
}
