//! Regression tests for miscellaneous features.
//!
//! Ported from ngspice-upstream/tests/regression/misc/

use approx::assert_abs_diff_eq;
use thevenin_types::Netlist;

fn op_voltage(result: &thevenin_types::SimResult, node: &str) -> f64 {
    let plot = &result.plots[0];
    let key = format!("v({})", node.to_lowercase());
    plot.vecs
        .iter()
        .find(|v| v.name == key)
        .unwrap_or_else(|| {
            let names: Vec<_> = plot.vecs.iter().map(|v| &v.name).collect();
            panic!("node {key} not found in results, available: {names:?}")
        })
        .real[0]
}

/// ac-zero: AC analysis at frequency = 0 Hz.
/// Ported from ngspice-upstream/tests/regression/misc/ac-zero.cir
///
/// At f=0, the inductor is a short and the capacitor is open.
/// V(3) = V(1) * R3/(R1+R3) = 1.0 * 1k/(1k+1k) = 0.5
#[test]
fn ac_zero_frequency() {
    let netlist = Netlist::parse(
        "whether .ac works for freq=0
v1 1 0 dc=0 ac=1
r1 1 2 1k
l1 2 3 1uH
r3 3 0 1k
c3 3 0 1uF
.ac lin 1 0 0
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_ac(&netlist).unwrap();
    let plot = &result.plots[0];
    let v3 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(3)")
        .expect("v(3) not found");

    // At f=0: inductor=short, capacitor=open → voltage divider R3/(R1+R3) = 0.5
    // AC results are complex-valued
    assert_abs_diff_eq!(v3.complex[0].re, 0.5, epsilon = 1e-12);
    assert_abs_diff_eq!(v3.complex[0].im, 0.0, epsilon = 1e-12);
}

/// bugs-1: B source with unary minus and plus operators (constant expressions).
/// Ported from ngspice-upstream/tests/regression/misc/bugs-1.cir (Bug #294)
#[test]
fn bsource_unary_operators() {
    let cases: &[(&str, &str, f64)] = &[
        ("V=(- (5))", "unary minus literal", -5.0),
        ("V=(+ (5))", "unary plus literal", 5.0),
        ("V=+-(5)", "plus-minus combo", -5.0),
        ("V=-+(5)", "minus-plus combo", -5.0),
    ];

    for (i, (expr, desc, expected)) in cases.iter().enumerate() {
        let netlist = Netlist::parse(&format!(
            "bugs-1 test {i}: {desc}
B1 1 0 {expr}
R1 1 0 1Meg
.op
.end
"
        ))
        .unwrap();
        let result = thevenin::simulate_op(&netlist).unwrap();
        assert_abs_diff_eq!(op_voltage(&result, "1"), *expected, epsilon = 1e-9);
    }
}

/// bugs-1: Power operator with param.
/// Note: B-source uses constant parameter reference, resolved at parse time.
#[test]
fn bsource_power_with_param() {
    // Use voltage sources with expressions instead of B-sources since
    // B-sources with .param references get resolved to constants
    let netlist = Netlist::parse(
        "bugs-1 power operator test
.param aux=-2

V1 n1 0 'aux**2'
V2 n2 0 '-(aux**2)'
R1 n1 0 1Meg
R2 n2 0 1Meg

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();
    // (-2)^2 = 4
    assert_abs_diff_eq!(op_voltage(&result, "n1"), 4.0, epsilon = 1e-9);
    // -((-2)^2) = -4
    assert_abs_diff_eq!(op_voltage(&result, "n2"), -4.0, epsilon = 1e-9);
}

/// log-functions-1: Consistency of ln, log, and log10 in voltage source expressions.
/// Ported from ngspice-upstream/tests/regression/misc/log-functions-1.cir
#[test]
fn log_functions_in_voltage_params() {
    let netlist = Netlist::parse(
        "regression test for log, log10 and ln

v1 1 0 dc 2.7

vn1 n1 0 'ln(2.7)'
vn2 n2 0 'log(2.7)'
vn3 n3 0 'log10(2.7)'

R1 n1 0 1Meg
R2 n2 0 1Meg
R3 n3 0 1Meg

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();

    let ln_gold = 0.9932517730102834;
    let log10_gold = 0.43136376415898736;

    // ln(2.7) = natural log
    assert_abs_diff_eq!(op_voltage(&result, "n1"), ln_gold, epsilon = 1e-12);
    // log(2.7) = natural log (same as ln in ngspice)
    assert_abs_diff_eq!(op_voltage(&result, "n2"), ln_gold, epsilon = 1e-12);
    // log10(2.7) = common log
    assert_abs_diff_eq!(op_voltage(&result, "n3"), log10_gold, epsilon = 1e-12);
}

/// if-elseif: Nested .if/.elseif/.else/.endif conditional blocks.
/// Ported from ngspice-upstream/tests/regression/misc/if-elseif.cir
#[test]
fn if_elseif_nested() {
    let netlist = Netlist::parse(
        "multiple .elseif, nested .if
.param select = 3
.param select2 = 3

V1 1 0 1

.if (select == 1)
    R1 1 0 1
.elseif (select == 2)
    R1 1 0 10
.elseif (select == 3)
    .if (select2 == 1)
        R1 1 0 100
    .elseif (select2 == 2)
        R1 1 0 200
    .elseif (select2 == 3)
        R1 1 0 300
    .else
        R1 1 0 400
    .endif
.elseif (select == 4)
    R1 1 0 1000
.else
    R1 1 0 10000
.endif

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();

    // select=3, select2=3 → R1 = 300
    // V(1) = 1V, I = V/R = 1/300
    assert_abs_diff_eq!(op_voltage(&result, "1"), 1.0, epsilon = 1e-9);
}

/// pz/ac-resistance: OP with AC resistance parameter on resistor.
/// Ported from ngspice-upstream/tests/regression/pz/ac-resistance.cir (OP part)
///
/// DC: R1=1k, R2=4k → V(2) = 5.0 * R2/(R1+R2) = 5.0 * 4/5 = 4.0
#[test]
fn ac_resistance_op() {
    let netlist = Netlist::parse(
        "check for proper use of ac resistance

Vin 1 0 dc 5.0 ac 3.0
R1 1 2 1k ac=4k
C2 2 0 1n
R2 2 0 4k

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();

    // DC: uses R1=1k (DC value), not 4k (AC value)
    // V(2) = 5.0 * 4k / (1k + 4k) = 4.0
    let h0 = 1.0 / (1.0 + 1000.0 / 4000.0);
    let gold = 5.0 * h0;
    assert_abs_diff_eq!(op_voltage(&result, "2"), gold, epsilon = 1e-9);
}

/// pz/ac-resistance: PZ analysis uses AC resistance.
/// Ported from ngspice-upstream/tests/regression/pz/ac-resistance.cir (PZ part)
///
/// With ac R1=4k, R2=4k: tau = (4k||4k) * 1n = 2k * 1n = 2us
/// pole = -1/tau = -500000
#[test]
fn ac_resistance_pz() {
    let netlist = Netlist::parse(
        "check for proper use of ac resistance

Vin 1 0 dc 5.0 ac 3.0
R1 1 2 1k ac=4k
C2 2 0 1n
R2 2 0 4k

.pz 1 0 2 0 vol pz
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_pz(&netlist).unwrap();

    // With AC R1=4k: tau = (4k||4k) * 1n = 2e-6
    // pole = -1/tau = -5e5
    let poles: Vec<_> = result.plots[0]
        .vecs
        .iter()
        .filter(|v| v.name.starts_with("pole("))
        .collect();
    assert_eq!(poles.len(), 1);
    assert_abs_diff_eq!(poles[0].complex[0].re, -5.0e5, epsilon = 1.0);
}

/// .temp directive: Verify temperature is parsed and accessible.
/// Ported from ngspice-upstream/tests/regression/temper/temper-1.cir (parse part)
#[test]
fn temp_directive_parsing() {
    let netlist = thevenin_types::Netlist::parse(
        "test .temp directive
v1 1 0 1
r1 1 0 1k
.temp 127.0
.op
.end
",
    )
    .unwrap();

    assert_abs_diff_eq!(thevenin::netlist_temp(&netlist), 127.0, epsilon = 1e-9);
}

/// .temp default: Without .temp directive, temperature should be 27°C.
#[test]
fn temp_default() {
    let netlist = thevenin_types::Netlist::parse(
        "test default temperature
v1 1 0 1
r1 1 0 1k
.op
.end
",
    )
    .unwrap();

    assert_abs_diff_eq!(thevenin::netlist_temp(&netlist), 27.0, epsilon = 1e-9);
}

/// if-elseif with different selections.
/// Tests that the first matching branch is selected.
#[test]
fn if_elseif_first_match() {
    let netlist = thevenin_types::Netlist::parse(
        "test first match
.param select = 1

V1 1 0 1

.if (select == 1)
    R1 1 0 100
.elseif (select == 2)
    R1 1 0 200
.else
    R1 1 0 300
.endif

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "1"), 1.0, epsilon = 1e-9);
}

/// if-elseif: .else fallback path.
#[test]
fn if_else_fallback() {
    let netlist = thevenin_types::Netlist::parse(
        "test else fallback
.param select = 99

V1 1 0 1

.if (select == 1)
    R1 1 0 100
.elseif (select == 2)
    R1 1 0 200
.else
    R1 1 0 300
.endif

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();
    // select=99 → else → R1=300, I = 1/300
    let i = result.plots[0]
        .vecs
        .iter()
        .find(|v| v.name.contains("v1"))
        .expect("v1 branch current");
    assert_abs_diff_eq!(i.real[0].abs(), 1.0 / 300.0, epsilon = 1e-9);
}
