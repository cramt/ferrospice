//! Regression tests for .lib file processing and subcircuit scoping.
//!
//! Ported from ngspice-upstream/tests/regression/lib-processing/

use approx::assert_abs_diff_eq;
use std::path::Path;
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

fn simulate_with_libs(cir: &str, base_dir: &Path) -> thevenin_types::SimResult {
    let mut netlist = Netlist::parse(cir).unwrap();
    thevenin::libproc::process_libs(&mut netlist, base_dir).unwrap();
    thevenin::simulate_op(&netlist).unwrap()
}

/// ex1a: Basic .lib processing — load subcircuit from library file.
///
/// Circuit: I1=-1mA, X1 uses sub1 (contains X2=sub_in_lib(R=2k) + R2=2k),
/// sub_in_lib comes from ex1.lib RES section (R3=2k).
/// V(9) = I * R_parallel = -1mA * (2k || 2k) = -1mA * 1k = -1V
/// Vcheck = V(9) - V(check0) = 1V, so V(check0) should be V(9) - 1V
/// Test passes when V(check0) ≈ 0 → V(9) = 1V → 1mA * 1k = 1V ✓
#[test]
fn lib_ex1a() {
    let fixtures =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/regression/lib-processing");
    let cir = std::fs::read_to_string(fixtures.join("ex1a.cir")).unwrap();
    let result = simulate_with_libs(&cir, &fixtures);

    // V(check0) should be ≈ 0 (the test's success criterion)
    assert_abs_diff_eq!(op_voltage(&result, "check0"), 0.0, epsilon = 1e-6);
}

/// ex1b: Library scoping — .lib inside subcircuit, with top-level shadowing.
///
/// sub1 uses .lib 'ex1.lib' RES to get sub_in_lib (R=2k).
/// Top-level sub_in_lib has R=4k.
/// X3 (top-level) uses the 4k version.
/// R2 = 4k at top level.
/// Total parallel: sub1(sub_in_lib=2k) || R2=4k || X3(sub_in_lib=4k) = 2k||4k||4k = 1k
/// V(9) = 1mA * 1k = 1V, Vcheck = 1V → V(check0) = 0
#[test]
fn lib_ex1b() {
    let fixtures =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/regression/lib-processing");
    let cir = std::fs::read_to_string(fixtures.join("ex1b.cir")).unwrap();
    let result = simulate_with_libs(&cir, &fixtures);
    assert_abs_diff_eq!(op_voltage(&result, "check0"), 0.0, epsilon = 1e-6);
}

/// ex2a: Nested library includes.
#[test]
fn lib_ex2a() {
    let fixtures =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/regression/lib-processing");
    let cir = std::fs::read_to_string(fixtures.join("ex2a.cir")).unwrap();
    let result = simulate_with_libs(&cir, &fixtures);
    assert_abs_diff_eq!(op_voltage(&result, "check1"), 0.0, epsilon = 1e-6);
    assert_abs_diff_eq!(op_voltage(&result, "check2"), 0.0, epsilon = 1e-6);
}

/// ex3a: Library inclusion from external file.
#[test]
fn lib_ex3a() {
    let fixtures =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/regression/lib-processing");
    let cir = std::fs::read_to_string(fixtures.join("ex3a.cir")).unwrap();
    let result = simulate_with_libs(&cir, &fixtures);
    assert_abs_diff_eq!(op_voltage(&result, "check1"), 0.0, epsilon = 1e-6);
    assert_abs_diff_eq!(op_voltage(&result, "check2"), 0.0, epsilon = 1e-6);
}

/// scope-1: Subcircuit name scoping — nested subcircuits with same name
/// use the locally-defined version.
///
/// sub1: local sub has R=4k, X1+R1 = 4k||4k = 2k → V = 1mA * 2k = 2V
/// sub2: local sub has R=2k, X1+R1 = 2k||2k = 1k → V = 1mA * 1k = 1V
#[test]
fn scope_1() {
    let netlist = Netlist::parse(
        "scope-1, subckt scopes
i1001_t n1001_t 0 -1mA
x1001_t n1001_t 0 sub1

i1002_t n1002_t 0 -1mA
x1002_t n1002_t 0 sub2

.subckt sub1 n1 n2
.subckt sub n1 n2
R1 n1 n2 4k
.ends
X1 n1 n2 sub
R1 n1 n2 4k
.ends

.subckt sub2 n1 n2
.subckt sub n1 n2
R1 n1 n2 2k
.ends
X1 n1 n2 sub
R1 n1 n2 2k
.ends

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();
    // sub1: 4k||4k = 2k → V = 2V
    assert_abs_diff_eq!(op_voltage(&result, "n1001_t"), 2.0, epsilon = 1e-6);
    // sub2: 2k||2k = 1k → V = 1V
    assert_abs_diff_eq!(op_voltage(&result, "n1002_t"), 1.0, epsilon = 1e-6);
}

/// scope-2: Complex scoping with name shadowing.
///
/// sub1: local sub(4k)||R1(4k) = 2k → 2V
/// sub2: local sub(2k)||R1(2k) = 1k → 1V
/// top-level sub: R1(3k) → 3V (negative because current source is -1mA flowing into node)
#[test]
fn scope_2() {
    let netlist = Netlist::parse(
        "scope-2, subckt scopes
i1001_t n1001_t 0 -1mA
x1001_t n1001_t 0 sub1

i1002_t n1002_t 0 -1mA
x1002_t n1002_t 0 sub2

i1003_t n1003_t 0 -1mA
x1003_t n1003_t 0 sub

.subckt sub1 n1 n2
.subckt sub n1 n2
R1 n1 n2 4k
.ends
X1 n1 n2 sub
R1 n1 n2 4k
.ends

.subckt sub n1 n2
R1 n1 n2 3k
.ends

.subckt sub2 n1 n2
.subckt sub n1 n2
R1 n1 n2 2k
.ends
X1 n1 n2 sub
R1 n1 n2 2k
.ends

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "n1001_t"), 2.0, epsilon = 1e-6);
    assert_abs_diff_eq!(op_voltage(&result, "n1002_t"), 1.0, epsilon = 1e-6);
    assert_abs_diff_eq!(op_voltage(&result, "n1003_t"), 3.0, epsilon = 1e-6);
}

/// scope-3: Parameter scoping in nested subcircuits.
///
/// .param foo = 2k
/// x1001_t: sub1 foo='foo*3' = 6k
///   sub1 default foo=5k, but instance passes 6k
///   Inside sub1: X1 calls sub with foo='3*foo'=18k, R1='5*foo'=30k
///   sub default foo=10k, but instance passes 18k → R1=18k
///   Parallel: 18k || 30k = 11.25k → V = 1mA * 11.25k = 11.25V
///
/// x1002_t: sub2 foo='foo*4' = 8k
///   sub2 default foo=121k, but instance passes 8k
///   Inside sub2: X1 calls sub with foo='11*foo'=88k, R1='foo*13'=104k
///   sub default foo=117k, but instance passes 88k → R1=88k
///   Parallel: 88k || 104k = 47.6667k → V = 1mA * 47.6667k ≈ 47.667V
///
/// x1003_t: sub foo=3k → top-level sub: R1='foo*11'=33k → V = 33V
#[test]
fn scope_3() {
    let netlist = Netlist::parse(
        "scope-3, subckt scopes
.param foo = 2k

i1001_t n1001_t 0 -1mA
x1001_t n1001_t 0 sub1 foo='foo*3'

i1002_t n1002_t 0 -1mA
x1002_t n1002_t 0 sub2 foo='foo*4'

i1003_t n1003_t 0 -1mA
x1003_t n1003_t 0 sub foo=3k

.subckt sub1 n1 n2 foo=5k
.subckt sub n1 n2 foo=10k
R1 n1 n2 'foo'
.ends
X1 n1 n2 sub foo='3*foo'
R1 n1 n2 '5*foo'
.ends

.subckt sub n1 n2 foo=17k
R1 n1 n2 'foo*11'
.ends

.subckt sub2 n1 n2 foo=121k
.subckt sub n1 n2 foo=117k
R1 n1 n2 'foo'
.ends
X1 n1 n2 sub foo='11*foo'
R1 n1 n2 'foo*13'
.ends

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();
    assert_abs_diff_eq!(op_voltage(&result, "n1001_t"), 11.25, epsilon = 1e-6);
    assert_abs_diff_eq!(
        op_voltage(&result, "n1002_t"),
        47.666666666666666,
        epsilon = 1e-6
    );
    assert_abs_diff_eq!(op_voltage(&result, "n1003_t"), 33.0, epsilon = 1e-6);
}
