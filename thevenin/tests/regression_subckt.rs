//! Regression tests for subcircuit processing.
//!
//! Ported from ngspice-upstream/tests/regression/subckt-processing/

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

/// global-1: Treatment of multiple .global cards (accumulative behavior).
/// Ported from ngspice-upstream/tests/regression/subckt-processing/global-1.cir
///
/// Three global nodes (n1001_t, n1002_t, n1003_t) declared across two .global cards.
/// Each global node has a current source and parallel resistors (including one inside
/// a subcircuit that connects to the global node).
///
/// n1001_t: I=-3mA, R=1k||2k = 2k/3 → V = 3mA * 2k/3 = 2.0V
/// n1002_t: I=-5mA, R=1k||4k = 4k/5 → V = 5mA * 4k/5 = 4.0V
/// n1003_t: I=-9mA, R=1k||8k = 8k/9 → V = 9mA * 8k/9 = 8.0V
#[test]
fn global_multiple_cards() {
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

    let result = thevenin::simulate_op(&netlist).unwrap();

    // Global nodes should have their expected voltages
    assert_abs_diff_eq!(op_voltage(&result, "n1001_t"), 2.0, epsilon = 1e-6);
    assert_abs_diff_eq!(op_voltage(&result, "n1002_t"), 4.0, epsilon = 1e-6);
    assert_abs_diff_eq!(op_voltage(&result, "n1003_t"), 8.0, epsilon = 1e-6);
}

/// model-scope-5: Scoping of nested .model definitions.
/// Ported from ngspice-upstream/tests/regression/subckt-processing/model-scope-5.cir
///
/// Tests that model names are resolved correctly in nested subcircuit scopes,
/// with proper shadowing behavior.
#[test]
fn model_scope_nested() {
    let netlist = Netlist::parse(
        "check scoping of nested .model definitions

i1  n1001_t 0  dc=-1
i2  n1002_t 0  dc=-1
i3  n1003_t 0  dc=-1
i4  n1004_t 0  dc=-1
i5  n1005_t 0  dc=-1
i6  n1006_t 0  dc=-1
i7  n1007_t 0  dc=-1

x1    n1001_t  sub1

x2    n1002_t n1003_t n1004_t n1005_t n1006_t n1007_t  sub2

.subckt sub1 2
  .model my r r=2k
  r1  2 0  my
.ends

.subckt sub2 3 41a 41b 42a 42b 5
  r2  3 0   my

  x31  41a 41b  sub3
  x32  42a 42b  sub3

  .subckt sub3 4 5
    .model my  r r=8k
    .model any r r=42
    r5  4 0  1k
    r6  5 0  my
  .ends

  .model just r r=43
  r5  5 0  just
.ends

.model my r r=4k

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();

    // x1 uses sub1's local .model my (r=2k) → v = 1A * 2k = 2k
    assert_abs_diff_eq!(op_voltage(&result, "n1001_t"), 2000.0, epsilon = 1.0);
    // x2.r2 uses top-level .model my (r=4k) → v = 1A * 4k = 4k
    assert_abs_diff_eq!(op_voltage(&result, "n1002_t"), 4000.0, epsilon = 1.0);
    // x2.x31.r5 = 1k (literal) → v = 1A * 1k = 1k
    assert_abs_diff_eq!(op_voltage(&result, "n1003_t"), 1000.0, epsilon = 1.0);
    // x2.x31.r6 uses sub3's local .model my (r=8k) → v = 1A * 8k = 8k
    assert_abs_diff_eq!(op_voltage(&result, "n1004_t"), 8000.0, epsilon = 1.0);
    // x2.x32.r5 = 1k → v = 1k
    assert_abs_diff_eq!(op_voltage(&result, "n1005_t"), 1000.0, epsilon = 1.0);
    // x2.x32.r6 uses sub3's local .model my (r=8k) → v = 8k
    assert_abs_diff_eq!(op_voltage(&result, "n1006_t"), 8000.0, epsilon = 1.0);
    // x2.r5 uses sub2's local .model just (r=43) → v = 1A * 43 = 43
    assert_abs_diff_eq!(op_voltage(&result, "n1007_t"), 43.0, epsilon = 0.01);
}
