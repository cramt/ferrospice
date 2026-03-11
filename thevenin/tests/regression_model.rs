//! Regression tests for model handling.
//!
//! Ported from ngspice-upstream/tests/regression/model/

use approx::assert_abs_diff_eq;
use thevenin_types::Netlist;

fn op_current(result: &thevenin_types::SimResult, src: &str) -> f64 {
    let plot = &result.plots[0];
    let key = format!("{}#branch", src.to_lowercase());
    plot.vecs
        .iter()
        .find(|v| v.name == key)
        .unwrap_or_else(|| {
            let names: Vec<_> = plot.vecs.iter().map(|v| &v.name).collect();
            panic!("source current {key} not found in results, available: {names:?}")
        })
        .real[0]
}

/// special-names-1: Model names starting with a digit (e.g. 1n4002, 2sk456).
/// Ported from ngspice-upstream/tests/regression/model/special-names-1.cir
///
/// Verifies that models named like 1n4002, 2sk456, 1smb4148 are accepted and
/// that the saturation current ratios match expectations.
#[test]
fn special_names_diode_models() {
    let netlist = Netlist::parse(
        "check special modelnames starting with a digit
i1 1 0 dc -100u

vm1 1 t1 dc 0
d1 t1 0 dplain

vm2 1 t2 dc 0
d2 t2 0 1n4002

vm3 1 t3 dc 0
d3 t3 0 2sk456

vm4 1 t4 dc 0
d4 t4 0 1smb4148

.model dplain   d(is=1.0f)
.model 1n4002   d(is=1.3f)
.model 2sk456   d(is=1.5f)
.model 1smb4148 d(is=1.7f)

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();

    // Current ratios should match Is ratios (at same forward bias)
    let i1 = op_current(&result, "vm1");
    let i2 = op_current(&result, "vm2");
    let i3 = op_current(&result, "vm3");
    let i4 = op_current(&result, "vm4");

    assert_abs_diff_eq!(i2 / i1, 1.3, epsilon = 1e-6);
    assert_abs_diff_eq!(i3 / i1, 1.5, epsilon = 1e-6);
    assert_abs_diff_eq!(i4 / i1, 1.7, epsilon = 1e-6);
}

/// instance-defaults: .model accepts instance default parameters like resistance=2k.
/// Ported from ngspice-upstream/tests/regression/model/instance-defaults.cir
#[test]
fn model_instance_defaults() {
    let netlist = Netlist::parse(
        "check whether .model accepts instance defaults
v1 1 0 dc=1
r1 1 0 myres

.model myres r(resistance=2k)

.op
.end
",
    )
    .unwrap();

    let result = thevenin::simulate_op(&netlist).unwrap();

    // V=1V, R=2k → I = 1/2k = 0.5mA
    let i_v1 = op_current(&result, "v1");
    assert_abs_diff_eq!(i_v1, -0.5e-3, epsilon = 1e-9);
}
