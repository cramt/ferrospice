//! Integration tests ported from ngspice-upstream/tests/transient/
//!
//! The fourbitadder test requires BJT model support (US-013) and subcircuit
//! expansion (US-016), so it is marked #[ignore] until those features are
//! implemented.

/// Path to the ngspice transient test directory.
const TRANSIENT_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../ngspice-upstream/tests/transient"
);

#[test]
#[ignore] // Requires BJT model (US-013) and subcircuit expansion (US-016)
fn test_fourbitadder() {
    let path = format!("{TRANSIENT_DIR}/fourbitadder.cir");
    let src = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let netlist = ferrospice_netlist::Netlist::parse(&src)
        .unwrap_or_else(|e| panic!("cannot parse {path}: {e}"));

    // Once BJT and subcircuit support are implemented, this should:
    // 1. Run transient analysis
    // 2. Compare v(1) output against fourbitadder.out reference
    let result = ferrospice_core::simulate_tran(&netlist).unwrap();

    // Load reference output.
    let out_path = format!("{TRANSIENT_DIR}/fourbitadder.out");
    let _out_src = std::fs::read_to_string(&out_path)
        .unwrap_or_else(|e| panic!("cannot read {out_path}: {e}"));

    // Verify we got transient results.
    assert!(!result.plots.is_empty(), "should have at least one plot");
    let plot = &result.plots[0];
    assert_eq!(plot.name, "tran1");

    // Find v(1) in the output.
    let v1 = plot
        .vecs
        .iter()
        .find(|v| v.name == "v(1)")
        .expect("should have v(1) vector");
    assert!(!v1.real.is_empty(), "v(1) should have data points");
}
