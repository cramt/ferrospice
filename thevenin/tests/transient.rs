//! Integration tests ported from ngspice-upstream/tests/transient/
//!
//! The fourbitadder test requires BJT model support (US-013) and subcircuit
//! expansion (US-016), so it is marked #[ignore] until those features are
//! implemented.

const FOURBITADDER_CIR: &str =
    include_str!("fixtures/transient/fourbitadder.cir");
const FOURBITADDER_OUT: &str =
    include_str!("fixtures/transient/fourbitadder.out");

#[test]
#[ignore] // Requires BJT model (US-013) and subcircuit expansion (US-016)
fn test_fourbitadder() {
    let netlist = thevenin_types::Netlist::parse(FOURBITADDER_CIR)
        .unwrap_or_else(|e| panic!("cannot parse fourbitadder.cir: {e}"));

    // Once BJT and subcircuit support are implemented, this should:
    // 1. Run transient analysis
    // 2. Compare v(1) output against fourbitadder.out reference
    let result = thevenin::simulate_tran(&netlist).unwrap();

    // Reference output available as FOURBITADDER_OUT.
    let _ = FOURBITADDER_OUT;

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
