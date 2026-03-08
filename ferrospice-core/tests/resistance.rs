//! Integration tests ported from ngspice-upstream/tests/resistance/
//!
//! Loads actual .cir files, parses with ferrospice_netlist, simulates with
//! ferrospice_core, and compares output against expected .out reference values.

#[cfg(target_arch = "wasm32")]
wasm_bindgen_test::wasm_bindgen_test_configure!(run_in_browser);
use approx::assert_abs_diff_eq;
use ferrospice_core::simulate_op;
use ferrospice_netlist::{Netlist, SimResult};
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

/// Path to the ngspice resistance test directory.
const RESISTANCE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../ngspice-upstream/tests/resistance"
);

/// Load and parse a .cir file from the resistance test directory.
fn load_cir(name: &str) -> Netlist {
    let path = format!("{RESISTANCE_DIR}/{name}");
    let src = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    Netlist::parse(&src).unwrap_or_else(|e| panic!("cannot parse {path}: {e}"))
}

/// Extract a node voltage from a DC operating point result.
fn op_voltage(result: &SimResult, node: &str) -> f64 {
    let plot = &result.plots[0];
    let name = format!("v({})", node);
    plot.vecs
        .iter()
        .find(|v| v.name == name)
        .unwrap_or_else(|| panic!("no vector {name}"))
        .real[0]
}

/// Extract a branch current from a DC operating point result.
fn op_branch_current(result: &SimResult, vsrc: &str) -> f64 {
    let plot = &result.plots[0];
    let name = format!("{}#branch", vsrc.to_lowercase());
    plot.vecs
        .iter()
        .find(|v| v.name == name)
        .unwrap_or_else(|| panic!("no vector {name}"))
        .real[0]
}

/// Filter output text using the same rules as ngspice check.sh.
///
/// Removes lines containing metadata, analysis keywords, headers, etc.
/// Used for output comparison between ferrospice and ngspice reference output.
fn filter_output(text: &str) -> Vec<String> {
    let filter_patterns = [
        "SPARSE",
        "KLU",
        "CPU",
        "Dynamic",
        "Note",
        "Circuit",
        "Trying",
        "Reference",
        "Date",
        "Doing",
        "---",
        "v-sweep",
        "time",
        "est",
        "Error",
        "Warning",
        "Data",
        "Index",
        "trans",
        "acan",
        "oise",
        "nalysis",
        "ole",
        "Total",
        "memory",
        "urrent",
        "Got",
        "Added",
        "BSIM",
        "bsim",
        "B4SOI",
        "b4soi",
        "codemodel",
        "binary raw file",
        "ngspice",
        "Operating",
    ];

    text.lines()
        .filter(|line| {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                return false;
            }
            !filter_patterns.iter().any(|pat| trimmed.contains(pat))
        })
        .map(|l| {
            // Normalize whitespace for comparison (like diff -w)
            l.split_whitespace().collect::<Vec<_>>().join(" ")
        })
        .filter(|l| !l.is_empty())
        .collect()
}

/// Parse numerical values from ngspice .out reference file data rows.
///
/// Extracts tab-separated data rows (lines starting with a number index)
/// and returns (index, Vec<f64>) for each row.
fn extract_data_rows(text: &str) -> Vec<(usize, Vec<f64>)> {
    text.lines()
        .filter_map(|line| {
            let trimmed = line.trim();
            // Data rows start with an integer index
            let mut parts = trimmed.split('\t').filter(|s| !s.is_empty());
            let idx: usize = parts.next()?.trim().parse().ok()?;
            let values: Vec<f64> = parts
                .filter_map(|s| {
                    // Handle complex format "real,imag" by taking just the real part
                    let s = s.trim().trim_end_matches(',');
                    s.parse().ok()
                })
                .collect();
            if values.is_empty() {
                None
            } else {
                Some((idx, values))
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// res_simple.cir — A simple resistor with a voltage source
// ---------------------------------------------------------------------------
// Circuit: R1 1 0 10k, V1 1 0 DC 1
// Analysis: .TRAN 1ns 6ns, .PRINT TRAN I(V1)
// Expected: V(1) = 1V, I(V1) = -0.0001A (constant for all timesteps)

#[test]
fn res_simple_parses() {
    let netlist = load_cir("res_simple.cir");
    assert_eq!(netlist.title, "A simple resistor with a voltage source");
    assert!(netlist.elements().count() >= 2); // R1 and V1
}

#[test]
fn res_simple_dc_op() {
    // Even though .cir uses .TRAN, we can compute the DC operating point.
    // For a purely resistive DC circuit, this matches the "Initial Transient Solution"
    // shown in the ngspice .out file.
    let netlist = load_cir("res_simple.cir");
    let result = simulate_op(&netlist).unwrap();

    // From .out: Node 1 = 1V, v1#branch = -0.0001
    assert_abs_diff_eq!(op_voltage(&result, "1"), 1.0, epsilon = 1e-9);
    assert_abs_diff_eq!(op_branch_current(&result, "V1"), -1e-4, epsilon = 1e-9);
}

#[test]
fn res_simple_transient_data_constant() {
    // Verify that the expected .out data has constant I(V1) = -0.0001
    // at all timesteps (since it's a purely resistive DC circuit).
    let out_path = format!("{RESISTANCE_DIR}/res_simple.out");
    let out_text = std::fs::read_to_string(&out_path).unwrap();
    let rows = extract_data_rows(&out_text);

    assert!(!rows.is_empty(), "should have data rows in .out");
    for (idx, values) in &rows {
        assert!(
            values.len() >= 2,
            "row {idx} should have time and current columns"
        );
        // Second column is I(V1), should be -0.0001 everywhere
        assert_abs_diff_eq!(values[1], -1e-4, epsilon = 1e-12);
    }
}

#[ignore = "TODO: requires transient analysis (.tran) - US-009"]
#[test]
fn res_simple_full_output_comparison() {
    // Full transient output comparison against .out reference.
    // Requires .tran simulation to produce time-series data.
    let _netlist = load_cir("res_simple.cir");
}

// ---------------------------------------------------------------------------
// res_partition.cir — Resistive partition with different AC/DC ratios
// ---------------------------------------------------------------------------
// Circuit: vin 1 0 DC 1V AC 1V, r1 1 2 5K, r2 2 0 5K ac=15k
// Analysis: .OP, .AC DEC 10 1 10K
// Expected DC: V(1)=1, V(2)=0.5, I(vin)=-0.0001

#[test]
fn res_partition_parses() {
    let netlist = load_cir("res_partition.cir");
    assert_eq!(
        netlist.title,
        "* Resistive partition with different ratios for AC/DC (Print V(2))"
    );
}

#[test]
fn res_partition_dc_op() {
    let netlist = load_cir("res_partition.cir");
    let result = simulate_op(&netlist).unwrap();

    // From .out: V(1)=1.0, V(2)=0.5, vin#branch=-0.0001
    assert_abs_diff_eq!(op_voltage(&result, "1"), 1.0, epsilon = 1e-9);
    assert_abs_diff_eq!(op_voltage(&result, "2"), 0.5, epsilon = 1e-9);
    assert_abs_diff_eq!(op_branch_current(&result, "vin"), -1e-4, epsilon = 1e-9);
}

#[test]
fn res_partition_ac_expected_values() {
    // Verify the .out AC data: V(2) = 0.75 at all frequencies
    // (AC uses r2=15k: V(2) = 15k/(5k+15k) = 0.75)
    let out_path = format!("{RESISTANCE_DIR}/res_partition.out");
    let out_text = std::fs::read_to_string(&out_path).unwrap();
    let rows = extract_data_rows(&out_text);

    // AC data rows should have complex frequency and complex V(2)
    assert!(!rows.is_empty(), "should have AC data rows");
    for (_idx, values) in &rows {
        // Check real part of V(2) is 0.75 (purely resistive, no imaginary)
        if values.len() >= 3 {
            assert_abs_diff_eq!(values[2], 0.75, epsilon = 1e-6);
        }
    }
}

#[ignore = "TODO: requires AC analysis (.ac) - US-012"]
#[test]
fn res_partition_full_output_comparison() {
    // Full AC output comparison against .out reference.
    let _netlist = load_cir("res_partition.cir");
}

// ---------------------------------------------------------------------------
// res_array.cir — Parallel resistor array with models, geometry, multipliers
// ---------------------------------------------------------------------------
// Circuit: 5 parallel resistor branches with various features:
//   R1: 10k basic, R2: 10k ac=5k, R3: model(RSH=1000) L=11u W=2u ac=2.5k,
//   R4: 10k ac=5k m=2, R5: 10 ac=5 scale=1k
// Analysis: .OP, .TRAN 1ns 10ns, .AC DEC 100 1MEG 100MEG
// Expected DC OP from .out:
//   All nodes V(1)-V(6) = 1.0V (all connected to Vin=1V through 0V sources)
//   I(VR1)=1e-4, I(VR2)=1e-4, I(VR3)=9.09e-5, I(VR4)=2e-4, I(VR5)=1e-4

#[test]
fn res_array_parses() {
    let netlist = load_cir("res_array.cir");
    assert_eq!(netlist.title, "A paralled resistor array");
    // Should have Vin + 5 VR sources + 5 resistors + model = 11 elements
    assert!(netlist.elements().count() >= 11);
}

#[ignore = "TODO: requires resistor model support (RSH, geometry W/L, multiplier m, scale) and AC/transient analysis"]
#[test]
fn res_array_dc_op() {
    // Model-based resistor R3 uses rmodel1 with RSH=1000, which requires
    // model evaluation to compute effective resistance. Not yet supported.
    let netlist = load_cir("res_array.cir");
    let result = simulate_op(&netlist).unwrap();

    // Expected from .out:
    assert_abs_diff_eq!(op_voltage(&result, "1"), 1.0, epsilon = 1e-9);
    assert_abs_diff_eq!(op_branch_current(&result, "VR1"), 1e-4, epsilon = 1e-9);
    assert_abs_diff_eq!(op_branch_current(&result, "VR2"), 1e-4, epsilon = 1e-9);
    assert_abs_diff_eq!(
        op_branch_current(&result, "VR3"),
        9.090909e-5,
        epsilon = 1e-9
    );
    assert_abs_diff_eq!(op_branch_current(&result, "VR4"), 2e-4, epsilon = 1e-9);
    assert_abs_diff_eq!(op_branch_current(&result, "VR5"), 1e-4, epsilon = 1e-9);
}

#[ignore = "TODO: requires transient analysis (.tran) - US-009, AC analysis (.ac) - US-012, and resistor model support"]
#[test]
fn res_array_full_output_comparison() {
    let _netlist = load_cir("res_array.cir");
}

// ---------------------------------------------------------------------------
// Output filtering / comparison utilities test
// ---------------------------------------------------------------------------

#[test]
fn test_filter_output() {
    let text = r#"
Circuit: A simple resistor

Doing analysis at TEMP = 300.150000 and TNOM = 300.150000

Initial Transient Solution
--------------------------

Node                                   Voltage
----                                   -------
1                                            1
v1#branch                              -0.0001


No. of Data Rows : 60
"#;
    let filtered = filter_output(text);
    // Should keep: "Node Voltage", "1 1", "v1#branch -0.0001"
    // Should filter: "Circuit:", "Doing analysis", "---", "Data"
    assert!(
        !filtered.iter().any(|l| l.contains("Circuit")),
        "Circuit line should be filtered"
    );
    assert!(
        !filtered.iter().any(|l| l.contains("Doing")),
        "Doing line should be filtered"
    );
    assert!(
        !filtered.iter().any(|l| l.contains("Data")),
        "Data line should be filtered"
    );
    // Numerical content should survive
    assert!(
        filtered.iter().any(|l| l.contains("-0.0001")),
        "v1#branch value should survive filtering"
    );
}

#[test]
fn test_extract_data_rows() {
    let text = r#"Index   time            v1#branch
--------------------------------------------------------------------------------
0	0.000000e+00	-1.000000e-04
1	3.000000e-13	-1.000000e-04
2	6.000000e-13	-1.000000e-04
"#;
    let rows = extract_data_rows(text);
    assert_eq!(rows.len(), 3);
    assert_eq!(rows[0].0, 0);
    assert_abs_diff_eq!(rows[0].1[0], 0.0, epsilon = 1e-15);
    assert_abs_diff_eq!(rows[0].1[1], -1e-4, epsilon = 1e-15);
    assert_eq!(rows[2].0, 2);
}
