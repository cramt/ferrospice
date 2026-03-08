use ferrospice_netlist::{Netlist, SimPlot, SimResult, SimVector};

use crate::mna::{MnaError, assemble_mna};

/// Compute the DC operating point of a circuit.
///
/// Assembles the MNA system from the netlist, solves it, and returns
/// a `SimResult` containing node voltages and voltage source branch currents.
pub fn simulate_op(netlist: &Netlist) -> Result<SimResult, MnaError> {
    let mna = assemble_mna(netlist)?;
    let solution = mna.solve()?;

    let mut vecs = Vec::new();

    // Node voltages
    for (name, _idx) in mna.node_map.iter() {
        let v = solution.voltage(name).unwrap_or(0.0);
        vecs.push(SimVector {
            name: format!("v({})", name),
            real: vec![v],
            complex: vec![],
        });
    }

    // Voltage source branch currents
    for vsrc in &mna.vsource_names {
        let i = solution.branch_current(vsrc).unwrap_or(0.0);
        vecs.push(SimVector {
            name: format!("{}#branch", vsrc.to_lowercase()),
            real: vec![i],
            complex: vec![],
        });
    }

    Ok(SimResult {
        plots: vec![SimPlot {
            name: "op1".to_string(),
            vecs,
        }],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ferrospice_netlist::Netlist;

    /// Helper to get a node voltage from the OP result.
    fn op_voltage(result: &SimResult, node: &str) -> f64 {
        let plot = &result.plots[0];
        let name = format!("v({})", node);
        plot.vecs
            .iter()
            .find(|v| v.name == name)
            .unwrap_or_else(|| panic!("no vector {name}"))
            .real[0]
    }

    /// Helper to get a branch current from the OP result.
    fn op_branch_current(result: &SimResult, vsrc: &str) -> f64 {
        let plot = &result.plots[0];
        let name = format!("{}#branch", vsrc.to_lowercase());
        plot.vecs
            .iter()
            .find(|v| v.name == name)
            .unwrap_or_else(|| panic!("no vector {name}"))
            .real[0]
    }

    #[test]
    fn test_voltage_divider_op() {
        // V1=5V, R1=1k, R2=1k → V(mid) = 2.5V
        let netlist = Netlist::parse(
            "Voltage divider
V1 1 0 5
R1 1 mid 1k
R2 mid 0 1k
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        assert_eq!(result.plots.len(), 1);
        assert_eq!(result.plots[0].name, "op1");

        assert_abs_diff_eq!(op_voltage(&result, "1"), 5.0, epsilon = 1e-9);
        assert_abs_diff_eq!(op_voltage(&result, "mid"), 2.5, epsilon = 1e-9);

        // I(V1) = V/R_total = 5/2000 = 2.5mA (negative by convention)
        let i_v1 = op_branch_current(&result, "V1");
        assert_abs_diff_eq!(i_v1.abs(), 2.5e-3, epsilon = 1e-9);
    }

    #[test]
    fn test_current_source_with_series_resistors() {
        // I1=1mA from 0 to 1, R1=1k (1 to mid), R2=2k (mid to 0)
        // Total R = 3k, V(1) = 1e-3 * 3000 = 3.0V, V(mid) = 1e-3 * 2000 = 2.0V
        let netlist = Netlist::parse(
            "Series resistors with current source
I1 0 1 1m
R1 1 mid 1k
R2 mid 0 2k
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        assert_abs_diff_eq!(op_voltage(&result, "1"), 3.0, epsilon = 1e-9);
        assert_abs_diff_eq!(op_voltage(&result, "mid"), 2.0, epsilon = 1e-9);
    }

    #[test]
    fn test_res_simple_op() {
        // Port of ngspice-upstream/tests/resistance/res_simple.cir
        // R1=10k, V1=1V DC → I(V1) = -1/10000 = -0.0001A
        let netlist = Netlist::parse(
            "A simple resistor with a voltage source
R1 1 0 10k
V1 1 0 DC 1
.op
.end
",
        )
        .unwrap();

        let result = simulate_op(&netlist).unwrap();

        assert_abs_diff_eq!(op_voltage(&result, "1"), 1.0, epsilon = 1e-9);

        let i_v1 = op_branch_current(&result, "V1");
        assert_abs_diff_eq!(i_v1, -1e-4, epsilon = 1e-9);
    }
}
