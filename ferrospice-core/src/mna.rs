use std::collections::BTreeMap;

use ferrospice_netlist::{Element, ElementKind, Expr, Netlist};
use thiserror::Error;

use crate::LinearSystem;

/// Ground node name — the reference node excluded from the MNA matrix.
const GROUND: &str = "0";

#[derive(Error, Debug)]
pub enum MnaError {
    #[error("unsupported element for MNA assembly: {0}")]
    UnsupportedElement(String),
    #[error("non-numeric value in element {element}: parameter expressions not yet supported")]
    NonNumericValue { element: String },
    #[error("voltage source {0} has no DC value")]
    NoVoltageValue(String),
    #[error("failed to solve MNA system: {0}")]
    SolveError(#[from] crate::SparseMatrixError),
}

/// Maps node names to matrix indices. Ground node "0" is excluded.
#[derive(Debug, Clone)]
pub struct NodeMap {
    /// node name -> matrix index
    map: BTreeMap<String, usize>,
}

impl NodeMap {
    fn new() -> Self {
        Self {
            map: BTreeMap::new(),
        }
    }

    /// Get or assign an index for a node. Returns `None` for ground.
    fn index(&mut self, node: &str) -> Option<usize> {
        if node == GROUND {
            return None;
        }
        let next = self.map.len();
        Some(*self.map.entry(node.to_string()).or_insert(next))
    }

    /// Number of non-ground nodes.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Returns true if there are no non-ground nodes.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Look up a node's index (returns None for ground or unknown nodes).
    pub fn get(&self, node: &str) -> Option<usize> {
        self.map.get(node).copied()
    }

    /// Iterate over (node_name, index) pairs.
    pub fn iter(&self) -> impl Iterator<Item = (&str, usize)> {
        self.map.iter().map(|(k, &v)| (k.as_str(), v))
    }
}

/// The assembled MNA system ready for solving.
#[derive(Debug)]
pub struct MnaSystem {
    /// The linear system (matrix + RHS).
    pub system: LinearSystem,
    /// Mapping from node names to matrix indices (first N entries of solution).
    pub node_map: NodeMap,
    /// Names of voltage sources whose branch currents appear in the solution
    /// (entries N..N+M of solution vector).
    pub vsource_names: Vec<String>,
}

impl MnaSystem {
    /// Solve the MNA system, returning an `MnaSolution`.
    pub fn solve(&self) -> Result<MnaSolution<'_>, MnaError> {
        let x = self.system.solve()?;
        let n = self.node_map.len();
        Ok(MnaSolution {
            values: x,
            num_nodes: n,
            node_map: &self.node_map,
            vsource_names: &self.vsource_names,
        })
    }
}

/// Solution of an MNA system with named access to node voltages and branch currents.
#[derive(Debug)]
pub struct MnaSolution<'a> {
    values: Vec<f64>,
    num_nodes: usize,
    node_map: &'a NodeMap,
    vsource_names: &'a [String],
}

impl MnaSolution<'_> {
    /// Get the voltage at a node by name. Ground returns 0.0.
    pub fn voltage(&self, node: &str) -> Option<f64> {
        if node == GROUND {
            return Some(0.0);
        }
        self.node_map.get(node).map(|i| self.values[i])
    }

    /// Get the branch current through a voltage source by name.
    pub fn branch_current(&self, vsource: &str) -> Option<f64> {
        let vsource_lower = vsource.to_lowercase();
        self.vsource_names
            .iter()
            .position(|n| n.to_lowercase() == vsource_lower)
            .map(|i| self.values[self.num_nodes + i])
    }
}

/// Extract a numeric value from an `Expr`, or return an error.
fn expr_value(expr: &Expr, element_name: &str) -> Result<f64, MnaError> {
    match expr {
        Expr::Num(v) => Ok(*v),
        _ => Err(MnaError::NonNumericValue {
            element: element_name.to_string(),
        }),
    }
}

/// Assemble an MNA system from a parsed netlist.
///
/// Currently supports: resistors (R), independent voltage sources (V),
/// independent current sources (I), capacitors (C, open in DC),
/// and inductors (L, short in DC).
pub fn assemble_mna(netlist: &Netlist) -> Result<MnaSystem, MnaError> {
    // First pass: collect all nodes and count voltage sources to determine matrix size.
    let mut node_map = NodeMap::new();
    let mut vsource_count = 0usize;

    for element in netlist.elements() {
        match &element.kind {
            ElementKind::Resistor { pos, neg, .. } => {
                node_map.index(pos);
                node_map.index(neg);
            }
            ElementKind::VoltageSource { pos, neg, .. } => {
                node_map.index(pos);
                node_map.index(neg);
                vsource_count += 1;
            }
            ElementKind::CurrentSource { pos, neg, .. } => {
                node_map.index(pos);
                node_map.index(neg);
            }
            ElementKind::Capacitor { pos, neg, .. } => {
                // Capacitor is open circuit in DC — just register nodes.
                node_map.index(pos);
                node_map.index(neg);
            }
            ElementKind::Inductor { pos, neg, .. } => {
                // Inductor is short circuit in DC — modeled as 0V voltage source.
                node_map.index(pos);
                node_map.index(neg);
                vsource_count += 1;
            }
            _ => {
                // Skip unsupported elements for now
            }
        }
    }

    let n = node_map.len();
    let dim = n + vsource_count;
    let mut system = LinearSystem::new(dim);
    let mut vsource_names = Vec::with_capacity(vsource_count);
    let mut vsource_idx = 0usize;

    // Second pass: stamp each element.
    for element in netlist.elements() {
        stamp_element(
            element,
            &node_map,
            &mut system,
            &mut vsource_names,
            &mut vsource_idx,
            n,
        )?;
    }

    Ok(MnaSystem {
        system,
        node_map,
        vsource_names,
    })
}

/// Stamp a single element into the MNA system.
fn stamp_element(
    element: &Element,
    node_map: &NodeMap,
    system: &mut LinearSystem,
    vsource_names: &mut Vec<String>,
    vsource_idx: &mut usize,
    num_nodes: usize,
) -> Result<(), MnaError> {
    match &element.kind {
        ElementKind::Resistor {
            pos, neg, value, ..
        } => {
            let g = 1.0 / expr_value(value, &element.name)?;
            let ni = node_map.get(pos);
            let nj = node_map.get(neg);
            stamp_conductance(&mut system.matrix, ni, nj, g);
        }
        ElementKind::VoltageSource { pos, neg, source } => {
            let v = source
                .dc
                .as_ref()
                .map(|e| expr_value(e, &element.name))
                .transpose()?
                .unwrap_or(0.0);
            let ni = node_map.get(pos);
            let nj = node_map.get(neg);
            let branch = num_nodes + *vsource_idx;

            // KCL stamps: branch current enters pos, exits neg
            if let Some(i) = ni {
                system.matrix.add(i, branch, 1.0);
                system.matrix.add(branch, i, 1.0);
            }
            if let Some(j) = nj {
                system.matrix.add(j, branch, -1.0);
                system.matrix.add(branch, j, -1.0);
            }
            // Branch equation: V(pos) - V(neg) = v
            system.rhs[branch] = v;

            vsource_names.push(element.name.clone());
            *vsource_idx += 1;
        }
        ElementKind::CurrentSource { pos, neg, source } => {
            let i_val = source
                .dc
                .as_ref()
                .map(|e| expr_value(e, &element.name))
                .transpose()?
                .unwrap_or(0.0);
            let ni = node_map.get(pos);
            let nj = node_map.get(neg);

            // SPICE convention: current flows from n+ to n- through external circuit,
            // so current exits n+ (subtract from RHS) and enters n- (add to RHS).
            if let Some(i) = ni {
                system.rhs[i] -= i_val;
            }
            if let Some(j) = nj {
                system.rhs[j] += i_val;
            }
        }
        ElementKind::Capacitor { .. } => {
            // Capacitor is open circuit in DC — no contribution to DC matrix.
        }
        ElementKind::Inductor { pos, neg, .. } => {
            // Inductor is short circuit in DC — modeled as 0V voltage source.
            let ni = node_map.get(pos);
            let nj = node_map.get(neg);
            let branch = num_nodes + *vsource_idx;

            if let Some(i) = ni {
                system.matrix.add(i, branch, 1.0);
                system.matrix.add(branch, i, 1.0);
            }
            if let Some(j) = nj {
                system.matrix.add(j, branch, -1.0);
                system.matrix.add(branch, j, -1.0);
            }
            // V(pos) - V(neg) = 0 (short circuit)
            system.rhs[branch] = 0.0;

            vsource_names.push(element.name.clone());
            *vsource_idx += 1;
        }
        _ => {
            // Skip unsupported elements silently for now
        }
    }
    Ok(())
}

/// Stamp a conductance G between nodes ni and nj into the matrix.
/// `None` means ground (not in matrix).
fn stamp_conductance(
    matrix: &mut crate::SparseMatrix,
    ni: Option<usize>,
    nj: Option<usize>,
    g: f64,
) {
    if let Some(i) = ni {
        matrix.add(i, i, g);
    }
    if let Some(j) = nj {
        matrix.add(j, j, g);
    }
    if let (Some(i), Some(j)) = (ni, nj) {
        matrix.add(i, j, -g);
        matrix.add(j, i, -g);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ferrospice_netlist::Netlist;

    #[test]
    fn test_voltage_divider() {
        // V1=5V, R1=1k (between node 1 and mid), R2=1k (between mid and 0)
        // Expected: V(mid) = 2.5V, V(1) = 5V, I(V1) = -2.5mA
        let netlist = Netlist::parse(
            "Voltage divider test
V1 1 0 5
R1 1 mid 1k
R2 mid 0 1k
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();

        // 2 non-ground nodes (1, mid) + 1 voltage source = 3x3 system
        assert_eq!(mna.system.dim(), 3);
        assert_eq!(mna.node_map.len(), 2);
        assert_eq!(mna.vsource_names.len(), 1);

        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("1").unwrap(), 5.0, epsilon = 1e-12);
        assert_abs_diff_eq!(solution.voltage("mid").unwrap(), 2.5, epsilon = 1e-12);
        assert_abs_diff_eq!(solution.voltage("0").unwrap(), 0.0, epsilon = 1e-12);

        // Current through V1: I = V/R_total = 5/2000 = 2.5mA
        // Convention: current flows from pos to neg through source (into node 1),
        // so branch current is negative (current flows out of source into circuit).
        let i_v1 = solution.branch_current("V1").unwrap();
        assert_abs_diff_eq!(i_v1.abs(), 2.5e-3, epsilon = 1e-12);
    }

    #[test]
    fn test_current_source_with_resistors() {
        // I1=1mA from 0 to node 1, R1=1k between node 1 and 0
        // V(1) = I * R = 1e-3 * 1000 = 1.0V
        let netlist = Netlist::parse(
            "Current source test
I1 0 1 1m
R1 1 0 1k
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        assert_eq!(mna.system.dim(), 1); // 1 node, 0 voltage sources

        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("1").unwrap(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_series_resistors_with_voltage_source() {
        // V1=10V, R1=2k and R2=3k in series
        // V(mid) = 10 * 3k/(2k+3k) = 6V
        let netlist = Netlist::parse(
            "Series resistors
V1 in 0 10
R1 in mid 2k
R2 mid 0 3k
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        let solution = mna.solve().unwrap();

        assert_abs_diff_eq!(solution.voltage("in").unwrap(), 10.0, epsilon = 1e-12);
        assert_abs_diff_eq!(solution.voltage("mid").unwrap(), 6.0, epsilon = 1e-12);
    }

    #[test]
    fn test_ground_node_excluded() {
        // Single resistor from node 1 to ground with voltage source
        let netlist = Netlist::parse(
            "Ground test
V1 1 0 3.3
R1 1 0 330
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        // 1 node + 1 voltage source = 2x2
        assert_eq!(mna.system.dim(), 2);
        assert_eq!(mna.node_map.get("0"), None); // ground not in map
        assert!(mna.node_map.get("1").is_some());

        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("1").unwrap(), 3.3, epsilon = 1e-12);
        // I = V/R = 3.3/330 = 0.01A
        let i_v1 = solution.branch_current("V1").unwrap();
        assert_abs_diff_eq!(i_v1.abs(), 0.01, epsilon = 1e-12);
    }

    #[test]
    fn test_multiple_current_sources() {
        // Two current sources into the same node with a resistor
        // I1=2mA from 0 to 1, I2=3mA from 0 to 1, R1=1k from 1 to 0
        // V(1) = (2m + 3m) * 1k = 5V
        let netlist = Netlist::parse(
            "Multiple current sources
I1 0 1 2m
I2 0 1 3m
R1 1 0 1k
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("1").unwrap(), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn test_rc_circuit_dc_op() {
        // RC circuit: V1=10V, R1=1k, C1=1u between mid and 0
        // In DC, capacitor is open circuit, so no current flows through R1.
        // V(mid) = V1 = 10V (no voltage drop across R1)
        let netlist = Netlist::parse(
            "RC circuit DC test
V1 1 0 10
R1 1 mid 1k
C1 mid 0 1u
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        // 2 non-ground nodes (1, mid) + 1 voltage source = 3x3
        // Capacitor adds no branch equation
        assert_eq!(mna.system.dim(), 3);

        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("1").unwrap(), 10.0, epsilon = 1e-12);
        assert_abs_diff_eq!(solution.voltage("mid").unwrap(), 10.0, epsilon = 1e-12);
        // No current flows (capacitor is open)
        assert_abs_diff_eq!(solution.branch_current("V1").unwrap(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_rl_circuit_dc_op() {
        // RL circuit: V1=5V, R1=1k, L1=1m between mid and 0
        // In DC, inductor is short circuit, so full current flows.
        // V(mid) = 0V (inductor shorts mid to ground)
        // I = V1/R1 = 5/1000 = 5mA
        let netlist = Netlist::parse(
            "RL circuit DC test
V1 1 0 5
R1 1 mid 1k
L1 mid 0 1m
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        // 2 non-ground nodes (1, mid) + 1 voltage source + 1 inductor branch = 4x4
        assert_eq!(mna.system.dim(), 4);

        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("1").unwrap(), 5.0, epsilon = 1e-12);
        assert_abs_diff_eq!(solution.voltage("mid").unwrap(), 0.0, epsilon = 1e-12);
        // Current through V1 = -5mA (flows out of source into circuit)
        assert_abs_diff_eq!(
            solution.branch_current("V1").unwrap().abs(),
            5e-3,
            epsilon = 1e-12
        );
        // Inductor branch current = 5mA
        assert_abs_diff_eq!(
            solution.branch_current("L1").unwrap().abs(),
            5e-3,
            epsilon = 1e-12
        );
    }
}
