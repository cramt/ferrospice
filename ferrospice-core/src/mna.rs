use std::collections::BTreeMap;

use ferrospice_netlist::{Element, ElementKind, Expr, Netlist};
use thiserror::Error;

use crate::LinearSystem;
use crate::bjt::{BjtInstance, BjtModel};
use crate::diode::DiodeModel;
use crate::jfet::{JfetInstance, JfetModel};
use crate::mosfet::{MosfetInstance, MosfetModel};
use crate::subckt::flatten_netlist;

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
    #[error("subcircuit expansion error: {0}")]
    SubcktError(#[from] crate::subckt::SubcktError),
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

/// A resolved diode instance with matrix indices and model parameters.
#[derive(Debug, Clone)]
pub struct DiodeInstance {
    /// Diode element name.
    pub name: String,
    /// Anode node matrix index (None = ground).
    pub anode_idx: Option<usize>,
    /// Cathode node matrix index (None = ground).
    pub cathode_idx: Option<usize>,
    /// Internal node index when RS > 0 (between RS and junction).
    pub internal_idx: Option<usize>,
    /// Resolved diode model parameters.
    pub model: DiodeModel,
}

/// A resolved capacitor instance with matrix indices.
#[derive(Debug, Clone)]
pub struct CapacitorInstance {
    /// Positive node matrix index (None = ground).
    pub pos_idx: Option<usize>,
    /// Negative node matrix index (None = ground).
    pub neg_idx: Option<usize>,
    /// Capacitance value in Farads.
    pub capacitance: f64,
    /// Initial condition voltage (from IC= parameter), if specified.
    pub ic: Option<f64>,
}

/// A resolved inductor instance with matrix indices.
#[derive(Debug, Clone)]
pub struct InductorInstance {
    /// Positive node matrix index (None = ground).
    pub pos_idx: Option<usize>,
    /// Negative node matrix index (None = ground).
    pub neg_idx: Option<usize>,
    /// Index of this inductor's branch current in the solution vector.
    pub branch_idx: usize,
    /// Inductance value in Henrys.
    pub inductance: f64,
    /// Initial condition current (from IC= parameter), if specified.
    pub ic: Option<f64>,
}

/// A resolved voltage source instance with matrix indices and waveform.
#[derive(Debug, Clone)]
pub struct VoltageSourceInstance {
    /// Index of this source's branch equation in the RHS vector.
    pub branch_idx: usize,
    /// DC value.
    pub dc_value: f64,
    /// Transient waveform, if any.
    pub waveform: Option<ferrospice_netlist::Waveform>,
}

/// A resolved current source instance with matrix indices and waveform.
#[derive(Debug, Clone)]
pub struct CurrentSourceInstance {
    /// Positive node matrix index (None = ground).
    pub pos_idx: Option<usize>,
    /// Negative node matrix index (None = ground).
    pub neg_idx: Option<usize>,
    /// DC value.
    pub dc_value: f64,
    /// Transient waveform, if any.
    pub waveform: Option<ferrospice_netlist::Waveform>,
}

/// A resolved VCVS (E) instance with matrix indices.
#[derive(Debug, Clone)]
pub struct VcvsInstance {
    /// Branch equation index in the solution vector.
    pub branch_idx: usize,
    /// Output positive node index (None = ground).
    pub out_pos_idx: Option<usize>,
    /// Output negative node index (None = ground).
    pub out_neg_idx: Option<usize>,
    /// Control positive node index (None = ground).
    pub ctrl_pos_idx: Option<usize>,
    /// Control negative node index (None = ground).
    pub ctrl_neg_idx: Option<usize>,
    /// Voltage gain.
    pub gain: f64,
}

/// A resolved VCCS (G) instance with matrix indices.
#[derive(Debug, Clone)]
pub struct VccsInstance {
    /// Output positive node index (None = ground).
    pub out_pos_idx: Option<usize>,
    /// Output negative node index (None = ground).
    pub out_neg_idx: Option<usize>,
    /// Control positive node index (None = ground).
    pub ctrl_pos_idx: Option<usize>,
    /// Control negative node index (None = ground).
    pub ctrl_neg_idx: Option<usize>,
    /// Transconductance.
    pub gm: f64,
}

/// A resolved CCCS (F) instance with matrix indices.
#[derive(Debug, Clone)]
pub struct CccsInstance {
    /// Output positive node index (None = ground).
    pub out_pos_idx: Option<usize>,
    /// Output negative node index (None = ground).
    pub out_neg_idx: Option<usize>,
    /// Branch index of the controlling voltage source.
    pub ctrl_branch_idx: usize,
    /// Current gain.
    pub gain: f64,
}

/// A resolved CCVS (H) instance with matrix indices.
#[derive(Debug, Clone)]
pub struct CcvsInstance {
    /// Branch equation index in the solution vector.
    pub branch_idx: usize,
    /// Output positive node index (None = ground).
    pub out_pos_idx: Option<usize>,
    /// Output negative node index (None = ground).
    pub out_neg_idx: Option<usize>,
    /// Branch index of the controlling voltage source.
    pub ctrl_branch_idx: usize,
    /// Transresistance.
    pub rm: f64,
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
    /// Resolved diode instances for NR iteration.
    pub diodes: Vec<DiodeInstance>,
    /// Resolved capacitor instances for transient analysis.
    pub capacitors: Vec<CapacitorInstance>,
    /// Resolved inductor instances for transient analysis.
    pub inductors: Vec<InductorInstance>,
    /// Resolved BJT instances for NR iteration.
    pub bjts: Vec<BjtInstance>,
    /// Resolved MOSFET instances for NR iteration.
    pub mosfets: Vec<MosfetInstance>,
    /// Resolved JFET instances for NR iteration.
    pub jfets: Vec<JfetInstance>,
    /// Resolved voltage source instances (for transient waveform evaluation).
    pub voltage_sources: Vec<VoltageSourceInstance>,
    /// Resolved current source instances (for transient waveform evaluation).
    pub current_sources: Vec<CurrentSourceInstance>,
    /// Resolved VCVS (E) instances.
    pub vcvs_sources: Vec<VcvsInstance>,
    /// Resolved VCCS (G) instances.
    pub vccs_sources: Vec<VccsInstance>,
    /// Resolved CCCS (F) instances.
    pub cccs_sources: Vec<CccsInstance>,
    /// Resolved CCVS (H) instances.
    pub ccvs_sources: Vec<CcvsInstance>,
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
/// inductors (L, short in DC), and diodes (D, nonlinear).
pub fn assemble_mna(netlist: &Netlist) -> Result<MnaSystem, MnaError> {
    // Flatten subcircuit calls before assembly.
    let flat_netlist = flatten_netlist(netlist)?;
    assemble_mna_flat(&flat_netlist)
}

/// Assemble an MNA system from a flattened (no subcircuit calls) netlist.
fn assemble_mna_flat(netlist: &Netlist) -> Result<MnaSystem, MnaError> {
    // Build a map of model definitions for lookup.
    let models: BTreeMap<String, &ferrospice_netlist::ModelDef> = netlist
        .items
        .iter()
        .filter_map(|item| {
            if let ferrospice_netlist::Item::Model(m) = item {
                Some((m.name.to_uppercase(), m))
            } else {
                None
            }
        })
        .collect();

    // First pass: collect all nodes and count voltage sources to determine matrix size.
    // Also build a vsource name → offset map for F/H controlled source lookup.
    let mut node_map = NodeMap::new();
    let mut vsource_count = 0usize;
    let mut internal_node_count = 0usize;
    let mut vsource_offset_map: BTreeMap<String, usize> = BTreeMap::new();

    for element in netlist.elements() {
        match &element.kind {
            ElementKind::Resistor { pos, neg, .. } => {
                node_map.index(pos);
                node_map.index(neg);
            }
            ElementKind::VoltageSource { pos, neg, .. } => {
                node_map.index(pos);
                node_map.index(neg);
                vsource_offset_map.insert(element.name.to_lowercase(), vsource_count);
                vsource_count += 1;
            }
            ElementKind::CurrentSource { pos, neg, .. } => {
                node_map.index(pos);
                node_map.index(neg);
            }
            ElementKind::Capacitor { pos, neg, .. } => {
                node_map.index(pos);
                node_map.index(neg);
            }
            ElementKind::Inductor { pos, neg, .. } => {
                node_map.index(pos);
                node_map.index(neg);
                vsource_offset_map.insert(element.name.to_lowercase(), vsource_count);
                vsource_count += 1;
            }
            ElementKind::Diode {
                anode,
                cathode,
                model,
                params,
            } => {
                node_map.index(anode);
                node_map.index(cathode);
                // Check if RS > 0 — need an internal node.
                let mut dm = if let Some(mdef) = models.get(&model.to_uppercase()) {
                    DiodeModel::from_model_def(mdef)
                } else {
                    DiodeModel::default()
                };
                dm = dm.with_instance_params(params);
                if dm.has_series_resistance() {
                    internal_node_count += 1;
                }
            }
            ElementKind::Bjt {
                c,
                b,
                e,
                substrate,
                model,
                params,
            } => {
                node_map.index(c);
                node_map.index(b);
                node_map.index(e);
                if let Some(s) = substrate {
                    node_map.index(s);
                }
                let _ = params; // instance params handled in second pass
                let bm = if let Some(mdef) = models.get(&model.to_uppercase()) {
                    BjtModel::from_model_def(mdef)
                } else {
                    BjtModel::new(crate::bjt::BjtType::Npn)
                };
                internal_node_count += bm.internal_node_count();
            }
            ElementKind::Mosfet {
                d,
                g,
                s,
                bulk,
                model,
                params,
            } => {
                node_map.index(d);
                node_map.index(g);
                node_map.index(s);
                node_map.index(bulk);
                let _ = params;
                let mm = if let Some(mdef) = models.get(&model.to_uppercase()) {
                    MosfetModel::from_model_def(mdef)
                } else {
                    MosfetModel::new(crate::mosfet::MosfetType::Nmos)
                };
                internal_node_count += mm.internal_node_count();
            }
            ElementKind::Jfet {
                d,
                g,
                s,
                model,
                params,
            } => {
                node_map.index(d);
                node_map.index(g);
                node_map.index(s);
                let _ = params;
                let jm = if let Some(mdef) = models.get(&model.to_uppercase()) {
                    JfetModel::from_model_def(mdef)
                } else {
                    JfetModel::new(crate::jfet::JfetType::Njf)
                };
                internal_node_count += jm.internal_node_count();
            }
            ElementKind::Vcvs {
                out_pos,
                out_neg,
                in_pos,
                in_neg,
                ..
            } => {
                node_map.index(out_pos);
                node_map.index(out_neg);
                node_map.index(in_pos);
                node_map.index(in_neg);
                vsource_offset_map.insert(element.name.to_lowercase(), vsource_count);
                vsource_count += 1; // VCVS adds a branch equation
            }
            ElementKind::Vccs {
                out_pos,
                out_neg,
                in_pos,
                in_neg,
                ..
            } => {
                node_map.index(out_pos);
                node_map.index(out_neg);
                node_map.index(in_pos);
                node_map.index(in_neg);
                // No branch added — VCCS stamps directly into admittance matrix
            }
            ElementKind::Cccs {
                out_pos, out_neg, ..
            } => {
                node_map.index(out_pos);
                node_map.index(out_neg);
                // No branch added — references existing vsource branch
            }
            ElementKind::Ccvs {
                out_pos, out_neg, ..
            } => {
                node_map.index(out_pos);
                node_map.index(out_neg);
                vsource_offset_map.insert(element.name.to_lowercase(), vsource_count);
                vsource_count += 1; // CCVS adds a branch equation
            }
            _ => {}
        }
    }

    let n = node_map.len() + internal_node_count;
    let dim = n + vsource_count;
    let mut system = LinearSystem::new(dim);
    let mut vsource_names = Vec::with_capacity(vsource_count);
    let mut vsource_idx = 0usize;
    let mut diodes = Vec::new();
    let mut bjts = Vec::new();
    let mut mosfets = Vec::new();
    let mut jfets = Vec::new();
    let mut capacitors = Vec::new();
    let mut inductors = Vec::new();
    let mut voltage_sources = Vec::new();
    let mut current_sources = Vec::new();
    let mut vcvs_sources = Vec::new();
    let mut vccs_sources = Vec::new();
    let mut cccs_sources = Vec::new();
    let mut ccvs_sources = Vec::new();
    let mut internal_idx = node_map.len(); // internal nodes start after external nodes

    // Second pass: stamp each element.
    for element in netlist.elements() {
        match &element.kind {
            ElementKind::Diode {
                anode,
                cathode,
                model,
                params,
            } => {
                let mut dm = if let Some(mdef) = models.get(&model.to_uppercase()) {
                    DiodeModel::from_model_def(mdef)
                } else {
                    DiodeModel::default()
                };
                dm = dm.with_instance_params(params);

                let anode_idx = node_map.get(anode);
                let cathode_idx = node_map.get(cathode);

                let int_idx = if dm.has_series_resistance() {
                    let idx = internal_idx;
                    internal_idx += 1;
                    Some(idx)
                } else {
                    None
                };

                diodes.push(DiodeInstance {
                    name: element.name.clone(),
                    anode_idx,
                    cathode_idx,
                    internal_idx: int_idx,
                    model: dm,
                });
                // Diode stamps are applied during NR iteration, not here.
            }
            ElementKind::Bjt {
                c,
                b,
                e,
                model,
                params,
                ..
            } => {
                let bm = if let Some(mdef) = models.get(&model.to_uppercase()) {
                    BjtModel::from_model_def(mdef)
                } else {
                    BjtModel::new(crate::bjt::BjtType::Npn)
                };
                let bm = bm.with_instance_params(params);

                let col_idx = node_map.get(c);
                let base_idx = node_map.get(b);
                let emit_idx = node_map.get(e);

                // Create internal nodes for series resistances
                let base_prime_idx = if bm.rb > 0.0 {
                    let idx = internal_idx;
                    internal_idx += 1;
                    Some(idx)
                } else {
                    base_idx
                };
                let col_prime_idx = if bm.rc > 0.0 {
                    let idx = internal_idx;
                    internal_idx += 1;
                    Some(idx)
                } else {
                    col_idx
                };
                let emit_prime_idx = if bm.re > 0.0 {
                    let idx = internal_idx;
                    internal_idx += 1;
                    Some(idx)
                } else {
                    emit_idx
                };

                // Extract area and M from instance params
                let mut area = 1.0;
                let mut m_mult = 1.0;
                for p in params {
                    if let Expr::Num(v) = &p.value {
                        match p.name.to_uppercase().as_str() {
                            "AREA" => area = *v,
                            "M" => m_mult = *v,
                            _ => {}
                        }
                    }
                }

                bjts.push(BjtInstance {
                    name: element.name.clone(),
                    col_idx,
                    base_idx,
                    emit_idx,
                    base_prime_idx,
                    col_prime_idx,
                    emit_prime_idx,
                    model: bm,
                    area,
                    m: m_mult,
                });
                // BJT stamps are applied during NR iteration, not here.
            }
            ElementKind::Mosfet {
                d,
                g,
                s,
                bulk,
                model,
                params,
            } => {
                let mm = if let Some(mdef) = models.get(&model.to_uppercase()) {
                    MosfetModel::from_model_def(mdef)
                } else {
                    MosfetModel::new(crate::mosfet::MosfetType::Nmos)
                };

                let drain_idx = node_map.get(d);
                let gate_idx = node_map.get(g);
                let source_idx = node_map.get(s);
                let bulk_idx = node_map.get(bulk);

                // Create internal nodes for series resistances
                let drain_prime_idx = if mm.rd > 0.0 {
                    let idx = internal_idx;
                    internal_idx += 1;
                    Some(idx)
                } else {
                    drain_idx
                };
                let source_prime_idx = if mm.rs > 0.0 {
                    let idx = internal_idx;
                    internal_idx += 1;
                    Some(idx)
                } else {
                    source_idx
                };

                // Extract W, L, AD, AS, PD, PS, M from instance params
                let mut w = 1e-4; // default 100um
                let mut l = 1e-4;
                let mut ad = 0.0;
                let mut as_ = 0.0;
                let mut pd = 0.0;
                let mut ps = 0.0;
                let mut m_mult = 1.0;
                for p in params {
                    if let Expr::Num(v) = &p.value {
                        match p.name.to_uppercase().as_str() {
                            "W" => w = *v,
                            "L" => l = *v,
                            "AD" => ad = *v,
                            "AS" => as_ = *v,
                            "PD" => pd = *v,
                            "PS" => ps = *v,
                            "M" => m_mult = *v,
                            _ => {}
                        }
                    }
                }

                mosfets.push(MosfetInstance {
                    name: element.name.clone(),
                    drain_idx,
                    gate_idx,
                    source_idx,
                    bulk_idx,
                    drain_prime_idx,
                    source_prime_idx,
                    model: mm,
                    w,
                    l,
                    ad,
                    as_,
                    pd,
                    ps,
                    m: m_mult,
                });
                // MOSFET stamps are applied during NR iteration, not here.
            }
            ElementKind::Jfet {
                d,
                g,
                s,
                model,
                params,
            } => {
                let jm = if let Some(mdef) = models.get(&model.to_uppercase()) {
                    JfetModel::from_model_def(mdef)
                } else {
                    JfetModel::new(crate::jfet::JfetType::Njf)
                };

                let drain_idx = node_map.get(d);
                let gate_idx = node_map.get(g);
                let source_idx = node_map.get(s);

                // Create internal nodes for series resistances
                let drain_prime_idx = if jm.rd > 0.0 {
                    let idx = internal_idx;
                    internal_idx += 1;
                    Some(idx)
                } else {
                    drain_idx
                };
                let source_prime_idx = if jm.rs > 0.0 {
                    let idx = internal_idx;
                    internal_idx += 1;
                    Some(idx)
                } else {
                    source_idx
                };

                // Extract AREA and M from instance params
                let mut area = 1.0;
                let mut m_mult = 1.0;
                for p in params {
                    if let Expr::Num(v) = &p.value {
                        match p.name.to_uppercase().as_str() {
                            "AREA" => area = *v,
                            "M" => m_mult = *v,
                            _ => {}
                        }
                    }
                }

                jfets.push(JfetInstance {
                    name: element.name.clone(),
                    drain_idx,
                    gate_idx,
                    source_idx,
                    drain_prime_idx,
                    source_prime_idx,
                    model: jm,
                    area,
                    m: m_mult,
                });
                // JFET stamps are applied during NR iteration, not here.
            }
            ElementKind::Capacitor {
                pos,
                neg,
                value,
                params,
            } => {
                let cap_val = expr_value(value, &element.name)?;
                let pos_idx = node_map.get(pos);
                let neg_idx = node_map.get(neg);
                let ic = extract_ic_param(params);
                capacitors.push(CapacitorInstance {
                    pos_idx,
                    neg_idx,
                    capacitance: cap_val,
                    ic,
                });
                // No DC stamp for capacitor (open circuit in DC).
            }
            ElementKind::Inductor {
                pos,
                neg,
                value,
                params,
            } => {
                let ind_val = expr_value(value, &element.name)?;
                let pos_idx = node_map.get(pos);
                let neg_idx = node_map.get(neg);
                let branch = n + vsource_idx;
                let ic = extract_ic_param(params);
                inductors.push(InductorInstance {
                    pos_idx,
                    neg_idx,
                    branch_idx: branch,
                    inductance: ind_val,
                    ic,
                });
                // Stamp inductor as 0V voltage source (short circuit in DC).
                stamp_inductor_topology(&mut system, pos_idx, neg_idx, branch);
                vsource_names.push(element.name.clone());
                vsource_idx += 1;
            }
            ElementKind::Vcvs {
                out_pos,
                out_neg,
                in_pos,
                in_neg,
                gain,
            } => {
                let gain_val = expr_value(gain, &element.name)?;
                let out_pos_idx = node_map.get(out_pos);
                let out_neg_idx = node_map.get(out_neg);
                let ctrl_pos_idx = node_map.get(in_pos);
                let ctrl_neg_idx = node_map.get(in_neg);
                let branch = n + vsource_idx;

                stamp_vcvs(
                    &mut system,
                    out_pos_idx,
                    out_neg_idx,
                    ctrl_pos_idx,
                    ctrl_neg_idx,
                    branch,
                    gain_val,
                );

                vcvs_sources.push(VcvsInstance {
                    branch_idx: branch,
                    out_pos_idx,
                    out_neg_idx,
                    ctrl_pos_idx,
                    ctrl_neg_idx,
                    gain: gain_val,
                });

                vsource_names.push(element.name.clone());
                vsource_idx += 1;
            }
            ElementKind::Vccs {
                out_pos,
                out_neg,
                in_pos,
                in_neg,
                gm,
            } => {
                let gm_val = expr_value(gm, &element.name)?;
                let out_pos_idx = node_map.get(out_pos);
                let out_neg_idx = node_map.get(out_neg);
                let ctrl_pos_idx = node_map.get(in_pos);
                let ctrl_neg_idx = node_map.get(in_neg);

                stamp_vccs(
                    &mut system,
                    out_pos_idx,
                    out_neg_idx,
                    ctrl_pos_idx,
                    ctrl_neg_idx,
                    gm_val,
                );

                vccs_sources.push(VccsInstance {
                    out_pos_idx,
                    out_neg_idx,
                    ctrl_pos_idx,
                    ctrl_neg_idx,
                    gm: gm_val,
                });
            }
            ElementKind::Cccs {
                out_pos,
                out_neg,
                vsrc,
                gain,
            } => {
                let gain_val = expr_value(gain, &element.name)?;
                let out_pos_idx = node_map.get(out_pos);
                let out_neg_idx = node_map.get(out_neg);
                let ctrl_offset =
                    vsource_offset_map
                        .get(&vsrc.to_lowercase())
                        .ok_or_else(|| {
                            MnaError::UnsupportedElement(format!(
                                "controlling voltage source '{}' not found for {}",
                                vsrc, element.name
                            ))
                        })?;
                let ctrl_branch_idx = n + ctrl_offset;

                stamp_cccs(
                    &mut system,
                    out_pos_idx,
                    out_neg_idx,
                    ctrl_branch_idx,
                    gain_val,
                );

                cccs_sources.push(CccsInstance {
                    out_pos_idx,
                    out_neg_idx,
                    ctrl_branch_idx,
                    gain: gain_val,
                });
            }
            ElementKind::Ccvs {
                out_pos,
                out_neg,
                vsrc,
                rm,
            } => {
                let rm_val = expr_value(rm, &element.name)?;
                let out_pos_idx = node_map.get(out_pos);
                let out_neg_idx = node_map.get(out_neg);
                let ctrl_offset =
                    vsource_offset_map
                        .get(&vsrc.to_lowercase())
                        .ok_or_else(|| {
                            MnaError::UnsupportedElement(format!(
                                "controlling voltage source '{}' not found for {}",
                                vsrc, element.name
                            ))
                        })?;
                let ctrl_branch_idx = n + ctrl_offset;
                let branch = n + vsource_idx;

                stamp_ccvs(
                    &mut system,
                    out_pos_idx,
                    out_neg_idx,
                    ctrl_branch_idx,
                    branch,
                    rm_val,
                );

                ccvs_sources.push(CcvsInstance {
                    branch_idx: branch,
                    out_pos_idx,
                    out_neg_idx,
                    ctrl_branch_idx,
                    rm: rm_val,
                });

                vsource_names.push(element.name.clone());
                vsource_idx += 1;
            }
            _ => {
                stamp_element(
                    element,
                    &node_map,
                    &mut system,
                    &mut vsource_names,
                    &mut vsource_idx,
                    n,
                    &mut voltage_sources,
                    &mut current_sources,
                )?;
            }
        }
    }

    Ok(MnaSystem {
        system,
        node_map,
        vsource_names,
        diodes,
        bjts,
        mosfets,
        jfets,
        capacitors,
        inductors,
        voltage_sources,
        current_sources,
        vcvs_sources,
        vccs_sources,
        cccs_sources,
        ccvs_sources,
    })
}

/// Stamp a single element into the MNA system.
#[expect(clippy::too_many_arguments)]
fn stamp_element(
    element: &Element,
    node_map: &NodeMap,
    system: &mut LinearSystem,
    vsource_names: &mut Vec<String>,
    vsource_idx: &mut usize,
    num_nodes: usize,
    voltage_sources: &mut Vec<VoltageSourceInstance>,
    current_sources: &mut Vec<CurrentSourceInstance>,
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

            voltage_sources.push(VoltageSourceInstance {
                branch_idx: branch,
                dc_value: v,
                waveform: source.waveform.clone(),
            });

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

            current_sources.push(CurrentSourceInstance {
                pos_idx: ni,
                neg_idx: nj,
                dc_value: i_val,
                waveform: source.waveform.clone(),
            });
        }
        ElementKind::Capacitor { .. } | ElementKind::Inductor { .. } => {
            // Handled directly in assemble_mna (need instance info for transient).
        }
        _ => {
            // Skip unsupported elements silently for now
        }
    }
    Ok(())
}

/// Stamp the topology of an inductor branch equation (same pattern as voltage source).
fn stamp_inductor_topology(
    system: &mut LinearSystem,
    ni: Option<usize>,
    nj: Option<usize>,
    branch: usize,
) {
    if let Some(i) = ni {
        system.matrix.add(i, branch, 1.0);
        system.matrix.add(branch, i, 1.0);
    }
    if let Some(j) = nj {
        system.matrix.add(j, branch, -1.0);
        system.matrix.add(branch, j, -1.0);
    }
    // V(pos) - V(neg) = 0 (short circuit in DC)
    system.rhs[branch] = 0.0;
}

/// Extract IC= parameter value from element params.
fn extract_ic_param(params: &[ferrospice_netlist::Param]) -> Option<f64> {
    for p in params {
        if p.name.to_uppercase() == "IC"
            && let Expr::Num(v) = &p.value
        {
            return Some(*v);
        }
    }
    None
}

/// Stamp a VCVS (E) element into the MNA system.
///
/// Branch equation: V(out+) - V(out-) = gain * (V(ctrl+) - V(ctrl-))
/// KCL: branch current enters out+ and exits out-.
fn stamp_vcvs(
    system: &mut LinearSystem,
    out_pos: Option<usize>,
    out_neg: Option<usize>,
    ctrl_pos: Option<usize>,
    ctrl_neg: Option<usize>,
    branch: usize,
    gain: f64,
) {
    // KCL stamps (same as voltage source topology)
    if let Some(i) = out_pos {
        system.matrix.add(i, branch, 1.0);
        system.matrix.add(branch, i, 1.0);
    }
    if let Some(j) = out_neg {
        system.matrix.add(j, branch, -1.0);
        system.matrix.add(branch, j, -1.0);
    }
    // Control voltage contribution to branch equation:
    // V(out+) - V(out-) - gain * (V(ctrl+) - V(ctrl-)) = 0
    if let Some(cp) = ctrl_pos {
        system.matrix.add(branch, cp, -gain);
    }
    if let Some(cn) = ctrl_neg {
        system.matrix.add(branch, cn, gain);
    }
    // RHS = 0 (no independent source component)
}

/// Stamp a VCCS (G) element into the MNA system.
///
/// I = gm * (V(ctrl+) - V(ctrl-)), current flows from out+ to out-.
fn stamp_vccs(
    system: &mut LinearSystem,
    out_pos: Option<usize>,
    out_neg: Option<usize>,
    ctrl_pos: Option<usize>,
    ctrl_neg: Option<usize>,
    gm: f64,
) {
    if let Some(i) = out_pos {
        if let Some(cp) = ctrl_pos {
            system.matrix.add(i, cp, gm);
        }
        if let Some(cn) = ctrl_neg {
            system.matrix.add(i, cn, -gm);
        }
    }
    if let Some(j) = out_neg {
        if let Some(cp) = ctrl_pos {
            system.matrix.add(j, cp, -gm);
        }
        if let Some(cn) = ctrl_neg {
            system.matrix.add(j, cn, gm);
        }
    }
}

/// Stamp a CCCS (F) element into the MNA system.
///
/// I_out = gain * I_ctrl, where I_ctrl is the branch current of the controlling source.
fn stamp_cccs(
    system: &mut LinearSystem,
    out_pos: Option<usize>,
    out_neg: Option<usize>,
    ctrl_branch: usize,
    gain: f64,
) {
    if let Some(i) = out_pos {
        system.matrix.add(i, ctrl_branch, gain);
    }
    if let Some(j) = out_neg {
        system.matrix.add(j, ctrl_branch, -gain);
    }
}

/// Stamp a CCVS (H) element into the MNA system.
///
/// Branch equation: V(out+) - V(out-) = rm * I_ctrl
/// KCL: branch current enters out+ and exits out-.
fn stamp_ccvs(
    system: &mut LinearSystem,
    out_pos: Option<usize>,
    out_neg: Option<usize>,
    ctrl_branch: usize,
    branch: usize,
    rm: f64,
) {
    // KCL stamps (same as voltage source topology)
    if let Some(i) = out_pos {
        system.matrix.add(i, branch, 1.0);
        system.matrix.add(branch, i, 1.0);
    }
    if let Some(j) = out_neg {
        system.matrix.add(j, branch, -1.0);
        system.matrix.add(branch, j, -1.0);
    }
    // Control current contribution to branch equation:
    // V(out+) - V(out-) - rm * I_ctrl = 0
    system.matrix.add(branch, ctrl_branch, -rm);
    // RHS = 0
}

/// Stamp a conductance G between nodes ni and nj into the matrix.
/// `None` means ground (not in matrix).
pub fn stamp_conductance(
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

    #[test]
    fn test_vcvs_voltage_follower() {
        // VCVS with gain=1 (voltage follower): V(out) = V(in)
        // V1=5V at node "in", E1 copies to "out", R1=1k load on "out"
        let netlist = Netlist::parse(
            "VCVS voltage follower
V1 in 0 5
R1 in 0 10k
E1 out 0 in 0 1
R2 out 0 1k
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        // 2 nodes (in, out) + 2 branches (V1, E1)
        assert_eq!(mna.vsource_names.len(), 2);

        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("in").unwrap(), 5.0, epsilon = 1e-12);
        assert_abs_diff_eq!(solution.voltage("out").unwrap(), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn test_vcvs_amplifier() {
        // VCVS with gain=10: V(out) = 10 * V(in)
        // V1=1V, E1 gain=10
        let netlist = Netlist::parse(
            "VCVS amplifier
V1 in 0 1
R1 in 0 10k
E1 out 0 in 0 10
R2 out 0 1k
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("in").unwrap(), 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(solution.voltage("out").unwrap(), 10.0, epsilon = 1e-12);
    }

    #[test]
    fn test_vcvs_inverting_amplifier() {
        // Op-amp as high-gain VCVS in inverting configuration:
        // V1=1V → R1=1k → inv node → R2=2k → out
        // E1 out 0 0 inv 100000 (gain from 0 to inv, i.e., V(out) = -100000 * V(inv))
        // With ideal gain: V(out) ≈ -R2/R1 * V(in) = -2V
        let netlist = Netlist::parse(
            "Inverting amplifier
V1 in 0 1
R1 in inv 1k
R2 inv out 2k
E1 out 0 0 inv 100000
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        let solution = mna.solve().unwrap();
        // V(out) ≈ -2.0 (with high gain, very close to ideal)
        assert_abs_diff_eq!(solution.voltage("out").unwrap(), -2.0, epsilon = 1e-3);
        // V(inv) ≈ 0 (virtual ground)
        assert_abs_diff_eq!(solution.voltage("inv").unwrap(), 0.0, epsilon = 1e-3);
    }

    #[test]
    fn test_vccs_transconductance() {
        // VCCS: G1 drives current gm * V(in) into node "out" through R2
        // SPICE convention: current flows from out- to out+ externally,
        // so G1 0 out means current enters "out" from ground.
        // V1=2V at "in", G1 gm=1m, R2=1k load
        // V(out) = gm * V(in) * R2 = 1e-3 * 2 * 1000 = 2V
        let netlist = Netlist::parse(
            "VCCS test
V1 in 0 2
R1 in 0 10k
G1 0 out in 0 1m
R2 out 0 1k
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        // No branch for G — only V1 has a branch
        assert_eq!(mna.vsource_names.len(), 1);

        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("in").unwrap(), 2.0, epsilon = 1e-12);
        assert_abs_diff_eq!(solution.voltage("out").unwrap(), 2.0, epsilon = 1e-12);
    }

    #[test]
    fn test_cccs_current_amplifier() {
        // CCCS: F1 outputs gain * I(Vsense)
        // Vsense is a 0V sensing source in series with R1
        // V1=10V, R1=10k → I(Vsense) = 1mA (in MNA branch convention)
        // F1 0 out: current enters "out" from ground
        // F1 gain=5 → I_out = 5mA into R2=1k → V(out) = 5V
        let netlist = Netlist::parse(
            "CCCS test
V1 1 0 10
R1 1 sense 10k
Vsense sense 0 0
F1 0 out Vsense 5
R2 out 0 1k
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        let solution = mna.solve().unwrap();
        // I(Vsense) = 1mA, F1 gain=5 → 5mA into R2=1k
        assert_abs_diff_eq!(solution.voltage("out").unwrap(), 5.0, epsilon = 1e-9);
    }

    #[test]
    fn test_ccvs_transresistance() {
        // CCVS: H1 output voltage = rm * I(Vsense)
        // V1=5V, R1=5k → I(Vsense) = 1mA
        // H1 rm=2k → V(out) = 2k * 1mA = 2V
        let netlist = Netlist::parse(
            "CCVS test
V1 1 0 5
R1 1 sense 5k
Vsense sense 0 0
H1 out 0 Vsense 2k
R2 out 0 10k
.op
.end
",
        )
        .unwrap();

        let mna = assemble_mna(&netlist).unwrap();
        // V1, Vsense, H1 all have branches
        assert_eq!(mna.vsource_names.len(), 3);

        let solution = mna.solve().unwrap();
        assert_abs_diff_eq!(solution.voltage("out").unwrap(), 2.0, epsilon = 1e-9);
    }
}
