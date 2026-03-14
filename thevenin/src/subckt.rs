//! Subcircuit expansion (flattening).
//!
//! Transforms a hierarchical netlist with `.subckt` definitions and `X` instance
//! calls into a flat netlist where every element is at the top level.
//!
//! # Algorithm
//!
//! 1. Collect all `.subckt` definitions from the netlist (including nested ones).
//! 2. For each `X` element, look up the subcircuit definition.
//! 3. Map subcircuit port nodes to the caller's connection nodes.
//! 4. Internal nodes (non-port, non-ground) get prefixed with the instance path.
//! 5. Element names get prefixed with the instance path.
//! 6. Recursively expand any nested `X` calls.
//! 7. Parameter substitution: instance params override subcircuit defaults.

use std::collections::{BTreeMap, HashMap};

use thevenin_types::{Element, ElementKind, Expr, Item, Netlist, Param, SubcktDef};

/// Remap node names inside `v(...)` calls in a B-source spec string.
fn remap_v_calls_in_spec(spec: &str, remap: &impl Fn(&str) -> String) -> String {
    let mut result = String::with_capacity(spec.len() + 32);
    let mut pos = 0;

    while pos < spec.len() {
        // Find 'v(' case-insensitively
        let remaining = &spec[pos..];
        let found = remaining.char_indices().find(|&(i, c)| {
            (c == 'v' || c == 'V') && remaining[i + 1..].starts_with('(')
        });

        let (rel, ch) = match found {
            Some(x) => x,
            None => {
                result.push_str(&spec[pos..]);
                break;
            }
        };

        let v_start = pos + rel;

        // Check word boundary
        let is_word_start = v_start == 0 || {
            let prev = spec[..v_start].chars().last().unwrap_or(' ');
            !prev.is_alphanumeric() && prev != '_'
        };

        if is_word_start {
            let after_paren = v_start + 2;
            if let Some(rel_close) = spec[after_paren..].find(')') {
                let close = after_paren + rel_close;
                let content = &spec[after_paren..close];
                let parts: Vec<&str> = content.split(',').collect();
                let remapped: Vec<String> = parts.iter().map(|p| remap(p.trim())).collect();

                result.push_str(&spec[pos..v_start]);
                // Preserve the 'v' or 'V' case
                result.push(ch);
                result.push('(');
                result.push_str(&remapped.join(", "));
                result.push(')');
                pos = close + 1;
                continue;
            }
        }

        result.push_str(&spec[pos..v_start + ch.len_utf8()]);
        pos = v_start + ch.len_utf8();
    }

    result
}
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SubcktError {
    #[error("undefined subcircuit: {0}")]
    UndefinedSubckt(String),
    #[error("port count mismatch for subcircuit {subckt}: expected {expected}, got {got}")]
    PortMismatch {
        subckt: String,
        expected: usize,
        got: usize,
    },
    #[error("recursion limit exceeded expanding subcircuit {0}")]
    RecursionLimit(String),
}

const MAX_RECURSION_DEPTH: usize = 64;

/// Flatten a hierarchical netlist by expanding all subcircuit calls.
///
/// Returns a new netlist with no `X` elements and no `.subckt` definitions —
/// all elements are at the top level. Model definitions from subcircuits are
/// hoisted to the top level with prefixed names.
pub fn flatten_netlist(netlist: &Netlist) -> Result<Netlist, SubcktError> {
    // Collect only top-level subcircuit definitions.
    let mut subckt_defs: BTreeMap<String, SubcktDef> = BTreeMap::new();
    collect_subckt_defs(&netlist.items, &mut subckt_defs);

    // Collect global nodes.
    let global_nodes: Vec<String> = netlist
        .items
        .iter()
        .filter_map(|item| {
            if let Item::Global(nodes) = item {
                Some(nodes.clone())
            } else {
                None
            }
        })
        .flatten()
        .collect();

    // Build flattened items.
    let mut flat_items: Vec<Item> = Vec::new();

    for item in &netlist.items {
        match item {
            Item::Element(elem) => {
                if let ElementKind::SubcktCall {
                    ports,
                    subckt,
                    params,
                } = &elem.kind
                {
                    expand_subckt_call(
                        &elem.name,
                        ports,
                        subckt,
                        params,
                        &subckt_defs,
                        &global_nodes,
                        &mut flat_items,
                        0,
                    )?;
                } else {
                    flat_items.push(item.clone());
                }
            }
            // Drop subcircuit definitions — they've been collected
            Item::Subckt(_) => {}
            _ => flat_items.push(item.clone()),
        }
    }

    Ok(Netlist {
        title: netlist.title.clone(),
        source: netlist.source.clone(),
        items: flat_items,
    })
}

/// Collect subcircuit definitions from items (top-level only, not nested).
fn collect_subckt_defs(items: &[Item], defs: &mut BTreeMap<String, SubcktDef>) {
    for item in items {
        if let Item::Subckt(subckt) = item {
            defs.insert(subckt.name.to_uppercase(), subckt.clone());
        }
    }
}

/// Expand a single subcircuit call, recursively expanding nested calls.
#[expect(clippy::too_many_arguments)]
fn expand_subckt_call(
    instance_name: &str,
    ports: &[String],
    subckt_name: &str,
    instance_params: &[Param],
    subckt_defs: &BTreeMap<String, SubcktDef>,
    global_nodes: &[String],
    output: &mut Vec<Item>,
    depth: usize,
) -> Result<(), SubcktError> {
    if depth >= MAX_RECURSION_DEPTH {
        return Err(SubcktError::RecursionLimit(subckt_name.to_string()));
    }

    let subckt = subckt_defs
        .get(&subckt_name.to_uppercase())
        .ok_or_else(|| SubcktError::UndefinedSubckt(subckt_name.to_string()))?;

    // Collect locally-defined subcircuits within this subcircuit definition.
    // These take priority over the global defs for nested X calls.
    let mut local_defs = subckt_defs.clone();
    collect_subckt_defs(&subckt.items, &mut local_defs);

    if ports.len() != subckt.ports.len() {
        return Err(SubcktError::PortMismatch {
            subckt: subckt_name.to_string(),
            expected: subckt.ports.len(),
            got: ports.len(),
        });
    }

    // Build port mapping: subckt port name -> caller node name
    let mut node_map: HashMap<String, String> = HashMap::new();
    for (port, caller_node) in subckt.ports.iter().zip(ports.iter()) {
        node_map.insert(port.to_lowercase(), caller_node.clone());
    }

    // Build parameter map: subckt default params, then instance overrides
    let mut param_map: BTreeMap<String, f64> = BTreeMap::new();
    for p in &subckt.params {
        if let Expr::Num(v) = &p.value {
            param_map.insert(p.name.to_uppercase(), *v);
        }
    }
    for p in instance_params {
        match &p.value {
            Expr::Num(v) => {
                param_map.insert(p.name.to_uppercase(), *v);
            }
            Expr::Brace(s) => {
                // Evaluate expression in context of current params
                let mut ctx = crate::expr::EvalContext::default();
                for (k, &v) in &param_map {
                    ctx.params.insert(k.clone(), v);
                }
                if let Ok(val) = ctx.eval_str(s) {
                    param_map.insert(p.name.to_uppercase(), val);
                }
            }
            Expr::Param(name) => {
                if let Some(&v) = param_map.get(&name.to_uppercase()) {
                    param_map.insert(p.name.to_uppercase(), v);
                }
            }
        }
    }

    // Also collect .param definitions inside the subcircuit
    for item in &subckt.items {
        if let Item::Param(params) = item {
            for p in params {
                if let Expr::Num(v) = &p.value {
                    param_map.insert(p.name.to_uppercase(), *v);
                }
            }
        }
    }

    // Instance prefix for unique naming (lowercase, matching ngspice convention)
    let prefix = instance_name.to_lowercase();

    // Collect local model names for this subcircuit scope so we can
    // prefix them to avoid collisions with same-named models in other scopes.
    let local_model_names: std::collections::HashSet<String> = subckt
        .items
        .iter()
        .filter_map(|item| {
            if let Item::Model(m) = item {
                Some(m.name.to_uppercase())
            } else {
                None
            }
        })
        .collect();

    // Process each item in the subcircuit body
    for item in &subckt.items {
        match item {
            Item::Element(elem) => {
                if let ElementKind::SubcktCall {
                    ports: nested_ports,
                    subckt: nested_subckt,
                    params: nested_params,
                } = &elem.kind
                {
                    // Remap port nodes for the nested call
                    let remapped_ports: Vec<String> = nested_ports
                        .iter()
                        .map(|n| remap_node(n, &prefix, &node_map, global_nodes))
                        .collect();
                    let nested_instance_name = format!("{}.{}", prefix, elem.name.to_lowercase());
                    let resolved_params = resolve_params(nested_params, &param_map);
                    expand_subckt_call(
                        &nested_instance_name,
                        &remapped_ports,
                        nested_subckt,
                        &resolved_params,
                        &local_defs,
                        global_nodes,
                        output,
                        depth + 1,
                    )?;
                } else {
                    // Remap element nodes and name, including model references
                    let mut new_elem =
                        remap_element(elem, &prefix, &node_map, global_nodes, &param_map);
                    remap_model_ref(&mut new_elem.kind, &prefix, &local_model_names);
                    output.push(Item::Element(new_elem));
                }
            }
            Item::Model(model_def) => {
                // Hoist model definitions to top level with prefixed name
                // to maintain proper scoping
                let mut prefixed_model = model_def.clone();
                prefixed_model.name = format!("{}.{}", prefix, model_def.name);
                output.push(Item::Model(prefixed_model));
            }
            Item::Subckt(_) => {
                // Nested subcircuit definitions are already collected globally
            }
            Item::Param(_) => {
                // Parameters are consumed into param_map above
            }
            _ => {
                // Other items (analyses, comments) are dropped from subcircuit body
            }
        }
    }

    Ok(())
}

/// Remap a node name within a subcircuit expansion.
///
/// - Ground ("0") is never remapped.
/// - Port nodes are mapped to the caller's connection nodes.
/// - Global nodes are kept as-is.
/// - Internal nodes are prefixed with the instance path.
fn remap_node(
    node: &str,
    prefix: &str,
    port_map: &HashMap<String, String>,
    global_nodes: &[String],
) -> String {
    // Ground is universal
    if node == "0" {
        return "0".to_string();
    }

    // Check port mapping (case-insensitive)
    if let Some(caller_node) = port_map.get(&node.to_lowercase()) {
        return caller_node.clone();
    }

    // Check global nodes (case-insensitive)
    if global_nodes.iter().any(|g| g.eq_ignore_ascii_case(node)) {
        return node.to_string();
    }

    // Internal node: prefix with instance path
    format!("{}.{}", prefix, node.to_lowercase())
}

/// Remap all node references in an element.
fn remap_element(
    elem: &Element,
    prefix: &str,
    port_map: &HashMap<String, String>,
    global_nodes: &[String],
    param_map: &BTreeMap<String, f64>,
) -> Element {
    let remap = |n: &str| remap_node(n, prefix, port_map, global_nodes);
    let resolve_expr = |e: &Expr| resolve_expr(e, param_map);

    let new_name = format!("{}.{}", prefix, elem.name.to_lowercase());

    let new_kind = match &elem.kind {
        ElementKind::Resistor {
            pos,
            neg,
            value,
            params,
        } => ElementKind::Resistor {
            pos: remap(pos),
            neg: remap(neg),
            value: resolve_expr(value),
            params: resolve_params(params, param_map),
        },
        ElementKind::Capacitor {
            pos,
            neg,
            value,
            params,
        } => ElementKind::Capacitor {
            pos: remap(pos),
            neg: remap(neg),
            value: resolve_expr(value),
            params: resolve_params(params, param_map),
        },
        ElementKind::Inductor {
            pos,
            neg,
            value,
            params,
        } => ElementKind::Inductor {
            pos: remap(pos),
            neg: remap(neg),
            value: resolve_expr(value),
            params: resolve_params(params, param_map),
        },
        ElementKind::VoltageSource { pos, neg, source } => ElementKind::VoltageSource {
            pos: remap(pos),
            neg: remap(neg),
            source: resolve_source(source, param_map),
        },
        ElementKind::CurrentSource { pos, neg, source } => ElementKind::CurrentSource {
            pos: remap(pos),
            neg: remap(neg),
            source: resolve_source(source, param_map),
        },
        ElementKind::Diode {
            anode,
            cathode,
            model,
            params,
        } => ElementKind::Diode {
            anode: remap(anode),
            cathode: remap(cathode),
            model: model.clone(),
            params: resolve_params(params, param_map),
        },
        ElementKind::Bjt {
            c,
            b,
            e,
            substrate,
            model,
            params,
        } => ElementKind::Bjt {
            c: remap(c),
            b: remap(b),
            e: remap(e),
            substrate: substrate.as_ref().map(|s| remap(s)),
            model: model.clone(),
            params: resolve_params(params, param_map),
        },
        ElementKind::Mosfet {
            d,
            g,
            s,
            bulk,
            body,
            model,
            params,
        } => ElementKind::Mosfet {
            d: remap(d),
            g: remap(g),
            s: remap(s),
            bulk: remap(bulk),
            body: body.as_ref().map(|b| remap(b)),
            model: model.clone(),
            params: resolve_params(params, param_map),
        },
        ElementKind::Jfet {
            d,
            g,
            s,
            model,
            params,
        } => ElementKind::Jfet {
            d: remap(d),
            g: remap(g),
            s: remap(s),
            model: model.clone(),
            params: resolve_params(params, param_map),
        },
        ElementKind::Mesa {
            d,
            g,
            s,
            model,
            params,
        } => ElementKind::Mesa {
            d: remap(d),
            g: remap(g),
            s: remap(s),
            model: model.clone(),
            params: resolve_params(params, param_map),
        },
        ElementKind::MutualCoupling { l1, l2, coupling } => ElementKind::MutualCoupling {
            l1: format!("{}.{}", prefix, l1.to_lowercase()),
            l2: format!("{}.{}", prefix, l2.to_lowercase()),
            coupling: resolve_expr(coupling),
        },
        ElementKind::Vcvs {
            out_pos,
            out_neg,
            in_pos,
            in_neg,
            gain,
        } => ElementKind::Vcvs {
            out_pos: remap(out_pos),
            out_neg: remap(out_neg),
            in_pos: remap(in_pos),
            in_neg: remap(in_neg),
            gain: resolve_expr(gain),
        },
        ElementKind::Cccs {
            out_pos,
            out_neg,
            vsrc,
            gain,
        } => ElementKind::Cccs {
            out_pos: remap(out_pos),
            out_neg: remap(out_neg),
            vsrc: format!("{}.{}", prefix, vsrc.to_lowercase()),
            gain: resolve_expr(gain),
        },
        ElementKind::Vccs {
            out_pos,
            out_neg,
            in_pos,
            in_neg,
            gm,
        } => ElementKind::Vccs {
            out_pos: remap(out_pos),
            out_neg: remap(out_neg),
            in_pos: remap(in_pos),
            in_neg: remap(in_neg),
            gm: resolve_expr(gm),
        },
        ElementKind::Ccvs {
            out_pos,
            out_neg,
            vsrc,
            rm,
        } => ElementKind::Ccvs {
            out_pos: remap(out_pos),
            out_neg: remap(out_neg),
            vsrc: format!("{}.{}", prefix, vsrc.to_lowercase()),
            rm: resolve_expr(rm),
        },
        ElementKind::BehavioralSource { pos, neg, spec } => ElementKind::BehavioralSource {
            pos: remap(pos),
            neg: remap(neg),
            spec: remap_v_calls_in_spec(spec, &remap),
        },
        // SubcktCall should have been handled before calling this
        ElementKind::SubcktCall { .. } => {
            unreachable!("SubcktCall should be expanded before remap")
        }
        ElementKind::Ltra {
            pos1,
            neg1,
            pos2,
            neg2,
            model,
            params,
        } => ElementKind::Ltra {
            pos1: remap(pos1),
            neg1: remap(neg1),
            pos2: remap(pos2),
            neg2: remap(neg2),
            model: model.clone(),
            params: resolve_params(params, param_map),
        },
        ElementKind::Txl {
            pos1,
            neg1,
            pos2,
            neg2,
            model,
            params,
        } => ElementKind::Txl {
            pos1: remap(pos1),
            neg1: remap(neg1),
            pos2: remap(pos2),
            neg2: remap(neg2),
            model: model.clone(),
            params: resolve_params(params, param_map),
        },
        ElementKind::Cpl {
            in_nodes,
            out_nodes,
            gnd,
            model,
            params,
        } => ElementKind::Cpl {
            in_nodes: in_nodes.iter().map(|n| remap(n)).collect(),
            out_nodes: out_nodes.iter().map(|n| remap(n)).collect(),
            gnd: gnd.clone(),
            model: model.clone(),
            params: resolve_params(params, param_map),
        },
        ElementKind::Xspice { connections, model } => ElementKind::Xspice {
            connections: connections
                .iter()
                .map(|c| match c {
                    thevenin_types::XspiceConnection::Scalar(s) => {
                        thevenin_types::XspiceConnection::Scalar(remap(s))
                    }
                    thevenin_types::XspiceConnection::Array(ports) => {
                        thevenin_types::XspiceConnection::Array(
                            ports.iter().map(|p| remap(p)).collect(),
                        )
                    }
                })
                .collect(),
            model: model.clone(),
        },
        ElementKind::Raw(s) => ElementKind::Raw(s.clone()),
    };

    Element {
        name: new_name,
        kind: new_kind,
    }
}

/// Remap model references in an element if they match a locally-scoped model.
fn remap_model_ref(
    kind: &mut ElementKind,
    prefix: &str,
    local_model_names: &std::collections::HashSet<String>,
) {
    let remap = |model: &mut String| {
        if local_model_names.contains(&model.to_uppercase()) {
            *model = format!("{}.{}", prefix, model);
        }
    };

    match kind {
        ElementKind::Diode { model, .. } => remap(model),
        ElementKind::Bjt { model, .. } => remap(model),
        ElementKind::Mosfet { model, .. } => remap(model),
        ElementKind::Jfet { model, .. } => remap(model),
        ElementKind::Mesa { model, .. } => remap(model),
        // For resistors with model-referenced values (Expr::Param), prefix the param name
        ElementKind::Resistor { value, .. } => {
            if let Expr::Param(name) = value
                && local_model_names.contains(&name.to_uppercase())
            {
                *value = Expr::Param(format!("{}.{}", prefix, name));
            }
        }
        _ => {}
    }
}

/// Resolve a parameter expression using the parameter map.
fn resolve_expr(expr: &Expr, param_map: &BTreeMap<String, f64>) -> Expr {
    match expr {
        Expr::Num(_) => expr.clone(),
        Expr::Param(name) => {
            if let Some(&val) = param_map.get(&name.to_uppercase()) {
                Expr::Num(val)
            } else {
                expr.clone()
            }
        }
        Expr::Brace(s) => {
            // Try to evaluate using the expression evaluator with current params.
            let mut ctx = crate::expr::EvalContext::default();
            for (k, &v) in param_map {
                ctx.params.insert(k.clone(), v);
            }
            if let Ok(val) = ctx.eval_str(s) {
                Expr::Num(val)
            } else {
                expr.clone()
            }
        }
    }
}

/// Resolve parameters in a list.
fn resolve_params(params: &[Param], param_map: &BTreeMap<String, f64>) -> Vec<Param> {
    params
        .iter()
        .map(|p| Param {
            name: p.name.clone(),
            value: resolve_expr(&p.value, param_map),
        })
        .collect()
}

/// Resolve a Source's expressions.
fn resolve_source(
    source: &thevenin_types::Source,
    param_map: &BTreeMap<String, f64>,
) -> thevenin_types::Source {
    thevenin_types::Source {
        dc: source.dc.as_ref().map(|e| resolve_expr(e, param_map)),
        ac: source.ac.as_ref().map(|ac| thevenin_types::AcSpec {
            mag: resolve_expr(&ac.mag, param_map),
            phase: ac.phase.as_ref().map(|p| resolve_expr(p, param_map)),
        }),
        waveform: source.waveform.clone(), // waveform params could also be resolved
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use thevenin_types::Netlist;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    #[test]
    fn test_simple_subcircuit_expansion() {
        let netlist = Netlist::parse(
            "Subcircuit test
.subckt DIVIDER in out
R1 in mid 1k
R2 mid out 1k
.ends DIVIDER
V1 1 0 10
X1 1 0 DIVIDER
.op
.end
",
        )
        .unwrap();

        let flat = flatten_netlist(&netlist).unwrap();

        // Should have: V1, x1.r1, x1.r2, .op
        let elements: Vec<&Element> = flat
            .items
            .iter()
            .filter_map(|i| {
                if let Item::Element(e) = i {
                    Some(e)
                } else {
                    None
                }
            })
            .collect();

        assert_eq!(elements.len(), 3); // V1 + 2 from subcircuit
        assert_eq!(elements[0].name, "V1");
        assert_eq!(elements[1].name, "x1.r1");
        assert_eq!(elements[2].name, "x1.r2");

        // Verify node remapping
        if let ElementKind::Resistor { pos, neg, .. } = &elements[1].kind {
            assert_eq!(pos, "1"); // port "in" -> caller "1"
            assert_eq!(neg, "x1.mid"); // internal node
        } else {
            panic!("expected resistor");
        }

        if let ElementKind::Resistor { pos, neg, .. } = &elements[2].kind {
            assert_eq!(pos, "x1.mid"); // internal node
            assert_eq!(neg, "0"); // port "out" -> caller "0" (ground)
        } else {
            panic!("expected resistor");
        }
    }

    #[test]
    fn test_subcircuit_with_params() {
        let netlist = Netlist::parse(
            "Param test
.subckt RLOAD p n PARAMS: rval=1k
R1 p n rval
.ends RLOAD
V1 1 0 5
X1 1 0 RLOAD PARAMS: rval=2k
.op
.end
",
        )
        .unwrap();

        let flat = flatten_netlist(&netlist).unwrap();

        let elements: Vec<&Element> = flat
            .items
            .iter()
            .filter_map(|i| {
                if let Item::Element(e) = i {
                    Some(e)
                } else {
                    None
                }
            })
            .collect();

        // R1 inside subcircuit should have value 2k (overridden by instance)
        if let ElementKind::Resistor { value, .. } = &elements[1].kind {
            if let Expr::Num(v) = value {
                assert!((*v - 2000.0).abs() < 1e-9, "expected 2k, got {}", v);
            } else {
                panic!("expected numeric value for resistor");
            }
        } else {
            panic!("expected resistor");
        }
    }

    #[test]
    fn test_nested_subcircuit() {
        let netlist = Netlist::parse(
            "Nested subcircuit test
.subckt INNER a b
R1 a b 1k
.ends INNER
.subckt OUTER a b
X1 a mid INNER
X2 mid b INNER
.ends OUTER
V1 in 0 10
X1 in 0 OUTER
.op
.end
",
        )
        .unwrap();

        let flat = flatten_netlist(&netlist).unwrap();

        let elements: Vec<&Element> = flat
            .items
            .iter()
            .filter_map(|i| {
                if let Item::Element(e) = i {
                    Some(e)
                } else {
                    None
                }
            })
            .collect();

        // V1 + x1.x1.r1 + x1.x2.r1
        assert_eq!(elements.len(), 3);
        assert_eq!(elements[1].name, "x1.x1.r1");
        assert_eq!(elements[2].name, "x1.x2.r1");

        // x1.x1.r1: INNER called with (in, x1.mid) -> R1 from "in" to "x1.mid"
        if let ElementKind::Resistor { pos, neg, .. } = &elements[1].kind {
            assert_eq!(pos, "in");
            assert_eq!(neg, "x1.mid");
        } else {
            panic!("expected resistor");
        }

        // x1.x2.r1: INNER called with (x1.mid, 0) -> R1 from "x1.mid" to "0"
        if let ElementKind::Resistor { pos, neg, .. } = &elements[2].kind {
            assert_eq!(pos, "x1.mid");
            assert_eq!(neg, "0");
        } else {
            panic!("expected resistor");
        }
    }

    #[test]
    fn test_multiple_instances() {
        let netlist = Netlist::parse(
            "Multiple instances
.subckt BUF in out
R1 in out 100
.ends BUF
V1 1 0 5
X1 1 2 BUF
X2 2 3 BUF
R1 3 0 1k
.op
.end
",
        )
        .unwrap();

        let flat = flatten_netlist(&netlist).unwrap();

        let elements: Vec<&Element> = flat
            .items
            .iter()
            .filter_map(|i| {
                if let Item::Element(e) = i {
                    Some(e)
                } else {
                    None
                }
            })
            .collect();

        // V1 + x1.r1 + x2.r1 + R1
        assert_eq!(elements.len(), 4);

        // x1.r1: in=1, out=2
        if let ElementKind::Resistor { pos, neg, .. } = &elements[1].kind {
            assert_eq!(pos, "1");
            assert_eq!(neg, "2");
        }
        // x2.r1: in=2, out=3
        if let ElementKind::Resistor { pos, neg, .. } = &elements[2].kind {
            assert_eq!(pos, "2");
            assert_eq!(neg, "3");
        }
    }

    #[test]
    fn test_undefined_subcircuit_error() {
        let netlist = Netlist::parse(
            "Error test
X1 1 0 NONEXISTENT
.op
.end
",
        )
        .unwrap();

        let result = flatten_netlist(&netlist);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SubcktError::UndefinedSubckt(_)
        ));
    }

    #[test]
    fn test_port_count_mismatch_error() {
        let netlist = Netlist::parse(
            "Port mismatch test
.subckt BUF in out
R1 in out 100
.ends BUF
X1 1 2 3 BUF
.op
.end
",
        )
        .unwrap();

        let result = flatten_netlist(&netlist);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SubcktError::PortMismatch { .. }
        ));
    }

    #[test]
    fn test_subcircuit_with_model() {
        let netlist = Netlist::parse(
            "Subcircuit with model
.subckt DCLAMP a k
.model DMOD D IS=1e-14
D1 a k DMOD
.ends DCLAMP
V1 1 0 5
X1 1 0 DCLAMP
.op
.end
",
        )
        .unwrap();

        let flat = flatten_netlist(&netlist).unwrap();

        // Check that model definition was hoisted
        let models: Vec<_> = flat
            .items
            .iter()
            .filter_map(|i| {
                if let Item::Model(m) = i {
                    Some(m)
                } else {
                    None
                }
            })
            .collect();
        assert_eq!(models.len(), 1);
        assert_eq!(models[0].name, "x1.DMOD");

        // Check diode element was expanded
        let elements: Vec<&Element> = flat
            .items
            .iter()
            .filter_map(|i| {
                if let Item::Element(e) = i {
                    Some(e)
                } else {
                    None
                }
            })
            .collect();
        assert_eq!(elements.len(), 2); // V1 + x1.d1

        if let ElementKind::Diode {
            anode,
            cathode,
            model,
            ..
        } = &elements[1].kind
        {
            assert_eq!(anode, "1"); // port "a" -> "1"
            assert_eq!(cathode, "0"); // port "k" -> "0"
            assert_eq!(model, "x1.DMOD"); // model name prefixed with instance path
        } else {
            panic!("expected diode");
        }
    }

    #[test]
    fn test_no_subcircuits_passthrough() {
        let netlist = Netlist::parse(
            "No subcircuits
V1 1 0 5
R1 1 0 1k
.op
.end
",
        )
        .unwrap();

        let flat = flatten_netlist(&netlist).unwrap();

        let elements: Vec<&Element> = flat
            .items
            .iter()
            .filter_map(|i| {
                if let Item::Element(e) = i {
                    Some(e)
                } else {
                    None
                }
            })
            .collect();

        assert_eq!(elements.len(), 2);
        assert_eq!(elements[0].name, "V1");
        assert_eq!(elements[1].name, "R1");
    }
}
