//! SPICE netlist parser.
//!
//! # Pre-processing
//! Before any structural parsing the raw text is normalised:
//! 1. Lines starting with `+` are appended to the previous logical line.
//! 2. Inline comments (`$` or `;` outside parentheses) are stripped.
//! 3. The first non-empty logical line becomes the title.
//!
//! # Value / suffix rules
//! | Suffix (case-insensitive) | Multiplier |
//! |---------------------------|-----------|
//! | `T`                       | 1 × 10¹²  |
//! | `G`                       | 1 × 10⁹   |
//! | `MEG`                     | 1 × 10⁶   |
//! | `K`                       | 1 × 10³   |
//! | `M`                       | 1 × 10⁻³  ← milli, NOT mega |
//! | `U`                       | 1 × 10⁻⁶  |
//! | `N`                       | 1 × 10⁻⁹  |
//! | `P`                       | 1 × 10⁻¹² |
//! | `F`                       | 1 × 10⁻¹⁵ |
//! | `A`                       | 1 × 10⁻¹⁸ |
//! | `MIL`                     | 25.4 × 10⁻⁶ |

use crate::{
    AcSpec, AcVariation, Analysis, DcSweep, Element, ElementKind, Expr, Item, ModelDef, Netlist,
    Param, PwlPoint, PzAnalysisType, PzInputType, Source, SubcktDef, Waveform, XspiceConnection,
};

// ---------------------------------------------------------------------------
// Error type
// ---------------------------------------------------------------------------

#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("netlist is empty")]
    Empty,
    #[error("line {line}: {msg}")]
    Syntax { line: usize, msg: String },
}

fn syntax(line: usize, msg: impl Into<String>) -> ParseError {
    ParseError::Syntax {
        line,
        msg: msg.into(),
    }
}

// ---------------------------------------------------------------------------
// Pre-processing: fold continuations, strip inline comments
// ---------------------------------------------------------------------------

/// Strip everything from the first `$` or `;` that appears outside
/// parentheses (those are inline SPICE comments).
fn strip_inline_comment(s: &str) -> &str {
    let mut depth: i32 = 0;
    for (i, c) in s.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => depth -= 1,
            '$' | ';' if depth == 0 => return &s[..i],
            _ => {}
        }
    }
    s
}

/// Convert raw text into logical lines (continuations folded, comments
/// stripped, blank lines discarded).
fn preprocess(input: &str) -> Vec<(usize, String)> {
    // (original_line_number, logical_line_text)
    let mut out: Vec<(usize, String)> = Vec::new();

    // Index into `out` of the last non-comment logical line, for attaching
    // continuation lines past comment lines.
    let mut last_noncomment: Option<usize> = None;

    for (lineno, raw) in input.lines().enumerate() {
        let line = strip_inline_comment(raw).trim();
        if line.is_empty() {
            continue;
        }
        if let Some(stripped) = line.strip_prefix('+') {
            // Continuation: append to the last non-comment logical line.
            // Standard SPICE allows comment lines between a statement and its
            // continuation lines (e.g. `*+ commented-out param` before `+ real=val`).
            if let Some(idx) = last_noncomment {
                out[idx].1.push(' ');
                out[idx].1.push_str(stripped.trim());
            }
            // If there's no previous non-comment line, silently drop it.
        } else if line.starts_with('*') {
            // Full-line comment: keep in output (produces Item::Comment) but
            // do NOT update last_noncomment so that subsequent '+' lines
            // still attach to the previous non-comment line.
            out.push((lineno + 1, line.to_string()));
        } else {
            last_noncomment = Some(out.len());
            out.push((lineno + 1, line.to_string()));
        }
    }

    out
}

/// Evaluate a simple numeric condition used in `.if` / `.elseif` directives.
///
/// Supports: `(expr == expr)`, `(expr != expr)`, `(expr > expr)`, etc.
/// where each expr is either a numeric literal or a parameter name.
fn eval_condition(cond: &str, params: &std::collections::HashMap<String, f64>) -> bool {
    let cond = cond.trim();
    // Strip outer parentheses
    let cond = if cond.starts_with('(') && cond.ends_with(')') {
        &cond[1..cond.len() - 1]
    } else {
        cond
    };
    let cond = cond.trim();

    // Find comparison operator
    let (lhs, op, rhs) = if let Some(i) = cond.find("==") {
        (&cond[..i], "==", &cond[i + 2..])
    } else if let Some(i) = cond.find("!=") {
        (&cond[..i], "!=", &cond[i + 2..])
    } else if let Some(i) = cond.find(">=") {
        (&cond[..i], ">=", &cond[i + 2..])
    } else if let Some(i) = cond.find("<=") {
        (&cond[..i], "<=", &cond[i + 2..])
    } else if let Some(i) = cond.find('>') {
        (&cond[..i], ">", &cond[i + 1..])
    } else if let Some(i) = cond.find('<') {
        (&cond[..i], "<", &cond[i + 1..])
    } else if let Some(i) = cond.find('=') {
        // Single = also treated as == (ngspice compat)
        (&cond[..i], "==", &cond[i + 1..])
    } else {
        // If no operator, treat as truthy: non-zero = true
        let val = resolve_param_val(cond.trim(), params);
        return val != 0.0;
    };

    let l = resolve_param_val(lhs.trim(), params);
    let r = resolve_param_val(rhs.trim(), params);

    match op {
        "==" => (l - r).abs() < 1e-15,
        "!=" => (l - r).abs() >= 1e-15,
        ">" => l > r,
        "<" => l < r,
        ">=" => l >= r,
        "<=" => l <= r,
        _ => false,
    }
}

fn resolve_param_val(tok: &str, params: &std::collections::HashMap<String, f64>) -> f64 {
    if let Some(v) = parse_spice_number(tok) {
        v
    } else if let Some(&v) = params.get(&tok.to_lowercase()) {
        v
    } else {
        0.0
    }
}

/// Process `.if` / `.elseif` / `.else` / `.endif` conditional blocks.
///
/// This runs on preprocessed logical lines, collecting `.param` definitions
/// as it goes, then evaluating conditions to include/exclude lines.
fn process_conditionals(lines: Vec<(usize, String)>) -> Vec<(usize, String)> {
    let mut params: std::collections::HashMap<String, f64> = std::collections::HashMap::new();
    let mut result = Vec::new();

    // Stack of (branch_taken, any_branch_taken) for nested .if blocks
    let mut stack: Vec<(bool, bool)> = Vec::new();

    for (lineno, line) in lines {
        let upper = line.trim().to_uppercase();
        let tokens: Vec<&str> = upper.split_whitespace().collect();

        if tokens.is_empty() {
            continue;
        }

        match tokens[0] {
            ".IF" => {
                let cond_str = line.trim()[3..].trim();
                let active = stack.iter().all(|(taken, _)| *taken);
                let cond_result = if active {
                    eval_condition(cond_str, &params)
                } else {
                    false
                };
                stack.push((active && cond_result, active && cond_result));
            }
            ".ELSEIF" => {
                let parent_active =
                    stack.len() <= 1 || stack[..stack.len() - 1].iter().all(|(t, _)| *t);
                if let Some((taken, any_taken)) = stack.last_mut() {
                    if parent_active && !*any_taken {
                        let cond_str = line.trim()[7..].trim();
                        let cond_result = eval_condition(cond_str, &params);
                        *taken = cond_result;
                        if cond_result {
                            *any_taken = true;
                        }
                    } else {
                        *taken = false;
                    }
                }
            }
            ".ELSE" => {
                let parent_active =
                    stack.len() <= 1 || stack[..stack.len() - 1].iter().all(|(t, _)| *t);
                if let Some((taken, any_taken)) = stack.last_mut() {
                    if parent_active && !*any_taken {
                        *taken = true;
                        *any_taken = true;
                    } else {
                        *taken = false;
                    }
                }
            }
            ".ENDIF" => {
                stack.pop();
            }
            _ => {
                let active = stack.is_empty() || stack.iter().all(|(taken, _)| *taken);
                if active {
                    // Collect .param definitions for future condition evaluation
                    if tokens[0] == ".PARAM" || tokens[0] == ".PARAMETERS" {
                        for tok in line.split_whitespace().skip(1) {
                            if let Some((k, v)) = tok.split_once('=') {
                                let k = k.trim().to_lowercase();
                                let v = v.trim();
                                if let Some(val) = parse_spice_number(v) {
                                    params.insert(k, val);
                                }
                            }
                        }
                    }
                    result.push((lineno, line));
                }
            }
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Tokeniser
// ---------------------------------------------------------------------------

/// Split a logical line into tokens on whitespace, but keep parenthesised
/// groups (waveform args like `PULSE(0 1 1n)`) as single tokens.
fn tokenize(line: &str) -> Vec<String> {
    let mut tokens: Vec<String> = Vec::new();
    let mut current = String::new();
    let mut depth: i32 = 0;
    let mut in_squote = false;
    let mut in_brace = false;

    for c in line.chars() {
        match c {
            '\'' if !in_brace => {
                in_squote = !in_squote;
                current.push(c);
            }
            '{' if !in_squote => {
                in_brace = true;
                current.push(c);
            }
            '}' if in_brace => {
                in_brace = false;
                current.push(c);
            }
            '(' if !in_squote && !in_brace => {
                depth += 1;
                current.push(c);
            }
            ')' if !in_squote && !in_brace => {
                depth -= 1;
                current.push(c);
            }
            ' ' | '\t' if depth == 0 && !in_squote && !in_brace => {
                if !current.is_empty() {
                    tokens.push(current.clone());
                    current.clear();
                }
            }
            _ => current.push(c),
        }
    }
    if !current.is_empty() {
        tokens.push(current);
    }
    tokens
}

// ---------------------------------------------------------------------------
// Value / expression parsing
// ---------------------------------------------------------------------------

/// Parse a SPICE numeric literal (possibly with SI suffix) into `f64`.
///
/// Returns `None` if the token does not begin with a digit, sign, or `.`.
pub fn parse_spice_number(s: &str) -> Option<f64> {
    let s = s.trim();
    if s.is_empty() {
        return None;
    }
    let bytes = s.as_bytes();
    let mut i = 0;

    // Optional sign
    if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
        i += 1;
    }

    // Must have at least one digit or start with '.'
    if i >= bytes.len() {
        return None;
    }
    let starts_with_digit = bytes[i].is_ascii_digit() || bytes[i] == b'.';
    if !starts_with_digit {
        return None;
    }

    // Integer digits
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        i += 1;
    }
    // Decimal part
    if i < bytes.len() && bytes[i] == b'.' {
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
    }
    // Exponent: only consume if followed by optional sign + digit(s)
    let exp_start = i;
    if i < bytes.len() && (bytes[i] == b'e' || bytes[i] == b'E') {
        let j = i + 1;
        let k = if j < bytes.len() && (bytes[j] == b'+' || bytes[j] == b'-') {
            j + 1
        } else {
            j
        };
        if k < bytes.len() && bytes[k].is_ascii_digit() {
            i = k;
            while i < bytes.len() && bytes[i].is_ascii_digit() {
                i += 1;
            }
        } else {
            i = exp_start; // 'e' is a suffix here (atto? unlikely — be safe)
        }
    }

    let base: f64 = s[..i].parse().ok()?;
    let suffix = s[i..].trim_start();
    let multiplier = parse_si_suffix(suffix);
    Some(base * multiplier)
}

fn parse_si_suffix(s: &str) -> f64 {
    let u = s.to_uppercase();
    let u = u.trim_start();
    // Check longer prefixes first to avoid `M` shadowing `MEG`/`MIL`
    if u.starts_with("MEG") {
        1e6
    } else if u.starts_with("MIL") {
        25.4e-6
    } else if u.starts_with('T') {
        1e12
    } else if u.starts_with('G') {
        1e9
    } else if u.starts_with('K') {
        1e3
    } else if u.starts_with('M') {
        1e-3 // milli — NOT mega
    } else if u.starts_with('U') {
        1e-6
    } else if u.starts_with('N') {
        1e-9
    } else if u.starts_with('P') {
        1e-12
    } else if u.starts_with('F') {
        1e-15
    } else if u.starts_with('A') {
        1e-18
    } else {
        1.0
    }
}

/// Parse a token as an [`Expr`].
///
/// * `{...}` → `Expr::Brace`
/// * `'...'` → `Expr::Brace` (ngspice single-quote expressions)
/// * starts with digit / sign / `.` → `Expr::Num`
/// * otherwise → `Expr::Param` (a parameter name or node)
fn parse_expr(tok: &str) -> Expr {
    if let Some(inner) = tok.strip_prefix('{').and_then(|s| s.strip_suffix('}')) {
        return Expr::Brace(inner.to_string());
    }
    if let Some(inner) = tok.strip_prefix('\'').and_then(|s| s.strip_suffix('\'')) {
        return Expr::Brace(inner.to_string());
    }
    if let Some(v) = parse_spice_number(tok) {
        return Expr::Num(v);
    }
    Expr::Param(tok.to_string())
}

/// Split a `"key=value"` token into a [`Param`].
fn parse_kv(tok: &str) -> Option<Param> {
    let (k, v) = tok.split_once('=')?;
    Some(Param {
        name: k.trim().to_string(),
        value: parse_expr(v.trim()),
    })
}

/// Collect parameters from tokens: any token containing `=` (or following a
/// bare `PARAMS:` / `PARAM:` keyword) is treated as a key=value pair.
fn collect_params(tokens: &[String]) -> Vec<Param> {
    let mut params = Vec::new();
    let mut i = 0;
    while i < tokens.len() {
        let upper = tokens[i].to_uppercase();
        if upper == "PARAMS:" || upper == "PARAM:" {
            i += 1;
            continue;
        }
        if tokens[i].contains('=') {
            if let Some(kv) = parse_kv(&tokens[i]) {
                params.push(kv);
            }
        } else if i + 1 < tokens.len() && tokens[i + 1].starts_with('=') {
            // Handle `key = value` split across tokens
            let key = tokens[i].clone();
            let val_tok = if tokens[i + 1] == "=" && i + 2 < tokens.len() {
                i += 2;
                &tokens[i]
            } else {
                // `key =value`
                i += 1;
                &tokens[i][1..]
            };
            params.push(Param {
                name: key,
                value: parse_expr(val_tok),
            });
        }
        i += 1;
    }
    params
}

// ---------------------------------------------------------------------------
// Waveform parsing
// ---------------------------------------------------------------------------

fn parse_waveform(tok: &str) -> Option<Waveform> {
    // tok looks like `PULSE(0 1 1n ...)` or `SIN(...)` etc.
    let upper = tok.to_uppercase();

    let find_inner = |prefix: &str| -> Option<&str> {
        if upper.starts_with(prefix) && tok.ends_with(')') {
            Some(&tok[prefix.len()..tok.len() - 1])
        } else {
            None
        }
    };

    let args_of = |prefix: &str| -> Option<Vec<Expr>> {
        let inner = find_inner(prefix)?;
        Some(inner.split_whitespace().map(parse_expr).collect())
    };

    if let Some(args) = args_of("PULSE(") {
        let g = |i| args.get(i).cloned();
        Some(Waveform::Pulse {
            v1: g(0)?,
            v2: g(1)?,
            td: g(2),
            tr: g(3),
            tf: g(4),
            pw: g(5),
            per: g(6),
        })
    } else if let Some(args) = args_of("SIN(") {
        let g = |i| args.get(i).cloned();
        Some(Waveform::Sin {
            v0: g(0)?,
            va: g(1)?,
            freq: g(2),
            td: g(3),
            theta: g(4),
            phi: g(5),
        })
    } else if let Some(args) = args_of("EXP(") {
        let g = |i| args.get(i).cloned();
        Some(Waveform::Exp {
            v1: g(0)?,
            v2: g(1)?,
            td1: g(2),
            tau1: g(3),
            td2: g(4),
            tau2: g(5),
        })
    } else if let Some(inner) = find_inner("PWL(") {
        let nums: Vec<Expr> = inner.split_whitespace().map(parse_expr).collect();
        let points = nums
            .chunks(2)
            .filter_map(|c| {
                if c.len() == 2 {
                    Some(PwlPoint {
                        time: c[0].clone(),
                        value: c[1].clone(),
                    })
                } else {
                    None
                }
            })
            .collect();
        Some(Waveform::Pwl(points))
    } else if let Some(args) = args_of("SFFM(") {
        let g = |i| args.get(i).cloned();
        Some(Waveform::Sffm {
            v0: g(0)?,
            va: g(1)?,
            fc: g(2),
            fs: g(3),
            md: g(4),
        })
    } else if let Some(args) = args_of("AM(") {
        let g = |i| args.get(i).cloned();
        Some(Waveform::Am {
            va: g(0)?,
            vo: g(1)?,
            fc: g(2)?,
            fs: g(3)?,
            td: g(4),
        })
    } else {
        None
    }
}

// ---------------------------------------------------------------------------
// Source parsing (V / I sources)
// ---------------------------------------------------------------------------

/// Parse everything after `name n+ n-` for a V or I source.
fn parse_source(tokens: &[String]) -> Source {
    let mut src = Source::default();
    let mut i = 0;
    while i < tokens.len() {
        let upper = tokens[i].to_uppercase();
        match upper.as_str() {
            "DC" => {
                i += 1;
                if i < tokens.len() {
                    src.dc = Some(parse_expr(&tokens[i]));
                }
            }
            "AC" => {
                i += 1;
                let mag = if i < tokens.len() {
                    parse_expr(&tokens[i])
                } else {
                    Expr::Num(1.0)
                };
                let phase = if i + 1 < tokens.len() && parse_spice_number(&tokens[i + 1]).is_some()
                {
                    i += 1;
                    Some(parse_expr(&tokens[i]))
                } else {
                    None
                };
                src.ac = Some(AcSpec { mag, phase });
            }
            _ => {
                // Handle dc=VALUE or DC=VALUE format
                if upper.starts_with("DC=") {
                    let val_str = &tokens[i][3..]; // skip "dc=" or "DC="
                    if let Some(v) = parse_spice_number(val_str) {
                        src.dc = Some(Expr::Num(v));
                    } else {
                        src.dc = Some(parse_expr(val_str));
                    }
                } else if upper.starts_with("AC=") {
                    let val_str = &tokens[i][3..]; // skip "ac=" or "AC="
                    let mag = parse_expr(val_str);
                    src.ac = Some(AcSpec { mag, phase: None });
                } else if let Some(w) = parse_waveform(&tokens[i]) {
                    // Try as a waveform
                    src.waveform = Some(w);
                } else if src.dc.is_none() {
                    // Bare value with no keyword: treat as DC
                    if let Some(v) = parse_spice_number(&tokens[i]) {
                        src.dc = Some(Expr::Num(v));
                    } else if parse_expr(&tokens[i]) != Expr::Param(tokens[i].to_string()) {
                        // It's a brace expr
                        src.dc = Some(parse_expr(&tokens[i]));
                    }
                }
            }
        }
        i += 1;
    }
    src
}

/// Collect model parameters: key=value pairs and bare flags.
/// Bare tokens (without `=`) are treated as flag parameters with value 1.0.
fn collect_model_params(tokens: &[String]) -> Vec<Param> {
    let mut params = Vec::new();
    let mut current_key = String::new();
    let mut i = 0;
    while i < tokens.len() {
        let upper = tokens[i].to_uppercase();
        if upper == "PARAMS:" || upper == "PARAM:" {
            i += 1;
            continue;
        }
        if tokens[i].contains('=') {
            if let Some(kv) = parse_kv(&tokens[i]) {
                current_key = kv.name.clone();
                params.push(kv);
            }
        } else if i + 1 < tokens.len() && tokens[i + 1].starts_with('=') {
            // Handle `key = value` split across tokens
            let key = tokens[i].clone();
            let val_tok = if tokens[i + 1] == "=" && i + 2 < tokens.len() {
                i += 2;
                &tokens[i]
            } else {
                // `key =value`
                i += 1;
                &tokens[i][1..]
            };
            current_key = key.clone();
            params.push(Param {
                name: key,
                value: parse_expr(val_tok),
            });
        } else if !current_key.is_empty() && parse_spice_number(&tokens[i]).is_some() {
            // Bare numeric token after a keyed param — treat as continuation value
            // (used by CPL model matrix params: R=0.3 0 0 0.3 ...)
            params.push(Param {
                name: current_key.clone(),
                value: parse_expr(&tokens[i]),
            });
        } else {
            // Bare non-numeric token or no prior key — treat as flag param with value 1.0
            current_key.clear();
            params.push(Param {
                name: tokens[i].clone(),
                value: Expr::Num(1.0),
            });
        }
        i += 1;
    }
    params
}

// ---------------------------------------------------------------------------
// XSPICE connection parsing
// ---------------------------------------------------------------------------

/// Parse XSPICE connections from tokens after the element name.
///
/// Tokens may contain bracket groups like `[a1 a2 a3]` which become
/// `XspiceConnection::Array`. Everything else is a scalar port. The last
/// scalar token is the model name and is returned separately.
fn parse_xspice_connections(
    lineno: usize,
    tokens: &[String],
) -> Result<(Vec<XspiceConnection>, String), ParseError> {
    let mut connections: Vec<XspiceConnection> = Vec::new();
    let mut i = 0;

    while i < tokens.len() {
        let tok = &tokens[i];
        if tok.starts_with('[') && tok.ends_with(']') && tok.len() > 2 {
            // Single-token bracket group: [cs] or [a1 a2 a3] wouldn't be single
            // but [port] with one element is.
            let inner = &tok[1..tok.len() - 1];
            let ports: Vec<String> = inner.split_whitespace().map(String::from).collect();
            connections.push(XspiceConnection::Array(ports));
            i += 1;
        } else if let Some(after_bracket) = tok.strip_prefix('[') {
            // Multi-token bracket group: collect until we find a token ending with `]`
            let first = if !after_bracket.is_empty() {
                vec![after_bracket.to_string()]
            } else {
                vec![]
            };
            let mut ports = first;
            i += 1;
            loop {
                if i >= tokens.len() {
                    return Err(syntax(lineno, "A: unclosed bracket in XSPICE connections"));
                }
                let t = &tokens[i];
                if t.ends_with(']') {
                    let inner = &t[..t.len() - 1];
                    if !inner.is_empty() {
                        ports.push(inner.to_string());
                    }
                    i += 1;
                    break;
                } else {
                    ports.push(t.clone());
                    i += 1;
                }
            }
            connections.push(XspiceConnection::Array(ports));
        } else {
            // Scalar port (or model name — we'll figure out which at the end)
            connections.push(XspiceConnection::Scalar(tok.clone()));
            i += 1;
        }
    }

    // The last element must be a scalar — that's the model name.
    match connections.pop() {
        Some(XspiceConnection::Scalar(model)) => Ok((connections, model)),
        Some(other) => {
            // Put it back, missing model
            connections.push(other);
            Err(syntax(
                lineno,
                "A: last XSPICE token must be scalar model name",
            ))
        }
        None => Err(syntax(lineno, "A: missing XSPICE connections and model")),
    }
}

// ---------------------------------------------------------------------------
// Element line parsing
// ---------------------------------------------------------------------------

fn parse_element(lineno: usize, line: &str) -> Result<Element, ParseError> {
    let tokens = tokenize(line);
    if tokens.is_empty() {
        return Err(syntax(lineno, "empty element line"));
    }
    let name = tokens[0].clone();
    let letter = name.chars().next().unwrap().to_ascii_uppercase();
    let rest = &tokens[1..];

    macro_rules! need {
        ($idx:expr, $what:expr) => {
            rest.get($idx)
                .ok_or_else(|| syntax(lineno, format!("missing {}", $what)))?
                .as_str()
        };
    }

    let kind = match letter {
        'R' | 'C' | 'L' => {
            let pos = need!(0, "n+").to_string();
            let neg = need!(1, "n-").to_string();
            let value = parse_expr(need!(2, "value"));
            let params: Vec<Param> = rest[3..]
                .iter()
                .filter(|t| t.contains('='))
                .filter_map(|t| parse_kv(t))
                .collect();
            match letter {
                'R' => ElementKind::Resistor {
                    pos,
                    neg,
                    value,
                    params,
                },
                'C' => ElementKind::Capacitor {
                    pos,
                    neg,
                    value,
                    params,
                },
                _ => ElementKind::Inductor {
                    pos,
                    neg,
                    value,
                    params,
                },
            }
        }
        'V' => {
            let pos = need!(0, "n+").to_string();
            let neg = need!(1, "n-").to_string();
            let source = parse_source(&rest[2..]);
            ElementKind::VoltageSource { pos, neg, source }
        }
        'I' => {
            let pos = need!(0, "n+").to_string();
            let neg = need!(1, "n-").to_string();
            let source = parse_source(&rest[2..]);
            ElementKind::CurrentSource { pos, neg, source }
        }
        'D' => {
            let anode = need!(0, "anode").to_string();
            let cathode = need!(1, "cathode").to_string();
            let model = need!(2, "model").to_string();
            let params: Vec<_> = rest[3..].iter().filter_map(|t| parse_kv(t)).collect();
            ElementKind::Diode {
                anode,
                cathode,
                model,
                params,
            }
        }
        'Q' => {
            // Q c b e [substrate] model [params]
            // We detect the model by looking for the first token that is not
            // a known node pattern — heuristic: if there are 5 positional
            // tokens before params, the 4th is substrate.
            let c = need!(0, "collector").to_string();
            let b_node = need!(1, "base").to_string();
            let e = need!(2, "emitter").to_string();
            let positional: Vec<&str> = rest[3..]
                .iter()
                .filter(|t| !t.contains('='))
                .map(|s| s.as_str())
                .collect();
            let params: Vec<_> = rest[3..].iter().filter_map(|t| parse_kv(t)).collect();
            let (substrate, model) = if positional.len() >= 2 {
                (Some(positional[0].to_string()), positional[1].to_string())
            } else {
                (
                    None,
                    positional
                        .first()
                        .map(|s| s.to_string())
                        .ok_or_else(|| syntax(lineno, "Q: missing model"))?,
                )
            };
            ElementKind::Bjt {
                c,
                b: b_node,
                e,
                substrate,
                model,
                params,
            }
        }
        'M' => {
            // M d g s bulk [body] model [params]
            // Detect 4 vs 5 terminal: count positional (non-kv) tokens after d,g,s.
            let d = need!(0, "drain").to_string();
            let g = need!(1, "gate").to_string();
            let s = need!(2, "source").to_string();
            let positional: Vec<&str> = rest[3..]
                .iter()
                .filter(|t| !t.contains('='))
                .map(|s| s.as_str())
                .collect();
            let params: Vec<_> = rest[3..].iter().filter_map(|t| parse_kv(t)).collect();
            let (bulk, body, model) = if positional.len() >= 3 {
                // 5-terminal: bulk body model
                (
                    positional[0].to_string(),
                    Some(positional[1].to_string()),
                    positional[2].to_string(),
                )
            } else if positional.len() >= 2 {
                // 4-terminal: bulk model
                (positional[0].to_string(), None, positional[1].to_string())
            } else {
                return Err(syntax(lineno, "M: missing bulk/model"));
            };
            ElementKind::Mosfet {
                d,
                g,
                s,
                bulk,
                body,
                model,
                params,
            }
        }
        'J' | 'Z' => {
            let d = need!(0, "drain").to_string();
            let g = need!(1, "gate").to_string();
            let s = need!(2, "source").to_string();
            let model = need!(3, "model").to_string();
            let params: Vec<_> = rest[4..].iter().filter_map(|t| parse_kv(t)).collect();
            if letter == 'Z' {
                ElementKind::Mesa {
                    d,
                    g,
                    s,
                    model,
                    params,
                }
            } else {
                ElementKind::Jfet {
                    d,
                    g,
                    s,
                    model,
                    params,
                }
            }
        }
        'K' => {
            let l1 = need!(0, "L1").to_string();
            let l2 = need!(1, "L2").to_string();
            let coupling = parse_expr(need!(2, "coupling"));
            ElementKind::MutualCoupling { l1, l2, coupling }
        }
        'E' => {
            let out_pos = need!(0, "out+").to_string();
            let out_neg = need!(1, "out-").to_string();
            let in_pos = need!(2, "in+").to_string();
            let in_neg = need!(3, "in-").to_string();
            let gain = parse_expr(need!(4, "gain"));
            ElementKind::Vcvs {
                out_pos,
                out_neg,
                in_pos,
                in_neg,
                gain,
            }
        }
        'F' => {
            let out_pos = need!(0, "out+").to_string();
            let out_neg = need!(1, "out-").to_string();
            let vsrc = need!(2, "vsrc").to_string();
            let gain = parse_expr(need!(3, "gain"));
            ElementKind::Cccs {
                out_pos,
                out_neg,
                vsrc,
                gain,
            }
        }
        'G' => {
            let out_pos = need!(0, "out+").to_string();
            let out_neg = need!(1, "out-").to_string();
            let in_pos = need!(2, "in+").to_string();
            let in_neg = need!(3, "in-").to_string();
            let gm = parse_expr(need!(4, "gm"));
            ElementKind::Vccs {
                out_pos,
                out_neg,
                in_pos,
                in_neg,
                gm,
            }
        }
        'H' => {
            let out_pos = need!(0, "out+").to_string();
            let out_neg = need!(1, "out-").to_string();
            let vsrc = need!(2, "vsrc").to_string();
            let rm = parse_expr(need!(3, "transresistance"));
            ElementKind::Ccvs {
                out_pos,
                out_neg,
                vsrc,
                rm,
            }
        }
        'B' => {
            let pos = need!(0, "n+").to_string();
            let neg = need!(1, "n-").to_string();
            // Rest is the V= or I= expression, preserve verbatim
            let spec = rest[2..].join(" ");
            ElementKind::BehavioralSource { pos, neg, spec }
        }
        'X' => {
            // Xname port... subckt_name [PARAMS: key=val ...]
            // The subckt name is the LAST positional (non key=val) token.
            let positional: Vec<String> = rest
                .iter()
                .filter(|t| {
                    let u = t.to_uppercase();
                    u != "PARAMS:" && u != "PARAM:" && !t.contains('=')
                })
                .cloned()
                .collect();
            if positional.is_empty() {
                return Err(syntax(lineno, "X: missing subcircuit name"));
            }
            let subckt = positional.last().unwrap().clone();
            let ports = positional[..positional.len() - 1].to_vec();
            let params = collect_params(rest);
            ElementKind::SubcktCall {
                ports,
                subckt,
                params,
            }
        }
        'O' => {
            // O<name> n1+ n1- n2+ n2- model [IC=v1,i1,v2,i2]
            let pos1 = need!(0, "n1+").to_string();
            let neg1 = need!(1, "n1-").to_string();
            let pos2 = need!(2, "n2+").to_string();
            let neg2 = need!(3, "n2-").to_string();
            let model = need!(4, "model").to_string();
            let params = collect_params(&rest[5..]);
            ElementKind::Ltra {
                pos1,
                neg1,
                pos2,
                neg2,
                model,
                params,
            }
        }
        'Y' => {
            // Y<name> n1+ n1- n2+ n2- model [LEN=val]
            let pos1 = need!(0, "n1+").to_string();
            let neg1 = need!(1, "n1-").to_string();
            let pos2 = need!(2, "n2+").to_string();
            let neg2 = need!(3, "n2-").to_string();
            let model = need!(4, "model").to_string();
            let params = collect_params(&rest[5..]);
            ElementKind::Txl {
                pos1,
                neg1,
                pos2,
                neg2,
                model,
                params,
            }
        }
        'P' => {
            // P<name> <in_nodes...> gnd <out_nodes...> gnd model
            let positional: Vec<&str> = rest
                .iter()
                .filter(|t| !t.contains('='))
                .map(|s| s.as_str())
                .collect();
            // Last positional token that isn't "0" is the model name
            let model_pos = positional
                .iter()
                .rposition(|t| *t != "0")
                .ok_or_else(|| syntax(lineno, "P: missing model"))?;
            let model = positional[model_pos].to_string();
            let node_tokens = &positional[..model_pos];
            // Split node_tokens by the first "0" to get in_nodes and out_nodes
            let gnd_pos = node_tokens
                .iter()
                .position(|t| *t == "0")
                .ok_or_else(|| syntax(lineno, "P: missing ground separator"))?;
            let in_nodes: Vec<String> = node_tokens[..gnd_pos]
                .iter()
                .map(|s| s.to_string())
                .collect();
            let gnd = "0".to_string();
            let remaining = &node_tokens[gnd_pos + 1..];
            let out_end = remaining
                .iter()
                .position(|t| *t == "0")
                .unwrap_or(remaining.len());
            let out_nodes: Vec<String> =
                remaining[..out_end].iter().map(|s| s.to_string()).collect();
            let params = collect_params(rest);
            ElementKind::Cpl {
                in_nodes,
                out_nodes,
                gnd,
                model,
                params,
            }
        }
        'A' => {
            let (connections, model) = parse_xspice_connections(lineno, rest)?;
            ElementKind::Xspice { connections, model }
        }
        _ => {
            // Unknown element — store everything after the name verbatim
            ElementKind::Raw(rest.join(" "))
        }
    };

    Ok(Element { name, kind })
}

// ---------------------------------------------------------------------------
// .func parsing
// ---------------------------------------------------------------------------

/// Parse `.func name(args) body` into `Item::Func`.
fn parse_func_def(line: &str) -> Result<Item, ParseError> {
    // Find function name and argument list
    // Format: .func name(arg1[, arg2, ...]) body
    let rest = line
        .strip_prefix('.')
        .and_then(|s| {
            let s = s.trim_start();
            s.strip_prefix("func").or_else(|| s.strip_prefix("FUNC"))
        })
        .map(|s| s.trim_start())
        .unwrap_or("");

    // Find the opening paren for args
    let lparen = rest.find('(').ok_or_else(|| ParseError::Syntax {
        line: 0,
        msg: "expected '(' in .func definition".into(),
    })?;
    let name = rest[..lparen].trim().to_string();

    let rparen = rest.find(')').ok_or_else(|| ParseError::Syntax {
        line: 0,
        msg: "expected ')' in .func definition".into(),
    })?;

    let args_str = &rest[lparen + 1..rparen];
    let args: Vec<String> = if args_str.trim().is_empty() {
        Vec::new()
    } else {
        args_str.split(',').map(|a| a.trim().to_string()).collect()
    };

    // Body is everything after the closing paren
    let body_raw = rest[rparen + 1..].trim();
    // Strip surrounding quotes or braces
    let body = if let Some(inner) = body_raw.strip_prefix('{').and_then(|s| s.strip_suffix('}')) {
        inner.trim().to_string()
    } else if let Some(inner) = body_raw
        .strip_prefix('\'')
        .and_then(|s| s.strip_suffix('\''))
    {
        inner.trim().to_string()
    } else {
        body_raw.to_string()
    };

    Ok(Item::Func { name, args, body })
}

// ---------------------------------------------------------------------------
// Dot-command parsing
// ---------------------------------------------------------------------------

fn parse_dot(
    lineno: usize,
    line: &str,
    rest_lines: &mut std::iter::Peekable<impl Iterator<Item = (usize, String)>>,
) -> Result<Item, ParseError> {
    let tokens = tokenize(line);
    let keyword = tokens[0].to_uppercase();
    let kw = keyword.as_str();

    match kw {
        ".OP" => Ok(Item::Analysis(Analysis::Op)),

        ".DC" => {
            // .dc src start stop step [src2 start2 stop2 step2]
            let src = tokens
                .get(1)
                .ok_or_else(|| syntax(lineno, ".dc: missing source"))?
                .clone();
            let start = parse_expr(tokens.get(2).map(|s| s.as_str()).unwrap_or("0"));
            let stop = parse_expr(tokens.get(3).map(|s| s.as_str()).unwrap_or("0"));
            let step = parse_expr(tokens.get(4).map(|s| s.as_str()).unwrap_or("0"));
            let src2 = if tokens.len() > 8 {
                Some(DcSweep {
                    src: tokens[5].clone(),
                    start: parse_expr(&tokens[6]),
                    stop: parse_expr(&tokens[7]),
                    step: parse_expr(&tokens[8]),
                })
            } else {
                None
            };
            Ok(Item::Analysis(Analysis::Dc {
                src,
                start,
                stop,
                step,
                src2,
            }))
        }

        ".TRAN" => {
            let tstep = parse_expr(tokens.get(1).map(|s| s.as_str()).unwrap_or("0"));
            let tstop = parse_expr(tokens.get(2).map(|s| s.as_str()).unwrap_or("0"));
            let tstart = tokens.get(3).map(|s| parse_expr(s));
            let tmax = tokens.get(4).map(|s| parse_expr(s));
            Ok(Item::Analysis(Analysis::Tran {
                tstep,
                tstop,
                tstart,
                tmax,
            }))
        }

        ".AC" => {
            let variation = match tokens.get(1).map(|s| s.to_uppercase()).as_deref() {
                Some("DEC") => AcVariation::Dec,
                Some("OCT") => AcVariation::Oct,
                Some("LIN") => AcVariation::Lin,
                other => return Err(syntax(lineno, format!(".ac: unknown variation {other:?}"))),
            };
            let n: u32 = tokens.get(2).and_then(|s| s.parse().ok()).unwrap_or(1);
            let fstart = parse_expr(tokens.get(3).map(|s| s.as_str()).unwrap_or("0"));
            let fstop = parse_expr(tokens.get(4).map(|s| s.as_str()).unwrap_or("0"));
            Ok(Item::Analysis(Analysis::Ac {
                variation,
                n,
                fstart,
                fstop,
            }))
        }

        ".NOISE" => {
            // .noise output [ref] src DEC|OCT|LIN n fstart fstop
            // Detect variation keyword to split output/ref from src
            let var_idx = tokens[1..]
                .iter()
                .position(|t| matches!(t.to_uppercase().as_str(), "DEC" | "OCT" | "LIN"))
                .map(|i| i + 1); // adjust for slice offset
            let var_idx = var_idx.ok_or_else(|| syntax(lineno, ".noise: missing variation"))?;
            let before_var = &tokens[1..var_idx]; // [output, maybe ref, src]
            let (output, ref_node, src) = match before_var.len() {
                2 => (before_var[0].clone(), None, before_var[1].clone()),
                3 => (
                    before_var[0].clone(),
                    Some(before_var[1].clone()),
                    before_var[2].clone(),
                ),
                _ => return Err(syntax(lineno, ".noise: unexpected argument count")),
            };
            let variation = match tokens[var_idx].to_uppercase().as_str() {
                "DEC" => AcVariation::Dec,
                "OCT" => AcVariation::Oct,
                _ => AcVariation::Lin,
            };
            let n: u32 = tokens
                .get(var_idx + 1)
                .and_then(|s| s.parse().ok())
                .unwrap_or(1);
            let fstart = parse_expr(tokens.get(var_idx + 2).map(|s| s.as_str()).unwrap_or("0"));
            let fstop = parse_expr(tokens.get(var_idx + 3).map(|s| s.as_str()).unwrap_or("0"));
            Ok(Item::Analysis(Analysis::Noise {
                output,
                ref_node,
                src,
                variation,
                n,
                fstart,
                fstop,
            }))
        }

        ".TF" => {
            let output = tokens
                .get(1)
                .ok_or_else(|| syntax(lineno, ".tf: missing output"))?
                .clone();
            let input = tokens
                .get(2)
                .ok_or_else(|| syntax(lineno, ".tf: missing input"))?
                .clone();
            Ok(Item::Analysis(Analysis::Tf { output, input }))
        }

        ".SENS" => {
            let output = tokens[1..].to_vec();
            Ok(Item::Analysis(Analysis::Sens { output }))
        }

        ".PZ" => {
            let node_i = tokens
                .get(1)
                .ok_or_else(|| syntax(lineno, ".pz: missing nodeI"))?
                .clone();
            let node_g = tokens
                .get(2)
                .ok_or_else(|| syntax(lineno, ".pz: missing nodeG"))?
                .clone();
            let node_j = tokens
                .get(3)
                .ok_or_else(|| syntax(lineno, ".pz: missing nodeJ"))?
                .clone();
            let node_k = tokens
                .get(4)
                .ok_or_else(|| syntax(lineno, ".pz: missing nodeK"))?
                .clone();
            let input_type = match tokens.get(5).map(|s| s.to_uppercase()).as_deref() {
                Some("VOL") => PzInputType::Vol,
                Some("CUR") => PzInputType::Cur,
                _ => {
                    return Err(syntax(
                        lineno,
                        ".pz: invalid input type, expected VOL or CUR",
                    ));
                }
            };
            let analysis_type = match tokens.get(6).map(|s| s.to_uppercase()).as_deref() {
                Some("POL") => PzAnalysisType::Pol,
                Some("ZER") => PzAnalysisType::Zer,
                Some("PZ") => PzAnalysisType::Pz,
                _ => {
                    return Err(syntax(
                        lineno,
                        ".pz: invalid analysis type, expected POL, ZER, or PZ",
                    ));
                }
            };
            Ok(Item::Analysis(Analysis::Pz {
                node_i,
                node_g,
                node_j,
                node_k,
                input_type,
                analysis_type,
            }))
        }

        ".MODEL" => {
            let name = tokens
                .get(1)
                .ok_or_else(|| syntax(lineno, ".model: missing name"))?
                .clone();
            let raw_kind = tokens
                .get(2)
                .ok_or_else(|| syntax(lineno, ".model: missing type"))?
                .clone();

            // Handle parenthesized params: "NPN(BF=80 RB=100 ...)"
            // The tokenizer keeps everything inside parens as one token.
            let (kind, params) = if let Some(paren_start) = raw_kind.find('(') {
                let type_part = raw_kind[..paren_start].to_uppercase();
                let inner = raw_kind[paren_start + 1..].trim_end_matches(')');
                // Re-tokenize the inner params (they are space-separated key=value pairs)
                let inner_tokens: Vec<String> =
                    inner.split_whitespace().map(String::from).collect();
                let mut params = collect_params(&inner_tokens);
                // Also collect any params after the parenthesized token
                params.extend(collect_params(&tokens[3..]));
                (type_part, params)
            } else {
                (raw_kind.to_uppercase(), collect_model_params(&tokens[3..]))
            };
            Ok(Item::Model(ModelDef { name, kind, params }))
        }

        ".SUBCKT" => {
            let name = tokens
                .get(1)
                .ok_or_else(|| syntax(lineno, ".subckt: missing name"))?
                .clone();
            // Ports are positional tokens before the first key=val or PARAMS:
            let mut ports = Vec::new();
            let mut params_start = tokens.len();
            for (i, tok) in tokens[2..].iter().enumerate() {
                let u = tok.to_uppercase();
                if u == "PARAMS:" || u == "PARAM:" || tok.contains('=') {
                    params_start = i + 2; // +2 because we're indexing tokens[2..]
                    break;
                }
                ports.push(tok.clone());
            }
            let params = collect_params(&tokens[params_start..]);
            let items = parse_subckt_body(rest_lines)?;
            Ok(Item::Subckt(SubcktDef {
                name,
                ports,
                params,
                items,
            }))
        }

        ".ENDS" => {
            // Should be consumed by parse_subckt_body; hitting it at top level
            // is a format quirk — drop it silently.
            Ok(Item::Raw(line.to_string()))
        }

        ".PARAM" | ".PARAMETERS" => {
            let params = collect_params(&tokens[1..]);
            Ok(Item::Param(params))
        }

        ".FUNC" => {
            // .func name(arg1, arg2, ...) body
            // or .func name(arg1, arg2, ...) {body}
            // or .func name(arg1, arg2, ...) 'body'
            parse_func_def(line)
        }

        ".INCLUDE" => {
            let file = tokens
                .get(1)
                .map(|s| s.trim_matches('"').to_string())
                .unwrap_or_default();
            Ok(Item::Include(file))
        }

        ".LIB" => {
            let file = tokens
                .get(1)
                .map(|s| s.trim_matches(|c| c == '"' || c == '\'').to_string())
                .unwrap_or_default();
            let entry = tokens.get(2).cloned();
            Ok(Item::Lib { file, entry })
        }

        ".GLOBAL" => {
            let nodes = tokens[1..].to_vec();
            Ok(Item::Global(nodes))
        }

        ".OPTIONS" | ".OPTION" => {
            let opts = collect_params(&tokens[1..]);
            Ok(Item::Options(opts))
        }

        ".SAVE" | ".PROBE" => {
            let vecs = tokens[1..].to_vec();
            Ok(Item::Save(vecs))
        }

        ".CONTROL" => {
            // Skip everything until .endc
            while let Some((_, l)) = rest_lines.peek() {
                if l.to_uppercase().trim_start().starts_with(".ENDC") {
                    rest_lines.next();
                    break;
                }
                rest_lines.next();
            }
            Ok(Item::Raw(String::new()))
        }

        ".ENDL" => {
            // Library section end marker — handled by .lib processing
            Ok(Item::Raw(line.to_string()))
        }

        ".END" => {
            // Sentinel — signals end of netlist; the caller should stop.
            Ok(Item::Raw(".end".to_string()))
        }

        ".TEMP" => {
            // Temperature specification: .temp <value_in_celsius>
            let val = tokens
                .get(1)
                .and_then(|t| parse_spice_number(t))
                .unwrap_or(27.0);
            Ok(Item::Temp(val))
        }

        _ => Ok(Item::Raw(line.to_string())),
    }
}

// ---------------------------------------------------------------------------
// Subcircuit body parsing (recursive)
// ---------------------------------------------------------------------------

fn parse_subckt_body(
    lines: &mut std::iter::Peekable<impl Iterator<Item = (usize, String)>>,
) -> Result<Vec<Item>, ParseError> {
    let mut items = Vec::new();
    loop {
        match lines.peek() {
            None => break,
            Some((_, l)) if l.to_uppercase().trim_start().starts_with(".ENDS") => {
                lines.next(); // consume .ends
                break;
            }
            _ => {}
        }
        let (lineno, line) = lines.next().unwrap();
        items.push(parse_line(lineno, &line, lines)?);
    }
    Ok(items)
}

// ---------------------------------------------------------------------------
// Single logical line dispatcher
// ---------------------------------------------------------------------------

fn parse_line(
    lineno: usize,
    line: &str,
    rest: &mut std::iter::Peekable<impl Iterator<Item = (usize, String)>>,
) -> Result<Item, ParseError> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Ok(Item::Raw(String::new()));
    }

    let first = trimmed.chars().next().unwrap();

    if first == '*' {
        return Ok(Item::Comment(trimmed[1..].trim().to_string()));
    }

    if first == '.' {
        return parse_dot(lineno, trimmed, rest);
    }

    // Circuit element
    parse_element(lineno, trimmed).map(Item::Element)
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

pub fn parse(input: &str) -> Result<Netlist, ParseError> {
    let logical = preprocess(input);
    let logical = process_conditionals(logical);
    let mut iter = logical.into_iter().peekable();

    // First logical line is always the title
    let (_, title) = iter.next().ok_or(ParseError::Empty)?;

    let mut items = Vec::new();
    while let Some((lineno, line)) = iter.next() {
        // .end (without 's') terminates the top-level netlist
        if line.to_uppercase().trim_start() == ".END" {
            break;
        }
        items.push(parse_line(lineno, &line, &mut iter)?);
    }

    Ok(Netlist { title, items, source: input.to_string() })
}

// ---------------------------------------------------------------------------
// Parser tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ElementKind, Item};
    use approx::assert_abs_diff_eq;
    #[cfg(target_arch = "wasm32")]
    use wasm_bindgen_test::wasm_bindgen_test as test;

    fn parse_ok(src: &str) -> Netlist {
        Netlist::parse(src).expect("parse failed")
    }

    // ---- pre-processing ------------------------------------------------

    #[test]
    fn continuation_lines_folded() {
        let n = parse_ok("Title\n.model NMOD NMOS\n+ VTO=1.5 KP=2e-5\n.end");
        assert_eq!(n.items.len(), 1);
        if let Item::Model(m) = &n.items[0] {
            assert_eq!(m.params.len(), 2);
        } else {
            panic!("expected Model, got {:?}", n.items[0]);
        }
    }

    #[test]
    fn inline_comment_stripped() {
        let n = parse_ok("Title\nR1 a b 1k $ this is a comment\n.end");
        let e = n.elements().next().unwrap();
        if let ElementKind::Resistor { value, .. } = &e.kind {
            assert_abs_diff_eq!(
                match value {
                    crate::Expr::Num(v) => *v,
                    _ => panic!(),
                },
                1e3
            );
        }
    }

    #[test]
    fn inline_comment_inside_parens_not_stripped() {
        // The $ inside PULSE(...) must NOT be treated as a comment delimiter
        let n = parse_ok("Title\nV1 a 0 PULSE(0 1 1n 1n 1n 5n 10n)\n.end");
        let e = n.elements().next().unwrap();
        assert!(matches!(
            &e.kind,
            ElementKind::VoltageSource {
                source: crate::Source {
                    waveform: Some(_),
                    ..
                },
                ..
            }
        ));
    }

    #[test]
    fn title_is_first_line_even_if_star_comment_follows() {
        let n = parse_ok("My Circuit\n* this is a comment\nR1 a b 1k\n.end");
        assert_eq!(n.title, "My Circuit");
        assert_eq!(n.items.len(), 2); // comment + element
    }

    // ---- values --------------------------------------------------------

    #[test]
    fn m_is_milli_not_mega() {
        assert_abs_diff_eq!(parse_spice_number("1m").unwrap(), 1e-3, epsilon = 1e-15);
        assert_abs_diff_eq!(parse_spice_number("1M").unwrap(), 1e-3, epsilon = 1e-15);
        assert_abs_diff_eq!(parse_spice_number("1Meg").unwrap(), 1e6, epsilon = 1e0);
        assert_abs_diff_eq!(parse_spice_number("1MEG").unwrap(), 1e6, epsilon = 1e0);
    }

    #[test]
    fn mil_suffix() {
        assert_abs_diff_eq!(
            parse_spice_number("1mil").unwrap(),
            25.4e-6,
            epsilon = 1e-15
        );
    }

    #[test]
    fn scientific_notation_no_suffix() {
        assert_abs_diff_eq!(
            parse_spice_number("1.5e-9").unwrap(),
            1.5e-9,
            epsilon = 1e-18
        );
        assert_abs_diff_eq!(parse_spice_number("2E3").unwrap(), 2000.0, epsilon = 1e-9);
    }

    #[test]
    fn value_with_unit_suffix_ignored() {
        // "1kohm" → suffix K → 1e3, trailing "ohm" ignored
        assert_abs_diff_eq!(parse_spice_number("1kohm").unwrap(), 1e3, epsilon = 1e-9);
    }

    // ---- elements -------------------------------------------------------

    #[test]
    fn resistor_capacitor_inductor() {
        let n = parse_ok("T\nR1 a b 1k\nC2 c d 100n\nL3 e f 1u\n.end");
        let elems: Vec<_> = n.elements().collect();
        assert_eq!(elems.len(), 3);
        assert!(matches!(&elems[0].kind, ElementKind::Resistor { .. }));
        assert!(matches!(&elems[1].kind, ElementKind::Capacitor { .. }));
        assert!(matches!(&elems[2].kind, ElementKind::Inductor { .. }));
    }

    #[test]
    fn voltage_source_dc() {
        let n = parse_ok("T\nV1 vcc 0 DC 5\n.end");
        if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
            assert!(matches!(source.dc, Some(Expr::Num(_))));
        }
    }

    #[test]
    fn voltage_source_pulse() {
        let n = parse_ok("T\nV1 a 0 PULSE(0 5 0 1n 1n 5n 10n)\n.end");
        if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
            assert!(matches!(source.waveform, Some(Waveform::Pulse { .. })));
        }
    }

    #[test]
    fn voltage_source_dc_and_ac() {
        let n = parse_ok("T\nV1 a 0 DC 0 AC 1 0\n.end");
        if let ElementKind::VoltageSource { source, .. } = &n.elements().next().unwrap().kind {
            assert!(source.dc.is_some());
            assert!(source.ac.is_some());
        }
    }

    #[test]
    fn mosfet_with_params() {
        let n = parse_ok("T\nM1 d g s b NMOD W=10u L=1u\n.end");
        if let ElementKind::Mosfet { model, params, .. } = &n.elements().next().unwrap().kind {
            assert_eq!(model, "NMOD");
            assert_eq!(params.len(), 2);
            assert_eq!(params[0].name, "W");
            assert_eq!(params[1].name, "L");
        } else {
            panic!();
        }
    }

    #[test]
    fn subcircuit_call_x() {
        let n = parse_ok("T\nX1 in out gnd OPAMP gain=100\n.end");
        if let ElementKind::SubcktCall {
            ports,
            subckt,
            params,
        } = &n.elements().next().unwrap().kind
        {
            assert_eq!(subckt, "OPAMP");
            assert_eq!(ports, &["in", "out", "gnd"]);
            assert_eq!(params.len(), 1);
            assert_eq!(params[0].name, "gain");
        } else {
            panic!();
        }
    }

    #[test]
    fn bjt_with_substrate() {
        let n = parse_ok("T\nQ1 c b e sub BC547\n.end");
        if let ElementKind::Bjt {
            substrate, model, ..
        } = &n.elements().next().unwrap().kind
        {
            assert_eq!(substrate.as_deref(), Some("sub"));
            assert_eq!(model, "BC547");
        } else {
            panic!();
        }
    }

    // ---- dot commands ---------------------------------------------------

    #[test]
    fn dot_op() {
        let n = parse_ok("T\n.op\n.end");
        assert!(matches!(n.items[0], Item::Analysis(Analysis::Op)));
    }

    #[test]
    fn dot_tran() {
        let n = parse_ok("T\n.tran 0.1m 10m\n.end");
        if let Item::Analysis(Analysis::Tran {
            tstep,
            tstop,
            tstart,
            tmax,
        }) = &n.items[0]
        {
            assert_abs_diff_eq!(
                match tstep {
                    Expr::Num(v) => *v,
                    _ => panic!(),
                },
                1e-4
            );
            assert_abs_diff_eq!(
                match tstop {
                    Expr::Num(v) => *v,
                    _ => panic!(),
                },
                1e-2
            );
            assert!(tstart.is_none());
            assert!(tmax.is_none());
        } else {
            panic!("{:?}", n.items[0]);
        }
    }

    #[test]
    fn dot_dc() {
        let n = parse_ok("T\n.dc V1 0 5 1\n.end");
        if let Item::Analysis(Analysis::Dc {
            src,
            start: _,
            stop,
            step: _,
            src2,
        }) = &n.items[0]
        {
            assert_eq!(src, "V1");
            assert_abs_diff_eq!(
                match stop {
                    Expr::Num(v) => *v,
                    _ => panic!(),
                },
                5.0
            );
            assert!(src2.is_none());
        } else {
            panic!();
        }
    }

    #[test]
    fn dot_ac() {
        let n = parse_ok("T\n.ac DEC 10 1 1Meg\n.end");
        if let Item::Analysis(Analysis::Ac {
            variation,
            n: pts,
            fstart: _,
            fstop,
        }) = &n.items[0]
        {
            assert_eq!(*variation, AcVariation::Dec);
            assert_eq!(*pts, 10);
            assert_abs_diff_eq!(
                match fstop {
                    Expr::Num(v) => *v,
                    _ => panic!(),
                },
                1e6
            );
        } else {
            panic!();
        }
    }

    #[test]
    fn dot_model() {
        let n = parse_ok("T\n.model NMOD NMOS VTO=1.5 KP=2e-5\n.end");
        if let Item::Model(m) = &n.items[0] {
            assert_eq!(m.name, "NMOD");
            assert_eq!(m.kind, "NMOS");
            assert_eq!(m.params.len(), 2);
        } else {
            panic!();
        }
    }

    #[test]
    fn dot_subckt_roundtrip() {
        let src = "Buffer\n.subckt BUF in out\nR1 in out 0\n.ends BUF\n.end";
        let n = parse_ok(src);
        if let Item::Subckt(s) = &n.items[0] {
            assert_eq!(s.name, "BUF");
            assert_eq!(s.ports, ["in", "out"]);
            assert_eq!(s.items.len(), 1);
        } else {
            panic!();
        }
    }

    #[test]
    fn dot_param() {
        let n = parse_ok("T\n.param Rval=1k Cval=100n\n.end");
        if let Item::Param(ps) = &n.items[0] {
            assert_eq!(ps.len(), 2);
            assert_eq!(ps[0].name, "Rval");
        } else {
            panic!();
        }
    }

    // ---- round-trip ---------------------------------------------------

    #[test]
    fn roundtrip_simple_divider() {
        let src = "\
Voltage divider
V1 in 0 DC 5
R1 in mid 1k
R2 mid 0 1k
.op
.end";
        let n = parse_ok(src);
        let generated = n.to_string();
        // Re-parse and check key properties survive
        let n2 = Netlist::parse(&generated).unwrap();
        assert_eq!(n2.title, "Voltage divider");
        assert_eq!(n2.elements().count(), 3);
    }

    #[test]
    fn roundtrip_tran_pulse_source() {
        let src = "\
RC step
V1 in 0 PULSE(0 5 0 1n 1n 5n 10n)
R1 in out 1k
C1 out 0 100n
.tran 0.1n 50n
.end";
        let n = parse_ok(src);
        let out = n.to_string();
        let n2 = Netlist::parse(&out).unwrap();
        if let ElementKind::VoltageSource { source, .. } = &n2.elements().next().unwrap().kind {
            assert!(matches!(source.waveform, Some(Waveform::Pulse { .. })));
        }
    }
}
