//! Expression evaluator for SPICE parameter expressions.
//!
//! Supports arithmetic, comparison, boolean, ternary operators and
//! standard math functions matching ngspice's expression syntax.

use std::collections::BTreeMap;

use thevenin_types::Expr;
use thiserror::Error;

use crate::mna::MnaError;

#[derive(Error, Debug, Clone)]
pub enum ExprError {
    #[error("unknown variable: {0}")]
    UnknownVariable(String),
    #[error("unknown function: {0}")]
    UnknownFunction(String),
    #[error("wrong number of arguments for {name}: expected {expected}, got {got}")]
    WrongArgCount {
        name: String,
        expected: usize,
        got: usize,
    },
    #[error("parse error in expression: {0}")]
    ParseError(String),
}

/// Context for expression evaluation: parameters and user-defined functions.
#[derive(Debug, Clone, Default)]
pub struct EvalContext {
    /// Parameter values (case-insensitive, stored uppercase).
    pub params: BTreeMap<String, f64>,
    /// User-defined functions: name -> (arg_names, body_expression).
    pub funcs: BTreeMap<String, (Vec<String>, String)>,
}

impl EvalContext {
    /// Evaluate a `thevenin_types::Expr` to a numeric value.
    pub fn eval_expr(&self, expr: &Expr) -> Result<f64, ExprError> {
        match expr {
            Expr::Num(v) => Ok(*v),
            Expr::Param(name) => {
                let key = name.to_uppercase();
                self.params
                    .get(&key)
                    .copied()
                    .ok_or_else(|| ExprError::UnknownVariable(name.clone()))
            }
            Expr::Brace(s) => self.eval_str(s),
        }
    }

    /// Evaluate an expression string to a numeric value.
    pub fn eval_str(&self, s: &str) -> Result<f64, ExprError> {
        let tokens = tokenize(s)?;
        let mut pos = 0;
        let result = parse_ternary(&tokens, &mut pos, self)?;
        if pos < tokens.len() {
            return Err(ExprError::ParseError(format!(
                "unexpected token at position {pos}: {s}"
            )));
        }
        Ok(result)
    }

    /// Evaluate a `thevenin_types::Expr`, returning MnaError on failure.
    pub fn eval_expr_mna(&self, expr: &Expr, context: &str) -> Result<f64, MnaError> {
        self.eval_expr(expr).map_err(|e| MnaError::ExprError {
            element: context.to_string(),
            detail: e.to_string(),
        })
    }
}

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Num(f64),
    Ident(String),
    Plus,
    Minus,
    Star,
    Slash,
    Percent,
    Power,
    Lparen,
    Rparen,
    Comma,
    Question,
    Colon,
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    And,
    Or,
    Not,
}

fn tokenize(input: &str) -> Result<Vec<Token>, ExprError> {
    let bytes = input.as_bytes();
    let mut tokens = Vec::new();
    let mut i = 0;

    while i < bytes.len() {
        let b = bytes[i];

        // Skip whitespace
        if b.is_ascii_whitespace() {
            i += 1;
            continue;
        }

        // Number: digit or decimal point followed by digit
        if b.is_ascii_digit() || (b == b'.' && i + 1 < bytes.len() && bytes[i + 1].is_ascii_digit())
        {
            let start = i;
            // Consume digits and dots
            while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
                i += 1;
            }
            // Scientific notation
            if i < bytes.len() && (bytes[i] == b'e' || bytes[i] == b'E') {
                i += 1;
                if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
                    i += 1;
                }
                while i < bytes.len() && bytes[i].is_ascii_digit() {
                    i += 1;
                }
            }
            // SPICE SI suffix
            let num_end = i;
            let suffix_start = i;
            if i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                // Collect suffix
                let suf_start = i;
                while i < bytes.len() && bytes[i].is_ascii_alphabetic() {
                    i += 1;
                }
                let suffix = &input[suf_start..i];
                let multiplier = spice_suffix(suffix);
                if let Some(mult) = multiplier {
                    let base: f64 = input[start..num_end].parse().map_err(|_| {
                        ExprError::ParseError(format!("bad number: {}", &input[start..i]))
                    })?;
                    tokens.push(Token::Num(base * mult));
                    continue;
                }
                // Not a known suffix — rewind
                i = suffix_start;
            }
            let s = &input[start..num_end];
            let v: f64 = s
                .parse()
                .map_err(|_| ExprError::ParseError(format!("bad number: {s}")))?;
            tokens.push(Token::Num(v));
            continue;
        }

        // Identifier (starts with letter or underscore)
        if b.is_ascii_alphabetic() || b == b'_' {
            let start = i;
            while i < bytes.len()
                && (bytes[i].is_ascii_alphanumeric() || bytes[i] == b'_' || bytes[i] == b'.')
            {
                i += 1;
            }
            let name = &input[start..i];
            // Check for named constants
            match name.to_uppercase().as_str() {
                "PI" => tokens.push(Token::Num(std::f64::consts::PI)),
                "E" if !matches!(tokens.last(), Some(Token::Num(_))) => {
                    tokens.push(Token::Num(std::f64::consts::E));
                }
                "TRUE" | "YES" => tokens.push(Token::Num(1.0)),
                "FALSE" | "NO" => tokens.push(Token::Num(0.0)),
                _ => tokens.push(Token::Ident(name.to_string())),
            }
            continue;
        }

        // Two-character operators
        if i + 1 < bytes.len() {
            let two = &input[i..i + 2];
            match two {
                "**" => {
                    tokens.push(Token::Power);
                    i += 2;
                    continue;
                }
                "==" => {
                    tokens.push(Token::Eq);
                    i += 2;
                    continue;
                }
                "!=" => {
                    tokens.push(Token::Ne);
                    i += 2;
                    continue;
                }
                "<=" => {
                    tokens.push(Token::Le);
                    i += 2;
                    continue;
                }
                ">=" => {
                    tokens.push(Token::Ge);
                    i += 2;
                    continue;
                }
                "&&" => {
                    tokens.push(Token::And);
                    i += 2;
                    continue;
                }
                "||" => {
                    tokens.push(Token::Or);
                    i += 2;
                    continue;
                }
                _ => {}
            }
        }

        // Single-character operators
        match b {
            b'+' => tokens.push(Token::Plus),
            b'-' => tokens.push(Token::Minus),
            b'*' => tokens.push(Token::Star),
            b'/' => tokens.push(Token::Slash),
            b'%' => tokens.push(Token::Percent),
            b'^' => tokens.push(Token::Power),
            b'(' => tokens.push(Token::Lparen),
            b')' => tokens.push(Token::Rparen),
            b',' => tokens.push(Token::Comma),
            b'?' => tokens.push(Token::Question),
            b':' => tokens.push(Token::Colon),
            b'<' => tokens.push(Token::Lt),
            b'>' => tokens.push(Token::Gt),
            b'!' => tokens.push(Token::Not),
            b'~' => tokens.push(Token::Not), // bitwise not treated as logical not
            _ => {
                return Err(ExprError::ParseError(format!(
                    "unexpected character: '{}'",
                    b as char
                )));
            }
        }
        i += 1;
    }

    Ok(tokens)
}

fn spice_suffix(s: &str) -> Option<f64> {
    let su = s.to_uppercase();
    match su.as_str() {
        "T" => Some(1e12),
        "G" => Some(1e9),
        "MEG" => Some(1e6),
        "K" => Some(1e3),
        "M" | "MIL" => {
            // M alone is milli in SPICE, MIL is 25.4e-6
            if su == "MIL" {
                Some(25.4e-6)
            } else {
                Some(1e-3)
            }
        }
        "U" => Some(1e-6),
        "N" => Some(1e-9),
        "P" => Some(1e-12),
        "F" => Some(1e-15),
        "A" => Some(1e-18),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Recursive descent parser/evaluator
// ---------------------------------------------------------------------------

fn peek(tokens: &[Token], pos: usize) -> Option<&Token> {
    tokens.get(pos)
}

fn is_unary_context(tokens: &[Token], pos: usize) -> bool {
    if pos == 0 {
        return true;
    }
    matches!(
        tokens[pos - 1],
        Token::Plus
            | Token::Minus
            | Token::Star
            | Token::Slash
            | Token::Percent
            | Token::Power
            | Token::Lparen
            | Token::Comma
            | Token::Question
            | Token::Colon
            | Token::Eq
            | Token::Ne
            | Token::Lt
            | Token::Le
            | Token::Gt
            | Token::Ge
            | Token::And
            | Token::Or
            | Token::Not
    )
}

// Precedence 1 (lowest): ternary ? :
fn parse_ternary(tokens: &[Token], pos: &mut usize, ctx: &EvalContext) -> Result<f64, ExprError> {
    let cond = parse_or(tokens, pos, ctx)?;
    if matches!(peek(tokens, *pos), Some(Token::Question)) {
        *pos += 1;
        let then_val = parse_ternary(tokens, pos, ctx)?;
        if !matches!(peek(tokens, *pos), Some(Token::Colon)) {
            return Err(ExprError::ParseError(
                "expected ':' in ternary expression".into(),
            ));
        }
        *pos += 1;
        let else_val = parse_ternary(tokens, pos, ctx)?;
        Ok(if cond != 0.0 { then_val } else { else_val })
    } else {
        Ok(cond)
    }
}

// Precedence 2: ||
fn parse_or(tokens: &[Token], pos: &mut usize, ctx: &EvalContext) -> Result<f64, ExprError> {
    let mut left = parse_and(tokens, pos, ctx)?;
    while matches!(peek(tokens, *pos), Some(Token::Or)) {
        *pos += 1;
        let right = parse_and(tokens, pos, ctx)?;
        left = if left != 0.0 || right != 0.0 {
            1.0
        } else {
            0.0
        };
    }
    Ok(left)
}

// Precedence 3: &&
fn parse_and(tokens: &[Token], pos: &mut usize, ctx: &EvalContext) -> Result<f64, ExprError> {
    let mut left = parse_equality(tokens, pos, ctx)?;
    while matches!(peek(tokens, *pos), Some(Token::And)) {
        *pos += 1;
        let right = parse_equality(tokens, pos, ctx)?;
        left = if left != 0.0 && right != 0.0 {
            1.0
        } else {
            0.0
        };
    }
    Ok(left)
}

// Precedence 4: == !=
fn parse_equality(tokens: &[Token], pos: &mut usize, ctx: &EvalContext) -> Result<f64, ExprError> {
    let mut left = parse_comparison(tokens, pos, ctx)?;
    loop {
        match peek(tokens, *pos) {
            Some(Token::Eq) => {
                *pos += 1;
                let right = parse_comparison(tokens, pos, ctx)?;
                left = if (left - right).abs() < f64::EPSILON {
                    1.0
                } else {
                    0.0
                };
            }
            Some(Token::Ne) => {
                *pos += 1;
                let right = parse_comparison(tokens, pos, ctx)?;
                left = if (left - right).abs() >= f64::EPSILON {
                    1.0
                } else {
                    0.0
                };
            }
            _ => break,
        }
    }
    Ok(left)
}

// Precedence 5: < <= > >=
fn parse_comparison(
    tokens: &[Token],
    pos: &mut usize,
    ctx: &EvalContext,
) -> Result<f64, ExprError> {
    let mut left = parse_additive(tokens, pos, ctx)?;
    loop {
        match peek(tokens, *pos) {
            Some(Token::Lt) => {
                *pos += 1;
                let right = parse_additive(tokens, pos, ctx)?;
                left = if left < right { 1.0 } else { 0.0 };
            }
            Some(Token::Le) => {
                *pos += 1;
                let right = parse_additive(tokens, pos, ctx)?;
                left = if left <= right { 1.0 } else { 0.0 };
            }
            Some(Token::Gt) => {
                *pos += 1;
                let right = parse_additive(tokens, pos, ctx)?;
                left = if left > right { 1.0 } else { 0.0 };
            }
            Some(Token::Ge) => {
                *pos += 1;
                let right = parse_additive(tokens, pos, ctx)?;
                left = if left >= right { 1.0 } else { 0.0 };
            }
            _ => break,
        }
    }
    Ok(left)
}

// Precedence 6: + -
fn parse_additive(tokens: &[Token], pos: &mut usize, ctx: &EvalContext) -> Result<f64, ExprError> {
    let mut left = parse_multiplicative(tokens, pos, ctx)?;
    loop {
        match peek(tokens, *pos) {
            Some(Token::Plus) if !is_unary_context(tokens, *pos) => {
                *pos += 1;
                let right = parse_multiplicative(tokens, pos, ctx)?;
                left += right;
            }
            Some(Token::Minus) if !is_unary_context(tokens, *pos) => {
                *pos += 1;
                let right = parse_multiplicative(tokens, pos, ctx)?;
                left -= right;
            }
            _ => break,
        }
    }
    Ok(left)
}

// Precedence 7: * / %
fn parse_multiplicative(
    tokens: &[Token],
    pos: &mut usize,
    ctx: &EvalContext,
) -> Result<f64, ExprError> {
    let mut left = parse_power(tokens, pos, ctx)?;
    loop {
        match peek(tokens, *pos) {
            Some(Token::Star) => {
                *pos += 1;
                let right = parse_power(tokens, pos, ctx)?;
                left *= right;
            }
            Some(Token::Slash) => {
                *pos += 1;
                let right = parse_power(tokens, pos, ctx)?;
                left /= right;
            }
            Some(Token::Percent) => {
                *pos += 1;
                let right = parse_power(tokens, pos, ctx)?;
                left %= right;
            }
            _ => break,
        }
    }
    Ok(left)
}

// Precedence 8: ** ^ (right-associative)
fn parse_power(tokens: &[Token], pos: &mut usize, ctx: &EvalContext) -> Result<f64, ExprError> {
    let base = parse_unary(tokens, pos, ctx)?;
    if matches!(peek(tokens, *pos), Some(Token::Power)) {
        *pos += 1;
        let exp = parse_power(tokens, pos, ctx)?; // right-associative
        Ok(base.powf(exp))
    } else {
        Ok(base)
    }
}

// Precedence 9: unary + - !
fn parse_unary(tokens: &[Token], pos: &mut usize, ctx: &EvalContext) -> Result<f64, ExprError> {
    match peek(tokens, *pos) {
        Some(Token::Plus) if is_unary_context(tokens, *pos) => {
            *pos += 1;
            parse_unary(tokens, pos, ctx)
        }
        Some(Token::Minus) if is_unary_context(tokens, *pos) => {
            *pos += 1;
            let val = parse_unary(tokens, pos, ctx)?;
            Ok(-val)
        }
        Some(Token::Not) => {
            *pos += 1;
            let val = parse_unary(tokens, pos, ctx)?;
            Ok(if val == 0.0 { 1.0 } else { 0.0 })
        }
        _ => parse_primary(tokens, pos, ctx),
    }
}

// Precedence 10: atoms, function calls, parenthesized expressions
fn parse_primary(tokens: &[Token], pos: &mut usize, ctx: &EvalContext) -> Result<f64, ExprError> {
    match peek(tokens, *pos) {
        Some(Token::Num(v)) => {
            let v = *v;
            *pos += 1;
            Ok(v)
        }
        Some(Token::Lparen) => {
            *pos += 1;
            let val = parse_ternary(tokens, pos, ctx)?;
            if !matches!(peek(tokens, *pos), Some(Token::Rparen)) {
                return Err(ExprError::ParseError("expected ')'".into()));
            }
            *pos += 1;
            Ok(val)
        }
        Some(Token::Ident(name)) => {
            let name = name.clone();
            *pos += 1;
            // Function call?
            if matches!(peek(tokens, *pos), Some(Token::Lparen)) {
                *pos += 1;
                let mut args = Vec::new();
                if !matches!(peek(tokens, *pos), Some(Token::Rparen)) {
                    args.push(parse_ternary(tokens, pos, ctx)?);
                    while matches!(peek(tokens, *pos), Some(Token::Comma)) {
                        *pos += 1;
                        args.push(parse_ternary(tokens, pos, ctx)?);
                    }
                }
                if !matches!(peek(tokens, *pos), Some(Token::Rparen)) {
                    return Err(ExprError::ParseError(
                        "expected ')' after function args".into(),
                    ));
                }
                *pos += 1;
                eval_function(&name, &args, ctx)
            } else {
                // Variable lookup
                let key = name.to_uppercase();
                ctx.params
                    .get(&key)
                    .copied()
                    .ok_or(ExprError::UnknownVariable(name))
            }
        }
        Some(tok) => Err(ExprError::ParseError(format!("unexpected token: {tok:?}"))),
        None => Err(ExprError::ParseError("unexpected end of expression".into())),
    }
}

// ---------------------------------------------------------------------------
// Built-in functions
// ---------------------------------------------------------------------------

fn eval_function(name: &str, args: &[f64], ctx: &EvalContext) -> Result<f64, ExprError> {
    let name_upper = name.to_uppercase();

    // Check user-defined functions first
    if let Some((param_names, body)) = ctx.funcs.get(&name_upper) {
        if args.len() != param_names.len() {
            return Err(ExprError::WrongArgCount {
                name: name.to_string(),
                expected: param_names.len(),
                got: args.len(),
            });
        }
        // Create a child context with function parameters
        let mut child_ctx = ctx.clone();
        for (pname, &val) in param_names.iter().zip(args.iter()) {
            child_ctx.params.insert(pname.to_uppercase(), val);
        }
        return child_ctx.eval_str(body);
    }

    // Built-in functions
    match name_upper.as_str() {
        // Trigonometric
        "SIN" => require_args(name, args, 1).map(|_| args[0].sin()),
        "COS" => require_args(name, args, 1).map(|_| args[0].cos()),
        "TAN" => require_args(name, args, 1).map(|_| args[0].tan()),
        "ASIN" => require_args(name, args, 1).map(|_| args[0].asin()),
        "ACOS" => require_args(name, args, 1).map(|_| args[0].acos()),
        "ATAN" | "ARCTAN" => require_args(name, args, 1).map(|_| args[0].atan()),
        "ATAN2" => require_args(name, args, 2).map(|_| args[0].atan2(args[1])),

        // Hyperbolic
        "SINH" => require_args(name, args, 1).map(|_| args[0].sinh()),
        "COSH" => require_args(name, args, 1).map(|_| args[0].cosh()),
        "TANH" => require_args(name, args, 1).map(|_| args[0].tanh()),
        "ASINH" => require_args(name, args, 1).map(|_| args[0].asinh()),
        "ACOSH" => require_args(name, args, 1).map(|_| args[0].acosh()),
        "ATANH" => require_args(name, args, 1).map(|_| args[0].atanh()),

        // Exponential / logarithmic
        "EXP" => require_args(name, args, 1).map(|_| args[0].exp()),
        "LOG" | "LN" => require_args(name, args, 1).map(|_| args[0].ln()),
        "LOG10" => require_args(name, args, 1).map(|_| args[0].log10()),
        "SQRT" => require_args(name, args, 1).map(|_| args[0].sqrt()),
        "SQR" => require_args(name, args, 1).map(|_| args[0] * args[0]),

        // Power
        "POW" | "PWR" => require_args(name, args, 2).map(|_| args[0].abs().powf(args[1])),

        // Rounding / integer
        "ABS" => require_args(name, args, 1).map(|_| args[0].abs()),
        "SGN" | "SIGN" => require_args(name, args, 1).map(|_| {
            if args[0] > 0.0 {
                1.0
            } else if args[0] < 0.0 {
                -1.0
            } else {
                0.0
            }
        }),
        "INT" => require_args(name, args, 1).map(|_| args[0].trunc()),
        "NINT" => require_args(name, args, 1).map(|_| args[0].round()),
        "FLOOR" => require_args(name, args, 1).map(|_| args[0].floor()),
        "CEIL" | "CEILING" => require_args(name, args, 1).map(|_| args[0].ceil()),

        // Min/max
        "MIN" => require_args(name, args, 2).map(|_| args[0].min(args[1])),
        "MAX" => require_args(name, args, 2).map(|_| args[0].max(args[1])),

        // Step/ramp functions (B-source)
        "U" => require_args(name, args, 1).map(|_| if args[0] > 0.0 { 1.0 } else { 0.0 }),
        "U2" => require_args(name, args, 1).map(|_| {
            if args[0] <= 0.0 {
                0.0
            } else if args[0] < 1.0 {
                args[0]
            } else {
                1.0
            }
        }),
        "URAMP" => require_args(name, args, 1).map(|_| if args[0] > 0.0 { args[0] } else { 0.0 }),

        // Predicate functions (return 0 or 1)
        "EQ0" => require_args(name, args, 1).map(|_| if args[0] == 0.0 { 1.0 } else { 0.0 }),
        "NE0" => require_args(name, args, 1).map(|_| if args[0] != 0.0 { 1.0 } else { 0.0 }),
        "GT0" => require_args(name, args, 1).map(|_| if args[0] > 0.0 { 1.0 } else { 0.0 }),
        "LT0" => require_args(name, args, 1).map(|_| if args[0] < 0.0 { 1.0 } else { 0.0 }),
        "GE0" => require_args(name, args, 1).map(|_| if args[0] >= 0.0 { 1.0 } else { 0.0 }),
        "LE0" => require_args(name, args, 1).map(|_| if args[0] <= 0.0 { 1.0 } else { 0.0 }),

        // PWL function: pwl(x, x1, y1, x2, y2, ...)
        "PWL" => {
            if args.len() < 3 || args.len().is_multiple_of(2) {
                return Err(ExprError::WrongArgCount {
                    name: name.to_string(),
                    expected: 3, // at least x, x1, y1
                    got: args.len(),
                });
            }
            let x = args[0];
            let pairs: Vec<(f64, f64)> = args[1..].chunks(2).map(|c| (c[0], c[1])).collect();
            Ok(pwl_interp(x, &pairs))
        }

        _ => Err(ExprError::UnknownFunction(name.to_string())),
    }
}

fn require_args(name: &str, args: &[f64], expected: usize) -> Result<(), ExprError> {
    if args.len() != expected {
        Err(ExprError::WrongArgCount {
            name: name.to_string(),
            expected,
            got: args.len(),
        })
    } else {
        Ok(())
    }
}

/// Piecewise-linear interpolation.
fn pwl_interp(x: f64, pairs: &[(f64, f64)]) -> f64 {
    if pairs.is_empty() {
        return 0.0;
    }
    if x <= pairs[0].0 {
        return pairs[0].1;
    }
    if x >= pairs[pairs.len() - 1].0 {
        return pairs[pairs.len() - 1].1;
    }
    for i in 1..pairs.len() {
        if x <= pairs[i].0 {
            let (x0, y0) = pairs[i - 1];
            let (x1, y1) = pairs[i];
            let t = (x - x0) / (x1 - x0);
            return y0 + t * (y1 - y0);
        }
    }
    pairs[pairs.len() - 1].1
}

// ---------------------------------------------------------------------------
// Netlist expression resolution
// ---------------------------------------------------------------------------

/// Build an `EvalContext` from `.param` and `.func` items in a netlist.
pub fn build_context(items: &[thevenin_types::Item]) -> EvalContext {
    let mut ctx = EvalContext::default();
    collect_context(items, &mut ctx);
    ctx
}

/// Collect .param and .func definitions into a context, resolving params in order.
fn collect_context(items: &[thevenin_types::Item], ctx: &mut EvalContext) {
    for item in items {
        match item {
            thevenin_types::Item::Param(params) => {
                for p in params {
                    if let Ok(val) = ctx.eval_expr(&p.value) {
                        ctx.params.insert(p.name.to_uppercase(), val);
                    }
                }
            }
            thevenin_types::Item::Func { name, args, body } => {
                ctx.funcs.insert(
                    name.to_uppercase(),
                    (
                        args.iter().map(|a| a.to_uppercase()).collect(),
                        body.clone(),
                    ),
                );
            }
            _ => {}
        }
    }
}

/// Resolve all `Expr::Param` and `Expr::Brace` in a netlist to `Expr::Num`.
/// Also converts constant B-source expressions to V/I sources.
pub fn resolve_netlist_exprs(
    netlist: &mut thevenin_types::Netlist,
) -> Result<EvalContext, ExprError> {
    let ctx = build_context(&netlist.items);
    resolve_items(&mut netlist.items, &ctx)?;
    resolve_bsources(&mut netlist.items, &ctx)?;
    Ok(ctx)
}

fn resolve_items(items: &mut [thevenin_types::Item], ctx: &EvalContext) -> Result<(), ExprError> {
    for item in items.iter_mut() {
        match item {
            thevenin_types::Item::Element(el) => resolve_element(&mut el.kind, ctx)?,
            thevenin_types::Item::Subckt(_) => {
                // Don't resolve inside subcircuit definitions — the subcircuit
                // expander handles parameter substitution with instance params.
            }
            thevenin_types::Item::Analysis(analysis) => resolve_analysis(analysis, ctx)?,
            _ => {}
        }
    }
    Ok(())
}

fn try_resolve_expr(expr: &mut Expr, ctx: &EvalContext) {
    if let Ok(val) = ctx.eval_expr(expr) {
        *expr = Expr::Num(val);
    }
}

fn resolve_source(source: &mut thevenin_types::Source, ctx: &EvalContext) -> Result<(), ExprError> {
    if let Some(dc) = &mut source.dc {
        try_resolve_expr(dc, ctx);
    }
    if let Some(ac) = &mut source.ac {
        try_resolve_expr(&mut ac.mag, ctx);
        if let Some(phase) = &mut ac.phase {
            try_resolve_expr(phase, ctx);
        }
    }
    if let Some(wf) = &mut source.waveform {
        resolve_waveform(wf, ctx);
    }
    Ok(())
}

fn resolve_waveform(wf: &mut thevenin_types::Waveform, ctx: &EvalContext) {
    match wf {
        thevenin_types::Waveform::Pulse {
            v1,
            v2,
            td,
            tr,
            tf,
            pw,
            per,
        } => {
            try_resolve_expr(v1, ctx);
            try_resolve_expr(v2, ctx);
            for e in [td, tr, tf, pw, per].into_iter().flatten() {
                try_resolve_expr(e, ctx);
            }
        }
        thevenin_types::Waveform::Sin {
            v0,
            va,
            freq,
            td,
            theta,
            phi,
        } => {
            try_resolve_expr(v0, ctx);
            try_resolve_expr(va, ctx);
            for e in [freq, td, theta, phi].into_iter().flatten() {
                try_resolve_expr(e, ctx);
            }
        }
        thevenin_types::Waveform::Exp {
            v1,
            v2,
            td1,
            tau1,
            td2,
            tau2,
        } => {
            try_resolve_expr(v1, ctx);
            try_resolve_expr(v2, ctx);
            for e in [td1, tau1, td2, tau2].into_iter().flatten() {
                try_resolve_expr(e, ctx);
            }
        }
        thevenin_types::Waveform::Pwl(points) => {
            for p in points {
                try_resolve_expr(&mut p.time, ctx);
                try_resolve_expr(&mut p.value, ctx);
            }
        }
        thevenin_types::Waveform::Sffm { v0, va, fc, fs, md } => {
            try_resolve_expr(v0, ctx);
            try_resolve_expr(va, ctx);
            for e in [fc, fs, md].into_iter().flatten() {
                try_resolve_expr(e, ctx);
            }
        }
        thevenin_types::Waveform::Am { va, vo, fc, fs, td } => {
            try_resolve_expr(va, ctx);
            try_resolve_expr(vo, ctx);
            try_resolve_expr(fc, ctx);
            try_resolve_expr(fs, ctx);
            if let Some(e) = td {
                try_resolve_expr(e, ctx);
            }
        }
    }
}

fn resolve_params(params: &mut [thevenin_types::Param], ctx: &EvalContext) {
    for p in params {
        try_resolve_expr(&mut p.value, ctx);
    }
}

fn resolve_element(
    kind: &mut thevenin_types::ElementKind,
    ctx: &EvalContext,
) -> Result<(), ExprError> {
    use thevenin_types::ElementKind;
    match kind {
        ElementKind::Resistor { value, params, .. } => {
            try_resolve_expr(value, ctx);
            resolve_params(params, ctx);
        }
        ElementKind::Capacitor { value, params, .. } => {
            try_resolve_expr(value, ctx);
            resolve_params(params, ctx);
        }
        ElementKind::Inductor { value, params, .. } => {
            try_resolve_expr(value, ctx);
            resolve_params(params, ctx);
        }
        ElementKind::VoltageSource { source, .. } => {
            resolve_source(source, ctx)?;
        }
        ElementKind::CurrentSource { source, .. } => {
            resolve_source(source, ctx)?;
        }
        ElementKind::Diode { params, .. } => {
            resolve_params(params, ctx);
        }
        ElementKind::Bjt { params, .. } => {
            resolve_params(params, ctx);
        }
        ElementKind::Mosfet { params, .. } => {
            resolve_params(params, ctx);
        }
        ElementKind::Jfet { params, .. } => {
            resolve_params(params, ctx);
        }
        ElementKind::Mesa { params, .. } => {
            resolve_params(params, ctx);
        }
        ElementKind::MutualCoupling { coupling, .. } => {
            try_resolve_expr(coupling, ctx);
        }
        ElementKind::Vcvs { gain, .. } => {
            try_resolve_expr(gain, ctx);
        }
        ElementKind::Cccs { gain, .. } => {
            try_resolve_expr(gain, ctx);
        }
        ElementKind::Vccs { gm, .. } => {
            try_resolve_expr(gm, ctx);
        }
        ElementKind::Ccvs { rm, .. } => {
            try_resolve_expr(rm, ctx);
        }
        ElementKind::SubcktCall { params, .. } => {
            resolve_params(params, ctx);
        }
        ElementKind::Ltra { params, .. } => {
            resolve_params(params, ctx);
        }
        ElementKind::Txl { params, .. } => {
            resolve_params(params, ctx);
        }
        ElementKind::Cpl { params, .. } => {
            resolve_params(params, ctx);
        }
        ElementKind::BehavioralSource { .. } | ElementKind::Raw(_) => {}
    }
    Ok(())
}

/// Convert constant B-source expressions to V/I sources.
fn resolve_bsources(
    items: &mut [thevenin_types::Item],
    ctx: &EvalContext,
) -> Result<(), ExprError> {
    for item in items.iter_mut() {
        if let thevenin_types::Item::Element(el) = item
            && let thevenin_types::ElementKind::BehavioralSource { pos, neg, spec } = &el.kind
        {
            // Parse spec: "V=expr" or "I=expr" or "V = expr" etc.
            let spec_trimmed = spec.trim();
            let (is_voltage, expr_str) = if let Some(rest) = spec_trimmed
                .strip_prefix("V=")
                .or_else(|| spec_trimmed.strip_prefix("v="))
            {
                (true, rest.trim())
            } else if let Some(rest) = spec_trimmed
                .strip_prefix("V =")
                .or_else(|| spec_trimmed.strip_prefix("v ="))
            {
                (true, rest.trim())
            } else if let Some(rest) = spec_trimmed
                .strip_prefix("I=")
                .or_else(|| spec_trimmed.strip_prefix("i="))
            {
                (false, rest.trim())
            } else if let Some(rest) = spec_trimmed
                .strip_prefix("I =")
                .or_else(|| spec_trimmed.strip_prefix("i ="))
            {
                (false, rest.trim())
            } else {
                continue;
            };

            // Strip surrounding quotes/braces from expression
            let expr_clean =
                if let Some(inner) = expr_str.strip_prefix('{').and_then(|s| s.strip_suffix('}')) {
                    inner.trim()
                } else if let Some(inner) = expr_str
                    .strip_prefix('\'')
                    .and_then(|s| s.strip_suffix('\''))
                {
                    inner.trim()
                } else {
                    expr_str
                };

            // Try to evaluate as a constant expression
            if let Ok(val) = ctx.eval_str(expr_clean) {
                let pos = pos.clone();
                let neg = neg.clone();
                let source = thevenin_types::Source {
                    dc: Some(Expr::Num(val)),
                    ac: None,
                    waveform: None,
                };
                if is_voltage {
                    el.kind = thevenin_types::ElementKind::VoltageSource { pos, neg, source };
                } else {
                    el.kind = thevenin_types::ElementKind::CurrentSource { pos, neg, source };
                }
            }
            // If evaluation fails (has circuit variable references), leave as B-source
        }
    }
    Ok(())
}

fn resolve_analysis(
    analysis: &mut thevenin_types::Analysis,
    ctx: &EvalContext,
) -> Result<(), ExprError> {
    use thevenin_types::Analysis;
    match analysis {
        Analysis::Op => {}
        Analysis::Dc {
            start,
            stop,
            step,
            src2,
            ..
        } => {
            try_resolve_expr(start, ctx);
            try_resolve_expr(stop, ctx);
            try_resolve_expr(step, ctx);
            if let Some(s2) = src2 {
                try_resolve_expr(&mut s2.start, ctx);
                try_resolve_expr(&mut s2.stop, ctx);
                try_resolve_expr(&mut s2.step, ctx);
            }
        }
        Analysis::Tran {
            tstep,
            tstop,
            tstart,
            tmax,
        } => {
            try_resolve_expr(tstep, ctx);
            try_resolve_expr(tstop, ctx);
            if let Some(e) = tstart {
                try_resolve_expr(e, ctx);
            }
            if let Some(e) = tmax {
                try_resolve_expr(e, ctx);
            }
        }
        Analysis::Ac { fstart, fstop, .. } => {
            try_resolve_expr(fstart, ctx);
            try_resolve_expr(fstop, ctx);
        }
        Analysis::Noise { fstart, fstop, .. } => {
            try_resolve_expr(fstart, ctx);
            try_resolve_expr(fstop, ctx);
        }
        Analysis::Tf { .. } | Analysis::Sens { .. } | Analysis::Pz { .. } => {}
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn eval(s: &str) -> f64 {
        let ctx = EvalContext::default();
        ctx.eval_str(s).unwrap()
    }

    fn eval_ctx(s: &str, params: &[(&str, f64)]) -> f64 {
        let mut ctx = EvalContext::default();
        for &(k, v) in params {
            ctx.params.insert(k.to_uppercase(), v);
        }
        ctx.eval_str(s).unwrap()
    }

    #[test]
    fn basic_arithmetic() {
        assert_eq!(eval("1+2"), 3.0);
        assert_eq!(eval("1+2*3"), 7.0);
        assert_eq!(eval("(1+2)*3"), 9.0);
        assert_eq!(eval("10/4"), 2.5);
        assert_eq!(eval("10%3"), 1.0);
    }

    #[test]
    fn double_minus() {
        assert_eq!(eval("2--3"), 5.0);
        assert_eq!(eval("1--1"), 2.0);
    }

    #[test]
    fn unary() {
        assert_eq!(eval("-3"), -3.0);
        assert_eq!(eval("+3"), 3.0);
        assert_eq!(eval("-(1+2)"), -3.0);
    }

    #[test]
    fn power() {
        assert_eq!(eval("2**3"), 8.0);
        assert_eq!(eval("2^3"), 8.0);
        // Right-associative: 2^3^2 = 2^(3^2) = 2^9 = 512
        assert_eq!(eval("2^3^2"), 512.0);
    }

    #[test]
    fn comparison() {
        assert_eq!(eval("1<2"), 1.0);
        assert_eq!(eval("2<1"), 0.0);
        assert_eq!(eval("1<=1"), 1.0);
        assert_eq!(eval("1==1"), 1.0);
        assert_eq!(eval("1!=2"), 1.0);
    }

    #[test]
    fn boolean_ops() {
        assert_eq!(eval("1&&1"), 1.0);
        assert_eq!(eval("1&&0"), 0.0);
        assert_eq!(eval("0||1"), 1.0);
        assert_eq!(eval("!0"), 1.0);
        assert_eq!(eval("!1"), 0.0);
    }

    #[test]
    fn ternary() {
        assert_eq!(eval("1 ? 10 : 20"), 10.0);
        assert_eq!(eval("0 ? 10 : 20"), 20.0);
    }

    #[test]
    fn math_functions() {
        let eps = 1e-12;
        assert!((eval("sin(0)") - 0.0).abs() < eps);
        assert!((eval("cos(0)") - 1.0).abs() < eps);
        assert!((eval("exp(1)") - std::f64::consts::E).abs() < eps);
        assert!((eval("log(1)") - 0.0).abs() < eps);
        assert!((eval("sqrt(4)") - 2.0).abs() < eps);
        assert!((eval("abs(-5)") - 5.0).abs() < eps);
    }

    #[test]
    fn rounding_functions() {
        assert_eq!(eval("floor(1.7)"), 1.0);
        assert_eq!(eval("ceil(1.2)"), 2.0);
        assert_eq!(eval("int(1.9)"), 1.0);
        assert_eq!(eval("nint(1.5)"), 2.0);
        assert_eq!(eval("nint(2.5)"), 3.0); // round half away from zero
    }

    #[test]
    fn parameters() {
        assert_eq!(eval_ctx("x+1", &[("x", 5.0)]), 6.0);
        assert_eq!(eval_ctx("x*y", &[("x", 3.0), ("y", 4.0)]), 12.0);
    }

    #[test]
    fn user_function() {
        let mut ctx = EvalContext::default();
        ctx.funcs.insert(
            "DOUBLE".to_string(),
            (vec!["X".to_string()], "x*2".to_string()),
        );
        assert_eq!(ctx.eval_str("double(5)").unwrap(), 10.0);
    }

    #[test]
    fn spice_suffixes_in_expr() {
        assert_eq!(eval("1k"), 1000.0);
        assert_eq!(eval("1k + 500"), 1500.0);
        assert_eq!(eval("2.5n"), 2.5e-9);
    }

    #[test]
    fn step_functions() {
        assert_eq!(eval("u(1)"), 1.0);
        assert_eq!(eval("u(-1)"), 0.0);
        assert_eq!(eval("uramp(2)"), 2.0);
        assert_eq!(eval("uramp(-1)"), 0.0);
    }

    #[test]
    fn predicate_functions() {
        assert_eq!(eval("eq0(0)"), 1.0);
        assert_eq!(eval("eq0(1)"), 0.0);
        assert_eq!(eval("gt0(1)"), 1.0);
        assert_eq!(eval("lt0(-1)"), 1.0);
    }
}
