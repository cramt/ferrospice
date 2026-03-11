//! Batch-mode output formatter matching ngspice `--batch` text output.
//!
//! Produces text output that, after filtering with the check.sh FILTER regex,
//! matches the reference `.out` files from `ngspice-upstream/tests/`.

use thevenin_types::{Analysis, Expr, Item, Netlist, SimPlot, SimResult, SimVector};

/// Parsed `.print` directive: analysis type keyword + list of output variable names.
#[derive(Debug, Clone)]
pub struct PrintDirective {
    pub analysis: String,
    pub vars: Vec<String>,
}

/// Format a simulation result in ngspice batch-mode text format.
pub fn format_batch_output(netlist: &Netlist, result: &SimResult) -> String {
    let mut out = String::new();
    let title = &netlist.title;
    let (temp, tnom) = netlist_temp_tnom(netlist);

    // Circuit header
    out.push_str(&format!("\nCircuit: {title}\n"));

    out.push_str(&format!(
        "\nDoing analysis at TEMP = {temp:.6} and TNOM = {tnom:.6}\n"
    ));

    let prints = parse_print_directives(netlist);
    let analyses = collect_analyses(netlist);

    // Initial transient solution (for .tran analyses)
    let has_tran = analyses.iter().any(|a| matches!(a, Analysis::Tran { .. }));
    if has_tran {
        format_initial_tran_solution(&mut out, result);
    }

    // Transfer function output (for .tf)
    for analysis in &analyses {
        if let Analysis::Tf { output, input } = analysis {
            format_tf_output(&mut out, result, output, input);
        }
    }

    // Data tables for each plot (skip OP plots - they're only used for initial solution)
    for plot in &result.plots {
        let plot_type = plot_analysis_type(&plot.name);
        if plot_type == "op" {
            continue;
        }

        let matching_prints: Vec<&PrintDirective> = prints
            .iter()
            .filter(|p| p.analysis.eq_ignore_ascii_case(&plot_type))
            .collect();

        if !matching_prints.is_empty() {
            for p in &matching_prints {
                format_print_table(&mut out, title, plot, p);
            }
        } else if plot_type == "pz" || plot_type == "sens" {
            format_print_all_table(&mut out, title, plot, &plot_type);
        }
    }

    out
}

/// Get TEMP and TNOM from netlist options.
fn netlist_temp_tnom(netlist: &Netlist) -> (f64, f64) {
    let mut temp = 27.0_f64;
    let mut tnom = 27.0_f64;
    for item in &netlist.items {
        match item {
            Item::Temp(t) => temp = *t,
            Item::Options(params) => {
                for p in params {
                    if let Expr::Num(v) = &p.value
                        && p.name.eq_ignore_ascii_case("TNOM")
                    {
                        tnom = *v;
                    }
                }
            }
            _ => {}
        }
    }
    (temp + 273.15, tnom + 273.15)
}

/// Parse `.print` and `.plot` directives from Raw items.
fn parse_print_directives(netlist: &Netlist) -> Vec<PrintDirective> {
    let mut directives = Vec::new();
    for item in &netlist.items {
        if let Item::Raw(line) = item {
            let trimmed = line.trim();
            let lower = trimmed.to_lowercase();
            if (lower.starts_with(".print") || lower.starts_with(".plot"))
                && let Some(d) = parse_single_print(trimmed)
            {
                directives.push(d);
            }
        }
    }
    directives
}

/// Parse a single `.print` or `.plot` line.
fn parse_single_print(line: &str) -> Option<PrintDirective> {
    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 2 {
        return None;
    }
    let analysis = parts.get(1)?.to_lowercase();

    if parts.len() < 3 {
        return Some(PrintDirective {
            analysis,
            vars: vec!["all".to_string()],
        });
    }

    let vars: Vec<String> = parts[2..]
        .iter()
        .flat_map(|s| s.split(','))
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_lowercase())
        .collect();

    Some(PrintDirective { analysis, vars })
}

/// Collect all Analysis items from the netlist.
fn collect_analyses(netlist: &Netlist) -> Vec<Analysis> {
    netlist
        .items
        .iter()
        .filter_map(|item| {
            if let Item::Analysis(a) = item {
                Some(a.clone())
            } else {
                None
            }
        })
        .collect()
}

/// Determine the analysis type from a plot name like "op1", "tran1", "ac1", "dc1".
fn plot_analysis_type(name: &str) -> String {
    name.trim_end_matches(|c: char| c.is_ascii_digit())
        .to_lowercase()
}

/// Format the initial transient solution (OP at t=0).
fn format_initial_tran_solution(out: &mut String, result: &SimResult) {
    let op_plot = result.plots.iter().find(|p| p.name.starts_with("op"));
    if let Some(plot) = op_plot {
        out.push_str("\nInitial Transient Solution\n");
        out.push_str("--------------------------\n\n");
        out.push_str("Node                                   Voltage\n");
        out.push_str("----                                   -------\n");

        for vec in &plot.vecs {
            if !vec.real.is_empty() {
                let val = vec.real[0];
                // Use bare node name (strip "v()" wrapper)
                let display_name = strip_v_wrapper(&vec.name);
                out.push_str(&format!(
                    "{:<40}{:>12}\n",
                    display_name,
                    format_value_compact(val)
                ));
            }
        }
        out.push_str("\n\n");
    }
}

/// Strip "v()" wrapper from node name: "v(1)" -> "1", "v1#branch" unchanged.
fn strip_v_wrapper(name: &str) -> String {
    if name.starts_with("v(") && name.ends_with(')') {
        name[2..name.len() - 1].to_string()
    } else {
        name.to_string()
    }
}

/// Format a value as integer if round, otherwise as compact decimal.
fn format_value_compact(val: f64) -> String {
    if val == 0.0 {
        "0".to_string()
    } else if val.fract() == 0.0 && val.abs() < 1e15 {
        format!("{}", val as i64)
    } else {
        // Use enough precision to distinguish values
        format!("{val}")
    }
}

/// Format a float in ngspice-compatible scientific notation: always sign + 2-digit exponent.
/// E.g., "1.234567e+00", "-5.678900e-04"
fn format_sci(val: f64) -> String {
    if val == 0.0 {
        return "0.000000e+00".to_string();
    }
    let abs_val = val.abs();
    let exp = abs_val.log10().floor() as i32;
    let mantissa = val / 10.0_f64.powi(exp);

    // Handle edge case where mantissa rounds to 10.0
    if mantissa.abs() >= 9.9999995 {
        let exp = exp + 1;
        let mantissa = val / 10.0_f64.powi(exp);
        return format!("{:.6e}", mantissa * 10.0_f64.powi(exp)).replace(|_: char| false, ""); // fallback
    }

    let sign = if exp >= 0 { '+' } else { '-' };
    let exp_abs = exp.unsigned_abs();
    format!("{mantissa:.6}e{sign}{exp_abs:02}")
}

/// Format transfer function output.
fn format_tf_output(out: &mut String, result: &SimResult, _output: &str, _input: &str) {
    let tf_plot = result.plots.iter().find(|p| p.name.starts_with("tf"));
    if let Some(plot) = tf_plot {
        out.push_str("Transfer function information:\n");
        for vec in &plot.vecs {
            if !vec.real.is_empty() {
                let val = vec.real[0];
                out.push_str(&format!("{} = {}\n", vec.name, format_sci(val)));
            }
        }
        out.push('\n');
    }
}

/// Format a .print data table.
fn format_print_table(out: &mut String, title: &str, plot: &SimPlot, print_dir: &PrintDirective) {
    let resolved_vars = resolve_print_vars(plot, &print_dir.vars);
    if resolved_vars.is_empty() {
        return;
    }

    let n_rows = resolved_vars[0].1.len();
    if n_rows == 0 {
        return;
    }

    out.push_str(&format!("\nNo. of Data Rows : {n_rows}\n"));

    // Title line (centered in 80-char field)
    let title_padded = format!("{:^80}", title);
    out.push_str(&title_padded);
    out.push('\n');

    // Column headers + separators
    format_table_header(out, &resolved_vars);

    // Data rows with pagination
    let page_size = 55; // ngspice paginates around 55 rows
    for i in 0..n_rows {
        if i > 0 && i % page_size == 0 {
            out.push('\n');
            // Repeat column headers for new page
            format_table_header(out, &resolved_vars);
        }
        out.push_str(&format!("{i}\t"));
        for (_, values) in &resolved_vars {
            if i < values.len() {
                out.push_str(&format!("{}\t", format_sci(values[i])));
            }
        }
        out.push('\n');
    }
    out.push('\n');
}

/// Format column header and separators for a data table.
fn format_table_header(out: &mut String, vars: &[(String, Vec<f64>)]) {
    out.push_str(&"-".repeat(80));
    out.push('\n');
    out.push_str("Index\t");
    for (name, _) in vars {
        out.push_str(&format!("{name}\t"));
    }
    out.push('\n');
    out.push_str(&"-".repeat(80));
    out.push('\n');
}

/// Format a table with all vectors (for pz/sens "all").
fn format_print_all_table(out: &mut String, title: &str, plot: &SimPlot, _plot_type: &str) {
    let vecs_with_data: Vec<(&str, &[f64])> = plot
        .vecs
        .iter()
        .filter(|v| !v.real.is_empty())
        .map(|v| (v.name.as_str(), v.real.as_slice()))
        .collect();

    if vecs_with_data.is_empty() {
        return;
    }

    let n_rows = vecs_with_data[0].1.len();
    let cols_per_page = 4;

    for chunk_start in (0..vecs_with_data.len()).step_by(cols_per_page) {
        let chunk_end = (chunk_start + cols_per_page).min(vecs_with_data.len());
        let chunk = &vecs_with_data[chunk_start..chunk_end];

        let title_padded = format!("{:^80}", title);
        out.push_str(&title_padded);
        out.push('\n');

        out.push_str(&"-".repeat(80));
        out.push('\n');
        out.push_str("Index\t");
        for (name, _) in chunk {
            out.push_str(&format!("{name}\t"));
        }
        out.push('\n');
        out.push_str(&"-".repeat(80));
        out.push('\n');

        for i in 0..n_rows {
            out.push_str(&format!("{i}\t"));
            for (_, values) in chunk {
                if i < values.len() {
                    out.push_str(&format!("{}\t", format_sci(values[i])));
                }
            }
            out.push('\n');
        }
        out.push('\n');
    }
}

/// Resolve print variable names to actual data vectors.
fn resolve_print_vars(plot: &SimPlot, var_names: &[String]) -> Vec<(String, Vec<f64>)> {
    let mut result = Vec::new();

    let sweep_vec = find_sweep_vector(plot);
    if let Some(sv) = &sweep_vec {
        result.push((sv.name.clone(), sv.real.clone()));
    }

    for var in var_names {
        if var == "all" {
            for v in &plot.vecs {
                if sweep_vec.as_ref().is_some_and(|sv| sv.name == v.name) {
                    continue;
                }
                if !v.real.is_empty() {
                    result.push((v.name.clone(), v.real.clone()));
                } else if !v.complex.is_empty() {
                    let mags: Vec<f64> = v
                        .complex
                        .iter()
                        .map(|c| (c.re * c.re + c.im * c.im).sqrt())
                        .collect();
                    result.push((v.name.clone(), mags));
                }
            }
            return result;
        }

        if let Some((name, data)) = resolve_single_var(plot, var) {
            result.push((name, data));
        }
    }

    result
}

/// Find the sweep/independent variable in a plot.
fn find_sweep_vector(plot: &SimPlot) -> Option<&SimVector> {
    let sweep_names = ["time", "frequency", "v-sweep", "i-sweep"];
    for name in &sweep_names {
        if let Some(v) = plot.vecs.iter().find(|v| v.name == *name) {
            return Some(v);
        }
    }
    None
}

/// Resolve a single print variable to data.
fn resolve_single_var(plot: &SimPlot, var: &str) -> Option<(String, Vec<f64>)> {
    let var_lower = var.to_lowercase();

    if let Some(inner) = strip_func(&var_lower, "db") {
        let vec = find_vec_by_name(plot, &inner)?;
        if !vec.complex.is_empty() {
            let db_data: Vec<f64> = vec
                .complex
                .iter()
                .map(|c| 20.0 * (c.re * c.re + c.im * c.im).sqrt().log10())
                .collect();
            return Some((format!("db({inner})"), db_data));
        }
        let base = resolve_base_var(plot, &inner)?;
        let db_data: Vec<f64> = base.iter().map(|v| 20.0 * v.abs().log10()).collect();
        return Some((format!("db({inner})"), db_data));
    }
    if let Some(inner) = strip_func(&var_lower, "ph") {
        let vec = find_vec_by_name(plot, &inner)?;
        if !vec.complex.is_empty() {
            let ph_data: Vec<f64> = vec
                .complex
                .iter()
                .map(|c| c.im.atan2(c.re).to_degrees())
                .collect();
            return Some((format!("ph({inner})"), ph_data));
        }
        return None;
    }
    if let Some(inner) = strip_func(&var_lower, "abs") {
        let base = resolve_base_var(plot, &inner)?;
        let abs_data: Vec<f64> = base.iter().map(|v| v.abs()).collect();
        return Some((format!("abs({inner})"), abs_data));
    }

    if let Some(inner) = var_lower.strip_prefix('-') {
        let base = resolve_base_var(plot, inner)?;
        let neg_data: Vec<f64> = base.iter().map(|v| -v).collect();
        return Some((format!("-{inner}"), neg_data));
    }

    let data = resolve_base_var(plot, &var_lower)?;
    Some((var_lower, data))
}

/// Resolve a base variable name to data.
fn resolve_base_var(plot: &SimPlot, name: &str) -> Option<Vec<f64>> {
    let vec = find_vec_by_name(plot, name)?;
    if !vec.real.is_empty() {
        Some(vec.real.clone())
    } else if !vec.complex.is_empty() {
        Some(
            vec.complex
                .iter()
                .map(|c| (c.re * c.re + c.im * c.im).sqrt())
                .collect(),
        )
    } else {
        None
    }
}

/// Find a SimVector by name, handling naming conventions.
fn find_vec_by_name<'a>(plot: &'a SimPlot, name: &str) -> Option<&'a SimVector> {
    let name_lower = name.to_lowercase();

    // Direct match
    if let Some(v) = plot
        .vecs
        .iter()
        .find(|v| v.name.to_lowercase() == name_lower)
    {
        return Some(v);
    }

    // i(source) -> source#branch
    if let Some(inner) = strip_func(&name_lower, "i") {
        let branch_name = format!("{inner}#branch");
        if let Some(v) = plot
            .vecs
            .iter()
            .find(|v| v.name.to_lowercase() == branch_name)
        {
            return Some(v);
        }
    }

    None
}

/// Strip a function wrapper: "db(v(2))" with func="db" returns "v(2)".
fn strip_func(s: &str, func: &str) -> Option<String> {
    let prefix = format!("{func}(");
    if s.starts_with(&prefix) && s.ends_with(')') {
        Some(s[prefix.len()..s.len() - 1].to_string())
    } else {
        None
    }
}

/// Apply the check.sh filter to text output.
pub fn filter_output(text: &str) -> String {
    let normalized = normalize_exponents(text);

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
    ];

    let mut lines: Vec<&str> = Vec::new();
    for line in normalized.lines() {
        if line.starts_with("binary raw file")
            || (line.starts_with("ngspice") && line.contains("done"))
        {
            continue;
        }
        if line.contains("Operating") {
            continue;
        }
        if filter_patterns.iter().any(|pat| line.contains(pat)) {
            continue;
        }
        lines.push(line);
    }

    lines.join("\n")
}

/// Normalize Windows-style 3-digit exponents to 2-digit.
fn normalize_exponents(text: &str) -> String {
    let mut result = String::with_capacity(text.len());
    let bytes = text.as_bytes();
    let len = bytes.len();
    let mut i = 0;

    while i < len {
        if i + 4 < len
            && (bytes[i] == b'.' || bytes[i].is_ascii_digit())
            && (bytes[i + 1] == b'e' || bytes[i + 1] == b'E')
        {
            let sign_offset = if i + 2 < len && (bytes[i + 2] == b'+' || bytes[i + 2] == b'-') {
                1
            } else {
                0
            };
            let zero_pos = i + 2 + sign_offset;
            if zero_pos + 2 < len
                && bytes[zero_pos] == b'0'
                && bytes[zero_pos + 1].is_ascii_digit()
                && bytes[zero_pos + 2].is_ascii_digit()
            {
                result.push(bytes[i] as char);
                result.push(bytes[i + 1] as char);
                if sign_offset == 1 {
                    result.push(bytes[i + 2] as char);
                }
                result.push(bytes[zero_pos + 1] as char);
                result.push(bytes[zero_pos + 2] as char);
                i = zero_pos + 3;
                continue;
            }
        }
        result.push(bytes[i] as char);
        i += 1;
    }

    result
}

/// Compare two outputs after filtering, using diff -B -w semantics with numeric tolerance.
///
/// Non-numeric tokens are compared exactly. Numeric tokens (parseable as f64) are
/// compared with a relative tolerance of 1e-4 and absolute tolerance of 1e-15.
///
/// When line counts differ, falls back to interpolation-aware comparison for
/// time-series data sections (transient, DC sweep, etc.), linearly interpolating
/// the actual data at the expected data's independent variable points.
pub fn compare_filtered(expected: &str, actual: &str) -> Result<(), String> {
    let expected_filtered = filter_output(expected);
    let actual_filtered = filter_output(actual);

    let norm_expected = normalize_for_diff(&expected_filtered);
    let norm_actual = normalize_for_diff(&actual_filtered);

    // Try exact line-by-line comparison first.
    if norm_expected.len() == norm_actual.len() {
        let all_match = norm_expected
            .iter()
            .zip(norm_actual.iter())
            .all(|(e, a)| lines_match_approx(e, a));
        if all_match {
            return Ok(());
        }
    }

    // Fall back to interpolation-aware comparison.
    compare_with_interpolation(&norm_expected, &norm_actual)
}

/// Check if a normalized line looks like a data row: starts with an integer index
/// followed by numeric values (e.g. "5 1.234e-09 -2.000e+00").
fn is_data_row(line: &str) -> bool {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 2 {
        return false;
    }
    // First token must be a non-negative integer (row index).
    if tokens[0].parse::<u64>().is_err() {
        return false;
    }
    // Second token must be a float (time/sweep variable).
    tokens[1].parse::<f64>().is_ok()
}

/// Parse a data row into (index, independent_var, dependent_values).
fn parse_data_row(line: &str) -> Option<(u64, f64, Vec<f64>)> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < 2 {
        return None;
    }
    let idx = tokens[0].parse::<u64>().ok()?;
    let indep = tokens[1].parse::<f64>().ok()?;
    let deps: Vec<f64> = tokens[2..]
        .iter()
        .filter_map(|t| t.parse::<f64>().ok())
        .collect();
    if deps.len() + 2 != tokens.len() {
        return None; // Some token wasn't numeric.
    }
    Some((idx, indep, deps))
}

/// Linearly interpolate value at `x` given sorted data points `xs` and `ys`.
fn lerp_at(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    debug_assert_eq!(xs.len(), ys.len());
    if xs.is_empty() {
        return 0.0;
    }
    if x <= xs[0] {
        return ys[0];
    }
    if x >= xs[xs.len() - 1] {
        return ys[ys.len() - 1];
    }
    // Binary search for the interval.
    let i = match xs.binary_search_by(|v| v.partial_cmp(&x).unwrap()) {
        Ok(i) => return ys[i], // Exact match.
        Err(i) => i,
    };
    let x0 = xs[i - 1];
    let x1 = xs[i];
    let y0 = ys[i - 1];
    let y1 = ys[i];
    let frac = (x - x0) / (x1 - x0);
    y0 + frac * (y1 - y0)
}

/// Compare with interpolation: separate text and data lines, compare text exactly,
/// interpolate data for comparison.
fn compare_with_interpolation(
    norm_expected: &[String],
    norm_actual: &[String],
) -> Result<(), String> {
    // Separate text lines and data rows.
    let exp_text: Vec<&str> = norm_expected
        .iter()
        .filter(|l| !is_data_row(l))
        .map(|s| s.as_str())
        .collect();
    let act_text: Vec<&str> = norm_actual
        .iter()
        .filter(|l| !is_data_row(l))
        .map(|s| s.as_str())
        .collect();

    // Compare text lines.
    if exp_text.len() != act_text.len() {
        return Err(format_diff(norm_expected, norm_actual));
    }
    for (e, a) in exp_text.iter().zip(act_text.iter()) {
        if !lines_match_approx(e, a) {
            return Err(format_diff(norm_expected, norm_actual));
        }
    }

    // Parse data rows.
    let exp_data: Vec<(u64, f64, Vec<f64>)> = norm_expected
        .iter()
        .filter_map(|l| parse_data_row(l))
        .collect();
    let act_data: Vec<(u64, f64, Vec<f64>)> = norm_actual
        .iter()
        .filter_map(|l| parse_data_row(l))
        .collect();

    if exp_data.is_empty() && act_data.is_empty() {
        return Ok(());
    }
    if exp_data.is_empty() || act_data.is_empty() {
        return Err(format_diff(norm_expected, norm_actual));
    }

    // Check that both have the same number of dependent columns.
    let n_deps = exp_data[0].2.len();
    if act_data[0].2.len() != n_deps {
        return Err(format_diff(norm_expected, norm_actual));
    }

    // Extract independent variable arrays.
    let exp_x: Vec<f64> = exp_data.iter().map(|(_, x, _)| *x).collect();
    let act_x: Vec<f64> = act_data.iter().map(|(_, x, _)| *x).collect();

    // For each dependent column, build actual arrays and interpolate at expected points.
    for col in 0..n_deps {
        let act_y: Vec<f64> = act_data.iter().map(|(_, _, deps)| deps[col]).collect();

        for (i, exp_row) in exp_data.iter().enumerate() {
            let exp_val = exp_row.2[col];
            let interp_val = lerp_at(exp_x[i], &act_x, &act_y);

            let abs_diff = (exp_val - interp_val).abs();
            let rel_tol = 1e-4 * exp_val.abs().max(interp_val.abs());
            if abs_diff > rel_tol.max(1e-15) {
                return Err(format!(
                    "Interpolation mismatch at x={:.6e}, col {}: expected {:.6e}, got {:.6e} (diff={:.6e})\n{}",
                    exp_x[i],
                    col,
                    exp_val,
                    interp_val,
                    abs_diff,
                    format_diff(norm_expected, norm_actual),
                ));
            }
        }
    }

    Ok(())
}

/// Check if two lines match, allowing numeric tolerance.
fn lines_match_approx(expected: &str, actual: &str) -> bool {
    let exp_tokens: Vec<&str> = expected.split_whitespace().collect();
    let act_tokens: Vec<&str> = actual.split_whitespace().collect();

    if exp_tokens.len() != act_tokens.len() {
        return false;
    }

    for (e, a) in exp_tokens.iter().zip(act_tokens.iter()) {
        if e == a {
            continue;
        }
        match (e.parse::<f64>(), a.parse::<f64>()) {
            (Ok(ev), Ok(av)) => {
                let abs_diff = (ev - av).abs();
                let rel_tol = 1e-4 * ev.abs().max(av.abs());
                if abs_diff > rel_tol.max(1e-15) {
                    return false;
                }
            }
            _ => return false,
        }
    }
    true
}

/// Format a diff between expected and actual for error reporting.
fn format_diff(expected: &[String], actual: &[String]) -> String {
    let mut diff = String::new();
    diff.push_str("=== Expected (filtered) ===\n");
    for line in expected.iter() {
        diff.push_str(line);
        diff.push('\n');
    }
    diff.push_str("=== Actual (filtered) ===\n");
    for line in actual.iter() {
        diff.push_str(line);
        diff.push('\n');
    }
    diff
}

/// Normalize text for diff comparison: collapse whitespace, remove blank lines.
fn normalize_for_diff(text: &str) -> Vec<String> {
    text.lines()
        .map(|line| line.split_whitespace().collect::<Vec<_>>().join(" "))
        .filter(|line| !line.is_empty())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_removes_expected_patterns() {
        let input = "Circuit: test\nDoing analysis\n0\t1.000e+00\n---\nIndex\ttime\n";
        let filtered = filter_output(input);
        assert!(filtered.contains("0\t1.000e+00"));
        assert!(!filtered.contains("Circuit"));
        assert!(!filtered.contains("Doing"));
        assert!(!filtered.contains("---"));
        assert!(!filtered.contains("Index"));
    }

    #[test]
    fn test_normalize_exponents() {
        assert_eq!(normalize_exponents("1.0e+001"), "1.0e+01");
        assert_eq!(normalize_exponents("2.5e-003"), "2.5e-03");
        assert_eq!(normalize_exponents("3.0E+010"), "3.0E+10");
        assert_eq!(normalize_exponents("normal text"), "normal text");
    }

    #[test]
    fn test_compare_filtered_matching() {
        let a = "Circuit: test\n0\t1.000000e+00\t";
        let b = "Circuit: test\n0\t1.000000e+00\t";
        assert!(compare_filtered(a, b).is_ok());
    }

    #[test]
    fn test_compare_filtered_ignores_blanks_and_whitespace() {
        let a = "0\t1.000000e+00\n\n";
        let b = "0  1.000000e+00\n";
        assert!(compare_filtered(a, b).is_ok());
    }

    #[test]
    fn test_parse_print_directive() {
        let d = parse_single_print(".print dc v(2) i(vs)").unwrap();
        assert_eq!(d.analysis, "dc");
        assert_eq!(d.vars, vec!["v(2)", "i(vs)"]);
    }

    #[test]
    fn test_parse_print_with_commas() {
        let d = parse_single_print(".PRINT AC V(2), I(VR1)").unwrap();
        assert_eq!(d.analysis, "ac");
        assert_eq!(d.vars, vec!["v(2)", "i(vr1)"]);
    }

    #[test]
    fn test_strip_func() {
        assert_eq!(strip_func("db(v(2))", "db"), Some("v(2)".to_string()));
        assert_eq!(
            strip_func("ph(i(vmeas))", "ph"),
            Some("i(vmeas)".to_string())
        );
        assert_eq!(strip_func("v(2)", "db"), None);
    }

    #[test]
    fn test_format_sci() {
        assert_eq!(format_sci(0.0), "0.000000e+00");
        assert_eq!(format_sci(1.0), "1.000000e+00");
        assert_eq!(format_sci(-1e-4), "-1.000000e-04");
        assert_eq!(format_sci(3e-13), "3.000000e-13");
        assert_eq!(format_sci(1.536e-10), "1.536000e-10");
    }

    #[test]
    fn test_strip_v_wrapper() {
        assert_eq!(strip_v_wrapper("v(1)"), "1");
        assert_eq!(strip_v_wrapper("v1#branch"), "v1#branch");
        assert_eq!(strip_v_wrapper("v(out)"), "out");
    }
}
