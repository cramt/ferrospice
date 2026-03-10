//! `.lib` file processing — loads library files and extracts named sections.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use thevenin_types::{Item, Netlist};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum LibError {
    #[error("library file not found: {0}")]
    FileNotFound(PathBuf),
    #[error("library section not found: {section} in {file}")]
    SectionNotFound { file: String, section: String },
    #[error("I/O error reading {path}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("parse error in library file {file}: {source}")]
    Parse {
        file: String,
        source: thevenin_types::ParseError,
    },
}

/// Process all `.lib` directives in a netlist, loading files from `base_dir`.
///
/// Recursively processes library includes. `.lib 'file' section` loads the
/// named section from the file. `.LIB section` / `.ENDL section` markers
/// within library files delimit sections.
pub fn process_libs(netlist: &mut Netlist, base_dir: &Path) -> Result<(), LibError> {
    let mut cache: BTreeMap<PathBuf, String> = BTreeMap::new();
    process_items(&mut netlist.items, base_dir, &mut cache, 0)
}

fn process_items(
    items: &mut Vec<Item>,
    base_dir: &Path,
    cache: &mut BTreeMap<PathBuf, String>,
    depth: usize,
) -> Result<(), LibError> {
    if depth > 20 {
        return Ok(()); // prevent infinite recursion
    }

    let mut i = 0;
    while i < items.len() {
        match &items[i] {
            Item::Lib {
                file,
                entry: Some(entry),
            } => {
                let file = file.clone();
                let entry = entry.clone();
                let lib_path = base_dir.join(&file);

                // Load file content (from cache or disk)
                let content = if let Some(cached) = cache.get(&lib_path) {
                    cached.clone()
                } else {
                    let content = std::fs::read_to_string(&lib_path).map_err(|e| LibError::Io {
                        path: lib_path.clone(),
                        source: e,
                    })?;
                    cache.insert(lib_path.clone(), content.clone());
                    content
                };

                // Extract the named section
                let section_items = extract_section(&content, &entry, &file)?;

                // Remove the .lib directive and insert section items
                items.remove(i);
                let section_len = section_items.len();
                for (j, item) in section_items.into_iter().enumerate() {
                    items.insert(i + j, item);
                }

                // Determine the lib file's directory for nested includes
                let lib_dir = lib_path.parent().unwrap_or(base_dir);

                // Recursively process the inserted items
                let end = i + section_len;
                let mut sub_items: Vec<Item> = items.drain(i..end).collect();
                process_items(&mut sub_items, lib_dir, cache, depth + 1)?;
                for (j, item) in sub_items.into_iter().enumerate() {
                    items.insert(i + j, item);
                }
                // Don't increment i — process the newly inserted items
            }
            Item::Subckt(sub) => {
                // Process .lib directives inside subcircuit definitions
                let mut sub_clone = sub.clone();
                process_items(&mut sub_clone.items, base_dir, cache, depth + 1)?;
                items[i] = Item::Subckt(sub_clone);
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }
    Ok(())
}

/// Extract items from a named `.LIB section` / `.ENDL section` block.
fn extract_section(content: &str, section: &str, file: &str) -> Result<Vec<Item>, LibError> {
    let section_upper = section.to_uppercase();

    // Parse the file as a temporary netlist to get items
    // We wrap it with a dummy title and .end
    let wrapped = format!("lib_wrapper\n{content}\n.end\n");
    let parsed = Netlist::parse(&wrapped).map_err(|e| LibError::Parse {
        file: file.to_string(),
        source: e,
    })?;

    // Find the section boundaries
    let mut in_section = false;
    let mut section_items = Vec::new();
    let mut section_depth = 0;

    for item in &parsed.items {
        if let Item::Raw(s) = item {
            let upper = s.trim().to_uppercase();
            // Check for .LIB section_name (section start marker)
            if let Some(rest) = upper.strip_prefix(".LIB") {
                let name = rest.trim();
                if !name.is_empty()
                    && !name.contains('\'')
                    && !name.contains('"')
                    && name == section_upper
                    && !in_section
                {
                    in_section = true;
                    section_depth = 0;
                    continue;
                } else if in_section {
                    section_depth += 1;
                }
            }
            // Check for .ENDL
            if upper.starts_with(".ENDL") {
                if in_section && section_depth == 0 {
                    in_section = false;
                    continue;
                } else if in_section {
                    section_depth -= 1;
                }
            }
        }

        if in_section {
            section_items.push(item.clone());
        }
    }

    if section_items.is_empty() {
        // Try finding the section by scanning raw text (more robust for nested .LIB)
        return extract_section_raw(content, section, file);
    }

    Ok(section_items)
}

/// Fallback: extract section by raw text scanning.
fn extract_section_raw(content: &str, section: &str, file: &str) -> Result<Vec<Item>, LibError> {
    let section_upper = section.to_uppercase();
    let mut in_section = false;
    let mut section_lines = Vec::new();

    for line in content.lines() {
        let trimmed = line.trim();
        let upper = trimmed.to_uppercase();

        if !in_section {
            // Look for `.LIB section_name` (without quotes — section marker, not file include)
            if let Some(rest) = upper.strip_prefix(".LIB") {
                let name = rest.trim();
                if name == section_upper && !name.contains('\'') && !name.contains('"') {
                    in_section = true;
                    continue;
                }
            }
        } else {
            // Check for .ENDL
            if upper.starts_with(".ENDL") {
                in_section = false;
                continue;
            }
            section_lines.push(line);
        }
    }

    if section_lines.is_empty() {
        return Err(LibError::SectionNotFound {
            file: file.to_string(),
            section: section.to_string(),
        });
    }

    // Parse the extracted lines
    let section_text = format!("lib_section\n{}\n.end\n", section_lines.join("\n"));
    let parsed = Netlist::parse(&section_text).map_err(|e| LibError::Parse {
        file: file.to_string(),
        source: e,
    })?;

    Ok(parsed.items)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_simple_section() {
        let content = r"
* comment
.LIB RES
.subckt sub1 n1 n2
R1 n1 n2 2k
.ends
.ENDL RES
";
        let items = extract_section(content, "RES", "test.lib").unwrap();
        assert!(!items.is_empty());
        // Should have the subcircuit
        assert!(items.iter().any(|i| matches!(i, Item::Subckt(_))));
    }

    #[test]
    fn extract_section_not_found() {
        let content = ".LIB FOO\n.ENDL FOO\n";
        let result = extract_section(content, "BAR", "test.lib");
        assert!(result.is_err());
    }
}
