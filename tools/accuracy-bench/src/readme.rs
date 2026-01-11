//! README accuracy table generation and validation.

use crate::FunctionResult;
use std::fmt::Write;

const MARKER_START: &str = "<!-- ACCURACY_START -->";
const MARKER_END: &str = "<!-- ACCURACY_END -->";

/// Generate the accuracy section content (without markers).
pub fn generate_accuracy_section(results: &[FunctionResult]) -> String {
    let mut out = String::new();

    writeln!(out, "### Accuracy\n").unwrap();
    writeln!(
        out,
        "Relative error statistics measured against MPFR reference implementations.\n"
    )
    .unwrap();

    // Combined table with both I16F16 and I32F32
    writeln!(
        out,
        "| Function | I16F16 Mean | I16F16 Median | I16F16 P95 | I32F32 Mean | I32F32 Median | I32F32 P95 |"
    )
    .unwrap();
    writeln!(
        out,
        "|----------|-------------|---------------|------------|-------------|---------------|------------|"
    )
    .unwrap();
    for r in results {
        writeln!(
            out,
            "| {} | {:.2e} | {:.2e} | {:.2e} | {:.2e} | {:.2e} | {:.2e} |",
            r.name,
            r.i16f16.rel_mean,
            r.i16f16.rel_p50,
            r.i16f16.rel_p95,
            r.i32f32.rel_mean,
            r.i32f32.rel_p50,
            r.i32f32.rel_p95
        )
        .unwrap();
    }

    out
}

/// Update the README file with new accuracy data.
/// Returns Ok(true) if changes were made, Ok(false) if already up-to-date.
pub fn update_readme(readme_path: &str, results: &[FunctionResult]) -> Result<bool, String> {
    let content = std::fs::read_to_string(readme_path)
        .map_err(|e| format!("Failed to read {readme_path}: {e}"))?;

    let new_section = generate_accuracy_section(results);
    let new_content = replace_section(&content, &new_section)?;

    if new_content == content {
        return Ok(false);
    }

    std::fs::write(readme_path, &new_content)
        .map_err(|e| format!("Failed to write {readme_path}: {e}"))?;

    Ok(true)
}

/// Verify the README accuracy section matches current results.
/// Returns Ok(()) if valid, Err with details if mismatched.
pub fn verify_readme(readme_path: &str, results: &[FunctionResult]) -> Result<(), String> {
    let content = std::fs::read_to_string(readme_path)
        .map_err(|e| format!("Failed to read {readme_path}: {e}"))?;

    let current = extract_section(&content)?;
    let expected = generate_accuracy_section(results);

    let current_values = parse_table_values(&current)?;
    let expected_values = parse_table_values(&expected)?;

    let mut mismatches = Vec::new();

    for (key, exp_val) in &expected_values {
        match current_values.get(key) {
            Some(cur_val) => {
                // Allow 1% tolerance for floating-point formatting differences
                let rel_diff = if *exp_val != 0.0 {
                    (cur_val - exp_val).abs() / exp_val.abs()
                } else {
                    cur_val.abs()
                };
                if rel_diff > 0.01 {
                    mismatches.push(format!(
                        "  {}: README has {:.2e}, expected {:.2e}",
                        key, cur_val, exp_val
                    ));
                }
            }
            None => {
                mismatches.push(format!("  {}: missing from README", key));
            }
        }
    }

    // Check for extra entries in README
    for key in current_values.keys() {
        if !expected_values.contains_key(key) {
            mismatches.push(format!("  {}: unexpected entry in README", key));
        }
    }

    if mismatches.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "README accuracy section is out of date:\n{}\n\nRun `cargo run --release` in tools/accuracy-bench to update.",
            mismatches.join("\n")
        ))
    }
}

fn replace_section(content: &str, new_section: &str) -> Result<String, String> {
    let start_idx = content.find(MARKER_START)
        .ok_or("README missing accuracy section start marker. Add <!-- ACCURACY_START --> where you want the table.")?;
    let end_idx = content.find(MARKER_END)
        .ok_or("README missing accuracy section end marker. Add <!-- ACCURACY_END --> after the start marker.")?;

    if end_idx <= start_idx {
        return Err("ACCURACY_END marker must come after ACCURACY_START".to_string());
    }

    let mut result = String::new();
    result.push_str(&content[..start_idx]);
    result.push_str(MARKER_START);
    result.push('\n');
    result.push_str(new_section);
    result.push_str(MARKER_END);
    result.push_str(&content[end_idx + MARKER_END.len()..]);

    Ok(result)
}

fn extract_section(content: &str) -> Result<String, String> {
    let start_idx = content
        .find(MARKER_START)
        .ok_or("README missing accuracy section start marker")?;
    let end_idx = content
        .find(MARKER_END)
        .ok_or("README missing accuracy section end marker")?;

    let section_start = start_idx + MARKER_START.len();
    Ok(content[section_start..end_idx].to_string())
}

/// Parse all numeric values from the accuracy tables into a map.
/// Keys are like "sin/I16F16/mean", "cos/I32F32/p95", etc.
fn parse_table_values(section: &str) -> Result<std::collections::HashMap<String, f64>, String> {
    let mut values = std::collections::HashMap::new();

    for line in section.lines() {
        let line = line.trim();

        // Skip non-data lines
        if !line.starts_with('|') || line.contains("---") || line.contains("Function") {
            continue;
        }

        // Parse table row: | func | i16f16_mean | i16f16_median | i16f16_p95 | i32f32_mean | i32f32_median | i32f32_p95 |
        let parts: Vec<&str> = line
            .split('|')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .collect();

        if parts.len() >= 7 {
            let func = parts[0];
            if let Ok(v) = parts[1].parse::<f64>() {
                values.insert(format!("{}/I16F16/mean", func), v);
            }
            if let Ok(v) = parts[2].parse::<f64>() {
                values.insert(format!("{}/I16F16/median", func), v);
            }
            if let Ok(v) = parts[3].parse::<f64>() {
                values.insert(format!("{}/I16F16/p95", func), v);
            }
            if let Ok(v) = parts[4].parse::<f64>() {
                values.insert(format!("{}/I32F32/mean", func), v);
            }
            if let Ok(v) = parts[5].parse::<f64>() {
                values.insert(format!("{}/I32F32/median", func), v);
            }
            if let Ok(v) = parts[6].parse::<f64>() {
                values.insert(format!("{}/I32F32/p95", func), v);
            }
        }
    }

    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_table_values() {
        let section = r#"
### Accuracy

Relative error statistics measured against MPFR reference implementations.

| Function | I16F16 Mean | I16F16 Median | I16F16 P95 | I32F32 Mean | I32F32 Median | I32F32 P95 |
|----------|-------------|---------------|------------|-------------|---------------|------------|
| sin | 7.30e-5 | 6.05e-5 | 1.80e-4 | 1.41e-9 | 1.16e-9 | 3.49e-9 |
| cos | 7.96e-5 | 6.44e-5 | 2.03e-4 | 1.50e-9 | 1.20e-9 | 3.60e-9 |
"#;
        let values = parse_table_values(section).unwrap();

        assert!((values["sin/I16F16/mean"] - 7.30e-5).abs() < 1e-10);
        assert!((values["sin/I32F32/mean"] - 1.41e-9).abs() < 1e-14);
        assert!((values["cos/I16F16/p95"] - 2.03e-4).abs() < 1e-10);
    }
}
