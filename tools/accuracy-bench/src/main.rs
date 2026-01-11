//! Accuracy benchmark for fixed_analytics.
//!
//! Run with: cargo run --release
//! Compare: cargo run --release -- --baseline path/to/baseline.json

use accuracy_bench::{
    build_registry, readme, report::Report, sampling::SampleStrategy, test_function,
};
use rayon::prelude::*;
use std::{env, fs, path::Path, process};

const README_PATH: &str = "../../README.md";

fn main() {
    let args: Vec<String> = env::args().collect();
    let baseline_path = args
        .iter()
        .position(|a| a == "--baseline")
        .and_then(|i| args.get(i + 1))
        .map(String::as_str);

    let strategy = SampleStrategy::thorough();

    eprintln!("--- fixed_analytics accuracy benchmark ---");
    eprintln!(
        "Points per function: ~{}",
        strategy.grid_points + strategy.random_points + strategy.boundary_points * 2
    );

    let registry = build_registry();
    eprintln!("Testing {} functions...\n", registry.len());

    let results: Vec<_> = registry
        .par_iter()
        .map(|f| {
            eprintln!("  {}", f.name());
            test_function(f.as_ref(), &strategy)
        })
        .collect();

    let report = Report::new(results);

    // Save JSON report
    fs::create_dir_all("reports").ok();
    let json_path = format!("reports/accuracy-{}.json", report.timestamp);
    fs::write(&json_path, report.to_json()).expect("Failed to write report");
    eprintln!("Report saved: {json_path}");

    // Determine README path (handle running from different directories)
    let readme_path = find_readme_path();

    if let Some(baseline_path) = baseline_path {
        // CI mode: verify README and compare to baseline
        let mut all_passed = true;

        // Verify README is up-to-date
        if let Some(ref path) = readme_path {
            eprintln!("\nVerifying README accuracy section...");
            match readme::verify_readme(path, &report.results) {
                Ok(()) => eprintln!("README: OK"),
                Err(e) => {
                    eprintln!("README: FAILED\n{e}");
                    all_passed = false;
                }
            }
        } else {
            eprintln!("\nWarning: Could not find README.md to verify");
        }

        // Compare to baseline
        let baseline_passed = compare_and_report(&report, baseline_path);
        if !baseline_passed {
            all_passed = false;
        }

        process::exit(if all_passed { 0 } else { 1 });
    } else {
        // Local mode: update README and print table
        if let Some(ref path) = readme_path {
            match readme::update_readme(path, &report.results) {
                Ok(true) => eprintln!("README.md updated with latest accuracy data"),
                Ok(false) => eprintln!("README.md already up-to-date"),
                Err(e) => eprintln!("Warning: Could not update README: {e}"),
            }
        }

        report.print_table();
    }
}

/// Find the README.md file, checking multiple possible locations.
fn find_readme_path() -> Option<String> {
    let candidates = [
        README_PATH,
        "README.md",
        "../README.md",
        "../../README.md",
        "../../../README.md",
    ];

    for candidate in candidates {
        if Path::new(candidate).exists()
            && let Ok(content) = fs::read_to_string(candidate)
            && content.contains("fixed_analytics")
        {
            return Some(candidate.to_string());
        }
    }

    None
}

fn compare_and_report(current: &Report, baseline_path: &str) -> bool {
    let baseline_json = match fs::read_to_string(baseline_path) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("Failed to read baseline: {e}");
            return false;
        }
    };

    let baseline: Report = match serde_json::from_str(&baseline_json) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Failed to parse baseline: {e}");
            return false;
        }
    };

    println!("\n================================================================================");
    println!("  ACCURACY COMPARISON");
    println!("================================================================================\n");

    println!(
        "{:<12} {:>14} {:>14} {:>14} {:>8}",
        "Function", "Baseline", "Current", "Delta", "Status"
    );
    println!(
        "{:<12} {:>14} {:>14} {:>14} {:>8}",
        "", "(rel_mean)", "(rel_mean)", "", ""
    );
    println!("{}", "-".repeat(70));

    let mut all_passed = true;

    for current_fn in &current.results {
        let baseline_fn = baseline.results.iter().find(|b| b.name == current_fn.name);

        let Some(baseline_fn) = baseline_fn else {
            println!(
                "{:<12} {:>14} {:>14.6e} {:>14} {:>8}",
                current_fn.name, "NEW", current_fn.i16f16.rel_mean, "-", "?"
            );
            continue;
        };

        // Check I16F16
        let (passed_16, status_16) =
            check_regression(baseline_fn.i16f16.rel_mean, current_fn.i16f16.rel_mean);
        if !passed_16 {
            all_passed = false;
        }

        let delta_16 = current_fn.i16f16.rel_mean - baseline_fn.i16f16.rel_mean;
        println!(
            "{:<12} {:>14.6e} {:>14.6e} {:>+14.6e} {:>8}",
            format!("{} I16", current_fn.name),
            baseline_fn.i16f16.rel_mean,
            current_fn.i16f16.rel_mean,
            delta_16,
            status_16
        );

        // Check I32F32
        let (passed_32, status_32) =
            check_regression(baseline_fn.i32f32.rel_mean, current_fn.i32f32.rel_mean);
        if !passed_32 {
            all_passed = false;
        }

        let delta_32 = current_fn.i32f32.rel_mean - baseline_fn.i32f32.rel_mean;
        println!(
            "{:<12} {:>14.6e} {:>14.6e} {:>+14.6e} {:>8}",
            format!("{} I32", current_fn.name),
            baseline_fn.i32f32.rel_mean,
            current_fn.i32f32.rel_mean,
            delta_32,
            status_32
        );
    }

    println!("{}", "-".repeat(70));

    if all_passed {
        println!("\nResult: PASSED (no regressions)\n");
    } else {
        println!("\nResult: FAILED (regression detected)\n");
    }

    all_passed
}

fn check_regression(baseline: f64, current: f64) -> (bool, &'static str) {
    // Allow 0.1% tolerance to avoid floating-point noise triggering false regressions
    let tolerance = baseline * 0.001;
    if current > baseline + tolerance {
        (false, "REGRESS")
    } else if current < baseline - tolerance {
        (true, "IMPROVE")
    } else {
        (true, "SAME")
    }
}
