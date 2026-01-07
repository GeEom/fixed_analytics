//! Accuracy benchmark for fixed_analytics.
//!
//! Run with: cargo run --release
//! Compare: cargo run --release -- --baseline path/to/baseline.json

use accuracy_bench::{build_registry, report::Report, sampling::SampleStrategy, test_function};
use rayon::prelude::*;
use std::{env, fs, process};

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

    let results = registry
        .par_iter()
        .map(|f| {
            eprintln!("  {}", f.name());
            test_function(f.as_ref(), &strategy)
        })
        .collect();

    let report = Report::new(results);

    fs::create_dir_all("reports").ok();
    let json_path = format!("reports/accuracy-{}.json", report.timestamp);
    fs::write(&json_path, report.to_json()).expect("Failed to write report");
    eprintln!("Report saved: {json_path}");

    if let Some(path) = baseline_path {
        let passed = compare_and_report(&report, path);
        process::exit(if passed { 0 } else { 1 });
    } else {
        report.print_table();
    }
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
        let (passed_16, status_16) = check_regression(
            baseline_fn.i16f16.rel_mean,
            current_fn.i16f16.rel_mean,
        );
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
        let (passed_32, status_32) = check_regression(
            baseline_fn.i32f32.rel_mean,
            current_fn.i32f32.rel_mean,
        );
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
