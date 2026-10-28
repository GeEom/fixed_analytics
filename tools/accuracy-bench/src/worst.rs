//! The worst-case table: running the sweep, writing it as text, comparing it
//! with a baseline, and its section of the README.

use crate::sut::types;
use crate::sweep::{ALL, FunStats};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::fmt::Write;

pub const WORST_PATH: &str = "worst.txt";
const MARKER_START: &str = "<!-- WORST_START -->";
const MARKER_END: &str = "<!-- WORST_END -->";
const README_TYPES: [&str; 3] = ["I16F16", "I32F32", "I64F64"];

/// Runs every function on every type, in the order of `ALL` then `types()`.
pub fn sweep() -> Vec<FunStats> {
    let types = types();
    let jobs: Vec<(usize, usize)> = (0..ALL.len())
        .flat_map(|f| (0..types.len()).map(move |t| (f, t)))
        .collect();
    let mut out: Vec<((usize, usize), FunStats)> = jobs
        .par_iter()
        .map(|&(f, t)| ((f, t), (types[t].run)(ALL[f])))
        .collect();
    out.sort_by_key(|(k, _)| *k);
    out.into_iter().map(|(_, s)| s).collect()
}

const HEADER: &str = "# function\ttype\tsamples\tmax_ulps\tmax_rel_ulps\tparity_bad\tparity_n\tmono_bad\tmono_n\terrors\tworst\tworst_rel\texample";

pub fn to_text(stats: &[FunStats]) -> String {
    let mut out = String::new();
    writeln!(out, "{HEADER}").unwrap();
    for s in stats {
        writeln!(
            out,
            "{}\t{}\t{}\t{:e}\t{:e}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            s.fun,
            s.ty,
            s.samples,
            s.max_ulps,
            s.max_rel_ulps,
            s.parity_bad,
            s.parity_n,
            s.mono_bad,
            s.mono_n,
            s.errors,
            s.worst,
            s.worst_rel,
            s.example
        )
        .unwrap();
    }
    out
}

pub fn parse(text: &str) -> BTreeMap<(String, String), FunStats> {
    let mut map = BTreeMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#') && !l.is_empty()) {
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 10 {
            continue;
        }
        let s = FunStats {
            fun: c[0].to_string(),
            ty: c[1].to_string(),
            samples: c[2].parse().unwrap_or(0),
            max_ulps: c[3].parse().unwrap_or(f64::INFINITY),
            max_rel_ulps: c[4].parse().unwrap_or(f64::INFINITY),
            parity_bad: c[5].parse().unwrap_or(u64::MAX),
            parity_n: c[6].parse().unwrap_or(0),
            mono_bad: c[7].parse().unwrap_or(u64::MAX),
            mono_n: c[8].parse().unwrap_or(0),
            errors: c[9].parse().unwrap_or(u64::MAX),
            ..FunStats::default()
        };
        map.insert((s.fun.clone(), s.ty.clone()), s);
    }
    map
}

fn worse(cur: f64, base: f64) -> bool {
    cur > base * (1.0 + 1e-9) + 1e-12
}

/// Prints every row past the crate's stated bound, or that breaks a symmetry
/// or monotonicity, or errs; true if none.
pub fn check_bounds(stats: &[FunStats]) -> bool {
    let limit = f64::from(fixed_analytics::ERROR_ULPS);
    let mut ok = true;
    for s in stats {
        let fun = crate::sweep::ALL.iter().find(|f| f.name() == s.fun);
        let exact = fun.is_some_and(|f| f.rounds_once());
        // A directed form's error reaches toward a whole ulp, which the sweep checks exactly.
        let directed = fun.is_some_and(|f| f.direction().is_some());
        let bound = if exact { 0.5 } else { limit };
        let over = if exact || directed { s.max_ulps > bound } else { s.max_ulps >= bound };
        if over || s.parity_bad > 0 || s.mono_bad > 0 || s.errors > 0 {
            ok = false;
            println!(
                "BOUND {} {}: {:.3e} ulps (bound {bound}), symmetry {}, monotonicity {}, errors {}: {}",
                s.fun, s.ty, s.max_ulps, s.parity_bad, s.mono_bad, s.errors, s.example
            );
        }
    }
    println!("Bounds: {}", if ok { "every function within its bound" } else { "FAILED" });
    ok
}

/// Prints every row where `current` is worse than `baseline`; true if none.
pub fn compare(current: &[FunStats], baseline: &str) -> bool {
    let base = parse(baseline);
    let mut ok = true;
    let mut improved = 0;
    for s in current {
        let Some(b) = base.get(&(s.fun.clone(), s.ty.clone())) else {
            continue;
        };
        let mut why = Vec::new();
        if worse(s.max_ulps, b.max_ulps) {
            why.push(format!("max {:.3e} -> {:.3e} ulps", b.max_ulps, s.max_ulps));
        }
        if worse(s.max_rel_ulps, b.max_rel_ulps) {
            why.push(format!(
                "relative {:.3e} -> {:.3e} ulps",
                b.max_rel_ulps, s.max_rel_ulps
            ));
        }
        if s.parity_bad > b.parity_bad {
            why.push(format!("symmetry {} -> {}", b.parity_bad, s.parity_bad));
        }
        if s.mono_bad > b.mono_bad {
            why.push(format!("monotonicity {} -> {}", b.mono_bad, s.mono_bad));
        }
        if s.errors > b.errors {
            why.push(format!("errors {} -> {}", b.errors, s.errors));
        }
        if why.is_empty() {
            if s.max_ulps < b.max_ulps || s.parity_bad < b.parity_bad || s.mono_bad < b.mono_bad {
                improved += 1;
            }
        } else {
            ok = false;
            println!("REGRESS {} {}: {}", s.fun, s.ty, why.join(", "));
        }
    }
    println!("Worst case: {improved} function-type pairs improve");
    ok
}

fn fmt_ulps(u: f64) -> String {
    if u < 100.0 {
        format!("{u:.2}")
    } else {
        format!("{u:.2e}")
    }
}

pub fn readme_section(stats: &[FunStats]) -> String {
    let n_types = types().len();
    let mut out = String::new();
    writeln!(out, "### Worst-case error\n").unwrap();
    writeln!(
        out,
        "The largest error seen, in ulps of the type, against MPFR on exact inputs. The sweep covers {n_types} layouts from 8 to 128 bits, every value of each 8- and 16-bit one and {} seeded samples of each wider one; tools/accuracy-bench/worst.txt holds every layout.\n",
        crate::sweep::sample_count()
    )
    .unwrap();
    writeln!(out, "| Function | {} |", README_TYPES.join(" | ")).unwrap();
    writeln!(out, "|----------|{}", "------|".repeat(README_TYPES.len())).unwrap();
    let mut by_fun: BTreeMap<usize, Vec<String>> = BTreeMap::new();
    for s in stats {
        if let Some(col) = README_TYPES.iter().position(|t| *t == s.ty) {
            let row = ALL.iter().position(|f| f.name() == s.fun).unwrap_or(usize::MAX);
            let cells = by_fun.entry(row).or_insert_with(|| vec![String::new(); README_TYPES.len()]);
            cells[col] = fmt_ulps(s.max_ulps);
        }
    }
    for (row, cells) in by_fun {
        writeln!(out, "| {} | {} |", ALL[row].name(), cells.join(" | ")).unwrap();
    }
    out
}

fn section_bounds(content: &str) -> Result<(usize, usize), String> {
    let start = content
        .find(MARKER_START)
        .ok_or("README has no <!-- WORST_START --> marker")?;
    let end = content
        .find(MARKER_END)
        .ok_or("README has no <!-- WORST_END --> marker")?;
    if end <= start {
        return Err("WORST_END must follow WORST_START".to_string());
    }
    Ok((start + MARKER_START.len(), end))
}

pub fn update_readme(path: &str, stats: &[FunStats]) -> Result<bool, String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let (a, b) = section_bounds(&content)?;
    let new = format!("{}\n{}{}", &content[..a], readme_section(stats), &content[b..]);
    if new == content {
        return Ok(false);
    }
    std::fs::write(path, new).map_err(|e| e.to_string())?;
    Ok(true)
}

pub fn verify_readme(path: &str, stats: &[FunStats]) -> Result<(), String> {
    let content = std::fs::read_to_string(path).map_err(|e| e.to_string())?;
    let (a, b) = section_bounds(&content)?;
    if content[a..b].trim() == readme_section(stats).trim() {
        Ok(())
    } else {
        Err("README worst-case section is out of date; run `cargo run --release` in tools/accuracy-bench".to_string())
    }
}
