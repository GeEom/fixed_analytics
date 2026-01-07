//! Error metrics and statistical analysis.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy)]
pub struct ErrorMeasurement {
    pub absolute: f64,
    pub relative: Option<f64>,
}

pub fn compute_error(computed: f64, reference: f64) -> Option<ErrorMeasurement> {
    if !computed.is_finite() || !reference.is_finite() {
        return None;
    }
    let absolute = (computed - reference).abs();
    let relative = if reference.abs() > 1e-15 {
        Some(absolute / reference.abs())
    } else {
        None
    };
    Some(ErrorMeasurement { absolute, relative })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorStats {
    pub count: usize,
    pub abs_max: f64,
    pub abs_mean: f64,
    pub abs_p50: f64,
    pub abs_p95: f64,
    pub abs_p99: f64,
    pub rel_max: f64,
    pub rel_mean: f64,
    pub rel_p50: f64,
    pub rel_p95: f64,
    pub rel_p99: f64,
}

impl ErrorStats {
    pub fn from_errors(errors: &[ErrorMeasurement]) -> Self {
        if errors.is_empty() {
            return Self::empty();
        }

        let mut abs_vals: Vec<f64> = errors.iter().map(|e| e.absolute).collect();
        abs_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mut rel_vals: Vec<f64> = errors.iter().filter_map(|e| e.relative).collect();
        rel_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let abs_max = *abs_vals.last().unwrap_or(&0.0);
        let abs_mean = mean(&abs_vals);
        let abs_p50 = percentile(&abs_vals, 0.50);
        let abs_p95 = percentile(&abs_vals, 0.95);
        let abs_p99 = percentile(&abs_vals, 0.99);

        let (rel_max, rel_mean, rel_p50, rel_p95, rel_p99) = if rel_vals.is_empty() {
            (0.0, 0.0, 0.0, 0.0, 0.0)
        } else {
            (
                *rel_vals.last().unwrap_or(&0.0),
                mean(&rel_vals),
                percentile(&rel_vals, 0.50),
                percentile(&rel_vals, 0.95),
                percentile(&rel_vals, 0.99),
            )
        };

        Self {
            count: abs_vals.len(),
            abs_max, abs_mean, abs_p50, abs_p95, abs_p99,
            rel_max, rel_mean, rel_p50, rel_p95, rel_p99,
        }
    }

    pub fn empty() -> Self {
        Self {
            count: 0,
            abs_max: 0.0, abs_mean: 0.0, abs_p50: 0.0, abs_p95: 0.0, abs_p99: 0.0,
            rel_max: 0.0, rel_mean: 0.0, rel_p50: 0.0, rel_p95: 0.0, rel_p99: 0.0,
        }
    }
}

fn mean(vals: &[f64]) -> f64 {
    if vals.is_empty() { 0.0 } else { vals.iter().sum::<f64>() / vals.len() as f64 }
}

fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() { return 0.0; }
    let idx = ((sorted.len() - 1) as f64 * p).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}
