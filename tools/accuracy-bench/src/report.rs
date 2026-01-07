//! Report generation.

use crate::FunctionResult;
use comfy_table::{ContentArrangement, Table};
use serde::{Deserialize, Serialize};
use std::time::{SystemTime, UNIX_EPOCH};

#[derive(Debug, Serialize, Deserialize)]
pub struct Report {
    pub timestamp: u64,
    pub results: Vec<FunctionResult>,
}

impl Report {
    pub fn new(results: Vec<FunctionResult>) -> Self {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0);
        Self { timestamp, results }
    }

    pub fn print_table(&self) {
        println!("\n================================================================================");
        println!("  ACCURACY REPORT");
        println!("================================================================================\n");

        let mut table = Table::new();
        table.set_content_arrangement(ContentArrangement::Dynamic);
        table.set_header(vec![
            "Function",
            "I16F16 rel_mean",
            "I16F16 rel_max",
            "I32F32 rel_mean",
            "I32F32 rel_max",
            "Samples",
        ]);

        for r in &self.results {
            table.add_row(vec![
                r.name.clone(),
                format!("{:.6e}", r.i16f16.rel_mean),
                format!("{:.6e}", r.i16f16.rel_max),
                format!("{:.6e}", r.i32f32.rel_mean),
                format!("{:.6e}", r.i32f32.rel_max),
                r.samples_tested.to_string(),
            ]);
        }
        println!("{table}\n");
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap_or_else(|_| "{}".to_string())
    }
}
