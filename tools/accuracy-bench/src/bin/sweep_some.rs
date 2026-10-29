//! The worst-case sweep for the named functions only, on every type:
//! `cargo run --release --bin sweep_some -- [--quick] name…`.

use accuracy_bench::{sut::types, sweep, worst};
use rayon::prelude::*;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let quick = args.iter().any(|a| a == "--quick");
    sweep::QUICK.store(quick, std::sync::atomic::Ordering::Relaxed);
    let funs: Vec<sweep::Fun> = sweep::ALL
        .iter()
        .copied()
        .filter(|f| args.iter().any(|a| a == f.name()))
        .collect();
    let types = types();
    let jobs: Vec<(sweep::Fun, usize)> = funs
        .iter()
        .flat_map(|&f| (0..types.len()).map(move |t| (f, t)))
        .collect();
    let stats: Vec<_> = jobs.par_iter().map(|&(f, t)| (types[t].run)(f)).collect();
    print!("{}", worst::to_text(&stats));
    worst::check_bounds(&stats);
}
