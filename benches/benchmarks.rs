//! Benchmarks for CORDIC functions, on `I16F16` and on `I64F64`.

#![allow(missing_docs, reason = "benchmark code does not need documentation")]

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use fixed::types::{I16F16, I64F64};
use fixed_analytics::{
    CordicNumber, acos, acosh, acoth, asin, asinh, atan, atan2, atanh, cos, cosh, coth, exp, ln,
    log2, log10, pow2, sin, sin_cos, sinh, sinh_cosh, sqrt, tan, tanh,
};

fn bench_type<T: CordicNumber>(c: &mut Criterion, name: &str) {
    let angle = T::from_num(0.5);
    let large_angle = T::from_num(1000.0);
    let x = T::from_num(0.5);
    let large_x = T::from_num(1.5);
    let pos_x = T::from_num(2.0);
    let big = T::from_num(50.0);
    let small = T::from_num(0.005);

    {
        let mut g = c.benchmark_group(format!("{name}/circular"));
        g.bench_function("sin", |b| b.iter(|| sin(black_box(angle))));
        g.bench_function("cos", |b| b.iter(|| cos(black_box(angle))));
        g.bench_function("tan", |b| b.iter(|| tan(black_box(angle))));
        g.bench_function("sin_cos", |b| b.iter(|| sin_cos(black_box(angle))));
        g.bench_function("sin_cos_large", |b| {
            b.iter(|| sin_cos(black_box(large_angle)));
        });
        g.bench_function("asin", |b| b.iter(|| asin(black_box(x))));
        g.bench_function("acos", |b| b.iter(|| acos(black_box(x))));
        g.bench_function("atan", |b| b.iter(|| atan(black_box(x))));
        g.bench_function("atan2", |b| {
            b.iter(|| atan2(black_box(x), black_box(T::one())));
        });
        g.finish();
    }
    {
        let mut g = c.benchmark_group(format!("{name}/hyperbolic"));
        g.bench_function("sinh", |b| b.iter(|| sinh(black_box(x))));
        g.bench_function("cosh", |b| b.iter(|| cosh(black_box(x))));
        g.bench_function("tanh", |b| b.iter(|| tanh(black_box(x))));
        g.bench_function("coth", |b| b.iter(|| coth(black_box(x))));
        g.bench_function("sinh_cosh", |b| b.iter(|| sinh_cosh(black_box(x))));
        g.bench_function("asinh", |b| b.iter(|| asinh(black_box(x))));
        g.bench_function("asinh_large", |b| b.iter(|| asinh(black_box(big))));
        g.bench_function("acosh", |b| b.iter(|| acosh(black_box(large_x))));
        g.bench_function("atanh", |b| b.iter(|| atanh(black_box(x))));
        g.bench_function("acoth", |b| b.iter(|| acoth(black_box(large_x))));
        g.finish();
    }
    {
        let mut g = c.benchmark_group(format!("{name}/exponential"));
        g.bench_function("exp", |b| b.iter(|| exp(black_box(x))));
        g.bench_function("pow2", |b| b.iter(|| pow2(black_box(x))));
        g.bench_function("ln", |b| b.iter(|| ln(black_box(pos_x))));
        g.bench_function("log2", |b| b.iter(|| log2(black_box(pos_x))));
        g.bench_function("log10", |b| b.iter(|| log10(black_box(pos_x))));
        g.finish();
    }
    {
        let mut g = c.benchmark_group(format!("{name}/algebraic"));
        g.bench_function("sqrt", |b| b.iter(|| sqrt(black_box(pos_x))));
        g.bench_function("sqrt_large", |b| b.iter(|| sqrt(black_box(big))));
        g.bench_function("sqrt_small", |b| b.iter(|| sqrt(black_box(small))));
        g.finish();
    }
}

fn bench_i16f16(c: &mut Criterion) {
    bench_type::<I16F16>(c, "I16F16");
}

fn bench_i64f64(c: &mut Criterion) {
    bench_type::<I64F64>(c, "I64F64");
}

criterion_group!(benches, bench_i16f16, bench_i64f64);
criterion_main!(benches);
