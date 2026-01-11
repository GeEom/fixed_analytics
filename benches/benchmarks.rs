//! Benchmarks for CORDIC functions.

#![allow(missing_docs, reason = "benchmark code does not need documentation")]

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use fixed::types::I16F16;
use fixed_analytics::{
    acos, acosh, acoth, asin, asinh, atan, atan2, atanh, cos, cosh, coth, exp, ln, log2, log10,
    sin, sin_cos, sinh, sinh_cosh, sqrt, tan, tanh,
};

fn bench_circular(c: &mut Criterion) {
    let angle = I16F16::from_num(0.5);
    let x = I16F16::from_num(0.5);

    c.bench_function("sin", |b| b.iter(|| sin(black_box(angle))));
    c.bench_function("cos", |b| b.iter(|| cos(black_box(angle))));
    c.bench_function("tan", |b| b.iter(|| tan(black_box(angle))));
    c.bench_function("sin_cos", |b| b.iter(|| sin_cos(black_box(angle))));
    c.bench_function("asin", |b| b.iter(|| asin(black_box(x))));
    c.bench_function("acos", |b| b.iter(|| acos(black_box(x))));
    c.bench_function("atan", |b| b.iter(|| atan(black_box(x))));
    c.bench_function("atan2", |b| {
        b.iter(|| atan2(black_box(x), black_box(I16F16::ONE)));
    });
}

fn bench_hyperbolic(c: &mut Criterion) {
    let x = I16F16::from_num(0.5);
    let large_x = I16F16::from_num(1.5);

    c.bench_function("sinh", |b| b.iter(|| sinh(black_box(x))));
    c.bench_function("cosh", |b| b.iter(|| cosh(black_box(x))));
    c.bench_function("tanh", |b| b.iter(|| tanh(black_box(x))));
    c.bench_function("coth", |b| b.iter(|| coth(black_box(x))));
    c.bench_function("sinh_cosh", |b| b.iter(|| sinh_cosh(black_box(x))));
    c.bench_function("asinh", |b| b.iter(|| asinh(black_box(x))));
    c.bench_function("acosh", |b| b.iter(|| acosh(black_box(large_x))));
    c.bench_function("atanh", |b| b.iter(|| atanh(black_box(x))));
    c.bench_function("acoth", |b| b.iter(|| acoth(black_box(large_x))));
}

fn bench_exponential(c: &mut Criterion) {
    let x = I16F16::from_num(0.5);
    let pos_x = I16F16::from_num(2.0);

    c.bench_function("exp", |b| b.iter(|| exp(black_box(x))));
    c.bench_function("ln", |b| b.iter(|| ln(black_box(pos_x))));
    c.bench_function("log2", |b| b.iter(|| log2(black_box(pos_x))));
    c.bench_function("log10", |b| b.iter(|| log10(black_box(pos_x))));
}

fn bench_algebraic(c: &mut Criterion) {
    let x = I16F16::from_num(2.0);

    c.bench_function("sqrt", |b| b.iter(|| sqrt(black_box(x))));
}

criterion_group!(
    benches,
    bench_circular,
    bench_hyperbolic,
    bench_exponential,
    bench_algebraic
);
criterion_main!(benches);
