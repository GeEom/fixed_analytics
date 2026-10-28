//! Prints the engine's constants as 384 normalized bits, for src/consts.rs.

use rug::float::Constant;
use rug::{Float, Integer};

const BITS: u32 = 384;

fn hex(v: u64) -> String {
    format!(
        "0x{:04X}_{:04X}_{:04X}_{:04X}",
        v >> 48,
        (v >> 32) & 0xFFFF,
        (v >> 16) & 0xFFFF,
        v & 0xFFFF
    )
}

fn parts(x: &Float) -> String {
    if x.is_zero() {
        return "(false, 0, [0; 6])".to_string();
    }
    let neg = x.is_sign_negative();
    let x = &Float::with_val(x.prec(), x.abs_ref());
    let e = x.get_exp().unwrap_or(0);
    let m: Integer = (Float::with_val(BITS + 64, x) << (BITS as i32 - e))
        .to_integer()
        .unwrap();
    let limbs: Vec<String> = (0..6)
        .map(|i| hex((m.clone() >> (64 * i) as u32).to_u64_wrapping()))
        .collect();
    format!("({neg}, {e}, [{}])", limbs.join(", "))
}

fn table(name: &str, xs: &[Float]) {
    println!("pub const {name}: [C; {}] = [", xs.len());
    for x in xs {
        println!("    {},", parts(x));
    }
    println!("];");
}

fn emit(name: &str, x: Float) {
    let e = x.get_exp().unwrap_or(0);
    let m: Integer = (Float::with_val(BITS + 64, &x) << (BITS as i32 - e))
        .to_integer()
        .unwrap();
    let limbs: Vec<String> = (0..6)
        .map(|i| hex((m.clone() >> (64 * i) as u32).to_u64_wrapping()))
        .collect();
    println!("pub const {name}: C = (false, {e}, [{}]);", limbs.join(", "));
}

fn main() {
    let p = BITS + 128;
    let pi = Float::with_val(p, Constant::Pi);
    let ln2 = Float::with_val(p, Constant::Log2);
    let ln10 = Float::with_val(p, 10).ln();
    let one = Float::with_val(p, 1);
    emit("PI", pi.clone());
    emit("LN_2", ln2.clone());
    emit("LN_10", ln10.clone());
    emit("FRAC_1_LN_2", one.clone() / &ln2);
    emit("FRAC_1_LN_10", one.clone() / &ln10);
    emit("FRAC_1_2PI", one.clone() / (pi.clone() * 2u32));
    emit("EULER_GAMMA", Float::with_val(p, Constant::Euler));
    emit("FRAC_2_SQRT_PI", Float::with_val(p, 2) / pi.clone().sqrt());
    emit("SQRT_PI", pi.clone().sqrt());
    emit("SQRT_2", Float::with_val(p, 2).sqrt());
    let mut fact = Float::with_val(p, 1);
    let mut inv_fact = vec![];
    for n in 0..49u32 {
        if n > 0 {
            fact *= n;
        }
        inv_fact.push(Float::with_val(p, 1) / &fact);
    }
    table("INV_FACT", &inv_fact);
    let inv_odd: Vec<Float> = (0..49u32).map(|k| Float::with_val(p, 1) / (2 * k + 1)).collect();
    table("INV_ODD", &inv_odd);
    let atan16: Vec<Float> = (0..17u32).map(|j| Float::with_val(p, j) / 16u32).map(|x| x.atan()).collect();
    table("ATAN_SIXTEENTHS", &atan16);
    let ln16: Vec<Float> = (-5..8i32).map(|j| (Float::with_val(p, j) / 16u32 + 1u32).ln()).collect();
    table("LN_SIXTEENTHS", &ln16);
    let angle = |j: u32| pi.clone() * j / 128u32;
    table("SIN_PI_128THS", &(0..33u32).map(|j| angle(j).sin()).collect::<Vec<_>>());
    table("COS_PI_128THS", &(0..33u32).map(|j| angle(j).cos()).collect::<Vec<_>>());
    table("EXP2_32NDS", &(0..32u32).map(|j| (Float::with_val(p, j) / 32u32).exp2()).collect::<Vec<_>>());
}
