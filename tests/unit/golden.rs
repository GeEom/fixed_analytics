//! Pins every function's output bits in tests/golden.txt; CHANGELOG.md
//! explains how a change is re-pinned and recorded.

#[cfg(test)]
mod tests {
    use crate::unit::support::Lcg;
    use fixed::types::extra::{U5, U61};
    use fixed::types::{I16F16, I24F8, I32F32, I64F64};
    use fixed::{FixedI8, FixedI64};
    use fixed_analytics::{
        CordicNumber, acos, acosh, acoth, asin, asin_from_gap, asinh, atan, atan_range, atan2,
        atanh, cos, cosh, coth, div_or_saturate, dot, exp, exp_range, hypot, ln, log2, log10,
        mul_add, mul_div, norm, normalize, pow, pow2, rsqrt, scale_by, scale_exponent, sin,
        sin_cos, sin_cos_range, sinh, sinh_cosh, sqrt, sqrt_or_zero, sqrt_product, tan, tanh,
        tanh_range,
    };
    use std::collections::BTreeMap;
    use std::fmt::Write;

    const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden.txt");
    const CHANGELOG: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/CHANGELOG.md");
    const INPUTS: usize = 600;

    const FUNCTIONS: &[&str] = &[
        "sin",
        "cos",
        "tan",
        "sin_cos",
        "asin",
        "acos",
        "atan",
        "atan2",
        "sinh",
        "cosh",
        "tanh",
        "sinh_cosh",
        "coth",
        "asinh",
        "acosh",
        "atanh",
        "acoth",
        "exp",
        "pow2",
        "ln",
        "log2",
        "log10",
        "pow",
        "sqrt",
        "sqrt_or_zero",
        "mul_add",
        "mul_div",
        "div_or_saturate",
        "dot",
        "hypot",
        "norm",
        "normalize",
        "scale_exponent",
        "scale_by",
        "sqrt_product",
        "rsqrt",
        "asin_from_gap",
        "sin_cos_range",
        "exp_range",
        "atan_range",
        "tanh_range",
    ];

    trait Raw: CordicNumber {
        const NAME: &'static str;
        fn raw(self) -> i128;
        fn from_raw(raw: i128) -> Self;
    }

    macro_rules! raw {
        ($($t:ty, $bits:ty, $name:literal);+ $(;)?) => {$(
            impl Raw for $t {
                const NAME: &'static str = $name;
                fn raw(self) -> i128 {
                    i128::from(self.to_bits())
                }
                fn from_raw(raw: i128) -> Self {
                    let clamped = raw.clamp(<$bits>::MIN.into(), <$bits>::MAX.into());
                    Self::from_bits(<$bits>::try_from(clamped).unwrap_or(0))
                }
            }
        )+};
    }

    raw! {
        FixedI8<U5>, i8, "I3F5";
        I16F16, i32, "I16F16";
        I24F8, i32, "I24F8";
        I32F32, i64, "I32F32";
        FixedI64<U61>, i64, "I3F61";
        I64F64, i128, "I64F64";
    }

    const fn pair<T: Copy>(range: (T, T)) -> [Option<T>; 2] {
        [Some(range.0), Some(range.1)]
    }

    fn eval<T: CordicNumber>(name: &str, x: T, y: T) -> [Option<T>; 2] {
        let one = |v: T| [Some(v), None];
        let res = |v: fixed_analytics::Result<T>| [v.ok(), None];
        match name {
            "sin" => one(sin(x)),
            "cos" => one(cos(x)),
            "tan" => one(tan(x)),
            "sin_cos" => {
                let (s, c) = sin_cos(x);
                [Some(s), Some(c)]
            }
            "asin" => res(asin(x)),
            "acos" => res(acos(x)),
            "atan" => one(atan(x)),
            "atan2" => one(atan2(x, y)),
            "sinh" => one(sinh(x)),
            "cosh" => one(cosh(x)),
            "tanh" => one(tanh(x)),
            "sinh_cosh" => {
                let (s, c) = sinh_cosh(x);
                [Some(s), Some(c)]
            }
            "coth" => res(coth(x)),
            "asinh" => one(asinh(x)),
            "acosh" => res(acosh(x)),
            "atanh" => res(atanh(x)),
            "acoth" => res(acoth(x)),
            "exp" => one(exp(x)),
            "pow2" => one(pow2(x)),
            "ln" => res(ln(x)),
            "log2" => res(log2(x)),
            "log10" => res(log10(x)),
            "pow" => res(pow(x, y)),
            "sqrt" => res(sqrt(x)),
            "sqrt_or_zero" => one(sqrt_or_zero(x)),
            "mul_add" => one(mul_add(x, y, x)),
            "mul_div" => one(mul_div(x, y, y)),
            "div_or_saturate" => one(div_or_saturate(x, y)),
            "dot" => one(dot(&[x, y, x], &[y, y, x])),
            "hypot" => one(hypot(x, y)),
            "norm" => one(norm(&[x, y, x])),
            "normalize" => normalize([x, y]).map_or([None, None], |v| [Some(v[0]), Some(v[1])]),
            "scale_exponent" => {
                let e = scale_exponent(x);
                [
                    e.map(|v| T::from_num(v.rem_euclid(4))),
                    e.map(|v| T::from_num(v.div_euclid(16).clamp(-3, 3))),
                ]
            }
            "scale_by" => one(scale_by(x, y.to_i32() % 40)),
            "sqrt_product" => res(sqrt_product(x, y)),
            "rsqrt" => res(rsqrt(x)),
            "asin_from_gap" => res(asin_from_gap(x)),
            "sin_cos_range" => {
                let ((s_lo, s_hi), (c_lo, c_hi)) = sin_cos_range(x, y);
                [Some(s_lo + s_hi), Some(c_lo + c_hi)]
            }
            "exp_range" => pair(exp_range(x, y)),
            "atan_range" => pair(atan_range(x, y)),
            "tanh_range" => pair(tanh_range(x, y)),
            _ => [None, None],
        }
    }

    /// Fixed values, then raw bits whose length is uniform, so tiny and huge
    /// values come up as often as ordinary ones.
    fn inputs<T: Raw>(seed: u64) -> Vec<T> {
        let mut v: Vec<T> = [0.0, 1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 3.0, -3.0]
            .iter()
            .map(|&x| T::from_num(x))
            .collect();
        v.extend([
            T::min_value(),
            T::max_value(),
            T::from_raw(1),
            T::from_raw(-1),
        ]);
        let mut rng = Lcg(seed);
        let bits = T::total_bits();
        while v.len() < INPUTS {
            let len = u32::try_from(rng.next_u64() % u64::from(bits)).unwrap_or(0);
            let r = i128::from(rng.next_u64()) << 64 | i128::from(rng.next_u64());
            let mag = r.unsigned_abs().checked_shr(128 - len).unwrap_or(0);
            let raw = i128::try_from(mag).unwrap_or(0);
            v.push(T::from_raw(if rng.next_u64() & 1 == 0 {
                raw
            } else {
                -raw
            }));
        }
        v
    }

    fn fnv(hash: &mut u64, bytes: &[u8]) {
        for &b in bytes {
            *hash = (*hash ^ u64::from(b)).wrapping_mul(0x0100_0000_01B3);
        }
    }

    fn digest<T: Raw>(name: &str) -> u64 {
        let seed = name.bytes().fold(0x9E37_79B9_7F4A_7C15, |h: u64, b| {
            h.rotate_left(7) ^ u64::from(b)
        });
        let xs = inputs::<T>(seed);
        let ys = inputs::<T>(seed ^ 0x5555);
        let mut hash = 0xCBF2_9CE4_8422_2325;
        for (&x, &y) in xs.iter().zip(&ys) {
            for out in eval(name, x, y) {
                match out {
                    Some(v) => fnv(&mut hash, &v.raw().to_le_bytes()),
                    None => fnv(&mut hash, b"none"),
                }
            }
        }
        hash
    }

    fn current() -> BTreeMap<(String, String), u64> {
        let mut map = BTreeMap::new();
        for &name in FUNCTIONS {
            for (ty, hash) in [
                (<FixedI8<U5> as Raw>::NAME, digest::<FixedI8<U5>>(name)),
                (<I16F16 as Raw>::NAME, digest::<I16F16>(name)),
                (<I24F8 as Raw>::NAME, digest::<I24F8>(name)),
                (<I32F32 as Raw>::NAME, digest::<I32F32>(name)),
                (<FixedI64<U61> as Raw>::NAME, digest::<FixedI64<U61>>(name)),
                (<I64F64 as Raw>::NAME, digest::<I64F64>(name)),
            ] {
                map.insert((name.to_string(), ty.to_string()), hash);
            }
        }
        map
    }

    /// Each row: function, type, hash and the version that last changed it.
    type Pinned = BTreeMap<(String, String), (u64, String)>;

    fn pinned() -> Result<Pinned, String> {
        let text = std::fs::read_to_string(GOLDEN).map_err(|e| format!("{GOLDEN}: {e}"))?;
        let mut map = BTreeMap::new();
        for line in text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
        {
            let cols: Vec<&str> = line.split_whitespace().collect();
            let [name, ty, hash, changed] = cols.as_slice() else {
                return Err(format!("malformed golden row: {line}"));
            };
            let hash = u64::from_str_radix(hash, 16).map_err(|e| format!("{line}: {e}"))?;
            map.insert(
                ((*name).to_string(), (*ty).to_string()),
                (hash, (*changed).to_string()),
            );
        }
        Ok(map)
    }

    #[test]
    fn golden_outputs() -> Result<(), String> {
        let now = current();
        let old = pinned().unwrap_or_default();
        let changed: Vec<String> = now
            .iter()
            .filter(|(k, h)| old.get(*k).is_some_and(|(o, _)| o != *h))
            .map(|((name, ty), _)| format!("{name} on {ty}"))
            .collect();
        let missing = now.keys().any(|k| !old.contains_key(k));
        if std::env::var_os("GOLDEN_UPDATE").is_some() {
            let version = env!("CARGO_PKG_VERSION");
            let mut out = String::from("# function type hash changed_in\n");
            for ((name, ty), hash) in &now {
                let stamp = match old.get(&(name.clone(), ty.clone())) {
                    Some((o, v)) if o == hash => v.clone(),
                    Some(_) => version.to_string(),
                    None => "-".to_string(),
                };
                writeln!(out, "{name} {ty} {hash:016x} {stamp}").map_err(|e| e.to_string())?;
            }
            return std::fs::write(GOLDEN, out).map_err(|e| e.to_string());
        }
        if changed.is_empty() && !missing {
            return Ok(());
        }
        Err(format!(
            "outputs changed for {}{}; rerun with GOLDEN_UPDATE=1 and name each changed function in CHANGELOG.md",
            changed.join(", "),
            if missing { " (and new rows)" } else { "" }
        ))
    }

    #[test]
    fn golden_changes_in_changelog() -> Result<(), String> {
        let log = std::fs::read_to_string(CHANGELOG).map_err(|e| format!("{CHANGELOG}: {e}"))?;
        let mut sections: BTreeMap<String, String> = BTreeMap::new();
        let mut current: Option<String> = None;
        for line in log.lines() {
            if let Some(head) = line.strip_prefix("## ") {
                let version = head.trim_matches(|c: char| c == '[' || c == ']' || c == ' ');
                let version = version.split_whitespace().next().unwrap_or("");
                current = Some(version.trim_end_matches(']').to_string());
            } else if let Some(v) = &current {
                let body = sections.entry(v.clone()).or_default();
                body.push_str(line);
                body.push('\n');
            }
        }
        let mut missing = Vec::new();
        for ((name, ty), (_, version)) in pinned()? {
            if version == "-" {
                continue;
            }
            let named = sections
                .get(&version)
                .is_some_and(|body| body.contains(&format!("`{name}`")));
            if !named {
                missing.push(format!("{name} ({ty}) in {version}"));
            }
        }
        if missing.is_empty() {
            Ok(())
        } else {
            Err(format!(
                "CHANGELOG.md doesn't name these output changes: {}",
                missing.join(", ")
            ))
        }
    }
}
