//! Coefficient tables as fixed point, built at compile time from the 384-bit
//! constants, for the working widths' Horner evaluation.

use crate::consts::{C, INV_FACT, INV_ODD};

/// A table of coefficients in `[0, 1]` as fixed point with 62, 126 and 190
/// fractional bits, for the one-, two- and three-limb working widths.
pub struct Series {
    pub q62: &'static [i64],
    pub q126: &'static [i128],
    pub q190: &'static [[u64; 3]],
}

/// `c · 2^126` rounded, for `|c| ≤ 1`.
const fn q126(entry: &C) -> i128 {
    let (neg, exp, mant) = *entry;
    let top = ((mant[5] as u128) << 64) | mant[4] as u128;
    let shift = 2 - exp;
    let value = if shift <= 0 || shift >= 128 {
        0
    } else {
        ((top >> shift) + ((top >> (shift - 1)) & 1)).cast_signed()
    };
    if neg { -value } else { value }
}

const fn q62(entry: &C) -> i64 {
    let (neg, exp, mant) = *entry;
    let shift = 2 - exp;
    let value = if shift <= 0 || shift >= 64 {
        0
    } else {
        ((mant[5] >> shift) + ((mant[5] >> (shift - 1)) & 1)).cast_signed()
    };
    if neg { -value } else { value }
}

/// `|c| · 2^190` rounded, as three limbs, least significant first.
const fn q190(entry: &C) -> [u64; 3] {
    let (_, exp, mant) = *entry;
    let shift = 2 - exp;
    if shift < 0 || shift >= 190 {
        return [0; 3];
    }
    let shift = shift.cast_unsigned();
    let wide = [mant[2], mant[3], mant[4], mant[5]];
    let mut out = [0u64; 3];
    let mut rest: &mut [u64] = &mut out;
    let mut index = 0;
    while let [slot, tail @ ..] = rest {
        let bit = 64 + 64 * index + shift;
        *slot = word(&wide, bit);
        rest = tail;
        index += 1;
    }
    if (word(&wide, 63 + shift) & 1) == 1 {
        let (low, c0) = out[0].overflowing_add(1);
        let (mid, c1) = out[1].overflowing_add(c0 as u64);
        out = [low, mid, out[2] + c1 as u64];
    }
    out
}

/// The 64 bits of `w` from bit `bit` up, reading zero past its end.
const fn word(w: &[u64; 4], bit: u32) -> u64 {
    let lo = limb(w, bit / 64);
    let hi = limb(w, bit / 64 + 1);
    let shift = bit % 64;
    if shift == 0 {
        lo
    } else {
        (lo >> shift) | (hi << (64 - shift))
    }
}

const fn limb(w: &[u64; 4], i: u32) -> u64 {
    match i {
        0 => w[0],
        1 => w[1],
        2 => w[2],
        3 => w[3],
        _ => 0,
    }
}

const fn tables<const L: usize>(table: &[C; L]) -> ([i64; L], [i128; L], [[u64; 3]; L]) {
    let mut q62s = [0; L];
    let mut q126s = [0; L];
    let mut q190s = [[0; 3]; L];
    let mut rest: &[C] = table;
    let mut xs: &mut [i64] = &mut q62s;
    let mut ys: &mut [i128] = &mut q126s;
    let mut zs: &mut [[u64; 3]] = &mut q190s;
    while let ([head, tail @ ..], [q62_slot, xt @ ..], [q126_slot, yt @ ..], [q190_slot, zt @ ..]) =
        (rest, xs, ys, zs)
    {
        *q62_slot = q62(head);
        *q126_slot = q126(head);
        *q190_slot = q190(head);
        rest = tail;
        xs = xt;
        ys = yt;
        zs = zt;
    }
    (q62s, q126s, q190s)
}

pub const FACT: Series = Series {
    q62: &tables(&INV_FACT).0,
    q126: &tables(&INV_FACT).1,
    q190: &tables(&INV_FACT).2,
};

pub const ODD: Series = Series {
    q62: &tables(&INV_ODD).0,
    q126: &tables(&INV_ODD).1,
    q190: &tables(&INV_ODD).2,
};
