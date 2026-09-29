//! Lehmer's Euclidean algorithm for the gcds and half gcds of rational reconstruction.
//! Runs of quotients come from the leading words, so the big numbers see one 2x2
//! cofactor matrix per run instead of one long division per quotient.

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

// The 63 bits of the little-endian `digits` starting at bit `shift`.
fn window(digits: impl Iterator<Item = u64>, shift: u64) -> i128 {
    let mut digits = digits.skip((shift / 64) as usize);
    let lo = u128::from(digits.next().unwrap_or(0));
    let hi = u128::from(digits.next().unwrap_or(0));
    ((hi << 64 | lo) >> (shift % 64)) as i128 & i128::from(i64::MAX)
}

fn run(u: &BigUint, v: &BigUint, limit: i128) -> Option<[i64; 4]> {
    let shift = u.bits().saturating_sub(63);
    let (x, y) = (
        window(u.iter_u64_digits(), shift),
        window(v.iter_u64_digits(), shift),
    );
    cofactors(x, y, limit)
}

// Knuth's Algorithm L: the cofactors `[a, b, c, d]` of the quotients that the leading
// words `x >= y` of `u >= v` determine, mapping `(u, v)` to `(a u + b v, c u + d v)`.
// Stops before `|d|` would exceed `limit`, and returns `None` when no quotient is
// determined.
fn cofactors(mut x: i128, mut y: i128, limit: i128) -> Option<[i64; 4]> {
    let (mut a, mut b, mut c, mut d) = (1i128, 0i128, 0i128, 1i128);
    // Every operand is below 2^64, so the quotients take the hardware 64-bit divide.
    let quotient = |n: i128, m: i128| (n >= 0 && m > 0).then(|| i128::from(n as u64 / m as u64));
    while let (Some(q), Some(r)) = (quotient(x + a, y + c), quotient(x + b, y + d)) {
        if q != r || (b - q * d).abs() > limit {
            break;
        }
        (a, c) = (c, a - q * c);
        (b, d) = (d, b - q * d);
        (x, y) = (y, x - q * y);
    }
    (b != 0).then(|| [a, b, c, d].map(|t| t as i64))
}

fn apply(m: [i64; 4], u: &BigInt, v: &BigInt) -> (BigInt, BigInt) {
    (u * m[0] + v * m[1], u * m[2] + v * m[3])
}

fn nonnegative(x: BigInt) -> BigUint {
    x.into_parts().1
}

// Operands of at most this many words take the allocation-free gcd.
const WORDS: usize = 8;
type Words = [u64; WORDS];

fn bits(x: &Words) -> u64 {
    x.iter()
        .rposition(|&w| w != 0)
        .map_or(0, |i| 64 * i as u64 + u64::from(64 - x[i].leading_zeros()))
}

// `(a u + b v, c u + d v)` for Lehmer cofactors, whose results are nonnegative. The
// signs of `a` and `b` differ, so each limb sum stays within an `i128`.
fn apply_words(m: [i64; 4], u: &Words, v: &Words) -> (Words, Words) {
    let combine = |a: i64, b: i64| {
        let (mut out, mut carry) = ([0; WORDS], 0i128);
        for i in 0..WORDS {
            let s = i128::from(a) * i128::from(u[i]) + i128::from(b) * i128::from(v[i]) + carry;
            out[i] = s as u64;
            carry = s >> 64;
        }
        out
    };
    (combine(m[0], m[1]), combine(m[2], m[3]))
}

fn words(x: &BigUint) -> Option<Words> {
    let mut out = [0; WORDS];
    for (i, d) in x.iter_u64_digits().enumerate() {
        *out.get_mut(i)? = d;
    }
    Some(out)
}

// Lehmer's gcd of `u >= v` on fixed words, finishing in a single word.
fn gcd_words(mut u: Words, mut v: Words) -> BigUint {
    while bits(&v) > 64 {
        let shift = bits(&u).saturating_sub(63);
        let m = cofactors(
            window(u.iter().copied(), shift),
            window(v.iter().copied(), shift),
            i128::from(i64::MAX),
        );
        if let Some(m) = m {
            (u, v) = apply_words(m, &u, &v);
        } else {
            // A quotient too large for a word, which is rare after the first step.
            let (a, b) = (
                BigUint::from_slice(&to_u32(&u)),
                BigUint::from_slice(&to_u32(&v)),
            );
            let r = a % &b;
            (u, v) = (v, words(&r).unwrap_or_else(|| unreachable!("r < v")));
        }
    }
    let y = v[0];
    if y == 0 {
        return BigUint::from_slice(&to_u32(&u));
    }
    let top = bits(&u).div_ceil(64) as usize;
    let x = u[..top].iter().rev().fold(0, |r, &w| {
        ((u128::from(r) << 64 | u128::from(w)) % u128::from(y)) as u64
    });
    euclid(x, y).into()
}

fn to_u32(x: &Words) -> Vec<u32> {
    x.iter()
        .flat_map(|&w| [w as u32, (w >> 32) as u32])
        .collect()
}

fn euclid(mut x: u64, mut y: u64) -> u64 {
    while x != 0 {
        (x, y) = (y % x, x);
    }
    y
}

// `gcd(u, y)` for a word `y <= u`.
fn finish(u: BigUint, y: u64) -> BigUint {
    if y == 0 {
        return u;
    }
    euclid((u % y).iter_u64_digits().next().unwrap_or(0), y).into()
}

/// Greatest common divisor of two nonnegative integers.
pub(crate) fn gcd(a: &BigUint, b: &BigUint) -> BigUint {
    if let (Some(x), Some(y)) = (words(a), words(b)) {
        return if a >= b {
            gcd_words(x, y)
        } else {
            gcd_words(y, x)
        };
    }
    let (mut u, mut v) = if a >= b {
        (a.clone(), b.clone())
    } else {
        (b.clone(), a.clone())
    };
    while v.bits() > 64 {
        if let Some(m) = run(&u, &v, i128::from(i64::MAX)) {
            let (s, t) = apply(m, &u.into(), &v.into());
            (u, v) = (nonnegative(s), nonnegative(t));
        } else {
            (u, v) = (v.clone(), u % v);
        }
    }
    finish(u, v.iter_u64_digits().next().unwrap_or(0))
}

/// `n / d` in lowest terms, for `d > 0`.
pub(crate) fn fraction(n: BigInt, d: &BigUint) -> BigRational {
    let g = gcd(n.magnitude(), d);
    let (sign, n) = n.into_parts();
    BigRational::new_raw(BigInt::from_biguint(sign, n / &g), BigInt::from(d / g))
}

/// Wang's rational reconstruction of `x mod m`: the fraction `r / s` with `|r|` and
/// `0 < s` at most `bound`, found at the first Euclidean remainder of `(m, x)` that is at
/// most `bound`.
pub(crate) fn wang(x: &BigInt, m: &BigInt, bound: &BigInt) -> Option<BigRational> {
    let (mut r0, mut r1) = (m.clone(), x.mod_floor(m));
    let (mut s0, mut s1) = (BigInt::zero(), BigInt::from(1u32));
    let bits = bound.bits();
    while &r1 > bound {
        // `r0 <= 2 |d| r0'`, so this keeps `r0' > bound` and the first remainder at most
        // `bound` cannot be stepped over.
        let limit = r0.bits().saturating_sub(bits + 2).min(62);
        if let Some(m) = (limit > 0)
            .then(|| run(r0.magnitude(), r1.magnitude(), 1 << limit))
            .flatten()
        {
            (r0, r1) = apply(m, &r0, &r1);
            (s0, s1) = apply(m, &s0, &s1);
        } else {
            let (q, r) = r0.div_rem(&r1);
            r0 = std::mem::replace(&mut r1, r);
            let s = &s0 - q * &s1;
            s0 = std::mem::replace(&mut s1, s);
        }
    }
    if s1.is_zero() || s1.magnitude() > bound.magnitude() {
        return None;
    }
    if !gcd(r1.magnitude(), s1.magnitude()).is_one() {
        return None;
    }
    let sign = if s1.sign() == Sign::Minus { -r1 } else { r1 };
    Some(BigRational::new_raw(sign, s1.abs()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use polycore::crt;

    // A deterministic stream of big integers with mixed sizes and structure.
    fn samples() -> impl Iterator<Item = BigUint> {
        let mut s = 0x9e37_79b9_7f4a_7c15u64;
        (0..400).map(move |i| {
            let mut next = || {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                s
            };
            let words = 1 + i % 24;
            let mut x = BigUint::zero();
            for _ in 0..words {
                x = (x << 64u32) + next();
            }
            x >> (next() % 64)
        })
    }

    #[test]
    fn gcd_matches_euclid() {
        let xs: Vec<_> = samples().collect();
        for (a, b) in xs.iter().zip(xs.iter().skip(1)) {
            let common = &xs[a.bits() as usize % xs.len()];
            let (a, b) = (a * common, b * common);
            assert_eq!(gcd(&a, &b), a.gcd(&b));
            assert_eq!(gcd(&a, &BigUint::zero()), a);
        }
    }

    #[test]
    fn wang_matches_plain_reconstruction() {
        let m: BigInt = samples()
            .take(8)
            .fold(BigUint::from(1u32), |m, x| m * (x | BigUint::from(1u32)))
            .into();
        let context = crt::WangContext::new(&m).unwrap();
        let bound = (&m / 2u32).sqrt();
        for (n, d) in samples().zip(samples().skip(7)) {
            let (n, d) = (BigInt::from(n) % &bound, BigInt::from(d) % &bound + 1u32);
            // A reduced fraction and an unrelated residue, which usually fails.
            let Some(inv) = d.modinv(&m) else { continue };
            for x in [(&n * inv).mod_floor(&m), &n * 7 + &d] {
                assert_eq!(wang(&x, &m, &bound), context.reconstruct(&x));
            }
        }
    }
}
