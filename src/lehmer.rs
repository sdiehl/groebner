//! Lehmer's Euclidean algorithm for the gcds and half gcds of rational reconstruction.
//! Runs of quotients come from the leading words, so the big numbers see one 2x2
//! cofactor matrix per run instead of one long division per quotient.

use num_bigint::{BigInt, BigUint, Sign};
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

// The 63 bits of `a` starting at bit `shift`.
fn window(a: &BigUint, shift: u64) -> i128 {
    let mut digits = a.iter_u64_digits().skip((shift / 64) as usize);
    let lo = u128::from(digits.next().unwrap_or(0));
    let hi = u128::from(digits.next().unwrap_or(0));
    ((hi << 64 | lo) >> (shift % 64)) as i128 & i128::from(i64::MAX)
}

// Knuth's Algorithm L: the cofactors `[a, b, c, d]` of the quotients that the leading
// words of `u >= v` determine, mapping `(u, v)` to `(a u + b v, c u + d v)`. Stops before
// `|d|` would exceed `limit`, and returns `None` when no quotient is determined.
fn run(u: &BigUint, v: &BigUint, limit: i128) -> Option<[i64; 4]> {
    let shift = u.bits().saturating_sub(63);
    let (mut x, mut y) = (window(u, shift), window(v, shift));
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

/// Greatest common divisor of two nonnegative integers.
pub(crate) fn gcd(a: &BigUint, b: &BigUint) -> BigUint {
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
    let Some(mut y) = v.iter_u64_digits().next() else {
        return u;
    };
    let mut x = (u % y).iter_u64_digits().next().unwrap_or(0);
    while x != 0 {
        (x, y) = (y % x, x);
    }
    y.into()
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
