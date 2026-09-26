//! Loader and checker for the public benchmark corpus in `tests/corpus`.
//!
//! Each system is a `<name>.txt` file (header lines `char:` and `vars:`, then one expanded
//! polynomial per line) next to a `<name>.ref` file computed by an independent implementation.
//! The reference holds the reduced GRevLex basis over `prime`, one row per polynomial: its
//! leading exponent vector and an FNV-1a hash of its monic terms sorted by exponent vector.

#![allow(dead_code, unreachable_pub)]

use groebner::{
    F4Field, GroebnerError, MonomialOrder, ParseCoefficient, Polynomial, PolynomialRing, Zp,
    groebner_basis, groebner_basis_f4,
};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, ToPrimitive};
use std::fs;
use std::path::{Path, PathBuf};

pub struct System {
    pub name: String,
    pub path: PathBuf,
    pub char: u64,
    pub vars: Vec<String>,
    pub polys: Vec<String>,
    pub reference: Option<Reference>,
}

pub struct Reference {
    pub prime: u64,
    pub time: f64,
    pub rows: Vec<Row>,
}

pub type Row = (Vec<u32>, u64);

pub fn load_dir(root: &Path) -> Vec<System> {
    let mut paths = Vec::new();
    let Ok(sources) = fs::read_dir(root) else {
        return Vec::new();
    };
    for source in sources.flatten() {
        let Ok(entries) = fs::read_dir(source.path()) else {
            continue;
        };
        for entry in entries.flatten() {
            if entry.path().extension().is_some_and(|e| e == "txt") {
                paths.push(entry.path());
            }
        }
    }
    paths.sort();
    paths.iter().filter_map(|p| load(p).ok()).collect()
}

pub fn load(path: &Path) -> Result<System, String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut char = 0;
    let mut vars = Vec::new();
    let mut polys = Vec::new();
    for line in text.lines().map(str::trim).filter(|l| !l.is_empty()) {
        if line.starts_with('#') {
        } else if let Some(c) = line.strip_prefix("char:") {
            char = c
                .trim()
                .parse()
                .map_err(|_| format!("bad char in {path:?}"))?;
        } else if let Some(v) = line.strip_prefix("vars:") {
            vars = v.split(',').map(|s| s.trim().to_string()).collect();
        } else {
            polys.push(line.to_string());
        }
    }
    let source = path
        .parent()
        .and_then(Path::file_name)
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let stem = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    let reference = load_reference(&path.with_extension("ref"));
    Ok(System {
        name: format!("{source}/{stem}"),
        path: path.to_path_buf(),
        char,
        vars,
        polys,
        reference,
    })
}

fn load_reference(path: &Path) -> Option<Reference> {
    let text = fs::read_to_string(path).ok()?;
    let mut prime = 0;
    let mut time = 0.0;
    let mut rows = Vec::new();
    for line in text.lines() {
        if let Some(p) = line.strip_prefix("prime:") {
            prime = p.trim().parse().ok()?;
        } else if let Some(t) = line.strip_prefix("time:") {
            time = t.trim().parse().ok()?;
        } else if let Some((exps, hash)) = line.split_once('|') {
            let exps = exps
                .split_whitespace()
                .map(str::parse)
                .collect::<Result<_, _>>()
                .ok()?;
            rows.push((exps, u64::from_str_radix(hash.trim(), 16).ok()?));
        }
    }
    Some(Reference { prime, time, rows })
}

const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

fn fnv(mut h: u64, word: u64) -> u64 {
    for byte in word.to_le_bytes() {
        h = (h ^ u64::from(byte)).wrapping_mul(FNV_PRIME);
    }
    h
}

/// Canonical rows for a basis given each coefficient's residue modulo the reference prime.
pub fn fingerprint<F>(
    basis: &[Polynomial<F>],
    residue: impl Fn(&F) -> Option<u64>,
) -> Option<Vec<Row>> {
    let mut rows = Vec::with_capacity(basis.len());
    for poly in basis {
        let lead = poly.terms.first()?.monomial.exponents().to_vec();
        let mut terms: Vec<(&[u32], u64)> = poly
            .terms
            .iter()
            .map(|t| Some((t.monomial.exponents(), residue(&t.coefficient)?)))
            .collect::<Option<_>>()?;
        terms.sort_by(|a, b| b.0.cmp(a.0));
        let mut h = fnv(FNV_OFFSET, terms.len() as u64);
        for (exps, c) in terms {
            h = exps.iter().fold(h, |h, &e| fnv(h, u64::from(e)));
            h = fnv(h, c);
        }
        rows.push((lead, h));
    }
    rows.sort();
    Some(rows)
}

pub fn compare(actual: &[Row], expected: &[Row]) -> Result<(), String> {
    if actual == expected {
        return Ok(());
    }
    let lead = |rows: &[Row]| rows.iter().map(|r| r.0.clone()).collect::<Vec<_>>();
    if lead(actual) != lead(expected) {
        let missing = expected
            .iter()
            .filter(|r| !actual.iter().any(|a| a.0 == r.0))
            .count();
        let extra = actual
            .iter()
            .filter(|r| !expected.iter().any(|e| e.0 == r.0))
            .count();
        return Err(format!(
            "leading monomials differ: {} vs {} polys, {missing} missing, {extra} extra",
            actual.len(),
            expected.len()
        ));
    }
    let bad = actual
        .iter()
        .zip(expected)
        .filter(|(a, e)| a.1 != e.1)
        .count();
    Err(format!(
        "{bad} of {} polynomials have wrong tails",
        actual.len()
    ))
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Field {
    /// Runtime modulus `Zp` over the reference prime.
    Zp,
    /// Const generic `PrimeField<32003>`, when the reference prime is 32003.
    Gf32003,
    /// Rationals, when the system has characteristic zero; checked by its image mod p.
    Qq,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Algorithm {
    F4,
    Buchberger,
}

impl Field {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "zp" => Some(Self::Zp),
            "gf32003" => Some(Self::Gf32003),
            "qq" => Some(Self::Qq),
            _ => None,
        }
    }

    pub fn applies(self, system: &System) -> bool {
        let prime = system.reference.as_ref().map_or(0, |r| r.prime);
        match self {
            Self::Zp => prime != 0,
            Self::Gf32003 => prime == 32003,
            Self::Qq => prime != 0 && system.char == 0,
        }
    }
}

impl Algorithm {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "f4" => Some(Self::F4),
            "buchberger" => Some(Self::Buchberger),
            _ => None,
        }
    }
}

type F32003 = groebner::PrimeField<32003>;

/// Compute the basis of `system` and return its fingerprint modulo the reference prime.
pub fn run(system: &System, field: Field, algorithm: Algorithm) -> Result<Vec<Row>, String> {
    let prime = system.reference.as_ref().map_or(0, |r| r.prime);
    match field {
        Field::Zp => {
            let ring =
                PolynomialRing::<Zp>::with_modulus(&system.vars, MonomialOrder::GRevLex, prime)
                    .map_err(|e| e.to_string())?;
            let basis = compute(&ring, system, algorithm)?;
            fingerprint(&basis, |c| Some(c.value())).ok_or_else(|| "empty polynomial".into())
        }
        Field::Gf32003 => {
            let ring = PolynomialRing::<F32003>::new(&system.vars, MonomialOrder::GRevLex)
                .map_err(|e| e.to_string())?;
            let basis = compute(&ring, system, algorithm)?;
            fingerprint(&basis, |c| Some(u64::from(c.value())))
                .ok_or_else(|| "empty polynomial".into())
        }
        Field::Qq => {
            let ring = PolynomialRing::<BigRational>::new(&system.vars, MonomialOrder::GRevLex)
                .map_err(|e| e.to_string())?;
            let basis = compute(&ring, system, algorithm)?;
            fingerprint(&basis, |c| rational_residue(c, prime))
                .ok_or_else(|| "unlucky prime for rational basis".into())
        }
    }
}

fn compute<F: F4Field + ParseCoefficient>(
    ring: &PolynomialRing<F>,
    system: &System,
    algorithm: Algorithm,
) -> Result<Vec<Polynomial<F>>, String> {
    let polys = system
        .polys
        .iter()
        .map(|p| ring.parse(p).map_err(|e| format!("parse: {e}")))
        .collect::<Result<Vec<_>, _>>()?;
    let polys: Vec<_> = polys.into_iter().filter(|p| !p.is_zero()).collect();
    let result: Result<_, GroebnerError> = match algorithm {
        Algorithm::F4 => groebner_basis_f4(polys, true),
        Algorithm::Buchberger => groebner_basis(polys, true),
    };
    result.map_err(|e| format!("{e:?}"))
}

fn rational_residue(c: &BigRational, p: u64) -> Option<u64> {
    let modp = |n: &BigInt| {
        let r = (n.abs() % p).to_u64()?;
        Some(if n.is_negative() { (p - r) % p } else { r })
    };
    let num = modp(c.numer())?;
    let den = modp(c.denom())?;
    Some(mul_mod(num, inverse_mod(den, p)?, p))
}

fn mul_mod(a: u64, b: u64, p: u64) -> u64 {
    ((u128::from(a) * u128::from(b)) % u128::from(p)) as u64
}

fn inverse_mod(a: u64, p: u64) -> Option<u64> {
    if a == 0 {
        return None;
    }
    let (mut base, mut exp, mut acc) = (a, p - 2, 1);
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base, p);
        }
        base = mul_mod(base, base, p);
        exp >>= 1;
    }
    Some(acc)
}
