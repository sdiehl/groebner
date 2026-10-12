//! C ABI over a ring handle, with polynomials passed as strings in the [`PolynomialRing`]
//! syntax.
//!
//! Every entry point returns a [`GbStatus`]. On failure it writes an error message to
//! `out`; on success the computing functions write their result there instead. Strings
//! written to `out` are freed with [`gb_string_free`] and rings with [`gb_ring_free`].
#![allow(unsafe_code)]

use std::ffi::{CStr, CString, c_char, c_int};
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::slice;

use num_rational::BigRational;

use crate::{
    F4Field, GroebnerError, Ideal, MonomialOrder, ParseCoefficient, ParsePolynomialError,
    Polynomial, PolynomialRing, Zp, groebner_basis_f4,
};

/// Status code returned by every entry point.
#[repr(C)]
#[derive(Clone, Copy)]
pub enum GbStatus {
    Success = 0,
    InvalidArgument = 1,
    ParseError = 2,
    ComputeError = 3,
    Panic = 4,
}

/// Monomial order codes accepted by [`gb_ring_new`].
#[repr(C)]
pub enum GbOrder {
    Lex = 0,
    GrLex = 1,
    GRevLex = 2,
}

/// Opaque polynomial ring over Q or a prime field.
pub enum GbRing {
    Rational(PolynomialRing<BigRational>),
    Prime(PolynomialRing<Zp>),
}

struct Failure(GbStatus, String);

impl From<ParsePolynomialError> for Failure {
    fn from(e: ParsePolynomialError) -> Self {
        Self(GbStatus::ParseError, e.to_string())
    }
}

impl From<GroebnerError> for Failure {
    fn from(e: GroebnerError) -> Self {
        Self(GbStatus::ComputeError, e.to_string())
    }
}

fn invalid(message: impl Into<String>) -> Failure {
    Failure(GbStatus::InvalidArgument, message.into())
}

fn order(code: c_int) -> Result<MonomialOrder, Failure> {
    match code {
        c if c == GbOrder::Lex as c_int => Ok(MonomialOrder::Lex),
        c if c == GbOrder::GrLex as c_int => Ok(MonomialOrder::GrLex),
        c if c == GbOrder::GRevLex as c_int => Ok(MonomialOrder::GRevLex),
        _ => Err(invalid(format!("unknown order {code}"))),
    }
}

fn parse_all<F: ParseCoefficient>(
    ring: &PolynomialRing<F>,
    polys: &[&str],
) -> Result<Vec<Polynomial<F>>, Failure> {
    Ok(polys
        .iter()
        .map(|p| ring.parse(p))
        .collect::<Result<_, _>>()?)
}

fn basis<F: F4Field + ParseCoefficient>(
    ring: &PolynomialRing<F>,
    polys: &[&str],
) -> Result<String, Failure> {
    let basis = groebner_basis_f4(parse_all(ring, polys)?, true)?;
    let lines = basis
        .iter()
        .map(|p| ring.format(p))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(lines.join("\n"))
}

fn reduce<F: F4Field + ParseCoefficient>(
    ring: &PolynomialRing<F>,
    polys: &[&str],
    f: &str,
) -> Result<String, Failure> {
    let ideal = Ideal::new(parse_all(ring, polys)?)?;
    Ok(ring.format(&ideal.normal_form(&ring.parse(f)?))?)
}

unsafe fn str_arg<'a>(s: *const c_char) -> Result<&'a str, Failure> {
    if s.is_null() {
        return Err(invalid("null string"));
    }
    unsafe { CStr::from_ptr(s) }
        .to_str()
        .map_err(|e| invalid(e.to_string()))
}

unsafe fn str_args<'a>(ptr: *const *const c_char, len: usize) -> Result<Vec<&'a str>, Failure> {
    if len == 0 {
        return Ok(Vec::new());
    }
    if ptr.is_null() {
        return Err(invalid("null array"));
    }
    let ptrs = unsafe { slice::from_raw_parts(ptr, len) };
    ptrs.iter().map(|&s| unsafe { str_arg(s) }).collect()
}

unsafe fn ring_arg<'a>(ring: *const GbRing) -> Result<&'a GbRing, Failure> {
    unsafe { ring.as_ref() }.ok_or_else(|| invalid("null ring"))
}

unsafe fn write_string(out: *mut *mut c_char, text: String) {
    if !out.is_null() {
        unsafe { *out = CString::new(text).unwrap_or_default().into_raw() };
    }
}

/// Run `f`, catching panics. Errors are written to `out` and returned as a status.
unsafe fn guard<T>(
    out: *mut *mut c_char,
    f: impl FnOnce() -> Result<T, Failure>,
) -> Result<T, c_int> {
    let failure = match catch_unwind(AssertUnwindSafe(f)) {
        Ok(Ok(value)) => return Ok(value),
        Ok(Err(failure)) => failure,
        Err(_) => Failure(GbStatus::Panic, "internal panic".into()),
    };
    unsafe { write_string(out, failure.1) };
    Err(failure.0 as c_int)
}

unsafe fn respond(out: *mut *mut c_char, f: impl FnOnce() -> Result<String, Failure>) -> c_int {
    match unsafe { guard(out, f) } {
        Ok(text) => {
            unsafe { write_string(out, text) };
            GbStatus::Success as c_int
        }
        Err(status) => status,
    }
}

/// Build a ring in `nvars` variables, over Q when `modulus` is 0 and GF(`modulus`) otherwise.
///
/// # Safety
/// `vars` must point to `nvars` valid null terminated UTF-8 strings. `ring` must be writable
/// and `out` must be null or writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gb_ring_new(
    vars: *const *const c_char,
    nvars: usize,
    order_code: c_int,
    modulus: u64,
    ring: *mut *mut GbRing,
    out: *mut *mut c_char,
) -> c_int {
    let built = unsafe {
        guard(out, || {
            if ring.is_null() {
                return Err(invalid("null ring output"));
            }
            let vars = str_args(vars, nvars)?;
            let order = order(order_code)?;
            Ok(if modulus == 0 {
                GbRing::Rational(PolynomialRing::new(vars, order)?)
            } else {
                GbRing::Prime(PolynomialRing::with_modulus(vars, order, modulus)?)
            })
        })
    };
    match built {
        Ok(built) => {
            unsafe { *ring = Box::into_raw(Box::new(built)) };
            GbStatus::Success as c_int
        }
        Err(status) => status,
    }
}

/// Free a ring built by [`gb_ring_new`].
///
/// # Safety
/// `ring` must be null or a handle from [`gb_ring_new`] not already freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gb_ring_free(ring: *mut GbRing) {
    if !ring.is_null() {
        drop(unsafe { Box::from_raw(ring) });
    }
}

/// Reduced Groebner basis of `npolys` polynomials, written to `out` one per line.
///
/// # Safety
/// `ring` must be a live handle, `polys` must point to `npolys` valid null terminated UTF-8
/// strings, and `out` must be null or writable.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gb_groebner_basis(
    ring: *const GbRing,
    polys: *const *const c_char,
    npolys: usize,
    out: *mut *mut c_char,
) -> c_int {
    unsafe {
        respond(out, || {
            let polys = str_args(polys, npolys)?;
            match ring_arg(ring)? {
                GbRing::Rational(r) => basis(r, &polys),
                GbRing::Prime(r) => basis(r, &polys),
            }
        })
    }
}

/// Normal form of `f` modulo the ideal generated by `npolys` polynomials, written to `out`.
///
/// # Safety
/// As [`gb_groebner_basis`], and `f` must be a valid null terminated UTF-8 string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gb_reduce(
    ring: *const GbRing,
    polys: *const *const c_char,
    npolys: usize,
    f: *const c_char,
    out: *mut *mut c_char,
) -> c_int {
    unsafe {
        respond(out, || {
            let polys = str_args(polys, npolys)?;
            let f = str_arg(f)?;
            match ring_arg(ring)? {
                GbRing::Rational(r) => reduce(r, &polys, f),
                GbRing::Prime(r) => reduce(r, &polys, f),
            }
        })
    }
}

/// Free a string written to `out` by this library.
///
/// # Safety
/// `s` must be null or a pointer written to `out` by this library and not already freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn gb_string_free(s: *mut c_char) {
    if !s.is_null() {
        drop(unsafe { CString::from_raw(s) });
    }
}
