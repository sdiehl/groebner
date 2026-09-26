# Groebner Basis

This is an implementation of the F4 and Buchberger algorithm for computing Groebner bases in Rust.

Examples:

- [Buchberger example](examples/buchberger.rs)
- [F4 example](examples/f4.rs)

## Usage

```bash
cargo add groebner
```

Parse polynomials with named variables under a monomial order and compute a reduced Groebner basis with F4 (over `PrimeField<P>`, `Zp` or `BigRational`):

```rust
use groebner::{groebner_basis_f4, MonomialOrder, PolynomialRing, PrimeField};

let ring = PolynomialRing::<PrimeField<32003>>::new(["x", "y"], MonomialOrder::GRevLex)?;
let polys = ring.parse_many("x^2 - y; x*y - 1")?;
let basis = groebner_basis_f4(polys, true)?;
for p in &basis {
    println!("{}", ring.format(p)?);
}
```

Buchberger over any `Field` is available as `groebner_basis` (and `groebner_basis_parallel` with the default `parallel` feature). Runtime primes up to 64 bits use `Zp`:

```rust
use groebner::{groebner_basis_f4, MonomialOrder, PolynomialRing, Zp};

let ring = PolynomialRing::<Zp>::with_modulus(["x", "y"], MonomialOrder::Lex, 1_000_003)?;
let basis = groebner_basis_f4(ring.parse_many("x^2 - y; x*y - 1")?, true)?;
```

`Ideal` wraps a reduced basis and answers the usual questions: membership and normal forms, elimination ideals, zero-dimensionality, the standard monomial basis and its dimension, radical membership, multiplication matrices, and change of order (FGLM for zero-dimensional ideals):

```rust
use groebner::{Ideal, MonomialOrder, PolynomialRing};
use num_rational::BigRational;

let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::GRevLex)?;
let ideal = Ideal::new(ring.parse_many("x^2 + y^2 - 1; x - y^3")?)?;
assert_eq!(ideal.vector_space_dimension(), Some(6));
let lex = ideal.change_order(MonomialOrder::Lex)?;
assert!(lex.contains(&ring.parse("y^6 + y^2 - 1")?)?);
```

`Ideal::lift` returns cofactors `h` with `f = sum h[k] * generators[k]` for members (as in Singular's `lift`), so an external checker can confirm membership with `verify_lift`, which uses only ring addition and multiplication:

```rust
use groebner::{Ideal, MonomialOrder, PolynomialRing, verify_lift};
use num_rational::BigRational;

let ring = PolynomialRing::<BigRational>::new(["x", "y"], MonomialOrder::GRevLex)?;
let ideal = Ideal::new(ring.parse_many("x^2 - y; y^2 - x")?)?;
let f = ring.parse("x^4 - x")?;
let h = ideal.lift(&f)?.ok_or("not a member")?;
assert!(verify_lift(ideal.generators(), &h, &f));
```

Monomial orders: `Lex`, `GrLex`, `GRevLex`, `MonomialOrder::weighted(weights, tie_break)`, and product orders via `MonomialOrder::block` or `MonomialOrder::elimination(k, rest)`.

Systems with one symbolic parameter work over `RationalFunction`, the field Q(a). The result is the basis for a generic value of the parameter, and `specialize` substitutes a number:

```rust
use groebner::{groebner_basis_f4, MonomialOrder, PolynomialRing, RationalFunction};

let ring = PolynomialRing::<RationalFunction>::with_parameter(["x", "y"], MonomialOrder::Lex, "a")?;
let basis = groebner_basis_f4(ring.parse_many("x^2 - a; x*y - 1")?, true)?;
assert_eq!(ring.format(&basis[1])?, "y^2 - 1/a");
```

## Test Suite

```bash
cargo test
cargo bench
```

The [test suite](SUITE.md) is the full list of known Groebner bases for a variety of large multivariate polynomial systems from several textbooks and some trusted Mathematica generated corpus. Both the Rust algo implementations have to correctly produce the same textbook outputs and Mathematica for all inputs, up to re-ordering.

## References

The classic papers on this topic:

1. Cox, D., Little, J., O'Shea, D. "Ideals, Varieties, and Algorithms"
1. Buchberger, B. "Gröbner Bases: An Algorithmic Method in Polynomial Ideal Theory"
1. Giovini, A., Mora, T., Niesi, G., Robbiano, L., & Traverso, C. (1991, June). “One sugar cube, please” or selection strategies in the Buchberger algorithm. In Proceedings of the 1991 international symposium on Symbolic and algebraic computation (pp. 49-54).
1. Gebauer, R., & Möller, H. M. (1988). On an installation of Buchberger's algorithm. Journal of Symbolic computation, 6(2-3), 275-286.
1. Roune, B. H., & Stillman, M. (2012, July). Practical Gröbner basis computation. In Proceedings of the 37th International Symposium on Symbolic and Algebraic Computation (pp. 203-210).
1. Faugère, J. C. (1999). A new efficient algorithm for computing Gröbner bases (F4). Journal of Pure and Applied Algebra, 139(1-3), 61-88.
1. Faugère, J. C., Gianni, P., Lazard, D., & Mora, T. (1993). Efficient computation of zero-dimensional Gröbner bases by change of ordering. Journal of Symbolic Computation, 16(4), 329-344.
1. Monagan, M., & Pearce, R. (2015). A compact parallel implementation of F4. In Proceedings of PASCO 2015 (pp. 95-100).

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.
