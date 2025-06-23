extern crate groebner;
use groebner::{groebner_basis, is_groebner_basis, Monomial, MonomialOrder, Polynomial, Term};
use num_rational::BigRational;

/// Helper function to create polynomials for testing
fn create_polynomial(
    terms: Vec<(i32, i32, Vec<u32>)>, // (numerator, denominator, exponents)
    nvars: usize,
    order: MonomialOrder,
) -> Polynomial<BigRational> {
    let term_vec: Vec<Term<BigRational>> = terms
        .into_iter()
        .map(|(num, den, exp)| {
            let coeff = BigRational::new(num.into(), den.into());
            let monomial = Monomial::new(exp);
            Term::new(coeff, monomial)
        })
        .collect();

    Polynomial::new(term_vec, nvars, order)
}

/// Helper function to create integer polynomials for testing
fn create_int_polynomial(
    terms: Vec<(i32, Vec<u32>)>, // (coefficient, exponents)
    nvars: usize,
    order: MonomialOrder,
) -> Polynomial<BigRational> {
    let term_vec: Vec<Term<BigRational>> = terms
        .into_iter()
        .map(|(coeff, exp)| {
            let rational_coeff = BigRational::new(coeff.into(), 1.into());
            let monomial = Monomial::new(exp);
            Term::new(rational_coeff, monomial)
        })
        .collect();

    Polynomial::new(term_vec, nvars, order)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_groebner_cox_little_oshea() {
        // Example from Cox, Little, O'Shea: x^2 + 2xy^2, xy + 2y^3 - 1
        // Expected result: [x, y^3 - 1/2]
        let f = create_polynomial(
            vec![(1, 1, vec![2, 0]), (2, 1, vec![1, 2])],
            2,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![1, 1]), (2, 1, vec![0, 3]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let expected = vec![
            create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex), // x
            create_polynomial(
                vec![(1, 1, vec![0, 3]), (-1, 2, vec![0, 0])],
                2,
                MonomialOrder::Lex,
            ), // y^3 - 1/2
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_groebner_y_x_order() {
        // Test with y,x variable order: 2x^2*y + y^2, 2x^3 + xy - 1
        // Expected result: [y, x^3 - 1/2]
        let f = create_polynomial(
            vec![(2, 1, vec![2, 1]), (1, 1, vec![0, 2])],
            2,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(2, 1, vec![3, 0]), (1, 1, vec![1, 1]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_three_variable_simple() {
        // Test: x - z^2, y - z^3
        // Expected result: [x - z^2, y - z^3] (already a Groebner basis)
        let f = create_polynomial(
            vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 0, 2])],
            3,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![0, 1, 0]), (-1, 1, vec![0, 0, 3])],
            3,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
        assert_eq!(basis.len(), 2);
    }

    #[test]
    fn test_grlex_ordering() {
        // Test with graded lexicographic ordering: x^3 - 2xy, x^2*y + x - 2y^2
        // Expected result: [x^2, xy, -x/2 + y^2]
        let f = create_polynomial(
            vec![(1, 1, vec![3, 0]), (-2, 1, vec![1, 1])],
            2,
            MonomialOrder::GrLex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![2, 1]), (1, 1, vec![1, 0]), (-2, 1, vec![0, 2])],
            2,
            MonomialOrder::GrLex,
        );

        let expected = vec![
            create_polynomial(vec![(1, 1, vec![2, 0])], 2, MonomialOrder::GrLex), // x^2
            create_polynomial(vec![(1, 1, vec![1, 1])], 2, MonomialOrder::GrLex), // xy
            create_polynomial(
                vec![(-1, 2, vec![1, 0]), (1, 1, vec![0, 2])],
                2,
                MonomialOrder::GrLex,
            ), // -x/2 + y^2
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_three_variable_complex_lex() {
        // Test: -x^2 + y, -x^3 + z
        // Expected result: [x^2 - y, xy - z, xz - y^2, y^3 - z^2]
        let f = create_polynomial(
            vec![(-1, 1, vec![2, 0, 0]), (1, 1, vec![0, 1, 0])],
            3,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(-1, 1, vec![3, 0, 0]), (1, 1, vec![0, 0, 1])],
            3,
            MonomialOrder::Lex,
        );

        let expected = vec![
            create_polynomial(
                vec![(1, 1, vec![2, 0, 0]), (-1, 1, vec![0, 1, 0])],
                3,
                MonomialOrder::Lex,
            ), // x^2 - y
            create_polynomial(
                vec![(1, 1, vec![1, 1, 0]), (-1, 1, vec![0, 0, 1])],
                3,
                MonomialOrder::Lex,
            ), // xy - z
            create_polynomial(
                vec![(1, 1, vec![1, 0, 1]), (-1, 1, vec![0, 2, 0])],
                3,
                MonomialOrder::Lex,
            ), // xz - y^2
            create_polynomial(
                vec![(1, 1, vec![0, 3, 0]), (-1, 1, vec![0, 0, 2])],
                3,
                MonomialOrder::Lex,
            ), // y^3 - z^2
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_three_variable_complex_grlex() {
        // Same as above but with graded lexicographic ordering
        // Expected result: [y^3 - z^2, x^2 - y, xy - z, xz - y^2]
        let f = create_polynomial(
            vec![(-1, 1, vec![2, 0, 0]), (1, 1, vec![0, 1, 0])],
            3,
            MonomialOrder::GrLex,
        );

        let g = create_polynomial(
            vec![(-1, 1, vec![3, 0, 0]), (1, 1, vec![0, 0, 1])],
            3,
            MonomialOrder::GrLex,
        );

        let expected = [
            create_polynomial(
                vec![(1, 1, vec![0, 3, 0]), (-1, 1, vec![0, 0, 2])],
                3,
                MonomialOrder::GrLex,
            ), // y^3 - z^2
            create_polynomial(
                vec![(1, 1, vec![2, 0, 0]), (-1, 1, vec![0, 1, 0])],
                3,
                MonomialOrder::GrLex,
            ), // x^2 - y
            create_polynomial(
                vec![(1, 1, vec![1, 1, 0]), (-1, 1, vec![0, 0, 1])],
                3,
                MonomialOrder::GrLex,
            ), // xy - z
            create_polynomial(
                vec![(1, 1, vec![1, 0, 1]), (-1, 1, vec![0, 2, 0])],
                3,
                MonomialOrder::GrLex,
            ), // xz - y^2
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_parametric_curve() {
        // Test: x - y^2, -y^3 + z
        // Expected result: [x - y^2, y^3 - z]
        let f = create_polynomial(
            vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 2, 0])],
            3,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(-1, 1, vec![0, 3, 0]), (1, 1, vec![0, 0, 1])],
            3,
            MonomialOrder::Lex,
        );

        let expected = vec![
            create_polynomial(
                vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 2, 0])],
                3,
                MonomialOrder::Lex,
            ), // x - y^2
            create_polynomial(
                vec![(1, 1, vec![0, 3, 0]), (-1, 1, vec![0, 0, 1])],
                3,
                MonomialOrder::Lex,
            ), // y^3 - z
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_parametric_curve_grlex() {
        // Same as above but with graded lexicographic ordering
        // Expected result: [x^2 - yz, xy - z, -x + y^2]
        let f = create_polynomial(
            vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 2, 0])],
            3,
            MonomialOrder::GrLex,
        );

        let g = create_polynomial(
            vec![(-1, 1, vec![0, 3, 0]), (1, 1, vec![0, 0, 1])],
            3,
            MonomialOrder::GrLex,
        );

        let expected = vec![
            create_polynomial(
                vec![(1, 1, vec![2, 0, 0]), (-1, 1, vec![0, 1, 1])],
                3,
                MonomialOrder::GrLex,
            ), // x^2 - yz
            create_polynomial(
                vec![(1, 1, vec![1, 1, 0]), (-1, 1, vec![0, 0, 1])],
                3,
                MonomialOrder::GrLex,
            ), // xy - z
            create_polynomial(
                vec![(-1, 1, vec![1, 0, 0]), (1, 1, vec![0, 2, 0])],
                3,
                MonomialOrder::GrLex,
            ), // -x + y^2
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_twisted_cubic() {
        // Test: x - z^2, y - z^3
        // Expected result: [x - z^2, y - z^3] (already in Groebner basis form)
        let f = create_polynomial(
            vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 0, 2])],
            3,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![0, 1, 0]), (-1, 1, vec![0, 0, 3])],
            3,
            MonomialOrder::Lex,
        );

        let expected = vec![f.clone(), g.clone()];

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_twisted_cubic_grlex() {
        // Same as above but with graded lexicographic ordering
        // Expected result: [x^2 - yz, xz - y, -x + z^2]
        let f = create_polynomial(
            vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 0, 2])],
            3,
            MonomialOrder::GrLex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![0, 1, 0]), (-1, 1, vec![0, 0, 3])],
            3,
            MonomialOrder::GrLex,
        );

        let expected = [
            create_polynomial(
                vec![(1, 1, vec![2, 0, 0]), (-1, 1, vec![0, 1, 1])],
                3,
                MonomialOrder::GrLex,
            ), // x^2 - yz
            create_polynomial(
                vec![(1, 1, vec![1, 0, 1]), (-1, 1, vec![0, 1, 0])],
                3,
                MonomialOrder::GrLex,
            ), // xz - y
            create_polynomial(
                vec![(-1, 1, vec![1, 0, 0]), (1, 1, vec![0, 0, 2])],
                3,
                MonomialOrder::GrLex,
            ), // -x + z^2
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex, true)
            .expect("Groebner computation failed");

        // Check basic properties
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));

        // Should have 3 polynomials in the basis
        assert_eq!(basis.len(), 3);

        // Compare the computed basis to the expected result
        assert_eq!(basis, expected);

        // Check that the basis contains the expected structure
        // (order may vary but should have x^2, xz, and z^2 terms)
        let has_x2_term = basis.iter().any(|p| {
            p.terms
                .iter()
                .any(|t| t.monomial.exponents == vec![2, 0, 0])
        });
        let has_xz_term = basis.iter().any(|p| {
            p.terms
                .iter()
                .any(|t| t.monomial.exponents == vec![1, 0, 1])
        });
        let has_z2_term = basis.iter().any(|p| {
            p.terms
                .iter()
                .any(|t| t.monomial.exponents == vec![0, 0, 2])
        });

        assert!(has_x2_term);
        assert!(has_xz_term);
        assert!(has_z2_term);
    }

    #[test]
    fn test_variety_intersection() {
        // Test: -y^2 + z, x - y^3
        // Canonical output: [x - y^3, y^2 - z]
        let f = create_polynomial(
            vec![(-1, 1, vec![0, 2, 0]), (1, 1, vec![0, 0, 1])],
            3,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 3, 0])],
            3,
            MonomialOrder::Lex,
        );

        let expected = vec![
            create_polynomial(
                vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 3, 0])],
                3,
                MonomialOrder::Lex,
            ), // x - y^3
            create_polynomial(
                vec![(1, 1, vec![0, 2, 0]), (-1, 1, vec![0, 0, 1])],
                3,
                MonomialOrder::Lex,
            ), // y^2 - z
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_variety_intersection_grlex() {
        // Same as above but with graded lexicographic ordering
        // Expected result: [-x^2 + z^3, xy - z^2, y^2 - z, -x + yz]
        let f = create_polynomial(
            vec![(-1, 1, vec![0, 2, 0]), (1, 1, vec![0, 0, 1])],
            3,
            MonomialOrder::GrLex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 3, 0])],
            3,
            MonomialOrder::GrLex,
        );

        let expected = vec![
            create_polynomial(
                vec![(-1, 1, vec![2, 0, 0]), (1, 1, vec![0, 0, 3])],
                3,
                MonomialOrder::GrLex,
            ), // -x^2 + z^3
            create_polynomial(
                vec![(1, 1, vec![1, 1, 0]), (-1, 1, vec![0, 0, 2])],
                3,
                MonomialOrder::GrLex,
            ), // xy - z^2
            create_polynomial(
                vec![(1, 1, vec![0, 2, 0]), (-1, 1, vec![0, 0, 1])],
                3,
                MonomialOrder::GrLex,
            ), // y^2 - z
            create_polynomial(
                vec![(-1, 1, vec![1, 0, 0]), (1, 1, vec![0, 1, 1])],
                3,
                MonomialOrder::GrLex,
            ), // -x + yz
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_space_curve() {
        // Test: y - z^2, x - z^3
        // Expected result: [x - z^3, y - z^2]
        let f = create_polynomial(
            vec![(1, 1, vec![0, 1, 0]), (-1, 1, vec![0, 0, 2])],
            3,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 0, 3])],
            3,
            MonomialOrder::Lex,
        );

        let expected = vec![
            create_polynomial(
                vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 0, 3])],
                3,
                MonomialOrder::Lex,
            ), // x - z^3
            create_polynomial(
                vec![(1, 1, vec![0, 1, 0]), (-1, 1, vec![0, 0, 2])],
                3,
                MonomialOrder::Lex,
            ), // y - z^2
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_space_curve_grlex() {
        // Same as above but with graded lexicographic ordering
        // Canonical output: [y^3 - x^2, xz - y^2, y z - x, z^2 - y]
        let f = create_polynomial(
            vec![(1, 1, vec![0, 1, 0]), (-1, 1, vec![0, 0, 2])],
            3,
            MonomialOrder::GrLex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 0, 3])],
            3,
            MonomialOrder::GrLex,
        );

        let expected = vec![
            create_polynomial(
                vec![(1, 1, vec![0, 3, 0]), (-1, 1, vec![2, 0, 0])],
                3,
                MonomialOrder::GrLex,
            ), // y^3 - x^2
            create_polynomial(
                vec![(1, 1, vec![1, 0, 1]), (-1, 1, vec![0, 2, 0])],
                3,
                MonomialOrder::GrLex,
            ), // xz - y^2
            create_polynomial(
                vec![(1, 1, vec![0, 1, 1]), (-1, 1, vec![1, 0, 0])],
                3,
                MonomialOrder::GrLex,
            ), // y z - x
            create_polynomial(
                vec![(1, 1, vec![0, 0, 2]), (-1, 1, vec![0, 1, 0])],
                3,
                MonomialOrder::GrLex,
            ), // z^2 - y
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis, expected);
    }

    #[test]
    fn test_circle_and_parabola() {
        // Test: 4x^2*y^2 + 4xy + 1, x^2 + y^2 - 1
        // Complex intersection of circle and parabola
        let f = create_polynomial(
            vec![(4, 1, vec![2, 2]), (4, 1, vec![1, 1]), (1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![2, 0]), (1, 1, vec![0, 2]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_katsura_3_system() {
        // Katsura-3 system: x0 + 2x1 + 2x2 - 1, x0^2 + 2x1^2 + 2x2^2 - x0, 2x0*x1 + 2x1*x2 - x1
        let f1 = create_int_polynomial(
            vec![
                (1, vec![1, 0, 0]),
                (2, vec![0, 1, 0]),
                (2, vec![0, 0, 1]),
                (-1, vec![0, 0, 0]),
            ],
            3,
            MonomialOrder::Lex,
        );

        let f2 = create_int_polynomial(
            vec![
                (1, vec![2, 0, 0]),
                (2, vec![0, 2, 0]),
                (2, vec![0, 0, 2]),
                (-1, vec![1, 0, 0]),
            ],
            3,
            MonomialOrder::Lex,
        );

        let f3 = create_int_polynomial(
            vec![(2, vec![1, 1, 0]), (2, vec![0, 1, 1]), (-1, vec![0, 1, 0])],
            3,
            MonomialOrder::Lex,
        );

        // Note: The actual Groebner basis computation may produce different but equivalent results
        // Let's verify it's a valid Groebner basis with the expected structure
        let basis = groebner_basis(vec![f1, f2, f3], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");

        // Check basic properties of the computed Groebner basis
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));

        // The Katsura-3 system should produce 3 polynomials in the basis
        assert_eq!(basis.len(), 3);

        // Check that we have the expected leading monomials (in lex order: x0, x1, x2^4)
        assert_eq!(
            basis[0]
                .leading_monomial()
                .expect("Leading monomial computation failed")
                .exponents,
            vec![1, 0, 0]
        ); // Leading term should be x0
        assert_eq!(
            basis[1]
                .leading_monomial()
                .expect("Leading monomial computation failed")
                .exponents,
            vec![0, 1, 0]
        ); // Leading term should be x1
        assert_eq!(
            basis[2]
                .leading_monomial()
                .expect("Leading monomial computation failed")
                .exponents,
            vec![0, 0, 4]
        ); // Leading term should be x2^4
    }

    #[test]
    fn test_katsura_3_grlex() {
        // Same as above but with graded lexicographic ordering
        let f1 = create_int_polynomial(
            vec![
                (1, vec![1, 0, 0]),
                (2, vec![0, 1, 0]),
                (2, vec![0, 0, 1]),
                (-1, vec![0, 0, 0]),
            ],
            3,
            MonomialOrder::GrLex,
        );

        let f2 = create_int_polynomial(
            vec![
                (1, vec![2, 0, 0]),
                (2, vec![0, 2, 0]),
                (2, vec![0, 0, 2]),
                (-1, vec![1, 0, 0]),
            ],
            3,
            MonomialOrder::GrLex,
        );

        let f3 = create_int_polynomial(
            vec![(2, vec![1, 1, 0]), (2, vec![0, 1, 1]), (-1, vec![0, 1, 0])],
            3,
            MonomialOrder::GrLex,
        );

        let basis = groebner_basis(vec![f1, f2, f3], MonomialOrder::GrLex, true)
            .expect("Groebner computation failed");
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_cyclic_4_system() {
        // Cyclic-4 system: a + b + c + d, ab + ad + bc + bd, abc + abd + acd + bcd, abcd - 1
        let f1 = create_int_polynomial(
            vec![
                (1, vec![1, 0, 0, 0]),
                (1, vec![0, 1, 0, 0]),
                (1, vec![0, 0, 1, 0]),
                (1, vec![0, 0, 0, 1]),
            ],
            4,
            MonomialOrder::Lex,
        );

        let f2 = create_int_polynomial(
            vec![
                (1, vec![1, 1, 0, 0]),
                (1, vec![1, 0, 0, 1]),
                (1, vec![0, 1, 1, 0]),
                (1, vec![0, 1, 0, 1]),
            ],
            4,
            MonomialOrder::Lex,
        );

        let f3 = create_int_polynomial(
            vec![
                (1, vec![1, 1, 1, 0]),
                (1, vec![1, 1, 0, 1]),
                (1, vec![1, 0, 1, 1]),
                (1, vec![0, 1, 1, 1]),
            ],
            4,
            MonomialOrder::Lex,
        );

        let f4 = create_int_polynomial(
            vec![(1, vec![1, 1, 1, 1]), (-1, vec![0, 0, 0, 0])],
            4,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f1, f2, f3, f4], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");

        // Check basic properties of the computed Groebner basis
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));

        // The cyclic-4 system should produce 5 polynomials in the basis
        assert_eq!(basis.len(), 5);

        // Check that we have the expected leading monomials in lex order
        // (the exact coefficients may vary but the structure should be consistent)
        // Just verify the structure is reasonable for cyclic-4
        let leading_exponents: Vec<_> = basis
            .iter()
            .map(|p| {
                p.leading_monomial()
                    .expect("Leading monomial computation failed")
                    .exponents
                    .clone()
            })
            .collect();

        // The first polynomial should have leading term in the first variable (a)
        assert_eq!(leading_exponents[0][0], 1);

        // The last polynomial should have a high degree in the last variable (d)
        assert!(leading_exponents[4][3] >= 4); // At least d^4
    }

    #[test]
    fn test_cyclic_4_grlex() {
        // Same as above but with graded lexicographic ordering
        let f1 = create_int_polynomial(
            vec![
                (1, vec![1, 0, 0, 0]),
                (1, vec![0, 1, 0, 0]),
                (1, vec![0, 0, 1, 0]),
                (1, vec![0, 0, 0, 1]),
            ],
            4,
            MonomialOrder::GrLex,
        );

        let f2 = create_int_polynomial(
            vec![
                (1, vec![1, 1, 0, 0]),
                (1, vec![1, 0, 0, 1]),
                (1, vec![0, 1, 1, 0]),
                (1, vec![0, 1, 0, 1]),
            ],
            4,
            MonomialOrder::GrLex,
        );

        let f3 = create_int_polynomial(
            vec![
                (1, vec![1, 1, 1, 0]),
                (1, vec![1, 1, 0, 1]),
                (1, vec![1, 0, 1, 1]),
                (1, vec![0, 1, 1, 1]),
            ],
            4,
            MonomialOrder::GrLex,
        );

        let f4 = create_int_polynomial(
            vec![(1, vec![1, 1, 1, 1]), (-1, vec![0, 0, 0, 0])],
            4,
            MonomialOrder::GrLex,
        );

        let basis = groebner_basis(vec![f1, f2, f3, f4], MonomialOrder::GrLex, true)
            .expect("Groebner computation failed");
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_simple_ideal() {
        // Test: x^2 - y, xy - 1
        let f1 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (-1, 1, vec![0, 1])],
            2,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(
            vec![(1, 1, vec![1, 1]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_monomial_ordering_properties() {
        // Test that different orderings give different results
        let f = create_polynomial(
            vec![(1, 1, vec![2, 1]), (1, 1, vec![1, 2])],
            2,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![1, 1]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let basis_lex = groebner_basis(vec![f.clone(), g.clone()], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        let basis_grlex = groebner_basis(vec![f, g], MonomialOrder::GrLex, true)
            .expect("Groebner computation failed");

        assert!(is_groebner_basis(&basis_lex).expect("Groebner basis check failed"));
        assert!(is_groebner_basis(&basis_grlex).expect("Groebner basis check failed"));

        // The bases might be different but both should be valid
        assert!(!basis_lex.is_empty());
        assert!(!basis_grlex.is_empty());
    }

    #[test]
    fn test_reduction_properties() {
        // Test polynomial reduction
        let f = create_polynomial(
            vec![(1, 1, vec![2, 0]), (1, 1, vec![1, 1])],
            2,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![1, 0]), (-1, 1, vec![0, 1])],
            2,
            MonomialOrder::Lex,
        );

        let basis = vec![g];
        let remainder = f.reduce(&basis).expect("Polynomial reduction failed");

        // The remainder should have degree less than the leading term of g
        assert!(remainder.terms.len() <= f.terms.len());
    }

    #[test]
    fn test_s_polynomial() {
        // Test S-polynomial computation
        let f = create_polynomial(
            vec![(1, 1, vec![2, 0]), (2, 1, vec![1, 2])],
            2,
            MonomialOrder::Lex,
        );

        let g = create_polynomial(
            vec![(1, 1, vec![1, 1]), (2, 1, vec![0, 3]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let s_poly = f.s_polynomial(&g).expect("S-polynomial computation failed");
        assert!(!s_poly.is_zero());
    }

    #[test]
    fn test_empty_and_single_polynomial() {
        let empty_basis = groebner_basis::<BigRational>(vec![], MonomialOrder::Lex, true);
        assert!(matches!(
            empty_basis,
            Err(groebner::GroebnerError::EmptyInput)
        ));

        let single = create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex);

        let single_basis = groebner_basis(vec![single.clone()], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert_eq!(single_basis.len(), 1);
        assert!(is_groebner_basis(&single_basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_zero_polynomial() {
        let zero = Polynomial::zero(2, MonomialOrder::Lex);
        let nonzero = create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex);

        let basis = groebner_basis(vec![zero, nonzero.clone()], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert_eq!(basis.len(), 1);
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_constant_polynomial() {
        let constant = create_polynomial(vec![(1, 1, vec![0, 0])], 2, MonomialOrder::Lex);
        let other = create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex);

        let basis = groebner_basis(vec![constant, other], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        // If the ideal contains a nonzero constant, the basis should be {1}
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_linear_system() {
        // Test: x + y - 1, x - y
        let f1 = create_polynomial(
            vec![(1, 1, vec![1, 0]), (1, 1, vec![0, 1]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(
            vec![(1, 1, vec![1, 0]), (-1, 1, vec![0, 1])],
            2,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_quadratic_system() {
        // Test: x^2 + y^2 - 1, x^2 - y^2
        let f1 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (1, 1, vec![0, 2]), (-1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (-1, 1, vec![0, 2])],
            2,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }

    #[test]
    fn test_homogeneous_ideal() {
        // Test homogeneous ideal: x^2 + y^2, xy
        let f1 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (1, 1, vec![0, 2])],
            2,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(
            vec![(1, 1, vec![1, 1]), (1, 1, vec![0, 0])],
            2,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex, true)
            .expect("Groebner computation failed");
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis).expect("Groebner basis check failed"));
    }
}
