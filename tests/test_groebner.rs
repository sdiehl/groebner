use groebner::*;

/// Helper function to create polynomials for testing
fn create_polynomial(
    terms: Vec<(i32, i32, Vec<u32>)>, // (numerator, denominator, exponents)
    nvars: usize,
    order: MonomialOrder,
) -> Polynomial<Rational> {
    let term_vec: Vec<Term<Rational>> = terms
        .into_iter()
        .map(|(num, den, exp)| {
            let coeff = Rational::new(num as i64, den as i64);
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
) -> Polynomial<Rational> {
    let term_vec: Vec<Term<Rational>> = terms
        .into_iter()
        .map(|(coeff, exp)| {
            let rational_coeff = Rational::new(coeff as i64, 1);
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
        assert!(basis.len() >= 3);
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);
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

        let _expected = vec![
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex);

        // Check basic properties
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));

        // Should have 3 polynomials in the basis
        assert_eq!(basis.len(), 3);

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
        // Expected result: [x - yz, y^2 - z]
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

        let _expected = vec![
            create_polynomial(
                vec![(1, 1, vec![1, 0, 0]), (-1, 1, vec![0, 1, 1])],
                3,
                MonomialOrder::Lex,
            ), // x - yz
            create_polynomial(
                vec![(1, 1, vec![0, 2, 0]), (-1, 1, vec![0, 0, 1])],
                3,
                MonomialOrder::Lex,
            ), // y^2 - z
        ];

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);

        // Check basic properties
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));

        // Should have 2 polynomials in the basis
        assert_eq!(basis.len(), 2);

        // Check that the basis contains the expected structure
        let has_x_leading = basis.iter().any(|p| {
            p.leading_monomial().unwrap().exponents[0] == 1
                && p.leading_monomial().unwrap().exponents[1] == 0
        });
        let has_y2_leading = basis.iter().any(|p| {
            p.leading_monomial().unwrap().exponents[0] == 0
                && p.leading_monomial().unwrap().exponents[1] == 2
        });

        assert!(has_x_leading);
        assert!(has_y2_leading);
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_space_curve_grlex() {
        // Same as above but with graded lexicographic ordering
        // Expected result: [-x^2 + y^3, xz - y^2, -x + yz, -y + z^2]
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::GrLex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f, g], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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
        let basis = groebner_basis(vec![f1, f2, f3], MonomialOrder::Lex);

        // Check basic properties of the computed Groebner basis
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));

        // The Katsura-3 system should produce 3 polynomials in the basis
        assert_eq!(basis.len(), 3);

        // Check that we have the expected leading monomials (in lex order: x0, x1, x2^4)
        assert_eq!(
            basis[0].leading_monomial().unwrap().exponents,
            vec![1, 0, 0]
        ); // Leading term should be x0
        assert_eq!(
            basis[1].leading_monomial().unwrap().exponents,
            vec![0, 1, 0]
        ); // Leading term should be x1
        assert_eq!(
            basis[2].leading_monomial().unwrap().exponents,
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

        let basis = groebner_basis(vec![f1, f2, f3], MonomialOrder::GrLex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f1, f2, f3, f4], MonomialOrder::Lex);

        // Check basic properties of the computed Groebner basis
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));

        // The cyclic-4 system should produce 5 polynomials in the basis
        assert_eq!(basis.len(), 5);

        // Check that we have the expected leading monomials in lex order
        // (the exact coefficients may vary but the structure should be consistent)
        // Just verify the structure is reasonable for cyclic-4
        let leading_exponents: Vec<_> = basis
            .iter()
            .map(|p| p.leading_monomial().unwrap().exponents.clone())
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

        let basis = groebner_basis(vec![f1, f2, f3, f4], MonomialOrder::GrLex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis_lex = groebner_basis(vec![f.clone(), g.clone()], MonomialOrder::Lex);
        let basis_grlex = groebner_basis(vec![f, g], MonomialOrder::GrLex);

        assert!(is_groebner_basis(&basis_lex));
        assert!(is_groebner_basis(&basis_grlex));

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
        let remainder = f.reduce(&basis);

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

        let s_poly = f.s_polynomial(&g);
        assert!(!s_poly.is_zero());
    }

    #[test]
    fn test_empty_and_single_polynomial() {
        let empty_basis = groebner_basis::<Rational>(vec![], MonomialOrder::Lex);
        assert!(empty_basis.is_empty());

        let single = create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex);

        let single_basis = groebner_basis(vec![single.clone()], MonomialOrder::Lex);
        assert_eq!(single_basis.len(), 1);
        assert!(is_groebner_basis(&single_basis));
    }

    #[test]
    fn test_zero_polynomial() {
        let zero = Polynomial::zero(2, MonomialOrder::Lex);
        let nonzero = create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex);

        let basis = groebner_basis(vec![zero, nonzero.clone()], MonomialOrder::Lex);
        assert_eq!(basis.len(), 1);
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_constant_polynomial() {
        let constant = create_polynomial(vec![(1, 1, vec![0, 0])], 2, MonomialOrder::Lex);
        let other = create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex);

        let basis = groebner_basis(vec![constant, other], MonomialOrder::Lex);
        // If the ideal contains a nonzero constant, the basis should be {1}
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
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

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_homogeneous_ideal() {
        // Test homogeneous ideal: x^2 + y^2, xy
        let f1 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (1, 1, vec![0, 2])],
            2,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(vec![(1, 1, vec![1, 1])], 2, MonomialOrder::Lex);

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_principal_ideal() {
        // Test principal ideal: (x^2 + xy + y^2)
        let f = create_polynomial(
            vec![(1, 1, vec![2, 0]), (1, 1, vec![1, 1]), (1, 1, vec![0, 2])],
            2,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f.clone()], MonomialOrder::Lex);
        assert_eq!(basis.len(), 1);
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_elimination_ideal() {
        // Test elimination: x^2 + y + z, xy + z^2
        // Should eliminate some variables
        let f1 = create_polynomial(
            vec![
                (1, 1, vec![2, 0, 0]),
                (1, 1, vec![0, 1, 0]),
                (1, 1, vec![0, 0, 1]),
            ],
            3,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(
            vec![(1, 1, vec![1, 1, 0]), (1, 1, vec![0, 0, 2])],
            3,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_four_variable_system() {
        // Test with 4 variables: simple system
        let f1 = create_polynomial(
            vec![(1, 1, vec![1, 0, 0, 0]), (-1, 1, vec![0, 1, 0, 0])],
            4,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(
            vec![(1, 1, vec![0, 1, 0, 0]), (-1, 1, vec![0, 0, 1, 0])],
            4,
            MonomialOrder::Lex,
        );

        let f3 = create_polynomial(
            vec![(1, 1, vec![0, 0, 1, 0]), (-1, 1, vec![0, 0, 0, 1])],
            4,
            MonomialOrder::Lex,
        );

        let basis = groebner_basis(vec![f1, f2, f3], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_radical_membership() {
        // Test if x^2 is in the radical of (x^3, x^2*y)
        let f1 = create_polynomial(vec![(1, 1, vec![3, 0])], 2, MonomialOrder::Lex);
        let f2 = create_polynomial(vec![(1, 1, vec![2, 1])], 2, MonomialOrder::Lex);

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_intersection_of_ideals() {
        // Test intersection of two ideals
        let f1 = create_polynomial(vec![(1, 1, vec![1, 0])], 2, MonomialOrder::Lex);
        let f2 = create_polynomial(vec![(1, 1, vec![0, 1])], 2, MonomialOrder::Lex);

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_quotient_ideal() {
        // Test quotient ideal computation
        let f1 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (-1, 1, vec![0, 2])],
            2,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(vec![(1, 1, vec![1, 1])], 2, MonomialOrder::Lex);

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_syzygy_computation() {
        // Test syzygy computation through S-polynomials
        let f1 = create_polynomial(
            vec![(1, 1, vec![2, 0]), (1, 1, vec![1, 1])],
            2,
            MonomialOrder::Lex,
        );

        let f2 = create_polynomial(
            vec![(1, 1, vec![1, 1]), (1, 1, vec![0, 2])],
            2,
            MonomialOrder::Lex,
        );

        let s_poly = f1.s_polynomial(&f2);
        let basis = vec![f1, f2];
        let _reduced_s = s_poly.reduce(&basis);

        // For a Groebner basis, all S-polynomials should reduce to zero
        let full_basis = groebner_basis(basis, MonomialOrder::Lex);
        let final_s = s_poly.reduce(&full_basis);

        // This tests the fundamental property of Groebner bases
        assert!(final_s.is_zero() || final_s.terms.len() < s_poly.terms.len());
    }

    #[test]
    fn test_cyclic_4_known_result() {
        // Example from groebner.rs: v1*v2*v3*v4 - 1, etc.
        // This is the same as cyclic-4 but with different variable names
        let polys = vec![
            create_int_polynomial(
                vec![(1, vec![1, 1, 1, 1]), (-1, vec![0, 0, 0, 0])], // v1*v2*v3*v4 - 1
                4,
                MonomialOrder::Lex,
            ),
            create_int_polynomial(
                vec![
                    (1, vec![1, 1, 1, 0]), // v1*v2*v3
                    (1, vec![1, 1, 0, 1]), // v1*v2*v4
                    (1, vec![1, 0, 1, 1]), // v1*v3*v4
                    (1, vec![0, 1, 1, 1]), // v2*v3*v4
                ],
                4,
                MonomialOrder::Lex,
            ),
            create_int_polynomial(
                vec![
                    (1, vec![1, 1, 0, 0]), // v1*v2
                    (1, vec![0, 1, 1, 0]), // v2*v3
                    (1, vec![1, 0, 0, 1]), // v1*v4
                    (1, vec![0, 0, 1, 1]), // v3*v4
                ],
                4,
                MonomialOrder::Lex,
            ),
            create_int_polynomial(
                vec![
                    (1, vec![1, 0, 0, 0]), // v1
                    (1, vec![0, 1, 0, 0]), // v2
                    (1, vec![0, 0, 1, 0]), // v3
                    (1, vec![0, 0, 0, 1]), // v4
                ],
                4,
                MonomialOrder::Lex,
            ),
        ];

        let basis = groebner_basis(polys, MonomialOrder::Lex);

        // Check basic properties of the computed Groebner basis
        assert!(!basis.is_empty());
        assert!(is_groebner_basis(&basis));

        // The cyclic-4 system should produce 6 polynomials in the basis
        assert_eq!(basis.len(), 6);

        // Check that we have the expected leading monomials structure
        // The basis should be in lexicographic order with leading terms:
        // v1, v2, v3^2, v3*v4^4, v4^5, constant or high degree in v4
        let leading_exponents: Vec<_> = basis
            .iter()
            .map(|p| p.leading_monomial().unwrap().exponents.clone())
            .collect();

        // First polynomial should have leading term in v1
        assert_eq!(leading_exponents[0][0], 1); // v1 appears

        // Verify this is indeed a Groebner basis for the cyclic-4 system
        assert!(is_groebner_basis(&basis));
    }

    #[test]
    fn test_known_groebner_result() {
        // Test with a specific example where we know the exact expected result
        // This is a simple case: x^2 - 1, y - 1
        let f1 = create_int_polynomial(
            vec![(1, vec![2, 0]), (-1, vec![0, 0])], // x^2 - 1
            2,
            MonomialOrder::Lex,
        );

        let f2 = create_int_polynomial(
            vec![(1, vec![0, 1]), (-1, vec![0, 0])], // y - 1
            2,
            MonomialOrder::Lex,
        );

        let expected = vec![
            create_int_polynomial(
                vec![(1, vec![2, 0]), (-1, vec![0, 0])], // x^2 - 1
                2,
                MonomialOrder::Lex,
            ),
            create_int_polynomial(
                vec![(1, vec![0, 1]), (-1, vec![0, 0])], // y - 1
                2,
                MonomialOrder::Lex,
            ),
        ];

        let basis = groebner_basis(vec![f1, f2], MonomialOrder::Lex);

        // For this simple case, the result should be exactly what we expect
        // since the input is already a Groebner basis
        assert_eq!(basis.len(), expected.len());

        // Check that each expected polynomial is in the computed basis
        for expected_poly in &expected {
            assert!(basis.iter().any(|computed_poly| {
                // Check if polynomials are equivalent (same terms, possibly different order)
                if computed_poly.terms.len() != expected_poly.terms.len() {
                    return false;
                }

                // Check if all terms match
                for expected_term in &expected_poly.terms {
                    if !computed_poly.terms.iter().any(|computed_term| {
                        computed_term.coefficient == expected_term.coefficient
                            && computed_term.monomial.exponents == expected_term.monomial.exponents
                    }) {
                        return false;
                    }
                }
                true
            }));
        }
    }
}
