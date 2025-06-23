# Groebner Basis Test Suite

| Test Name | System (Input Polynomials) | Canonical Groebner Basis (Output) |
|-----------|---------------------------|-----------------------------------|
| Basic Cox-Little-O'Shea | x² + 2xy², xy + 2y³ - 1 | x, y³ - 1/2 |
| Groebner y,x Order | 2x²y + y², 2x³ + xy - 1 | y, x³ - 1/2 |
| Three Variable Simple | x - z², y - z³ | x - z², y - z³ |
| Grlex Ordering | x³ - 2xy, x²y + x - 2y² | x², xy, -x/2 + y² |
| Three Variable Complex Lex | -x² + y, -x³ + z | x² - y, xy - z, xz - y², y³ - z² |
| Three Variable Complex Grlex | -x² + y, -x³ + z | y³ - z², x² - y, xy - z, xz - y² |
| Parametric Curve | x - y², -y³ + z | x - y², y³ - z |
| Parametric Curve Grlex | x - y², -y³ + z | x² - yz, xy - z, -x + y² |
| Twisted Cubic | x - z², y - z³ | x - z², y - z³ |
| Twisted Cubic Grlex | x - z², y - z³ | x² - yz, xz - y, -x + z² |
| Variety Intersection | -y² + z, x - y³ | x - yz, y² - z |
| Variety Intersection Grlex | -y² + z, x - y³ | -x² + z³, xy - z², y² - z, -x + yz |
| Space Curve | y - z², x - z³ | x - z³, y - z² |
| Space Curve Grlex | y - z², x - z³ | -x² + y³, xz - y², -x + yz, -y + z² |
| Circle-Parabola Intersection | 4x²y² + 4xy + 1, x² + y² - 1 | x² + y² - 1, 4x²y² + 4xy + 1 |
| Katsura-3 System | x₀ + 2x₁ + 2x₂ - 1, x₀² + 2x₁² + 2x₂² - x₀, 2x₀x₁ + 2x₁x₂ - x₁ | x₀ + 2x₁ + 2x₂ - 1, x₁ + 2x₂² - 1/2, x₂⁴ - x₂² + 1/4 |
| Katsura-3 Grlex | x₀ + 2x₁ + 2x₂ - 1, x₀² + 2x₁² + 2x₂² - x₀, 2x₀x₁ + 2x₁x₂ - x₁ | x₀ + 2x₁ + 2x₂ - 1, x₁ + 2x₂² - 1/2, x₂⁴ - x₂² + 1/4 |
| Cyclic-4 System | a + b + c + d, ab + ad + bc + bd, abc + abd + acd + bcd, abcd - 1 | (see code for full basis) |
| Cyclic-4 Grlex | a + b + c + d, ab + ad + bc + bd, abc + abd + acd + bcd, abcd - 1 | (see code for full basis) |
| Simple Ideal | x² - y, xy - 1 | x² - y, xy - 1, y² - y |
| Monomial Ordering Properties | x²y + x², xy - 1 | (see code for both orderings) |
| Reduction Properties | x² + xy, x - y | (see code for reduction) |
| S-Polynomial | x² + 2xy², xy + 2y³ - 1 | (see code for S-polynomial) |
| Empty and Single Polynomial | (empty), x | x |
| Zero Polynomial | 0, x | x |
| Constant Polynomial | 1, x | 1 |
| Linear System | x + y - 1, x - y | x - 1, y |
| Quadratic System | x² + y² - 1, x² - y² | x² - y², y² - 1 |
| Homogeneous Ideal | x² + y², xy | x² + y², xy |
| Principal Ideal | x² + xy + y² | x² + xy + y² |
| Elimination Ideal | x² + y + z, xy + z² | (see code for basis) |
| Four Variable System | x - y, y - z, z - w | x - y, y - z, z - w |
| Radical Membership | x³, x²y | x³, x²y |
| Intersection of Ideals | x, y | x, y |
| Quotient Ideal | x² - y², xy | x² - y², xy |
| Syzygy Computation | x² + xy, xy + y² | (see code for S-polynomial reduction) |
| Cyclic-4 Known Result | v₁v₂v₃v₄ - 1, v₁v₂v₃ + v₁v₂v₄ + v₁v₃v₄ + v₂v₃v₄, v₁v₂ + v₂v₃ + v₁v₄ + v₃v₄, v₁ + v₂ + v₃ + v₄ | (see code for full basis) |
| Known Groebner Result | x² - 1, y - 1 | x² - 1, y - 1 |
| Katsura-2 System | x + 2y - 1, x² - x + 2y² | x + 2y - 1, y² - (1/3)y |
| Polynomial Division Example | x³ - 2xy, x²y + x - 2y² | x - 2y², y³ |
| System with Common Factor | x(x + y - 1), y(x + y - 1) | x² + xy - x, xy + y² - y |
| Three Spheres Intersection | x² + y² - 1, y² + z² - 1, x² + z² - 1 | y² + z² - 1, z² - (1/2) |
| Inconsistent System | x² + y², x² + y² + 1 | 1 |
| Single Variable Polynomial | x³ - 2x + 1 | x³ - 2x + 1 |
| System with Rational Coefficients | x² + y/2 - 1, 3x/4 + y² - 1/3 | x + 4y²/3 - 4/9, y⁴ - 2y²/3 + 9y/32 - 65/144 |