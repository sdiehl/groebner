# Groebner Basis Test Suite

| Test Name | System (Input Polynomials) | Canonical Groebner Basis (Output) |
|-----------|---------------------------|-----------------------------------|
| Basic Cox-Little-O'Shea | x² + 2xy², xy + 2y³ - 1 | x, y³ - 1/2 |
| Three Variable Simple | x - z², y - z³ | x - z², y - z³ |
| Grlex Ordering | x³ - 2xy, x²y + x - 2y² | (see code, structure checked) |
| Three Variable Complex Lex | -x² + y, -x³ + z | (see code, structure checked) |
| Parametric Curve | x - y², -y³ + z | (see code, structure checked) |
| Twisted Cubic | x - z², y - z³ | x - z², y - z³ |
| Variety Intersection | -y² + z, x - y³ | (see code, structure checked) |
| Space Curve | y - z², x - z³ | (see code, structure checked) |
| Circle-Parabola Intersection | 4x²y² + 4xy + 1, x² + y² - 1 | (see code, structure checked) |
| Katsura-3 System | x₀ + 2x₁ + 2x₂ - 1, x₀² + 2x₁² + 2x₂² - x₀, 2x₀x₁ + 2x₁x₂ - x₁ | (see code, structure checked) |
| Cyclic-4 System | a + b + c + d, ab + ad + bc + bd, abc + abd + acd + bcd, abcd - 1 | (see code, structure checked) |
| Simple Ideal | x² - y, xy - 1 | (see code, structure checked) |
| Linear System | x + y - 1, x - y | (see code, structure checked) |
| Quadratic System | x² + y² - 1, x² - y² | (see code, structure checked) |
| Principal Ideal | x² + xy + y² | (see code, structure checked) |
| Four Variable System | x - y, y - z, z - w | (see code, structure checked) |
| Katsura-2 System | x + 2y - 1, x² - x + 2y² | x + 2y - 1, y² - (1/3)y |
| Polynomial Division Example | x³ - 2xy, x²y + x - 2y² | x - 2y², y³ |
| System with Common Factor | x(x + y - 1), y(x + y - 1) | x² + xy - x, xy + y² - y |
| Three Spheres Intersection | x² + y² - 1, y² + z² - 1, x² + z² - 1 | y² + z² - 1, z² - (1/2) |
| Inconsistent System | x² + y², x² + y² + 1 | 1 |
| Single Variable Polynomial | x³ - 2x + 1 | x³ - 2x + 1 |
| System with Rational Coefficients | x² + y/2 - 1, 3x/4 + y² - 1/3 | x + 4y²/3 - 4/9, y⁴ - 2y²/3 + 9y/32 - 65/144 |