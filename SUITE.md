# Groebner Basis Test Suite

| Test Name | System (Input Polynomials) | Canonical Groebner Basis (Output) |
|-----------|---------------------------|-----------------------------------|
| Basic Cox-Little-O'Shea | $x^2 + 2xy^2$, $xy + 2y^3 - 1$ | $x$, $y^3 - \frac{1}{2}$ |
| Groebner y,x Order | $2x^2y + y^2$, $2x^3 + xy - 1$ | $y$, $x^3 - \frac{1}{2}$ |
| Three Variable Simple | $x - z^2$, $y - z^3$ | $x - z^2$, $y - z^3$ |
| Grlex Ordering | $x^3 - 2xy$, $x^2y + x - 2y^2$ | $x^2$, $xy$, $-\frac{1}{2}x + y^2$ |
| Three Variable Complex Lex | $-x^2 + y$, $-x^3 + z$ | $x^2 - y$, $xy - z$, $xz - y^2$, $y^3 - z^2$ |
| Three Variable Complex Grlex | $-x^2 + y$, $-x^3 + z$ | $y^3 - z^2$, $x^2 - y$, $xy - z$, $xz - y^2$ |
| Parametric Curve | $x - y^2$, $-y^3 + z$ | $x - y^2$, $y^3 - z$ |
| Parametric Curve Grlex | $x - y^2$, $-y^3 + z$ | $x^2 - yz$, $xy - z$, $-x + y^2$ |
| Twisted Cubic | $x - z^2$, $y - z^3$ | $x - z^2$, $y - z^3$ |
| Twisted Cubic Grlex | $x - z^2$, $y - z^3$ | $x^2 - yz$, $xz - y$, $-x + z^2$ |
| Variety Intersection | $-y^2 + z$, $x - y^3$ | $x - y^3$, $y^2 - z$ |
| Variety Intersection Grlex | $-y^2 + z$, $x - y^3$ | $-x^2 + z^3$, $xy - z^2$, $y^2 - z$, $-x + yz$ |
| Space Curve | $y - z^2$, $x - z^3$ | $x - z^3$, $y - z^2$ |
| Space Curve Grlex | $y - z^2$, $x - z^3$ | $y^3 - x^2$, $xz - y^2$, $yz - x$, $z^2 - y$ |
| Circle-Parabola Intersection | $4x^2y^2 + 4xy + 1$, $x^2 + y^2 - 1$ | $x^2 + y^2 - 1$, $4x^2y^2 + 4xy + 1$ |
| Katsura-3 System | $x_0 + 2x_1 + 2x_2 - 1$, $x_0^2 + 2x_1^2 + 2x_2^2 - x_0$, $2x_0x_1 + 2x_1x_2 - x_1$ | $x_0 + 2x_1 + 2x_2 - 1$, $x_1 + 2x_2^2 - \frac{1}{2}$, $x_2^4 - x_2^2 + \frac{1}{4}$ |
| Katsura-3 Grlex | $x_0 + 2x_1 + 2x_2 - 1$, $x_0^2 + 2x_1^2 + 2x_2^2 - x_0$, $2x_0x_1 + 2x_1x_2 - x_1$ | $x_0 + 2x_1 + 2x_2 - 1$, $x_1 + 2x_2^2 - \frac{1}{2}$, $x_2^4 - x_2^2 + \frac{1}{4}$ |
| Cyclic-4 System | $a + b + c + d$, $ab + ad + bc + bd$, $abc + abd + acd + bcd$, $abcd - 1$ | (see test_cyclic_4_system) |
| Cyclic-4 Grlex | $a + b + c + d$, $ab + ad + bc + bd$, $abc + abd + acd + bcd$, $abcd - 1$ | (see test_cyclic_4_grlex) |
| Simple Ideal | $x^2 - y$, $xy - 1$ | $x^2 - y$, $xy - 1$, $y^2 - y$ |
| Monomial Ordering Properties | $x^2y + x^2$, $xy - 1$ | (see test_monomial_ordering_properties) |
| Reduction Properties | $x^2 + xy$, $x - y$ | (see test_reduction_properties) |
| S-Polynomial | $x^2 + 2xy^2$, $xy + 2y^3 - 1$ | (see test_s_polynomial) |
| Empty and Single Polynomial | (empty), $x$ | $x$ |
| Zero Polynomial | $0$, $x$ | $x$ |
| Constant Polynomial | $1$, $x$ | $1$ |
| Linear System | $x + y - 1$, $x - y$ | $x - 1$, $y$ |
| Quadratic System | $x^2 + y^2 - 1$, $x^2 - y^2$ | $x^2 - y^2$, $y^2 - 1$ |
| Homogeneous Ideal | $x^2 + y^2$, $xy$ | $x^2 + y^2$, $xy$ |
| Principal Ideal | $x^2 + xy + y^2$ | $x^2 + xy + y^2$ |
| Elimination Ideal | $x^2 + y + z$, $xy + z^2$ | (see test_elimination_ideal) |
| Four Variable System | $x - y$, $y - z$, $z - w$ | $x - y$, $y - z$, $z - w$ |
| Radical Membership | $x^3$, $x^2y$ | $x^3$, $x^2y$ |
| Intersection of Ideals | $x$, $y$ | $x$, $y$ |
| Quotient Ideal | $x^2 - y^2$, $xy$ | $x^2 - y^2$, $xy$ |
| Syzygy Computation | $x^2 + xy$, $xy + y^2$ | (see test_syzygy_computation) |
| Cyclic-4 Known Result | $v_1v_2v_3v_4 - 1$, $v_1v_2v_3 + v_1v_2v_4 + v_1v_3v_4 + v_2v_3v_4$, $v_1v_2 + v_2v_3 + v_1v_4 + v_3v_4$, $v_1 + v_2 + v_3 + v_4$ | (see test_cyclic_4_known_result) |
| Known Groebner Result | $x^2 - 1$, $y - 1$ | $x^2 - 1$, $y - 1$ |
| Katsura-2 System | $x + 2y - 1$, $x^2 - x + 2y^2$ | $x + 2y - 1$, $y^2 - \frac{1}{3}y$ |
| Polynomial Division Example | $x^3 - 2xy$, $x^2y + x - 2y^2$ | $x - 2y^2$, $y^3$ |
| System with Common Factor | $x(x + y - 1)$, $y(x + y - 1)$ | $x^2 + xy - x$, $xy + y^2 - y$ |
| Three Spheres Intersection | $x^2 + y^2 - 1$, $y^2 + z^2 - 1$, $x^2 + z^2 - 1$ | $y^2 + z^2 - 1$, $z^2 - \frac{1}{2}$ |
| Inconsistent System | $x^2 + y^2$, $x^2 + y^2 + 1$ | $1$ |
| Single Variable Polynomial | $x^3 - 2x + 1$ | $x^3 - 2x + 1$ |
| System with Rational Coefficients | $x^2 + \frac{y}{2} - 1$, $\frac{3x}{4} + y^2 - \frac{1}{3}$ | $x + \frac{4}{3}y^2 - \frac{4}{9}$, $y^4 - \frac{2}{3}y^2 + \frac{9}{32}y - \frac{65}{144}$ |