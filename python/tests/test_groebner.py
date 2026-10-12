import unittest

from groebner import GroebnerError, contains, groebner_basis, reduce


class TestGroebner(unittest.TestCase):
    def test_rational(self):
        basis = groebner_basis(["x^2 + y^2 - 1", "x - y"], order="lex")
        self.assertEqual(basis, ["x - y", "y^2 - 1/2"])

    def test_prime_field(self):
        basis = groebner_basis("x^2 - y; x*y - 1", "x, y", "lex", 32003)
        self.assertEqual(basis, ["x + 32002*y^2", "y^3 + 32002"])
        self.assertTrue(contains("x^2 - y", basis, "x, y", "lex", 32003))

    def test_reduce(self):
        polys = ["x^2 + y^2 - 1", "x - y"]
        self.assertEqual(reduce("x^2", polys, order="lex"), "1/2")
        self.assertTrue(contains("2*y^2 - 1", polys))
        self.assertFalse(contains("y", polys))

    def test_errors(self):
        with self.assertRaises(GroebnerError):
            groebner_basis("x^2 +", vars="x")
        with self.assertRaises(GroebnerError):
            groebner_basis("x^2", order="nope")


if __name__ == "__main__":
    unittest.main()
