import unittest

from groebner import GroebnerError, Order, Ring, contains, groebner_basis, reduce

CIRCLE = ["x^2 + y^2 - 1", "x - y"]


class TestGroebner(unittest.TestCase):
    def test_rational(self):
        basis = groebner_basis(CIRCLE, order=Order.LEX)
        self.assertEqual(basis, ["x - y", "y^2 - 1/2"])

    def test_prime_field(self):
        ring = Ring(["x", "y"], Order.LEX, modulus=32003)
        basis = ring.groebner_basis(["x^2 - y", "x*y - 1"])
        self.assertEqual(basis, ["x + 32002*y^2", "y^3 + 32002"])
        self.assertTrue(ring.contains("x^2 - y", basis))

    def test_reduce(self):
        self.assertEqual(reduce("x^2", CIRCLE, order=Order.LEX), "1/2")
        self.assertTrue(contains("2*y^2 - 1", CIRCLE))
        self.assertFalse(contains("y", CIRCLE))

    def test_errors(self):
        with self.assertRaises(GroebnerError):
            groebner_basis(["x^2 +"], ["x"])
        with self.assertRaises(GroebnerError):
            Ring(["x", "x"])
        with self.assertRaises(TypeError):
            groebner_basis("x^2 - y")


if __name__ == "__main__":
    unittest.main()
