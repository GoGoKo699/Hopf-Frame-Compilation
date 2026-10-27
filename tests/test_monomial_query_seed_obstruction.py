"""Exact seed-class identity, arbitrary phases, helpers and its scope."""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, CI, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_field_shear_intertwiner import _permutation_matrix, _shear


def _diagonal(values):
    return [[value if i == j else CZ for j in range(len(values))]
            for i, value in enumerate(values)]


def _orbit(queries, seed):
    return [multiply(multiply(q, seed), adjoint(q)) for q in queries]


def _plane(orbit, c, a, b):
    result = eye(len(orbit[0]))
    for factor in (orbit[c ^ a], adjoint(orbit[c]), orbit[c ^ b],
                   adjoint(orbit[c ^ a ^ b])):
        result = multiply(result, factor)
    return result


class MonomialQuerySeedObstructionTests(unittest.TestCase):
    def test_nonhermitian_phase_seed_on_every_shifted_three_bit_plane(self):
        queries = [_diagonal([CO * (-1) ** ((x & z).bit_count() % 2)
                              for z in range(8)]) for x in range(8)]
        permutation = _permutation_matrix([(z + 1) % 8 for z in range(8)])
        phases = _diagonal([CO, CI, -CO, -CI, CI, CO, -CI, -CO])
        seed = multiply(permutation, phases)
        self.assertNotEqual(seed, adjoint(seed))
        self.assertNotEqual(multiply(seed, seed), eye(8))
        orbit = _orbit(queries, seed)
        for c in range(8):
            for a in range(1, 8):
                for b in range(1, 8):
                    if a != b:
                        self.assertEqual(_plane(orbit, c, a, b), eye(8))
        # The actual-adjoint convention is essential for this seed.
        wrong = eye(8)
        for x in (1, 0, 2, 3):
            wrong = multiply(wrong, orbit[x])
        self.assertNotEqual(wrong, eye(8))

    def test_field_fourier_basis_with_arbitrary_helper_mixing(self):
        # Four field-bank wires plus one arbitrary helper.  Seed mixes the
        # helper with the bank in the diagonal query coordinates.
        queries = [_permutation_matrix([_shear(2, x, z & 15) | (z & 16)
                                        for z in range(32)]) for x in range(4)]
        fourier = native_columns(5, [('H', 2), ('H', 3)])
        diagonal_queries = [multiply(multiply(fourier, q), fourier)
                            for q in queries]
        for q in diagonal_queries:
            self.assertTrue(all(q[i][j] == CZ for i in range(32)
                                for j in range(32) if i != j))
        mono = native_columns(5, [('CX', 4, 0), ('T', 0), ('CX', 1, 4),
                                  ('S', 2), ('CX', 3, 1)])
        seed = multiply(multiply(fourier, mono), fourier)
        self.assertNotEqual(seed, adjoint(seed))
        orbit = _orbit(queries, seed)
        for c in range(4):
            self.assertEqual(_plane(orbit, c, 1, 2), eye(32))
        # Four distinct output labels give a nonidentity X word.  The
        # antisymmetric pair of computational columns is an eigenvector -1.
        output_word = _permutation_matrix([z ^ 15 for z in range(16)])
        vector = [CO if z == 0 else -CO if z == 15 else CZ for z in range(16)]
        image = [sum((output_word[i][j] * vector[j] for j in range(16)), CZ)
                 for i in range(16)]
        self.assertEqual(image, [-value for value in vector])

    def test_nonmonomial_seed_can_break_the_plane_identity(self):
        queries = [_diagonal([CO * (-1) ** ((x & z).bit_count() % 2)
                              for z in range(4)]) for x in range(4)]
        seed = native_columns(2, [('H', 0), ('CX', 0, 1), ('T', 1), ('H', 1)])
        orbit = _orbit(queries, seed)
        plane = _plane(orbit, 0, 1, 2)
        self.assertNotEqual(plane, eye(4))
        self.assertFalse(all(plane[i][j] == (plane[0][0] if i == j else CZ)
                             for i in range(4) for j in range(4)))
        transported = multiply(multiply(seed, queries[1]), adjoint(seed))
        self.assertNotEqual(multiply(queries[2], transported),
                            multiply(transported, queries[2]))


    def test_quantitative_cross_commutator_and_hermitian_part(self):
        queries = [_diagonal([CO * (-1) ** ((x & z).bit_count() % 2)
                              for z in range(4)]) for x in range(4)]
        xx = _permutation_matrix([z ^ 3 for z in range(4)])
        # This exact rational unitary is deliberately not a Hermitian
        # reflection. No native synthesis is asserted for the rational case.
        rational = add(scale(eye(4), Fraction(24, 25)),
                       scale(xx, CI * Fraction(7, 25)))
        # A literal native Pauli T rotation, with its scalar phase retained.
        boundary = native_columns(2, [('H', 0), ('H', 1), ('CX', 0, 1),
                                       ('T', 1), ('CX', 0, 1),
                                       ('H', 0), ('H', 1)])
        for seed, kappa_squared in (
                (rational, Fraction(672, 625) ** 2), (boundary, Fraction(2))):
            self.assertEqual(multiply(seed, adjoint(seed)), eye(4))
            self.assertNotEqual(seed, adjoint(seed))
            b_a = multiply(multiply(seed, queries[1]), adjoint(seed))
            orbit = _orbit(queries, seed)
            w = _plane(orbit, 0, 1, 2)
            formula = eye(4)
            for factor in (queries[1], b_a, queries[2], b_a,
                           queries[2], queries[1]):
                formula = multiply(formula, factor)
            self.assertEqual(w, formula)
            commutator = add(multiply(b_a, queries[2]),
                             scale(multiply(queries[2], b_a), -1))
            displacement = add(w, scale(eye(4), -1))
            self.assertEqual(multiply(adjoint(commutator), commutator),
                             scale(eye(4), kappa_squared))
            self.assertEqual(multiply(adjoint(displacement), displacement),
                             scale(eye(4), kappa_squared))
            real_part = scale(add(w, adjoint(w)), Fraction(1, 2))
            self.assertEqual(real_part,
                             scale(eye(4), 1 - kappa_squared / 2))
            if seed == rational:
                self.assertGreater((2 - kappa_squared) / 4, 0)
            else:
                self.assertEqual(real_part, scale(eye(4), 0))
