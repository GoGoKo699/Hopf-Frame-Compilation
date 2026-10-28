"""A finite physical Clifford seed with complete noncommuting copy grading."""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, CQ2, Q2, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_field_shear_intertwiner import _permutation_matrix, _shear


def _pauli(x, z):
    phase = (CO, CI, -CO, -CI)[(x & z).bit_count() % 4]
    return [[phase * (-1 if (z & column).bit_count() % 2 else 1)
             if row == column ^ x else CZ
             for column in range(16)] for row in range(16)]


def _trace(matrix):
    return sum((matrix[i][i] for i in range(len(matrix))), CZ)


def _expand_clifford_macros(word):
    result = []
    for gate in word:
        if gate[0] == 'Z':
            result += [('S', gate[1])] * 2
        elif gate[0] == 'CZ':
            _, source, target = gate
            result += [('H', target), ('CX', source, target), ('H', target)]
        else:
            result.append(gate)
    return result


def _seed_words():
    # Macro gates here are individually Hermitian, so reversal is the
    # actual adjoint. Expanding Z afterward gives its reduced SS word.
    ka = [('Z', 1), ('CZ', 0, 1), ('CX', 3, 1)]
    kb = [('CZ', 0, 1), ('CZ', 0, 2)]
    seed = list(reversed(kb)) + [('H', 0)] + kb
    seed += list(reversed(ka)) + [('H', 1)] + ka
    return (_expand_clifford_macros(seed),
            _expand_clifford_macros(list(reversed(seed))))


def _seed_and_orbit():
    seed, _ = _seed_words()
    r = native_columns(4, seed)
    queries = [_permutation_matrix([_shear(2, x, z) for z in range(16)])
               for x in range(4)]
    return r, [multiply(multiply(q, r), q) for q in queries]


class FieldCopyCliffordPilotTests(unittest.TestCase):
    def test_literal_seed_factors_and_complete_native_inverse(self):
        invsqrt2 = CQ2(Q2(0, Fraction(1, 2)))
        a = scale(add(scale(_pauli(2, 1), -CO), _pauli(0, 10)), invsqrt2)
        b = scale(add(_pauli(1, 6), _pauli(0, 1)), invsqrt2)
        self.assertEqual(multiply(a, a), eye(16))
        self.assertEqual(multiply(b, b), eye(16))
        self.assertEqual(multiply(a, b), multiply(b, a))
        seed, inverse = _seed_words()
        actual = native_columns(4, seed)
        self.assertEqual(actual, multiply(a, b))
        self.assertEqual(native_columns(4, inverse), adjoint(actual))
        self.assertEqual(actual, adjoint(actual))
        self.assertEqual(multiply(actual, actual), eye(16))
        for word in (seed, inverse):
            self.assertEqual(len(word), 26)
            self.assertEqual({g[0] for g in word}, {'H', 'S', 'CX'})
            self.assertEqual({w for g in word for w in g[1:]}, set(range(4)))

    def test_central_pauli_and_all_literal_commutators(self):
        _, orbit = _seed_and_orbit()
        p = _pauli(4, 1)
        expected = {(0, 1): scale(p, -CO), (0, 2): eye(16), (0, 3): p,
                    (1, 2): p, (1, 3): eye(16), (2, 3): scale(p, -CO)}
        self.assertEqual(_trace(p), CZ)
        for r in orbit:
            self.assertEqual(multiply(r, r), eye(16))
            self.assertEqual(r, adjoint(r))
            self.assertEqual(multiply(r, p), multiply(p, r))
        for (x, y), target in expected.items():
            word = eye(16)
            for index in (x, y, x, y):
                word = multiply(word, orbit[index])
            self.assertEqual(word, target)

    def test_sector_rank_and_all_parity_word_traces(self):
        _, orbit = _seed_and_orbit()
        p = _pauli(4, 1)
        for sign in (1, -1):
            pair_signs = {(0, 1): -sign, (0, 2): 1, (0, 3): sign,
                          (1, 2): sign, (1, 3): 1, (2, 3): -sign}
            rows = [sum((pair_signs[tuple(sorted((x, y)))] == -1) << y
                        for y in range(4) if y != x) for x in range(4)]
            span = {0}
            for row in rows:
                span |= {value ^ row for value in tuple(span)}
            self.assertEqual(len(span), 16)
        for mask in range(1, 16):
            word = eye(16)
            for index in range(4):
                if mask >> index & 1:
                    word = multiply(word, orbit[index])
            self.assertEqual(_trace(word), CZ)
            self.assertEqual(_trace(multiply(p, word)), CZ)
