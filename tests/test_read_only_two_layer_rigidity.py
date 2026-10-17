"""Exact bounded checks of the preserved-address two-layer proof.

The fixtures retain interlayer address phases and the actual zero-address
word. They also audit the native scope counterexample: unrestricted dirty
Toffoli has two T layers but uses transformed address wires as T targets.
No finite enumeration supplies a uniform-depth lower bound.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, adjoint, multiply, native_columns, scale,
)
from tests.test_dirty_indicator_depth import _toffoli_two_layers
from tests.test_t_depth import _word


def _cz(control, target):
    return [('H', target), ('CX', control, target), ('H', target)]


def _fibre(matrix, address_bits, x):
    width = len(matrix).bit_length() - 1 - address_bits
    indices = [x | (w << address_bits) for w in range(1 << width)]
    return [[matrix[i][j] for j in indices] for i in indices]


def _pauli(width, x_mask, z_mask):
    # Matrix X^x Z^z, with no unrecorded Hermitian phase.
    word = [('Z', j) for j in range(width) if (z_mask >> j) & 1]
    word += [('X', j) for j in range(width) if (x_mask >> j) & 1]
    return native_columns(width, word)


def _identify_projective_pauli(matrix):
    width = len(matrix).bit_length() - 1
    phases = (CO, CI, -CO, -CI)
    for x, z in product(range(1 << width), repeat=2):
        canonical = _pauli(width, x, z)
        for exponent, scalar in enumerate(phases):
            if matrix == scale(canonical, scalar):
                return x, z, exponent, canonical
    raise AssertionError('not a projective Pauli')


class ReadOnlyTwoLayerRigidityTests(unittest.TestCase):
    def test_clifford_interlayer_normal_form_keeps_address_phases(self):
        a0, a1, work = 0, 1, 2
        words = [
            [('H', work), ('S', a0), ('CX', a0, work)] + _cz(a1, work),
            [('S', work), ('CX', a1, work)] + _cz(a0, a1)
            + [('H', work), ('CX', a0, work)],
            [('SDG', a1), ('H', work)] + _cz(a0, work)
            + [('CX', a1, work)],
        ]
        for word in words:
            full = native_columns(3, word)
            for i, row in enumerate(full):
                for j, value in enumerate(row):
                    if (i & 3) != (j & 3):
                        self.assertEqual(value, CZ)
            c0 = _fibre(full, 2, 0)
            values = []
            for x in range(4):
                fibre = _fibre(full, 2, x)
                p = multiply(fibre, adjoint(c0))
                values.append(_identify_projective_pauli(p))
                _, _, exponent, canonical = values[-1]
                self.assertEqual(fibre,
                                 scale(multiply(canonical, c0),
                                       (CO, CI, -CO, -CI)[exponent]))
            for component in (0, 1):
                self.assertEqual(values[0][component], 0)
                self.assertEqual(values[3][component],
                                 values[1][component] ^ values[2][component])
            exponents = [item[2] for item in values]
            self.assertEqual(exponents[0], 0)
            self.assertEqual((exponents[3] - exponents[1] - exponents[2]) % 2, 0)

    def test_three_band_identity_uses_the_actual_zero_address_word(self):
        a0, a1, work = 0, 1, 2
        c_words = [
            [('H', work), ('S', a0), ('CX', a0, work)] + _cz(a1, work),
            [('S', work), ('CX', a1, work)] + _cz(a0, a1)
            + [('H', work), ('CX', a0, work)],
            [('SDG', a1), ('H', work)] + _cz(a0, work)
            + [('CX', a1, work)],
        ]
        full_word = c_words[0] + [('T', work)] + c_words[1]
        full_word += [('TDG', work)] + c_words[2]
        full = native_columns(3, full_word)
        c = [native_columns(3, word) for word in c_words]
        zero = [_fibre(matrix, 2, 0) for matrix in c]
        u0 = _fibre(full, 2, 0)
        d2 = native_columns(1, [('TDG', 0)])
        for x in range(4):
            p = [multiply(_fibre(matrix, 2, x), adjoint(z))
                 for matrix, z in zip(c, zero)]
            # The complete scalar is retained inside p; no phase quotient.
            m = multiply(multiply(multiply(multiply(zero[2], d2), p[1]),
                                  adjoint(d2)), adjoint(zero[2]))
            r = multiply(multiply(adjoint(zero[0]), p[0]), zero[0])
            reconstructed = multiply(multiply(multiply(p[2], m), u0), r)
            self.assertEqual(_fibre(full, 2, x), reconstructed)
        # A proof must not replace U_0 by a Pauli before its error argument.
        with self.assertRaises(AssertionError):
            _identify_projective_pauli(u0)

    def test_physical_t_layer_has_flat_pauli_expansion(self):
        width = 2
        layer = native_columns(width, [('T', 0), ('TDG', 1)])
        for x_mask, z_mask in product(range(4), repeat=2):
            p = _pauli(width, x_mask, z_mask)
            image = multiply(multiply(layer, p), adjoint(layer))
            coefficients = []
            for x, z in product(range(4), repeat=2):
                candidate = _pauli(width, x, z)
                relative = multiply(adjoint(candidate), image)
                coefficient = sum((relative[i][i] for i in range(4)), CZ) / 4
                if coefficient != CZ:
                    coefficients.append(coefficient)
            k = x_mask.bit_count()
            self.assertEqual(len(coefficients), 1 << k)
            self.assertTrue(all(value.abs2() == Q2(Fraction(1, 1 << k))
                                for value in coefficients))

    def test_unrestricted_two_layer_toffoli_violates_read_only_scope(self):
        schedule = _toffoli_two_layers(0, 1, 2, 3)
        layers = [gates for kind, gates in schedule if kind == 'T']
        self.assertEqual(len(layers), 2)
        self.assertTrue(any(gate[1] in (0, 1) for gates in layers for gate in gates))
        actual = native_columns(4, _word(schedule))
        for column in range(16):
            row = column ^ (4 if (column & 1) and (column & 2) else 0)
            self.assertEqual([actual[i][column] for i in range(16)],
                             [CO if i == row else CZ for i in range(16)])


if __name__ == '__main__':
    unittest.main()
