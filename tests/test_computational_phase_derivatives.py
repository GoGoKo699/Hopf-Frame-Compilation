"""Exact actual-code gauges and completed computational derivatives.

The occupied target is CZ, not a sharpness example for the quartic bound.
All three available C-active clean flips give the same derivative masks.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import (
    _clifford_columns, _conjugate, _symplectic, _tau,
)


def _block(matrix, output_clean, input_clean):
    size = len(matrix) // 2
    return [[matrix[2 * row + output_clean][2 * col + input_clean]
             for col in range(size)] for row in range(size)]


def _diagonal_blocks(top, bottom):
    size = len(top)
    result = [[CZ for _ in range(2 * size)] for _ in range(2 * size)]
    for row in range(size):
        for col in range(size):
            result[2 * row][2 * col] = top[row][col]
            result[2 * row + 1][2 * col + 1] = bottom[row][col]
    return result


class ComputationalPhaseDerivativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # The retained three-band r=1 identity split, followed by data CZ.
        k_word = [('H', 1), ('H', 0), ('CX', 0, 1), ('CX', 0, 2)]
        e_word = [('T', wire) for wire in range(3)]
        source = [('CX', 1, 0), ('H', 1)]
        e0 = _inverse(source) + [('TDG', 0), ('T', 1)] + source
        uc0 = k_word + e_word
        uab = [('CX', 0, 2)] + e0 + uc0
        g_word = [('H', 2), ('CX', 1, 2), ('H', 2)]
        cls.word = uab + _inverse(uc0) + g_word
        cls.v = native_columns(3, cls.word)
        cls.kc = native_columns(3, _inverse(g_word) + k_word)
        cls.uc = native_columns(3, _inverse(g_word) + uc0)
        cls.uab = native_columns(3, uab)
        cls.g = _block(cls.v, 0, 0)
        cls.clean_x, cls.clean_z = _tau(1, 3), _tau(8, 3)
        cls.s = _conjugate(cls.uc, cls.clean_z)
        cls.r = _tau((1 | 2) << 3, 3)
        cls.normalizations = []
        for wire in range(3):
            axis = _tau(1 << (3 + wire), 3)
            q_out = _conjugate(adjoint(cls.kc), axis)
            # q_out=[[0,B],[B†,0]], D=diag(I,B†) fixes clean0.
            b = _block(q_out, 0, 1)
            d = _diagonal_blocks(eye(4), adjoint(b))
            v_prime = multiply(adjoint(d), cls.v)
            logical_l = multiply(adjoint(cls.g), _block(v_prime, 1, 1))
            cls.normalizations.append((axis, q_out, d, v_prime, logical_l))

    def test_actual_native_word_and_all_clean_flip_clifford_gauges(self):
        self.assertEqual(native_columns(3, _inverse(self.word)), adjoint(self.v))
        self.assertEqual(self.v, multiply(adjoint(self.uc), self.uab))
        self.assertEqual(_conjugate(self.v, self.clean_z), self.clean_z)
        for value in range(4):
            expected = -CO if value == 3 else CO
            for row in range(8):
                self.assertEqual(self.v[row][2 * value],
                                 expected if row == 2 * value else CZ)
        all_masks = []
        for axis, q_out, d, v_prime, logical_l in self.normalizations:
            self.assertEqual(_conjugate(d, self.clean_x), q_out)
            self.assertEqual(_block(d, 0, 0), eye(4))
            self.assertEqual(v_prime, _diagonal_blocks(
                self.g, multiply(self.g, logical_l)))
            q_in = _conjugate(adjoint(self.uab), axis)
            self.assertEqual(_conjugate(adjoint(v_prime), self.clean_x), q_in)
            self.assertEqual(_block(q_in, 0, 1), logical_l)
            self.assertEqual(_block(q_in, 1, 0), adjoint(logical_l))
            columns = _clifford_columns(logical_l, 2)
            all_masks.append([columns[2] & 3, columns[3] & 3])
        # This space has dimension one; changing the actual clean flip adds none.
        self.assertEqual(all_masks, [[1, 0]] * 3)

    def test_literal_central_pullback_and_every_computational_derivative(self):
        plus = scale(add(eye(8), self.s), Fraction(1, 2))
        minus = scale(add(eye(8), scale(self.s, -1)), Fraction(1, 2))
        for _, _, d, _, logical_l in self.normalizations:
            uc_prime = multiply(self.uc, d)
            kc_prime = multiply(self.kc, d)
            for z in range(4):
                logical_z = _tau(z << 2, 2)
                physical_z = _tau((z << 1) << 3, 3)
                transformed_z = _conjugate(d, physical_z)
                chi = transformed_z != physical_z
                self.assertEqual(transformed_z,
                                 multiply(physical_z, self.clean_z)
                                 if chi else physical_z)
                alpha = scale(self.r, CI) if z & 1 else eye(8)
                if bool(z & 2) != chi:
                    alpha = scale(alpha, -1)
                c = _conjugate(uc_prime, physical_z)
                a = _conjugate(self.uab, physical_z)
                self.assertEqual(multiply(c, a),
                                 add(plus, multiply(alpha, minus)))
                central_pullback = _conjugate(adjoint(kc_prime), alpha)
                central_logical = _block(central_pullback, 0, 0)
                self.assertEqual(central_pullback,
                                 _diagonal_blocks(central_logical, central_logical))
                p = multiply(logical_z, _conjugate(logical_l, logical_z))
                self.assertEqual(_conjugate(self.g, p), central_logical)
                derivative = multiply(central_logical, adjoint(p))
                # Here E_s is empty, so the derived <=2-band word is Clifford.
                mask = z & 1
                for value in range(4):
                    phase = (-1 if ((value == 3) != ((value ^ mask) == 3))
                             else 1)
                    for row in range(4):
                        self.assertEqual(derivative[row][value],
                                         phase * CO
                                         if row == value else CZ)

    def test_dependent_axis_guard_for_boolean_reality(self):
        # Four commuting axes Z0,Z1,Z2,Z0Z1Z2 have one relation.
        # This is a two-band physical word, not one independent native band.
        parity = [('CX', 0, 2), ('CX', 1, 2)]
        word = [('T', wire) for wire in range(3)]
        word += parity + [('T', 2)] + _inverse(parity)
        e = native_columns(3, word)
        self.assertEqual(native_columns(3, _inverse(word)), adjoint(e))
        r = _tau(7, 3)  # X0 X1 X2 is active for every one of the four axes.
        derivative = multiply(multiply(r, adjoint(e)), multiply(r, e))
        expected = [1 if value in (0, 7) else 0 for value in range(8)]
        for value, bit in enumerate(expected):
            for row in range(8):
                self.assertEqual(derivative[row][value],
                                 (-CO if bit else CO)
                                 if row == value else CZ)
        # Its exact sign is quadratic, so it is not a computational Pauli.
        anf = expected[:]
        for wire in range(3):
            for value in range(8):
                if value & (1 << wire):
                    anf[value] ^= anf[value ^ (1 << wire)]
        self.assertEqual(anf, [1, 1, 1, 1, 1, 1, 1, 0])


    def test_nontrivial_active_basis_matching_keeps_central_cubic_sign(self):
        # f0=X_t, fi=X_t Z_i.  The independent b_i are products of
        # the OTHER three f axes.  Their Gram matrices agree, although
        # their individual axis sets do not.  This is not F/B cancellation.
        width = 4
        f_axes = [1] + [1 | (1 << (width + bit)) for bit in range(1, 4)]
        all_f = 0
        for axis in f_axes:
            all_f ^= axis
        b_axes = [all_f ^ axis for axis in f_axes]
        self.assertNotEqual(set(f_axes), set(b_axes))
        for axes in (f_axes, b_axes):
            self.assertEqual(_rank(axes), 4)
            self.assertTrue(all(not _symplectic(a, b, width)
                                for a in axes for b in axes))
        gram_f = [sum(((axis >> bit) & 1) * (1 << j)
                      for j, axis in enumerate(f_axes))
                  for bit in range(2 * width)]
        gram_b = [sum(((axis >> bit) & 1) * (1 << j)
                      for j, axis in enumerate(b_axes))
                  for bit in range(2 * width)]
        self.assertEqual([[(a & b).bit_count() % 2 for b in gram_f]
                          for a in gram_f],
                         [[(a & b).bit_count() % 2 for b in gram_b]
                          for a in gram_b])

        # Fixed CNOT basis change O=I+11^T, O^2=I, followed by K_F.
        # Every displayed four-T block acts on distinct physical wires.
        k_f_word = [('CX', 0, bit) for bit in range(1, 4)] + [('H', 0)]
        orthogonal_word = [
            ('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1),
            ('CX', 0, 2), ('CX', 0, 3),
            ('CX', 1, 2), ('CX', 1, 3),
            ('CX', 2, 3), ('CX', 3, 2), ('CX', 2, 3),
            ('CX', 2, 0), ('CX', 2, 1),
            ('CX', 3, 0), ('CX', 3, 1),
        ]
        k_b_word = orthogonal_word + k_f_word
        for word, axes in ((k_f_word, f_axes), (k_b_word, b_axes)):
            change = native_columns(width, word)
            for bit, axis in enumerate(axes):
                self.assertEqual(_conjugate(change, _tau(1 << (width + bit), width)),
                                 _tau(axis, width))
        f_word = (_inverse(k_f_word)
                  + [('T', bit) for bit in range(4)] + k_f_word)
        b_word = (_inverse(k_b_word)
                  + [('TDG', bit) for bit in range(4)] + k_b_word)
        c0_word = (_inverse(k_f_word)
                   + [('SDG', bit) for bit in range(4)] + k_f_word)
        c_parity = [('CX', bit, 0) for bit in range(1, 4)]
        cs_word = c_parity + [('T', 0)] + _inverse(c_parity)
        a_word = [('T', 0)]
        word = a_word + b_word + f_word + _inverse(cs_word) + _inverse(c0_word)
        self.assertEqual([len(band) for band in _native_t_bands(word)], [1, 4, 4, 1])
        actual = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))

        f = native_columns(width, f_word)
        b = native_columns(width, b_word)
        c0 = native_columns(width, c0_word)
        d_active = multiply(f, b)
        p, p_prime = _tau(1 << width, width), _tau(15 << width, width)
        self.assertTrue(all(_symplectic(axis, 1 << width, width)
                            for axis in f_axes + b_axes))
        self.assertEqual(_conjugate(d_active, p), p_prime)
        self.assertEqual(_conjugate(c0, p), p_prime)
        self.assertEqual(multiply(d_active, d_active), multiply(c0, c0))
        central = multiply(adjoint(c0), d_active)
        self.assertEqual(multiply(central, central), eye(16))
        self.assertEqual(central, actual)
        signs = []
        for value in range(16):
            z1, z2, z3 = [(value >> bit) & 1 for bit in range(1, 4)]
            bit = z1 * z2 ^ z1 * z3 ^ z2 * z3 ^ z1 * z2 * z3
            signs.append(bit)
            for row in range(16):
                self.assertEqual(actual[row][value],
                                 (-CO if bit else CO) if row == value else CZ)
        # Literal constant phase is +1, with a genuine central cubic term.
        anf = signs[:]
        for bit in range(4):
            for value in range(16):
                if value & (1 << bit):
                    anf[value] ^= anf[value ^ (1 << bit)]
        self.assertEqual([i for i, coefficient in enumerate(anf) if coefficient],
                         [6, 10, 12, 14])


if __name__ == '__main__':
    unittest.main()
