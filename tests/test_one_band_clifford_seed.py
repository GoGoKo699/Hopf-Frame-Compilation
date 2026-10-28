"""Exact fixtures for the arbitrary-Clifford-seed, one-band-copy screen.

The all-size degree/error theorem is in CLIFFORD_SEED_COPY_STRUCTURE.md.
These checks retain signed native phases and arbitrary inactive inputs;
they are not a search for, or an implementation of, an all-size copy.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import (
    CI, CO, Q2, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_field_shear_intertwiner import _pencil_columns
from tests.test_one_band_copy_feasibility import _inverse
from tests.test_pauli_orbit_copy import _pauli


def _coefficients(values, modulus):
    """Integer Boolean Mobius inversion, retaining modulo-four/eight carries."""
    result = list(values)
    for bit in range((len(result) - 1).bit_length()):
        for mask in range(len(result)):
            if mask >> bit & 1:
                result[mask] = (result[mask] - result[mask ^ (1 << bit)]) % modulus
    return result


def _degree(values, modulus):
    return max((mask.bit_count() for mask, value in
                enumerate(_coefficients(values, modulus)) if value), default=0)


def _rank(columns):
    pivots = {}
    for value in columns:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def _linear(columns, value):
    result = 0
    for bit, column in enumerate(columns):
        if value >> bit & 1:
            result ^= column
    return result


def _normal_image(x, z, phase, word):
    """Literal native adjoint on i**phase X**x Z**z, with no phase quotient."""
    for name, *wires in word:
        if name == 'H':
            bit = wires[0]
            a, b = x >> bit & 1, z >> bit & 1
            phase += 2 * a * b
            if a != b:
                x ^= 1 << bit
                z ^= 1 << bit
        elif name == 'S':
            bit = wires[0]
            a = x >> bit & 1
            phase += a
            z ^= a << bit
        elif name == 'CX':
            control, target = wires
            x ^= ((x >> control) & 1) << target
            z ^= ((z >> target) & 1) << control
        else:
            raise ValueError(name)
    return x, z, phase % 4


def _normal_matrix(width, x, z, phase):
    # Remove the Hermitian helper's phase to use ordered X^x Z^z.
    scalar = (CO, CI, -CO, -CI)[(phase - (x & z).bit_count()) % 4]
    return scale(_pauli(width, x | (z << width)), scalar)


class OneBandCliffordSeedTests(unittest.TestCase):
    def test_signed_affine_intersections_and_nonpermutation_cubic_witness(self):
        checked = 0
        for width in range(4):
            size = 1 << width
            for columns in product(range(1, size), repeat=width):
                if _rank(columns) != width:
                    continue
                permutation = all(column.bit_count() == 1 for column in columns)
                for negative in range(size):
                    epsilon = [(-1) ** (negative >> bit & 1) for bit in range(width)]
                    for offset in range(size):
                        phases = [sum(epsilon[bit] *
                                      ((z >> bit & 1) -
                                       ((_linear(columns, z) ^ offset) >> bit & 1))
                                      for bit in range(width)) % 8
                                  for z in range(size)]
                        coefficients = _coefficients(phases, 8)
                        clifford = all(
                            (value % 2 == 0 if mask.bit_count() == 1 else
                             value % 4 == 0 if mask.bit_count() == 2 else value == 0)
                            for mask, value in enumerate(coefficients) if mask)
                        self.assertEqual(clifford, permutation)
                        checked += 1
        self.assertEqual(checked, 10853)

        # A=I+J on four coordinates is orthogonal and invertible, but is not
        # a permutation. Pair overlaps cannot replace the cubic condition.
        columns = tuple(15 ^ (1 << bit) for bit in range(4))
        self.assertEqual(_rank(columns), 4)
        for i, left in enumerate(columns):
            for j, right in enumerate(columns):
                self.assertEqual((left & right).bit_count() % 2, int(i == j))
        phases = [(z.bit_count() - _linear(columns, z).bit_count()) % 8
                  for z in range(16)]
        coefficients = _coefficients(phases, 8)
        self.assertEqual([coefficients[mask] for mask in (7, 11, 13, 14)], [4] * 4)
        self.assertTrue(all(coefficients[1 << bit] % 2 == 0 for bit in range(4)))
        self.assertTrue(all(value % 4 == 0 for mask, value in enumerate(coefficients)
                            if mask.bit_count() == 2))

    def test_native_signed_matching_scalar_and_actual_inverse_with_helper(self):
        width = 3
        # Swap active axes, flip one, and mix an arbitrary inactive helper.
        # X0 is emitted as the literal HSSH word.
        k_word = [('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1),
                  ('H', 0), ('S', 0), ('S', 0), ('H', 0),
                  ('H', 2), ('CX', 0, 2), ('S', 1)]
        native_k = native_columns(width, k_word)
        permutation, signs = [], []
        for bit in range(2):
            x, z, phase = _normal_image(0, 1 << bit, 0, k_word)
            self.assertEqual(x, 0)
            self.assertEqual(z.bit_count(), 1)
            self.assertIn(phase, (0, 2))
            permutation.append(z.bit_length() - 1)
            signs.append(phase // 2)
        self.assertEqual(permutation, [1, 0])
        self.assertEqual(signs, [0, 1])

        # Independently validate the sparse ordered-Pauli phase evaluator
        # against every full native matrix column, including helper values.
        for x, z in product(range(8), repeat=2):
            source = _normal_matrix(width, x, z, 0)
            target = _normal_matrix(width, *_normal_image(x, z, 0, k_word))
            self.assertEqual(multiply(multiply(native_k, source), adjoint(native_k)),
                             target)

        omega = (CO + CI) * Q2(0, Fraction(1, 2))
        roots = [CO]
        for _ in range(7):
            roots.append(roots[-1] * omega)
        for negative in range(4):
            epsilon = [(-1) ** (negative >> bit & 1) for bit in range(2)]
            band = [('TDG' if epsilon[bit] < 0 else 'T', bit) for bit in range(2)]
            actual_word = _inverse(band) + k_word + band
            actual = native_columns(width, actual_word)
            correction = []
            for bit, target in enumerate(permutation):
                exponent = (epsilon[target] - epsilon[bit] * (-1) ** signs[bit]) // 2
                formula = ((negative >> bit & 1) - (negative >> target & 1)
                           + epsilon[bit] * signs[bit])
                self.assertEqual(exponent, formula)
                correction += [('S', target)] * (exponent % 4)
            scalar = roots[-sum(e * s for e, s in zip(epsilon, signs)) % 8]
            expected = scale(native_columns(width, k_word + correction), scalar)
            self.assertEqual(actual, expected)
            actual_inverse = native_columns(width, _inverse(actual_word))
            self.assertEqual(actual_inverse, adjoint(actual))
            self.assertEqual(multiply(actual_inverse, actual), eye(8))

    def test_exact_pauli_gap_and_eighth_root_separation(self):
        # The nonmatching K=H, D=T sends Z to (X+Y)/sqrt(2).
        # Its distance from the nearest signed Pauli has squared Gram
        # exactly (2-sqrt(2)) I, retaining the literal phase and all columns.
        word = [('TDG', 0), ('H', 0), ('T', 0)]
        unitary = native_columns(1, word)
        z, x = _pauli(1, 2), _pauli(1, 1)
        image = multiply(multiply(unitary, z), adjoint(unitary))
        defect = [[entry - other for entry, other in zip(row, target)]
                  for row, target in zip(image, x)]
        gamma_squared = Q2(2, -1)
        self.assertEqual(multiply(adjoint(defect), defect),
                         scale(eye(2), gamma_squared))
        omega = (CO + CI) * Q2(0, Fraction(1, 2))
        roots = [CO]
        for _ in range(7):
            roots.append(roots[-1] * omega)
        # Every nonzero root difference has one of these four exact norms;
        # gamma_squared is the smallest, as 0<sqrt(2)<2.
        squared = [(root - CO).conj() * (root - CO) for root in roots[1:]]
        self.assertEqual(squared,
                         [CO * Q2(2, -1), CO * Q2(2), CO * Q2(2, 1), CO * Q2(4),
                          CO * Q2(2, 1), CO * Q2(2), CO * Q2(2, -1)])

    def test_actual_field_clifford_phase_needs_degree_four_modulo_four(self):
        # A real p=4 field shear on 16 arbitrary carrier wires. Only sparse
        # Pauli labels are propagated; no clean-input or dense-state shortcut.
        p, half = 4, 8
        g_word = [('S', half + 3)]
        r_word = [('H', 0), ('H', 1)]
        phases, xs, zs = [], [], []
        for address in range(16):
            query = [('CX', source, half + target)
                     for source, column in enumerate(_pencil_columns(p, address))
                     for target in range(half) if column >> target & 1]
            self.assertEqual(_normal_image(0, 1 << (half + 1), 0, query)[2], 0)
            word = _inverse(g_word) + query + r_word + query + g_word
            x, z, phase = _normal_image(0, 1 << (half + 1), 0, word)
            phases.append(phase)
            xs.append(x)
            zs.append(z)
        expected = [((address & 5) == 5) ^ ((address & 10) == 10)
                    for address in range(16)]
        self.assertEqual(phases, expected)
        coefficients = _coefficients(phases, 4)
        self.assertEqual({mask: value for mask, value in enumerate(coefficients) if value},
                         {5: 1, 10: 1, 15: 2})
        self.assertEqual(_degree(phases, 4), 4)
        for labels in (xs, zs):
            for bit in range(16):
                self.assertLessEqual(_degree([value >> bit & 1 for value in labels], 2), 2)

    def test_modulo_four_lift_even_division_and_nine_bit_singleton(self):
        # A degree-four Boolean sign can need degree eight when lifted to
        # Z/4. Omitting that carry would give an unjustified stronger theorem.
        sign = [int((x & 15) == 15) ^ int((x & 240) == 240) for x in range(256)]
        self.assertEqual(_degree(sign, 2), 4)
        lifted = _coefficients(sign, 4)
        self.assertEqual({mask: value for mask, value in enumerate(lifted) if value},
                         {15: 1, 240: 1, 255: 2})
        self.assertEqual(_degree(sign, 4), 8)
        even = _coefficients([2 * value for value in sign], 4)
        self.assertTrue(all(value % 2 == 0 for value in even))
        self.assertEqual([value // 2 for value in even], _coefficients(sign, 2))
        for target in (0, 173, 511):
            singleton = [int(x == target) for x in range(512)]
            self.assertEqual(_degree(singleton, 2), 9)
            self.assertEqual(_coefficients(singleton, 2)[511], 1)
            self.assertEqual(_coefficients([2 * value for value in singleton], 4)[511], 2)


if __name__ == '__main__':
    unittest.main()
