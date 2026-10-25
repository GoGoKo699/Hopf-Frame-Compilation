"""Exact dirty cubic-phase assembly from the retained R91 derivative.

The compiler is a paid component; no complete lookup source is inferred.
All native words retain arbitrary helper inputs and their actual inverses.
"""
from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import CQ2, Q2, CI, CO, CZ, adjoint, native_columns
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _gray_phase_batch(controls, helpers, signs, arity):
    """Closed parity-cube walk, with every phase on a private dirty wire."""
    assert arity in (2, 3)
    assert len(controls) == len(helpers) == len(signs)
    assert all(len(term) == arity for term in controls)
    assert len(set(helpers)) == len(helpers)
    assert all(sign in (-1, 1) for sign in signs)
    assert not set(helpers).intersection(wire for term in controls for wire in term)
    word = []
    old = 0
    for step in range(1 << arity):
        vertex = step ^ (step >> 1)
        difference = old ^ vertex
        if difference:
            index = difference.bit_length() - 1
            word += [('CX', term[index], helper)
                     for term, helper in zip(controls, helpers)]
        coefficient = (-1) ** vertex.bit_count()
        word += [('T' if coefficient * sign == 1 else 'TDG', helper)
                 for helper, sign in zip(helpers, signs)]
        old = vertex
    # Reflected binary Gray code ends at the highest coordinate.
    word += [('CX', term[arity - 1], helper)
             for term, helper in zip(controls, helpers)]
    return word


def _relative_phase_batch(pairs, helpers, signs):
    """Product F(u,v;h)^sign, with four private-target physical bands."""
    return _gray_phase_batch(pairs, helpers, signs, 2)


def _dirty_cubic_phase(linear, quadratic, cubic):
    """Return (native word, total width) for uncollected weighted terms.

    Quadratics are (u,v,sign) with sign +/-1 and weight 2. Cubics are
    (u,v,w) with weight 4. Linear coefficients are arbitrary integers.
    """
    n, q, c = len(linear), len(quadratic), len(cubic)
    assert all(len({u, v}) == 2 and 0 <= u < n and 0 <= v < n
               and sign in (-1, 1) for u, v, sign in quadratic)
    assert all(len(set(term)) == 3 and all(0 <= wire < n for wire in term)
               for term in cubic)
    h = list(range(n, n + q))
    triples = [(u, v, helper) for (u, v, _), helper in zip(quadratic, h)]
    triples += list(cubic)
    k = list(range(n + q, n + 2 * q + c))

    # Keep all Clifford linear phases before the one linear T band.
    word, linear_t = [], []
    for wire, coefficient in enumerate(linear):
        coefficient %= 8
        if coefficient % 2:
            linear_t.append(('T', wire))
            coefficient -= 1
        word += [('S', wire)] * (coefficient // 2)
    word += linear_t
    word += _relative_phase_batch([(u, v) for u, v, _ in quadratic],
                                  h, [-sign for _, _, sign in quadratic])
    word += _gray_phase_batch(triples, k, [1] * len(k), 3)
    return word, n + 2 * q + c


def _parity_edge_layer(p, parity, address_signed=True):
    """Paid +pi/4 unsigned or -Z_x*pi/4 signed coefficient edge layer."""
    assert p >= 1 and parity in (0, 1)
    size = 1 << p
    linear = [0] * (p + size)
    quadratic, cubic = [], []
    for subset in range(size):
        if subset.bit_count() % 2 != parity:
            continue
        for bit in range(p):
            if (subset >> bit) & 1:
                continue
            a, b = p + subset, p + (subset | (1 << bit))
            if address_signed:
                linear[b] -= 1
                quadratic += [(bit, b, 1), (a, b, 1)]
                cubic.append((bit, a, b))
            else:
                linear[b] += 1
                quadratic.append((a, b, -1))
    word, width = _dirty_cubic_phase(linear, quadratic, cubic)
    wall = [('H', p + subset) for subset in range(size)
            if subset.bit_count() % 2 == parity]
    return wall + word + list(reversed(wall)), width


def _phase_exponent(value, linear, quadratic, cubic):
    bit = lambda wire: (value >> wire) & 1
    return (sum(coefficient * bit(wire) for wire, coefficient in enumerate(linear))
            + 2 * sum(sign * bit(u) * bit(v) for u, v, sign in quadratic)
            + 4 * sum(bit(u) * bit(v) * bit(w) for u, v, w in cubic)) % 8


_SQRT_HALF = Q2(0, Fraction(1, 2))
_ROOT = CQ2(_SQRT_HALF, _SQRT_HALF)
_ROOTS = [CO]
for _ in range(7):
    _ROOTS.append(_ROOTS[-1] * _ROOT)


class DirtyCubicPhaseCompilerTests(unittest.TestCase):
    def _check_complete_diagonal(self, linear, quadratic, cubic):
        word, width = _dirty_cubic_phase(linear, quadratic, cubic)
        matrix = native_columns(width, word)
        for column in range(1 << width):
            phase = _ROOTS[_phase_exponent(column, linear, quadratic, cubic)]
            for row in range(1 << width):
                self.assertEqual(matrix[row][column], phase if row == column else CZ)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(matrix))
        return word, width, matrix

    def test_relative_primitive_all_columns_and_actual_inverse(self):
        word = _relative_phase_batch([(0, 1)], [2], [1])
        actual = native_columns(3, word)
        for column in range(8):
            u, v, h = [(column >> wire) & 1 for wire in range(3)]
            expected = _ROOTS[((4 * h - 2) * u * v) % 8]
            for row in range(8):
                self.assertEqual(actual[row][column], expected if row == column else CZ)
        self.assertEqual(Counter(gate[0] for gate in word),
                         Counter({'CX': 4, 'T': 2, 'TDG': 2}))
        self.assertEqual(native_columns(3, _inverse(word)), adjoint(actual))

    def test_mixed_polynomial_all_dirty_columns_and_correlated_helpers(self):
        linear, quadratic, cubic = (5, -1, 2), [(0, 1, -1)], [(0, 1, 2)]
        _, width, actual = self._check_complete_diagonal(linear, quadratic, cubic)
        # Entangle a logical difference with every helper and retain relative i.
        left, right = 1, 7 | (((1 << (width - 3)) - 1) << 3)
        amplitudes = [CO * _SQRT_HALF, CI * _SQRT_HALF]
        expected = {
            column: amplitude * _ROOTS[_phase_exponent(column, linear, quadratic, cubic)]
            for column, amplitude in zip((left, right), amplitudes)
        }
        for row in range(1 << width):
            output = actual[row][left] * amplitudes[0] + actual[row][right] * amplitudes[1]
            self.assertEqual(output, expected.get(row, CZ))

    def test_signed_quadratics_share_controls_and_keep_even_occurrences(self):
        self._check_complete_diagonal((0, 0, 0),
                                      [(0, 1, 1), (0, 2, -1)], [])
        # Two equal CS occurrences are a CZ; their helper corrections must return.
        self._check_complete_diagonal((0, 0), [(0, 1, 1), (0, 1, 1)], [])

    def test_shared_cubic_controls_all_columns(self):
        self._check_complete_diagonal((0, 0, 0, 0), [],
                                      [(0, 1, 2), (0, 1, 3)])

    def test_uniform_native_ledgers_and_private_band_targets(self):
        for n in (3, 4, 9):
            linear = [7] * n
            quadratic = [(0, wire, (-1) ** wire) for wire in range(1, n)]
            cubic = [(0, 1, wire) for wire in range(2, n)]
            word, width = _dirty_cubic_phase(linear, quadratic, cubic)
            q, c = len(quadratic), len(cubic)
            self.assertEqual(width, n + 2 * q + c)
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                self.assertLessEqual(counts['T'] + counts['TDG'], n + 12 * q + 8 * c)
                self.assertLessEqual(len(native), 4 * n + 24 * q + 16 * c)
                bands = _bands(native)
                self.assertLessEqual(len(bands), 13)
                for band in bands:
                    self.assertEqual(len({gate[1] for gate in band}), len(band))

    def test_parity_edge_ledger_scales_without_initialized_helpers(self):
        for p in (2, 3, 5):
            size, edges = 1 << p, p * (1 << p) // 4
            for parity in (0, 1):
                for signed in (False, True):
                    word, width = _parity_edge_layer(p, parity, signed)
                    self.assertLessEqual(width - p - size, 5 * edges)
                    for native in (word, _inverse(word)):
                        counts = Counter(gate[0] for gate in native)
                        self.assertLessEqual(counts['T'] + counts['TDG'], size + 32 * edges)
                        self.assertLessEqual(len(native), 6 * size + 64 * edges)
                        bands = _bands(native)
                        self.assertLessEqual(len(bands), 13)
                        for band in bands:
                            self.assertEqual(len({gate[1] for gate in band}), len(band))


if __name__ == '__main__':
    unittest.main()
