"""Focused checks for native C3 source structure, not a source search.

The arbitrary-width classification, extension and overlap results are
analytic. These fixtures retain literal phases, arbitrary work columns,
actual inverses and the distinction between independent physical bands
and a dependent commuting phase-gadget list. No returned query is built.
"""
from __future__ import annotations

import itertools
import math
import unittest

import numpy as np

from tests.exact_matrices import CI, CO, CZ, native_columns
from tests.test_one_band_axis_matching import (
    _one_qubit_cliffords, _projective_distance,
)
from tests.test_outer_axis_normalizer import _image, _tableau, _tau
from tests.test_read_only_three_layer_rigidity import H, T, X, Z, _wire
from tests.test_returned_diagonal_echo_degree import _ccz


def _native_phase(axis, orientation=1):
    """Literal T/T-dagger phase on the -1 eigenspace of a signed axis."""
    identity = np.eye(len(axis), dtype=complex)
    return ((identity + axis) / 2
            + np.exp(orientation * 1j * math.pi / 4) * (identity - axis) / 2)


def _band(labels, width, orientations=None):
    orientations = [1] * len(labels) if orientations is None else orientations
    result = np.eye(1 << width, dtype=complex)
    for label, orientation in zip(labels, orientations):
        result = result @ _native_phase(_tau(label, width), orientation)
    return result


def _is_clifford(operator, width):
    try:
        _tableau(operator, width)
        return True
    except ValueError:
        return False


def _is_c3(operator, width):
    # Images of Pauli generators suffice because Cliffords form a group.
    return all(_is_clifford(operator @ _tau(1 << bit, width)
                            @ operator.conj().T, width)
               for bit in range(2 * width))


def _rank(labels):
    pivots = {}
    for value in labels:
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def _pair(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) & 1


def _coefficients(operator, width):
    return np.array([np.trace(_tau(p, width) @ operator) / (1 << width)
                     for p in range(1 << (2 * width))])


def _inverse(schedule):
    inverse_name = {'T': 'TDG', 'TDG': 'T', 'CX': 'CX', 'H': 'H'}
    return [(kind, [(inverse_name[name], *wires)
                    for name, *wires in reversed(gates)])
            for kind, gates in reversed(schedule)]


def _flatten(schedule):
    return [gate for _, gates in schedule for gate in gates]


def _ledger(schedule):
    word = _flatten(schedule)
    t_count = sum(gate[0] in ('T', 'TDG') for gate in word)
    return t_count, len(word) - t_count, sum(kind == 'T' for kind, _ in schedule)


class NativeBandSourceStructure(unittest.TestCase):
    def test_two_native_bands_require_cross_axis_commutation(self):
        width = 2
        # X labels precede Z labels. Each listed individual band is independent
        # and isotropic; the combined union may have dependencies.
        fixtures = [([], [4, 8], [], [1, -1]),
                    ([4], [12], [-1], [1]),
                    ([4, 2], [4, 6], [1, -1], [-1, 1]),
                    ([4], [1], [1], [1]),
                    ([4, 8], [1, 2], [-1, 1], [1, -1])]
        for a, b, signs_a, signs_b in fixtures:
            for axes in (a, b):
                self.assertEqual(_rank(axes), len(axes))
                self.assertTrue(all(_pair(p, q, width) == 0
                                    for p, q in itertools.combinations(axes, 2)))
            commute = all(_pair(p, q, width) == 0 for p in a for q in b)
            source = _band(b, width, signs_b) @ _band(a, width, signs_a)
            self.assertEqual(_is_c3(source, width), commute)
            np.testing.assert_allclose(source.conj().T @ source, np.eye(4), atol=3e-14)
        # The two-band source need not be C3 merely because each band is.
        outside = T @ H @ T @ H
        self.assertFalse(_is_c3(outside, 1))

    def test_half_angle_requires_a_pauli_axis(self):
        for axis in (X, -X, Z):
            rotation = (np.eye(2) + 1j * axis) / math.sqrt(2)
            self.assertTrue(_is_clifford(rotation, 1))
            coefficients = _coefficients(rotation, 1)
            self.assertEqual(sum(abs(coefficients) > 1e-10), 2)
            np.testing.assert_allclose(abs(coefficients[abs(coefficients) > 1e-10]),
                                       1 / math.sqrt(2), atol=2e-14)
        # A traceless Hermitian Clifford involution is not enough.
        np.testing.assert_allclose(H @ H, np.eye(2), atol=2e-14)
        self.assertTrue(_is_clifford(H, 1))
        self.assertFalse(_is_clifford((np.eye(2) + 1j * H) / math.sqrt(2), 1))

    def test_left_and_right_extensions_keep_literal_axis_signs(self):
        source = Z @ H @ T
        self.assertTrue(_is_c3(source, 1))
        np.testing.assert_allclose(source @ Z @ source.conj().T, -X, atol=2e-14)
        np.testing.assert_allclose(source.conj().T @ X @ source, -Z, atol=2e-14)
        right_band = _native_phase(Z, -1)
        left_band = _native_phase(X, -1)
        pulled_band = _native_phase(-Z, -1)
        self.assertTrue(_is_c3(source @ right_band, 1))
        self.assertTrue(_is_c3(left_band @ source, 1))
        np.testing.assert_allclose(left_band @ source, source @ pulled_band, atol=2e-14)
        # Losing the pulled-back sign loses a literal relative phase.
        self.assertGreater(np.linalg.norm(left_band @ source - source @ right_band), 0.5)

    def test_zero_overlap_allows_non_c3_inner_source_and_exact_cancellation(self):
        width = 2  # The second work qubit is arbitrary and inactive.
        a, b, c = [1], [4], [4]
        self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 0)
        self.assertEqual(_pair(a[0], b[0], width), 1)
        first = _wire(H, 0) @ _wire(T, 0) @ _wire(H, 0)
        middle, outer = _wire(T.conj().T, 0), _wire(T, 0)
        self.assertFalse(_is_c3(middle @ first, width))
        source = outer @ middle @ first
        self.assertTrue(_is_c3(source, width))
        np.testing.assert_allclose(source, first, atol=2e-14)
        np.testing.assert_allclose(source.conj().T, first.conj().T, atol=2e-14)
        np.testing.assert_allclose(source.conj().T @ source, np.eye(4), atol=2e-14)
        # Declared source versus its fair simplified native implementation.
        first_word = [('C', [('H', 0)]), ('T', [('T', 0)]), ('C', [('H', 0)])]
        displayed = first_word + [('T', [('TDG', 0)]), ('T', [('T', 0)])]
        self.assertEqual(_ledger(displayed), (3, 2, 3))
        self.assertEqual(_ledger(_inverse(displayed)), (3, 2, 3))
        self.assertEqual(_ledger(first_word), (1, 2, 1))
        self.assertEqual(native_columns(width, _flatten(displayed)),
                         native_columns(width, _flatten(first_word)))
        self.assertEqual(native_columns(width, _flatten(_inverse(displayed))),
                         native_columns(width, _flatten(_inverse(first_word))))

    def test_overlap_orbits_and_relation_budget_do_not_suffice_for_c3(self):
        width = 2
        a, b, c = [12], [4, 8], [3, 15]
        overlap = _rank(a) + _rank(b) - _rank(a + b)
        self.assertEqual(overlap, 1)
        # Centralizer probes activate exactly the even two-coordinate code.
        code = {sum(_pair(axis, p, width) << i for i, axis in enumerate(b))
                for p in range(16) if all(_pair(axis, p, width) == 0 for axis in a)}
        self.assertEqual(code, {0, 3})
        self.assertEqual(_rank(c), len(c))
        self.assertEqual(_pair(*c, width), 0)
        for activation in code:
            images = []
            for axis in c:
                image = axis
                for i, middle_axis in enumerate(b):
                    if activation >> i & 1 and _pair(middle_axis, axis, width):
                        image ^= middle_axis
                images.append(image)
            self.assertEqual(set(images), set(c))
        noncommuting = sum(any(_pair(p, q, width) for q in c) for p in b)
        self.assertEqual(noncommuting, 2 * overlap)
        # The necessary centralizer action can meet its 2t budget without
        # the complete source being C3: other Pauli probes still matter.
        source = _band(c, width) @ _band(b, width) @ _band(a, width)
        self.assertFalse(_is_c3(source, width))

    def test_strict_robust_source_promise_uses_every_pauli_probe(self):
        cliffords = _one_qubit_cliffords()
        threshold = math.sqrt(2) / 4
        for source, expected in [(T, True), (T @ H @ T @ H, False),
                                 (H @ T @ H @ T @ H @ T @ H, False)]:
            distances = []
            for label in range(4):
                image = source @ _tau(label, 1) @ source.conj().T
                distances.append(min(_projective_distance(image, candidate)
                                     for candidate in cliffords))
            if expected:
                self.assertLess(max(distances), 3e-8)
            else:
                self.assertGreaterEqual(max(distances), threshold - 2e-12)
            self.assertEqual(_is_c3(source, 1), expected)

    def test_dependent_toffoli_axes_defeat_left_extension_inference(self):
        width = 3
        dimension = 1 << width
        source = np.zeros((dimension, dimension), dtype=complex)
        for column in range(dimension):
            row = column ^ (((column & 1) & ((column >> 1) & 1)) << 2)
            source[row, column] = 1
        labels = [(4 | controls) << width for controls in range(4)]
        self.assertEqual(_rank(labels), 3)
        self.assertEqual(labels[0] ^ labels[1] ^ labels[2] ^ labels[3], 0)
        dependent_band = _band(labels, width)
        self.assertTrue(_is_c3(source, width))
        self.assertTrue(_is_c3(dependent_band @ source, width))
        # The whole Pauli-image frame permutes the dependent axes with the
        # exact signs from the proof, not merely their binary span.
        for p in range(64):
            alpha, beta, gamma = [(p >> i) & 1 for i in range(3)]
            dressed = source @ _tau(p, width) @ source.conj().T
            for u, v in itertools.product((0, 1), repeat=2):
                image_label, sign = _image(dressed, (4 | u | (v << 1)) << width, width)
                self.assertEqual(image_label, (4 | (u ^ beta) | ((v ^ alpha) << 1)) << width)
                self.assertEqual(sign, gamma ^ (u & alpha) ^ (v & beta) ^ (alpha & beta))
        pulled_axis = source.conj().T @ _tau(labels[0], width) @ source
        coefficients = _coefficients(pulled_axis, width)
        self.assertEqual(set(np.flatnonzero(abs(coefficients) > 1e-10)), set(labels))
        np.testing.assert_allclose(coefficients[labels], [0.5, 0.5, 0.5, -0.5], atol=2e-14)
        np.testing.assert_allclose(pulled_axis @ pulled_axis, np.eye(dimension), atol=2e-14)
        self.assertLess(max(abs(coefficients)), 1)

        # Fully charged native source and actual reversed inverse, exact on
        # every work column. The dependent gadget list is never one band.
        toffoli = [('C', [('H', 2)])] + _ccz(0, 1, 2) + [('C', [('H', 2)])]
        parity_phase = []
        for u, v in itertools.product((0, 1), repeat=2):
            compute = [('CX', control, 2) for control, use in ((0, u), (1, v)) if use]
            parity_phase += [('C', compute), ('T', [('T', 2)]), ('C', list(reversed(compute)))]
        self.assertEqual(_ledger(toffoli), (7, 12, 5))
        self.assertEqual(_ledger(parity_phase), (4, 8, 4))
        schedule = toffoli + parity_phase
        inverse = _inverse(schedule)
        self.assertEqual(_ledger(schedule), (11, 20, 9))
        self.assertEqual(_ledger(inverse), (11, 20, 9))
        expected = [[CZ for _ in range(dimension)] for _ in range(dimension)]
        phases = [CO, CI, -CO, -CI]
        for column in range(dimension):
            row = column ^ (((column & 1) & ((column >> 1) & 1)) << 2)
            phase = sum(((label >> width) & row).bit_count() & 1 for label in labels)
            self.assertEqual(phase % 2, 0)
            expected[row][column] = phases[(phase // 2) % 4]
        self.assertEqual(native_columns(width, _flatten(schedule)), expected)
        expected_inverse = [[expected[column][row].conj() for column in range(dimension)]
                            for row in range(dimension)]
        self.assertEqual(native_columns(width, _flatten(inverse)), expected_inverse)


if __name__ == '__main__':
    unittest.main()
