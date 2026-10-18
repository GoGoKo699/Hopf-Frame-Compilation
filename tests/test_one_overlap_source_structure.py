"""Bounded checks for the analytic overlap-one source theorem.

These are full-operator fixtures and finite checks of its local identities,
not an exhaustive source search or a returned-query construction.
"""
from __future__ import annotations

import itertools
import math
import unittest

import numpy as np

from tests.exact_matrices import adjoint, native_columns
from tests.test_native_band_source_structure import (
    _band, _coefficients, _flatten, _inverse, _is_c3, _is_clifford,
    _ledger, _rank,
)
from tests.test_outer_axis_normalizer import _tau
from tests.test_read_only_three_layer_rigidity import H


def _word_inverse(word):
    result = []
    for name, *wires in reversed(word):
        if name == 'S':
            result.extend([('S', *wires)] * 3)
        else:
            result.append(({'T': 'TDG', 'TDG': 'T', 'CX': 'CX', 'H': 'H'}[name],
                           *wires))
    return result


def _support(operator, width):
    return set(np.flatnonzero(abs(_coefficients(operator, width)) > 1e-10))


class OneOverlapSourceStructure(unittest.TestCase):
    def test_paired_half_angle_and_necessary_hypotheses(self):
        identity = np.eye(4)
        u = _band([3, 15], 2)  # Independent commuting XX and YY axes.
        for axes, expected_size in [([_tau(4, 2), _tau(8, 2)], 1),
                                    ([u @ _tau(p, 2) @ u.conj().T
                                      for p in (4, 8)], 4)]:
            h1, h2 = axes
            np.testing.assert_allclose(h1 @ h2, h2 @ h1, atol=3e-14)
            for h in (h1, h2, h1 @ h2):
                self.assertAlmostEqual(abs(np.trace(h)), 0)
                self.assertTrue(_is_clifford(h, 2))
            for e, f in itertools.product((-1, 1), repeat=2):
                k = (identity + 1j * e * h1) @ (identity + 1j * f * h2) / 2
                self.assertTrue(_is_clifford(k, 2))
                self.assertEqual(len(_support(k, 2)), 4)
                np.testing.assert_allclose(np.trace(k) / 4, 0.5, atol=3e-14)
                self.assertEqual(_support(k, 2), _support(k.conj().T, 2))
            self.assertEqual(len(_support(h1, 2)), expected_size)
            self.assertEqual(len(_support(h1 @ h2, 2)), 1)
            if expected_size == 4:
                self.assertEqual(_support(h1, 2), _support(h2, 2))
                self.assertEqual(len(_support(h1 + h2, 2)), 2)
        # Commuting, orthogonal Hermitian Cliffords alone do not suffice.
        h1 = np.kron(np.eye(2), H)
        h2 = _tau(8, 2) @ h1
        self.assertEqual(len(_support(h1, 2)), 2)
        self.assertEqual(len(_support(h2, 2)), 2)
        self.assertFalse(_is_clifford((identity + 1j*h1) @ (identity + 1j*h2)/2, 2))

    def test_weight_one_duplicate_and_composite_obstruction(self):
        a, b, c = [4], [4], [1]
        self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 1)
        source = _band(c, 2) @ _band(b, 2, [-1]) @ _band(a, 2)
        np.testing.assert_allclose(source, _band(c, 2), atol=3e-14)
        self.assertTrue(_is_c3(source, 2))
        # A composite weight-one B relation cannot retain nontrivial C spread.
        a, b, c = [4, 8], [12], [1]
        self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 1)
        self.assertFalse(_is_c3(_band(c, 2) @ _band(b, 2) @ _band(a, 2), 2))

    def test_dependent_four_axis_exception_preserves_another_frame(self):
        width = 3
        ccz = np.diag([1] * 7 + [-1]).astype(complex)
        axes = [4, 12, 20, 28]  # Xt, Xt Za, Xt Zb, Xt Za Zb.
        self.assertEqual(_rank(axes), 3)
        e4 = _band(axes, width)
        source = ccz @ e4
        self.assertTrue(_is_c3(source, width))
        images = [ccz @ _tau(a, width) @ ccz for a in axes]
        self.assertTrue(all(len(_support(h, width)) == 4 for h in images))
        self.assertTrue(all(_support(h, width) == _support(images[0], width)
                            for h in images))
        for label in (8, 16):
            np.testing.assert_allclose(source @ _tau(label, width) @ source.conj().T,
                                       _tau(label, width), atol=4e-14)
        # Yt completes a commuting Pauli frame; Xt itself does not.
        np.testing.assert_allclose(source @ _tau(36, width) @ source.conj().T,
                                   -_tau(60, width), atol=4e-14)
        self.assertEqual(len(_support(source @ _tau(4, width) @ source.conj().T,
                                     width)), 4)

    def test_order_three_core_cannot_couple_to_spectators(self):
        count = 0
        for spectator_width in range(3):
            k = spectator_width
            dimension = 2 + 2*k
            identity = np.eye(dimension, dtype=np.uint8)
            omega = np.zeros((dimension, dimension), dtype=np.uint8)
            omega[0, 1] = omega[1, 0] = 1
            for j in range(k):
                omega[2+j, 2+k+j] = omega[2+k+j, 2+j] = 1
            positions = [(i, j) for i in range(k) for j in range(i, k)]
            for coefficients in itertools.product((0, 1), repeat=len(positions)):
                d = identity.copy()
                for (i, j), bit in zip(positions, coefficients):
                    d[2+k+i, 2+j] = d[2+k+j, 2+i] = bit
                for target in ([[0, 1], [1, 1]], [[1, 1], [1, 0]]):
                    m = identity.copy()
                    m[:2, :2] = target
                    for label in range(1 << dimension):
                        v = np.array([(label >> j) & 1 for j in range(dimension)],
                                     dtype=np.uint8)
                        transvection = identity ^ np.outer(v, v @ omega % 2)
                        action = m @ transvection @ d % 2
                        involution = np.array_equal(action @ action % 2, identity)
                        self.assertEqual(involution, bool(label & 3) and label >> 2 == 0)
                        count += 1
        self.assertEqual(count, 1096)

    def test_signed_two_qubit_cores_and_spectator_obstruction(self):
        count = 0
        for a in ([2, 6], [10, 14]):
            b, c = [8, 12], [2, 6]
            self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 1)
            for signs in itertools.product((-1, 1), repeat=6):
                source = (_band(c, 2, signs[4:]) @ _band(b, 2, signs[2:4])
                          @ _band(a, 2, signs[:2]))
                self.assertTrue(_is_c3(source, 2))
                count += 1
        self.assertEqual(count, 128)
        a, b, c = [34, 42], [16, 24], [2, 10]
        self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 1)
        self.assertFalse(_is_c3(_band(c, 3) @ _band(b, 3) @ _band(a, 3), 3))

    def test_literal_native_source_actual_inverse_and_fair_baseline(self):
        b = [('C', [('CX', 1, 0)]), ('T', [('T', 0), ('T', 1)]),
             ('C', [('CX', 1, 0)])]
        c = [('C', [('H', 1)])] + b + [('C', [('H', 1)])]
        schedule = c + b + c
        self.assertEqual(_ledger(schedule), (6, 10, 3))
        self.assertEqual(_ledger(_inverse(schedule)), (6, 10, 3))
        exact = native_columns(2, _flatten(schedule))
        self.assertEqual(native_columns(2, _flatten(_inverse(schedule))), adjoint(exact))
        matrix = np.array([[v.display() for v in row] for row in exact])
        omega = np.exp(1j * math.pi / 4)
        expected = np.zeros((4, 4), dtype=complex)
        expected[np.ix_([0, 2], [0, 2])] = omega * H
        expected[np.ix_([1, 3], [1, 3])] = omega**3 * np.eye(2)
        np.testing.assert_allclose(matrix, expected, atol=3e-14)
        self.assertTrue(_is_c3(matrix, 2))
        self.assertFalse(_is_c3(_band([2, 6], 2) @ _band([8, 12], 2), 2))
        self.assertFalse(_is_c3(_band([8, 12], 2) @ _band([2, 6], 2), 2))
        # One physical Y-axis pair removes the scalar and second T layer.
        middle = [('H', 1), ('CX', 0, 1), ('H', 1), ('S', 1), ('S', 1)]
        y_pair = ([('C', [('S', 1)]*3 + [('H', 1)])] + b
                  + [('C', [('H', 1), ('S', 1)])])
        baseline_schedule = [('C', middle)] + y_pair + [('C', [('S', 0)])]
        baseline = _flatten(baseline_schedule)
        self.assertEqual(_ledger(baseline_schedule), (2, 14, 1))
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in baseline), 2)
        self.assertEqual(native_columns(2, baseline), exact)
        self.assertEqual(native_columns(2, _word_inverse(baseline)), adjoint(exact))
        # The inverse keeps the same coordinates and has its own Clifford cost.
        inverse_pair = [(kind, [('TDG', *gate[1:]) if gate[0] == 'T' else gate
                                for gate in gates]) for kind, gates in y_pair]
        baseline_inverse = [('C', [('S', 0)]*3)] + inverse_pair + [('C', middle)]
        self.assertEqual(_ledger(baseline_inverse), (2, 16, 1))
        self.assertEqual(native_columns(2, _flatten(baseline_inverse)), adjoint(exact))


if __name__ == '__main__':
    unittest.main()
