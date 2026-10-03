"""Exact finite certificates for the rectangular blocked-query allocation.

These checks use the normalized abstract indicator ledger c_I=1. They
audit power-of-two rounding, reservations, and integer consequences of
the analytic inequalities; they do not assert that a literal indicator
has unit gate/width constants. The native query word is unchanged and
is checked in test_blocked_bilinear_lookup.py. No fitted asymptotic
exponent, floating square root, or large circuit simulation is used.
"""
from __future__ import annotations

from math import isqrt
import unittest

try:
    from .test_operator_source_compiler import _basis_action
except ImportError:
    from test_operator_source_compiler import _basis_action


def _allocation(r, m, a, budget):
    """Test-only arithmetic transcription, checked against candidate search."""
    assert r >= 0 and m >= 1 and a >= 1 and budget >= 0
    q, overhead = 1 << r, (a + 2) ** 3
    if q == 1:
        return dict(q=q, p=0, h=1, j=1, k=1, pool=0, helpers=0, direct=True)
    assert budget >= 16 * (overhead + r + 1)
    cap = min(q, isqrt(q * m), budget // (8 * overhead))
    j = 1 << (cap.bit_length() - 1)
    h = min(j, q // j)
    k = q // (h * j)
    p = k.bit_length() - 1
    return dict(q=q, p=p, h=h, j=j, k=k, pool=overhead * (h + j),
                helpers=p + 1, direct=False)


def _cases():
    """Finite edges around every width-cap transition through r=12."""
    for r in range(1, 13):
        q = 1 << r
        precisions = sorted({1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17,
                             q - 1, q, q + 1, 4 * q + 3})
        for a in (1, 2, 3, 7):
            overhead = (a + 2) ** 3
            minimum = 16 * (overhead + r + 1)
            budgets = {minimum, minimum + 1}
            for exponent in range(r + 2):
                transition = 8 * overhead * (1 << exponent)
                budgets.update(b for b in (transition - 1, transition, transition + 1)
                               if b >= minimum)
            for m in precisions:
                for budget in sorted(budgets):
                    yield r, m, a, budget


class RectangularQueryAllocationTests(unittest.TestCase):
    def test_exact_floor_choice_partition_and_simultaneous_reservation(self):
        for r, m, a, budget in _cases():
            f = _allocation(r, m, a, budget)
            q, overhead = 1 << r, (a + 2) ** 3
            # Independent finite feasibility search avoids relying on the
            # floor/square-root implementation of the claimed allocation.
            feasible = [1 << exponent for exponent in range(r + 1)
                        if (1 << (2 * exponent)) <= q * m
                        and 8 * overhead * (1 << exponent) <= budget]
            self.assertEqual(f['j'], max(feasible))
            self.assertEqual(f['h'] * f['j'] * f['k'], q)
            self.assertEqual(sum(x.bit_length() - 1 for x in (f['h'], f['j'], f['k'])), r)
            self.assertTrue(all(x > 0 and x & (x - 1) == 0
                                for x in (f['h'], f['j'], f['k'])))
            self.assertLessEqual(f['h'], f['j'])
            self.assertLessEqual(4 * f['pool'], budget)
            self.assertLessEqual(16 * f['helpers'], budget)
            self.assertLessEqual(16 * (f['pool'] + f['helpers']), 5 * budget)

    def test_integer_block_count_selected_count_and_indicator_bounds(self):
        for r, m, a, budget in _cases():
            f = _allocation(r, m, a, budget)
            q, overhead = 1 << r, (a + 2) ** 3
            self.assertLessEqual(f['k'] * budget ** 2,
                                 max(4 * budget ** 2, 256 * q * overhead ** 2))
            selected = m * q // f['j']
            self.assertEqual(selected, f['k'] * m * f['h'])
            # S <= 2m + 2 sqrt(Qm) + 16 QmP/B, with no irrational
            # approximation: square only a strictly positive excess.
            excess = selected * budget - 2 * m * budget - 16 * q * m * overhead
            if excess > 0:
                self.assertLessEqual(excess ** 2, 4 * budget ** 2 * q * m)
            self.assertLessEqual((f['j'] * overhead) ** 2, q * m * overhead ** 2)
            # The Clifford coordinate-change price is unchanged by shape.
            self.assertEqual(f['k'] * m * f['h'] * f['j'], q * m)

    def test_odd_addresses_precision_cap_and_width_transition(self):
        # An odd address does not force a high-block bit after rectangular
        # allocation: the extra output precision permits H=4,J=8,K=1.
        f = _allocation(r=5, m=2, a=1, budget=10_000)
        self.assertEqual((f['h'], f['j'], f['k']), (4, 8, 1))
        # When Q<m, cap J at Q, giving H=1 and no traversal selectors.
        f = _allocation(r=3, m=19, a=1, budget=10_000)
        self.assertEqual((f['h'], f['j'], f['k'], f['p']), (1, 8, 1, 0))
        # Immediately below an exact power-of-two width threshold, the
        # larger indicator does not fit its declared one-eighth allowance.
        overhead, transition = 27, 8 * 27 * 8
        below = _allocation(r=12, m=1, a=1, budget=transition - 1)
        at = _allocation(r=12, m=1, a=1, budget=transition)
        self.assertEqual((below['j'], at['j']), (4, 8))
        self.assertGreater(8 * overhead * at['j'], transition - 1)
        # Nonintegral sqrt(Qm) is bounded exactly before power-two rounding.
        f = _allocation(r=5, m=7, a=1, budget=10_000)
        self.assertEqual(f['j'], 8)
        self.assertGreater((2 * f['j']) ** 2, f['q'] * 7)

    def test_single_row_uses_only_literal_x_gates_and_no_helpers(self):
        for m in range(1, 7):
            for a in (1, 7):
                f = _allocation(r=0, m=m, a=a, budget=0)
                self.assertTrue(f['direct'])
                self.assertEqual((f['pool'], f['helpers'], f['p']), (0, 0, 0))
            for row in range(1 << m):
                word = [('X', bit) for bit in range(m) if (row >> bit) & 1]
                self.assertTrue(all(gate[1] < m for gate in word))
                for output in range(1 << m):
                    self.assertEqual(_basis_action(output, word), output ^ row)
                    self.assertEqual(_basis_action(output ^ row, list(reversed(word))), output)

    def test_balanced_square_has_exact_sqrt_precision_loss_in_selected_cost(self):
        # This is a negative control for the selected-middle upper ledger,
        # not a lower bound for arbitrary lookup circuits or total T-count.
        for s in range(2, 7):
            for t in range(1, s + 1):
                q, m, overhead = 1 << (2 * s), 1 << (2 * t), 27
                budget = max(16 * (overhead + 2 * s + 1),
                             8 * overhead * (1 << (s + t)))
                f = _allocation(2 * s, m, 1, budget)
                self.assertEqual((f['h'], f['j'], f['k']),
                                 (1 << (s - t), 1 << (s + t), 1))
                rectangular_selected = m * q // f['j']
                balanced_selected = m * (1 << s)
                self.assertEqual(rectangular_selected ** 2, q * m)
                self.assertEqual(balanced_selected, rectangular_selected * (1 << t))
                self.assertEqual((balanced_selected // rectangular_selected) ** 2, m)


if __name__ == '__main__':
    unittest.main()
