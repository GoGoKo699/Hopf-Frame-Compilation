"""Bounded integer certificates for the QBP compiler cost comparison.

These checks exercise finite precision sums, disjoint workspace accounting,
and literal eligibility thresholds. They do not prove asymptotic complexity,
optimality of a gradient algorithm, or any native statevector identity.
"""
from __future__ import annotations

from math import isqrt
import unittest


def _ceil_div(numerator: int, denominator: int) -> int:
    return (numerator + denominator - 1) // denominator


def _ceil_sqrt(value: int) -> int:
    root = isqrt(value)
    return root + (root * root != value)


def _bank_count(rows: int, word_bits: int, bank_pool: int) -> int:
    """The documented power-of-two choice, evaluated without floating point."""
    # floor(sqrt(rows/word_bits)) = isqrt(rows//word_bits).
    cap = max(1, min(rows, isqrt(rows // word_bits), bank_pool // word_bits))
    return 1 << (cap.bit_length() - 1)


class QBPCostComparisonTests(unittest.TestCase):
    def test_exact_weighted_precision_sum(self):
        for n in (*range(1, 33), 48, 64, 96, 128):
            dimension = 1 << n
            for precision in {6, 7, max(6, n), max(6, n * n), dimension + 6}:
                with self.subTest(n=n, precision=precision):
                    weighted = sum((1 << depth) * (precision + n - depth)
                                   for depth in range(n))
                    self.assertEqual(weighted,
                                     (precision + 2) * dimension - precision - n - 2)

    def test_banked_state_workspace_and_query_bound(self):
        saw_small_table = saw_bank_limited = saw_square_root_limited = False
        for n in (*range(6, 25), 32, 48, 64, 96, 128):
            selectors = n - 6
            for precision in {n, n + 1, 2 * n, n * n, 1 << n}:
                core = precision + 11
                base = precision + n + 7
                predicate_helper = core + selectors
                borrowed_signal = predicate_helper + 1
                self.assertEqual(borrowed_signal + 1, base)
                self.assertEqual(base, core + selectors + 2)
                self.assertLessEqual(core, base)
                self.assertGreaterEqual(predicate_helper, core)
                self.assertGreater(borrowed_signal, predicate_helper)
                budgets = {2 * base, 2 * base + 1, 3 * base, 7 * base,
                           max(2 * base, 1 << n), max(2 * base, 1 << (2 * n))}
                address_widths = {0, selectors // 2, selectors}
                if selectors:
                    address_widths.add(1)
                for dirty in budgets:
                    bank_pool = dirty - base
                    self.assertGreaterEqual(2 * bank_pool, dirty)
                    self.assertGreaterEqual(bank_pool, core)
                    for address in address_widths:
                        rows = 1 << address
                        banks = _bank_count(rows, core, bank_pool)
                        with self.subTest(n=n, precision=precision, dirty=dirty, address=address):
                            self.assertEqual(banks & (banks - 1), 0)
                            self.assertGreaterEqual(banks, 1)
                            self.assertLessEqual(banks, rows)
                            self.assertLessEqual(banks * core, bank_pool)
                            high_selectors = address - (banks.bit_length() - 1)
                            self.assertGreaterEqual(high_selectors, 0)
                            self.assertLessEqual(high_selectors, selectors)
                            # Bank interval [base, base+banks*core) is disjoint
                            # from the occupied core, selectors, helper, signal.
                            self.assertLess(borrowed_signal, base)
                            self.assertLessEqual(base + banks * core, dirty)
                            query = rows // banks + banks * core
                            bound = 4 * (core + _ceil_sqrt(rows * core)
                                         + _ceil_div(rows * core, dirty))
                            self.assertLessEqual(query, bound)
                        if rows < core:
                            saw_small_table = True
                            self.assertEqual(banks, 1)
                        elif bank_pool * bank_pool < rows * core:
                            saw_bank_limited = True
                        else:
                            saw_square_root_limited = True
        self.assertTrue(saw_small_table)
        self.assertTrue(saw_bank_limited)
        self.assertTrue(saw_square_root_limited)

    def test_old_complex_bank_one_wire_gap_at_equal_precision(self):
        for n in (*range(1, 33), 64, 128):
            for precision in {max(6, n), max(6, n + 1), max(6, n * n)}:
                with self.subTest(n=n, precision=precision):
                    # Here P=K. The old compiler needs one clean flag and
                    # may legally borrow the second available clean flag.
                    state_base = precision + n + 7
                    old_base = precision + n + 8
                    self.assertEqual(state_base + 1, old_base)
                    state_bank_minimum = 2 * state_base
                    old_effective_dirty = state_bank_minimum + 1
                    self.assertEqual(2 * old_base - old_effective_dirty, 1)
                    self.assertLess(old_effective_dirty, 2 * old_base)
                    self.assertEqual(old_effective_dirty + 1, 2 * old_base)
                    # Reallocation preserves the same physical total.
                    self.assertEqual(2 + state_bank_minimum, 1 + old_effective_dirty)
            if n > 6:
                for precision in {6, n - 1}:
                    # When P=n>K, the dimension floor supplies the gap.
                    state_bank_minimum = 2 * (n + n + 7)
                    self.assertGreaterEqual(state_bank_minimum + 1,
                                            2 * (precision + n + 8))


if __name__ == "__main__":
    unittest.main()
