"""Finite all-input checks for lookup with a reused dirty indicator pool.

The bounded native fixtures have at most two high guard controls, so they
use literal X/CX/CCX words. They do not implement the general two-borrowed-
helper MCX theorem. Both arbitrary dirty helpers are nevertheless reserved,
disjoint, and included in the complete-space checks. These fixtures audit
the identity whenever their literal width fits, independently of the
stronger sufficient pool used to optimize the asymptotic theorem.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_operator_source_compiler import _basis_action, _word_matrix
    from .test_parallel_dirty_lookup import (
        _classical_router, _expected, _symbolic,
    )
    from .test_t_depth import _inverse, _native_schedule, _router, _word
except ImportError:
    from test_operator_source_compiler import _basis_action, _word_matrix
    from test_parallel_dirty_lookup import _classical_router, _expected, _symbolic
    from test_t_depth import _inverse, _native_schedule, _router, _word


ATOL = 5e-11


def _guarded_indicator(low, high, output, helpers, chunk, *, native=False):
    """Chronological route, high-address equality toggle, actual unroute."""
    wires = low + high + output + helpers
    assert len(helpers) == 2 and len(set(wires)) == len(wires)
    assert len(output) == 1 << len(low)
    assert len(high) <= 2 and 0 <= chunk < 1 << len(high)
    negatives = [('X', wire) for bit, wire in enumerate(high)
                 if not (chunk >> bit & 1)]
    gate = ('X', 'CX', 'CCX')[len(high)]
    toggle = negatives + [(gate, *high, output[0])] + list(reversed(negatives))
    if native:
        router = _router(low, [[wire] for wire in output])
        return router + _native_schedule(toggle) + _inverse(router)
    router = _classical_router(low, [[wire] for wire in output])
    return router + toggle + list(reversed(router))


def _batched_fixture(rows, word_bits, bank_count, capacity, *, native=False):
    bits = rows.bit_length() - 1
    bank_bits = bank_count.bit_length() - 1
    selector_bits = capacity.bit_length() - 1
    assert rows == 1 << bits and bank_count == 1 << bank_bits
    assert capacity == 1 << selector_bits and bank_count * capacity <= rows
    address = list(range(bits))
    output = list(range(bits, bits + word_bits))
    start = bits + word_bits
    banks = [list(range(start + j * word_bits, start + (j + 1) * word_bits))
             for j in range(bank_count)]
    start += bank_count * word_bits
    selector = list(range(start, start + capacity))
    helpers = [start + capacity, start + capacity + 1]
    width = helpers[-1] + 1
    low = address[bank_bits:bank_bits + selector_bits]
    high = address[bank_bits + selector_bits:]
    mask = (1 << word_bits) - 1
    table = [((3 * row) ^ (row >> 1) ^ 1) & mask
             if row < max(1, rows // 2) else 0 for row in range(rows)]
    reverse = _inverse if native else lambda word: list(reversed(word))
    loader = []
    for chunk in range(rows // (bank_count * capacity)):
        linear = [('CX', selector[h], banks[j][bit])
                  for h in range(capacity) for j in range(bank_count)
                  for bit in range(word_bits)
                  if table[(chunk * capacity + h) * bank_count + j] >> bit & 1]
        if native:
            linear = [('C', linear)]
        indicator = _guarded_indicator(low, high, selector, helpers, chunk,
                                       native=native)
        loader += linear + indicator + reverse(linear) + reverse(indicator)
    copy = [('CX', source, target) for source, target in zip(banks[0], output)]
    if native:
        router = _router(address[:bank_bits], banks)
        copy = [('C', copy)]
    else:
        router = _classical_router(address[:bank_bits], banks)
    selected_copy = router + copy + reverse(router)
    query = loader + selected_copy + reverse(loader) + selected_copy
    missing_copy = loader + selected_copy + reverse(loader)
    return width, address, output, table, query, helpers, banks, missing_copy


def _router_t_count(bank_count, word_bits=1):
    # At each level the shared-control CCZ stage has 6*p+(p mod 2) T gates.
    return sum(6 * pairs + pairs % 2
               for level in range(bank_count.bit_length() - 1)
               for pairs in [word_bits * (bank_count >> (level + 1))])


class BatchedDirtyLookupTests(unittest.TestCase):
    def test_native_guarded_indicator_negative_controls_and_actual_inverse(self):
        low, high, selector, helpers = [0], [1, 2], [3, 4], [5, 6]
        for chunk in range(4):
            schedule = _guarded_indicator(low, high, selector, helpers, chunk,
                                           native=True)
            images = [basis ^ (1 << selector[basis & 1])
                      if (basis >> 1) & 3 == chunk else basis
                      for basis in range(128)]
            expected = np.eye(128)[:, np.argsort(images)]
            actual = _word_matrix(7, _word(schedule))
            np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(7, _word(_inverse(schedule))),
                                       expected.conj().T, atol=ATOL, rtol=0)
        # A wrong negative-control pattern moves a different coherent sector.
        wrong = _guarded_indicator(low, high, selector, helpers, 2, native=True)
        basis = (1 << high[0]) | 1 | (1 << helpers[0])
        matrix = _word_matrix(7, _word(wrong))
        self.assertAlmostEqual(abs(matrix[basis, basis]), 1)
        self.assertAlmostEqual(abs(matrix[basis ^ (1 << selector[1]), basis]), 0)
        with self.assertRaises(AssertionError):
            _guarded_indicator(low, high, selector, [high[0], 6], 1)

    def test_completed_query_symbolic_multi_batch_and_all_dirty_return(self):
        cases = ((32, 2, 4, 2), (16, 3, 2, 2), (8, 2, 1, 2),
                 (8, 2, 2, 1), (8, 2, 2, 4), (4, 2, 4, 1),
                 (2, 5, 1, 2), (1, 3, 1, 1))
        for rows, word_bits, banks, capacity in cases:
            with self.subTest(rows=rows, word_bits=word_bits, banks=banks,
                              capacity=capacity):
                width, address, output, table, query, helpers, _, _ = (
                    _batched_fixture(rows, word_bits, banks, capacity))
                initial = [{1 << wire} for wire in range(width)]
                actual = _symbolic(query, initial)
                self.assertEqual(actual, _expected(initial, address, output, table))
                self.assertEqual(_symbolic(list(reversed(query)), actual), initial)
                for helper in helpers:
                    self.assertEqual(actual[helper], initial[helper])
                if address:
                    inactive = [set(value) for value in initial]
                    inactive[address[-1]] = {0}  # Highest address bit equals one.
                    self.assertEqual(_symbolic(query, inactive), inactive)

    def test_completed_native_query_literal_phase_and_two_reserved_helpers(self):
        # One fixture has a nontrivial indicator router; the other has a
        # nontrivial bank router. Both have two high guard controls, and all
        # 512 input columns include arbitrary dirty values on both helpers.
        for banks, capacity in ((1, 2), (2, 1)):
            width, address, output, table, query, helpers, _, _ = (
                _batched_fixture(8, 1, banks, capacity, native=True))
            self.assertEqual(width, 9)
            self.assertEqual(len(helpers), 2)
            images = [basis ^ (table[basis & 7] << output[0])
                      for basis in range(1 << width)]
            expected = np.eye(1 << width)[:, np.argsort(images)]
            np.testing.assert_allclose(_word_matrix(width, _word(query)), expected,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(width, _word(_inverse(query))),
                                       expected.conj().T, atol=ATOL, rtol=0)

    def test_emitted_native_depth_count_width_and_distinct_t_targets(self):
        for bits in range(6):
            rows = 1 << bits
            for bank_bits in range(bits + 1):
                for selector_bits in range(max(0, bits - bank_bits - 2),
                                           bits - bank_bits + 1):
                    for word_bits in (1, 2, 3):
                        banks, capacity = 1 << bank_bits, 1 << selector_bits
                        high_bits = bits - bank_bits - selector_bits
                        chunks = 1 << high_bits
                        width, _, _, _, query, helpers, _, _ = _batched_fixture(
                            rows, word_bits, banks, capacity, native=True)
                        self.assertEqual(width, bits + word_bits
                                         + banks * word_bits + capacity + 2)
                        indicator_depth = 8 * selector_bits + 4 * (high_bits == 2)
                        indicator_count = (2 * _router_t_count(capacity)
                                           + 7 * (high_bits == 2))
                        self.assertEqual(sum(kind == 'T' for kind, _ in query),
                                         4 * chunks * indicator_depth + 16 * bank_bits)
                        actual_count = sum(len(gates) for kind, gates in query if kind == 'T')
                        self.assertEqual(actual_count, 4 * chunks * indicator_count
                                         + 4 * _router_t_count(banks, word_bits))
                        self.assertLessEqual(actual_count,
                                             140 * (rows // banks + banks * word_bits))
                        for kind, gates in query:
                            used = [wire for gate in gates for wire in gate[1:]]
                            self.assertTrue(all(0 <= wire < width for wire in used))
                            self.assertTrue(set(helpers).isdisjoint(used))
                            if kind == 'T':
                                self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                                targets = [gate[1] for gate in gates]
                                self.assertEqual(len(targets), len(set(targets)))

    def test_second_copy_is_needed_for_arbitrary_initial_selected_bank(self):
        _, address, output, table, query, helpers, banks, missing_copy = (
            _batched_fixture(8, 2, 2, 2))
        # Address three selects bank one. Nonzero unknown bank data must be
        # canceled even with dirty indicators and helpers initially nonzero.
        selector_start = banks[-1][-1] + 1
        basis = (3 | (1 << banks[1][0]) | (1 << selector_start)
                 | (1 << helpers[0]) | (1 << helpers[1]))
        expected = basis ^ (table[3] << output[0])
        self.assertEqual(_basis_action(basis, query), expected)
        self.assertNotEqual(_basis_action(basis, missing_copy), expected)
        self.assertEqual(_basis_action(expected, list(reversed(query))), basis)


if __name__ == '__main__':
    unittest.main()
