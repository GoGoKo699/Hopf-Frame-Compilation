"""Literal native and symbolic audits of scratch-free dirty indicators.

Routing to a fixed X target and actually undoing the routing implements an
indicator XOR on arbitrary dirty bits. These small checks protect phases,
address order, work return and T-layer support; they do not prove the
asymptotic complete-frame theorem.
"""
from __future__ import annotations

import unittest
import numpy as np

try:
    from .test_operator_source_compiler import _word_matrix
    from .test_t_depth import _inverse, _router, _word
except ImportError:
    from test_operator_source_compiler import _word_matrix
    from test_t_depth import _inverse, _router, _word


def _indicator(address, output):
    assert len(output) == 1 << len(address)
    assert len(set(address + output)) == len(address + output)
    router = _router(address, [[wire] for wire in output])
    return router + [('C', [('X', output[0])])] + _inverse(router)


def _classical_router(address, banks):
    word = []
    for level, control in enumerate(address):
        stride = 1 << level
        for bank in range(0, len(banks), 2 * stride):
            for left, right in zip(banks[bank], banks[bank + stride]):
                word += [('CX', left, right), ('CCX', control, right, left),
                         ('CX', left, right)]
    return word


def _classical_indicator(address, output):
    router = _classical_router(address, [[wire] for wire in output])
    return router + [('X', output[0])] + list(reversed(router))


# Boolean polynomials over GF(2), with bit masks for square-free monomials.
def _multiply(left, right):
    result = set()
    for a in left:
        for b in right:
            result.symmetric_difference_update({a | b})
    return result


def _symbolic(word, values):
    values = [set(value) for value in values]
    for gate in word:
        if gate[0] == 'X':
            change = {0}
        elif gate[0] == 'CX':
            change = values[gate[1]]
        else:
            assert gate[0] == 'CCX'
            change = _multiply(values[gate[1]], values[gate[2]])
        values[gate[-1]] ^= change
    return values


def _expected(values, address, output, table):
    result = [set(value) for value in values]
    for row, data in enumerate(table):
        indicator = {0}
        for bit, wire in enumerate(address):
            factor = values[wire] if row >> bit & 1 else values[wire] ^ {0}
            indicator = _multiply(indicator, factor)
        for bit, target in enumerate(output):
            if data >> bit & 1:
                result[target] ^= indicator
    return result


def _query_fixture(rows, word_bits, bank_count, *, native=False):
    bits, low_bits = rows.bit_length() - 1, bank_count.bit_length() - 1
    assert rows == 1 << bits and bank_count == 1 << low_bits
    address = list(range(bits))
    output = list(range(bits, bits + word_bits))
    first_bank = bits + word_bits
    banks = [list(range(first_bank + j * word_bits,
                        first_bank + (j + 1) * word_bits))
             for j in range(bank_count)]
    first_selector = first_bank + bank_count * word_bits
    selector = list(range(first_selector, first_selector + rows // bank_count))
    width = first_selector + len(selector)
    mask = (1 << word_bits) - 1
    table = [(((3 * row) ^ (row >> 1) ^ 1) & mask) if row < rows // 2 else 0
             for row in range(rows)]
    linear = [('CX', selector[h], banks[j][bit])
              for h in range(len(selector)) for j in range(bank_count)
              for bit in range(word_bits) if table[h * bank_count + j] >> bit & 1]
    copy = [('CX', source, target) for source, target in zip(banks[0], output)]
    if native:
        indicator = _indicator(address[low_bits:], selector)
        linear, copy = [('C', linear)], [('C', copy)]
        reverse = _inverse
        router = _router(address[:low_bits], banks)
    else:
        indicator = _classical_indicator(address[low_bits:], selector)
        reverse = lambda word: list(reversed(word))
        router = _classical_router(address[:low_bits], banks)
    loader = linear + indicator + reverse(linear) + reverse(indicator)
    query = (loader + router + copy + reverse(router) + reverse(loader)
             + router + copy + reverse(router))
    return width, address, output, table, query


class ParallelDirtyLookupTests(unittest.TestCase):
    def test_indicator_native_phase_and_actual_inverse(self):
        for bits in range(3):
            rows = 1 << bits
            width = bits + rows
            address, output = list(range(bits)), list(range(bits, width))
            schedule = _indicator(address, output)
            images = [basis ^ (1 << (bits + (basis & (rows - 1))))
                      for basis in range(1 << width)]
            expected = np.eye(1 << width)[:, np.argsort(images)]
            actual = _word_matrix(width, _word(schedule))
            np.testing.assert_allclose(actual, expected, atol=5e-11, rtol=0)
            np.testing.assert_allclose(_word_matrix(width, _word(_inverse(schedule))),
                                       expected.conj().T, atol=5e-11, rtol=0)
        # The router is not an involution. At x=3 the reversed orientation
        # flips original Y1, whereas the intended circuit flips Y3.
        router = _router([0, 1], [[2], [3], [4], [5]])
        wrong = _inverse(router) + [('C', [('X', 2)])] + router
        wrong_matrix = _word_matrix(6, _word(wrong))
        self.assertAlmostEqual(abs(wrong_matrix[3 ^ (1 << 3), 3]), 1)
        self.assertAlmostEqual(abs(wrong_matrix[3 ^ (1 << 5), 3]), 0)

    def test_indicator_symbolically_preserves_all_other_dirty_bits(self):
        for bits in range(6):
            rows = 1 << bits
            width = bits + rows
            address, output = list(range(bits)), list(range(bits, width))
            initial = [{1 << wire} for wire in range(width)]
            word = _classical_indicator(address, output)
            expected = _expected(initial, address, output, [1 << x for x in range(rows)])
            actual = _symbolic(word, initial)
            self.assertEqual(actual, expected)
            self.assertEqual(_symbolic(list(reversed(word)), actual), initial)

    def test_indicator_native_depth_and_scratch_free_width(self):
        for bits in range(7):
            rows, width = 1 << bits, bits + (1 << bits)
            schedule = _indicator(list(range(bits)), list(range(bits, width)))
            self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 8 * bits)
            count = sum(len(gates) for kind, gates in schedule if kind == 'T')
            self.assertEqual(count, 0 if bits == 0 else 12 * rows - 10)
            for kind, gates in schedule:
                self.assertTrue(all(0 <= wire < width for gate in gates for wire in gate[1:]))
                if kind == 'T':
                    self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                    targets = [gate[1] for gate in gates]
                    self.assertEqual(len(targets), len(set(targets)))

    def test_completed_query_symbolically_returns_work_and_inactive_rows(self):
        for rows, word_bits, bank_count in ((8, 2, 2), (2, 2, 2), (4, 1, 1)):
            width, address, output, table, query = _query_fixture(rows, word_bits, bank_count)
            initial = [{1 << wire} for wire in range(width)]
            actual = _symbolic(query, initial)
            self.assertEqual(actual, _expected(initial, address, output, table))
            self.assertEqual(_symbolic(list(reversed(query)), actual), initial)
            inactive = [set(value) for value in initial]
            inactive[address[-1]] = {0}
            self.assertEqual(_symbolic(query, inactive), inactive)

    def test_completed_native_query_has_literal_phase_and_return(self):
        width, address, output, table, query = _query_fixture(4, 1, 2, native=True)
        images = [basis ^ (table[basis & 3] << output[0]) for basis in range(1 << width)]
        expected = np.eye(1 << width)[:, np.argsort(images)]
        np.testing.assert_allclose(_word_matrix(width, _word(query)), expected,
                                   atol=5e-11, rtol=0)
        np.testing.assert_allclose(_word_matrix(width, _word(_inverse(query))),
                                   expected.conj().T, atol=5e-11, rtol=0)


if __name__ == '__main__':
    unittest.main()
