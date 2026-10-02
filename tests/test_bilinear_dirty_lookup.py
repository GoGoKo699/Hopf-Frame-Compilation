"""Bounded audits of lookup by a bilinear form and two dirty indicators.

The rank reduction acts on two arbitrary dirty input words, so its left
basis is inverse-transposed rather than the output basis of a linear
shear. Full matrices check literal phases and actual inverses; symbolic
checks retain every dirty input variable. These fixtures are not a general
compiler or an asymptotic proof.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_amortized_dirty_lookup import _normal_form
    from .test_operator_source_compiler import _basis_action, _word_matrix
    from .test_parallel_dirty_lookup import (
        _classical_indicator, _expected, _indicator, _symbolic,
    )
    from .test_t_depth import _ccz_batch, _inverse, _word
except ImportError:
    from test_amortized_dirty_lookup import _normal_form
    from test_operator_source_compiler import _basis_action, _word_matrix
    from test_parallel_dirty_lookup import _classical_indicator, _expected, _indicator, _symbolic
    from test_t_depth import _ccz_batch, _inverse, _word


ATOL = 5e-11


def _bilinear(matrix, left, right, target, *, native=False):
    """z ^= Y^T D X, with literal preservation of arbitrary Y and X."""
    assert len({*left, *right, target}) == len(left) + len(right) + 1
    before, rank = _normal_form(matrix, right, left)
    if not rank:
        return [], 0
    # _normal_form has P D Q=diag(I_rank), with chronological bases P on
    # left and Q^{-1} on right. Transpose each left CNOT in the SAME order:
    # its resulting operator is P^{-T}, the basis required by a bilinear form.
    left_set = set(left)
    before = [('CX', gate[2], gate[1]) if gate[1] in left_set else gate
              for gate in before]
    if native:
        pre = [('C', before)]
        middle = ([('C', [('H', target)])]
                  + _ccz_batch(target, list(zip(left[:rank], right[:rank])))
                  + [('C', [('H', target)])])
        return pre + middle + _inverse(pre), rank
    middle = [('CCX', left[i], right[i], target) for i in range(rank)]
    return before + middle + list(reversed(before)), rank


def _bilinear_fixture(height, width, word_bits, *, native=False, table=None):
    row_bits, column_bits = height.bit_length() - 1, width.bit_length() - 1
    assert height == 1 << row_bits and width == 1 << column_bits
    bits = row_bits + column_bits
    address = list(range(bits))
    output = list(range(bits, bits + word_bits))
    left = list(range(bits + word_bits, bits + word_bits + height))
    right = list(range(left[-1] + 1, left[-1] + 1 + width))
    total_width = right[-1] + 1
    rows = height * width
    if table is None:
        table = [((5 * index) ^ (index >> 1) ^ 1) & ((1 << word_bits) - 1)
                 if index < max(1, rows // 2) else 0 for index in range(rows)]
    assert len(table) == rows and all(0 <= data < 1 << word_bits for data in table)
    bilinear, ranks = [], []
    for bit, target in enumerate(output):
        matrix = [[table[row + height * column] >> bit & 1 for column in range(width)]
                  for row in range(height)]
        block, rank = _bilinear(matrix, left, right, target, native=native)
        bilinear += block
        ranks.append(rank)
    indicator = _indicator if native else _classical_indicator
    row_indicator = indicator(address[:row_bits], left)
    column_indicator = indicator(address[row_bits:], right)
    reverse = _inverse if native else lambda word: list(reversed(word))
    # Chronological four-corner identity. Commuting the two independent
    # indicators removes the redundant calls of the expanded nested echo.
    query = (bilinear + row_indicator + reverse(bilinear) + column_indicator
             + bilinear + reverse(row_indicator) + reverse(bilinear)
             + reverse(column_indicator))
    return total_width, address, output, table, query, ranks, query[len(bilinear):]


class BilinearDirtyLookupTests(unittest.TestCase):
    def test_rectangular_rank_bases_literal_native_phase_and_inverse(self):
        matrices = (
            [[0, 0], [0, 0]],
            [[1, 1, 0], [1, 1, 0]],
            [[1, 1, 0], [1, 0, 1]],
            [[0, 1], [1, 1], [1, 0]],
            [[0, 0, 1]],
            [[1], [0], [1]],
            [[1, 0, 1], [0, 1, 1], [1, 1, 0]],
            [[1, 1, 0], [0, 1, 1], [0, 0, 1]],
        )
        for matrix in matrices:
            height, width = len(matrix), len(matrix[0])
            left, right = list(range(height)), list(range(height, height + width))
            target = height + width
            schedule, rank = _bilinear(matrix, left, right, target, native=True)
            column_images = {sum((sum(matrix[i][j] * (x >> j & 1) for j in range(width)) % 2) << i
                                 for i in range(height)) for x in range(1 << width)}
            self.assertEqual(rank, len(column_images).bit_length() - 1)
            self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 4 if rank else 0)
            self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'),
                             6 * rank + rank % 2)
            images = []
            for basis in range(1 << (target + 1)):
                delta = sum(matrix[i][j] * (basis >> left[i] & 1) * (basis >> right[j] & 1)
                            for i in range(height) for j in range(width)) % 2
                images.append(basis ^ (delta << target))
            expected = np.eye(len(images))[:, np.argsort(images)]
            np.testing.assert_allclose(_word_matrix(target + 1, _word(schedule)), expected,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(target + 1, _word(_inverse(schedule))),
                                       expected.conj().T, atol=ATOL, rtol=0)

    def test_four_corner_query_symbolic_dirty_return_and_inactive_rows(self):
        for height, width, word_bits in ((1, 1, 3), (1, 4, 2), (4, 1, 2),
                                         (2, 4, 3), (4, 2, 2), (4, 4, 2)):
            total, address, output, table, query, _, _ = _bilinear_fixture(height, width, word_bits)
            initial = [{1 << wire} for wire in range(total)]
            actual = _symbolic(query, initial)
            self.assertEqual(actual, _expected(initial, address, output, table))
            self.assertEqual(_symbolic(list(reversed(query)), actual), initial)
            if address:
                inactive = [set(value) for value in initial]
                inactive[address[-1]] = {0}
                self.assertEqual(_symbolic(query, inactive), inactive)

    def test_completed_native_queries_full_input_phase_and_actual_inverse(self):
        # The two-output case has one rank-two and one rank-one matrix.
        # Rectangular cases also exercise the actual two-level router inverse.
        for height, width, word_bits, table in ((2, 2, 2, [1, 0, 2, 3]),
                                               (1, 4, 1, None), (4, 1, 1, None)):
            total, address, output, table, query, _, _ = _bilinear_fixture(
                height, width, word_bits, native=True, table=table)
            self.assertEqual(total, 8)
            images = [basis ^ (table[basis & ((1 << len(address)) - 1)] << output[0])
                      for basis in range(1 << total)]
            expected = np.eye(1 << total)[:, np.argsort(images)]
            np.testing.assert_allclose(_word_matrix(total, _word(query)), expected,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(total, _word(_inverse(query))),
                                       expected.conj().T, atol=ATOL, rtol=0)

    def test_emitted_four_bilinear_two_indicator_resource_ledger(self):
        for row_bits in range(4):
            for column_bits in range(4):
                height, width = 1 << row_bits, 1 << column_bits
                for word_bits in (1, 2, 3):
                    total, _, _, _, query, ranks, _ = _bilinear_fixture(
                        height, width, word_bits, native=True)
                    self.assertEqual(total, row_bits + column_bits + word_bits + height + width)
                    indicator_count = sum(0 if size == 1 else 12 * size - 10
                                          for size in (height, width))
                    count = sum(len(gates) for kind, gates in query if kind == 'T')
                    depth = sum(kind == 'T' for kind, _ in query)
                    self.assertEqual(count, 4 * sum(6 * rank + rank % 2 for rank in ranks)
                                     + 2 * indicator_count)
                    self.assertEqual(depth, 16 * sum(rank > 0 for rank in ranks)
                                     + 16 * (row_bits + column_bits))
                    self.assertLessEqual(count, 28 * word_bits * min(height, width)
                                         + 28 * (height + width - 2))
                    for kind, gates in query:
                        self.assertTrue(all(0 <= wire < total for gate in gates for wire in gate[1:]))
                        if kind == 'T':
                            self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                            targets = [gate[1] for gate in gates]
                            self.assertEqual(len(targets), len(set(targets)))

    def test_wrong_left_basis_and_missing_bilinear_leave_detectable_errors(self):
        matrix, left, right, target = [[1, 1], [1, 0]], [0, 1], [2, 3], 4
        before, rank = _normal_form(matrix, right, left)
        wrong = before + [('CCX', left[i], right[i], target) for i in range(rank)] + list(reversed(before))
        correct, _ = _bilinear(matrix, left, right, target)
        self.assertTrue(any(_basis_action(basis, wrong) != _basis_action(basis, correct)
                            for basis in range(1 << target)))
        _, _, output, _, query, _, missing = _bilinear_fixture(2, 2, 1, table=[1, 0, 0, 0])
        # Address (1,1) selects zero, while unknown Y0=X0=1 contaminates
        # the result if the first corner's bilinear call is omitted.
        basis = 3 | (1 << 3) | (1 << 5)
        self.assertEqual(_basis_action(basis, query), basis)
        self.assertEqual(_basis_action(basis, missing), basis ^ (1 << output[0]))


if __name__ == '__main__':
    unittest.main()
