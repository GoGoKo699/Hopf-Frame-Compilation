"""Finite native audits of a two-pass dirty traversal and amortized lookup.

Rectangular controlled shears are compiled by binary row/column elimination,
with one shared-control four-T-layer stage per nonzero rank. The traversal
uses one arbitrary dirty selector per high-address bit and its literal
inverse. No initialized indicator, native table oracle, or general API is
assumed. Small full matrices and Boolean-polynomial checks certify these
fixtures; the asymptotic result is proved separately.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_batched_dirty_lookup import _router_t_count
    from .test_operator_source_compiler import _basis_action, _word_matrix
    from .test_parallel_dirty_lookup import (
        _classical_indicator, _classical_router, _expected, _indicator,
        _multiply, _symbolic,
    )
    from .test_t_depth import (
        _ccz_batch, _inverse, _native_schedule, _router, _word,
    )
except ImportError:
    from test_batched_dirty_lookup import _router_t_count
    from test_operator_source_compiler import _basis_action, _word_matrix
    from test_parallel_dirty_lookup import (
        _classical_indicator, _classical_router, _expected, _indicator,
        _multiply, _symbolic,
    )
    from test_t_depth import _ccz_batch, _inverse, _native_schedule, _router, _word


ATOL = 5e-11


def _swap(left, right):
    return [('CX', left, right), ('CX', right, left), ('CX', left, right)]


def _normal_form(matrix, inputs, outputs):
    """Return pre-conjugation gates and rank for D=E A F.

    Chronological column-operation gates implement F^{-1} on the inputs;
    chronological row-operation gates implement E on the outputs. Undoing
    both after the diagonal shear gives the original rectangular matrix A.
    """
    work = [list(row) for row in matrix]
    height, width = len(outputs), len(inputs)
    assert len(work) == height and all(len(row) == width for row in work)
    assert all(bit in (0, 1) for row in work for bit in row)
    row_word, column_word, rank = [], [], 0
    while rank < min(height, width):
        pivot = next(((i, j) for i in range(rank, height)
                      for j in range(rank, width) if work[i][j]), None)
        if pivot is None:
            break
        i, j = pivot
        if i != rank:
            work[i], work[rank] = work[rank], work[i]
            row_word += _swap(outputs[i], outputs[rank])
        if j != rank:
            for row in work:
                row[j], row[rank] = row[rank], row[j]
            column_word += _swap(inputs[j], inputs[rank])
        for i in range(height):
            if i != rank and work[i][rank]:
                work[i] = [left ^ right for left, right in zip(work[i], work[rank])]
                row_word.append(('CX', outputs[rank], outputs[i]))
        for j in range(width):
            if j != rank and work[rank][j]:
                for row in work:
                    row[j] ^= row[rank]
                # Adding column rank into column j acts on input coordinate
                # rank, controlled by input j, in the right matrix factor.
                column_word.append(('CX', inputs[j], inputs[rank]))
        rank += 1
    assert work == [[int(i == j and i < rank) for j in range(width)]
                    for i in range(height)]
    return column_word + row_word, rank


def _controlled_shear(matrix, control, inputs, outputs, *, native=False):
    assert len({control, *inputs, *outputs}) == 1 + len(inputs) + len(outputs)
    before, rank = _normal_form(matrix, inputs, outputs)
    if not rank:
        return [], 0
    if native:
        hadamards = [('H', outputs[i]) for i in range(rank)]
        middle = ([('C', hadamards)]
                  + _ccz_batch(control, list(zip(inputs[:rank], outputs[:rank])))
                  + [('C', hadamards)])
        pre = [('C', before)]
        return pre + middle + _inverse(pre), rank
    middle = [('CCX', control, inputs[i], outputs[i]) for i in range(rank)]
    return before + middle + list(reversed(before)), rank


def _selected_shear(address, inputs, outputs, selectors, matrices, *, native=False):
    """First dirty DFS followed by the actual inverse of root-disabled DFS."""
    wires = address + inputs + outputs + selectors
    assert len(wires) == len(set(wires)) and len(selectors) == len(address)
    assert len(matrices) == 1 << len(address)
    if not address:
        linear = [('CX', inputs[j], outputs[i])
                  for i, row in enumerate(matrices[0]) for j, bit in enumerate(row) if bit]
        word = [('C', linear)] if native else linear
        return word, [_normal_form(matrices[0], inputs, outputs)[1]], word

    leaves = [_controlled_shear(matrix, selectors[-1], inputs, outputs, native=native)
              for matrix in matrices]
    reverse = _inverse if native else lambda word: list(reversed(word))

    def traversal(root_enabled):
        def descend(level, chunk):
            if level == len(address):
                return leaves[chunk][0]
            result = []
            for value in (0, 1):
                wire, target = address[level], selectors[level]
                negatives = [('X', wire)] if value == 0 else []
                toggle = ('CX', wire, target) if level == 0 else (
                    'CCX', selectors[level - 1], wire, target)
                update = negatives + [toggle] + list(reversed(negatives))
                if level == 0 and not root_enabled:
                    update = []
                if native:
                    update = _native_schedule(update)
                result += (update + descend(level + 1, chunk | (value << level))
                           + reverse(update))
            return result
        return descend(0, 0)

    first, second = traversal(True), traversal(False)
    return first + reverse(second), [rank for _, rank in leaves], first


def _expected_shear(initial, address, inputs, outputs, matrices):
    expected = [set(value) for value in initial]
    for chunk, matrix in enumerate(matrices):
        selected = {0}
        for bit, wire in enumerate(address):
            literal = initial[wire] if chunk >> bit & 1 else initial[wire] ^ {0}
            selected = _multiply(selected, literal)
        for i, row in enumerate(matrix):
            for j, bit in enumerate(row):
                if bit:
                    expected[outputs[i]] ^= _multiply(selected, initial[inputs[j]])
    return expected


def _amortized_fixture(rows, word_bits, bank_count, capacity, *, native=False):
    bits, bank_bits, selector_bits = (rows.bit_length() - 1,
                                     bank_count.bit_length() - 1,
                                     capacity.bit_length() - 1)
    assert rows == 1 << bits and bank_count == 1 << bank_bits
    assert capacity == 1 << selector_bits and bank_count * capacity <= rows
    address = list(range(bits))
    output = list(range(bits, bits + word_bits))
    start = bits + word_bits
    banks = [list(range(start + j * word_bits, start + (j + 1) * word_bits))
             for j in range(bank_count)]
    start += bank_count * word_bits
    indicator = list(range(start, start + capacity))
    high = address[bank_bits + selector_bits:]
    selectors = list(range(start + capacity, start + capacity + len(high)))
    width = start + capacity + len(high)
    table = [((5 * row) ^ (row >> 1) ^ 1) & ((1 << word_bits) - 1)
             if row < max(1, rows // 2) else 0 for row in range(rows)]
    matrices = [
        [[(table[(chunk * capacity + h) * bank_count + bank] >> bit) & 1
          for h in range(capacity)]
         for bank in range(bank_count) for bit in range(word_bits)]
        for chunk in range(rows // (bank_count * capacity))]
    selected, ranks, _ = _selected_shear(
        high, indicator, [wire for bank in banks for wire in bank], selectors,
        matrices, native=native)
    low = address[bank_bits:bank_bits + selector_bits]
    reverse = _inverse if native else lambda word: list(reversed(word))
    indicator_word = (_indicator if native else _classical_indicator)(low, indicator)
    loader = selected + indicator_word + reverse(selected) + reverse(indicator_word)
    copy = [('CX', source, target) for source, target in zip(banks[0], output)]
    if native:
        router = _router(address[:bank_bits], banks)
        copy = [('C', copy)]
    else:
        router = _classical_router(address[:bank_bits], banks)
    selected_copy = router + copy + reverse(router)
    query = loader + selected_copy + reverse(loader) + selected_copy
    return width, address, output, table, query, selected, ranks, selectors


class AmortizedDirtyLookupTests(unittest.TestCase):
    def test_controlled_rectangular_shear_rank_normalization_native_phases(self):
        for height, width in ((1, 1), (1, 3), (2, 2), (2, 3), (3, 2)):
            inputs = list(range(1, width + 1))
            outputs = list(range(width + 1, width + height + 1))
            for encoded in range(1 << (height * width)):
                matrix = [[(encoded >> (i * width + j)) & 1 for j in range(width)]
                          for i in range(height)]
                images_of_inputs = {sum(((sum(matrix[i][j] * (y >> j & 1)
                                              for j in range(width))) % 2) << i
                                       for i in range(height))
                                    for y in range(1 << width)}
                rank = len(images_of_inputs).bit_length() - 1
                schedule, claimed_rank = _controlled_shear(matrix, 0, inputs, outputs,
                                                            native=True)
                self.assertEqual(claimed_rank, rank)
                self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 4 if rank else 0)
                self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'),
                                 6 * rank + rank % 2)
                images = []
                for basis in range(1 << (1 + width + height)):
                    change = sum((sum(matrix[i][j] * (basis >> inputs[j] & 1)
                                      for j in range(width)) % 2) << outputs[i]
                                 for i in range(height)) if basis & 1 else 0
                    images.append(basis ^ change)
                expected = np.eye(len(images))[:, np.argsort(images)]
                np.testing.assert_allclose(_word_matrix(1 + width + height, _word(schedule)),
                                           expected, atol=ATOL, rtol=0)
                np.testing.assert_allclose(
                    _word_matrix(1 + width + height, _word(_inverse(schedule))),
                    expected.conj().T, atol=ATOL, rtol=0)

    def test_two_pass_dfs_symbolic_all_input_return_and_completed_lookup(self):
        for bits in range(4):
            address = list(range(bits))
            inputs, outputs = [bits, bits + 1], list(range(bits + 2, bits + 5))
            selectors = list(range(bits + 5, 2 * bits + 5))
            matrices = []
            for chunk in range(1 << bits):
                encoded = ((chunk + 1) * 53) ^ (chunk << 2) ^ (chunk >> 1)
                matrices.append([[(encoded >> (2 * i + j)) & 1 for j in range(2)]
                                 for i in range(3)])
            word, _, _ = _selected_shear(address, inputs, outputs, selectors, matrices)
            initial = [{1 << wire} for wire in range(2 * bits + 5)]
            for nonzero_selectors in (False, True):
                values = [set(value) for value in initial]
                if nonzero_selectors:
                    for wire in selectors:
                        values[wire] = {0}
                actual = _symbolic(word, values)
                self.assertEqual(actual, _expected_shear(values, address, inputs, outputs, matrices))
                self.assertEqual(_symbolic(list(reversed(word)), actual), values)
        for parameters in ((32, 2, 2, 2), (8, 3, 2, 1), (2, 5, 1, 2), (1, 3, 1, 1)):
            width, address, output, table, query, _, _, _ = _amortized_fixture(*parameters)
            initial = [{1 << wire} for wire in range(width)]
            self.assertEqual(_symbolic(query, initial), _expected(initial, address, output, table))
            if address:
                inactive = [set(value) for value in initial]
                inactive[address[-1]] = {0}
                self.assertEqual(_symbolic(query, inactive), inactive)

    def test_completed_nine_wire_queries_literal_phase_and_actual_inverse(self):
        for banks, capacity in ((1, 2), (2, 1)):
            width, address, output, table, query, _, _, selectors = (
                _amortized_fixture(8, 1, banks, capacity, native=True))
            self.assertEqual(width, 9)
            self.assertEqual(len(selectors), 2)
            images = [basis ^ (table[basis & 7] << output[0]) for basis in range(1 << width)]
            expected = np.eye(1 << width)[:, np.argsort(images)]
            np.testing.assert_allclose(_word_matrix(width, _word(query)), expected,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(width, _word(_inverse(query))),
                                       expected.conj().T, atol=ATOL, rtol=0)

    def test_emitted_resources_rank_zero_no_high_bits_and_distinct_t_targets(self):
        for bits in range(6):
            rows = 1 << bits
            for bank_bits in range(bits + 1):
                for low_bits in range(bits - bank_bits + 1):
                    for word_bits in (1, 2):
                        banks, capacity = 1 << bank_bits, 1 << low_bits
                        high_bits, chunks = bits - bank_bits - low_bits, rows // (banks * capacity)
                        width, _, _, _, query, selected, ranks, selectors = (
                            _amortized_fixture(rows, word_bits, banks, capacity, native=True))
                        self.assertEqual(width, bits + word_bits + banks * word_bits
                                         + capacity + high_bits)
                        self.assertEqual(len(selectors), high_bits)
                        internal = 8 * chunks - 16 if high_bits else 0
                        selected_count = (7 * internal + 2 * sum(6 * rank + rank % 2 for rank in ranks)
                                          if high_bits else 0)
                        selected_depth = (4 * internal + 8 * sum(rank > 0 for rank in ranks)
                                          if high_bits else 0)
                        self.assertEqual(sum(len(gates) for kind, gates in selected if kind == 'T'),
                                         selected_count)
                        self.assertEqual(sum(kind == 'T' for kind, _ in selected), selected_depth)
                        count = sum(len(gates) for kind, gates in query if kind == 'T')
                        self.assertEqual(count, 4 * selected_count
                                         + 8 * _router_t_count(capacity)
                                         + 4 * _router_t_count(banks, word_bits))
                        self.assertEqual(sum(kind == 'T' for kind, _ in query),
                                         4 * selected_depth + 32 * low_bits + 16 * bank_bits)
                        self.assertLessEqual(count, 400 * (rows // banks + banks * word_bits))
                        for kind, gates in query:
                            self.assertTrue(all(0 <= wire < width for gate in gates for wire in gate[1:]))
                            if kind == 'T':
                                self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                                targets = [gate[1] for gate in gates]
                                self.assertEqual(len(targets), len(set(targets)))
        zero, rank = _controlled_shear([[0, 0], [0, 0]], 0, [1, 2], [3, 4], native=True)
        self.assertEqual((zero, rank), ([], 0))

    def test_omitting_root_disabled_traversal_leaks_unknown_selector(self):
        # The selected row is zero, but dirty selector one causes the first
        # traversal alone to apply the other row's shear.
        word, _, first_only = _selected_shear([0], [1], [2], [3], [[[1]], [[0]]])
        basis = (1 << 0) | (1 << 1) | (1 << 3)
        self.assertEqual(_basis_action(basis, word), basis)
        self.assertEqual(_basis_action(basis, first_only), basis ^ (1 << 2))
        self.assertEqual(_basis_action(basis, list(reversed(word))), basis)


if __name__ == '__main__':
    unittest.main()
