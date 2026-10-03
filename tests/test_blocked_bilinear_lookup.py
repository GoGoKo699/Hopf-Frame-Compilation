"""Bounded audits of a dirty-controlled bilinear block lookup.

Native small matrices retain literal phases, arbitrary dirty helpers, and
actual inverses. Larger complete queries use exact Boolean-polynomial
gates, with a reduced C3X for the already-audited controlled bilinear
action. Every input wire, including both indicator words and traversal
selectors, remains a symbolic variable. These checks do not emit the
asymptotic chunked indicators or prove the width-sensitive frame bounds.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_amortized_dirty_lookup import _normal_form
    from .test_bilinear_dirty_lookup import _bilinear
    from .test_operator_source_compiler import _word_matrix
    from .test_parallel_dirty_lookup import _classical_indicator, _expected, _multiply
    from .test_t_depth import _ccz_batch, _inverse, _native_schedule, _word
except ImportError:
    from test_amortized_dirty_lookup import _normal_form
    from test_bilinear_dirty_lookup import _bilinear
    from test_operator_source_compiler import _word_matrix
    from test_parallel_dirty_lookup import _classical_indicator, _expected, _multiply
    from test_t_depth import _ccz_batch, _inverse, _native_schedule, _word


ATOL = 8e-11


def _controlled_bilinear(matrix, control, left, right, output, helper, *, native=False):
    wires = [control, *left, *right, output, helper]
    assert len(wires) == len(set(wires))
    before, rank = _normal_form(matrix, right, left)
    if not rank:
        return [], 0
    left_set = set(left)
    before = [('CX', gate[2], gate[1]) if gate[1] in left_set else gate
              for gate in before]
    if not native:
        # Semantic all-input interface; literal phases/helper return are
        # checked independently below using the D,E,D†,E† native word.
        middle = [('C3X', control, left[j], right[j], output) for j in range(rank)]
        return before + middle + list(reversed(before)), rank
    pre = [('C', before + [('H', output)])]
    diagonal = _ccz_batch(helper, list(zip(left[:rank], right[:rank])))
    flip = _native_schedule([('CCX', control, output, helper)])
    return (pre + diagonal + flip + _inverse(diagonal) + _inverse(flip)
            + _inverse(pre)), rank


def _selected_bilinear(address, left, right, output, selectors, helper, blocks,
                       *, native=False):
    """Enabled DFS followed by the actual inverse of root-disabled DFS."""
    wires = address + left + right + output + selectors + [helper]
    assert len(wires) == len(set(wires)) and len(selectors) == len(address)
    assert len(blocks) == 1 << len(address)
    reverse = _inverse if native else lambda word: list(reversed(word))
    if not address:
        word, ranks = [], []
        for target, matrix in zip(output, blocks[0]):
            leaf, rank = _bilinear(matrix, left, right, target, native=native)
            word += leaf
            ranks.append(rank)
        return word, [ranks], word
    leaves, ranks = [], []
    for block in blocks:
        leaf, leaf_ranks = [], []
        for target, matrix in zip(output, block):
            part, rank = _controlled_bilinear(
                matrix, selectors[-1], left, right, target, helper, native=native)
            leaf += part
            leaf_ranks.append(rank)
        leaves.append(leaf)
        ranks.append(leaf_ranks)

    def traversal(root_enabled):
        def descend(level, prefix):
            if level == len(address):
                return leaves[prefix]
            result = []
            for value in (0, 1):
                control, target = address[level], selectors[level]
                negative = [('X', control)] if value == 0 else []
                toggle = ('CX', control, target) if level == 0 else (
                    'CCX', selectors[level - 1], control, target)
                update = negative + [toggle] + list(reversed(negative))
                if level == 0 and not root_enabled:
                    update = []
                if native:
                    update = _native_schedule(update)
                result += (update + descend(level + 1, prefix | (value << level))
                           + reverse(update))
            return result
        return descend(0, 0)

    first, second = traversal(True), traversal(False)
    return first + reverse(second), ranks, first


def _symbolic(word, initial):
    values = [set(value) for value in initial]
    for name, *wires in word:
        assert name in ('X', 'CX', 'CCX', 'C3X')
        product = {0}
        for control in wires[:-1]:
            product = _multiply(product, values[control])
        values[wires[-1]] ^= product
    return values


def _expected_selected(initial, address, left, right, output, blocks):
    result = [set(value) for value in initial]
    for block, matrices in enumerate(blocks):
        select = {0}
        for bit, wire in enumerate(address):
            literal = initial[wire] if (block >> bit) & 1 else initial[wire] ^ {0}
            select = _multiply(select, literal)
        for target, matrix in zip(output, matrices):
            for i, row in enumerate(matrix):
                for j, bit in enumerate(row):
                    if bit:
                        result[target] ^= _multiply(
                            select, _multiply(initial[left[i]], initial[right[j]]))
    return result


def _fixture(a_bits, b_bits, high_bits, word_bits):
    bits = a_bits + b_bits + high_bits
    height, width = 1 << a_bits, 1 << b_bits
    address = list(range(bits))
    output = list(range(bits, bits + word_bits))
    left = list(range(bits + word_bits, bits + word_bits + height))
    right = list(range(left[-1] + 1, left[-1] + 1 + width))
    selectors = list(range(right[-1] + 1, right[-1] + 1 + high_bits))
    helper = right[-1] + 1 + high_bits
    total = helper + 1
    # Unequal binary block matrices, with zero rows in the upper half.
    table = [(((13 * row) ^ (row >> 1) ^ (row >> 2) ^ 5)
              & ((1 << word_bits) - 1)) if row < max(1, 1 << (bits - 1)) else 0
             for row in range(1 << bits)]
    blocks = [
        [[[table[i + height * j + height * width * c] >> bit & 1
           for j in range(width)] for i in range(height)] for bit in range(word_bits)]
        for c in range(1 << high_bits)
    ]
    high = address[a_bits + b_bits:]
    selected, ranks, first = _selected_bilinear(
        high, left, right, output, selectors, helper, blocks)
    ia = _classical_indicator(address[:a_bits], left)
    ib = _classical_indicator(address[a_bits:a_bits + b_bits], right)

    def echo(middle):
        reverse = lambda word: list(reversed(word))
        return (middle + ia + reverse(middle) + ib + middle + reverse(ia)
                + reverse(middle) + reverse(ib))

    return dict(total=total, address=address, output=output, left=left, right=right,
                selectors=selectors, helper=helper, high=high, blocks=blocks,
                table=table, query=echo(selected), selected=selected, ranks=ranks,
                first_only=echo(first), echo=echo)


class BlockedBilinearLookupTests(unittest.TestCase):
    def test_native_controlled_bilinear_literal_phase_dirty_return_and_counts(self):
        matrices = ([[0]], [[1]], [[1, 1], [1, 1]],
                    [[1, 1], [1, 0]], [[1, 1, 0], [1, 0, 1]])
        for matrix in matrices:
            height, width = len(matrix), len(matrix[0])
            control, target = 0, 1
            left = list(range(2, 2 + height))
            right = list(range(2 + height, 2 + height + width))
            helper, total = 2 + height + width, 3 + height + width
            schedule, rank = _controlled_bilinear(
                matrix, control, left, right, target, helper, native=True)
            images = []
            for basis in range(1 << total):
                bit = sum(matrix[i][j] * (basis >> left[i] & 1)
                          * (basis >> right[j] & 1)
                          for i in range(height) for j in range(width)) % 2
                images.append(basis ^ (((basis & 1) * bit) << target))
            expected = np.eye(1 << total)[:, np.argsort(images)]
            np.testing.assert_allclose(_word_matrix(total, _word(schedule)), expected,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(_word_matrix(total, _word(_inverse(schedule))),
                                       expected.conj().T, atol=ATOL, rtol=0)
            self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 16 if rank else 0)
            self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'),
                             12 * rank + 2 * (rank % 2) + 14 if rank else 0)
            for kind, gates in schedule:
                if kind == 'T':
                    targets = [gate[1] for gate in gates]
                    self.assertEqual(len(targets), len(set(targets)))

    def test_native_selected_block_on_every_dirty_selector_and_helper_input(self):
        # Eight-wire fixture: one high address, two input bits on each side,
        # one output, one selector, one returned dirty helper.
        address, left, right, output, selectors, helper = [0], [1, 2], [3, 4], [5], [6], 7
        blocks = [[[[1, 1], [1, 0]]], [[[0, 1], [0, 1]]]]
        schedule, _, _ = _selected_bilinear(
            address, left, right, output, selectors, helper, blocks, native=True)
        images = []
        for basis in range(256):
            matrix = blocks[basis & 1][0]
            bit = sum(matrix[i][j] * (basis >> left[i] & 1) * (basis >> right[j] & 1)
                      for i in range(2) for j in range(2)) % 2
            images.append(basis ^ (bit << output[0]))
        expected = np.eye(256)[:, np.argsort(images)]
        np.testing.assert_allclose(_word_matrix(8, _word(schedule)), expected,
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(_word_matrix(8, _word(_inverse(schedule))),
                                   expected.conj().T, atol=ATOL, rtol=0)

    def test_complete_query_symbolic_all_input_return_and_actual_inverse(self):
        for parameters in ((1, 1, 1, 1), (1, 1, 2, 2), (0, 1, 2, 1),
                           (1, 0, 2, 2), (1, 1, 0, 2)):
            f = _fixture(*parameters)
            initial = [{1 << wire} for wire in range(f['total'])]
            self.assertEqual(_symbolic(f['selected'], initial), _expected_selected(
                initial, f['high'], f['left'], f['right'], f['output'], f['blocks']))
            actual = _symbolic(f['query'], initial)
            self.assertEqual(actual, _expected(initial, f['address'], f['output'], f['table']))
            self.assertEqual(_symbolic(list(reversed(f['query'])), actual), initial)
            # The designated upper-half rows are zero, with all work arbitrary.
            inactive = [set(value) for value in initial]
            inactive[f['address'][-1]] = {0}
            self.assertEqual(_symbolic(f['query'], inactive), inactive)

    def test_two_pass_emitted_native_ledger_and_disjoint_t_targets(self):
        for high_bits in (1, 2, 3):
            for word_bits in (1, 2):
                f = _fixture(1, 1, high_bits, word_bits)
                schedule, ranks, _ = _selected_bilinear(
                    f['high'], f['left'], f['right'], f['output'], f['selectors'],
                    f['helper'], f['blocks'], native=True)
                positive = [rank for row in ranks for rank in row if rank]
                traversal_toffolis = 8 * ((1 << high_bits) - 2)
                self.assertEqual(sum(kind == 'T' for kind, _ in schedule),
                                 32 * len(positive) + 4 * traversal_toffolis)
                self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'),
                                 2 * sum(12 * r + 2 * (r % 2) + 14 for r in positive)
                                 + 7 * traversal_toffolis)
                for kind, gates in schedule:
                    self.assertTrue(all(0 <= wire < f['total']
                                        for gate in gates for wire in gate[1:]))
                    if kind == 'T':
                        targets = [gate[1] for gate in gates]
                        self.assertEqual(len(targets), len(set(targets)))

    def test_missing_helper_inverse_and_incorrect_block_selection_are_detected(self):
        # Omitting the final E† fails to return the arbitrary helper.
        diagonal = _ccz_batch(4, [(2, 3)])
        flip = _native_schedule([('CCX', 0, 1, 4)])
        wrong = diagonal + flip + _inverse(diagonal)
        correct = wrong + _inverse(flip)
        self.assertGreater(np.linalg.norm(_word_matrix(5, _word(wrong))
                                          - _word_matrix(5, _word(correct)), 2), 1.9)
        f = _fixture(1, 1, 2, 2)
        initial = [{1 << wire} for wire in range(f['total'])]
        expected = _expected(initial, f['address'], f['output'], f['table'])
        self.assertNotEqual(_symbolic(f['first_only'], initial), expected)
        wrong_selected, _, _ = _selected_bilinear(
            f['high'], f['left'], f['right'], f['output'], f['selectors'],
            f['helper'], list(reversed(f['blocks'])))
        self.assertNotEqual(_symbolic(f['echo'](wrong_selected), initial), expected)


if __name__ == '__main__':
    unittest.main()
