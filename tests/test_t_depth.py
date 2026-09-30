"""Finite native-circuit audits of shared-control T-depth schedules.

T-depth here permits arbitrary Clifford circuits between T layers; it is
not total circuit depth. These full-matrix fixtures verify literal phases
and arbitrary dirty-input return, but do not prove asymptotic lookup or
whole-frame bounds.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_operator_source_compiler import (
        _adjoint, _dirty_mask_lookup, _word_matrix,
    )
except ImportError:
    from test_operator_source_compiler import (
        _adjoint, _dirty_mask_lookup, _word_matrix,
    )


ATOL = 5e-11


def _ccz_batch(control, pairs):
    """Four T layers for prod_i CCZ(control,a_i,b_i), without helpers."""
    wires = [control] + [q for pair in pairs for q in pair]
    assert pairs and len(wires) == len(set(wires))
    count = len(pairs)
    first = [('T', q) for pair in pairs for q in pair]
    if count % 2:
        first.append(('T', control))
    common = [('S', control)] * (count // 2)
    shared_parities = [('CX', control, q) for pair in pairs for q in pair]
    pair_parities = [('CX', a, b) for a, b in pairs]
    triple_parities = [('CX', control, b) for _, b in pairs] + pair_parities
    return [
        ('C', common),
        ('T', first),
        ('C', shared_parities),
        ('T', [('TDG', q) for pair in pairs for q in pair]),
        ('C', _adjoint(shared_parities) + pair_parities),
        ('T', [('TDG', b) for _, b in pairs]),
        ('C', _adjoint(pair_parities) + triple_parities),
        ('T', [('T', b) for _, b in pairs]),
        ('C', _adjoint(triple_parities)),
    ]


def _fredkin_batch(control, pairs, control_value=1):
    outer = [('CX', a, b) for a, b in pairs]
    hadamards = [('H', a) for a, _ in pairs]
    negative = [('X', control)] if control_value == 0 else []
    return ([('C', negative + outer + hadamards)]
            + _ccz_batch(control, pairs)
            + [('C', hadamards + outer + negative)])


def _native_schedule(word):
    """Expand any miniature loader Toffolis into exact four-layer words."""
    schedule = []
    for gate in word:
        if gate[0] == 'CCX':
            _, left, right, target = gate
            schedule += ([('C', [('H', target)])]
                         + _ccz_batch(left, [(target, right)])
                         + [('C', [('H', target)])])
        else:
            assert gate[0] not in ('T', 'TDG')
            schedule.append(('C', [gate]))
    return schedule


def _inverse(schedule):
    return [(kind, _adjoint(word)) for kind, word in reversed(schedule)]


def _word(schedule):
    return [gate for _, stage in schedule for gate in stage]


def _router(address, banks):
    schedule = []
    for level in range(len(address)):
        stride = 1 << level
        pairs = [pair
                 for bank in range(0, len(banks), 2 * stride)
                 for pair in zip(banks[bank], banks[bank + stride])]
        schedule += _fredkin_batch(address[level], pairs)
    return schedule


def _dirty_query(rows, word_bits, bank_count):
    address_bits = rows.bit_length() - 1
    low_bits = bank_count.bit_length() - 1
    address = list(range(address_bits))
    output = list(range(address_bits, address_bits + word_bits))
    bank_start = address_bits + word_bits
    banks = [list(range(bank_start + j * word_bits,
                        bank_start + (j + 1) * word_bits))
             for j in range(bank_count)]
    bank_targets = [q for bank in banks for q in bank]
    high = address[low_bits:]
    selector = bank_start + bank_count * word_bits
    width = selector + (len(high) == 2)
    mask = (1 << word_bits) - 1
    table = [((3 * row) ^ (row >> 1) ^ 1) & mask for row in range(rows)]
    chunks = [sum(table[chunk * bank_count + j] << (j * word_bits)
                  for j in range(bank_count))
              for chunk in range(rows // bank_count)]
    if len(high) == 2:
        loader_word = _dirty_mask_lookup(high, bank_targets, selector, chunks)
    elif len(high) == 1:
        loader_word = []
        for chunk, data in enumerate(chunks):
            if chunk == 0:
                loader_word.append(('X', high[0]))
            loader_word += [('CX', high[0], q)
                            for j, q in enumerate(bank_targets) if (data >> j) & 1]
            if chunk == 0:
                loader_word.append(('X', high[0]))
    else:
        assert not high
        loader_word = [('X', q) for j, q in enumerate(bank_targets)
                       if (chunks[0] >> j) & 1]
    loader = _native_schedule(loader_word)
    router = _router(address[:low_bits], banks)
    copy = [('C', [('CX', source, target) for source, target in zip(banks[0], output)])]
    # The second copy cancels the initially unknown selected bank word.
    query = (loader + router + copy + _inverse(router) + _inverse(loader)
             + router + copy + _inverse(router))
    return width, address_bits, table, loader, router, query


class TDepthTests(unittest.TestCase):
    def assert_schedule(self, schedule, width):
        for kind, gates in schedule:
            if kind == 'T':
                self.assertTrue(gates)
                self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                targets = [gate[1] for gate in gates]
                self.assertEqual(len(targets), len(set(targets)))
            else:
                self.assertTrue(all(gate[0] in ('X', 'CX', 'H', 'S', 'SDG', 'Z')
                                    for gate in gates))
            self.assertTrue(all(0 <= q < width for gate in gates for q in gate[1:]))

    def test_shared_control_ccz_phase_polynomial_and_four_layers(self):
        for count in (1, 2, 3, 4):
            width = 1 + 2 * count
            pairs = [(1 + 2 * j, 2 + 2 * j) for j in range(count)]
            schedule = _ccz_batch(0, pairs)
            self.assert_schedule(schedule, width)
            self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 4)
            self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'),
                             6 * count + count % 2)
            expected = []
            for basis in range(1 << width):
                control = basis & 1
                exponent = count * control
                desired = 0
                for a, b in pairs:
                    left, right = (basis >> a) & 1, (basis >> b) & 1
                    exponent += (left + right - (control ^ left) - (control ^ right)
                                 - (left ^ right) + (control ^ left ^ right))
                    desired += 4 * control * left * right
                self.assertEqual((exponent - desired) % 8, 0)
                expected.append((-1) ** (desired // 4))
            actual = _word_matrix(width, _word(schedule))
            np.testing.assert_allclose(actual, np.diag(expected), atol=ATOL, rtol=0)

    def test_fredkin_batches_are_literal_all_input_permutations(self):
        for count in (1, 2, 3, 4):
            width = 1 + 2 * count
            pairs = [(1 + 2 * j, 2 + 2 * j) for j in range(count)]
            for control_value in (0, 1):
                schedule = _fredkin_batch(0, pairs, control_value)
                self.assert_schedule(schedule, width)
                actual = _word_matrix(width, _word(schedule))
                images = []
                for basis in range(1 << width):
                    output = basis
                    if (basis & 1) == control_value:
                        for left, right in pairs:
                            if ((basis >> left) ^ (basis >> right)) & 1:
                                output ^= (1 << left) | (1 << right)
                    images.append(output)
                expected = np.eye(1 << width)[:, np.argsort(images)]
                np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
                np.testing.assert_allclose(
                    _word_matrix(width, _word(_inverse(schedule))), expected.conj().T,
                    atol=ATOL, rtol=0)

    def test_dirty_selectswap_queries_keep_router_depth_logarithmic(self):
        # Cover two router levels, multi-bit words, and an arbitrary dirty
        # selector for the high-address loader, without a large simulation.
        for rows, word_bits, bank_count in ((4, 1, 4), (4, 2, 2), (8, 1, 2)):
            width, address_bits, table, loader, router, query = _dirty_query(
                rows, word_bits, bank_count)
            self.assertLessEqual(width, 8)
            self.assert_schedule(query, width)
            router_depth = sum(kind == 'T' for kind, _ in router)
            loader_depth = sum(kind == 'T' for kind, _ in loader)
            self.assertEqual(router_depth, 4 * (bank_count.bit_length() - 1))
            self.assertEqual(sum(kind == 'T' for kind, _ in query),
                             2 * loader_depth + 4 * router_depth)
            actual = _word_matrix(width, _word(query))
            images = [basis ^ (table[basis & (rows - 1)] << address_bits)
                      for basis in range(1 << width)]
            expected = np.eye(1 << width)[:, np.argsort(images)]
            np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
