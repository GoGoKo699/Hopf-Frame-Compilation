"""Literal finite checks of the parallel dirty-indicator construction.

The new indicator uses four-corner Boolean differences, exact Toffolis,
and returned dirty scratch. Polynomial checks keep every input symbolic;
the native phase audit uses at most eight wires. These tests verify fragile
circuit identities and schedules, not an asymptotic theorem or a complete
frame compiler.
"""
from __future__ import annotations

from itertools import zip_longest
import unittest

import numpy as np

try:
    from .test_operator_source_compiler import _word_matrix
    from .test_t_depth import _native_schedule, _inverse, _word
except ImportError:
    from test_operator_source_compiler import _word_matrix
    from test_t_depth import _native_schedule, _inverse, _word


def _reverse(stages):
    """Actual inverse of a chronological X/CX/CCX circuit."""
    return [list(reversed(stage)) for stage in reversed(stages)]


def _parallel(circuits):
    return [sum(row, []) for row in zip_longest(*circuits, fillvalue=[])]


def _fanout(control, targets):
    """Exact dirty-target fanout: inverse tree, root CNOT, forward tree."""
    if not targets:
        return []
    tree, covered = [], 1
    while covered < len(targets):
        count = min(covered, len(targets) - covered)
        tree.append([('CX', targets[j], targets[covered + j])
                     for j in range(count)])
        covered += count
    return _reverse(tree) + [[('CX', control, targets[0])]] + tree


def _masked_outer(edges, scratch):
    """XOR each u*v into its own target using two dirty copies per edge."""
    count = len(edges)
    assert len(scratch) >= 2 * count
    left, right = scratch[:count], scratch[count:2 * count]
    controls = set(q for u, v, _ in edges for q in (u, v))
    targets = [target for _, _, target in edges]
    assert len(set(targets)) == count
    assert not controls.intersection(targets)
    assert len(set(scratch[:2 * count])) == 2 * count
    assert not set(scratch[:2 * count]).intersection(controls | set(targets))
    fan_left = _parallel([
        _fanout(control, [left[j] for j, (u, _, _) in enumerate(edges)
                          if u == control])
        for control in sorted(set(u for u, _, _ in edges))])
    fan_right = _parallel([
        _fanout(control, [right[j] for j, (_, v, _) in enumerate(edges)
                          if v == control])
        for control in sorted(set(v for _, v, _ in edges))])
    batch = [[('CCX', left[j], right[j], target)
              for j, (_, _, target) in enumerate(edges)]]
    return (batch + fan_left + batch + fan_right + batch
            + _reverse(fan_left) + batch + _reverse(fan_right))


def _scratch_size(bits):
    if bits <= 1:
        return 0
    low, high = bits // 2, (bits + 1) // 2
    return (1 << low) + (1 << high) + 2 * (1 << bits)


def _indicator(address, output, scratch):
    """Recursive indicator XOR, with scratch reused across exact subcalls."""
    bits = len(address)
    assert len(output) == 1 << bits
    assert len(scratch) >= _scratch_size(bits)
    assert len(set(address + output + scratch)) == len(address + output + scratch)
    if bits == 0:
        return [[('X', output[0])]]
    if bits == 1:
        return [[('X', output[0])], [('CX', address[0], output[0])],
                [('CX', address[0], output[1])]]
    low = bits // 2
    p, q = 1 << low, 1 << (bits - low)
    left, right = scratch[:p], scratch[p:p + q]
    pool = scratch[p + q:]
    a = _indicator(address[:low], left, pool[:_scratch_size(low)])
    b = _indicator(address[low:], right, pool[:_scratch_size(bits - low)])
    edges = [(left[i], right[j], output[i + p * j])
             for j in range(q) for i in range(p)]
    outer = _masked_outer(edges, pool)
    return outer + a + outer + b + outer + _reverse(a) + outer + _reverse(b)


def _native_parallel(stages):
    """Expand disjoint Toffoli batches into disjoint exact native T layers."""
    result = []
    for stage in stages:
        if stage[0][0] != 'CCX':
            assert all(gate[0] in ('X', 'CX') for gate in stage)
            result.append(('C', stage))
            continue
        assert all(gate[0] == 'CCX' for gate in stage)
        used = [wire for gate in stage for wire in gate[1:]]
        assert len(used) == len(set(used))
        local = [_native_schedule([gate]) for gate in stage]
        for row in zip(*local):
            assert len(set(kind for kind, _ in row)) == 1
            result.append((row[0][0], sum((word for _, word in row), [])))
    return result


# A Boolean polynomial is its set of square-free monomials over GF(2).
# Integer bit masks encode monomials; zero encodes the constant one.
def _multiply(left, right):
    result = set()
    for a in left:
        for b in right:
            term = a | b
            if term in result:
                result.remove(term)
            else:
                result.add(term)
    return result


def _symbolic(stages, values):
    values = [set(value) for value in values]
    for stage in stages:
        for gate in stage:
            if gate[0] == 'X':
                change = {0}
            elif gate[0] == 'CX':
                change = values[gate[1]]
            else:
                assert gate[0] == 'CCX'
                change = _multiply(values[gate[1]], values[gate[2]])
            values[gate[-1]] = values[gate[-1]] ^ change
    return values


def _indicator_query_fixture(rows, word_bits, bank_count):
    """Actual loader/router/output echoes, including all returned scratch."""
    bits, low_bits = rows.bit_length() - 1, bank_count.bit_length() - 1
    assert rows == 1 << bits and bank_count == 1 << low_bits
    address = list(range(bits))
    output = list(range(bits, bits + word_bits))
    first_bank = bits + word_bits
    banks = [list(range(first_bank + j * word_bits,
                        first_bank + (j + 1) * word_bits))
             for j in range(bank_count)]
    high = address[low_bits:]
    high_rows = rows // bank_count
    first_selector = first_bank + bank_count * word_bits
    selector = list(range(first_selector, first_selector + high_rows))
    scratch = list(range(first_selector + high_rows,
                         first_selector + high_rows + _scratch_size(len(high))))
    width = first_selector + high_rows + len(scratch)
    mask = (1 << word_bits) - 1
    # The upper half is an inactive sector with an identically zero word.
    table = [(((3 * row) ^ (row >> 1) ^ 1) & mask) if row < rows // 2 else 0
             for row in range(rows)]
    linear = [[('CX', selector[h], banks[j][bit])]
              for h in range(high_rows) for j in range(bank_count)
              for bit in range(word_bits) if table[h * bank_count + j] >> bit & 1]
    indicator = _indicator(high, selector, scratch)
    loader = linear + indicator + _reverse(linear) + _reverse(indicator)
    router = []
    for level in range(low_bits):
        stride = 1 << level
        for j in range(0, bank_count, 2 * stride):
            for left, right in zip(banks[j], banks[j + stride]):
                router += [[('CX', left, right)],
                           [('CCX', address[level], right, left)],
                           [('CX', left, right)]]
    copy = [[('CX', source, target)] for source, target in zip(banks[0], output)]
    query = (loader + router + copy + _reverse(router) + _reverse(loader)
             + router + copy + _reverse(router))
    return width, address, output, table, query


class ParallelDirtyLookupTests(unittest.TestCase):
    def test_dirty_fanout_has_exact_action_and_logarithmic_layers(self):
        for count in range(1, 18):
            circuit = _fanout(0, list(range(1, count + 1)))
            initial = [{1 << wire} for wire in range(count + 1)]
            expected = [initial[0]] + [value ^ initial[0] for value in initial[1:]]
            self.assertEqual(_symbolic(circuit, initial), expected)
            self.assertEqual(sum(map(len, circuit)), 2 * count - 1)
            self.assertEqual(len(circuit), 2 * (count - 1).bit_length() + 1)
            for stage in circuit:
                wires = [wire for gate in stage for wire in gate[1:]]
                self.assertEqual(len(wires), len(set(wires)))

    def test_masked_outer_is_literal_native_permutation_and_inverse(self):
        # Two targets of the same product exercise both shared-control
        # fanouts and a parallel native Toffoli layer on all 256 inputs.
        width = 8
        edges = [(0, 1, 2), (0, 1, 3)]
        circuit = _masked_outer(edges, list(range(4, 8)))
        schedule = _native_parallel(circuit)
        self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 16)
        self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'), 56)
        actual = _word_matrix(width, _word(schedule))
        images = [basis ^ (12 if basis & 3 == 3 else 0) for basis in range(1 << width)]
        expected = np.eye(1 << width)[:, np.argsort(images)]
        np.testing.assert_allclose(actual, expected, atol=5e-11, rtol=0)
        np.testing.assert_allclose(
            _word_matrix(width, _word(_inverse(schedule))), expected.conj().T,
            atol=5e-11, rtol=0)

    def test_recursive_indicator_symbolically_restores_every_dirty_input(self):
        # No truth-table or dense-matrix growth: every address, output and
        # scratch bit is an independent symbolic input, including reused work.
        for bits in (1, 2, 3, 4):
            rows = 1 << bits
            address = list(range(bits))
            output = list(range(bits, bits + rows))
            scratch = list(range(bits + rows, bits + rows + _scratch_size(bits)))
            width = bits + rows + len(scratch)
            initial = [{1 << wire} for wire in range(width)]
            circuit = _indicator(address, output, scratch)
            expected = [set(value) for value in initial]
            for row, target in enumerate(output):
                indicator = {0}
                for bit, wire in enumerate(address):
                    factor = initial[wire] if row >> bit & 1 else initial[wire] ^ {0}
                    indicator = _multiply(indicator, factor)
                expected[target] ^= indicator
            actual = _symbolic(circuit, initial)
            self.assertEqual(actual, expected)
            self.assertEqual(_symbolic(_reverse(circuit), actual), initial)
        base = _indicator([0], [1, 2], [])
        expected = np.eye(8)[:, np.argsort([basis ^ (1 << (1 + (basis & 1)))
                                         for basis in range(8)])]
        np.testing.assert_allclose(_word_matrix(3, sum(base, [])), expected,
                                   atol=5e-11, rtol=0)

    def test_indicator_native_layers_and_reused_workspace_ledger(self):
        counts, depths = {1: 0}, {1: 0}
        for bits in range(1, 6):
            rows = 1 << bits
            size = _scratch_size(bits)
            self.assertLessEqual(size, 3 * rows)
            address = list(range(bits))
            output = list(range(bits, bits + rows))
            scratch = list(range(bits + rows, bits + rows + size))
            circuit = _indicator(address, output, scratch)
            schedule = _native_parallel(circuit)
            if bits > 1:
                low, high = bits // 2, (bits + 1) // 2
                self.assertLessEqual(max(_scratch_size(low), _scratch_size(high)), 2 * rows)
                counts[bits] = 2 * counts[low] + 2 * counts[high] + 16 * rows
                depths[bits] = 2 * depths[low] + 2 * depths[high] + 64
            self.assertEqual(sum(gate[0] == 'CCX' for stage in circuit for gate in stage),
                             counts[bits])
            self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'),
                             7 * counts[bits])
            self.assertEqual(sum(kind == 'T' for kind, _ in schedule), depths[bits])
            for kind, gates in schedule:
                self.assertTrue(all(0 <= wire < bits + rows + size
                                    for gate in gates for wire in gate[1:]))
                if kind == 'T':
                    targets = [gate[1] for gate in gates]
                    self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                    self.assertEqual(len(targets), len(set(targets)))

    def test_completed_query_returns_all_work_and_preserves_inactive_rows(self):
        # The first case composes a recursive decoder, multi-bit banks,
        # routing and both dirty echoes. The second has no high address.
        for rows, word_bits, bank_count in ((8, 2, 2), (2, 2, 2)):
            width, address, output, table, query = _indicator_query_fixture(
                rows, word_bits, bank_count)
            initial = [{1 << wire} for wire in range(width)]
            expected = [set(value) for value in initial]
            for row, value in enumerate(table):
                indicator = {0}
                for bit, wire in enumerate(address):
                    factor = initial[wire] if row >> bit & 1 else initial[wire] ^ {0}
                    indicator = _multiply(indicator, factor)
                for bit, target in enumerate(output):
                    if value >> bit & 1:
                        expected[target] ^= indicator
            actual = _symbolic(query, initial)
            self.assertEqual(actual, expected)
            self.assertEqual(_symbolic(_reverse(query), actual), initial)
            inactive = [set(value) for value in initial]
            inactive[address[-1]] = {0}  # Highest address bit fixed to one.
            self.assertEqual(_symbolic(query, inactive), inactive)


if __name__ == '__main__':
    unittest.main()
