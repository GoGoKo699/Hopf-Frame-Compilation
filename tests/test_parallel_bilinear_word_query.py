"""Bounded exact audits of the R91 parallel bilinear word query.

Phase/permutation words are evaluated on every basis input using integer
phases modulo eight, without dense large matrices. A separate five-wire
native H/CNOT/T check validates the shared-control Toffoli substitution
used by the larger controlled fixture. Hadamard conjugation of a complete
output-diagonal identity then gives the full XOR operator, including all
output columns, arbitrary dirty inputs, and literal phase.

These tests do not supply an asymptotic emitter or a query lower bound.
"""
from __future__ import annotations

from math import isqrt
import unittest

import numpy as np


def _inverse(word):
    return [(gate[0], gate[1], -gate[2]) if gate[0] == 'T' else gate
            for gate in reversed(word)]


def _parity_phase(target, controls, power):
    compute = [('CX', control, target) for control in controls]
    return compute + [('T', target, power)] + list(reversed(compute))


def _relative_grid(y, z, dirty):
    """Four physical T layers, with each dirty[i][j] a distinct target."""
    layers = [[], [], [], []]
    for i, yi in enumerate(y):
        for j, zj in enumerate(z):
            d = dirty[i][j]
            for layer, controls, power in ((0, (), 1), (1, (yi,), -1),
                                           (2, (zj,), -1),
                                           (3, (yi, zj), 1)):
                layers[layer].extend(_parity_phase(d, controls, power))
    return [gate for layer in layers for gate in layer], layers


def _linear_rows(x, dirty, matrices):
    return [('CX', x[k], dirty[i][j])
            for j, matrix in enumerate(matrices)
            for i, row in enumerate(matrix)
            for k, entry in enumerate(row) if entry]


def _bilinear_phase(y, x, z, dirty, matrices):
    f, _ = _relative_grid(y, z, dirty)
    linear = _linear_rows(x, dirty, matrices)
    return _inverse(f) + linear + f + _inverse(linear)


def _monomial(word, width):
    images = np.arange(1 << width, dtype=np.int64)
    phases = np.zeros(1 << width, dtype=np.int64)
    for gate in word:
        kind = gate[0]
        if kind == 'T':
            phases += gate[2] * ((images >> gate[1]) & 1)
        elif kind == 'X':
            images ^= 1 << gate[1]
        elif kind == 'CX':
            images ^= ((images >> gate[1]) & 1) << gate[2]
        elif kind == 'CCX':
            images ^= (((images >> gate[1]) & 1)
                       & ((images >> gate[2]) & 1)) << gate[3]
        else:
            raise ValueError(kind)
    return images, phases % 8


def _native_matrix(word, width):
    """Only used for a five-wire native substitution, never a large word."""
    matrix = np.eye(1 << width, dtype=complex)
    indices = np.arange(1 << width)
    for gate in word:
        if gate[0] == 'H':
            low = indices[(indices >> gate[1]) & 1 == 0]
            high = low | (1 << gate[1])
            lo, hi = matrix[low].copy(), matrix[high].copy()
            matrix[low] = (lo + hi) / np.sqrt(2)
            matrix[high] = (lo - hi) / np.sqrt(2)
        elif gate[0] == 'T':
            phase = np.exp(1j * np.pi / 4 * gate[2]
                           * ((indices >> gate[1]) & 1))
            matrix *= phase[:, None]
        elif gate[0] == 'CX':
            permutation = indices ^ (((indices >> gate[1]) & 1) << gate[2])
            matrix = matrix[permutation]
        else:
            raise ValueError(gate[0])
    return matrix


def _common_control_toffolis(control, pairs):
    layers = [[('T', control, len(pairs))], [], [], []]
    for y, e in pairs:
        layers[0].extend([('T', y, 1), ('T', e, 1)])
        layers[1].extend(_parity_phase(y, (control,), -1))
        layers[1].extend(_parity_phase(e, (control,), -1))
        layers[2].extend(_parity_phase(e, (y,), -1))
        layers[3].extend(_parity_phase(e, (control, y), 1))
    hadamards = [('H', e) for _, e in pairs]
    return hadamards + [g for layer in layers for g in layer] + hadamards


def _bit(values, wire):
    return (values >> wire) & 1


def _bilinear_exponent(values, y, x, z, matrices, control=None):
    """Independent target phase from the original table matrices."""
    parity = np.zeros_like(values)
    for j, matrix in enumerate(matrices):
        contribution = np.zeros_like(values)
        for i, row in enumerate(matrix):
            for k, entry in enumerate(row):
                if entry:
                    contribution ^= _bit(values, y[i]) & _bit(values, x[k])
        parity ^= contribution & _bit(values, z[j])
    if control is not None:
        parity &= _bit(values, control)
    return 4 * parity


def _pow2floor(value):
    return 1 << (value.bit_length() - 1)


def _allocation(rows, outputs, capacity):
    """Integer-only finite probe of the proof's normalized allocation."""
    if outputs >= rows and rows <= capacity:
        return None
    j = _pow2floor(min(rows, isqrt(rows * outputs), capacity))
    group = min(outputs, capacity)
    h = _pow2floor(min(rows // j, capacity // group))
    return h, j, rows // (h * j), group


class ParallelBilinearWordQueryTests(unittest.TestCase):
    MATRICES = (((1, 1), (0, 1)), ((1, 0), (1, 1)))

    def test_relative_phase_word_on_every_dirty_basis_input(self):
        word, layers = _relative_grid([0], [1], [[2]])
        values = np.arange(8)
        images, phases = _monomial(word, 3)
        np.testing.assert_array_equal(images, values)
        expected = (4 * _bit(values, 2) - 2) * _bit(values, 0) * _bit(values, 1)
        np.testing.assert_array_equal(phases, expected % 8)
        self.assertEqual(len(layers), 4)
        for layer in layers:
            self.assertEqual([g[1] for g in layer if g[0] == 'T'], [2])
        # F is not the CCZ word: its relative CS-dagger phase is literal.
        self.assertTrue(np.any(phases != (4 * _bit(values, 0)
                                         * _bit(values, 1) * _bit(values, 2)) % 8))
        back, inverse_phase = _monomial(word + _inverse(word), 3)
        np.testing.assert_array_equal(back, values)
        np.testing.assert_array_equal(inverse_phase, np.zeros(8, dtype=int))

    def test_shared_grid_uses_distinct_private_phase_targets(self):
        y, x, z, dirty = [0, 1], [2, 3], [4, 5], [[6, 7], [8, 9]]
        f, layers = _relative_grid(y, z, dirty)
        for layer in layers:
            targets = [g[1] for g in layer if g[0] == 'T']
            self.assertEqual(set(targets), {6, 7, 8, 9})
            self.assertEqual(len(targets), len(set(targets)))
            self.assertTrue(all(g[-1] not in y + x + z
                                for g in layer if g[0] == 'CX'))
        word = _bilinear_phase(y, x, z, dirty, self.MATRICES)
        values = np.arange(1 << 10)
        images, phases = _monomial(word, 10)
        np.testing.assert_array_equal(images, values)
        np.testing.assert_array_equal(phases, _bilinear_exponent(
            values, y, x, z, self.MATRICES))
        # Using F instead of its actual inverse leaves a detectable phase.
        linear = _linear_rows(x, dirty, self.MATRICES)
        _, wrong = _monomial(f + linear + f + _inverse(linear), 10)
        self.assertTrue(np.any(wrong != phases))

    def test_native_common_control_substitution_has_literal_phase(self):
        word = _common_control_toffolis(0, [(1, 3), (2, 4)])
        actual = _native_matrix(word, 5)
        images, phases = _monomial([('CCX', 0, 1, 3), ('CCX', 0, 2, 4)], 5)
        expected = np.zeros_like(actual)
        expected[images, np.arange(32)] = 1
        np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=0)
        np.testing.assert_array_equal(phases, np.zeros(32, dtype=int))
        np.testing.assert_allclose(_native_matrix(_inverse(word), 5),
                                   expected.conj().T, atol=2e-14, rtol=0)

    def test_controlled_grid_returns_all_private_dirty_fields(self):
        h, y, x, z = 0, [1, 2], [3, 4], [5, 6]
        e, dirty = [7, 8], [[9, 10], [11, 12]]
        phase = _bilinear_phase(e, x, z, dirty, self.MATRICES)
        toggle = [('CCX', h, yi, ei) for yi, ei in zip(y, e)]
        word = phase + toggle + _inverse(phase) + _inverse(toggle)
        values = np.arange(1 << 13)
        images, phases = _monomial(word, 13)
        np.testing.assert_array_equal(images, values)
        np.testing.assert_array_equal(phases, _bilinear_exponent(
            values, y, x, z, self.MATRICES, h))
        inactive = _bit(values, h) == 0
        np.testing.assert_array_equal(phases[inactive], np.zeros(sum(inactive), dtype=int))
        # The e return is an operator requirement, not a marginal promise.
        broken_images, _ = _monomial(word[:-len(toggle)], 13)
        self.assertTrue(np.any(broken_images != values))

    def test_four_corner_query_on_every_output_and_dirty_input(self):
        # Two address bits, two indicators, two outputs, four private grid
        # wires: 4096 basis inputs, evaluated exactly without dense matrices.
        address, y, x, z = [0, 1], [2, 3], [4, 5], [6, 7]
        dirty = [[8, 9], [10, 11]]
        f = _bilinear_phase(y, x, z, dirty, self.MATRICES)
        ia = [('X', y[0]), ('CX', address[0], y[0]), ('CX', address[0], y[1])]
        ib = [('X', x[0]), ('CX', address[1], x[0]), ('CX', address[1], x[1])]
        word = f + ia + _inverse(f) + ib + f + _inverse(ia) + _inverse(f) + _inverse(ib)
        values = np.arange(1 << 12)
        images, phases = _monomial(word, 12)
        expected = np.zeros_like(values)
        for j, matrix in enumerate(self.MATRICES):
            table_bit = np.array([matrix[a][b] for a, b in zip(
                _bit(values, address[0]), _bit(values, address[1]))])
            expected ^= table_bit & _bit(values, z[j])
        np.testing.assert_array_equal(images, values)
        np.testing.assert_array_equal(phases, 4 * expected)
        # The proved phase is diagonal on every output basis value; H on
        # both outputs conjugates it to exactly the requested XOR table.
        returned, total_phase = _monomial(word + _inverse(word), 12)
        np.testing.assert_array_equal(returned, values)
        np.testing.assert_array_equal(total_phase, np.zeros_like(values))

    def test_paired_router_returns_auxiliary_only_in_complete_selected_swap(self):
        # First certify the literal native shear core on all 2048 inputs.
        # A remains quantum data: only its dirty Y indicator is echoed.
        r, y, a, e = 0, [1, 2], [3, 4], [5, 6]
        dirty = [[7, 8], [9, 10]]
        matrices = (((1, 0), (0, 1)), ((0, 1), (1, 0)))
        bilinear = _bilinear_phase(y, a, e, dirty, matrices)
        indicator = [('X', y[0]), ('CX', r, y[0]), ('CX', r, y[1])]
        shear_core = bilinear + indicator + _inverse(bilinear) + _inverse(indicator)
        values = np.arange(1 << 11)
        images, phases = _monomial(shear_core, 11)
        selected = np.where(_bit(values, r), _bit(values, a[1]), _bit(values, a[0]))
        other = np.where(_bit(values, r), _bit(values, a[0]), _bit(values, a[1]))
        expected = 4 * ((selected & _bit(values, e[0]))
                        ^ (other & _bit(values, e[1])))
        np.testing.assert_array_equal(images, values)
        np.testing.assert_array_equal(phases, expected)

        # H on E conjugates this exact diagonal core to its XOR shear.
        # Compose those independently verified full-space shears below.
        # All 8192 inputs include P, the unknown initial program R, target,
        # both arbitrary arrays, both indicator bits and all four grid bits.
        p, r, target, a, e = 0, 1, 2, [3, 4], [5, 6]
        values = np.arange(1 << 13)

        def run(word, inputs):
            images = inputs.copy()
            for gate in word:
                if gate[0] == 'SHEAR':
                    source, destination = gate[1:]
                    old = images.copy()
                    for slot, destination_bit in enumerate(destination):
                        value = np.where(_bit(old, r),
                                         _bit(old, source[slot ^ 1]),
                                         _bit(old, source[slot]))
                        images ^= value << destination_bit
                elif gate[0] == 'SWAP':
                    difference = _bit(images, gate[1]) ^ _bit(images, gate[2])
                    images ^= (difference << gate[1]) | (difference << gate[2])
                elif gate[0] == 'QUERY':
                    images ^= _bit(images, p) << r
                else:
                    raise ValueError(gate[0])
            return images

        paired = [('SHEAR', a, e), ('SHEAR', e, a), ('SHEAR', a, e)]
        paired += [('SWAP', a[0], e[0]), ('SWAP', a[1], e[1])]
        # Reverse chronological words, not a claimed shortcut inverse.
        translated = [('QUERY',)] + paired + [('QUERY',)] + list(reversed(paired))
        selected_swap = translated + [('SWAP', target, a[0])] + list(reversed(translated))
        actual = run(selected_swap, values)
        chosen_a = np.where(_bit(values, p), _bit(values, a[1]), _bit(values, a[0]))
        difference = chosen_a ^ _bit(values, target)
        expected = values ^ (difference << target)
        expected ^= np.where(_bit(values, p), difference << a[1], difference << a[0])
        np.testing.assert_array_equal(actual, expected)
        np.testing.assert_array_equal(run(list(reversed(selected_swap)), actual), values)
        # The inner paired translation genuinely moves E; the outer word
        # returns E exactly. Omitting this distinction would be incorrect.
        inner = run(paired, values)
        self.assertTrue(np.any(_bit(inner, e[0]) != _bit(values, e[0])))
        for wire in e + [p, r, 7, 8, 9, 10, 11, 12]:
            np.testing.assert_array_equal(_bit(actual, wire), _bit(values, wire))

    def test_integer_allocations_cover_balanced_blocked_and_batched_regimes(self):
        seen = set()
        for address_bits in list(range(17)) + [31, 63, 127]:
            rows = 1 << address_bits
            outputs_to_check = set(range(1, 65))
            outputs_to_check.update(max(1, rows + delta) for delta in (-1, 0, 1))
            outputs_to_check.update([1 << 20, (1 << 31) + 1])
            capacities = {1 << k for k in range(min(address_bits + 4, 22))}
            capacities.update([rows, rows * 2, max(1, rows // 2)])
            for outputs in outputs_to_check:
                root = isqrt(rows * outputs)
                capacities.update([_pow2floor(root), _pow2floor(root) * 2])
                for capacity in capacities:
                    allocation = _allocation(rows, outputs, capacity)
                    if allocation is None:
                        seen.add('direct')
                        self.assertLessEqual(rows * rows, rows * outputs)
                        self.assertLessEqual(rows, capacity)
                        continue
                    h, j, blocks, group = allocation
                    self.assertEqual(h * j * blocks, rows)
                    for size in (h, j, blocks):
                        self.assertEqual(size & (size - 1), 0)
                    self.assertLessEqual(group * h, capacity)
                    self.assertLessEqual(h + j + group * h + h, 4 * capacity)
                    batches = (outputs + group - 1) // group
                    # Independent same-circuit ledgers: roundings may lose
                    # constants, but never an output-length depth term.
                    self.assertLessEqual(blocks * batches * capacity**2,
                                         8 * (capacity**2 + rows * outputs))
                    count = blocks * outputs * h + h + j
                    self.assertLessEqual(count * capacity,
                                         8 * ((root + 1) * capacity + rows * outputs))
                    self.assertLessEqual(h + j, 2 * rows)
                    if capacity**2 >= rows * outputs:
                        seen.add('balanced')
                    elif outputs <= capacity:
                        seen.add('blocked')
                    else:
                        seen.add('batched')
        self.assertEqual(seen, {'direct', 'balanced', 'blocked', 'batched'})


if __name__ == '__main__':
    unittest.main()
