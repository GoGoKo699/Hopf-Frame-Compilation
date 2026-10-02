"""Bounded native and classical audits of masked dirty-sum compression.

Local compressors use two parity CNOTs and two controlled increments on a
fresh arbitrary upper word. Classical words test all local inputs and full
small sum echoes. One native matrix of order 64 audits literal phase and
returned borrowed work. Its slow MCX expansion does not implement the
imported logarithmic-depth incrementer; fast-depth bounds remain analytic.
"""
from __future__ import annotations

import random
import unittest

import numpy as np

try:
    from .test_counter_dirty_indicator import _add, _increment, _primitive, _reverse, _serial
    from .test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from .test_t_depth import _inverse, _native_schedule, _word
except ImportError:
    from test_counter_dirty_indicator import _add, _increment, _primitive, _reverse, _serial
    from test_operator_source_compiler import _adjoint, _basis_action, _word_matrix
    from test_t_depth import _inverse, _native_schedule, _word


ATOL = 5e-11


def _read(basis, wires):
    return sum((basis >> wire & 1) << bit for bit, wire in enumerate(wires))


def _replace(basis, wires, value):
    for bit, wire in enumerate(wires):
        basis = (basis & ~(1 << wire)) | ((value >> bit & 1) << wire)
    return basis


def _controlled_increment(conditions, counter):
    """Literal high-to-low MCX word; counter and conditions are disjoint."""
    controls = [wire for wire, _ in conditions]
    assert len(set(controls + counter)) == len(controls) + len(counter)
    negative = [('X', wire) for wire, value in conditions if not value]
    increment = [('MCX', *controls, *counter[:bit], counter[bit])
                 for bit in reversed(range(len(counter)))]
    return negative + increment + _adjoint(negative) if counter else []


def _compressor(triple, upper):
    """Retain a,p=a xor b and emit s=a xor b xor c plus D+majority."""
    a, b, c = triple
    assert len(set(triple + upper)) == len(triple) + len(upper)
    return ([('CX', a, b), ('CX', b, c)]
            + _controlled_increment([(a, 1), (b, 0)], upper)
            + _controlled_increment([(b, 1), (c, 0)], upper))


def _small_native(word, helper):
    """Expand at most three controls, borrowing one returned dirty wire."""
    expanded = []
    for gate in word:
        if gate[0] != 'MCX':
            expanded.append(gate)
            continue
        controls, target = list(gate[1:-1]), gate[-1]
        assert helper not in gate[1:]
        if len(controls) == 2:
            expanded.append(('CCX', *controls, target))
        else:
            assert len(controls) == 3
            left = ('CCX', controls[0], controls[1], helper)
            right = ('CCX', helper, controls[2], target)
            expanded += [left, right, left, right]
    return _native_schedule(expanded)


def _compress_columns(columns, first_free):
    """One fresh upper word per triple; all triples see the old round."""
    bits = len(columns)
    columns = [list(column) for column in columns]
    word, fresh, history, compressions = [], [], [tuple(map(len, columns))], [0] * bits
    while any(len(column) > 2 for column in columns):
        following = [[] for _ in range(bits)]
        for column, wires in enumerate(columns):
            triples = len(wires) // 3
            compressions[column] += triples
            for offset in range(triples):
                triple = wires[3 * offset:3 * offset + 3]
                upper = list(range(first_free, first_free + bits - column - 1))
                first_free += len(upper)
                fresh += [(wire, higher) for higher, wire in enumerate(upper, column + 1)]
                word += _compressor(triple, upper)
                following[column].append(triple[2])
                for higher, wire in enumerate(upper, column + 1):
                    following[higher].append(wire)
            following[column] += wires[3 * triples:]
        columns = following
        history.append(tuple(map(len, columns)))
    for column, wires in enumerate(columns):
        while len(wires) < 2:
            wires.append(first_free)
            fresh.append((first_free, column))
            first_free += 1
    outputs = [[column[row] for column in columns] for row in range(2)]
    return word, outputs, fresh, history, compressions, first_free


def _sum_fixture(count):
    bits = count.bit_length()
    address = list(range(count))
    counters = [list(range(count + row * bits, count + (row + 1) * bits))
                for row in range(count)]
    accumulator = list(range(count + count * bits, count + (count + 1) * bits))
    columns = [[counter[column] for counter in counters] for column in range(bits)]
    compression, outputs, fresh, history, counts, first_free = _compress_columns(
        columns, accumulator[-1] + 1)
    loader = (compression + _add(outputs[0], accumulator)[0]
              + _add(outputs[1], accumulator)[0] + _adjoint(compression))
    pattern = [index % 2 for index in range(count)]
    increment = []
    for control, counter, value in zip(address, counters, pattern):
        helpers = list(range(first_free, first_free + bits))
        first_free += bits
        increment += _increment(control, counter, helpers, negative=not value)[0]
    echo = _adjoint(loader) + increment + loader + _adjoint(increment)
    return dict(bits=bits, address=address, counters=counters, accumulator=accumulator,
                compression=compression, outputs=outputs, fresh=fresh,
                history=history, counts=counts, loader=loader, echo=echo,
                pattern=pattern, width=first_free)


class DirtySumInterfaceTests(unittest.TestCase):
    def test_local_compressors_preserve_weight_on_every_input_and_actual_inverse(self):
        for bits in range(2, 7):
            for column in range(bits):
                upper = list(range(3, 3 + bits - column - 1))
                word = _compressor([0, 1, 2], upper)
                modulus = 1 << (bits - column)
                for basis in range(1 << (3 + len(upper))):
                    a, b, c = (basis >> wire & 1 for wire in range(3))
                    majority = int(a + b + c >= 2)
                    expected = _replace(basis, [1, 2], (a ^ b) | ((a ^ b ^ c) << 1))
                    expected = _replace(expected, upper,
                                        (_read(basis, upper) + majority) % (1 << len(upper)))
                    actual = _basis_action(basis, word)
                    self.assertEqual(actual, expected)
                    self.assertEqual(_basis_action(actual, _adjoint(word)), basis)
                    before = a + b + c + 2 * _read(basis, upper)
                    after = (actual >> 2 & 1) + 2 * _read(actual, upper)
                    self.assertEqual(before % modulus, after % modulus)
        # A canonical parity/carry extraction cannot replace the dirty upper
        # word: three distinct weight-one inputs all demand the pair (1,0).
        canonical_pairs = [(value.bit_count() % 2, int(value.bit_count() >= 2))
                           for value in range(8)]
        self.assertEqual(canonical_pairs.count((1, 0)), 3)
        self.assertGreater(canonical_pairs.count((1, 0)), 2)  # one garbage bit

    def test_native_compressor_has_literal_phase_and_returns_borrowed_helper(self):
        upper = [3, 4]
        word = _compressor([0, 1, 2], upper)
        schedule = _small_native(word, 5)
        images = [_basis_action(basis, word) for basis in range(64)]
        expected = np.eye(64)[:, images]
        actual = _word_matrix(6, _word(schedule))
        inverse = _word_matrix(6, _word(_inverse(schedule)))
        self.assertEqual(_word(_inverse(schedule)), _adjoint(_word(schedule)))
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(inverse, expected.conj().T, atol=ATOL, rtol=0)
        np.testing.assert_allclose(inverse @ actual, np.eye(64), atol=ATOL, rtol=0)
        self.assertTrue(all((image >> 5 & 1) == (basis >> 5 & 1)
                            for basis, image in enumerate(images)))
        # Deleting a native phase must fail even though the classical word is
        # unchanged. This checks literal phase, not merely a basis permutation.
        missing = list(_word(schedule))
        missing.pop(next(index for index, gate in enumerate(missing)
                         if gate[0] in ('T', 'TDG')))
        self.assertGreater(np.linalg.norm(_word_matrix(6, missing) - expected), 0.5)

    def test_full_masked_sum_echo_cancels_offsets_and_restores_all_work(self):
        generator = random.Random(0xD175)
        for count in (3, 4, 5):
            fixture = _sum_fixture(count)
            modulus = 1 << fixture['bits']
            backgrounds = [0, (1 << fixture['width']) - 1]
            backgrounds += [generator.getrandbits(fixture['width']) for _ in range(8)]
            for background in backgrounds:
                initial_sum = sum(_read(background, counter) for counter in fixture['counters'])
                offset = sum((background >> wire & 1) << column
                             for wire, column in fixture['fresh'])
                compressed = _basis_action(background, fixture['compression'])
                self.assertEqual(sum(_read(compressed, output) for output in fixture['outputs'])
                                 % modulus, (initial_sum + offset) % modulus)
                expected_loader = _replace(background, fixture['accumulator'],
                                           (_read(background, fixture['accumulator'])
                                            + initial_sum + offset) % modulus)
                self.assertEqual(_basis_action(background, fixture['loader']), expected_loader)
                self.assertEqual(_basis_action(expected_loader, _adjoint(fixture['loader'])),
                                 background)
                for address in range(1 << count):
                    basis = _replace(background, fixture['address'], address)
                    added = sum((address >> index & 1) == value
                                for index, value in enumerate(fixture['pattern']))
                    expected = _replace(basis, fixture['accumulator'],
                                        (_read(basis, fixture['accumulator']) + added) % modulus)
                    self.assertEqual(_basis_action(basis, fixture['echo']), expected)
                    self.assertEqual(_basis_action(expected, _adjoint(fixture['echo'])), basis)

    def test_column_recurrence_and_weighted_gate_ledger(self):
        profiles = [(count,) * count.bit_length() for count in (3, 4, 7, 8, 15, 31, 63)]
        profiles += [(3, 0, 7, 1), (0, 8, 2, 3, 0), (2, 2, 2, 4)]
        for profile in profiles:
            bits, count = len(profile), max(profile)
            columns = [list(range(sum(profile[:column]), sum(profile[:column + 1])))
                       for column in range(bits)]
            _, _, fresh, history, compressions, _ = _compress_columns(columns, sum(profile))
            self.assertEqual(len(fresh), len({wire for wire, _ in fresh}))
            for previous, following in zip(history, history[1:]):
                triples = [size // 3 for size in previous]
                expected = tuple(previous[column] - 2 * triples[column] + sum(triples[:column])
                                 for column in range(bits))
                self.assertEqual(following, expected)
                # Exact integer version of the proved 4^-j potential.
                scale = 4 ** (bits - 1)
                potential = lambda sizes: sum(max(size - 2, 0) * (scale // 4 ** column)
                                             for column, size in enumerate(sizes))
                self.assertLessEqual(3 * potential(following), 2 * potential(previous))
            self.assertTrue(all(size <= 2 for size in history[-1]))
            for column, used in enumerate(compressions):
                self.assertLessEqual(2 * used, profile[column] + sum(compressions[:column]))
            # Each compressor's two imported INC calls cost O(bits-column).
            # The integer recurrence is checked above; this finite ceiling
            # audits the claimed weighted summation rather than a depth fit.
            weighted = sum(used * (bits - column) for column, used in enumerate(compressions))
            self.assertLessEqual(weighted * 2 ** bits, 3 * count * 3 ** bits)

    def test_native_signed_increment_cycle_needs_the_sign_correction(self):
        counter, helpers = [0, 1], [2, 3]
        addition = _add(helpers, counter)
        flip = _primitive([('X', helpers[0])])
        raw = _serial([_reverse(addition), flip, addition, _reverse(flip)])
        complement = _primitive([('CX', helpers[0], wire) for wire in counter])
        corrected = _serial([complement, raw, _reverse(complement)])
        for basis in range(16):
            sign = -1 if basis >> helpers[0] & 1 else 1
            self.assertEqual(_basis_action(basis, raw[0]),
                             (basis & ~3) | (((basis & 3) + sign) % 4))
            self.assertEqual(_basis_action(basis, corrected[0]),
                             (basis & ~3) | (((basis & 3) + 1) % 4))
        self.assertNotEqual(_basis_action(4, raw[0]), _basis_action(4, corrected[0]))
        images = [(basis & ~3) | (((basis & 3) + 1) % 4) for basis in range(16)]
        actual = _word_matrix(4, _word(corrected[1]))
        np.testing.assert_allclose(actual, np.eye(16)[:, images], atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(actual @ actual - np.eye(16)), 1)
        np.testing.assert_allclose(np.linalg.matrix_power(actual, 4), np.eye(16),
                                   atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
