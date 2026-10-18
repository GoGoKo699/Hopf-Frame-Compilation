"""Bounded evidence for the noncommutative cycle-character construction.

The algebraic claims are analytic. These checks independently compare path
algebra multiplication, in-place factorization, integer phases, and literal
eight-wire H/CNOT/T circuits on every physical input column. Commuting
inversion pivots are also compared with independent increasing-path sums.
No asymptotic native lower bound or arbitrary-table query follows from
these fixtures.
"""
from __future__ import annotations

from collections import Counter
from math import comb
from random import Random
import unittest

import numpy as np

from tests.test_nonlinear_character_transport import _anf, _native_toffoli
from tests.test_parallel_bilinear_word_query import (
    _inverse, _monomial, _native_matrix, _parity_phase, _relative_grid,
)


class _CycleAlgebra:
    """Radical coordinates; integer zero denotes the unit 1 when appropriate."""

    def __init__(self, degree):
        self.degree = degree
        self.basis = [(i, k) for k in range(1, degree + 1)
                      for i in range(degree)]
        self.index = {path: bit for bit, path in enumerate(self.basis)}
        self.width = degree * degree
        self.products = [
            [(self.index[((i + k) % degree, length)],
              self.index[(i, k + length)])
             for length in range(1, degree - k + 1)]
            for i, k in self.basis
        ]
        self.derivative_mask = sum(
            1 << bit for bit, (i, k) in enumerate(self.basis)
            if any((i + step) % degree == 0 for step in range(k))
        )
        self.top_mask = ((1 << degree) - 1) << (degree * (degree - 1))

    def multiply(self, first, second):
        result = 0
        while first:
            bit = (first & -first).bit_length() - 1
            first &= first - 1
            for source, target in self.products[bit]:
                if second >> source & 1:
                    result ^= 1 << target
        return result

    def unit_product(self, first, second):
        return first ^ second ^ self.multiply(first, second)

    def inverse(self, radical):
        term = result = radical
        for _ in range(1, self.degree):
            term = self.multiply(term, radical)
            result ^= term
        return result

    def derivative(self, radical):
        return radical & self.derivative_mask

    def trace(self, radical):
        return (radical & self.top_mask).bit_count() % 2

    def character(self, radical):
        # delta(1+a)(1+a)^-1; the omitted scalar inverse term is handled
        # by initializing result with delta(a).
        term = result = self.derivative(radical)
        for _ in range(1, self.degree):
            term = self.multiply(term, radical)
            result ^= term
        return self.trace(result)


def _factor_gates(algebra, order=None, skip_processed=True):
    order = list(range(algebra.width)) if order is None else list(order)
    rank = {pivot: position for position, pivot in enumerate(order)}
    return [(pivot, source, target) for pivot in order
            for source, target in algebra.products[pivot]
            if not skip_processed or rank[source] > rank[pivot]]


def _apply_factor(value, gates):
    for first, second, target in gates:
        value ^= ((value >> first & 1) & (value >> second & 1)) << target
    return value


def _reconstruct(algebra, factors, order):
    result = 0
    for bit in order:
        if factors >> bit & 1:
            result = algebra.unit_product(result, 1 << bit)
    return result


def _factor_count(degree):
    half = degree // 2
    return (degree * ((degree - 1) ** 2 // 4)
            + half * degree - half * (half + 1) // 2)


def _colored_pivot_batches(algebra):
    """Color each permutation cycle's edges, including parallel 2-cycle edges."""
    degree = algebra.degree
    batches = []
    for length in range(1, degree):
        colors = [[], [], []]
        seen = set()
        for initial in range(degree):
            if initial in seen:
                continue
            cycle = []
            current = initial
            while current not in seen:
                seen.add(current)
                cycle.append(current)
                current = (current + length) % degree
            for position, start in enumerate(cycle):
                color = (2 if len(cycle) % 2 and position == len(cycle) - 1
                         else position % 2)
                colors[color].append(algebra.index[(start, length)])
        batches.extend(batch for batch in colors if batch)
    batches.append([algebra.index[(i, degree)] for i in range(degree)])
    return batches


def _transport_gates(algebra):
    degree = algebra.degree
    gates = [(i, algebra.index[((i + 1) % degree, k - 1)],
              algebra.index[(i, k)])
             for k in range(degree, 1, -1) for i in range(degree)]
    return gates + [(i, None, algebra.index[(i, 1)]) for i in range(degree)]


def _apply_transport(address, work, gates):
    for control, source, target in gates:
        value = address >> control & 1
        if source is not None:
            value &= work >> source & 1
        work ^= value << target
    return work


def _native_source(algebra, work, output):
    factor = [gate for first, second, target in _factor_gates(algebra)
              for gate in _native_toffoli(work[first], work[second], work[target])]
    middle = [gate for i in range(algebra.degree)
              for gate in [('H', work[algebra.index[(i, algebra.degree)]]),
                           ('CX', output, work[algebra.index[(i, algebra.degree)]]),
                           ('H', work[algebra.index[(i, algebra.degree)]])]]
    return factor + middle + _inverse(factor)


def _native_transport(algebra, address, work, helper):
    word = []
    for control, source, target in _transport_gates(algebra):
        if source is None:
            word.append(('CX', address[control], work[target]))
            continue
        relative, _ = _relative_grid([work[source]], [work[target]], [[helper]])
        linear = [('CX', address[control], helper)]
        word += ([('H', work[target])] + _inverse(relative) + linear + relative
                 + _inverse(linear) + [('H', work[target])])
    return word


def _serial_t_depth(word):
    """A valid schedule retaining every intervening Clifford as a barrier."""
    depth = 0
    run = Counter()
    for gate in word + [('BARRIER',)]:
        if gate[0] == 'T':
            run[gate[1]] += 1
        elif run:
            depth += max(run.values())
            run.clear()
    return depth


def _shared_ccz_layers(control, pairs):
    # T^m on the common control is one phase power: at most one T and
    # a Clifford power. Integer-phase evaluation also retains even m.
    layers = [[('T', control, len(pairs))], [], [], []]
    for source, target in pairs:
        layers[0] += [('T', source, 1), ('T', target, 1)]
        layers[1] += _parity_phase(source, (control,), -1)
        layers[1] += _parity_phase(target, (control,), -1)
        layers[2] += _parity_phase(target, (source,), -1)
        layers[3] += _parity_phase(target, (control, source), 1)
    return layers


def _triangular_pairs(degree):
    return ([(i, j) for i in range(degree) for j in range(i + 1, degree)],
            [(i, j) for i in range(degree) for j in range(i + 1)])


def _cycle_from_triangular(algebra, upper, lower):
    upper, lower = int(upper), int(lower)
    upper_pairs, lower_pairs = _triangular_pairs(algebra.degree)
    result = sum(((upper >> bit) & 1) << algebra.index[(i, j - i)]
                 for bit, (i, j) in enumerate(upper_pairs))
    result |= sum(((lower >> bit) & 1) << algebra.index[(i, algebra.degree + j - i)]
                  for bit, (i, j) in enumerate(lower_pairs))
    return result


def _upper_inverse(degree, upper):
    upper_pairs, _ = _triangular_pairs(degree)
    matrix = np.eye(degree, dtype=np.int64)
    for bit, (i, j) in enumerate(upper_pairs):
        matrix[i, j] = upper >> bit & 1
    inverse = np.eye(degree, dtype=np.int64)
    # Independent binary back substitution, not cycle-algebra factorization.
    for i in range(degree - 1, -1, -1):
        for j in range(i + 1, degree):
            inverse[i, j] = sum(matrix[i, k] * inverse[k, j]
                                for k in range(i + 1, j + 1)) % 2
    return matrix, inverse


def _inverse_lower_word(inverse):
    _, lower_pairs = _triangular_pairs(len(inverse))
    return sum(int(inverse[j, i]) << bit for bit, (i, j) in enumerate(lower_pairs))


def _parallel_product_shear(products, helpers):
    """Private-helper schedule for disjoint product targets and source sets."""
    individual = [_relative_grid([second], [target], [[helper]])[1]
                  for (_, second, target), helper in zip(products, helpers)]
    layers = []
    for index in range(4):
        before, phases, after = [], [], []
        for schedule in individual:
            layer = schedule[index]
            position = next(i for i, gate in enumerate(layer) if gate[0] == 'T')
            before += layer[:position]
            phases.append(layer[position])
            after += layer[position + 1:]
        layers.append(before + phases + after)
    relative = [gate for layer in layers for gate in layer]
    linear = [('CX', first, helper)
              for (first, _, _), helper in zip(products, helpers)]
    phase = _inverse(relative) + linear + relative + _inverse(linear)
    hadamards = [('H', target) for _, _, target in products]
    return hadamards + phase + hadamards, phase, layers


def _parallel_lower_shear(control, lower, extra, helpers):
    products = [(control, source, target) for source, target in zip(lower, extra)]
    return _parallel_product_shear(products, helpers)


def _upper_pivot_gates(degree, pivot):
    upper_pairs, _ = _triangular_pairs(degree)
    index = {pair: bit for bit, pair in enumerate(upper_pairs)}
    return [(index[i, pivot], index[pivot, j], index[i, j])
            for i in range(pivot) for j in range(pivot + 1, degree)]


def _upper_path_word(degree, upper, allowed):
    """Enumerate increasing paths directly, independently of pivot updates."""
    upper_pairs, _ = _triangular_pairs(degree)
    index = {pair: bit for bit, pair in enumerate(upper_pairs)}
    result = 0
    for bit, (i, j) in enumerate(upper_pairs):
        internal = [k for k in range(i + 1, j) if k in allowed]
        parity = 0
        for subset in range(1 << len(internal)):
            vertices = [i] + [k for position, k in enumerate(internal)
                              if subset >> position & 1] + [j]
            product = 1
            for first, second in zip(vertices, vertices[1:]):
                product &= upper >> index[first, second] & 1
            parity ^= product
        result |= parity << bit
    return result


class NoncommutativeCycleCharacterTests(unittest.TestCase):
    def test_path_algebra_derivation_cyclic_trace_and_nilpotence(self):
        rng = Random(4107)
        for degree in range(1, 8):
            algebra = _CycleAlgebra(degree)
            basis = [1 << bit for bit in range(algebra.width)]
            for first in basis:
                self.assertEqual(algebra.multiply(first, first), 0)
                for second in basis:
                    product = algebra.multiply(first, second)
                    self.assertEqual(algebra.derivative(product),
                                     algebra.multiply(algebra.derivative(first), second)
                                     ^ algebra.multiply(first, algebra.derivative(second)))
                    self.assertEqual(algebra.trace(product),
                                     algebra.trace(algebra.multiply(second, first)))
            for _ in range(100):
                first, second, third = [rng.getrandbits(algebra.width) for _ in range(3)]
                self.assertEqual(algebra.multiply(algebra.multiply(first, second), third),
                                 algebra.multiply(first, algebra.multiply(second, third)))
            path = 1 << algebra.index[(0, 1)]
            for i in range(1, degree):
                path = algebra.multiply(path, 1 << algebra.index[(i, 1)])
            self.assertEqual(path, 1 << algebra.index[(0, degree)])
            for element in basis:
                self.assertEqual(algebra.multiply(path, element), 0)
                self.assertEqual(algebra.multiply(element, path), 0)

    def test_character_homomorphism_and_two_sided_inverse(self):
        rng = Random(9813)
        for degree in range(1, 10):
            algebra = _CycleAlgebra(degree)
            words = (range(1 << algebra.width) if degree <= 3
                     else [rng.getrandbits(algebra.width) for _ in range(100)])
            for first in words:
                inverse = algebra.inverse(first)
                self.assertEqual(algebra.unit_product(first, inverse), 0)
                self.assertEqual(algebra.unit_product(inverse, first), 0)
                seconds = (range(1 << algebra.width) if degree <= 2
                           else [rng.getrandbits(algebra.width)])
                for second in seconds:
                    self.assertEqual(algebra.character(algebra.unit_product(first, second)),
                                     algebra.character(first) ^ algebra.character(second))

    def test_and_character_has_full_degree_including_even_cycles(self):
        for degree in range(1, 11):
            algebra = _CycleAlgebra(degree)
            # The first d coefficient bits are exactly the d edge coefficients.
            table = [algebra.character(address) for address in range(1 << degree)]
            self.assertEqual(table, [0] * ((1 << degree) - 1) + [1])
            coefficients = _anf(table)
            self.assertEqual([mask for mask, value in enumerate(coefficients) if value],
                             [(1 << degree) - 1])
            self.assertEqual(max(mask.bit_count() for mask, value in enumerate(coefficients)
                                 if value), degree)
        # Both the nilindex and logarithmic coefficient-width extensions fail.
        self.assertGreater(3, (4 - 1).bit_length())
        self.assertGreater(7, (7 * 7).bit_length())

    def test_lexicographic_factorization_exact_count_and_actual_inverse(self):
        rng = Random(2309)
        for degree in range(1, 10):
            algebra = _CycleAlgebra(degree)
            order = list(range(algebra.width))
            gates = _factor_gates(algebra)
            self.assertEqual(len(gates), _factor_count(degree))
            self.assertTrue(all(len({first, second, target}) == 3
                                for first, second, target in gates))
            self.assertTrue(all(second > first and target > first
                                for first, second, target in gates))
            words = (range(1 << algebra.width) if degree <= 3
                     else [rng.getrandbits(algebra.width) for _ in range(100)])
            for word in words:
                factors = _apply_factor(word, gates)
                self.assertEqual(_apply_factor(factors, list(reversed(gates))), word)
                self.assertEqual(_reconstruct(algebra, factors, order), word)
                self.assertEqual(algebra.trace(factors), algebra.character(word))
        algebra = _CycleAlgebra(2)
        wrong = _factor_gates(algebra, skip_processed=False)
        # Stored earlier factor bits must not be read as residual coefficients.
        self.assertTrue(any(_reconstruct(algebra, _apply_factor(word, wrong), range(4)) != word
                            for word in range(16)))

    def test_transport_is_left_multiplication_and_uses_actual_inverse(self):
        for degree in range(1, 4):
            algebra = _CycleAlgebra(degree)
            gates = _transport_gates(algebra)
            self.assertEqual(sum(source is not None for _, source, _ in gates),
                             degree * (degree - 1))
            self.assertEqual(sum(source is None for _, source, _ in gates), degree)
            differs_from_forward = differs_from_right = False
            for address in range(1 << degree):
                inverse = algebra.inverse(address)
                for word in range(1 << algebra.width):
                    actual = _apply_transport(address, word, gates)
                    self.assertEqual(actual, algebra.unit_product(address, word))
                    backward = _apply_transport(address, word, list(reversed(gates)))
                    self.assertEqual(backward, algebra.unit_product(inverse, word))
                    self.assertEqual(_apply_transport(address, actual, list(reversed(gates))), word)
                    differs_from_forward |= actual != backward
                    differs_from_right |= actual != algebra.unit_product(word, address)
            if degree >= 2:
                self.assertTrue(differs_from_forward)
                self.assertTrue(differs_from_right)

    def test_controlled_character_commutator_on_every_small_input(self):
        for degree in range(1, 4):
            algebra = _CycleAlgebra(degree)
            gates = _transport_gates(algebra)
            factors = _factor_gates(algebra)
            phases = [algebra.trace(_apply_factor(word, factors))
                      for word in range(1 << algebra.width)]
            for address in range(1 << degree):
                for work in range(1 << algebra.width):
                    for output in range(2):
                        current = _apply_transport(address, work, list(reversed(gates)))
                        phase = -4 * output * phases[current]
                        current = _apply_transport(address, current, gates)
                        phase += 4 * output * phases[current]
                        self.assertEqual(current, work)
                        self.assertEqual(phase % 8,
                                         4 * output * int(address == (1 << degree) - 1))

    def test_complete_native_eight_wire_query_and_exact_ledger(self):
        algebra = _CycleAlgebra(2)
        address, work, output, helper = [0, 1], [2, 3, 4, 5], 6, 7
        source = _native_source(algebra, work, output)
        transport = _native_transport(algebra, address, work, helper)
        hadamards = [('H', output)]
        word = (hadamards + _inverse(transport) + _inverse(source)
                + transport + source + hadamards)
        self.assertEqual((sum(g[0] == 'T' for g in source),
                          sum(g[0] != 'T' for g in source)), (14, 30))
        self.assertEqual((sum(g[0] == 'T' for g in transport),
                          sum(g[0] != 'T' for g in transport)), (16, 42))
        self.assertEqual(sum(g[0] == 'T' for g in word), 60)
        self.assertEqual(sum(g[0] != 'T' for g in word), 146)
        self.assertEqual(_serial_t_depth(source), 10)
        self.assertEqual(_serial_t_depth(transport), 16)
        self.assertLessEqual(_serial_t_depth(word), 52)
        self.assertTrue(all(g[1] not in address + [output] for g in word if g[0] == 'T'))
        self.assertTrue(all(g[-1] not in address for g in word if g[0] == 'CX'))
        expected = np.zeros((256, 256), dtype=complex)
        for column in range(256):
            result = column ^ (int(column & 3 == 3) << output)
            expected[result, column] = 1
        np.testing.assert_allclose(_native_matrix(word, 8), expected, rtol=0, atol=5e-14)
        np.testing.assert_allclose(_native_matrix(_inverse(word), 8),
                                   expected.conj().T, rtol=0, atol=5e-14)

    def test_four_layer_shared_control_ccz_on_every_column(self):
        for count in range(1, 6):
            pairs = [(1 + 2 * i, 2 + 2 * i) for i in range(count)]
            layers = _shared_ccz_layers(0, pairs)
            self.assertEqual(len(layers), 4)
            for layer in layers:
                targets = [gate[1] for gate in layer if gate[0] == 'T']
                self.assertEqual(len(targets), len(set(targets)))
            self.assertEqual(layers[0][0], ('T', 0, count))
            word = [gate for layer in layers for gate in layer]
            width = 2 * count + 1
            values = np.arange(1 << width)
            images, phases = _monomial(word, width)
            expected = np.zeros_like(values)
            for first, second in pairs:
                expected += 4 * (values & 1) * ((values >> first) & 1) * ((values >> second) & 1)
            np.testing.assert_array_equal(images, values)
            np.testing.assert_array_equal(phases, expected % 8)
            returned, inverse_phase = _monomial(word + _inverse(word), width)
            np.testing.assert_array_equal(returned, values)
            np.testing.assert_array_equal(inverse_phase, np.zeros_like(values))
            if count == 2:
                hadamards = [('H', target) for _, target in pairs]
                actual = _native_matrix(hadamards + word + hadamards, width)
                targets, _ = _monomial([('CCX', 0, source, target)
                                        for source, target in pairs], width)
                expected_matrix = np.zeros_like(actual)
                expected_matrix[targets, values] = 1
                np.testing.assert_allclose(actual, expected_matrix, rtol=0, atol=2e-14)

    def test_colored_factorization_matching_schedule_and_reconstruction(self):
        rng = Random(7291)
        for degree in range(2, 10):
            algebra = _CycleAlgebra(degree)
            batches = _colored_pivot_batches(algebra)
            order = [pivot for batch in batches for pivot in batch]
            self.assertEqual(sorted(order), list(range(algebra.width)))
            lengths = [algebra.basis[pivot][1] for pivot in order]
            self.assertEqual(lengths, sorted(lengths))
            self.assertLessEqual(len(batches) - 1, 3 * (degree - 1))
            for batch in batches[:-1]:
                endpoints = []
                for pivot in batch:
                    start, length = algebra.basis[pivot]
                    endpoints.extend([start, (start + length) % degree])
                self.assertEqual(len(endpoints), len(set(endpoints)))
            gates = _factor_gates(algebra, order)
            self.assertLessEqual(len(gates), degree * ((degree - 1) ** 2 // 4)
                                 + degree * (degree // 2))
            self.assertTrue(all(len(set(gate)) == 3 for gate in gates))
            for _ in range(100):
                word = rng.getrandbits(algebra.width)
                factors = _apply_factor(word, gates)
                self.assertEqual(_reconstruct(algebra, factors, order), word)
                self.assertEqual(_apply_factor(factors, list(reversed(gates))), word)
                self.assertEqual(algebra.trace(factors), algebra.character(word))

    def test_source_is_upper_inverse_lower_pairing(self):
        rng = Random(3901)
        for degree in range(1, 9):
            algebra = _CycleAlgebra(degree)
            upper_pairs, lower_pairs = _triangular_pairs(degree)
            upper_width, lower_width = len(upper_pairs), len(lower_pairs)
            samples = (((upper, lower) for upper in range(1 << upper_width)
                        for lower in range(1 << lower_width)) if degree <= 3
                       else [(rng.getrandbits(upper_width), rng.getrandbits(lower_width))
                             for _ in range(200)])
            for upper, lower in samples:
                matrix, inverse = _upper_inverse(degree, upper)
                identity = np.eye(degree, dtype=np.int64)
                np.testing.assert_array_equal(matrix @ inverse % 2, identity)
                np.testing.assert_array_equal(inverse @ matrix % 2, identity)
                lower_matrix = np.zeros_like(matrix)
                for bit, (i, j) in enumerate(lower_pairs):
                    lower_matrix[i, j] = lower >> bit & 1
                expected = int(np.trace(inverse @ lower_matrix) % 2)
                word = _cycle_from_triangular(algebra, upper, lower)
                self.assertEqual(algebra.character(word), expected)
                self.assertEqual(expected, (lower & _inverse_lower_word(inverse)).bit_count() % 2)

    def test_fourier_source_is_complete_inverse_word_query(self):
        # U is preserved by the completed diagonal source. This says nothing
        # about gate-level read-only U access inside a native implementation.
        hadamard = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
        for degree in range(1, 4):
            algebra = _CycleAlgebra(degree)
            upper_pairs, lower_pairs = _triangular_pairs(degree)
            lower_width = len(lower_pairs)
            fourier = np.array([[1.]])
            for _ in range(lower_width):
                fourier = np.kron(fourier, hadamard)
            for upper in range(1 << len(upper_pairs)):
                _, inverse = _upper_inverse(degree, upper)
                inverse_word = _inverse_lower_word(inverse)
                source = np.array([algebra.character(_cycle_from_triangular(algebra, upper, lower))
                                   for lower in range(1 << lower_width)])
                diagonal = 1 - 2 * source
                actual = (fourier * diagonal) @ fourier
                expected = np.zeros_like(actual)
                for column in range(1 << lower_width):
                    expected[column ^ inverse_word, column] = 1
                np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-14)
                # H D^dagger H is the actual reversed ideal word. The literal
                # native inverse is checked separately in the eight-wire test.
                actual_inverse = (fourier * diagonal.conj()) @ fourier
                np.testing.assert_allclose(actual_inverse, expected.conj().T,
                                           rtol=0, atol=2e-14)

    def test_dirty_lower_array_echo_controls_the_uncontrolled_source(self):
        for degree in range(1, 4):
            algebra = _CycleAlgebra(degree)
            upper_pairs, lower_pairs = _triangular_pairs(degree)
            lower_width = len(lower_pairs)
            lower_values = np.arange(1 << lower_width, dtype=np.int64)
            offsets = lower_values.copy()
            missing_return_detected = missing_offset_cancellation_detected = False
            for upper in range(1 << len(upper_pairs)):
                source = np.array([algebra.character(_cycle_from_triangular(algebra, upper, extra))
                                   for extra in offsets])
                for lower in lower_values:
                    expected_phase = source[lower]
                    for output in range(2):
                        # Chronological D(U,E), S_h, D(U,E)^dagger, S_h^dagger.
                        phase = 4 * source[offsets]
                        current = offsets ^ (lower if output else 0)
                        phase -= 4 * source[current]
                        before_return = current.copy()
                        current ^= lower if output else 0
                        np.testing.assert_array_equal(current, offsets)
                        np.testing.assert_array_equal(phase % 8,
                                                       np.full_like(offsets, 4 * output * expected_phase))
                        missing_return_detected |= bool(np.any(before_return != offsets))
                        # Omitting the first source leaves the unknown E offset
                        # in the phase; omitting the last shear leaves E changed.
                        wrong_phase = (-4 * source[before_return]) % 8
                        missing_offset_cancellation_detected |= bool(np.any(
                            wrong_phase != 4 * output * expected_phase))
            self.assertTrue(missing_return_detected)
            self.assertTrue(missing_offset_cancellation_detected)

    def test_single_inverse_entry_xor_needs_no_controlled_source(self):
        hadamard = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
        for degree in range(1, 4):
            algebra = _CycleAlgebra(degree)
            upper_pairs, lower_pairs = _triangular_pairs(degree)
            for upper in range(1 << len(upper_pairs)):
                _, inverse = _upper_inverse(degree, upper)
                source = [algebra.character(_cycle_from_triangular(algebra, upper, lower))
                          for lower in range(1 << len(lower_pairs))]
                for bit, (row, column) in enumerate(lower_pairs):
                    entry = int(inverse[column, row])
                    expected = np.eye(2)[:, [entry, 1 ^ entry]]
                    for lower in range(1 << len(lower_pairs)):
                        phases = []
                        for output in range(2):
                            # D, CX(h->L_bit), D^dagger, CX(h->L_bit).
                            # Only the translation is controlled: no extra E.
                            phase = 4 * source[lower]
                            current = lower ^ (output << bit)
                            phase -= 4 * source[current]
                            current ^= output << bit
                            self.assertEqual(current, lower)
                            self.assertEqual(phase % 8, 4 * output * entry)
                            phases.append(1 if phase % 8 == 0 else -1)
                        actual = (hadamard * np.array(phases)) @ hadamard
                        np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-14)
        # The added word is two CX gates and two output Hadamards. Source
        # and actual-source-inverse costs must be charged separately in general.

    def test_commuting_pivots_match_each_partial_path_set(self):
        # Exhaust every upper word and pivot subset through dimension four.
        for degree in range(1, 5):
            pairs, _ = _triangular_pairs(degree)
            pivots = list(range(1, degree - 1))
            for upper in range(1 << len(pairs)):
                images = {}
                for subset in range(1 << len(pivots)):
                    chosen = [pivot for bit, pivot in enumerate(pivots)
                              if subset >> bit & 1]
                    gates = [gate for pivot in chosen
                             for gate in _upper_pivot_gates(degree, pivot)]
                    image = _apply_factor(upper, gates)
                    self.assertEqual(image, _upper_path_word(degree, upper, chosen))
                    reverse_order = [gate for pivot in reversed(chosen)
                                     for gate in _upper_pivot_gates(degree, pivot)]
                    self.assertEqual(_apply_factor(upper, reverse_order), image)
                    images[subset] = image
                for first in images:
                    for second in images:
                        gates = [gate for bit, pivot in enumerate(pivots)
                                 if second >> bit & 1
                                 for gate in _upper_pivot_gates(degree, pivot)]
                        self.assertEqual(_apply_factor(images[first], gates),
                                         images[first ^ second])

    def test_commuting_pivot_inversion_and_actual_reverse(self):
        rng = Random(68271)
        for degree in range(1, 9):
            pairs, _ = _triangular_pairs(degree)
            batches = [_upper_pivot_gates(degree, pivot)
                       for pivot in range(1, degree - 1)]
            gates = [gate for batch in batches for gate in batch]
            self.assertEqual(len(gates), comb(degree, 3))
            self.assertEqual(max(map(len, batches), default=0),
                             (degree - 1) ** 2 // 4)
            for batch in batches:
                sources = {wire for first, second, _ in batch
                           for wire in (first, second)}
                targets = [target for _, _, target in batch]
                self.assertEqual(len(targets), len(set(targets)))
                self.assertTrue(sources.isdisjoint(targets))
            samples = (range(1 << len(pairs)) if degree <= 5
                       else [rng.getrandbits(len(pairs)) for _ in range(100)])
            for upper in samples:
                image = _apply_factor(upper, gates)
                _, inverse = _upper_inverse(degree, upper)
                expected = sum(int(inverse[i, j]) << bit
                               for bit, (i, j) in enumerate(pairs))
                self.assertEqual(image, expected)
                self.assertEqual(_apply_factor(image, list(reversed(gates))), upper)
                self.assertEqual(_apply_factor(image, gates), upper)
                for first in batches:
                    self.assertEqual(_apply_factor(_apply_factor(upper, first), first), upper)
                    for second in batches:
                        self.assertEqual(_apply_factor(upper, first + second),
                                         _apply_factor(upper, second + first))

    def test_simultaneous_quadratic_updates_miss_a_cubic_path(self):
        degree = 4
        pairs, _ = _triangular_pairs(degree)
        index = {pair: bit for bit, pair in enumerate(pairs)}
        chain = sum(1 << index[i, i + 1] for i in range(degree - 1))
        gates = [gate for pivot in range(1, degree - 1)
                 for gate in _upper_pivot_gates(degree, pivot)]
        simultaneous = chain
        for first, second, target in gates:
            simultaneous ^= ((chain >> first & 1) & (chain >> second & 1)) << target
        actual = _apply_factor(chain, gates)
        self.assertEqual(simultaneous >> index[0, 3] & 1, 0)
        self.assertEqual(actual >> index[0, 3] & 1, 1)
        self.assertNotEqual(simultaneous, actual)
        # A cross-pivot target becomes a control, so the individual updates
        # do not commute despite commutation of the completed pivots.
        first = (index[0, 1], index[1, 2], index[0, 2])
        second = (index[0, 2], index[2, 3], index[0, 3])
        self.assertNotEqual(_apply_factor(chain, [first, second]),
                            _apply_factor(chain, [second, first]))

    def test_native_pivot_parallel_phase_and_full_inversion_columns(self):
        for degree in range(3, 6):
            upper_pairs, _ = _triangular_pairs(degree)
            for pivot in range(1, degree - 1):
                products = _upper_pivot_gates(degree, pivot)
                count = len(products)
                helpers = list(range(len(upper_pairs), len(upper_pairs) + count))
                word, phase, layers = _parallel_product_shear(products, helpers)
                self.assertEqual(sum(gate[0] == 'T' for gate in word), 8 * count)
                self.assertEqual(sum(gate[0] != 'T' for gate in word), 20 * count)
                self.assertEqual(_serial_t_depth(word), 8)
                self.assertTrue(all(gate[1] in helpers for gate in word if gate[0] == 'T'))
                for layer in layers:
                    targets = [gate[1] for gate in layer if gate[0] == 'T']
                    self.assertEqual(set(targets), set(helpers))
                    self.assertEqual(len(targets), len(set(targets)))
                width = len(upper_pairs) + count
                values = np.arange(1 << width, dtype=np.int64)
                images, exponents = _monomial(phase, width)
                expected = np.zeros_like(values)
                for first, second, target in products:
                    expected ^= ((values >> first & 1) & (values >> second & 1)
                                 & (values >> target & 1))
                np.testing.assert_array_equal(images, values)
                np.testing.assert_array_equal(exponents, 4 * expected)

        # The smallest two-pivot inversion has six arbitrary matrix bits and
        # two reused arbitrary helpers. Check all 256 physical columns,
        # literal phases, both pivot orders and the gate-by-gate inverse.
        degree, upper_width, width = 4, 6, 8
        pairs, _ = _triangular_pairs(degree)
        expected = np.zeros((1 << width, 1 << width), dtype=complex)
        for value in range(1 << width):
            _, inverse = _upper_inverse(degree, value & ((1 << upper_width) - 1))
            image = sum(int(inverse[i, j]) << bit
                        for bit, (i, j) in enumerate(pairs))
            image |= value & ~((1 << upper_width) - 1)
            expected[image, value] = 1
        for order in ((1, 2), (2, 1)):
            word = [gate for pivot in order
                    for gate in _parallel_product_shear(
                        _upper_pivot_gates(degree, pivot), [6, 7])[0]]
            self.assertEqual(sum(gate[0] == 'T' for gate in word), 32)
            self.assertEqual(sum(gate[0] != 'T' for gate in word), 80)
            self.assertEqual(_serial_t_depth(word), 16)
            np.testing.assert_allclose(_native_matrix(word, width), expected,
                                       rtol=0, atol=3e-14)
            np.testing.assert_allclose(_native_matrix(_inverse(word), width),
                                       expected.conj().T, rtol=0, atol=3e-14)

    def test_pivot_compute_pair_uncompute_is_cycle_source(self):
        rng = Random(2807)
        for degree in range(1, 8):
            algebra = _CycleAlgebra(degree)
            upper_pairs, lower_pairs = _triangular_pairs(degree)
            lower_index = {pair: bit for bit, pair in enumerate(lower_pairs)}
            gates = [gate for pivot in range(1, degree - 1)
                     for gate in _upper_pivot_gates(degree, pivot)]
            samples = (((upper, lower) for upper in range(1 << len(upper_pairs))
                        for lower in range(1 << len(lower_pairs))) if degree <= 3
                       else [(rng.getrandbits(len(upper_pairs)),
                              rng.getrandbits(len(lower_pairs))) for _ in range(100)])
            for upper, lower in samples:
                image = _apply_factor(upper, gates)
                phase = sum(lower >> lower_index[i, i] & 1 for i in range(degree)) % 2
                for bit, (i, j) in enumerate(upper_pairs):
                    phase ^= (image >> bit & 1) & (lower >> lower_index[j, i] & 1)
                self.assertEqual(phase, algebra.character(
                    _cycle_from_triangular(algebra, upper, lower)))
                self.assertEqual(_apply_factor(image, list(reversed(gates))), upper)

    def test_private_lower_shear_native_word_parallel_depth_and_actual_inverse(self):
        for count in range(1, 4):
            control = 0
            lower = list(range(1, count + 1))
            extra = list(range(count + 1, 2 * count + 1))
            helpers = list(range(2 * count + 1, 3 * count + 1))
            word, phase, layers = _parallel_lower_shear(control, lower, extra, helpers)
            self.assertEqual(sum(gate[0] == 'T' for gate in word), 8 * count)
            self.assertEqual(sum(gate[0] != 'T' for gate in word), 20 * count)
            self.assertEqual(_serial_t_depth(word), 8)
            self.assertTrue(all(gate[1] in helpers for gate in word if gate[0] == 'T'))
            for layer in layers:
                targets = [gate[1] for gate in layer if gate[0] == 'T']
                self.assertEqual(set(targets), set(helpers))
                self.assertEqual(len(targets), len(set(targets)))
            width = 3 * count + 1
            values = np.arange(1 << width)
            images, exponents = _monomial(phase, width)
            expected = np.zeros_like(values)
            for first, second in zip(lower, extra):
                expected ^= ((values >> first) & 1) & ((values >> second) & 1)
            expected *= 4 * (values & 1)
            np.testing.assert_array_equal(images, values)
            np.testing.assert_array_equal(exponents, expected)
            returned, cancelled = _monomial(phase + _inverse(phase), width)
            np.testing.assert_array_equal(returned, values)
            np.testing.assert_array_equal(cancelled, np.zeros_like(values))
            if count == 1:
                native = _native_matrix(word, width)
                expected_matrix = np.zeros_like(native)
                for value in values:
                    target = value ^ (((value & 1) & ((value >> lower[0]) & 1)) << extra[0])
                    expected_matrix[target, value] = 1
                np.testing.assert_allclose(native, expected_matrix, rtol=0, atol=2e-14)
                np.testing.assert_allclose(_native_matrix(_inverse(word), width),
                                           expected_matrix.conj().T, rtol=0, atol=2e-14)


if __name__ == '__main__':
    unittest.main()
