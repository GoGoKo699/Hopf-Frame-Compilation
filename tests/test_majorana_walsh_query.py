"""Exact fixtures for a paid Walsh/Majorana probe and its restricted ledger.

Native checks retain all carrier/address columns, phases and actual
inverses. Binary routing checks do not claim an unrestricted depth bound;
the proof supplies the arbitrary-size Majorana-linear lower bounds.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CO, CZ, add, adjoint, eye, native_columns, scale, zeros
from tests.test_two_layer_diagonal_rigidity import _primitive_inverse, _reduce_adjacent_s
from tests.test_two_layer_galois_return import (
    _coefficients, _conjugate, _product, _symplectic, _tau,
)


def _inverse(word):
    return _reduce_adjacent_s(_primitive_inverse(word))


def _majoranas(width):
    return [(1 << wire) | (((1 << wire) - 1 | (kind << wire)) << width)
            for wire in range(width) for kind in (0, 1)]


def _pauli_word(label, width):
    word = []
    for wire in range(width):
        x, z = (label >> wire) & 1, (label >> (width + wire)) & 1
        flip = [('H', wire), ('S', wire), ('S', wire), ('H', wire)]
        if x and z:
            word += [('S', wire)] * 3 + flip + [('S', wire)]
        elif x:
            word += flip
        elif z:
            word += [('S', wire), ('S', wire)]
    return word


def _mask_label(p, bit):
    label = 0
    for mode, gamma in enumerate(_majoranas(1 << (p - 1))):
        if mode >> bit & 1:
            label ^= gamma
    return label


def _adjacent_swap(mode):
    wire = mode // 2
    if mode % 2 == 0:
        return [('S', wire)]
    # This signed swap sends gamma_mode to gamma_(mode+1) and
    # gamma_(mode+1) to -gamma_mode, with its literal native scalar.
    return [('H', wire), ('H', wire + 1), ('CX', wire, wire + 1),
            ('S', wire + 1), ('CX', wire, wire + 1),
            ('H', wire), ('H', wire + 1)]


def _matching_route(p, bit):
    m, distance = 1 << p, 1 << bit
    desired = [mode for start in range(0, m, 2 * distance)
               for left in range(start, start + distance)
               for mode in (left, left + distance)]
    modes, signs, swaps, word = list(range(m)), [1] * m, [], []
    for destination, mode in enumerate(desired):
        current = modes.index(mode)
        while current > destination:
            position = current - 1
            modes[position], modes[current] = modes[current], modes[position]
            signs[position], signs[current] = -signs[current], signs[position]
            swaps.append(position)
            word += _adjacent_swap(position)
            current -= 1
    return word, modes, signs, swaps


def _walsh_schedule(p):
    width = 1 << (p - 1)
    schedule = []
    for bit in range(p):
        route, _, signs, _ = _matching_route(p, bit)
        band = [('T' if signs[2 * wire] * signs[2 * wire + 1] == 1 else 'TDG', wire)
                for wire in range(width)]
        schedule += [('C', _pauli_word(_mask_label(p, bit), width) + route),
                     ('T', band), ('C', _inverse(route))]
    return schedule


def _flat_schedule(p):
    schedule = []
    for bit in range(p):
        route, _, signs, _ = _matching_route(p, bit)
        band = [('T' if signs[2 * wire] * signs[2 * wire + 1] == 1 else 'TDG', wire)
                for wire in range(1 << bit)]
        schedule += [('C', route), ('T', band), ('C', _inverse(route))]
    return schedule


def _word(schedule):
    return [gate for _, gates in schedule for gate in gates]


def _shift(word, offset):
    return [(name, *(wire + offset for wire in wires)) for name, *wires in word]


def _controlled_pauli(control, label, width, offset):
    word = []
    for wire in range(width):
        x, z = (label >> wire) & 1, (label >> (width + wire)) & 1
        target = offset + wire
        if x and z:
            word += [('S', target)] * 3 + [('CX', control, target), ('S', target)]
        elif x:
            word += [('CX', control, target)]
        elif z:
            word += [('H', target), ('CX', control, target), ('H', target)]
    return word


class MajoranaWalshQueryTests(unittest.TestCase):
    def test_routed_matchings_have_independent_native_bands_and_quadratic_ledger(self):
        for p in range(2, 7):
            m, width = 1 << p, 1 << (p - 1)
            integer_matrix = [[int(row == column) for column in range(m)]
                              for row in range(m)]
            total_swaps = 0
            for bit in range(p):
                distance = 1 << bit
                route, modes, signs, swaps = _matching_route(p, bit)
                total_swaps += len(swaps)
                self.assertEqual(len(swaps), m * (distance - 1) // 4)
                self.assertLessEqual(len(route), 9 * len(swaps))
                self.assertLessEqual(len(_inverse(route)), 9 * len(swaps))
                positions = {mode: position for position, mode in enumerate(modes)}
                for mode in range(m):
                    position = positions[mode]
                    orientation = signs[position & ~1] * signs[(position & ~1) + 1]
                    input_sign = signs[position] * (-1 if mode & distance else 1)
                    # Unscaled T/TDG plane, then the actual signed unroute.
                    image = {modes[position]: input_sign * signs[position],
                             modes[position ^ 1]: input_sign * signs[position ^ 1]
                             * orientation * (-1 if position & 1 else 1)}
                    self.assertEqual(image, {mode: -1 if mode & distance else 1,
                                             mode ^ distance: 1})
                old = integer_matrix
                integer_matrix = [[0] * m for _ in range(m)]
                for start in range(0, m, 2 * distance):
                    for left in range(start, start + distance):
                        right = left + distance
                        integer_matrix[left] = [a + b for a, b in zip(old[left], old[right])]
                        integer_matrix[right] = [a - b for a, b in zip(old[left], old[right])]
            self.assertEqual(total_swaps, m * (m - 1 - p) // 4)
            self.assertEqual(integer_matrix,
                             [[(-1) ** ((row & column).bit_count())
                               for column in range(m)] for row in range(m)])
            schedule = _walsh_schedule(p)
            bands = [gates for kind, gates in schedule if kind == 'T']
            self.assertEqual(len(bands), p)
            self.assertTrue(all(len(band) == width
                                and {gate[1] for gate in band} == set(range(width))
                                for band in bands))
            word = _word(schedule)
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), m * p // 2)
            clifford_bound = 18 * total_swaps + 4 * m * p
            self.assertLess(2 * clifford_bound, 9 * m * m)
            for orientation in (word, _inverse(word)):
                self.assertLessEqual(sum(g[0] not in ('T', 'TDG') for g in orientation),
                                     clifford_bound)
            flat_schedule = _flat_schedule(p)
            flat_word = _word(flat_schedule)
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in flat_word), m - 1)
            self.assertEqual([len(gates) for kind, gates in flat_schedule if kind == 'T'],
                             [1 << bit for bit in range(p)])
            for orientation in (flat_word, _inverse(flat_word)):
                self.assertLessEqual(sum(g[0] not in ('T', 'TDG') for g in orientation),
                                     18 * total_swaps)

    def test_two_qubit_walsh_has_literal_native_majorana_action_and_inverse(self):
        word = _word(_walsh_schedule(2))
        actual = native_columns(2, word)
        gamma = _majoranas(2)
        self.assertEqual(gamma, [1, 5, 6, 14])
        for column, label in enumerate(gamma):
            coefficients = _coefficients(_conjugate(actual, _tau(label, 2)), 2)
            expected = [CZ] * 16
            for row, output in enumerate(gamma):
                expected[output] = CO * Fraction((-1) ** ((row & column).bit_count()), 2)
            self.assertEqual(coefficients, expected)
        inverse = _inverse(word)
        self.assertEqual(native_columns(2, inverse), adjoint(actual))
        for orientation in (word, inverse):
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in orientation), 4)
            self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in orientation), 30)

    def test_complete_probe_preserves_all_inputs_and_needs_actual_mask_inverse(self):
        p, width = 2, 2
        gamma = _majoranas(width)
        first, second = _mask_label(p, 0), _mask_label(p, 1)
        self.assertEqual(_symplectic(first, second, width), 1)
        self.assertEqual(_product(_tau(first, width), _tau(second, width)),
                         scale(_product(_tau(second, width), _tau(first, width)), -1))
        walsh = _shift(_word(_walsh_schedule(p)), p)
        walsh_inverse = _inverse(walsh)
        mask = [_gate for bit in range(p)
                for _gate in _controlled_pauli(bit, _mask_label(p, bit), width, p)]
        seed = _shift(_pauli_word(gamma[0], width), p)
        word = walsh_inverse + _inverse(mask) + walsh + seed + walsh_inverse + mask + walsh
        actual = native_columns(p + width, word)
        for column in range(16):
            address, carrier = column & 3, column >> p
            probe = _tau(gamma[address], width)
            for row in range(16):
                expected = probe[row >> p][carrier] if row & 3 == address else CZ
                self.assertEqual(actual[row][column], expected)
        inverse = _inverse(word)
        self.assertEqual(native_columns(p + width, inverse), adjoint(actual))
        for orientation in (word, inverse):
            self.assertEqual(sum(g[0] in ('T', 'TDG') for g in orientation), 16)
            self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in orientation), 142)
        wrong = walsh_inverse + mask + walsh + seed + walsh_inverse + mask + walsh
        wrong_matrix = native_columns(p + width, wrong)
        for column in range(16):
            for row in range(16):
                self.assertEqual(wrong_matrix[row][column],
                                 actual[row][column] * (-1 if column & 3 == 3 else 1))
        self.assertNotEqual(wrong_matrix, actual)

    def test_flat_frame_equality_phase_and_reflection_keep_non_pauli_duals(self):
        width, p = 2, 2
        labels = _majoranas(width)
        gammas = [_tau(label, width) for label in labels]
        flat = zeros(4)
        for gamma in gammas:
            flat = add(flat, scale(gamma, Fraction(1, 2)))
        self.assertEqual(_product(flat, flat), eye(4))
        preparation_word = _word(_flat_schedule(p))
        preparation = native_columns(width, preparation_word)
        self.assertEqual(_conjugate(preparation, gammas[0]), flat)
        preparation_inverse = _inverse(preparation_word)
        self.assertEqual(native_columns(width, preparation_inverse), adjoint(preparation))
        flat_word = preparation_inverse + _pauli_word(labels[0], width) + preparation_word
        self.assertEqual(native_columns(width, flat_word), flat)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in flat_word), 6)
        self.assertLessEqual(sum(g[0] not in ('T', 'TDG') for g in flat_word), 40)
        walsh = native_columns(width, _word(_walsh_schedule(p)))
        frame = []
        for address in range(4):
            mask_word = [_gate for bit in range(p) if address >> bit & 1
                         for _gate in _pauli_word(_mask_label(p, bit), width)]
            member = _conjugate(native_columns(width, mask_word), flat)
            self.assertEqual(member, _conjugate(walsh, gammas[address]))
            self.assertEqual(sum(value != CZ for value in _coefficients(member, width)), 4)
            frame.append(member)
        for left in range(4):
            for right in range(4):
                anticommutator = add(_product(frame[left], frame[right]),
                                     _product(frame[right], frame[left]))
                self.assertEqual(anticommutator, scale(eye(4), 2 if left == right else 0))
                self.assertEqual(_product(frame[left], frame[right], frame[left], frame[right]),
                                 scale(eye(4), 1 if left == right else -1))

        all_labels = 0
        for label in labels:
            all_labels ^= label
        for index, gamma in enumerate(gammas):
            # At m=4, 2/sqrt(m)=1. The difference has operator norm one,
            # witnessed exactly by its squared absolute value.
            self.assertEqual(_product(flat, gamma, flat), add(scale(gamma, -1), flat))
            pauli = gammas[(index + 1) % 4]
            changed = add(_conjugate(_product(pauli, flat), gamma), scale(gamma, -1))
            self.assertEqual(changed, _conjugate(pauli, flat))
            self.assertEqual(_product(adjoint(changed), changed), eye(4))

            dual = _conjugate(walsh, _tau(all_labels ^ labels[index], width))
            self.assertEqual(sum(value != CZ for value in _coefficients(dual, width)), 4)
            for other, member in enumerate(frame):
                self.assertEqual(_conjugate(dual, member),
                                 scale(member, -1 if index == other else 1))


if __name__ == '__main__':
    unittest.main()
