"""Bounded evidence for uniform transport through controlled Pauli sources.

The proof supplies the width-independent degree theorem. These fixtures
retain Hermitian-section carries, the known Gottesman--Mochon example,
its actual native inverse, and paid tensor-copy allocation. They do not
assert a three-native-band realization or construct a returned query.
Seven-qubit Pauli-domain checks use monomial bitmaps, not a dense Pauli
basis; the only full matrices are the 128 native input columns.
"""
from __future__ import annotations

from collections import Counter
import itertools
import unittest

import numpy as np

from tests.test_nonlinear_character_transport import _anf
from tests.test_outer_axis_normalizer import _tau
from tests.test_parallel_bilinear_word_query import _inverse as _native_inverse
from tests.test_parallel_bilinear_word_query import _native_matrix
from tests.test_read_only_three_layer_rigidity import H, _wire
from tests.test_returned_diagonal_echo_degree import _ccz


def _degree(values):
    return max((mask.bit_count() for mask, bit in enumerate(_anf(values)) if bit),
               default=-1)


def _fg(control):
    z = [(control >> j) & 1 for j in range(4)]
    first = z[0] * z[1] | (z[0] * z[2] << 1) | (z[0] * z[3] << 2)
    second = z[2] * z[3] | (z[1] * z[3] << 1) | (z[1] * z[2] << 2)
    return first, second


def _gm():
    """Literal X^f Z^g, controls on wires 0..3 and targets on 4..6."""
    permutation, phase = [], []
    for value in range(128):
        f, g = _fg(value & 15)
        permutation.append(value ^ (f << 4))
        phase.append(4 * ((g & (value >> 4)).bit_count() % 2))
    return np.array(permutation), np.array(phase)


def _compose(left, right):
    """Composition left after right; phases are integers modulo eight."""
    return left[0][right[0]], (right[1] + left[1][right[0]]) % 8


def _adjoint(operator):
    inverse = np.argsort(operator[0])
    return inverse, -operator[1][inverse] % 8


def _pauli(label, width):
    mask = (1 << width) - 1
    x, z = label & mask, label >> width
    values = np.arange(1 << width)
    phase = np.array([(2 * (x & z).bit_count()
                       + 4 * (z & int(value)).bit_count()) % 8 for value in values])
    return values ^ x, phase


def _pauli_label(operator, width):
    """Return the actual signed Hermitian label, or None."""
    values = np.arange(1 << width)
    x = int(operator[0][0])
    if not np.array_equal(operator[0], values ^ x):
        return None
    relative = (operator[1] - operator[1][0]) % 8
    if np.any((relative != 0) & (relative != 4)):
        return None
    z = sum(int(relative[1 << j] == 4) << j for j in range(width))
    expected = np.array([4 * ((z & int(value)).bit_count() % 2) for value in values])
    if not np.array_equal(relative, expected):
        return None
    scalar = int(operator[1][0] - 2 * (x & z).bit_count()) % 8
    if scalar not in (0, 4):
        return None
    return x | (z << width), scalar // 4


def _conjugate(operator, probe):
    return _compose(_compose(operator, probe), _adjoint(operator))


def _is_clifford(operator, width):
    return all(_pauli_label(_conjugate(operator, _pauli(1 << j, width)), width)
               is not None for j in range(2 * width))


def _domain(operator, width):
    """Prune translations before checking the possible phase labels."""
    values = np.arange(1 << width)
    inverse = np.argsort(operator[0])
    result = {}
    for x in range(1 << width):
        image = operator[0][inverse ^ x]
        if not np.array_equal(image, values ^ image[0]):
            continue
        for z in range(1 << width):
            label = x | (z << width)
            output = _pauli_label(_conjugate(operator, _pauli(label, width)), width)
            if output is not None:
                result[label] = output
    return result


def _product_phase(first, second, width):
    mask = (1 << width) - 1
    x, y = first & mask, first >> width
    s, t = second & mask, second >> width
    return (-(x & t).bit_count() + (s & y).bit_count()
            + 2 * ((x & y & t).bit_count() + (s & y & t).bit_count()
                   + (x & s & y).bit_count() + (x & s & t).bit_count())) % 4


def _native_schedule():
    schedule = []
    for wires in ((2, 3, 4), (1, 3, 5), (1, 2, 6)):
        schedule.extend(_ccz(*wires))
    for first, second, target in ((0, 1, 4), (0, 2, 5), (0, 3, 6)):
        schedule += ([('C', [('H', target)])] + _ccz(first, second, target)
                     + [('C', [('H', target)])])
    return schedule


def _flatten(schedule):
    return [('T', gate[1], 1 if gate[0] == 'T' else -1)
            if gate[0] in ('T', 'TDG') else gate
            for _, gates in schedule for gate in gates]


class ControlledPauliTransport(unittest.TestCase):
    def test_hermitian_multiplication_and_quartic_scalar(self):
        for width in (1, 2):
            for first, second in itertools.product(range(1 << (2 * width)), repeat=2):
                phase = _product_phase(first, second, width)
                np.testing.assert_allclose(_tau(first, width) @ _tau(second, width),
                                           1j ** phase * _tau(first ^ second, width),
                                           rtol=0, atol=1e-15)
        phases = [_product_phase(value & 3, value >> 2, 1) for value in range(16)]
        self.assertEqual(_degree([phase & 1 for phase in phases]), 2)
        self.assertEqual(_degree([phase >> 1 for phase in phases]), 4)
        self.assertTrue(any(phase % 2 for phase in phases))

        hadamards = _wire(H, 0, 3) @ _wire(H, 1, 3) @ _wire(H, 2, 3)
        raw_phase, hermitian_middle = [], []
        omega = np.exp(1j * np.pi / 4)
        for control in range(16):
            f, g = _fg(control)
            dot = (f & g).bit_count()
            original = 1j ** (-dot) * _tau(f | (g << 3), 3)
            h = -2 * dot % 8
            np.testing.assert_allclose(original, omega ** h * _tau(f | (g << 3), 3),
                                       rtol=0, atol=2e-15)
            # Fixed target Hadamards swap f,g but retain a quartic raw sign.
            swapped = 1j ** (-dot) * _tau(g | (f << 3), 3)
            np.testing.assert_allclose(hadamards @ original @ hadamards,
                                       (-1) ** dot * swapped, rtol=0, atol=2e-15)
            raw_phase.append(dot % 2)
            hermitian_middle.append((h >> 1) & 1)
        self.assertEqual(_degree(raw_phase), 4)
        self.assertEqual(_degree(hermitian_middle), 4)
        self.assertEqual(raw_phase, [int(control == 15) for control in range(16)])

    def test_gm_exact_domain_c3_inverse_and_literal_coordinates(self):
        source, width = _gm(), 7
        for j in range(14):
            image = _conjugate(source, _pauli(1 << j, width))
            self.assertTrue(_is_clifford(image, width), msg=f"generator {j}")

        domain = _domain(source, width)
        self.assertEqual(domain, {z << width: (z << width, 0) for z in range(16)})
        self.assertLess(4, width)  # No seven-dimensional Pauli Lagrangian.

        square = _compose(source, source)
        np.testing.assert_array_equal(square[0], np.arange(128))
        np.testing.assert_array_equal(square[1], [4 * int(value & 15 == 15)
                                                 for value in range(128)])
        reverse_image = _conjugate(_adjoint(source), _pauli(1, width))
        expected_permutation, expected_phase = [], []
        for value in range(128):
            control = value & 15
            f, _ = _fg(control)
            shifted_f, _ = _fg(control ^ 1)
            expected_permutation.append(value ^ 1 ^ ((f ^ shifted_f) << 4))
            expected_phase.append(4 * int(control & 14 == 14))
        np.testing.assert_array_equal(reverse_image[0], expected_permutation)
        np.testing.assert_array_equal(reverse_image[1], expected_phase)
        self.assertFalse(_is_clifford(reverse_image, width))

        # Beigi--Shor Section 1.1: three R-controlled swaps after four CCZs.
        # Original order: A1,A2,A3,B1,B2,B3,R. The new target choices are
        # A1,B2,B3, which remove any extra cubic control phase.
        coordinate_permutation, original_permutation, original_phase = [], [], []
        for value in range(128):
            a1, a2, a3, b1, b2, b3, control = [(value >> j) & 1 for j in range(7)]
            new = (control | ((a1 ^ b1) << 1) | ((a2 ^ b2) << 2)
                   | ((a3 ^ b3) << 3) | (a1 << 4) | (b2 << 5) | (b3 << 6))
            coordinate_permutation.append(new)
            output = value
            if control:
                for j in range(3):
                    if ((value >> j) ^ (value >> (j + 3))) & 1:
                        output ^= (1 << j) | (1 << (j + 3))
            original_permutation.append(output)
            phase = a1 * a2 * a3 ^ a1 * b2 * b3 ^ b1 * a2 * b3 ^ b1 * b2 * a3
            original_phase.append(4 * phase)
        coordinates = np.array(coordinate_permutation), np.zeros(128, dtype=int)
        original = np.array(original_permutation), np.array(original_phase)
        transformed = _conjugate(coordinates, original)
        np.testing.assert_array_equal(transformed[0], source[0])
        np.testing.assert_array_equal(transformed[1], source[1])

    def test_actual_gm_endpoint_four_block_signed_recovery(self):
        source, width = _gm(), 7
        target_paulis = np.array([_tau(label, 3) for label in range(64)])
        source_blocks = []
        for control in range(16):
            f, g = _fg(control)
            source_blocks.append(1j ** (-(f & g).bit_count())
                                 * target_paulis[f | (g << 3)])
        roots = np.array([1, 1j, -1, -1j])
        saw_imaginary_alpha = saw_essential_carry = False
        largest_integer_dot = checked = 0

        for translation in (1, 15):
            # Actual Pauli inner words give exact Clifford endpoints by C3.
            inner = _pauli(translation, width)
            endpoint = _conjugate(source, inner)
            self.assertTrue(_is_clifford(endpoint, width))
            full_j = _compose(_adjoint(inner), endpoint)
            blocks = [source_blocks[z ^ translation] @ source_blocks[z].conj().T
                      for z in range(16)]
            for z, block in enumerate(blocks):
                expected = np.zeros((8, 8), dtype=complex)
                for target in range(8):
                    column = z | (target << 4)
                    self.assertEqual(int(full_j[0][column]) & 15, z)
                    expected[int(full_j[0][column]) >> 4, target] = np.exp(
                        1j * np.pi / 4 * full_j[1][column])
                np.testing.assert_allclose(block, expected, rtol=0, atol=2e-15)

            # All 64 target probes and two fixed control probes: 128 per
            # endpoint, rather than all 4^7 global Paulis or dense operators.
            for core, (u, v) in itertools.product(range(64), ((14, 0), (15, 5))):
                def transition(z):
                    return blocks[z ^ u] @ target_paulis[core] @ blocks[z].conj().T

                at_zero = transition(0)
                coefficients = np.einsum('pij,ji->p', target_paulis, at_zero) / 8
                target_label = int(np.argmax(np.abs(coefficients)))
                mu = int(np.argmin(np.abs(roots - coefficients[target_label])))
                np.testing.assert_allclose(at_zero, roots[mu] * target_paulis[target_label],
                                           rtol=0, atol=2e-15)
                saw_imaginary_alpha |= bool(mu % 2)

                output_v = v
                for j in range(4):
                    # Each output bit uses only blocks 0,u,e_j,e_j+u.
                    ratio = transition(1 << j) @ at_zero.conj().T
                    bit = int(np.trace(ratio).real < 0)
                    np.testing.assert_allclose(ratio, (-1) ** bit * np.eye(8),
                                               rtol=0, atol=2e-15)
                    output_v ^= bit << j

                output_dot, input_dot = (u & output_v).bit_count(), (u & v).bit_count()
                even_phase = (mu - output_dot + input_dot) % 4
                self.assertEqual(even_phase % 2, 0)
                sign = even_phase // 2
                parity_only = (mu - (output_dot % 2) + (input_dot % 2)) % 4 // 2
                saw_essential_carry |= parity_only != sign
                largest_integer_dot = max(largest_integer_dot, output_dot)
                probe = u | ((core & 7) << 4) | (v << 7) | ((core >> 3) << 11)
                output = (u | ((target_label & 7) << 4) | (output_v << 7)
                          | ((target_label >> 3) << 11))
                direct = _pauli_label(_conjugate(full_j, _pauli(probe, width)), width)
                self.assertEqual(direct, (output, sign))
                checked += 1

        self.assertEqual(checked, 256)
        self.assertTrue(saw_imaginary_alpha)
        self.assertGreaterEqual(largest_integer_dot, 2)
        self.assertTrue(saw_essential_carry)

    def test_literal_native_forward_inverse_and_paid_ledger(self):
        schedule = _native_schedule()
        word = _flatten(schedule)
        counts = Counter(gate[0] for gate in word)
        self.assertEqual(counts, {'T': 42, 'CX': 60, 'H': 6})
        self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 30)
        for kind, gates in schedule:
            if kind == 'T':
                targets = [gate[1] for gate in gates]
                self.assertEqual(len(targets), len(set(targets)))
        wires = {wire for gate in word for wire in (gate[1:2] if gate[0] == 'T' else gate[1:])}
        self.assertEqual(wires, set(range(7)))  # Seven data wires; no helpers.
        self.assertTrue(any(gate[0] == 'T' and gate[1] < 4 for gate in word))

        permutation, phase = _gm()
        expected = np.zeros((128, 128), dtype=complex)
        expected[permutation, np.arange(128)] = np.exp(1j * np.pi / 4 * phase)
        actual = _native_matrix(word, 7)
        inverse_word = _native_inverse(word)
        self.assertEqual(Counter(gate[0] for gate in inverse_word), counts)
        actual_inverse = _native_matrix(inverse_word, 7)
        np.testing.assert_allclose(actual, expected, rtol=0, atol=6e-14)
        np.testing.assert_allclose(actual_inverse, expected.conj().T, rtol=0, atol=6e-14)
        np.testing.assert_allclose(actual_inverse @ actual, np.eye(128), rtol=0, atol=8e-14)

    def test_two_copy_domain_and_allocation_without_dense_tensor(self):
        controls, targets, width = 8, 6, 14
        coefficients = {}
        for copy in range(2):
            c, t = 4 * copy, 3 * copy
            for left, right, target, is_z in (
                    (0, 1, 0, False), (0, 2, 1, False), (0, 3, 2, False),
                    (2, 3, 0, True), (1, 3, 1, True), (1, 2, 2, True)):
                coefficients[(c + left, c + right)] = 1 << (t + target + targets * is_z)

        # Translation derivatives are constant iff every linear coefficient
        # vanishes. This checks the tensor source's control radical directly.
        control_radical = []
        for translation in range(1 << controls):
            derivative = [0] * controls
            for (left, right), label in coefficients.items():
                if (translation >> right) & 1:
                    derivative[left] ^= label
                if (translation >> left) & 1:
                    derivative[right] ^= label
            if not any(derivative):
                control_radical.append(translation)
        self.assertEqual(control_radical, [0])

        def pair(first, second):
            mask = (1 << targets) - 1
            return (((first & mask) & (second >> targets)).bit_count()
                    + ((second & mask) & (first >> targets)).bit_count()) % 2

        target_kernel = [label for label in range(1 << (2 * targets))
                         if all(pair(label, coefficient) == 0
                                for coefficient in coefficients.values())]
        self.assertEqual(target_kernel, [0])
        # Exactly the eight control Zs remain: E_(F tensor F)=E_F direct-sum E_F.
        self.assertEqual(1 << controls, len(_domain(_gm(), 7)) ** 2)
        self.assertLess(controls, width)

        first_word = _flatten(_native_schedule())
        second_word = [(gate[0], gate[1] + 7, gate[2]) if gate[0] == 'T'
                       else (gate[0], *(wire + 7 for wire in gate[1:]))
                       for gate in first_word]
        paid_word = first_word + second_word
        counts = Counter(gate[0] for gate in paid_word)
        self.assertEqual((counts['T'], counts['CX'] + counts['H']), (84, 132))
        used = {wire for gate in paid_word
                for wire in (gate[1:2] if gate[0] == 'T' else gate[1:])}
        self.assertEqual(used, set(range(width)))
        # Serial copies have a paid depth bound 60 and zero extra helpers;
        # width independence of the transport bound does not erase this cost.
        self.assertEqual(2 * sum(kind == 'T' for kind, _ in _native_schedule()), 60)


if __name__ == '__main__':
    unittest.main()
