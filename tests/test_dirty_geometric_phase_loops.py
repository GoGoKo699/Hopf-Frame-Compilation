"""Bounded checks for dirty geometric scalar loops and their paid word.

The exact Clifford-angle digital fixture is deliberately coarse and is
not an accurate indicator. It verifies every physical column, the
Hamiltonian decomposition and the actual inverse. No finite check is
used as an adiabatic or asymptotic compiler theorem.
"""
from functools import lru_cache
from itertools import product
import unittest

import numpy as np

from tests.exact_matrices import CO, CZ, adjoint, eye, native_columns, scale
from tests.test_majorana_walsh_query import _inverse


def _omega_word():
    return [('S', 0), ('H', 0)] * 3


def _rotation_word(k, axis='Z'):
    # exp(-i k*pi*Z/4), including its literal SU(2) scalar.
    z = _omega_word() * ((-k) % 8) + [('S', 0)] * (k % 4)
    return z if axis == 'Z' else [('H', 0)] + z + [('H', 0)]


def _y_word():
    return [('S', 0)] * 3 + [('H', 0), ('S', 0), ('S', 0), ('H', 0), ('S', 0)]


@lru_cache(None)
def _local_table(k, axis='Z', pair=False, inverse=False):
    word = _rotation_word(k, axis)
    if pair:
        word = [('CX', 0, 1)] + [(name, 1) for name, _ in word] + [('CX', 0, 1)]
    if inverse:
        word = _inverse(word)
    matrix = native_columns(2 if pair else 1, word)
    omega = native_columns(1, _omega_word())[0][0]
    roots = [CO]
    for _ in range(7):
        roots.append(roots[-1] * omega)
    images, phases = [], []
    for col in range(len(matrix)):
        support = [(row, matrix[row][col]) for row in range(len(matrix))
                   if matrix[row][col] != CZ]
        assert len(support) == 1
        row, entry = support[0]
        images.append(row)
        phases.append(roots.index(entry))
    return tuple(images), tuple(phases), tuple(word)


def _apply_table(states, phases, wires, table):
    images, angles, _ = table
    local = sum(((states >> wire) & 1) << j for j, wire in enumerate(wires))
    for j, wire in enumerate(wires):
        changed = ((np.asarray(images)[local] ^ local) >> j) & 1
        states ^= changed << wire
    phases += np.asarray(angles)[local]


class DirtyGeometricPhaseLoopsTests(unittest.TestCase):
    def test_literal_native_half_echo_retains_the_central_minus_sign(self):
        omega = native_columns(1, _omega_word())[0][0]
        self.assertEqual(native_columns(1, _omega_word()), scale(eye(2), omega))
        dynamic = _rotation_word(1)
        for winding in (0, 1):
            middle = (_omega_word() * 6 + _y_word()) if winding else []
            end = _omega_word() * 4 if winding else []
            # Chronological U_first, Y, U_second, Y.
            word = dynamic + middle + _y_word() + _inverse(middle) + dynamic + end + _y_word()
            actual = native_columns(1, word)
            self.assertEqual(actual, scale(eye(2), -CO if winding else CO))
            self.assertEqual(native_columns(1, _inverse(word)), adjoint(actual))
            self.assertTrue(all(gate[0] in ('H', 'S') for gate in word))
        # Omitting the SU(2) scalar changes a nonwinding echo by i.
        bad = [('S', 0)] + _y_word() + [('S', 0)] + _y_word()
        self.assertNotEqual(native_columns(1, bad), eye(2))

    def test_all_sector_digital_decomposition_and_actual_inverse(self):
        p, m, width = 2, 4, 10
        initial = np.arange(1 << width, dtype=np.int64)
        instructions = []
        # t=1, J=2: mu=t*pi/J=pi/2, theta=0,pi/2,pi,3pi/2.
        # At these exact Clifford angles each local matrix is monomial.
        for step, (cosine, sine) in enumerate(((1, 0), (0, 1), (-1, 0), (0, -1))):
            row_instructions = []
            for u in range(m):
                signal = p + m + u
                row_instructions.append(((signal,), _local_table(4 + 2 * cosine)))
                row_instructions.append(((p + u, signal), _local_table(1, pair=True)))
                for j in range(p):
                    k = 1 if (u >> j) & 1 else -1
                    row_instructions.append(((j, signal), _local_table(k, pair=True)))
                row_instructions.append(((signal,), _local_table(2 * sine, axis='X')))
            states, phases = initial.copy(), np.zeros_like(initial)
            for wires, table in row_instructions:
                _apply_table(states, phases, wires, table)
            expected, expected_phases = initial.copy(), np.zeros_like(initial)
            for index in range(len(initial)):
                x = index & 3
                bank = (index >> p) & 15
                for u in range(m):
                    signal = p + m + u
                    bit = (index >> signal) & 1
                    twice_a = 2 * (x ^ u).bit_count() + 3 - 2 * ((bank >> u) & 1)
                    z = _local_table(twice_a + 2 * cosine)
                    xrot = _local_table(2 * sine, axis='X')
                    after_z = z[0][bit]
                    after_x = xrot[0][after_z]
                    expected[index] ^= (bit ^ after_x) << signal
                    expected_phases[index] += z[1][bit] + xrot[1][after_z]
            np.testing.assert_array_equal(states, expected)
            np.testing.assert_array_equal(phases % 8, expected_phases % 8)
            instructions.extend(row_instructions)
            if step in (1, 3):
                # Y = i X Z, implemented as two native axial rotations:
                # exp(-i*pi*X/2) exp(-i*pi*Z/2) = iY.
                # There are m=4 copies, so their common i^m scalar is one.
                for u in range(m):
                    signal = p + m + u
                    instructions.extend((((signal,), _local_table(2)),
                                         ((signal,), _local_table(2, axis='X'))))
        states, phases = initial.copy(), np.zeros_like(initial)
        for wires, table in instructions:
            _apply_table(states, phases, wires, table)
        target_phases = np.array([4 * (((int(i) >> p) >> (int(i) & 3)) & 1)
                                  for i in initial])
        # This deliberately coarse native word is not the target indicator phase.
        self.assertTrue(np.any(states != initial) or np.any(phases % 8 != target_phases))
        checked_words = set()
        for wires, table in reversed(instructions):
            word = list(table[2])
            if table[2] not in checked_words:
                inverse_matrix = native_columns(len(wires), _inverse(word))
                matrix = native_columns(len(wires), word)
                self.assertEqual(inverse_matrix, adjoint(matrix))
                checked_words.add(table[2])
            inverse_images = np.argsort(np.asarray(table[0]))
            inverse_phases = -np.asarray(table[1])[inverse_images]
            _apply_table(states, phases, wires,
                         (tuple(inverse_images), tuple(inverse_phases), ()))
        np.testing.assert_array_equal(states, initial)
        np.testing.assert_array_equal(phases % 8, np.zeros_like(initial))

    def test_tensor_operator_error_uses_sum_and_near_minus_identity_is_excluded(self):
        tau = 0.5
        alpha = 2 * np.arcsin(tau / 2)
        for angles in ((0.10, 0.12, 0.13), (0.20, 0.20, 0.20),
                       (0.49 * np.pi, 0.49 * np.pi), (0.01,) * 8):
            spectral = [sum(sign * angle for sign, angle in zip(signs, angles))
                        for signs in product((-1, 1), repeat=len(angles))]
            error = max(abs(np.exp(1j * phase) - 1) for phase in spectral)
            self.assertEqual(error <= tau, sum(angles) <= alpha)
        self.assertLess(2 * np.sin(0.20 / 2), tau)
        self.assertGreater(2 * np.sin(0.60 / 2), tau)
        # Two local -I factors cancel globally; the local convergence
        # hypothesis delta<pi/2 is indispensable to the stated criterion.
        near_minus = (np.pi - 0.01, np.pi - 0.01)
        error = max(abs(np.exp(1j * sum(s * a for s, a in zip(signs, near_minus))) - 1)
                    for signs in product((-1, 1), repeat=2))
        self.assertLess(error, tau)
        self.assertGreater(sum(near_minus), alpha)


if __name__ == '__main__':
    unittest.main()
