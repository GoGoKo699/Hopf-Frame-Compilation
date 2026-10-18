"""Bounded evidence for the baseline-only word obstruction.

The universal quantifier over words and widths is proved analytically.
These checks verify its native source, all small columns, actual inverse,
symplectic order, incidence argument and central-factor conclusion.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_read_only_three_layer_rigidity import H, T


def _gate_matrix(gate, width):
    name, *wires = gate
    result = np.zeros((1 << width, 1 << width), dtype=complex)
    for column in range(1 << width):
        if name == "CX":
            control, target = wires
            row = column ^ ((1 << target) if column >> control & 1 else 0)
            result[row, column] = 1
        else:
            target = wires[0]
            matrix = {"H": H, "T": T, "Tdg": T.conj().T,
                      "X": np.array([[0, 1], [1, 0]]),
                      "Z": np.diag([1, -1])}[name]
            for bit in range(2):
                row = (column & ~(1 << target)) | (bit << target)
                result[row, column] = matrix[bit, column >> target & 1]
    return result


def _word_matrix(word, width):
    result = np.eye(1 << width, dtype=complex)
    for gate in word:
        result = _gate_matrix(gate, width) @ result
    return result


def _inverse(word):
    names = {"H": "H", "CX": "CX", "T": "Tdg", "Tdg": "T"}
    return [(names[gate[0]], *gate[1:]) for gate in reversed(word)]


def _coordinate_word(work, offset=0):
    return ([('H', offset + j) for j in range(work)]
            + [('CX', offset + j, offset + j + 1)
               for j in reversed(range(work - 1))])


def _source_word(work, address, form):
    coordinate = _coordinate_word(work, address)
    forward = [('T', address + j) for j in range(work)]
    backward = [('Tdg', address + j) for j in range(work)]
    controlled = [('CX', bit, address + j)
                  for bit in range(address) if form >> bit & 1
                  for j in range(work)]
    # Right-to-left A_x N_x: execute N_x, then A_x.
    return (_inverse(coordinate) + backward + controlled + forward
            + coordinate + backward + controlled + forward)


def _symplectic_source(work):
    width = 2 * work
    identity = np.eye(width, dtype=np.int64)
    pairing = np.block([[np.zeros((work, work), dtype=np.int64),
                         np.eye(work, dtype=np.int64)],
                        [np.eye(work, dtype=np.int64),
                         np.zeros((work, work), dtype=np.int64)]])
    odd, even = identity.copy(), identity.copy()
    axes = []
    for j in range(work):
        z = np.zeros(width, dtype=np.int64)
        z[work + j] = 1
        x = np.zeros(width, dtype=np.int64)
        x[j] = 1
        if j + 1 < work:
            x[j + 1] = 1
        axes.extend((z, x))
        odd = (identity + np.outer(z, pairing @ z)) @ odd % 2
        even = (identity + np.outer(x, pairing @ x)) @ even % 2
    return odd @ even % 2, np.array(axes), pairing


class BaselineWordObstructionTests(unittest.TestCase):
    def test_native_source_axis_images_full_columns_and_actual_inverse(self):
        for work in range(1, 4):
            coordinate = _word_matrix(_coordinate_word(work), work)
            diagonal = _word_matrix([('T', j) for j in range(work)], work)
            x_all = _word_matrix([('X', j) for j in range(work)], work)
            a = diagonal @ x_all @ diagonal.conj().T
            n = coordinate @ a @ coordinate.conj().T
            np.testing.assert_allclose(a @ a, np.eye(1 << work), atol=2e-14)
            np.testing.assert_allclose(n @ n, np.eye(1 << work), atol=2e-14)
            for j in range(work):
                expected_word = [('X', j)]
                if j + 1 < work:
                    expected_word.append(('X', j + 1))
                expected = _word_matrix(expected_word, work)
                actual = coordinate @ _gate_matrix(('Z', j), work) @ coordinate.conj().T
                np.testing.assert_allclose(actual, expected, atol=2e-14)
            address = 2
            for form in range(1, 1 << address):
                word = _source_word(work, address, form)
                actual = _word_matrix(word, address + work)
                expected = np.zeros_like(actual)
                for x in range(1 << address):
                    fibre = a @ n if (x & form).bit_count() % 2 else np.eye(1 << work)
                    for column in range(1 << work):
                        for row in range(1 << work):
                            expected[x | (row << address), x | (column << address)] = fibre[row, column]
                np.testing.assert_allclose(actual, expected, atol=4e-14)
                inverse = _word_matrix(_inverse(word), address + work)
                np.testing.assert_allclose(inverse, actual.conj().T, atol=4e-14)
                np.testing.assert_allclose(inverse @ actual, np.eye(len(actual)), atol=5e-14)
                self.assertEqual(sum(g[0] in ('T', 'Tdg') for g in word), 4 * work)
                self.assertEqual(sum(g[0] in ('H', 'CX') for g in word),
                                 4 * work - 2 + 2 * work * form.bit_count())
                # Four distinct physical bands, each with disjoint targets.
                bands = []
                for gate in word:
                    if gate[0] in ('T', 'Tdg'):
                        if not bands or len(bands[-1]) == work:
                            bands.append([])
                        bands[-1].append(gate[1])
                self.assertEqual(len(bands), 4)
                self.assertTrue(all(len(set(band)) == work for band in bands))

    def test_chain_pairings_and_exact_odd_symplectic_orders(self):
        for work in range(1, 17):
            source, axes, pairing = _symplectic_source(work)
            expected = np.zeros((2 * work, 2 * work), dtype=np.int64)
            for j in range(2 * work - 1):
                expected[j, j + 1] = expected[j + 1, j] = 1
            np.testing.assert_array_equal(axes @ pairing @ axes.T % 2, expected)
            power = np.eye(2 * work, dtype=np.int64)
            for exponent in range(1, 2 * work + 2):
                power = power @ source % 2
                self.assertEqual(np.array_equal(power, np.eye(2 * work)),
                                 exponent == 2 * work + 1)

    def test_incidence_gram_is_nonsingular(self):
        for address in range(1, 7):
            labels = range(1, 1 << address)
            matrix = np.array([[(linear & z).bit_count() % 2 for z in labels]
                               for linear in labels], dtype=np.int64)
            if address == 1:
                np.testing.assert_array_equal(matrix, [[1]])
            else:
                expected = (1 << (address - 2)) * (
                    np.eye(len(matrix), dtype=np.int64)
                    + np.ones_like(matrix))
                np.testing.assert_array_equal(matrix @ matrix.T, expected)
            self.assertEqual(np.linalg.matrix_rank(matrix), len(matrix))

    def test_cube_words_fail_on_literal_native_baselines(self):
        # Any ordering has this same result because this baseline is cyclic.
        # Signed and all-positive cube products both fail; a target is absent.
        source, _, _ = _symplectic_source(1)
        for dimension in range(1, 9):
            count = 1 << dimension
            parity = [(x.bit_count() % 2) for x in range(count)]
            for signs in ([1] * count,
                          [(-1) ** x.bit_count() for x in range(count)]):
                exponent = sum(sign * value for sign, value in zip(signs, parity))
                self.assertNotEqual(exponent % 3, 0)
                self.assertFalse(np.array_equal(
                    np.linalg.matrix_power(source, exponent % 3) % 2,
                    np.eye(2, dtype=np.int64)))

        # The actual one-qubit Clifford witness also stays outside the
        # robust scalar radius. Minimize over a scalar by its eigenvalue arc.
        x = np.array([[0, 1], [1, 0]], dtype=complex)
        a = T @ x @ T.conj().T
        n = H @ a @ H
        for exponent in (1, 2):
            unitary = np.linalg.matrix_power(a @ n, exponent)
            angles = np.sort(np.angle(np.linalg.eigvals(unitary)))
            gaps = np.diff(np.append(angles, angles[0] + 2 * np.pi))
            projective_distance = 2 * np.sin((2 * np.pi - max(gaps)) / 4)
            self.assertGreaterEqual(projective_distance + 2e-14, 1 / np.sqrt(2))

    def test_balanced_word_cannot_detect_disjoint_target(self):
        # This tests the exact final cancellation, including inverse order.
        # It does not assert that a balanced word cancels the baseline.
        baselines = [np.eye(2), H, T @ H @ T.conj().T,
                     H @ T @ H @ T.conj().T]
        word = [(1, 1), (2, 1), (1, -1), (3, 1), (2, -1), (3, -1)]
        target = np.diag([1, -1])
        for bits in range(8):
            baseline_word = np.eye(4, dtype=complex)
            target_word = np.eye(4, dtype=complex)
            for label, sign in word:
                baseline = np.kron(baselines[label], np.eye(2))
                value = np.kron(baselines[label],
                                target if bits >> (label - 1) & 1 else np.eye(2))
                if sign == -1:
                    baseline = baseline.conj().T
                    value = value.conj().T
                baseline_word = baseline_word @ baseline
                target_word = target_word @ value
            np.testing.assert_allclose(target_word, baseline_word, atol=2e-14)


if __name__ == '__main__':
    unittest.main()
