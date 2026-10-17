"""Bounded exact checks of the Gaussian scalar-loop proof.

Integer Gauss sums retain literal signs. The path fixture is a real closed
Gaussian circuit with cubic address phase and complete work return. No test
claims an unrestricted native T-depth bound.
"""
from __future__ import annotations

import itertools
import random
import unittest

import numpy as np


def _inverse(rows):
    n = len(rows)
    a = [rows[i] | (1 << (n + i)) for i in range(n)]
    for i in range(n):
        pivot = next((j for j in range(i, n) if a[j] >> i & 1), None)
        if pivot is None:
            return None
        a[i], a[pivot] = a[pivot], a[i]
        for j in range(n):
            if j != i and (a[j] >> i & 1):
                a[j] ^= a[i]
    return [row >> n for row in a]


def _apply(rows, vector):
    return sum(((row & vector).bit_count() & 1) << i
               for i, row in enumerate(rows))


def _multiply(a, b):
    out = []
    for row in a:
        value = 0
        for i, second in enumerate(b):
            if row >> i & 1:
                value ^= second
        out.append(value)
    return out


def _polar(q):
    return [sum((((q[i] >> j) ^ (q[j] >> i)) & 1) << j
                for j in range(len(q))) for i in range(len(q))]


def _qvalue(q, vector):
    return sum((vector >> i & 1) * ((q[i] & vector).bit_count() & 1)
               for i in range(len(q))) & 1


def _c2(a):
    return sum(((a[i] >> i & 1) * (a[j] >> j & 1)
                + (a[i] >> j & 1) * (a[j] >> i & 1))
               for i in range(len(a)) for j in range(i + 1, len(a))) & 1


def _upper_form(n, mask):
    rows = [0] * n
    for bit, (i, j) in enumerate((i, j) for i in range(n)
                                for j in range(i, n)):
        rows[i] |= ((mask >> bit) & 1) << j
    return rows


def _gauss(q):
    return sum(1 - 2 * _qvalue(q, v) for v in range(1 << len(q)))


def _degree(values):
    a = list(values)
    for i in range((len(a) - 1).bit_length()):
        for mask in range(len(a)):
            if mask >> i & 1:
                a[mask] ^= a[mask ^ (1 << i)]
    return max((mask.bit_count() for mask, value in enumerate(a) if value),
               default=-1)


def _cubic_loop(x):
    """(F_B D_q)^3 with A=[[1,0],[x0,1]], Arf=x0*x1*x2.

    Gaussian D entries are lists (coefficient, first wire, second wire);
    repeated wire indices mean a linear phase. The in-place B map is
    explicitly factored into swaps and one parallel Gaussian phase round.
    Adjacent H gates in C_B followed by H_all are already canceled.
    """
    q = [('D', [(1, 0, 2), (1, 1, 3), (x[0], 1, 2),
                (x[1], 0, 0), (x[2], 3, 3)])]
    cb_then_fourier = [('SWAP', 0, 2), ('SWAP', 1, 3),
                      ('H', 1), ('H', 2),
                      ('D', [(x[0], 0, 1), (x[0], 3, 2)]),
                      ('H', 0), ('H', 3)]
    return (q + cb_then_fourier) * 3


def _path_phase(word, path):
    wires, phase, used = 0, 0, 0
    for gate in word:
        if gate[0] == 'H':
            wire = gate[1]
            value = path >> used & 1
            phase ^= (wires >> wire & 1) & value
            wires = (wires & ~(1 << wire)) | (value << wire)
            used += 1
        elif gate[0] == 'SWAP':
            i, j = gate[1:]
            if ((wires >> i) ^ (wires >> j)) & 1:
                wires ^= (1 << i) | (1 << j)
        elif gate[0] == 'D':
            for coefficient, i, j in gate[1]:
                phase ^= coefficient & (wires >> i & 1) & (wires >> j & 1)
        else:
            raise ValueError(gate)
    return phase


def _path_form(word):
    h = sum(gate[0] == 'H' for gate in word)
    constant = _path_phase(word, 0)
    q = [0] * h
    for i in range(h):
        q[i] |= (_path_phase(word, 1 << i) ^ constant) << i
        for j in range(i + 1, h):
            coefficient = (_path_phase(word, (1 << i) | (1 << j))
                           ^ _path_phase(word, 1 << i)
                           ^ _path_phase(word, 1 << j) ^ constant)
            q[i] |= coefficient << j
    return constant, q


def _suffix_labels(word):
    labels = []
    for position, gate in enumerate(word):
        if gate[0] != 'H':
            continue
        p, z = 0, 1 << gate[1]
        for later in word[position + 1:]:
            if later[0] == 'H':
                bit = 1 << later[1]
                if bool(p & bit) != bool(z & bit):
                    p ^= bit
                    z ^= bit
            elif later[0] == 'SWAP':
                i, j = later[1:]
                for which in (0, 1):
                    v = p if which == 0 else z
                    if ((v >> i) ^ (v >> j)) & 1:
                        v ^= (1 << i) | (1 << j)
                    if which == 0:
                        p = v
                    else:
                        z = v
            else:
                for coefficient, i, j in later[1]:
                    if coefficient and i != j:
                        z ^= (p >> j & 1) << i
                        z ^= (p >> i & 1) << j
        labels.append((p, z))
    return labels


def _native_matrix(word, width):
    matrix = np.eye(1 << width, dtype=complex)
    indices = np.arange(1 << width)
    for gate in word:
        if gate[0] == 'H':
            low = indices[(indices >> gate[1]) & 1 == 0]
            high = low | (1 << gate[1])
            old_low, old_high = matrix[low].copy(), matrix[high].copy()
            matrix[low] = (old_low + old_high) / np.sqrt(2)
            matrix[high] = (old_low - old_high) / np.sqrt(2)
        elif gate[0] == 'SWAP':
            i, j = gate[1:]
            change = ((indices >> i) ^ (indices >> j)) & 1
            matrix = matrix[indices ^ (change << i) ^ (change << j)]
        else:
            phase = np.zeros(1 << width, dtype=np.int64)
            for coefficient, i, j in gate[1]:
                phase ^= coefficient * ((indices >> i) & 1) * ((indices >> j) & 1)
            matrix *= (1 - 2 * phase)[:, None]
    return matrix


class GaussianScalarLoopTests(unittest.TestCase):
    def test_arf_formula_exhaustive_small_and_sampled_larger_forms(self):
        rng = random.Random(18423)
        for n, masks in ((2, range(1 << 3)),
                         (4, range(1 << 10)),
                         (6, [rng.randrange(1 << 21) for _ in range(300)]),
                         (8, [rng.randrange(1 << 36) for _ in range(150)])):
            for mask in masks:
                q = _upper_form(n, mask)
                inverse = _inverse(_polar(q))
                if inverse is None:
                    continue
                sign = _c2(_multiply(inverse, q)) ^ ((n // 2) * (n // 2 - 1) // 2 & 1)
                self.assertEqual(_gauss(q), (1 - 2 * sign) * (1 << (n // 2)))
                # Changing the matrix representative by an alternating matrix
                # leaves q and the formula unchanged.
                changed = q.copy()
                changed[0] ^= 1 << 1
                changed[1] ^= 1
                self.assertEqual(_gauss(changed), _gauss(q))
                self.assertEqual(_c2(_multiply(inverse, changed)),
                                 _c2(_multiply(inverse, q)))

    def test_literal_generalized_fourier_cube_for_all_four_bit_forms(self):
        n = 4
        for mask in range(1 << 10):
            q = _upper_form(n, mask)
            b = _polar(q)
            if _inverse(b) is None:
                continue
            # Integer matrix is 2^(n/2) F_B D_q; no floating signs.
            matrix = np.array([[1 - 2 * (((u & _apply(b, v)).bit_count()
                                          ^ _qvalue(q, v)) & 1)
                                for v in range(1 << n)]
                               for u in range(1 << n)], dtype=np.int64)
            sign = 1 if _gauss(q) > 0 else -1
            expected = sign * (1 << (3 * n // 2)) * np.eye(1 << n, dtype=np.int64)
            np.testing.assert_array_equal(matrix @ matrix @ matrix, expected)

    def test_cubic_scalar_loop_paths_markers_and_inverse_polar(self):
        scalar_table = []
        inverse_entries = []
        for x in itertools.product((0, 1), repeat=3):
            word = _cubic_loop(x)
            constant, q = _path_form(word)
            h = len(q)
            self.assertEqual(h, 12)
            inverse = _inverse(_polar(q))
            self.assertIsNotNone(inverse)
            scalar = x[0] & x[1] & x[2]
            scalar_table.append(scalar)
            inverse_entries.append(inverse)
            np.testing.assert_allclose(_native_matrix(word, 4),
                                       (1 - 2 * scalar) * np.eye(16), atol=2e-15)
            # Every path agrees with the derived quadratic, not just basis probes.
            phases = [_path_phase(word, path) for path in range(1 << h)]
            self.assertEqual(phases, [constant ^ _qvalue(q, path)
                                      for path in range(1 << h)])
            self.assertEqual(sum(1 - 2 * v for v in phases),
                             (1 - 2 * scalar) * (1 << (h // 2)))
            labels = _suffix_labels(word)
            for i in range(h):
                self.assertEqual(inverse[i] >> i & 1, 0)
                for j in range(i + 1, h):
                    # Later markers are on the left in the output Pauli product.
                    self.assertEqual(inverse[i] >> j & 1,
                                     (labels[j][1] & labels[i][0]).bit_count() & 1)
            markers = [0] + [1 << i for i in range(h)]
            markers += [(1 << i) | (1 << j) for i in range(h) for j in range(i + 1, h)]
            baseline = sum(1 - 2 * v for v in phases)
            for marker in markers:
                shifted = sum(1 - 2 * (value ^ ((path & marker).bit_count() & 1))
                              for path, value in enumerate(phases))
                expected = baseline * (1 - 2 * _qvalue(q, _apply(inverse, marker)))
                self.assertEqual(shifted, expected)
        self.assertEqual(_degree(scalar_table), 3)
        for i in range(12):
            for j in range(12):
                self.assertLessEqual(_degree([row[i] >> j & 1 for row in inverse_entries]), 2)

    def test_realification_of_s_is_a_real_clifford(self):
        # Basis order is (rebit, data); S data controls J=XZ on the rebit.
        s = np.diag([1, 1j])
        realified = np.block([[s.real, -s.imag], [s.imag, s.real]])
        actual = np.zeros((4, 4))
        for rebit in (0, 1):
            for data in (0, 1):
                column = 2 * rebit + data
                row = 2 * (rebit ^ data) + data
                actual[row, column] = (-1) ** (rebit * data)
        np.testing.assert_array_equal(actual, realified)


if __name__ == '__main__':
    unittest.main()
