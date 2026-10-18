"""A returned affine query can have an odd-order formal outer baseline.

Every small operator column and actual inverse is retained. The target is
affine; this is not a counterexample to any Boolean degree theorem.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_baseline_word_obstruction import (
    _coordinate_word, _gate_matrix, _inverse, _word_matrix,
)
from tests.test_outer_axis_normalizer import _tableau


def _returned_affine_word(work, address, form):
    fanout = [('CX', address, address + j) for j in range(1, work)]
    coordinate = _coordinate_word(work, address)
    forward = [('T', address + j) for j in range(work)]
    backward = [('Tdg', address + j) for j in range(work)]
    pauli = [('CX', bit, address + j)
             for bit in range(address) if form >> bit & 1
             for j in range(work)]
    output = [('CX', bit, address) for bit in range(address)
              if form >> bit & 1]
    return (_inverse(fanout) + forward + _inverse(coordinate) + backward
            + pauli + pauli + forward + coordinate + backward + fanout + output)


class ReturnedAffineAmbient(unittest.TestCase):
    def test_full_literal_return_actual_inverse_and_ledger(self):
        for work in range(1, 4):
            address = 2
            for form in range(1, 1 << address):
                word = _returned_affine_word(work, address, form)
                actual = _word_matrix(word, work + address)
                simplified = _word_matrix(
                    [('CX', bit, address) for bit in range(address)
                     if form >> bit & 1], work + address)
                np.testing.assert_allclose(actual, simplified, atol=5e-14)
                native_inverse = _word_matrix(_inverse(word), work + address)
                np.testing.assert_allclose(native_inverse @ actual,
                                           np.eye(len(actual)), atol=7e-14)
                self.assertEqual(sum(gate[0] in ('T', 'Tdg') for gate in word),
                                 4 * work)
                self.assertEqual(sum(gate[0] in ('H', 'CX') for gate in word),
                                 6 * work - 4 + (2 * work + 1) * form.bit_count())
                self.assertTrue(all(gate[1] >= address for gate in word
                                    if gate[0] in ('T', 'Tdg')))

    def test_actual_comparison_and_odd_formal_baseline(self):
        for work in range(1, 4):
            identity = np.eye(1 << work, dtype=complex)
            fanout = _word_matrix([('CX', 0, j) for j in range(1, work)], work)
            coordinate = _word_matrix(_coordinate_word(work), work)
            diagonal = _word_matrix([('T', j) for j in range(work)], work)
            x_all = _word_matrix([('X', j) for j in range(work)], work)
            output = _gate_matrix(('X', 0), work)
            outer = fanout @ diagonal.conj().T
            q = diagonal @ x_all @ diagonal.conj().T
            a = coordinate @ q @ coordinate.conj().T
            np.testing.assert_allclose(outer.conj().T @ output @ outer, q, atol=3e-14)
            for bit in (0, 1):
                aa = a if bit else identity
                qq = q if bit else identity
                p5 = output if bit else identity
                target = p5
                actual = aa.conj().T @ outer.conj().T @ p5.conj().T @ target @ outer
                formal = aa.conj().T @ outer.conj().T @ p5.conj().T @ outer
                np.testing.assert_allclose(actual, aa, atol=4e-14)
                np.testing.assert_allclose(formal, aa @ qq, atol=4e-14)
                np.testing.assert_allclose(formal @ qq, actual, atol=5e-14)
            symplectic, _ = _tableau(a @ q, work)
            power = np.eye(2 * work, dtype=int)
            for exponent in range(1, 2 * work + 2):
                power = power @ symplectic % 2
                self.assertEqual(np.array_equal(power, np.eye(2 * work, dtype=int)),
                                 exponent == 2 * work + 1)


if __name__ == '__main__':
    unittest.main()
