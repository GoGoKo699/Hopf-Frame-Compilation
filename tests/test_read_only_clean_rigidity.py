"""Bounded audits of the initialized three-Pauli-band invariant.

The analytic lower bound is in READ_ONLY_TWO_LAYER_RIGIDITY.md. The
native fixture below has general fixed conjugating bands, not two physical
T layers: it checks the clean-syndrome/slopes argument and literal phases.
It makes no construction or asymptotic claim for two-layer queries.
"""
import unittest

import numpy as np

from compiler_robust_hopf.native_coarse_fixture import commutator_word, toffoli_word
from tests.exact_matrices import CO, CZ, native_columns
from tests.test_operator_source_compiler import _adjoint, _word_matrix


def _flat(word):
    return [(name, *wires) for name, wires in word]


def _cz(control, target):
    return [('H', target), ('CX', control, target), ('H', target)]


def _conjugator(target, dirty, clean):
    # Matrix V = H_t Toffoli(c,t;d), so V Z_d V† = Z_d CNOT(c;t).
    return _flat(toffoli_word(clean, target, dirty)) + [('H', target)]


def _anf(table):
    result = list(table)
    for bit in range((len(table) - 1).bit_length()):
        for mask in range(len(table)):
            if mask & (1 << bit):
                result[mask] ^= result[mask ^ (1 << bit)]
    return result


def _degree(table):
    return max((mask.bit_count() for mask, value in enumerate(_anf(table))
                if value), default=0)


def _binary_rank(rows):
    rows = list(rows)
    rank = 0
    for bit in range(max(rows, default=0).bit_length() - 1, -1, -1):
        pivot = next((index for index in range(rank, len(rows))
                      if rows[index] & (1 << bit)), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for index in range(len(rows)):
            if index != rank and rows[index] & (1 << bit):
                rows[index] ^= rows[rank]
        rank += 1
    return rank


def _polar_rank(table):
    assert _degree(table) <= 2
    variables = (len(table) - 1).bit_length()
    rows = []
    for left in range(variables):
        row = 0
        for right in range(variables):
            value = (table[(1 << left) ^ (1 << right)]
                     ^ table[1 << left] ^ table[1 << right] ^ table[0])
            row |= value << right
        rows.append(row)
    return _binary_rank(rows)


class ReadOnlyCleanRigidityTests(unittest.TestCase):
    def test_literal_native_general_band_fixture_returns_all_dirty_columns(self):
        # Addresses a,b; arbitrary target/dirty; initialized clean.
        a, b, target, dirty, clean, width = 0, 1, 2, 3, 4, 5
        v = _conjugator(target, dirty, clean)
        # U_ab = X_c^a Z_d^b V Z_d^b V† X_c^a.
        word = [('CX', a, clean)] + _adjoint(v) + _cz(b, dirty)
        word += v + [('CX', a, clean)] + _cz(b, dirty)
        columns = range(1 << clean)
        actual = native_columns(width, word, columns)
        expected = [[CZ for _ in columns] for _ in range(1 << width)]
        for column in columns:
            f = ((column >> a) & 1) & ((column >> b) & 1)
            expected[column ^ (f << target)][column] = CO
        self.assertEqual(actual, expected)
        self.assertEqual(sum(gate[0] in {'T', 'TDG'} for gate in word), 14)
        self.assertEqual(sum(gate[0] not in {'T', 'TDG'} for gate in word), 26)
        # The address is a Clifford control throughout. Full initialized
        # columns imply return even with an arbitrary entangled reference.
        for name, *wires in word:
            for address in (a, b):
                if address in wires:
                    self.assertEqual((name, wires[0]), ('CX', address))

    def test_actual_inverse_and_two_error_hybrid_allow_clean_leakage(self):
        target, dirty, clean, width = 0, 1, 2, 3
        identity = np.eye(1 << width, dtype=complex)
        j = identity[:, :1 << clean]
        v = _word_matrix(width, _conjugator(target, dirty, clean))
        x_clean = _word_matrix(width, [('X', clean)])
        z_dirty = _word_matrix(width, [('Z', dirty)])
        x_target = _word_matrix(width, [('X', target)])
        # This fully native small perturbation is not exactly clean-returning.
        perturbation = [('CX', dirty, clean)]
        perturbation += _flat(commutator_word(3, clean))
        perturbation += [('CX', dirty, clean)]
        w = _word_matrix(width, perturbation)
        a_one = v @ z_dirty @ v.conj().T
        fibres, eta = {}, 0.0
        leakage = 0.0
        for a in (0, 1):
            for b in (0, 1):
                p = (x_clean if a else identity) @ (z_dirty if b else identity)
                r = x_clean if a else identity
                u = p @ (a_one if b else identity) @ w @ r
                desired = x_target if a and b else identity
                error = np.linalg.norm(u @ j - desired @ j, 2)
                eta = max(eta, error)
                leakage = max(leakage, np.linalg.norm((u @ j)[1 << clean:], 2))
                inverse_error = np.linalg.norm(u.conj().T @ j - desired.conj().T @ j, 2)
                self.assertLessEqual(inverse_error, error + 3e-11)
                fibres[a, b] = u
        self.assertGreater(leakage, 1e-8)
        self.assertLess(eta, 1 / 64)
        for a in (0, 1):
            # H consists of the b direction. Equation (12), including the
            # actual inverse, cancels the fixed perturbation on full space.
            base_p = x_clean if a else identity
            b_family = base_p @ a_one @ base_p.conj().T
            extracted = z_dirty @ fibres[a, 1] @ fibres[a, 0].conj().T
            np.testing.assert_allclose(extracted, b_family, atol=5e-11, rtol=0)
            expected = z_dirty @ (x_target if a else identity)
            self.assertLessEqual(np.linalg.norm(extracted @ j - expected @ j, 2),
                                 2 * eta + 5e-11)
        # Opposite left clean syndromes can have different affine slopes.
        # Ignoring this distinction would falsely strengthen the abstract
        # three-band theorem from quadratic to affine.
        self.assertGreater(np.linalg.norm((z_dirty @ x_target - z_dirty) @ j, 2), 1.9)

    def test_projective_pauli_label_gap(self):
        paulis = [_word_matrix(2, [('X', 0)] * (code & 1)
                              + [('X', 1)] * ((code >> 1) & 1)
                              + [('Z', 0)] * ((code >> 2) & 1)
                              + [('Z', 1)] * ((code >> 3) & 1))
                  for code in range(16)]
        for left in range(16):
            for right in range(left):
                # Opposite-eigenvalue quotients attain their optimal phase
                # at one of these four scalars.
                distance = min(np.linalg.norm(paulis[left] - phase * paulis[right], 2)
                               for phase in (1, -1, 1j, -1j))
                self.assertAlmostEqual(distance, np.sqrt(2), places=12)

    def test_quadratic_polar_certificate_and_excluded_witnesses(self):
        # Every product of two linear forms has polar rank at most two,
        # including dependent/zero forms and non-coordinate directions.
        for left in range(16):
            for right in range(16):
                table = [((left & x).bit_count() & 1)
                         & ((right & x).bit_count() & 1) for x in range(16)]
                self.assertLessEqual(_degree(table), 2)
                self.assertLessEqual(_polar_rank(table), 2)
        cubic = [int((x & 7) == 7) for x in range(8)]
        independent_pairs = [((x & 1) & ((x >> 1) & 1))
                             ^ (((x >> 2) & 1) & ((x >> 3) & 1)) for x in range(16)]
        self.assertEqual(_degree(cubic), 3)
        self.assertEqual(_degree(independent_pairs), 2)
        self.assertEqual(_polar_rank(independent_pairs), 4)


if __name__ == '__main__':
    unittest.main()
