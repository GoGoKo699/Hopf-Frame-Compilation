"""Small full-space witnesses for one-query orbit gaps and their premises.

The dimension and covariance proofs are analytic. These fixtures retain
literal phases, both carrier inputs and the actual inverse; the two-query
example is not an all-size source compiler.
"""
import math
import unittest

import numpy as np

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_coherent_feedback_sources import _shared_ccz
from tests.test_joint_lookup_pauli_source import _inverse
from tests.test_minimal_tree_probe import _tree_forward


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
H = (X + Z) / math.sqrt(2)
GAMMA = math.sqrt(2 - math.sqrt(2))


def _rank(rows):
    pivots = {}
    for row in rows:
        while row:
            bit = row.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = row
                break
            row ^= pivots[bit]
    return len(pivots)


def _xor_permutation(rows, control_width, offset=0, spectators=0):
    """All basis columns, with little-endian control, target, spectator banks."""
    width = control_width + len(rows) + spectators
    result = []
    for column in range(1 << width):
        target_delta = offset
        for target, row in enumerate(rows):
            target_delta ^= ((row & column).bit_count() & 1) << target
        result.append(column ^ (target_delta << control_width))
    return result


def _matrix(permutation):
    return np.eye(len(permutation), dtype=complex)[:, permutation]


def _cnot(control, target):
    return _matrix([column ^ (((column >> control) & 1) << target)
                    for column in range(4)])


def _target_label(p, address):
    m = 1 << p
    v = sum(1 << subset for subset in range(m) if subset & address == subset)
    av = _tree_forward(address, v, p)
    atv = sum(((_tree_forward(address, 1 << j, p) & v).bit_count() & 1) << j
              for j in range(m))
    return v | (av << m), atv | (v << m)


def _pairing(left, right):
    a, b = left
    c, d = right
    return ((a & d).bit_count() + (b & c).bit_count()) & 1


class QueryOrbitSpectralGapTests(unittest.TestCase):
    def test_linear_shears_have_a_strict_majority_eigenspace_on_all_inputs(self):
        # Every 2-by-3 binary matrix, with an arbitrary spectator.
        for encoded in range(1 << 6):
            rows = (encoded & 7, encoded >> 3)
            permutation = _xor_permutation(rows, 3, spectators=1)
            dimension = len(permutation)
            fixed = sum(row == column for column, row in enumerate(permutation))
            self.assertEqual(fixed, dimension >> _rank(rows))
            self.assertTrue(all(permutation[row] == column
                                for column, row in enumerate(permutation)))
            plus_dimension = (dimension + fixed) // 2
            self.assertGreater(plus_dimension, dimension // 2)

    def test_majority_intersection_survives_nonclifford_probe_and_decoder(self):
        relative = _cnot(0, 1)
        # Wire zero is genuinely non-Clifford; neither input is initialized.
        axis = math.cos(math.pi / 8) * Z + math.sin(math.pi / 8) * X
        probe = np.kron(Z, axis)
        identity = np.eye(4)
        constraints = np.vstack((relative - identity, probe - identity))
        _, singular, vh = np.linalg.svd(constraints)
        self.assertLess(singular[-1], 1e-13)
        common = vh[-1].conj()
        decoder = np.diag(np.exp(1j * np.array([0, 1, 3, 7]) / 5)) @ np.kron(H, H)
        state = decoder @ common
        images = (decoder @ probe @ decoder.conj().T,
                  decoder @ relative @ probe @ relative @ decoder.conj().T)
        for image in images:
            np.testing.assert_allclose(image @ state, state, atol=1e-13)
        targets = (np.kron(I2, X), np.kron(I2, Z))
        errors = [np.linalg.norm(image - target, ord=2)
                  for image, target in zip(images, targets)]
        self.assertGreaterEqual(max(errors), GAMMA - 1e-13)

    def test_majority_lemma_threshold_is_sharp_in_operator_norm(self):
        # L=I gives identical source images; their best midpoint is H.
        for target in (X, Z):
            self.assertAlmostEqual(np.linalg.norm(H - target, ord=2), GAMMA)
        _, vectors = np.linalg.eigh(H)
        state = vectors[:, -1]
        expectations = [np.vdot(state, target @ state).real for target in (X, Z)]
        self.assertAlmostEqual(sum(value ** 2 for value in expectations), 1)
        self.assertGreater(2 * (1 - (GAMMA - 0.01) ** 2 / 2) ** 2, 1)

    def test_affine_offsets_and_projective_phases_keep_address_covariance(self):
        def query(address):
            first, second = address & 1, (address >> 1) & 1
            rows = (1 ^ first ^ (second << 1), 2 ^ (first << 1) ^ second)
            offset = 1 ^ first ^ (second << 1)
            return _matrix(_xor_permutation(rows, 2, offset))

        base = query(0)
        linear = [base @ query(address) for address in range(4)]
        # A scalar cochain makes the representatives projective without
        # changing any adjoint action.
        representatives = [np.exp(1j * address.bit_count() ** 2 / 7) * word
                           for address, word in enumerate(linear)]
        axis = math.cos(math.pi / 8) * Z + math.sin(math.pi / 8) * X
        probe = np.kron(np.eye(8), axis)
        orbit = [word @ probe @ word.conj().T for word in representatives]
        for address in range(4):
            for shift in range(4):
                np.testing.assert_array_equal(linear[address] @ linear[shift],
                                              linear[address ^ shift])
                translated = representatives[shift] @ orbit[address] @ representatives[shift].conj().T
                np.testing.assert_allclose(translated, orbit[address ^ shift], atol=1e-13)

    def test_tree_target_changes_commutation_under_address_translation(self):
        for p in range(2, 8):
            labels = [_target_label(p, address) for address in (0, 1, 2, 3)]
            self.assertEqual(_pairing(labels[0], labels[1]), 1)
            self.assertEqual(_pairing(labels[2], labels[3]), 0)
        # If covariance moved an anticommuting pair within epsilon of a
        # commuting pair, a common eigenvector would force epsilon>=gamma.
        delta = GAMMA / 2 - 0.001
        self.assertGreater(2 * (1 - (2 * delta) ** 2 / 2) ** 2, 1)

    def test_translation_offset_is_a_boundary_to_the_majority_premise(self):
        first, second = H, X @ H @ X
        np.testing.assert_allclose(first @ second + second @ first, np.zeros((2, 2)), atol=1e-14)
        self.assertEqual(np.trace(X), 0)
        self.assertEqual(np.linalg.matrix_rank(X - I2), 1)
        s = np.diag([1, 1j])
        t_dagger = np.diag([1, np.exp(-1j * math.pi / 4)])
        # Literal S H T-dagger H S-dagger, with its actual adjoint.
        decoder = s @ H @ t_dagger @ H @ s.conj().T
        np.testing.assert_allclose(decoder @ first @ decoder.conj().T, Z, atol=1e-14)
        np.testing.assert_allclose(decoder @ second @ decoder.conj().T, X, atol=1e-14)

    def test_two_opposite_queries_have_literal_order_three_and_actual_inverse(self):
        forward, backward = _cnot(0, 1), _cnot(1, 0)
        relative = backward @ forward
        np.testing.assert_array_equal(relative @ relative @ relative, np.eye(4))
        np.testing.assert_array_equal(relative.conj().T, relative @ relative)
        self.assertFalse(np.array_equal(relative, relative.conj().T))
        self.assertEqual(4 - np.linalg.matrix_rank(relative - np.eye(4)), 2)
        family = (np.kron(I2, Y), np.kron(Y, Z), np.kron(Y, X))
        for index, probe in enumerate(family):
            np.testing.assert_array_equal(relative @ probe @ relative.conj().T,
                                          family[(index + 1) % 3])
            for other in family[index + 1:]:
                np.testing.assert_array_equal(probe @ other + other @ probe, np.zeros((4, 4)))

        # Address wire 0 controls both XOR stages; all eight columns and
        # their literal Clifford+T inverses are checked, with no helper.
        word = []
        for control, target in ((1, 2), (2, 1)):
            word += ([('H', target)] + _shared_ccz(0, [(control, target)])
                     + [('H', target)])
        actual = native_columns(3, word)
        expected = [[CZ for _ in range(8)] for _ in range(8)]
        for column in range(8):
            row = column
            if column & 1:
                row ^= ((row >> 1) & 1) << 2
                row ^= ((row >> 2) & 1) << 1
            expected[row][column] = CO
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(3, _inverse(word)), adjoint(actual))


if __name__ == '__main__':
    unittest.main()
