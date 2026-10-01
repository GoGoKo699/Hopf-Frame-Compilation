"""Four real modes and a small native spectator-borrowing diagnostic.

The Clifford and Cayley checks retain literal phases. The q=2 circuit is
not a certified-precision instance: it audits actual native gates on all
dirty columns, exact borrowed-helper return, and shared-flag telescoping.
No matrix has dimension greater than 128.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_one_clean_compiler import (
    _amplification_word, _coefficients, _encoded_signs, _literal_mask,
    _paired_source_word, _paired_weights, _sandwich_word,
)
from tests.test_operator_source_compiler import (
    ATOL, I2, X, Y, Z, _adjoint, _expand_toffolis, _pauli, _word_matrix,
)


def _magic_word(high, low):
    """M=CX H_high S_high S_low CZ, in chronological native order."""
    return [('S', high), ('S', low), ('H', low), ('CX', high, low),
            ('H', low), ('H', high), ('CX', high, low)]


def _regular_factors(generator):
    k = generator
    first = np.array([-k[0, 1] - k[2, 3], -k[0, 3] - k[1, 2],
                      -k[0, 2] + k[1, 3]]) / 2
    second = np.array([-k[0, 1] + k[2, 3], k[0, 3] - k[1, 2],
                       -k[0, 2] - k[1, 3]]) / 2
    r2, s2 = first @ first, second @ second
    denominator = np.sqrt((1 + r2 + s2) ** 2 - 4 * r2 * s2)
    a = ((1 - r2 + s2) * I2
         + 2j * sum(value * axis for value, axis in zip(first, (X, Y, Z))))
    b = ((1 + r2 - s2) * I2
         + 2j * sum(value * axis for value, axis in zip(second, (X, Y, Z))))
    return a / denominator, b / denominator, denominator


def _tree_generator(root, left, right):
    a, b, c = root, left, right
    return np.array([[0, -b, -a, -a * c],
                     [b, 0, -a * b, -a * b * c],
                     [a, a * b, 0, -c],
                     [a * c, a * b * c, c, 0]], dtype=float)


def _plane_rotation(first, second, angle):
    result = np.eye(4)
    c, s = np.cos(angle), np.sin(angle)
    result[np.ix_([first, second], [first, second])] = [[c, -s], [s, c]]
    return result


def _addressed_z_word(q, target, helper, address, flag, angle):
    """Two conditional literal phases implement an addressed Rz.

    Each flagged source center has controls (target,address,flag), and its
    C3X echo really borrows helper. The opposite target is never an address
    or source-core wire. Both helper inputs and both flag inputs are valid.
    """
    word = []
    for target_value in (0, 1):
        theta = angle if target_value == 0 else -angle
        signs = _encoded_signs(q, theta)
        primitive = _sandwich_word(
            q, _literal_mask(signs), None, flag, phase=True,
            predicate=(target, address), helper=helper)
        flips = [('X', target)] if target_value == 0 else []
        word += flips + _amplification_word(primitive, flag) + flips
    return word


def _addressed_block(width, target, address, block):
    result = np.eye(1 << width, dtype=complex)
    for basis in range(1 << width):
        if (basis >> target) & 1 or not (basis >> address) & 1:
            continue
        indices = [basis, basis | (1 << target)]
        result[np.ix_(indices, indices)] = block
    return result


def _addressed_rotation(width, target, address, angle, axis):
    rotation = np.cos(angle) * I2 - 1j * np.sin(angle) * axis
    return _addressed_block(width, target, address, rotation)


def _encoded_accepted_phase(q, angle):
    """Closed-form compression, independent of the native gate evaluator."""
    weights, fixed = _paired_weights(q)
    diagonal = []
    for theta in (angle, -angle):
        s, p = _coefficients(weights, fixed, _encoded_signs(q, theta))
        radius2 = s * s + p * p
        diagonal.append((5 - 20 * radius2 + 16 * radius2 ** 2) * (s - 1j * p))
    return np.diag(diagonal)


class NativeCayleyTests(unittest.TestCase):
    def test_literal_seven_clifford_magic_basis_and_orientation(self):
        word = _magic_word(1, 0)
        self.assertEqual(len(word), 7)
        self.assertFalse(any(gate[0] in ('T', 'TDG') for gate in word))
        actual = _word_matrix(2, word)
        expected = np.array([[1, 0, 1j, 0], [0, 1j, 0, 1],
                             [0, 1j, 0, -1], [1, 0, -1j, 0]]) / np.sqrt(2)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.det(actual), 1, delta=ATOL)
        np.testing.assert_allclose(_word_matrix(2, word + _adjoint(word)),
                                   np.eye(4), atol=ATOL, rtol=0)
        # This table distinguishes M W M^dagger from its reversed convention.
        for source, desired in (
                (np.kron(Y, I2), -np.kron(I2, Z)),
                (np.kron(X, Y), np.kron(I2, Y)),
                (np.kron(Z, Y), -np.kron(I2, X)),
                (np.kron(I2, Y), -np.kron(X, I2)),
                (np.kron(Y, X), -np.kron(Y, I2)),
                (np.kron(Y, Z), -np.kron(Z, I2))):
            np.testing.assert_allclose(actual @ source @ actual.conj().T,
                                       desired, atol=ATOL, rtol=0)

    def test_regular_cayley_factors_cover_zero_axes_and_complete_tree(self):
        magic = _word_matrix(2, _magic_word(1, 0))
        generators = []
        # Zero generator, either zero factor axis, equal axis norms, and a
        # generic real generator; the formulas never divide by an axis norm.
        for upper in ((0, 0, 0, 0, 0, 0), (.2, 0, 0, 0, 0, .2),
                      (.2, 0, 0, 0, 0, -.2), (.4, 0, 0, 0, 0, 0),
                      (.31, -.2, .17, .43, -.11, -.28)):
            generator = np.zeros((4, 4))
            generator[np.triu_indices(4, 1)] = upper
            generators.append(generator - generator.T)
        for tangents in ((.13, -.24, .31), (0, 0, .5), (2, -3, 4)):
            generator = _tree_generator(*tangents)
            generators.append(generator)
            a, b, c = (2 * np.arctan(t) for t in tangents)
            tree = (_plane_rotation(2, 3, c) @ _plane_rotation(0, 1, b)
                    @ _plane_rotation(0, 2, a))
            cayley = (np.eye(4) + generator) @ np.linalg.inv(np.eye(4) - generator)
            np.testing.assert_allclose(cayley, tree, atol=ATOL, rtol=0)
        for generator in generators:
            with self.subTest(generator=generator.tolist()):
                a, b, denominator = _regular_factors(generator)
                self.assertGreaterEqual(denominator, 1 - ATOL)
                for factor in (a, b):
                    np.testing.assert_allclose(factor.conj().T @ factor, I2,
                                               atol=ATOL, rtol=0)
                    self.assertAlmostEqual(np.linalg.det(factor), 1, delta=ATOL)
                desired = ((np.eye(4) + generator)
                           @ np.linalg.inv(np.eye(4) - generator))
                np.testing.assert_allclose(magic.conj().T @ np.kron(a, b) @ magic,
                                           desired, atol=ATOL, rtol=0)

    def test_native_addressed_spectator_return_and_shared_flag_telescope(self):
        q, low, high, address, flag, width = 2, 3, 4, 5, 6, 7
        angles = (.37, -.62)
        raw_a = _addressed_z_word(q, high, low, address, flag, angles[0])
        raw_b = ([('H', low)]
                 + _addressed_z_word(q, low, high, address, flag, angles[1])
                 + [('H', low)])
        native_a, native_b = _expand_toffolis(raw_a), _expand_toffolis(raw_b)
        # Four amplified phase primitives: 45 source calls per primitive,
        # two loader appearances per source. Predicate echoes contribute
        # 135 exact Toffolis per primitive, at seven T gates apiece.
        source_calls = 4 * 45
        source_loader_t = sum(gate[0] in ('T', 'TDG')
                              for gate in _paired_source_word(q))
        tof_count = sum(gate[0] == 'CCX' for gate in raw_a + raw_b)
        self.assertEqual(tof_count, 4 * 135)
        self.assertEqual(sum(gate[0] in ('T', 'TDG')
                             for gate in native_a + native_b),
                         2 * source_calls * source_loader_t + 7 * tof_count)
        self.assertTrue(any(gate[0] == 'CCX' and gate[-1] == low for gate in raw_a))
        self.assertTrue(any(gate[0] == 'CCX' and gate[-1] == high for gate in raw_b))
        actual_a = _word_matrix(width, native_a)
        actual_b = _word_matrix(width, native_b)
        identity = np.eye(1 << width)
        inactive = [basis for basis in range(1 << width)
                    if not (basis >> address) & 1]
        for actual, helper in ((actual_a, low), (actual_b, high)):
            np.testing.assert_allclose(actual.conj().T @ actual, identity,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(actual[:, inactive], identity[:, inactive],
                                       atol=ATOL, rtol=0)
            # Commuting with both X and Z means identity on the entire helper
            # factor, including arbitrary coherence and occupied flag inputs.
            for axis in (X, Z):
                helper_pauli = _pauli(width, {helper: axis})
                np.testing.assert_allclose(actual @ helper_pauli,
                                           helper_pauli @ actual, atol=ATOL, rtol=0)

        embedding = identity[:, :1 << flag]
        desired_a = _addressed_rotation(flag, high, address, angles[0], Z)
        desired_b = _addressed_rotation(flag, low, address, angles[1], X)
        encoded_a = _encoded_accepted_phase(q, angles[0])
        encoded_b = _encoded_accepted_phase(q, angles[1])
        hadamard = (X + Z) / np.sqrt(2)
        compressed_a = _addressed_block(flag, high, address, encoded_a)
        compressed_b = _addressed_block(
            flag, low, address, hadamard @ encoded_b @ hadamard)
        # Unlike the telescope alone, these comparisons fix the encoded
        # phase sign and the intended two addressed factor actions. Rz(theta)
        # has target-sector phases exp(-i theta), exp(+i theta).
        np.testing.assert_allclose(embedding.conj().T @ actual_a @ embedding,
                                   compressed_a, atol=ATOL, rtol=0)
        np.testing.assert_allclose(embedding.conj().T @ actual_b @ embedding,
                                   compressed_b, atol=ATOL, rtol=0)
        error_a = actual_a @ embedding - embedding @ desired_a
        error_b = actual_b @ embedding - embedding @ desired_b
        for error, desired, compressed in ((error_a, desired_a, compressed_a),
                                            (error_b, desired_b, compressed_b)):
            expected_gram = (2 * np.eye(1 << flag) - desired.conj().T @ compressed
                             - compressed.conj().T @ desired)
            np.testing.assert_allclose(error.conj().T @ error, expected_gram,
                                       atol=ATOL, rtol=0)
            # These coarse q=2 errors are approximately 1.28536 and 1.28642.
            # Reversing either intended phase gives error greater than 1.32.
            # This is a sign-sensitive diagnostic, not a precision guarantee.
            self.assertLess(np.linalg.norm(error, ord=2), 1.30)
        self.assertGreater(np.linalg.norm((actual_a @ embedding)[1 << flag:]), .01)
        composed_error = (actual_b @ actual_a @ embedding
                          - embedding @ desired_b @ desired_a)
        np.testing.assert_allclose(composed_error,
                                   actual_b @ error_a + error_b @ desired_a,
                                   atol=ATOL, rtol=0)
        self.assertLessEqual(np.linalg.norm(composed_error, ord=2),
                             np.linalg.norm(error_a, ord=2)
                             + np.linalg.norm(error_b, ord=2) + ATOL)

        magic_word = _magic_word(high, low)
        magic = _word_matrix(width, magic_word)
        logical_magic = _word_matrix(flag, magic_word)
        output = magic.conj().T @ actual_b @ actual_a @ magic @ embedding
        desired = logical_magic.conj().T @ desired_b @ desired_a @ logical_magic
        full_error = output - embedding @ desired
        np.testing.assert_allclose(
            full_error, magic.conj().T @ composed_error @ logical_magic,
            atol=ATOL, rtol=0)
        np.testing.assert_allclose(output.conj().T @ output, np.eye(1 << flag),
                                   atol=ATOL, rtol=0)
        # The accepted logical input includes both address rows, both targets,
        # and all eight precision-core inputs. No flag reset or source reset
        # is inserted between the two changing borrowed-pool partitions.


if __name__ == '__main__':
    unittest.main()
