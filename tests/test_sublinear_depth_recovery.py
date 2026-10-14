"""Bounded all-column audits of the compressed early-frame recovery.

The two-target source uses independently written subtree vectors and numerical
unitary completions. QR is only a finite verification fixture, never a claimed
compiler primitive. This module emits neither the derived Boolean tables nor a
native GKW loader, dirty query, or scalable complete-frame circuit. It checks
the marker embedding, recovery algebra, actual inverses, and protected-sector
phase/return contract; their native resource bounds remain analytic.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.frames import real_frame_matrix


ATOL = 2e-11
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)


def _complete_column(vector):
    """A small verification-only unitary having the specified first column."""
    vector = np.asarray(vector, dtype=complex)
    unitary, triangular = np.linalg.qr(
        np.column_stack((vector, np.eye(len(vector), dtype=complex)))
    )
    unitary[:, 0] *= triangular[0, 0] / abs(triangular[0, 0])
    return unitary


def _block_diagonal(*blocks):
    total = sum(len(block) for block in blocks)
    result = np.zeros((total, total), dtype=complex)
    offset = 0
    for block in blocks:
        result[offset:offset + len(block), offset:offset + len(block)] = block
        offset += len(block)
    return result


def _source(theta, *, wrong_left_tag=False):
    """Tag-major order; four data modes and four tag values, including padding."""
    a, b, c = theta
    ca, sa = np.cos(a), np.sin(a)
    cb, sb = np.cos(b), np.sin(b)
    cc, sc = np.cos(c), np.sin(c)
    # Independent recursive state/complement formulas, in computational order.
    root = np.array([ca * cb, ca * sb, sa * cc, sa * sc])
    root_marker = np.array([-sa * cb, -sa * sb, ca * cc, ca * sc])
    at_depth_one = _block_diagonal(
        _complete_column([-sb, cb]), _complete_column([-sc, cc])
    )
    loader = _block_diagonal(
        _complete_column(root), _complete_column(root_marker),
        at_depth_one, np.eye(4),
    )
    tags = (0, 1 if wrong_left_tag else 2, 1, 2)
    masks = (0, 2, 1, 0)
    compression = np.zeros((16, 16), dtype=complex)
    for old in range(16):
        tag, data = divmod(old, 4)
        new_tag = tag ^ tags[data]
        new_data = data ^ masks[new_tag]
        compression[4 * new_tag + new_data, old] = 1
    return loader @ compression, compression


def _damped_query(source, rounds):
    # K=4 padded tags. Either rounds=2 (the sufficient formula) or a larger
    # value is legal, with the attenuation adjusted to the same exact angle.
    theta = np.pi / (4 * rounds + 2)
    gamma = 2 * np.sin(theta)
    other = np.sqrt(1 - gamma * gamma)
    attenuation = np.array([[gamma, -other], [other, gamma]])
    return np.kron(attenuation, np.kron(np.kron(H, H), np.eye(4)) @ source)


def _recover(query, rounds, *, inverse=None):
    reflection = np.eye(32, dtype=complex)
    reflection[:4, :4] *= -1
    inverse = query.conj().T if inverse is None else inverse
    iterate = -query @ reflection @ inverse @ reflection
    return np.linalg.matrix_power(iterate, rounds) @ query


def _protected(query, rounds, *, global_minus=False):
    """Guard, complete physical work bank, and the actual final guard erasure."""
    eye = np.eye(32, dtype=complex)
    guarded = _block_diagonal(eye, query)
    reflection = np.eye(64, dtype=complex)
    reflection[32:36, 32:36] *= -1
    sign = -np.eye(64) if global_minus else _block_diagonal(eye, -eye)
    iterate = sign @ guarded @ reflection @ guarded.conj().T @ reflection
    interior = np.linalg.matrix_power(iterate, rounds) @ guarded
    # H is the highest bit. B=0 means the other work bits (tag and damping)
    # are zero; the four data labels are arbitrary on that sector.
    predicate = np.eye(64, dtype=complex)
    for data in range(4):
        predicate[[data, 32 + data]] = predicate[[32 + data, data]]
    return predicate @ interior @ predicate


class SublinearDepthRecoveryTests(unittest.TestCase):
    def test_compressed_marker_map_and_all_columns(self):
        fixtures = (
            (0.37, -0.61, 1.12),  # Unequal signed angles.
            (0.0, 0.83, -1.21),  # Singular state chart; all columns still matter.
            (np.pi / 2, np.pi, 0.0),
        )
        for angles in fixtures:
            with self.subTest(angles=angles):
                source, compression = _source(angles)
                expected = real_frame_matrix(angles)
                np.testing.assert_allclose(source.conj().T @ source, np.eye(16), atol=ATOL)
                np.testing.assert_array_equal(
                    np.argmax(abs(compression[:, :4]), axis=0), [0, 8, 4, 10]
                )
                encoded = source[:, :4]
                for column, tag in enumerate((0, 2, 1, 2)):
                    reference = np.zeros(16, dtype=complex)
                    reference[4 * tag:4 * tag + 4] = expected[:, column]
                    np.testing.assert_allclose(encoded[:, column], reference, atol=ATOL)
                left = np.kron(np.ones((1, 4)) / 2, np.eye(4))
                np.testing.assert_allclose(left @ encoded, expected / 2, atol=ATOL)
                for rounds in (2, 3):
                    recovered = _recover(_damped_query(source, rounds), rounds)[:, :4]
                    target = np.zeros((32, 4), dtype=complex)
                    target[:4] = expected
                    np.testing.assert_allclose(recovered, target, atol=ATOL)

    def test_wrong_marker_tag_can_preserve_root_but_corrupt_frame(self):
        angles = (0.37, -0.61, 1.12)
        correct, _ = _source(angles)
        wrong, _ = _source(angles, wrong_left_tag=True)
        np.testing.assert_allclose(wrong[:, 0], correct[:, 0], atol=ATOL)
        self.assertGreater(np.linalg.norm(wrong[:, 1] - correct[:, 1]), 0.1)
        recovered = _recover(_damped_query(wrong, 2), 2)[:, :4]
        target = np.zeros((32, 4), dtype=complex)
        target[:4] = real_frame_matrix(angles)
        self.assertGreater(np.linalg.norm(recovered - target, 2), 0.1)

    def test_recovery_uses_actual_inverse_for_different_completions(self):
        angles = (0.41, 0.0, -0.87)
        source, _ = _source(angles)
        # Same source columns, deliberately very different occupied-input action.
        change = np.diag([1] * 4 + [np.exp(0.73j)] * 12)
        alternate = source @ change
        np.testing.assert_allclose(alternate[:, :4], source[:, :4], atol=ATOL)
        self.assertGreater(np.linalg.norm(alternate - source, 2), 0.7)
        query = _damped_query(alternate, 3)
        baseline = _damped_query(source, 3)
        target = np.zeros((32, 4), dtype=complex)
        target[:4] = real_frame_matrix(angles)
        np.testing.assert_allclose(_recover(query, 3)[:, :4], target, atol=ATOL)
        wrong = _recover(query, 3, inverse=baseline.conj().T)[:, :4]
        self.assertGreater(np.linalg.norm(wrong - target, 2), 0.1)

    def test_protected_recovery_keeps_inactive_relative_phase(self):
        angles = (0.37, -0.61, 1.12)
        source, _ = _source(angles)
        query = _damped_query(source, 3)  # Odd count exposes an omitted active-only minus.
        expected = np.zeros((64, 32), dtype=complex)
        expected[:32] = _block_diagonal(real_frame_matrix(angles), np.eye(28))
        actual = _protected(query, 3)[:, :32]
        np.testing.assert_allclose(actual, expected, atol=ATOL)
        wrong = _protected(query, 3, global_minus=True)[:, :32]
        np.testing.assert_allclose(wrong[:, :4], expected[:, :4], atol=ATOL)
        self.assertAlmostEqual(np.linalg.norm(wrong - expected, 2), 2.0, places=10)
        # A coherent superposition detects the relative phase, which a check
        # confined to either basis sector could mistakenly dismiss as global.
        vector = np.zeros(32, dtype=complex)
        vector[0], vector[4] = 1 / np.sqrt(2), 1j / np.sqrt(2)
        self.assertGreater(np.linalg.norm((wrong - expected) @ vector), 1.4)

    def test_source_error_includes_work_leakage_and_final_flag_return(self):
        angles = (0.0, 0.83, -1.21)
        source, _ = _source(angles)
        angle = 0.006
        c, s = np.cos(angle), np.sin(angle)
        mixing = np.eye(4)
        mixing[np.ix_([0, 3], [0, 3])] = [[c, -s], [s, c]]
        perturbed = np.kron(mixing, np.eye(4)) @ source
        delta = np.linalg.norm((perturbed - source)[:, :4], 2)
        self.assertAlmostEqual(delta, 2 * np.sin(angle / 2), places=12)
        actual = _protected(_damped_query(perturbed, 3), 3)[:, :32]
        target = np.zeros((64, 32), dtype=complex)
        target[:32] = _block_diagonal(real_frame_matrix(angles), np.eye(28))
        error = np.linalg.norm(actual - target, 2)
        self.assertGreater(error, 1e-4)
        # Here the chosen perturbation also has full-operator distance delta,
        # so direct seven-call telescoping is sharper than the general 6delta
        # near-isometry-completion bound used for an arbitrary loader.
        self.assertLessEqual(error, 7 * delta + ATOL)
        self.assertGreater(np.linalg.norm(actual[32:, :4], 2), 1e-4)
        np.testing.assert_allclose(actual[:, 4:], target[:, 4:], atol=ATOL)


if __name__ == "__main__":
    unittest.main()
