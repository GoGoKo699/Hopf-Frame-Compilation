"""Small full-input checks for a carried flag/source code.

The original chain loader, including its literal native scalar phases,
supplies this code. Its cheap physical width transition is distinct from
the star loader's eigenbasis bridge. Matrices have dimension at most 128;
the checks do not establish a reusable grouped compiler or a T lower bound.
"""
from __future__ import annotations

import unittest

import numpy as np


ATOL = 4e-12
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
S = np.diag([1, 1j])
T = np.diag([1, np.exp(1j * np.pi / 4)])


def _pauli(width, entries):
    result = np.ones((1, 1), dtype=complex)
    for index in range(width):
        result = np.kron(result, entries.get(index, I2))
    return result


def _gamma(width, index):
    return _pauli(width, {**{earlier: Z for earlier in range(index)}, index: X})


def _literal_split(width, index):
    """Clifford-conjugated T dagger, without deleting its scalar phase."""
    pauli = _pauli(width, {index: Y, index + 1: X})
    phase = np.exp(-1j * np.pi / 4)
    return ((1 + phase) * np.eye(1 << width) + (1 - phase) * pauli) / 2


def _chain_loader(width, active):
    result = np.eye(1 << width, dtype=complex)
    for index in range(active - 1):
        result = _literal_split(width, index) @ result
    return result


def _source(width, active):
    loader = _chain_loader(width, active)
    return loader @ _gamma(width, 0) @ loader.conj().T


def _boundary(width, active):
    j_x = np.vstack([np.eye(1 << width), _gamma(width, 0)]) / np.sqrt(2)
    return np.kron(I2, _chain_loader(width, active)) @ j_x


def _shrink_word(width, larger, smaller):
    """Literal tail removal: chronological highest-index inverse first."""
    result = np.eye(1 << width, dtype=complex)
    for index in range(larger - 2, smaller - 2, -1):
        result = _literal_split(width, index).conj().T @ result
    return np.kron(I2, result)


def _syndrome_word(source):
    stabilizer = np.kron(X, source)
    size = len(stabilizer)
    zero = np.zeros((size, size), dtype=complex)
    controlled = np.block([[np.eye(size), zero], [zero, stabilizer]])
    hadamard = np.kron(H, np.eye(size))
    return hadamard @ controlled @ hadamard, stabilizer


class CorrelatedPrecisionCarryTests(unittest.TestCase):
    def test_literal_loader_phase_and_code_projector(self):
        # This native two-wire word maps Z0 to Y0 X1. Comparing the whole
        # matrix catches the phase lost by using only exp(i pi P / 8).
        cnot_1_to_0 = np.zeros((4, 4), dtype=complex)
        for source in range(4):
            target = source ^ (2 if source & 1 else 0)
            cnot_1_to_0[target, source] = 1
        basis = np.kron(S @ H, H) @ cnot_1_to_0
        native = basis @ np.kron(T.conj().T, I2) @ basis.conj().T
        np.testing.assert_allclose(native, _literal_split(2, 0), atol=ATOL, rtol=0)

        width = 5
        identity = np.eye(1 << width)
        for active in (2, 3, 4, 5):
            with self.subTest(active=active):
                source = _source(width, active)
                expected = sum(2 ** (-(index + 1) / 2) * _gamma(width, index)
                               for index in range(active - 1))
                expected += 2 ** (-(active - 1) / 2) * _gamma(width, active - 1)
                np.testing.assert_allclose(source, expected, atol=ATOL, rtol=0)
                boundary = _boundary(width, active)
                np.testing.assert_allclose(boundary.conj().T @ boundary, identity,
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(boundary @ boundary.conj().T,
                                           (np.eye(2 << width) + np.kron(X, source)) / 2,
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(np.kron(I2, source) @ boundary,
                                           np.kron(X, identity) @ boundary,
                                           atol=ATOL, rtol=0)

    def test_three_unequal_widths_preserve_all_dirty_and_reference_inputs(self):
        width = 5
        first = _shrink_word(width, 5, 4)
        second = _shrink_word(width, 4, 2)
        start, middle, finish = (_boundary(width, active) for active in (5, 4, 2))
        np.testing.assert_allclose(first @ start, middle, atol=ATOL, rtol=0)
        np.testing.assert_allclose(second @ first @ start, finish, atol=ATOL, rtol=0)

        # The input is arbitrary on the entire maximum-width dirty pool,
        # jointly entangled with a two-dimensional reference; no tail is zero.
        rng = np.random.default_rng(707)
        state = rng.normal(size=(1 << width, 2)) + 1j * rng.normal(size=(1 << width, 2))
        state /= np.linalg.norm(state)
        actual = second @ first @ start @ state
        np.testing.assert_allclose(actual, finish @ state, atol=ATOL, rtol=0)
        np.testing.assert_allclose(first.conj().T @ second.conj().T @ actual,
                                   start @ state, atol=ATOL, rtol=0)

    def test_released_tail_is_arbitrary_dirty_workspace(self):
        width, active = 5, 3
        boundary = _boundary(width, active)
        projector = boundary @ boundary.conj().T
        for tail_operator in (_pauli(width, {3: X}), _pauli(width, {4: Z}),
                              _pauli(width, {3: Y, 4: X})):
            full = np.kron(I2, tail_operator)
            np.testing.assert_allclose(full @ boundary, boundary @ tail_operator,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(full @ projector, projector @ full,
                                       atol=ATOL, rtol=0)

        # The corresponding last source wire was not available before shrinking.
        old = _boundary(width, width)
        old_projector = old @ old.conj().T
        full = np.kron(I2, _pauli(width, {4: Z}))
        self.assertGreater(np.linalg.norm(full @ old_projector - old_projector @ full, 2),
                           0.1)

    def test_realizable_group_mask_has_uniform_code_leakage(self):
        width = 5
        source = _source(width, width)
        boundary = _boundary(width, width)
        identity = np.eye(1 << width)
        mask = _pauli(width, {1: Z, 2: Z, 3: Z})
        masked_source = mask @ source @ mask
        coefficient = 1 / 8
        np.testing.assert_allclose((source @ masked_source + masked_source @ source) / 2,
                                   coefficient * identity, atol=ATOL, rtol=0)

        # Actual one-level group, native C=I, and its special forward atom.
        labels = 32  # 2**ceil(log2(8*s+12)) at s=1.
        sine = 1 / (8 * labels * np.sqrt(2))
        cosine = np.sqrt(1 - sine * sine)
        target = np.array([[cosine, -sine], [sine, cosine]])
        coarse_allowance = 1 / (4 * labels * np.sqrt(2))
        self.assertLessEqual(np.linalg.norm(target - np.eye(2), 2), coarse_allowance)
        self.assertAlmostEqual(labels * np.sqrt(2) * target[1, 0], coefficient)

        projector = boundary @ boundary.conj().T
        queried = np.kron(I2, mask) @ boundary
        changed_projector = (np.eye(2 << width) + np.kron(X, masked_source)) / 2
        np.testing.assert_allclose(changed_projector @ queried, queried, atol=ATOL, rtol=0)
        leakage = (np.eye(2 << width) - projector) @ queried
        np.testing.assert_allclose(leakage.conj().T @ leakage, 7 / 16 * identity,
                                   atol=ATOL, rtol=0)
        substitution_error = (np.kron(I2, source) - np.kron(X, identity)) @ queried
        np.testing.assert_allclose(substitution_error.conj().T @ substitution_error,
                                   7 / 4 * identity, atol=ATOL, rtol=0)

    def test_complete_syndrome_extraction_and_actual_inverse_recover_source(self):
        source = _source(5, 5)
        word, stabilizer = _syndrome_word(source)
        size = len(stabilizer)
        identity = np.eye(size)
        injection = np.vstack([identity, np.zeros_like(identity)])
        syndrome = np.vstack([(identity + stabilizer) / 2, (identity - stabilizer) / 2])
        np.testing.assert_allclose(word @ injection, syndrome, atol=ATOL, rtol=0)
        np.testing.assert_allclose(word.conj().T @ word, np.eye(2 * size), atol=ATOL, rtol=0)
        reflected = word.conj().T @ np.kron(Z, identity) @ word @ injection
        np.testing.assert_allclose(reflected, injection @ stabilizer, atol=ATOL, rtol=0)

    def test_approximate_extractor_bound_includes_nonzero_flag_leakage(self):
        source = _source(4, 4)
        word, stabilizer = _syndrome_word(source)
        size = len(stabilizer)
        identity = np.eye(size)
        injection = np.vstack([identity, np.zeros_like(identity)])
        syndrome = word @ injection
        angle = 0.07
        perturbation = np.kron(np.cos(angle) * I2 + 1j * np.sin(angle) * X, identity)
        actual = perturbation @ word
        error = np.linalg.norm(actual @ injection - syndrome, 2)
        reflected = actual.conj().T @ np.kron(Z, identity) @ actual @ injection
        self.assertLessEqual(np.linalg.norm(reflected - injection @ stabilizer, 2),
                             2 * error + ATOL)
        self.assertGreater(np.linalg.norm(reflected[size:], 2), 0.1)


if __name__ == "__main__":
    unittest.main()
