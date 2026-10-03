"""Bounded full-bank audits of a protected unary source across groups.

The q=4 preparation is the matrix of a literal native word on the entire
six-wire bank; its actual inverse is that matrix's adjoint. Predicates, signed shifts, and logical group actions
are reduced exact operators, not a scalable native group or QROM emitter.
Parameterized source perturbations are analytic unitary fixtures. Both
flags start at zero; inactive interior identity is asserted only for
H=0,h=0, with arbitrary bank state. Final H erasure acts on actual leakage.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_operator_source_compiler import H, _word_matrix
    from .test_unary_phase_gradient import _ideal_group, _native_gradient
except ImportError:
    from test_operator_source_compiler import H, _word_matrix
    from test_unary_phase_gradient import _ideal_group, _native_gradient


ATOL = 1e-10


class _ProtectedFixture:
    """Three logical targets, a six-wire bank, and precisely two flags."""

    def __init__(self):
        width, unary, word = _native_gradient(4)
        assert width == 6
        self.bank_width, self.unary = width, unary
        self.bank_dimension, self.logical_dimension = 1 << width, 8
        self.h, self.H = width + 3, width + 4
        self.dimension = 1 << (width + 5)
        self.indices = np.arange(self.dimension)
        self.bank = self.indices & (self.bank_dimension - 1)
        self.logical = (self.indices >> width) & 7
        self.native = _word_matrix(width, word)
        self.rows = [[1], [0, 1], [3, 2, 1, 0]]

    def columns(self, basis):
        basis = list(basis)
        state = np.zeros((self.dimension, len(basis)), dtype=complex)
        state[basis, np.arange(len(basis))] = 1
        return state

    def prepare(self, state, unitary):
        tensor = state.reshape(2, 2, 8, self.bank_dimension, -1)
        return np.einsum('ij,habjc->habic', unitary, tensor,
                         optimize=True).reshape(state.shape)

    def toggle_H(self, state):
        permutation = self.indices ^ ((self.bank == 0).astype(int) << self.H)
        return state[permutation]

    def toggle_h(self, state, end, *, retest_bank=False):
        predicate = (((self.indices >> self.H) & 1) != 0) & ((self.logical >> end) == 0)
        if retest_bank:
            predicate &= self.bank == 0
        return state[self.indices ^ (predicate.astype(int) << self.h)]

    def target_basis(self, state, target, matrix):
        mask = 1 << (self.bank_width + target)
        low = self.indices[(self.indices & mask) == 0]
        high = low | mask
        result = state.copy()
        result[low] = matrix[0, 0] * state[low] + matrix[0, 1] * state[high]
        result[high] = matrix[1, 0] * state[low] + matrix[1, 1] * state[high]
        return result

    def stage(self, state, target, end, rows):
        basis = np.diag([1, 1j]) @ H
        result = self.target_basis(state, target, basis.conj().T)
        image = self.indices.copy()
        for index in self.indices:
            local = int(self.logical[index])
            if not ((index >> self.h) & 1):
                continue
            if (local >> (target + 1)) & ((1 << (end - target - 1)) - 1):
                continue
            a = rows[target][local & ((1 << target) - 1)]
            if (local >> target) & 1:
                a = -a
            old_bank = int(self.bank[index])
            new_bank = old_bank
            for bit in self.unary:
                new_bank &= ~(1 << bit)
            for j, bit in enumerate(self.unary):
                if (old_bank >> bit) & 1:
                    new_bank |= 1 << self.unary[(j + a) % 4]
            image[index] = int(index) ^ old_bank ^ new_bank
        result = result[np.argsort(image)]
        return self.target_basis(result, target, basis)

    def interior(self, state, groups=(1, 2), rows=None, *, retest_bank=False):
        rows = self.rows if rows is None else rows
        start, result = 0, state
        for height in groups:
            end = start + height
            result = self.toggle_h(result, end, retest_bank=retest_bank)
            for target in range(start, end):
                result = self.stage(result, target, end, rows)
            result = self.toggle_h(result, end, retest_bank=retest_bank)
            start = end
        assert start == 3
        return result

    def apply(self, state, preparation=None, groups=(1, 2), rows=None,
              *, retest_bank=False, ideal_inverse=False):
        preparation = self.native if preparation is None else preparation
        result = self.toggle_H(state)
        result = self.prepare(result, preparation)  # Unconditional actual U.
        result = self.interior(result, groups, rows, retest_bank=retest_bank)
        inverse = self.native.conj().T if ideal_inverse else preparation.conj().T
        result = self.prepare(result, inverse)
        return self.toggle_H(result)  # Apply to the live bank, including leakage.

    def expected(self, initial, rows=None):
        rows = self.rows if rows is None else rows
        result = initial.copy()
        bank_zero = [local << self.bank_width for local in range(8)]
        result[bank_zero] = _ideal_group(rows, 4) @ initial[bank_zero]
        return result

    def perturbed(self, angle):
        # Rotate a unary wire, intentionally producing vacuum/two-hot
        # components. Those rejected components survive the group boundary.
        rotation = np.eye(self.bank_dimension, dtype=complex)
        bit = self.unary[0]
        c, s = np.cos(angle), np.sin(angle)
        for low in range(self.bank_dimension):
            if (low >> bit) & 1:
                continue
            high = low | (1 << bit)
            rotation[np.ix_([low, high], [low, high])] = [[c, -s], [s, c]]
        return rotation @ self.native


class ProtectedUnarySourceTests(unittest.TestCase):
    def test_native_full_bank_preparation_and_initial_zero_predicate(self):
        f = _ProtectedFixture()
        np.testing.assert_allclose(f.native.conj().T @ f.native, np.eye(64),
                                   atol=ATOL, rtol=0)
        ideal = np.zeros(64, dtype=complex)
        for j, bit in enumerate(f.unary):
            ideal[1 << bit] = np.exp(2j * np.pi * j / 4) / 2
        np.testing.assert_allclose(f.native[:, 0], ideal, atol=ATOL, rtol=0)
        initial = f.columns(range(64))
        encoded = f.toggle_H(initial)
        expected = initial.copy()
        expected[0, 0] = 0
        expected[1 << f.H, 0] = 1
        np.testing.assert_array_equal(encoded, expected)
        np.testing.assert_array_equal(f.toggle_H(encoded), initial)

    def test_unequal_groups_full_initial_bank_and_exact_inactive_identity(self):
        f = _ProtectedFixture()
        # Bounded batches retain every logical/bank input with both flags zero.
        for start in range(0, 512, 64):
            initial = f.columns(range(start, start + 64))
            expected = f.expected(initial)
            for groups in ((1, 2), (2, 1)):
                np.testing.assert_allclose(f.apply(initial, groups=groups), expected,
                                           atol=ATOL, rtol=0)
        # The prepared inactive bank need not be zero or within the unary code.
        initial = f.columns(range(512))
        prepared = f.prepare(initial, f.native)
        np.testing.assert_allclose(f.interior(prepared), prepared, atol=ATOL, rtol=0)
        # This is H=0,h=0 identity; no identity for arbitrary h=1 is asserted.

    def test_retained_source_error_includes_final_H_leakage_and_reference(self):
        f = _ProtectedFixture()
        initial = f.columns(local << f.bank_width for local in range(8))
        expected = f.expected(initial)
        for angle in (.013, -.047):
            actual = f.perturbed(angle)
            delta = np.linalg.norm(actual[:, 0] - f.native[:, 0])
            result = f.apply(initial, preparation=actual)
            error = np.linalg.norm(result - expected, 2)
            self.assertLessEqual(error, 2 * delta + ATOL)
            self.assertGreater(error, delta / 2)
            # Group h is exactly erased, while final H contains genuine
            # leakage because the actual source has not returned exactly.
            np.testing.assert_allclose(result[((f.indices >> f.h) & 1) != 0],
                                       0, atol=ATOL, rtol=0)
            self.assertGreater(np.linalg.norm(result[((f.indices >> f.H) & 1) != 0]),
                               abs(angle) / 4)
            # A two-dimensional reference with coherent logical amplitudes.
            reference = np.arange(1, 17).reshape(8, 2).astype(complex)
            reference[:, 1] *= 1j
            reference /= np.linalg.norm(reference)
            self.assertLessEqual(np.linalg.norm((result - expected) @ reference),
                                 2 * delta + ATOL)
            # Nonzero initial bank: unconditional actual U/U† cancel exactly
            # even with this perturbation outside the source code.
            inactive = f.columns([1, 7, 63, 65, 255, 511])
            np.testing.assert_allclose(f.apply(inactive, preparation=actual), inactive,
                                       atol=ATOL, rtol=0)

    def test_physical_bank_retest_and_ideal_inverse_are_detected(self):
        f = _ProtectedFixture()
        initial = f.columns(local << f.bank_width for local in range(8))
        expected = f.expected(initial)
        wrong = f.apply(initial, retest_bank=True)
        # The exact prepared source has no physical all-zero bank component;
        # retesting that bank disables every intended group.
        np.testing.assert_allclose(wrong, initial, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(wrong - expected, 2), 1.9)
        actual = f.perturbed(.037)
        identity_rows = [[0], [0, 0], [0, 0, 0, 0]]
        np.testing.assert_allclose(f.apply(initial, preparation=actual, rows=identity_rows),
                                   initial, atol=ATOL, rtol=0)
        wrong = f.apply(initial, preparation=actual, rows=identity_rows, ideal_inverse=True)
        delta = np.linalg.norm(actual[:, 0] - f.native[:, 0])
        self.assertAlmostEqual(np.linalg.norm(wrong - initial, 2), delta, places=10)
        self.assertGreater(np.linalg.norm(wrong[((f.indices >> f.H) & 1) != 0]), 0.01)


if __name__ == '__main__':
    unittest.main()
