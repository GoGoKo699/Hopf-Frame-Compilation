"""Coherent branch-selected residual preparation on every dirty input.

At q=5 the initialized isometry has 2048 rows and 256 branch/dirty input
columns. Complete native sectors and independent Majorana algebra retain
all rejected amplitudes and literal branch phases. External-reference
columns also traverse the flattened native word and its actual inverse.
Fine-precision checks use exact promises and ledgers, not large simulation.
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.native_branched_residual_state import emit_branched_residual_state
from tests.test_native_residual_rotation import _algebra_rotation, _majorana_triplet, _native_matrix
from tests.test_native_residual_state import (
    H, ZERO, _complex, _ideal_q as _one_branch_q, _norm_squared,
    _operator_error, _rounded_fixture, _small_native,
)
from tests.test_native_residual_table import ALPHABET, _inverse, _restrict_sector
from tests.test_native_two_qubit_residual_state import _apply_q, _finish_amplification


ATOL = 6e-10
FIXTURE = (
    ((Fraction(32767 * 65535, 32769 * 65537), Fraction(32767 * 512, 32769 * 65537)),
     (Fraction(512 * 3, 32769 * 5), Fraction(512 * 4, 32769 * 5))),
    ((Fraction(131071 * 262143, 131073 * 262145), Fraction(-131071 * 1024, 131073 * 262145)),
     (Fraction(1024 * 5, 131073 * 13), Fraction(-1024 * 12, 131073 * 13))),
)
ROOT_REFERENCE = (((1, 0), ZERO), FIXTURE[1])
PHASE_WITNESS = (((1, 0), ZERO), ((0, 1), ZERO))


def _rounded_branches(branches, q):
    return tuple(_rounded_fixture(branch, q) for branch in branches)


def _emit_fixture(branches, q):
    return emit_branched_residual_state(_rounded_branches(branches, q), q)


def _embedding(q):
    dirty = 1 << (q + 2)
    result = np.zeros((16 * dirty, 2 * dirty), dtype=complex)
    for branch in (0, 1):
        result[4 * branch * dirty:4 * branch * dirty + dirty,
               branch * dirty:(branch + 1) * dirty] = np.eye(dirty)
    return result


def _target_logical(branches):
    result = np.zeros((16, 2), dtype=complex)
    for branch, (a, w) in enumerate(branches):
        result[4 * branch, branch] = _complex(a)
        result[4 * branch + 2, branch] = _complex(w) / np.sqrt(2)
    return result


def _target(branches, q):
    return np.kron(_target_logical(branches), np.eye(1 << (q + 2)))


def _ideal_q(branches):
    """Embed independently constructed one-branch Q blocks literally."""
    result = np.zeros((16, 16), dtype=complex)
    for branch, fixture in enumerate(branches):
        rows = np.r_[np.arange(4), np.arange(4) + 8] + 4 * branch
        result[np.ix_(rows, rows)] = _one_branch_q(fixture)
    return result


def _ideal_amplification(q, *, include_branch_in_initial=False):
    initial = np.eye(16)
    initial[0, 0] = -1
    if not include_branch_in_initial:
        initial[4, 4] = -1
    good = np.eye(16)
    good[[0, 2, 4, 6], [0, 2, 4, 6]] = -1
    return -q @ initial @ q.conj().T @ good @ q


class NativeBranchedResidualStateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.word = _emit_fixture(FIXTURE, 5)
        word = cls.word
        cls.initial = _embedding(5)
        cls.native_sectors, cls.algebra_sectors = {}, {}
        cls.sector_errors = []
        native_cache, algebra_cache = {}, {}
        for table in word.tables:
            factors = {}
            for factor, rotation in enumerate(table.rotations[:2]):
                axis = "z" if factor == 0 else "y"
                for row in range(4):
                    program = rotation.programs[row]
                    if (axis, program) not in algebra_cache:
                        algebra_cache[axis, program] = _algebra_rotation(
                            program, _majorana_triplet(program), axis)
                    for mode in (0, 1):
                        fixed = {word.system: row & 1, word.branch: row >> 1, word.mode: mode}
                        reduced, phase_power, width = _restrict_sector(rotation.gates, word.mode + 1, fixed)
                        key = reduced, phase_power
                        if key not in native_cache:
                            native_cache[key] = ((1 + 1j) / np.sqrt(2)) ** phase_power * _native_matrix(reduced, width)
                        actual = native_cache[key]
                        expected = algebra_cache[axis, program] if mode == table.enable_value else np.eye(256)
                        cls.sector_errors.append(np.max(abs(actual - expected)))
                        factors[factor, row, mode] = actual
            mode = table.enable_value
            for row in range(4):
                zword, yword = factors[0, row, mode], factors[1, row, mode]
                cls.native_sectors[mode, row] = zword @ yword @ zword
                zword = algebra_cache["z", table.programs[0][row]]
                yword = algebra_cache["y", table.programs[1][row]]
                cls.algebra_sectors[mode, row] = zword @ yword @ zword
        cls.first = _apply_q(word, cls.native_sectors, cls.initial)
        cls.actual = _finish_amplification(word, cls.native_sectors, cls.first)
        cls.observed_state_error = _operator_error(cls.actual - _target(FIXTURE, 5))
        cls.observed_dirty_motion = np.linalg.norm(cls.first[np.arange(2048) % 128 != 0, 0])

    def test_exact_branch_promises_root_reference_and_relative_phase_witness(self):
        expected_distances = (Fraction(131072, 715860651), Fraction(524288, 11453377195))
        for branch, expected_distance in zip(FIXTURE, expected_distances):
            a, w = branch
            self.assertEqual(_norm_squared(a) + _norm_squared(w) / 2, 1)
            self.assertTrue(0 < _norm_squared(w) < 1)
            self.assertEqual(2 - 2 * a[0], expected_distance)
            self.assertLess(expected_distance, Fraction(1, 64 ** 2))
        for branches in (FIXTURE, ROOT_REFERENCE, PHASE_WITNESS):
            with self.subTest(branches=branches):
                for a, w in branches:
                    self.assertEqual(_norm_squared(a) + _norm_squared(w) / 2, 1)
                qideal = _ideal_q(branches)
                initial = np.eye(16)[:, [0, 4]]
                target = _target_logical(branches)
                accepted = qideal @ initial
                accepted[[1, 3, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15]] = 0
                np.testing.assert_allclose(accepted, target / 2, atol=3e-15, rtol=0)
                actual = _ideal_amplification(qideal) @ initial
                np.testing.assert_allclose(actual, target, atol=4e-15, rtol=0)
                wrong = _ideal_amplification(qideal, include_branch_in_initial=True) @ initial
                self.assertAlmostEqual(_operator_error(wrong - target), 1, delta=3e-15)
                self.assertAlmostEqual(np.linalg.norm((wrong - target) @ np.ones(2) / np.sqrt(2)),
                                       1 / np.sqrt(2), delta=3e-15)
                self.assertGreater(_operator_error(-actual - target), 1.99)
        # This exact phase-only pair is outside the small coarse radius but
        # in the bounded API. Phase-aligning the two columns independently
        # would erase the branch-Y signal and produce a different state.
        plus_input = np.eye(16)[:, [0, 4]] @ np.ones(2) / np.sqrt(2)
        prepared = _ideal_amplification(_ideal_q(PHASE_WITNESS)) @ plus_input
        branch_y = np.kron(np.eye(2), np.kron(np.array([[0, -1j], [1j, 0]]), np.eye(4)))
        branch_x = np.kron(np.eye(2), np.kron(np.array([[0, 1], [1, 0]]), np.eye(4)))
        self.assertAlmostEqual(np.vdot(prepared, branch_y @ prepared).real, 1, delta=3e-15)
        self.assertAlmostEqual(np.vdot(prepared, branch_x @ prepared).real, 0, delta=3e-15)
        aligned = prepared.copy()
        aligned[[4, 5, 6, 7, 12, 13, 14, 15]] *= -1j
        self.assertAlmostEqual(np.linalg.norm(prepared - aligned), 1, delta=3e-15)
        self.assertAlmostEqual(np.vdot(aligned, branch_y @ aligned).real, 0, delta=3e-15)

    def test_exact_reflection_excludes_branch_on_all_logical_inputs(self):
        word = self.word
        wires = (word.target, word.system, word.branch, word.mode)
        actual = _small_native(word.initial_reflection_gates, wires)
        expected = np.eye(16)
        expected[0, 0] = expected[4, 4] = -1
        np.testing.assert_allclose(actual, expected, atol=5e-15, rtol=0)
        self.assertNotIn(word.branch, {wire for _, operands in word.initial_reflection_gates for wire in operands})
        self.assertEqual(sum(name in ("T", "TDG") for name, _ in word.initial_reflection_gates), 7)
        controlled_h = _small_native(word.controlled_h_gates, (word.system, word.mode))
        ideal = np.eye(4, dtype=complex)
        ideal[2:, 2:] = H
        np.testing.assert_allclose(controlled_h, ideal, atol=3e-15, rtol=0)
        self.assertEqual(sum(name in ("T", "TDG") for name, _ in word.controlled_h_gates), 2)
        good = _small_native(word.good_reflection_gates, (word.target, word.mode))
        np.testing.assert_allclose(good, np.diag([-1, 1, 1, 1]), atol=3e-15, rtol=0)
        minus_wires = sorted({wire for _, operands in word.global_minus_gates for wire in operands})
        np.testing.assert_allclose(_small_native(word.global_minus_gates, minus_wires),
                                   -np.eye(1 << len(minus_wires)), atol=3e-15, rtol=0)

    def test_all_branch_dirty_columns_native_oracle_and_retained_leakage(self):
        word = self.word
        self.assertEqual(self.initial.shape, (2048, 256))
        self.assertLess(max(self.sector_errors), ATOL)
        algebra_first = _apply_q(word, self.algebra_sectors, self.initial)
        algebra_final = _finish_amplification(word, self.algebra_sectors, algebra_first)
        np.testing.assert_allclose(self.first, algebra_first, atol=ATOL, rtol=0)
        np.testing.assert_allclose(self.actual, algebra_final, atol=4 * ATOL, rtol=0)
        np.testing.assert_allclose(self.actual.conj().T @ self.actual, np.eye(256), atol=4 * ATOL, rtol=0)
        self.assertLess(self.observed_state_error, float(word.certificate.state_error_bound))
        # The general q=5 certificate is loose; this numerical fixture also
        # has a nonvacuous full branch/dirty isometry regression bound.
        self.assertLess(self.observed_state_error, .5)
        self.assertGreater(self.observed_dirty_motion, .05)
        for branch in (0, 1):
            wrong_branch_rows = ((np.arange(2048) >> word.branch) & 1) != branch
            np.testing.assert_allclose(self.actual[wrong_branch_rows, 128 * branch:128 * (branch + 1)],
                                       0, atol=4 * ATOL, rtol=0)
        accepted_rows = np.concatenate([np.arange(128) + 256 * row for row in range(4)])
        rejected = self.first.copy()
        rejected[accepted_rows] = 0
        self.assertGreater(np.linalg.norm(rejected[:, 0]), .5)
        projected = _finish_amplification(word, self.native_sectors, self.first - rejected)
        self.assertGreater(_operator_error(self.actual - projected), .5)
        for row in (0, 2):
            self.assertGreater(np.linalg.norm((self.native_sectors[1, row] - np.eye(256))[:, 0]), .5)

    def test_flattened_native_coherent_branch_reference_columns_and_actual_inverse(self):
        word = self.word
        self.assertEqual(word.q_inverse_gates, _inverse(word.q_gates))
        self.assertEqual(word.gates, tuple(gate for stage in word.stages for gate in stage.gates))
        self.assertEqual(tuple(stage.name for stage in word.stages),
                         ("q", "good_reflection", "q_inverse", "initial_reflection", "q", "global_minus"))
        # Each external-reference column coherently occupies both protocol
        # branches and both signal values, with different dirty correlations.
        reference = np.zeros((256, 2), dtype=complex)
        reference[[0, 1, 65, 126, 128, 131, 192, 255], 0] = np.array([1, 1j, -1, -1j, 1j, 1, -1j, -1]) / np.sqrt(8)
        reference[[2, 31, 64, 127, 129, 160, 193, 254], 1] = np.array([1j, -1, 1, 1j, -1j, 1, -1, 1j]) / np.sqrt(8)
        inputs = self.initial @ reference
        actual = _native_matrix(word.gates, 11, inputs)
        np.testing.assert_allclose(actual, self.actual @ reference, atol=5 * ATOL, rtol=0)
        first = _native_matrix(word.q_gates, 11, inputs)
        np.testing.assert_allclose(_native_matrix(word.q_inverse_gates, 11, first), inputs,
                                   atol=4 * ATOL, rtol=0)

    def test_fine_precision_branch_certificates_layout_and_t_ledger(self):
        for q in (5, 16, 80, 129):
            with self.subTest(q=q):
                word = _emit_fixture(FIXTURE, q)
                certificate = word.certificate
                self.assertEqual(word.core, tuple(range(q + 1)))
                self.assertEqual((word.signal, word.target, word.system, word.branch, word.mode),
                                 (q + 1, q + 2, q + 3, q + 4, q + 5))
                self.assertEqual((certificate.nqubits, certificate.dirty_qubits), (q + 6, q + 2))
                delta = Fraction(1, 1 << (2 * q + 20))
                self.assertEqual(certificate.input_error_bound, delta)
                self.assertEqual(certificate.normalization_tolerance, 3 * delta + Fraction(3, 2) * delta ** 2)
                rounded = _rounded_branches(FIXTURE, q)
                discrepancies = tuple(abs(_norm_squared(a) + _norm_squared(w) / 2 - 1) for a, w in rounded)
                self.assertEqual(certificate.normalization_discrepancies, discrepancies)
                self.assertTrue(all(value <= certificate.normalization_tolerance for value in discrepancies))
                self.assertEqual(certificate.q_operator_error_bound,
                                 max(table.operator_error_bound for table in word.tables))
                self.assertEqual(certificate.state_error_bound, 3 * certificate.q_operator_error_bound)
                self.assertLess(certificate.state_error_bound, Fraction(390, 1 << q))
                self.assertEqual(tuple(table.enable_value for table in word.tables), (0, 1))
                self.assertTrue(all(rotation.quadratic_axis_count == 0 for rotation in word.tables[0].rotations))
                kz, ky = (rotation.quadratic_axis_count for rotation in word.tables[1].rotations[:2])
                expected_t = 3240 * q + 3793 + 210 * (2 * kz + ky)
                self.assertEqual(certificate.t_count, expected_t)
                self.assertLessEqual(certificate.t_count, 3240 * q + 5053)
                self.assertEqual(sum(name in ("T", "TDG") for name, _ in word.gates), expected_t)
                self.assertEqual(word.q_inverse_gates, _inverse(word.q_gates))
                for name, wires in word.gates:
                    self.assertIn(name, ALPHABET)
                    self.assertEqual(len(wires), 2 if name == "CX" else 1)
                    self.assertEqual(len(wires), len(set(wires)))
                    self.assertTrue(all(0 <= wire < certificate.nqubits for wire in wires))
        # Root-reference specialization is the common-coarse residual pair;
        # its root branch consumes no additional initialized compiler qubit.
        reference_word = _emit_fixture(ROOT_REFERENCE, 16)
        self.assertEqual(reference_word.coefficients[0], ((Fraction(1), Fraction(0)), ZERO))
        self.assertEqual((reference_word.nqubits, reference_word.dirty_qubits), (22, 18))
        self.assertLess(reference_word.state_error_bound, Fraction(390, 1 << 16))
        phase_word = _emit_fixture(PHASE_WITNESS, 5)
        self.assertEqual(phase_word.t_count, 3240 * 5 + 3793)

    def test_rejects_invalid_precision_nested_rows_and_branchwise_promises(self):
        valid_branch = ((1, 0), ZERO)
        for q in (True, False, 5., Fraction(5), 4, -1):
            with self.subTest(q=q), self.assertRaises(ValueError):
                emit_branched_residual_state((valid_branch,) * 2, q)
        for branches in ((), (valid_branch,), (valid_branch,) * 3,
                         ((valid_branch[0],), valid_branch), ((valid_branch[0],) * 3, valid_branch)):
            with self.subTest(branches=branches), self.assertRaises(ValueError):
                emit_branched_residual_state(branches, 5)
        for value in (True, False, .5, complex(1), "1"):
            with self.subTest(value=value), self.assertRaises(TypeError):
                emit_branched_residual_state((((value, 0), ZERO), valid_branch), 5)
        for invalid in (((Fraction(1, 3), 0), ZERO), ((1, 0), (Fraction(1, 3), 0)),
                        ((2, 0), ZERO), (ZERO, ZERO), ((1, 0), (1, 0))):
            for branch in (0, 1):
                branches = [valid_branch, valid_branch]
                branches[branch] = invalid
                with self.subTest(branch=branch, invalid=invalid), self.assertRaises(ValueError):
                    emit_branched_residual_state(tuple(branches), 5)
        # The average norm equals one, but neither branch is normalized.
        # A check on the two branches' sum would incorrectly accept this.
        with self.assertRaises(ValueError):
            emit_branched_residual_state((((1, 0), (1, 0)), (ZERO, (1, 0))), 5)


if __name__ == "__main__":
    unittest.main()
