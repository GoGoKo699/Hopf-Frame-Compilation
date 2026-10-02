"""Two-system-qubit native residual preparation on every dirty input.

The q=5 initialized isometry has 2048 rows and 128 columns. Native table
sectors and independent Majorana algebra retain arbitrary target, signal,
and core inputs. The four-zero reflection is checked on every logical and
borrowed-helper column, including leaked inputs; no intermediate reset or
phase alignment occurs. Fine precision uses exact ledgers, not simulation.
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.native_two_qubit_residual_state import emit_two_qubit_residual_state
from tests.test_native_residual_rotation import _algebra_rotation, _majorana_triplet, _native_matrix
from tests.test_native_residual_state import (
    H, ZERO, _coin, _complex, _norm_squared, _operator_error, _rounded_fixture, _small_native,
)
from tests.test_native_residual_table import ALPHABET, _inverse, _restrict_sector


ATOL = 6e-10
ROOT_RADIUS = Fraction(524281, 524295)
FIXTURE = (
    (ROOT_RADIUS * Fraction(65535, 65537), ROOT_RADIUS * Fraction(512, 65537)),
    (Fraction(2048 * 3, 524295 * 5), Fraction(2048 * 4, 524295 * 5)),
    (Fraction(4096 * 5, 524295 * 13), Fraction(4096 * 12, 524295 * 13)),
    (Fraction(6144 * 8, 524295 * 17), Fraction(-6144 * 15, 524295 * 17)),
)


def _emit_fixture(fixture, q):
    rounded = _rounded_fixture(fixture, q)
    return emit_two_qubit_residual_state(rounded[0], rounded[1:], q)


def _target_vector(fixture):
    target = np.zeros(16, dtype=complex)
    target[0] = _complex(fixture[0])
    target[[2, 4, 6]] = [_complex(tail) / 2 for tail in fixture[1:]]
    return target


def _target(fixture, q):
    return np.kron(_target_vector(fixture)[:, None], np.eye(1 << (q + 2)))


def _ideal_q(fixture, *, inactive_zero=False):
    """Independent logical Q, upper-wire order (s,x_high,x_low,t)."""
    a, *tails = map(_complex, fixture)
    table = np.zeros((16, 16), dtype=complex)
    for mode in (0, 1):
        for system in range(4):
            start = 8 * mode + 2 * system
            z = a if mode == 0 else (0 if system == 0 else tails[system - 1])
            table[start:start + 2, start:start + 2] = (
                np.eye(2) if inactive_zero and mode == 1 and system == 0 else _coin(z))
    controlled_walsh = np.eye(16, dtype=complex)
    controlled_walsh[8:, 8:] = np.kron(np.kron(H, H), np.eye(2))
    mode_h = np.kron(H, np.eye(8))
    return mode_h @ table @ controlled_walsh @ mode_h


def _ideal_amplification(q):
    initial = np.eye(16)
    initial[0, 0] = -1
    good = np.eye(16)
    good[[0, 2, 4, 6], [0, 2, 4, 6]] = -1
    return -q @ initial @ q.conj().T @ good @ q


def _embedding(q):
    dirty = 1 << (q + 2)
    return np.vstack((np.eye(dirty, dtype=complex), np.zeros((15 * dirty, dirty), dtype=complex)))


def _apply_sectors(sectors, columns, *, inverse=False):
    result = np.empty_like(columns)
    size = len(columns) // 8
    for mode in (0, 1):
        for system in range(4):
            start = (4 * mode + system) * size
            matrix = sectors[mode, system]
            result[start:start + size] = (matrix.conj().T if inverse else matrix) @ columns[start:start + size]
    return result


def _apply_q(word, sectors, columns, *, inverse=False):
    width = word.mode + 1
    mode_h = (("H", (word.mode,)),)
    result = _native_matrix(mode_h, width, columns)
    if inverse:
        result = _apply_sectors(sectors, result, inverse=True)
        result = _native_matrix(word.controlled_h_gates, width, result)
    else:
        result = _native_matrix(word.controlled_h_gates, width, result)
        result = _apply_sectors(sectors, result)
    return _native_matrix(mode_h, width, result)


def _finish_amplification(word, sectors, first):
    width = word.mode + 1
    result = _native_matrix(word.good_reflection_gates, width, first)
    result = _apply_q(word, sectors, result, inverse=True)
    result = _native_matrix(word.initial_reflection_gates, width, result)
    result = _apply_q(word, sectors, result)
    return _native_matrix(word.global_minus_gates, width, result)


class NativeTwoQubitResidualStateTests(unittest.TestCase):
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
                for system in range(4):
                    program = rotation.programs[system]
                    if (axis, program) not in algebra_cache:
                        algebra_cache[axis, program] = _algebra_rotation(
                            program, _majorana_triplet(program), axis)
                    for mode in (0, 1):
                        fixed = {word.systems[0]: system & 1, word.systems[1]: system >> 1,
                                 word.mode: mode}
                        reduced, phase_power, width = _restrict_sector(rotation.gates, word.mode + 1, fixed)
                        key = reduced, phase_power
                        if key not in native_cache:
                            native_cache[key] = ((1 + 1j) / np.sqrt(2)) ** phase_power * _native_matrix(reduced, width)
                        actual = native_cache[key]
                        expected = (algebra_cache[axis, program] if mode == table.enable_value else np.eye(256))
                        cls.sector_errors.append(np.max(abs(actual - expected)))
                        factors[factor, system, mode] = actual
            mode = table.enable_value
            for system in range(4):
                zword, yword = factors[0, system, mode], factors[1, system, mode]
                cls.native_sectors[mode, system] = zword @ yword @ zword
                zword = algebra_cache["z", table.programs[0][system]]
                yword = algebra_cache["y", table.programs[1][system]]
                cls.algebra_sectors[mode, system] = zword @ yword @ zword
        cls.first = _apply_q(word, cls.native_sectors, cls.initial)
        cls.actual = _finish_amplification(word, cls.native_sectors, cls.first)
        cls.observed_state_error = _operator_error(cls.actual - _target(FIXTURE, 5))
        cls.observed_dirty_motion = np.linalg.norm(cls.first[np.arange(2048) % 128 != 0, 0])

    def test_exact_fixture_promise_half_amplitude_and_literal_ideal_amplification(self):
        self.assertEqual(_norm_squared(FIXTURE[0]) + sum(_norm_squared(tail) for tail in FIXTURE[1:]) / 4, 1)
        self.assertTrue(all(0 < _norm_squared(tail) < 1 for tail in FIXTURE[1:]))
        distance_squared = 2 - 2 * FIXTURE[0][0]
        self.assertEqual(distance_squared, Fraction(262144, 2290714761))
        self.assertLess(distance_squared, Fraction(1, 64 ** 2))
        fixtures = (FIXTURE, ((0, 1), ZERO, ZERO, ZERO),
                    ((Fraction(1, 2), 0), (1, 0), (0, 1), (-1, 0)))
        for fixture in fixtures:
            with self.subTest(fixture=fixture):
                self.assertEqual(_norm_squared(fixture[0]) + sum(_norm_squared(tail) for tail in fixture[1:]) / 4, 1)
                for q in (5, 16, 80, 129):
                    _rounded_fixture(fixture, q)
                qideal = _ideal_q(fixture)
                target = _target_vector(fixture)
                accepted = np.zeros(16, dtype=complex)
                accepted[[0, 2, 4, 6]] = qideal[[0, 2, 4, 6], 0]
                np.testing.assert_allclose(accepted, target / 2, atol=3e-15, rtol=0)
                actual = _ideal_amplification(qideal)[:, 0]
                np.testing.assert_allclose(actual, target, atol=4e-15, rtol=0)
                self.assertGreater(np.linalg.norm(-actual - target), 1.99)
                wrong_zero = _ideal_amplification(_ideal_q(fixture, inactive_zero=True))[:, 0]
                self.assertGreater(np.linalg.norm(wrong_zero - target), .1)
        # The boundary fixture is a valid bounded-table input, although
        # outside the coarse neighborhood. Permute its distinct phases to
        # make omission of either system-control bit separately visible.
        boundary = fixtures[-1]
        self.assertLess(_emit_fixture(boundary, 5).state_error_bound, Fraction(390, 1 << 5))
        witness = (boundary[0], boundary[2], boundary[3], boundary[1])
        qideal = _ideal_q(witness)
        good = np.eye(16)
        good[[0, 2, 4, 6], [0, 2, 4, 6]] = -1
        for missing_system in (1, 2):
            initial = np.eye(16)
            initial[0, 0] = initial[1 << missing_system, 1 << missing_system] = -1
            wrong_initial = -qideal @ initial @ qideal.conj().T @ good @ qideal
            self.assertGreater(np.linalg.norm(wrong_initial[:, 0] - _target_vector(witness)), .1)

    def test_exact_four_zero_reflection_returns_arbitrary_borrowed_helper(self):
        word = self.word
        helper = word.reflection_helper
        self.assertEqual(helper, word.core[0])
        wires = (helper, word.target, *word.systems, word.mode)
        actual = _small_native(word.initial_reflection_gates, wires)
        expected = np.eye(32)
        # Helper occupies the low bit; both arbitrary helper values acquire
        # minus exactly when the four logical inputs are all zero.
        expected[0, 0] = expected[1, 1] = -1
        np.testing.assert_allclose(actual, expected, atol=8e-15, rtol=0)
        self.assertEqual(sum(name in ("T", "TDG") for name, _ in word.initial_reflection_gates), 28)
        self.assertEqual({wire for _, operands in word.initial_reflection_gates for wire in operands}, set(wires))
        # This arbitrary entangled input includes nonzero logical flags and
        # both helper values. The full matrix identity already covers every
        # external reference, without a zero-helper or no-leakage premise.
        entangled = np.zeros((32, 2), dtype=complex)
        entangled[[0, 3, 13, 30], 0] = np.array([1, 1j, -1, -1j]) / 2
        entangled[[1, 6, 20, 31], 1] = np.array([1j, -1, 1, 1j]) / 2
        np.testing.assert_allclose(actual @ entangled, expected @ entangled, atol=8e-15, rtol=0)
        controlled_h = _small_native(word.controlled_h_gates, (*word.systems, word.mode))
        ideal = np.eye(8, dtype=complex)
        ideal[4:, 4:] = np.kron(H, H)
        np.testing.assert_allclose(controlled_h, ideal, atol=4e-15, rtol=0)
        self.assertEqual(sum(name in ("T", "TDG") for name, _ in word.controlled_h_gates), 4)
        good = _small_native(word.good_reflection_gates, (word.target, word.mode))
        np.testing.assert_allclose(good, np.diag([-1, 1, 1, 1]), atol=3e-15, rtol=0)
        minus_wires = sorted({wire for _, operands in word.global_minus_gates for wire in operands})
        np.testing.assert_allclose(_small_native(word.global_minus_gates, minus_wires),
                                   -np.eye(1 << len(minus_wires)), atol=3e-15, rtol=0)

    def test_complete_native_dirty_isometry_and_retained_intermediate_leakage(self):
        word = self.word
        self.assertEqual(self.initial.shape, (2048, 128))
        self.assertLess(max(self.sector_errors), ATOL)
        algebra_first = _apply_q(word, self.algebra_sectors, self.initial)
        algebra_final = _finish_amplification(word, self.algebra_sectors, algebra_first)
        np.testing.assert_allclose(self.first, algebra_first, atol=ATOL, rtol=0)
        np.testing.assert_allclose(self.actual, algebra_final, atol=4 * ATOL, rtol=0)
        np.testing.assert_allclose(self.actual.conj().T @ self.actual, np.eye(128), atol=4 * ATOL, rtol=0)
        self.assertLess(self.observed_state_error, float(word.certificate.state_error_bound))
        # The observed error is 0.360099. This finite-fixture regression
        # is nonvacuous despite the loose general analytic q=5 bound.
        self.assertLess(self.observed_state_error, .4)
        self.assertGreater(self.observed_dirty_motion, .05)
        accepted_rows = np.concatenate([np.arange(128) + 256 * system for system in range(4)])
        rejected = self.first.copy()
        rejected[accepted_rows] = 0
        self.assertGreater(np.linalg.norm(rejected[:, 0]), .5)
        projected = _finish_amplification(word, self.native_sectors, self.first - rejected)
        self.assertGreater(_operator_error(self.actual - projected), .5)
        self.assertGreater(np.linalg.norm((self.native_sectors[1, 0] - np.eye(256))[:, 0]), .5)

    def test_flattened_native_word_and_actual_inverse_on_coherent_reference_inputs(self):
        word = self.word
        self.assertEqual(word.q_inverse_gates, _inverse(word.q_gates))
        self.assertEqual(word.gates, tuple(gate for stage in word.stages for gate in stage.gates))
        self.assertEqual(tuple(stage.name for stage in word.stages),
                         ("q", "good_reflection", "q_inverse", "initial_reflection", "q", "global_minus"))
        reference = np.zeros((128, 2), dtype=complex)
        reference[[0, 1, 65, 126], 0] = np.array([1, 1j, -1, -1j]) / 2
        reference[[2, 31, 64, 127], 1] = np.array([1j, -1, 1, 1j]) / 2
        inputs = self.initial @ reference
        actual = _native_matrix(word.gates, 11, inputs)
        np.testing.assert_allclose(actual, self.actual @ reference, atol=5 * ATOL, rtol=0)
        first = _native_matrix(word.q_gates, 11, inputs)
        np.testing.assert_allclose(_native_matrix(word.q_inverse_gates, 11, first), inputs,
                                   atol=4 * ATOL, rtol=0)

    def test_fine_precision_certificates_layout_and_literal_native_t_counts(self):
        for q in (5, 16, 80, 129):
            with self.subTest(q=q):
                word = _emit_fixture(FIXTURE, q)
                certificate = word.certificate
                self.assertEqual(word.core, tuple(range(q + 1)))
                self.assertEqual((word.signal, word.target, word.systems, word.mode),
                                 (q + 1, q + 2, (q + 3, q + 4), q + 5))
                self.assertEqual(word.reflection_helper, word.core[0])
                self.assertEqual((certificate.nqubits, certificate.dirty_qubits), (q + 6, q + 2))
                delta = Fraction(1, 1 << (2 * q + 20))
                self.assertEqual(certificate.input_error_bound, delta)
                self.assertEqual(certificate.normalization_tolerance, Fraction(7, 2) * delta + Fraction(7, 4) * delta ** 2)
                rounded = _rounded_fixture(FIXTURE, q)
                discrepancy = abs(_norm_squared(rounded[0]) + sum(_norm_squared(tail) for tail in rounded[1:]) / 4 - 1)
                self.assertEqual(certificate.normalization_discrepancy, discrepancy)
                self.assertLessEqual(discrepancy, certificate.normalization_tolerance)
                self.assertEqual(certificate.q_operator_error_bound,
                                 max(table.operator_error_bound for table in word.tables))
                self.assertEqual(certificate.state_error_bound, 3 * certificate.q_operator_error_bound)
                self.assertLess(certificate.state_error_bound, Fraction(390, 1 << q))
                self.assertEqual(tuple(table.enable_value for table in word.tables), (0, 1))
                self.assertTrue(all(rotation.quadratic_axis_count == 0 for rotation in word.tables[0].rotations))
                kz, ky = (rotation.quadratic_axis_count for rotation in word.tables[1].rotations[:2])
                expected_t = 3240 * q + 3820 + 210 * (2 * kz + ky)
                self.assertEqual(certificate.t_count, expected_t)
                self.assertLessEqual(certificate.t_count, 3240 * q + 5080)
                self.assertEqual(sum(name in ("T", "TDG") for name, _ in word.gates), expected_t)
                self.assertEqual(word.q_inverse_gates, _inverse(word.q_gates))
                for name, wires in word.gates:
                    self.assertIn(name, ALPHABET)
                    self.assertEqual(len(wires), 2 if name == "CX" else 1)
                    self.assertEqual(len(wires), len(set(wires)))
                    self.assertTrue(all(0 <= wire < certificate.nqubits for wire in wires))

    def test_rejects_invalid_precision_rows_coefficients_and_inconsistent_promises(self):
        zero_tails = (ZERO,) * 3
        for q in (True, False, 5., Fraction(5), 4, -1):
            with self.subTest(q=q), self.assertRaises(ValueError):
                emit_two_qubit_residual_state((1, 0), zero_tails, q)
        for value in (True, False, .5, complex(1), "1"):
            with self.subTest(value=value), self.assertRaises(TypeError):
                emit_two_qubit_residual_state((value, 0), zero_tails, 5)
            with self.subTest(tail=value), self.assertRaises(TypeError):
                emit_two_qubit_residual_state((1, 0), ((value, 0), ZERO, ZERO), 5)
        for a, tails in (((Fraction(1, 3), 0), zero_tails),
                         ((1, 0), ((Fraction(1, 3), 0), ZERO, ZERO)),
                         ((2, 0), zero_tails), (ZERO, zero_tails),
                         ((1, 0), ((1, 0), ZERO, ZERO)),
                         (ZERO, ((2, 0), ZERO, ZERO))):
            with self.subTest(a=a, tails=tails), self.assertRaises(ValueError):
                emit_two_qubit_residual_state(a, tails, 5)
        for tails in ((), (ZERO,), (ZERO,) * 2, (ZERO,) * 4):
            with self.subTest(tails=tails), self.assertRaises(ValueError):
                emit_two_qubit_residual_state((1, 0), tails, 5)
        for malformed in ((), (1,), (1, 0, 0)):
            with self.subTest(pair=malformed), self.assertRaises(ValueError):
                emit_two_qubit_residual_state(malformed, zero_tails, 5)


if __name__ == "__main__":
    unittest.main()
