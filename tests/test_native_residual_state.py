"""Bounded native residual preparation, including every dirty-input column.

At q=5 the initialized isometry has 1024 rows and 128 columns.  Native
table sectors retain all target, signal and core inputs; the amplification
never discards a row or resets intermediate leakage.  Independent ideal
coins and Majorana source algebra supply literal-phase oracles.  Fine-q
checks use exact rational promises and gate ledgers, not large simulation.
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.native_residual_state import emit_residual_state
from tests.test_native_residual_rotation import (
    _algebra_rotation,
    _majorana_triplet,
    _native_matrix,
)
from tests.test_native_residual_table import ALPHABET, _inverse, _restrict_sector


ATOL = 4e-10
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
ZERO = (Fraction(0), Fraction(0))
REAL_FIXTURE = ((Fraction(127, 129), Fraction(0)),
                (Fraction(32, 129), Fraction(0)))
COMPLEX_FIXTURE = ((Fraction(511 * 255, 513 * 257), Fraction(511 * 32, 513 * 257)),
                   (Fraction(64 * 3, 513 * 5), Fraction(64 * 4, 513 * 5)))


def _norm_squared(pair):
    return sum((value * value for value in pair), Fraction(0))


def _rounded_fixture(fixture, q):
    """Nearest dyadics with an exact L1 certificate for each complex entry."""
    scale = 1 << (2 * q + 22)
    output = []
    for pair in fixture:
        rounded = tuple(Fraction((2 * value.numerator * scale + value.denominator)
                                 // (2 * value.denominator), scale) for value in pair)
        # L1 controls Euclidean distance without a floating square root.
        assert sum(abs(left - right) for left, right in zip(pair, rounded)) <= Fraction(1, scale)
        assert Fraction(1, scale) < Fraction(1, 1 << (2 * q + 20))
        output.append(rounded)
    return tuple(output)


def _complex(pair):
    return complex(float(pair[0]), float(pair[1]))


def _coin(z):
    complement = np.sqrt(max(0., 1 - abs(z) ** 2))
    return np.array([[z, -complement], [complement, z.conjugate()]])


def _ideal_q(fixture, *, inactive_zero=False):
    """Independent eight-dimensional Q, upper-wire order (s,x,t)."""
    a, w = map(_complex, fixture)
    table = np.zeros((8, 8), dtype=complex)
    controlled_h = np.eye(8, dtype=complex)
    for mode in (0, 1):
        for system in (0, 1):
            start = 4 * mode + 2 * system
            z = a if mode == 0 else (0 if system == 0 else w)
            table[start:start + 2, start:start + 2] = (
                np.eye(2) if inactive_zero and mode == 1 and system == 0 else _coin(z))
    controlled_h[4:, 4:] = np.kron(H, np.eye(2))
    mode_h = np.kron(H, np.eye(4))
    return mode_h @ table @ controlled_h @ mode_h


def _ideal_amplification(q):
    initial = np.eye(8)
    initial[0, 0] = -1
    good = np.eye(8)
    good[0, 0] = good[2, 2] = -1
    return -q @ initial @ q.conj().T @ good @ q


def _embedding(q):
    dirty = 1 << (q + 2)
    return np.vstack((np.eye(dirty, dtype=complex), np.zeros((7 * dirty, dirty), dtype=complex)))


def _target_vector(fixture):
    target = np.zeros(8, dtype=complex)
    target[0], target[2] = _complex(fixture[0]), _complex(fixture[1]) / np.sqrt(2)
    return target


def _target(fixture, q):
    return np.kron(_target_vector(fixture)[:, None], np.eye(1 << (q + 2)))


def _operator_error(columns):
    return np.sqrt(max(0., np.linalg.eigvalsh(columns.conj().T @ columns)[-1]))


def _apply_sectors(sectors, columns, *, inverse=False):
    result = np.empty_like(columns)
    size = len(columns) // 4
    for mode in (0, 1):
        for system in (0, 1):
            start = (2 * mode + system) * size
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


def _apply_amplification(word, sectors, columns):
    width = word.mode + 1
    first = _apply_q(word, sectors, columns)
    result = _native_matrix(word.good_reflection_gates, width, first)
    result = _apply_q(word, sectors, result, inverse=True)
    result = _native_matrix(word.initial_reflection_gates, width, result)
    result = _apply_q(word, sectors, result)
    return first, _native_matrix(word.global_minus_gates, width, result)


def _small_native(gates, wires):
    labels = {wire: index for index, wire in enumerate(wires)}
    remapped = tuple((name, tuple(labels[wire] for wire in operands)) for name, operands in gates)
    return _native_matrix(remapped, len(wires))


class NativeResidualStateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = COMPLEX_FIXTURE
        cls.word = emit_residual_state(*_rounded_fixture(cls.fixture, 5), 5)
        word = cls.word
        cls.initial = _embedding(5)
        cls.native_sectors = {}
        cls.algebra_sectors = {}
        cls.sector_errors = []
        cache = {}
        identity = np.eye(256)
        for table in word.tables:
            for mode in (0, 1):
                for system in (0, 1):
                    reduced, phase_power, width = _restrict_sector(
                        table.gates, word.mode + 1, {word.system: system, word.mode: mode})
                    key = reduced, phase_power
                    if key not in cache:
                        cache[key] = ((1 + 1j) / np.sqrt(2)) ** phase_power * _native_matrix(reduced, width)
                    actual = cache[key]
                    if mode == table.enable_value:
                        zword, yword = (
                            _algebra_rotation(programs[system], _majorana_triplet(programs[system]), axis)
                            for programs, axis in zip(table.programs[:2], ("z", "y")))
                        expected = zword @ yword @ zword
                        cls.native_sectors[mode, system] = actual
                        cls.algebra_sectors[mode, system] = expected
                    else:
                        expected = identity
                    cls.sector_errors.append(np.max(abs(actual - expected)))
        cls.first, cls.actual = _apply_amplification(word, cls.native_sectors, cls.initial)

    def test_exact_normalization_input_promises_and_ideal_amplification(self):
        # The first fixtures satisfy the relaxed 1/(4*sqrt(2)) closeness
        # sufficient for tail loading, not the actual tighter 1/64 coarse
        # guarantee. The endpoint fixtures exercise the broader table API.
        fixtures = (REAL_FIXTURE, COMPLEX_FIXTURE, ((0, 1), ZERO),
                    ((Fraction(1, 2), Fraction(1, 2)), (1, 0)))
        for fixture in fixtures:
            with self.subTest(fixture=fixture):
                a, w = fixture
                self.assertEqual(_norm_squared(a) + _norm_squared(w) / 2, 1)
                self.assertLessEqual(_norm_squared(w), 1)
                for q in (5, 16, 80, 129):
                    _rounded_fixture(fixture, q)
                qideal = _ideal_q(fixture)
                target = _target_vector(fixture)
                # Only clean flags are selected; x=0 and x=1 both remain.
                accepted = np.zeros(8, dtype=complex)
                accepted[[0, 2]] = qideal[[0, 2], 0]
                np.testing.assert_allclose(accepted, target / 2, atol=2e-15, rtol=0)
                np.testing.assert_allclose(_ideal_amplification(qideal)[:, 0], target,
                                           atol=3e-15, rtol=0)
                self.assertGreater(np.linalg.norm(-_ideal_amplification(qideal)[:, 0] - target), 1.99)
                wrong = _ideal_amplification(_ideal_q(fixture, inactive_zero=True))[:, 0]
                self.assertGreater(np.linalg.norm(wrong - target), .1)
        for fixture in fixtures[:2]:
            self.assertLess(2 - 2 * fixture[0][0], Fraction(1, 32))
            # The initial reflection must test x as well as the flags.
            # Replacing it with the good reflection fails already ideally.
            qideal = _ideal_q(fixture)
            good = np.diag([-1, 1, -1, 1, 1, 1, 1, 1])
            wrong_initial = -qideal @ good @ qideal.conj().T @ good @ qideal
            self.assertGreater(np.linalg.norm(wrong_initial[:, 0] - _target_vector(fixture)), .05)

    def test_exact_controlled_h_reflections_and_literal_global_minus(self):
        word = self.word
        actual = _small_native(word.controlled_h_gates, (word.system, word.mode))
        expected = np.eye(4, dtype=complex)
        expected[2:, 2:] = H
        np.testing.assert_allclose(actual, expected, atol=3e-15, rtol=0)
        self.assertEqual(sum(name in ("T", "TDG") for name, _ in word.controlled_h_gates), 2)
        initial = _small_native(word.initial_reflection_gates, (word.target, word.system, word.mode))
        expected = np.eye(8)
        expected[0, 0] = -1
        np.testing.assert_allclose(initial, expected, atol=4e-15, rtol=0)
        self.assertEqual(sum(name in ("T", "TDG") for name, _ in word.initial_reflection_gates), 7)
        good = _small_native(word.good_reflection_gates, (word.target, word.mode))
        np.testing.assert_allclose(good, np.diag([-1, 1, 1, 1]), atol=2e-15, rtol=0)
        minus_wires = sorted({wire for _, wires in word.global_minus_gates for wire in wires})
        np.testing.assert_allclose(_small_native(word.global_minus_gates, minus_wires),
                                   -np.eye(1 << len(minus_wires)), atol=2e-15, rtol=0)

    def test_all_dirty_columns_native_algebra_and_unprojected_amplification(self):
        word = self.word
        self.assertEqual(self.initial.shape, (1024, 128))
        self.assertLess(max(self.sector_errors), ATOL)
        oracle_first, oracle_final = _apply_amplification(word, self.algebra_sectors, self.initial)
        np.testing.assert_allclose(self.first, oracle_first, atol=ATOL, rtol=0)
        np.testing.assert_allclose(self.actual, oracle_final, atol=3 * ATOL, rtol=0)
        np.testing.assert_allclose(self.actual.conj().T @ self.actual, np.eye(128), atol=3 * ATOL, rtol=0)
        target = _target(self.fixture, 5)
        error = _operator_error(self.actual - target)
        self.assertLess(error, float(word.certificate.state_error_bound))
        # The q=5 analytic constant exceeds the maximum isometry error.
        # This fixture additionally has measured error 0.294105, supplying
        # a nonvacuous regression bound without claiming a general rate.
        self.assertLess(error, 1 / 3)
        # The flags carry rejected amplitudes after Q.  Keeping those rows
        # through the literal inverse is essential, even for root inputs.
        dirty = 128
        accepted_rows = np.r_[np.arange(dirty), np.arange(dirty) + 2 * dirty]
        rejected = self.first.copy()
        rejected[accepted_rows] = 0
        self.assertGreater(np.linalg.norm(rejected[:, 0]), .5)
        # Native Q also moves arbitrary work; a stage boundary is not a
        # license to reset its synthesis signal or precision core.
        changed_dirty = np.arange(1024) % dirty != 0
        self.assertGreater(np.linalg.norm(self.first[changed_dirty, 0]), .1)
        reset = self.first - rejected
        reset = _native_matrix(word.good_reflection_gates, 10, reset)
        reset = _apply_q(word, self.native_sectors, reset, inverse=True)
        reset = _native_matrix(word.initial_reflection_gates, 10, reset)
        reset = _apply_q(word, self.native_sectors, reset)
        reset = _native_matrix(word.global_minus_gates, 10, reset)
        self.assertGreater(_operator_error(self.actual - reset), .5)
        # The active zero row is a full U(0) completion, not an identity.
        self.assertGreater(np.linalg.norm((self.native_sectors[1, 0] - np.eye(256))[:, 0]), .5)

    def test_flattened_native_word_on_coherent_dirty_reference_inputs(self):
        word = self.word
        self.assertEqual(word.q_inverse_gates, _inverse(word.q_gates))
        self.assertEqual(word.gates, tuple(gate for stage in word.stages for gate in stage.gates))
        self.assertEqual(tuple(stage.name for stage in word.stages),
                         ("q", "good_reflection", "q_inverse", "initial_reflection", "q", "global_minus"))
        # Two columns are external-reference branches with coherent dirty
        # amplitudes.  Distinct complex phases preclude phase-aligned checks.
        dirty_reference = np.zeros((128, 2), dtype=complex)
        dirty_reference[[0, 1, 65, 127], 0] = np.array([1, 1j, -1, -1j]) / 2
        dirty_reference[[2, 31, 64, 126], 1] = np.array([1j, 1, 1j, -1]) / 2
        inputs = self.initial @ dirty_reference
        actual = _native_matrix(word.gates, 10, inputs)
        np.testing.assert_allclose(actual, self.actual @ dirty_reference, atol=4 * ATOL, rtol=0)
        first = _native_matrix(word.q_gates, 10, inputs)
        np.testing.assert_allclose(_native_matrix(word.q_inverse_gates, 10, first), inputs,
                                   atol=3 * ATOL, rtol=0)

    def test_fine_precision_certificate_layout_and_literal_native_ledger(self):
        for q in (5, 16, 80, 129):
            word = emit_residual_state(*_rounded_fixture(COMPLEX_FIXTURE, q), q)
            with self.subTest(q=q):
                self.assertEqual(word.core, tuple(range(q + 1)))
                self.assertEqual((word.signal, word.target, word.system, word.mode),
                                 (q + 1, q + 2, q + 3, q + 4))
                certificate = word.certificate
                self.assertEqual(certificate.nqubits, q + 5)
                self.assertEqual(certificate.dirty_qubits, q + 2)
                # Disjoint enabled sectors take a maximum, then the three
                # actual Q appearances incur the telescoping factor three.
                self.assertEqual(certificate.q_operator_error_bound,
                                 max(table.operator_error_bound for table in word.tables))
                self.assertEqual(certificate.state_error_bound, 3 * certificate.q_operator_error_bound)
                self.assertLess(certificate.state_error_bound, Fraction(390, 1 << q))
                count = sum(name in ("T", "TDG") for name, _ in word.gates)
                self.assertEqual(count, 3240 * q + 3793)
                self.assertEqual(certificate.t_count, count)
                self.assertEqual(tuple(table.enable_value for table in word.tables), (0, 1))
                self.assertEqual(word.q_inverse_gates, _inverse(word.q_gates))
                for name, wires in word.gates:
                    self.assertIn(name, ALPHABET)
                    self.assertEqual(len(wires), 2 if name == "CX" else 1)
                    self.assertEqual(len(wires), len(set(wires)))
                    self.assertTrue(all(0 <= wire < certificate.nqubits for wire in wires))

    def test_rejects_invalid_precision_coefficients_and_inconsistent_promises(self):
        for q in (True, False, 5., Fraction(5), 4, -1):
            with self.subTest(q=q), self.assertRaises(ValueError):
                emit_residual_state((1, 0), ZERO, q)
        for value in (True, False, .5, complex(1), "1"):
            with self.subTest(value=value), self.assertRaises(TypeError):
                emit_residual_state((value, 0), ZERO, 5)
        for a, w in (((Fraction(1, 3), 0), ZERO), ((1, 0), (Fraction(1, 3), 0)),
                     ((2, 0), ZERO), ((1, 0), (2, 0)), (ZERO, ZERO), ((1, 0), (1, 0))):
            with self.subTest(a=a, w=w), self.assertRaises(ValueError):
                emit_residual_state(a, w, 5)
        for malformed in ((), (1,), (1, 0, 0)):
            with self.subTest(pair=malformed), self.assertRaises(ValueError):
                emit_residual_state(malformed, ZERO, 5)


if __name__ == "__main__":
    unittest.main()
