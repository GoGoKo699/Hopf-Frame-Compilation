"""Fine-residual QBP integration on selected coherent dirty/reference inputs.

Two 2048-row columns traverse complete q=5 native magnitude and phase
words. Independent rational geometry, Majorana algebra, and ideal gradient
means provide separate oracles. All work/flag outcomes contribute to scores;
there is no postselection, intermediate reset, or per-branch phase fit.
This bounded integration does not repeat the primitive all-input audits.
"""
from fractions import Fraction
from math import isqrt
import unittest

import numpy as np

from compiler_robust_hopf import native_residual_qbp as native
from compiler_robust_hopf.complex_coarse_decoder import decode_complex_coarse_frame_histograms
from tests.test_native_residual_rotation import _algebra_rotation, _majorana_triplet, _native_matrix
from tests.test_native_residual_state import H, _apply_amplification, _operator_error, _small_native
from tests.test_native_residual_table import ALPHABET, _inverse
from tests.test_native_two_qubit_residual_state import _apply_q, _finish_amplification


ATOL = 9e-10
ANGLE_DENOMINATOR, PHASE_DENOMINATOR = 1024, 256


def _mul(left, right):
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c


def _matvec(matrix, vector):
    return tuple(tuple(sum((_mul(entry, value)[part] for entry, value in zip(row, vector)), Fraction(0))
                       for part in (0, 1)) for row in matrix)


def _complex(pair):
    return complex(float(pair[0]), float(pair[1]))


def _norm_squared(pair):
    return sum((value * value for value in pair), Fraction(0))


def _oracle():
    # Rational unit-circle parametrization independently fixes the angles.
    t, p = Fraction(1, ANGLE_DENOMINATOR), Fraction(1, PHASE_DENOMINATOR)
    cu, su = (1 - t * t) / (1 + t * t), 2 * t / (1 + t * t)
    cv, sv = (1 - p * p) / (1 + p * p), 2 * p / (1 + p * p)
    a, b = (cu * cv, su * sv), (su * cv, cu * sv)
    half = Fraction(1, 2)
    c_exact = (((half, -half), (-half, half)), ((half, half), (half, half)))
    state = _matvec(c_exact, (a, b))
    derivative = _matvec(c_exact, ((-su * cv, cu * sv), (cu * cv, -su * sv)))
    coarse = np.array([[_complex(value) for value in row] for row in c_exact])
    psi, dtheta = np.array(list(map(_complex, state))), np.array(list(map(_complex, derivative)))
    magnitude = np.sqrt(2) * float(-(cu * cu - su * su) + 4 * cu * su * cv * sv)
    phase = float((cu * cu - su * su) * (cv * cv - sv * sv)) / np.sqrt(2)
    wx = (4 * (a[0] - b[0]), -4 * (a[0] + b[0]))
    wy = (4 * (b[1] - a[1]), 4 * (b[1] + a[1]))
    return dict(cu=cu, su=su, cv=cv, sv=sv, a=a, b=b, coarse_exact=c_exact,
                state_exact=state, derivative_exact=derivative, C=coarse, psi=psi,
                derivative=dtheta, magnitude=magnitude, phase=np.array([phase, -phase]),
                wx=wx, wy=wy)


def _logical_output(oracle, stream, *, inverse="actual", branch_phase=1):
    """Independent four-dimensional branch/system protocol, O fixed to H."""
    c, psi = oracle["C"], oracle["psi"]
    reference = c[:, 0] if stream != "phase" else psi
    state = np.concatenate((reference, branch_phase * psi)) / np.sqrt(2)
    observable = np.eye(4, dtype=complex)
    observable[2:, 2:] = H
    state = observable @ state
    if stream != "phase":
        chosen = c.conj().T if inverse == "actual" else c if inverse == "forward" else np.eye(2)
        state = np.kron(np.eye(2), H @ chosen) @ state
    elif inverse == "extra":
        state = np.kron(np.eye(2), c.conj().T) @ state
    readout = H if stream == "X" else H @ np.diag([1, -1j])
    return np.kron(readout, np.eye(2)) @ state


def _logical_score(state, weights):
    return np.dot(abs(state) ** 2, np.array([1, 1, -1, -1]) * np.tile(weights, 2)).real


def _physical_columns(logical, dirty):
    # Upper logical bits are (mode,branch,system,target); flags are zero.
    clean = np.zeros(16, dtype=complex)
    clean[[0, 2, 4, 6]] = logical
    return np.kron(clean[:, None], dirty)


def _score_operator(bundle, columns, weights):
    indices = np.arange(len(columns))
    signs = 1 - 2 * ((indices >> bundle.branch) & 1)
    leaf = (indices >> bundle.system) & 1
    scores = signs * np.asarray(weights, dtype=float)[leaf]
    return columns.conj().T @ (scores[:, None] * columns)


def _algebra_sectors(word, cache):
    sectors = {}
    for table in word.tables:
        for row in range(len(table.rotation_data)):
            matrices = []
            for axis, programs in zip(("z", "y"), table.programs[:2]):
                program = programs[row]
                if (axis, program) not in cache:
                    cache[axis, program] = _algebra_rotation(program, _majorana_triplet(program), axis)
                matrices.append(cache[axis, program])
            sectors[table.enable_value, row] = matrices[0] @ matrices[1] @ matrices[0]
    return sectors


class NativeResidualQBPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.oracle = _oracle()
        cls.bundle = native.emit_residual_qbp(5)
        bundle = cls.bundle
        cls.dirty = np.zeros((128, 2), dtype=complex)
        cls.dirty[[0, 1, 65, 126], 0] = np.array([1, 1j, -1, -1j]) / 2
        cls.dirty[[2, 31, 64, 127], 1] = np.array([1j, -1, 1, 1j]) / 2
        initial = np.zeros((2048, 2), dtype=complex)
        initial[:128] = cls.dirty
        cls.initial = initial
        branch_h = (("H", (bundle.branch,)),)
        cls.outputs = {
            "X": _native_matrix(bundle.magnitude_x.gates, 11, initial),
            "phase": _native_matrix(bundle.phase.gates, 11, initial),
        }
        # Undo the X readout and apply the actual Y basis change. A literal
        # word identity below verifies that the entire preceding word agrees.
        cls.outputs["Y"] = _native_matrix(
            branch_h + (("SDG", (bundle.branch,)),) + branch_h, 11, cls.outputs["X"])
        cache = {}
        pair = bundle.pair_residual
        pair_sectors = _algebra_sectors(pair, cache)
        pair_first = _apply_q(pair, pair_sectors, _native_matrix(branch_h, 11, initial))
        pair_output = _finish_amplification(pair, pair_sectors, pair_first)
        pair_prepared = _native_matrix(bundle.coarse_gates, 11, pair_output)
        cls.algebra_outputs = {}
        for key, record in (("X", bundle.magnitude_x), ("Y", bundle.magnitude_y)):
            suffix = record.gates[len(bundle.pair_preparation_gates):]
            cls.algebra_outputs[key] = _native_matrix(suffix, 11, pair_prepared)
        plain = bundle.phase_residual
        plain_initial = np.zeros((1024, 2), dtype=complex)
        plain_initial[:128] = cls.dirty
        _, plain_output = _apply_amplification(plain, _algebra_sectors(plain, cache), plain_initial)
        remapped = np.zeros((2048, 2), dtype=complex)
        old_rows = np.arange(1024)
        new_rows = (old_rows & 511) | ((old_rows >> 9) << 10)
        remapped[new_rows] = plain_output
        phase_suffix = bundle.phase.gates[len(bundle.phase_residual_gates):]
        cls.algebra_outputs["phase"] = _native_matrix(phase_suffix, 11, remapped)
        cls.ideal_outputs = {key: _physical_columns(_logical_output(cls.oracle, key), cls.dirty)
                             for key in ("X", "Y", "phase")}
        cls.magnitude_score = sum(_score_operator(bundle, cls.outputs[key], cls.oracle[weights])
                                  for key, weights in (("X", "wx"), ("Y", "wy"))) / 2
        cls.phase_scores = tuple(_score_operator(bundle, cls.outputs["phase"], 2 * np.eye(2)[leaf])
                                 for leaf in (0, 1))
        cls.observed_preparation_errors = tuple(_operator_error(cls.outputs[key] - cls.ideal_outputs[key])
                                                for key in ("X", "phase"))
        cls.observed_magnitude_bias = np.linalg.norm(cls.magnitude_score - cls.oracle["magnitude"] * np.eye(2), ord=2)
        cls.observed_phase_bias = max(np.linalg.norm(score - wanted * np.eye(2), ord=2)
                                      for score, wanted in zip(cls.phase_scores, cls.oracle["phase"]))

    def test_exact_fixture_and_independently_certified_dyadic_inputs_at_fine_precision(self):
        oracle = self.oracle
        self.assertEqual(_norm_squared(oracle["a"]) + _norm_squared(oracle["b"]), 1)
        self.assertLess(2 * _norm_squared(oracle["b"]), 1)
        self.assertEqual(sum(_norm_squared(value) for value in oracle["state_exact"]), 1)
        distance_squared = 2 * (1 - oracle["cu"] * oracle["cv"])
        self.assertEqual(distance_squared, Fraction(262144, 4042387697))
        self.assertLess(distance_squared, Fraction(1, 64 ** 2))
        for q in (5, 16, 80, 129):
            with self.subTest(q=q):
                data = native.residual_qbp_data(q)
                self.assertEqual((data.angle_denominator, data.phase_denominator), (1024, 256))
                self.assertEqual((data.cos_u, data.sin_u, data.cos_v, data.sin_v),
                                 tuple(oracle[key] for key in ("cu", "su", "cv", "sv")))
                self.assertEqual(data.residual_root, oracle["a"])
                self.assertEqual(data.residual_tail, oracle["b"])
                self.assertEqual(data.target_state, oracle["state_exact"])
                self.assertEqual(data.magnitude_derivative, oracle["derivative_exact"])
                self.assertEqual(data.coarse_matrix, oracle["coarse_exact"])
                self.assertEqual(data.coarse_distance_squared, distance_squared)
                self.assertEqual(data.coarse_distance_bound, Fraction(1, 64))
                self.assertEqual((data.magnitude_x_weights, data.magnitude_y_weights), (oracle["wx"], oracle["wy"]))
                delta = Fraction(1, 1 << (2 * q + 20))
                self.assertEqual(data.input_error_bound, delta)
                # Independent, finer isqrt enclosure; do not reuse the
                # coefficient evaluator's midpoint or rounding certificate.
                scale = 1 << (2 * q + 32)
                root = isqrt(2 * scale * scale)
                lower, upper = Fraction(root, scale), Fraction(root + 1, scale)
                self.assertLess(lower * lower, 2)
                self.assertGreater(upper * upper, 2)
                root_error = sum((rounded - exact) ** 2 for rounded, exact in zip(data.programmed_root, oracle["a"]))
                tail_error = sum(max(abs(rounded - lower * exact), abs(rounded - upper * exact)) ** 2
                                 for rounded, exact in zip(data.programmed_scaled_tail, oracle["b"]))
                self.assertLessEqual(root_error, delta ** 2)
                self.assertLessEqual(tail_error, delta ** 2)
                self.assertEqual(data.coefficient_error_squared_bounds[0], root_error)
                self.assertLessEqual(tail_error, data.coefficient_error_squared_bounds[1])
                self.assertTrue(all(error <= delta ** 2 for error in data.coefficient_error_squared_bounds))
                for value in data.programmed_root + data.programmed_scaled_tail:
                    self.assertEqual(value.denominator & (value.denominator - 1), 0)

    def test_literal_coarse_observable_actual_inverse_and_complete_stream_ledgers(self):
        bundle = self.bundle
        np.testing.assert_allclose(_small_native(bundle.coarse_gates, (bundle.system,)), self.oracle["C"], atol=3e-15, rtol=0)
        self.assertEqual(bundle.coarse_inverse_gates, _inverse(bundle.coarse_gates))
        np.testing.assert_allclose(_small_native(bundle.coarse_inverse_gates, (bundle.system,)),
                                   self.oracle["C"].conj().T, atol=3e-15, rtol=0)
        observable = np.eye(4, dtype=complex)
        observable[2:, 2:] = H
        np.testing.assert_allclose(_small_native(bundle.controlled_observable_gates, (bundle.system, bundle.branch)),
                                   observable, atol=3e-15, rtol=0)
        for q in (5, 16, 80, 129):
            current = bundle if q == 5 else native.emit_residual_qbp(q)
            with self.subTest(q=q):
                self.assertEqual((current.nqubits, current.dirty_qubits), (q + 6, q + 2))
                self.assertEqual(current.coarse_inverse_gates, _inverse(current.coarse_gates))
                self.assertEqual(sum(name in ("T", "TDG") for name, _ in current.coarse_gates), 2)
                self.assertEqual(sum(name in ("T", "TDG") for name, _ in current.controlled_observable_gates), 2)
                remapped = tuple((name, tuple(current.mode if wire == current.phase_residual.mode else wire
                                              for wire in wires)) for name, wires in current.phase_residual.gates)
                self.assertEqual(current.phase_residual_gates, remapped)
                self.assertNotIn(current.branch, {wire for _, wires in remapped for wire in wires})
                branch_h = (("H", (current.branch,)),)
                self.assertEqual(current.pair_preparation_gates, branch_h + current.pair_residual.gates + current.coarse_gates)
                self.assertEqual(current.phase_preparation_gates, remapped + current.coarse_gates + branch_h)
                shared = current.pair_preparation_gates + current.controlled_observable_gates + current.coarse_inverse_gates + (("H", (current.system,)),)
                self.assertEqual(current.magnitude_x.gates, shared + branch_h)
                self.assertEqual(current.magnitude_y.gates, shared + (("SDG", (current.branch,)),) + branch_h)
                self.assertEqual(current.phase.gates, current.phase_preparation_gates + current.controlled_observable_gates
                                 + (("SDG", (current.branch,)),) + branch_h)
                self.assertEqual(current.pair_preparation_t_count, current.pair_residual.t_count + 2)
                self.assertEqual(current.phase_preparation_t_count, current.phase_residual.t_count + 2)
                for stream in (current.magnitude_x, current.magnitude_y, current.phase):
                    self.assertEqual(stream.gates, tuple(gate for stage in stream.stages for gate in stage.gates))
                    self.assertEqual(stream.t_count, sum(name in ("T", "TDG") for name, _ in stream.gates))
                    self.assertLess(stream.preparation_error_bound, Fraction(390, 1 << q))
                self.assertEqual(current.magnitude_x.t_count, current.pair_residual.t_count + 6)
                self.assertEqual(current.magnitude_y.t_count, current.magnitude_x.t_count)
                self.assertEqual(current.phase.t_count, current.phase_residual.t_count + 4)
                # Three coarse appearances and both controlled observables
                # remain charged across a magnitude/phase execution pair.
                self.assertEqual(current.magnitude_x.t_count + current.phase.t_count,
                                 current.pair_residual.t_count + current.phase_residual.t_count + 10)
        for stream in (bundle.magnitude_x, bundle.magnitude_y, bundle.phase):
            for name, wires in stream.gates:
                self.assertIn(name, ALPHABET)
                self.assertEqual(len(wires), 2 if name == "CX" else 1)
                self.assertTrue(all(0 <= wire < bundle.nqubits for wire in wires))

    def test_ideal_gradient_means_original_mathematical_reference_and_wrong_wrappers(self):
        oracle = self.oracle
        response = H @ oracle["psi"]
        self.assertAlmostEqual(2 * np.vdot(oracle["derivative"], response).real, oracle["magnitude"], delta=2e-15)
        np.testing.assert_allclose(2 * np.imag(oracle["psi"].conj() * response), oracle["phase"], atol=2e-15, rtol=0)
        weights = {key: np.asarray(oracle[name], dtype=float) for key, name in (("X", "wx"), ("Y", "wy"))}
        means = {key: _logical_score(_logical_output(oracle, key), weights[key]) for key in ("X", "Y")}
        self.assertAlmostEqual((means["X"] + means["Y"]) / 2, oracle["magnitude"], delta=3e-15)
        self.assertGreater(abs(means["Y"] / 2), 1e-5)
        phase_output = _logical_output(oracle, "phase")
        np.testing.assert_allclose([_logical_score(phase_output, 2 * np.eye(2)[leaf]) for leaf in (0, 1)],
                                   oracle["phase"], atol=3e-15, rtol=0)
        for inverse in ("omitted", "forward"):
            wrong = sum(_logical_score(_logical_output(oracle, key, inverse=inverse), weights[key])
                        for key in ("X", "Y")) / 2
            self.assertGreater(abs(wrong - oracle["magnitude"]), 1)
        changed_phase = sum(_logical_score(_logical_output(oracle, key, branch_phase=1j), weights[key])
                            for key in ("X", "Y")) / 2
        self.assertGreater(abs(changed_phase - oracle["magnitude"]), 1)
        wrong_phase = _logical_output(oracle, "phase", inverse="extra")
        self.assertGreater(abs(_logical_score(wrong_phase, (2, 0)) - oracle["phase"][0]), .5)
        # This comparison is mathematical only: no native fine-frame word
        # or gate count is claimed for the non-native target W'. The same
        # calibrated observable H is used by both ideal protocols.
        psi = oracle["psi"]
        frame = np.column_stack((psi, [-psi[1].conjugate(), psi[0].conjugate()]))
        selected_h = np.eye(4, dtype=complex)
        selected_h[2:, 2:] = H
        original = np.kron(H, H @ frame.conj().T) @ selected_h @ np.tile(psi, 2) / np.sqrt(2)
        self.assertAlmostEqual(_logical_score(original, (2, -2)), oracle["magnitude"], delta=3e-15)

    def test_complete_native_streams_keep_all_outcomes_and_match_independent_algebra(self):
        bundle = self.bundle
        for key in ("X", "Y", "phase"):
            with self.subTest(stream=key):
                actual = self.outputs[key]
                np.testing.assert_allclose(actual, self.algebra_outputs[key], atol=4 * ATOL, rtol=0)
                np.testing.assert_allclose(actual.conj().T @ actual, np.eye(2), atol=4 * ATOL, rtol=0)
                np.testing.assert_allclose(np.sum(abs(actual) ** 2, axis=0), 1, atol=4 * ATOL, rtol=0)
        for key, bound in (("X", bundle.magnitude_x.preparation_error_bound), ("phase", bundle.phase.preparation_error_bound)):
            self.assertLess(_operator_error(self.outputs[key] - self.ideal_outputs[key]), float(bound))
        # Measured selected-column errors are 0.296773/0.329699, with
        # magnitude bias 0.036290 and phase-coordinate bias 0.080265.
        # These are finite-fixture regressions, not uniform-q estimates.
        self.assertLess(max(self.observed_preparation_errors), .4)
        self.assertLess(self.observed_magnitude_bias, .05)
        self.assertLess(self.observed_phase_bias, .1)
        expected_mag = sum(_score_operator(bundle, self.algebra_outputs[key], self.oracle[weights])
                           for key, weights in (("X", "wx"), ("Y", "wy"))) / 2
        np.testing.assert_allclose(self.magnitude_score, expected_mag, atol=10 * ATOL, rtol=0)
        for leaf, score in enumerate(self.phase_scores):
            np.testing.assert_allclose(score, _score_operator(bundle, self.algebra_outputs["phase"], 2 * np.eye(2)[leaf]),
                                       atol=8 * ATOL, rtol=0)
        np.testing.assert_allclose(sum(self.phase_scores), np.zeros((2, 2)), atol=8 * ATOL, rtol=0)
        # Flag leakage is present, and the normalized probabilities above
        # retain it instead of conditioning on the returned-zero event.
        indices = np.arange(2048)
        rejected = ((indices >> bundle.target) & 1) | ((indices >> bundle.mode) & 1)
        for key in ("X", "phase"):
            self.assertGreater(np.sum(abs(self.outputs[key][rejected != 0, 0]) ** 2), 1e-4)

    def test_exact_integer_histogram_decoders_and_unprojected_phase_coordinates(self):
        hx, hy, shots = [3, -2], [-1, 2], 13
        original = hx.copy(), hy.copy()
        expected = sum(weight * count for weight, count in zip(self.oracle["wx"], hx))
        expected += sum(weight * count for weight, count in zip(self.oracle["wy"], hy))
        self.assertEqual(native.decode_magnitude_histograms(hx, hy, shots), expected / shots)
        self.assertEqual((hx, hy), original)
        self.assertEqual(native.decode_magnitude_histograms((0, 0), (0, 0), 7), 0)
        huge = 1 << 90
        self.assertEqual(native.decode_magnitude_histograms((huge, 0), (0, 0), huge + 7),
                         self.oracle["wx"][0] * huge / (huge + 7))
        phase = native.decode_phase_histogram((2, -1), 5)
        self.assertEqual(phase, (Fraction(4, 5), Fraction(-2, 5)))
        self.assertEqual(sum(phase), Fraction(2, 5))
        theta = np.pi / 4 + 2 * np.arctan(1 / ANGLE_DENOMINATOR)
        phi = np.pi / 4 + 2 * np.arctan(1 / PHASE_DENOMINATOR)
        ry = np.array([[1, -1], [1, 1]]) / np.sqrt(2)
        rz = np.diag([(1 - 1j) / np.sqrt(2), (1 + 1j) / np.sqrt(2)])
        decoded = decode_complex_coarse_frame_histograms([theta], [ry], [rz], [-phi, phi], hx, hy, shots)
        np.testing.assert_allclose(decoded, [float(expected / shots)], atol=2e-15, rtol=0)

    def test_rejects_invalid_precision_shots_and_signed_histograms(self):
        for q in (True, False, 5., Fraction(5), 4, -1):
            with self.subTest(q=q), self.assertRaises(ValueError):
                native.residual_qbp_data(q)
            with self.subTest(q=q), self.assertRaises(ValueError):
                native.emit_residual_qbp(q)
        for shots in (0, -1, True, False, 1., Fraction(1)):
            with self.subTest(shots=shots), self.assertRaises(ValueError):
                native.decode_magnitude_histograms((0, 0), (0, 0), shots)
            with self.subTest(shots=shots), self.assertRaises(ValueError):
                native.decode_phase_histogram((0, 0), shots)
        for counts in ((), (1,), (1, 0, 0)):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                native.decode_magnitude_histograms(counts, (0, 0), 5)
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                native.decode_phase_histogram(counts, 5)
        for value in (True, False, 1., Fraction(1), "1"):
            with self.subTest(count=value), self.assertRaises(TypeError):
                native.decode_magnitude_histograms((value, 0), (0, 0), 5)
            with self.subTest(count=value), self.assertRaises(TypeError):
                native.decode_phase_histogram((value, 0), 5)
        with self.assertRaises(ValueError):
            native.decode_magnitude_histograms((2, -2), (1, -1), 5)
        with self.assertRaises(ValueError):
            native.decode_phase_histogram((3, -3), 5)


if __name__ == "__main__":
    unittest.main()
