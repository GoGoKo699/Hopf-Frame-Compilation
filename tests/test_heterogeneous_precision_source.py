"""Finite small-ratio source checks, including all-column failure history.

The dense dilations below are diagnostics, not native synthesis routines.
Endpoint checks manipulate exponents and small integers, never 2**(-2**n).
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.heterogeneous_precision_source import (
    cayley_kernel_coefficients,
    cayley_pair,
    common_source_kernel_coefficients,
    common_source_kernel_error_bound,
    common_source_parameters,
    common_source_ratio,
    dyadic_precision,
    exponent_metadata,
    kernel_source_error_bound,
    small_ratio_parameters,
    source_shift_residual_squared,
)


J = np.array([[0., -1.], [1., 0.]])
ATOL = 4e-13


def _source(beta, length):
    vector = float(beta) ** np.arange(length)
    return vector / np.linalg.norm(vector)


def _shift(length, shift):
    result = np.zeros((length, length))
    for column in range(shift, length):
        result[column - shift, column] = 1
    return result


def _rotation(parameter):
    # Independent trigonometric construction, rather than the rational helper.
    angle = 2 * np.arctan(float(parameter))
    return np.cos(angle) * np.eye(2) + np.sin(angle) * J


def _kernel(beta, mantissa, order=2, length=4, probability=Fraction(1, 16), bits=8):
    coefficients = cayley_kernel_coefficients(mantissa, order, probability, bits)
    result = (1 - float(probability)) * np.eye(2 * length)
    for shift, coefficient in enumerate(coefficients, 1):
        result += float(probability * coefficient / order) * np.kron(
            _shift(length, shift), np.linalg.matrix_power(J, shift))
    return result


def _direct_sum(matrices):
    size = sum(len(matrix) for matrix in matrices)
    result = np.zeros((size, size), dtype=complex)
    start = 0
    for matrix in matrices:
        stop = start + len(matrix)
        result[start:stop, start:stop] = matrix
        start = stop
    return result


def _unitary_dilation(contraction):
    left, singular, right = np.linalg.svd(contraction)
    complement = np.sqrt(np.maximum(0., 1 - singular ** 2))
    top = (left * complement) @ left.conj().T
    bottom = (right.conj().T * complement) @ right
    return np.block([[contraction, top], [bottom, -contraction.conj().T]])


def _selected_term(length, shift, coefficient, separate_flags=True, *,
                   target_power=None, compare_modulo=False):
    """Actual full permutation/scalar word, with a nontrivial dirty rejection."""
    flags = 4 if separate_flags else 2
    shape = (flags, length, 2, 2)
    size = int(np.prod(shape))
    scalar = np.zeros((size, size), dtype=complex)
    # The accepted block is f I_dirty, while the rejected block moves dirty.
    polar = np.array([[0., 1j], [1., 0.]])
    f = float(coefficient)
    rejection = np.sqrt(1 - f * f)
    local = np.block([[f * np.eye(2), -rejection * polar.conj().T],
                      [rejection * polar, f * np.eye(2)]])
    for underflow in range(2 if separate_flags else 1):
        for source in range(length):
            for target in range(2):
                indices = [np.ravel_multi_index((2 * underflow + flag if separate_flags else flag,
                                                 source, target, dirty), shape)
                           for flag in range(2) for dirty in range(2)]
                scalar[np.ix_(indices, indices)] = local
    permutation = np.zeros_like(scalar)
    target_word = np.linalg.matrix_power(J, shift if target_power is None else target_power)
    for column in range(size):
        flag, source, target, dirty = np.unravel_index(column, shape)
        threshold = shift % length if compare_modulo else shift
        flag ^= (2 if separate_flags else 1) * (source < threshold)
        destination = (source - shift) % length
        for output in range(2):
            row = np.ravel_multi_index((flag, destination, output, dirty), shape)
            permutation[row, column] = target_word[output, target]
    return permutation @ scalar


def _plane(size, left, right, angle):
    matrix = np.eye(size)
    matrix[np.ix_((left, right), (left, right))] = (
        np.cos(angle) * np.eye(2) + np.sin(angle) * J)
    return matrix


def _common_plane_kernel(size, left, right, beta, mantissa, multiplier, remainder,
                         *, order=2, length=4, probability=Fraction(1, 64), bits=10):
    result = (1 - float(probability)) * np.eye(size * length)
    coefficients = common_source_kernel_coefficients(
        mantissa, order, probability, bits, multiplier=multiplier, remainder=remainder)
    for degree, coefficient in enumerate(coefficients, 1):
        local = np.zeros((size, size))
        local[np.ix_((left, right), (left, right))] = np.linalg.matrix_power(J, degree)
        result += float(probability * coefficient / order) * np.kron(
            _shift(length, multiplier * degree), local)
    return result


class HeterogeneousPrecisionSourceTests(unittest.TestCase):
    def test_exact_finite_source_shift_error_including_cap(self):
        for beta in (Fraction(1, 2), Fraction(1, 4), Fraction(3, 8)):
            for length in (2, 4, 7):
                weights = [beta ** j for j in range(length)]
                norm_squared = sum((value * value for value in weights), Fraction())
                for shift in range(2 * length + 2):
                    residual = [(weights[j + shift] if j + shift < length else 0) -
                                beta ** shift * weights[j] for j in range(length)]
                    exact = sum((value * value for value in residual), Fraction()) / norm_squared
                    self.assertEqual(source_shift_residual_squared(beta, length, shift), exact)
                    source = _source(beta, length)
                    observed = _shift(length, shift) @ source - float(beta ** shift) * source
                    self.assertAlmostEqual(float(observed @ observed), float(exact), delta=2e-16)
                    self.assertLessEqual(exact, beta ** (2 * length) / (1 - beta ** (2 * length)))

    def test_exact_dyadic_coefficients_and_symbolic_large_precision_margins(self):
        for height in (2, 4, 8, 16, 32):
            parameters = small_ratio_parameters(height, suffix_constant=12)
            self.assertEqual(parameters.order, height)
            self.assertEqual(parameters.source_length, 2 * height)
            self.assertEqual(parameters.rare_probability, Fraction(1, 8 * height))
            self.assertEqual(parameters.rescaling_bits, 6 + 2 * (height.bit_length() - 1))
            for mantissa in (Fraction(0), Fraction(1), Fraction(-3, 4), Fraction(17, 32)):
                coefficients = cayley_kernel_coefficients(
                    mantissa, parameters.order, parameters.rare_probability, parameters.rescaling_bits)
                bound = parameters.coefficient_precision_bound(dyadic_precision(mantissa))
                beta = Fraction(1, 4)
                parameter = beta * mantissa / (1 << parameters.rescaling_bits)
                for degree, coefficient in enumerate(coefficients, 1):
                    self.assertLessEqual(abs(coefficient), Fraction(1, 4))
                    self.assertLessEqual(dyadic_precision(coefficient), bound)
                    self.assertEqual(parameters.rare_probability * coefficient * beta ** degree /
                                     parameters.order,
                                     2 * (1 - parameters.rare_probability) * parameter ** degree)
        for height in (16, 64, 256, 1 << 20, 1 << 80):
            parameters = small_ratio_parameters(height, suffix_constant=12)
            self.assertTrue(parameters.cap_exponent_at_least_three_precision_lengths)
            self.assertEqual(parameters.minimum_exponent_log2,
                             height + 1 - (height.bit_length() - 1))
            # E_min*P=2N, and M(E_min-B)>=3N. These compare exponents;
            # no huge integer N or underflowed endpoint amplitude is created.
            self.assertEqual(Fraction(2, height) * parameters.order, 2)
            self.assertLessEqual((2 * height * parameters.rescaling_bits - 1).bit_length(), height)
            self.assertLessEqual(height * parameters.rare_probability, Fraction(1, 8))

    def test_truncated_kernel_full_source_columns_and_inactive_baseline(self):
        for beta in (Fraction(1, 2), Fraction(1, 4)):
            for mantissa in (Fraction(0), Fraction(1), Fraction(-3, 4)):
                actual = _kernel(beta, mantissa)
                source_column = np.kron(_source(beta, 4)[:, None], np.eye(2))
                parameter = beta * mantissa / 256
                expected = Fraction(15, 16) * source_column @ _rotation(parameter)
                error = np.linalg.norm(actual @ source_column - np.asarray(expected, dtype=float), 2)
                bound = kernel_source_error_bound(beta, mantissa, 2, 4, Fraction(1, 16), 8)
                self.assertLessEqual(error, float(bound) + ATOL)
                self.assertLessEqual(np.linalg.norm(actual, 2), 1 + ATOL)
                cosine, sine = cayley_pair(parameter)
                self.assertEqual(cosine * cosine + sine * sine, 1)
                np.testing.assert_allclose(float(cosine) * np.eye(2) + float(sine) * J,
                                           _rotation(parameter), atol=ATOL, rtol=0)
                if mantissa == 0:
                    np.testing.assert_array_equal(actual, 15 / 16 * np.eye(8))
                    self.assertGreater(np.linalg.norm(actual - np.eye(8), 2), .06)
                else:
                    self.assertGreater(error, 1e-6)

    def test_shift_and_scalar_failures_need_distinct_retained_flags(self):
        length, shift, coefficient = 4, 1, Fraction(3, 5)
        word = _selected_term(length, shift, coefficient)
        np.testing.assert_allclose(word.conj().T @ word, np.eye(len(word)), atol=ATOL, rtol=0)
        accepted_size = 4 * length
        accepted = word[:accepted_size, :accepted_size]
        expected = float(coefficient) * np.kron(_shift(length, shift), np.kron(J, np.eye(2)))
        np.testing.assert_allclose(accepted, expected, atol=ATOL, rtol=0)
        # Source index zero underflows. Neither an occupied dirty input nor a
        # simultaneous scalar rejection may turn that into an accepted wrap.
        np.testing.assert_array_equal(accepted[:, :4], np.zeros((accepted_size, 4)))
        dirty_failure = word[accepted_size:2 * accepted_size, :accepted_size]
        self.assertGreater(np.linalg.norm(dirty_failure), .5)
        wrong = _selected_term(length, shift, coefficient, separate_flags=False)
        np.testing.assert_allclose(wrong.conj().T @ wrong, np.eye(len(wrong)), atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(wrong[:accepted_size, :4]), .7)

    def test_heterogeneous_history_product_preserves_one_uniform_gamma(self):
        betas = (Fraction(1, 2), Fraction(1, 4))
        rows = ((Fraction(3, 4), Fraction(-1, 2)),
                (Fraction(0), Fraction(3, 4)),
                (Fraction(1, 2), Fraction(0)))
        data_size, columns, history_size = 32, 8, 4
        # Rectangular source embeddings have separate input and output sizes.
        source = np.zeros((data_size, columns), dtype=complex)
        for branch, beta in enumerate(betas):
            source[16 * branch:16 * (branch + 1), 4 * branch:4 * (branch + 1)] = np.kron(
                _source(beta, 4)[:, None], np.eye(4))
        state = np.zeros((history_size, 2, data_size, columns), dtype=complex)
        state[0, 0] = source
        no_history = np.vstack((source, np.zeros_like(source)))
        compressed, target = source.copy(), np.eye(columns, dtype=complex)
        total_bound = Fraction()
        for stage, mantissas in enumerate(rows):
            contractions = [np.kron(_kernel(beta, mantissa), np.eye(2))
                            for beta, mantissa in zip(betas, mantissas)]
            contraction = _direct_sum(contractions)
            word = _unitary_dilation(contraction)
            np.testing.assert_allclose(word.conj().T @ word, np.eye(2 * data_size), atol=ATOL, rtol=0)
            for history in range(history_size):
                state[history] = (word @ state[history].reshape(2 * data_size, columns)).reshape(
                    2, data_size, columns)
            if stage + 1 < len(rows):
                state[:, 1] = np.roll(state[:, 1].copy(), 1, axis=0)
            no_history = word @ no_history
            compressed = contraction @ compressed
            target = _direct_sum([np.kron(_rotation(beta * mantissa / 256), np.eye(2))
                                  for beta, mantissa in zip(betas, mantissas)]) @ target
            total_bound += max(kernel_source_error_bound(beta, mantissa, 2, 4, Fraction(1, 16), 8)
                               for beta, mantissa in zip(betas, mantissas))
        np.testing.assert_allclose(state[0, 0], compressed, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(no_history[:data_size] - compressed, 2), .01)
        gamma = Fraction(15, 16) ** len(rows)
        self.assertGreaterEqual(gamma, 1 - len(rows) * Fraction(1, 16))
        self.assertLessEqual(np.linalg.norm(compressed - float(gamma) * source @ target, 2),
                             float(total_bound) + ATOL)
        self.assertAlmostEqual(float(np.linalg.norm(state.ravel())), np.sqrt(columns), delta=ATOL)

    def test_pruned_feature_layout_preserves_class_source_on_all_mode_columns(self):
        # Historical class-dependent source, not the newer common source.
        # Feature inputs are ordinary occupied data in this check.
        mode_count, length = 16, 4
        betas = (Fraction(1, 2), Fraction(1, 4))
        source = np.zeros((length * mode_count, mode_count))
        class_of = [(label % 8) // 4 for label in range(mode_count)]
        for label, branch in enumerate(class_of):
            source[label::mode_count, label] = _source(betas[branch], length)
        projector = np.diag([branch == 0 for branch in class_of]).astype(float)

        def plane(left, right, angle):
            matrix = np.eye(mode_count)
            matrix[np.ix_((left, right), (left, right))] = (
                np.cos(angle) * np.eye(2) + np.sin(angle) * J)
            return matrix

        stages = []
        for root in (0, 4):
            # Bottom-up carrier merge, followed by the parent/feature splitter.
            # No node crossing the two component roots is present.
            for left, right in ((root, root + 1), (root + 2, root + 3),
                                (root, root + 2)):
                stages.extend((plane(left, right, -.23), plane(left, 8 + right, .17)))
        product = np.eye(mode_count)
        for stage in stages:
            np.testing.assert_allclose(projector @ stage, stage @ projector, atol=ATOL, rtol=0)
            np.testing.assert_allclose(np.kron(np.eye(length), stage) @ source,
                                       source @ stage, atol=ATOL, rtol=0)
            product = stage @ product
        np.testing.assert_allclose(np.kron(np.eye(length), product) @ source,
                                   source @ product, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(product[:8, 8:]), .3)

        # For this CLASS-DEPENDENT source, the crossing changes the geometric
        # source required by the mode. The common-source test below permits it.
        crossing = plane(0, 8 + 4, .37)
        self.assertGreater(np.linalg.norm(projector @ crossing - crossing @ projector, 2), .3)
        wrong = np.kron(np.eye(length), crossing) @ source - source @ crossing
        self.assertGreater(np.linalg.norm(wrong, 2), .03)
        self.assertGreater(np.linalg.norm(wrong[:, 8 + 4]), .03)

    def test_common_exponent_remainders_and_symbolic_coefficient_ledger(self):
        spacing, bits, order, probability = 11, 10, 4, Fraction(1, 64)
        beta = common_source_ratio(spacing, bits)
        self.assertEqual(beta, Fraction(1, 2))
        for exponent in (11, 14, 25, 33, 37):
            multiplier, remainder = exponent_metadata(exponent, spacing)
            self.assertEqual(exponent, multiplier * spacing + remainder)
            self.assertTrue(0 <= remainder < spacing)
            for mantissa in (Fraction(0), Fraction(3, 4), Fraction(-5, 8)):
                coefficients = common_source_kernel_coefficients(
                    mantissa, order, probability, bits,
                    multiplier=multiplier, remainder=remainder)
                parameter = mantissa / (1 << exponent)
                for degree, coefficient in enumerate(coefficients, 1):
                    self.assertEqual(probability * coefficient * beta ** (multiplier * degree) /
                                     order, 2 * (1 - probability) * parameter ** degree)
                    self.assertLessEqual(abs(coefficient), Fraction(1, 2))
                    self.assertLessEqual(dyadic_precision(coefficient),
                                         degree * (dyadic_precision(mantissa) + remainder +
                                                   multiplier * bits))
                if remainder and mantissa:
                    wrong = common_source_kernel_coefficients(
                        mantissa, order, probability, bits,
                        multiplier=multiplier, remainder=0)
                    self.assertNotEqual(wrong, coefficients)
        for height in (16, 64, 256, 1 << 20, 1 << 80):
            for maximum in (None, height * (1 << ((height.bit_length() - 1) // 2)),
                            3 * height ** 2 // (height.bit_length() - 1)):
                parameters = common_source_parameters(
                    height, suffix_constant=20, maximum_multiplier=maximum)
                log_height = height.bit_length() - 1
                self.assertEqual(parameters.spacing_log2 + 2 * log_height, height)
                self.assertEqual(parameters.source_length, 4 * height * height)
                self.assertTrue(parameters.cap_exponent_at_least_three_precision_lengths)
                # q<=n*Delta+overhead, so n*q<=N+n*overhead. Compare the
                # polynomial overhead's bit length against log2(N)=n.
                overhead = parameters.coefficient_precision_overhead(3 * log_height + 7)
                self.assertEqual(overhead, height * (3 * log_height + 7 +
                                 parameters.maximum_multiplier * parameters.rescaling_bits))
                if height >= 64:
                    self.assertTrue(parameters.stage_precision_at_most_twice_endpoint(
                        3 * log_height + 7))
                self.assertEqual(parameters.suffix_bits, 20 * (height + 1).bit_length())

    def test_saturated_shifts_are_zero_and_target_power_is_degree(self):
        length, coefficient = 4, Fraction(3, 5)
        for multiplier, degree in ((2, 2), (3, 2), (4, 1)):
            shift = multiplier * degree
            word = _selected_term(length, shift, coefficient, target_power=degree)
            np.testing.assert_allclose(word.conj().T @ word, np.eye(len(word)), atol=ATOL, rtol=0)
            accepted_size = 4 * length
            np.testing.assert_array_equal(word[:accepted_size, :accepted_size],
                                          np.zeros((accepted_size, accepted_size)))
            wrong = _selected_term(length, shift, coefficient, target_power=degree,
                                   compare_modulo=True)
            self.assertGreater(np.linalg.norm(wrong[:accepted_size, :accepted_size], 2), .5)
        # Before saturation, J**nu and J**(a*nu) are already distinguishable.
        word = _selected_term(length, 2, coefficient, target_power=1)
        expected = float(coefficient) * np.kron(_shift(length, 2), np.kron(J, np.eye(2)))
        np.testing.assert_allclose(word[:16, :16], expected, atol=ATOL, rtol=0)
        wrong = _selected_term(length, 2, coefficient, target_power=2)
        self.assertGreater(np.linalg.norm(wrong[:16, :16] - expected, 2), .5)

    def test_pair_symmetric_metadata_uncomputes_on_occupied_inputs(self):
        # Axes: logical pair/target label, metadata, source, underflow. Metadata
        # storage is tested on BOTH input values, not silently initialized.
        shape = (4, 2, 4, 2)
        size = int(np.prod(shape))

        def word(symmetric):
            result = np.zeros((size, size))
            for column in range(size):
                label, metadata, source, flag = np.unravel_index(column, shape)
                tag = label // 2 if symmetric else label % 2
                loaded = metadata ^ tag
                shift = 1 + 2 * loaded
                for target in range(2):
                    output_label = 2 * (label // 2) + target
                    output_tag = output_label // 2 if symmetric else output_label % 2
                    output = (output_label, loaded ^ output_tag,
                              (source - shift) % 4, flag ^ (source < shift))
                    result[np.ravel_multi_index(output, shape), column] = J[target, label % 2]
            return result

        actual = word(True)
        np.testing.assert_array_equal(actual.T @ actual, np.eye(size))
        for metadata in range(2):
            columns = [index for index in range(size)
                       if np.unravel_index(index, shape)[1] == metadata]
            other = [index for index in range(size)
                     if np.unravel_index(index, shape)[1] != metadata]
            np.testing.assert_array_equal(actual[np.ix_(other, columns)],
                                          np.zeros((len(other), len(columns))))
        wrong = word(False)
        np.testing.assert_array_equal(wrong.T @ wrong, np.eye(size))
        empty = [index for index in range(size) if np.unravel_index(index, shape)[1] == 0]
        occupied = [index for index in range(size) if np.unravel_index(index, shape)[1] == 1]
        self.assertEqual(np.linalg.norm(wrong[np.ix_(occupied, empty)], 2), 1.)

    def test_common_source_interleaved_scales_and_zero_merge_splitter_history(self):
        # A connected four-label tree, occupied feature inputs, and a single
        # geometric source. Different exponents occur on overlapping root
        # planes. A zero right-child merge STILL receives its splitter.
        beta, spacing, bits = Fraction(1, 2), 11, 10
        probability, length, modes = Fraction(1, 64), 4, 8
        source = np.kron(_source(beta, length)[:, None], np.eye(modes))
        specifications = ((0, 1, 1, 3, Fraction(3, 4)),
                          (2, 3, 3, 2, Fraction(0)),
                          (0, 2, 2, 1, Fraction(-5, 8)))
        data_size, history_size = modes * length, 4
        state = np.zeros((history_size, 2, data_size, modes), dtype=complex)
        state[0, 0] = source
        compressed, ideal = source.copy(), np.eye(modes)
        no_history = np.vstack((source, np.zeros_like(source)))
        total_bound = Fraction()
        for stage, (left, right, multiplier, remainder, mantissa) in enumerate(specifications):
            contraction = _common_plane_kernel(modes, left, right, beta, mantissa,
                                                multiplier, remainder)
            dilation = _unitary_dilation(contraction)
            splitter = _plane(modes, left, 4 + right, .23)
            lifted = np.kron(np.eye(length), splitter)
            unitary = _direct_sum((lifted, lifted)) @ dilation
            np.testing.assert_allclose(unitary.conj().T @ unitary,
                                       np.eye(2 * data_size), atol=ATOL, rtol=0)
            # Full occupied-feature intertwining holds, even at zero merges.
            np.testing.assert_allclose(lifted @ source, source @ splitter, atol=ATOL, rtol=0)
            for history in range(history_size):
                state[history] = (unitary @ state[history].reshape(2 * data_size, modes)).reshape(
                    2, data_size, modes)
            if stage + 1 < len(specifications):
                state[:, 1] = np.roll(state[:, 1].copy(), 1, axis=0)
            no_history = unitary @ no_history
            compressed = lifted @ contraction @ compressed
            parameter = mantissa / (1 << (multiplier * spacing + remainder))
            ideal = splitter @ _plane(modes, left, right, 2 * np.arctan(float(parameter))) @ ideal
            total_bound += common_source_kernel_error_bound(
                beta, mantissa, 2, length, probability, bits,
                multiplier=multiplier, remainder=remainder)
        np.testing.assert_allclose(state[0, 0], compressed, atol=ATOL, rtol=0)
        self.assertLessEqual(np.linalg.norm(compressed - float((1 - probability) ** 3) *
                                           source @ ideal, 2), float(total_bound) + ATOL)
        self.assertGreater(np.linalg.norm(no_history[:data_size] - compressed, 2), .01)

        # Independently check why the identity-merge splitter must remain.
        # Q's complete leaf columns have top block Gamma W^T only with equal
        # child attenuation before the root merge.
        left_angle, root_angle, split_angle = .31, -.42, .23
        coarse = _plane(4, 0, 1, left_angle) @ _plane(4, 0, 2, root_angle)
        correct = np.eye(8)
        omitted = np.eye(8)
        exact_stages = (_plane(8, 0, 1, -left_angle), _plane(8, 0, 5, split_angle),
                        _plane(8, 2, 7, split_angle), _plane(8, 0, 2, -root_angle),
                        _plane(8, 0, 6, split_angle))
        for index, stage in enumerate(exact_stages):
            correct = stage @ correct
            if index != 2:
                omitted = stage @ omitted
            np.testing.assert_allclose(np.kron(np.eye(length), stage) @ source,
                                       source @ stage, atol=ATOL, rtol=0)
        attenuation = np.cos(split_angle)
        gamma_rows = np.diag((attenuation ** 2, 1., attenuation, 1.))
        np.testing.assert_allclose(correct[:4, :4], gamma_rows @ coarse.T, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(omitted[:4, :4] - gamma_rows @ coarse.T, 2), .02)

    def test_invalid_parameters_are_rejected_without_rounding(self):
        for height in (0, 1, 3, True, 2.0):
            with self.subTest(height=height), self.assertRaises(ValueError):
                small_ratio_parameters(height, suffix_constant=12)
        for value in (Fraction(1, 3), .5, True):
            with self.subTest(value=value), self.assertRaises(ValueError):
                dyadic_precision(value)
        for mantissa, probability, bits in ((Fraction(3, 2), Fraction(1, 16), 8),
                                           (Fraction(1, 2), Fraction(3, 16), 8),
                                           (Fraction(1), Fraction(1, 16), 0)):
            with self.assertRaises(ValueError):
                cayley_kernel_coefficients(mantissa, 2, probability, bits)
        with self.assertRaises(ValueError):
            source_shift_residual_squared(Fraction(1, 2), 4, -1)
        for order in ("2", True, 0):
            with self.subTest(order=order), self.assertRaises(ValueError):
                kernel_source_error_bound(Fraction(1, 4), Fraction(1, 2), order, 4,
                                          Fraction(1, 16), 8)
        for height in (2, 8, 32):
            with self.assertRaises(ValueError):
                common_source_parameters(height, suffix_constant=20)
        with self.assertRaises(ValueError):
            common_source_ratio(4, 4)
        with self.assertRaises(ValueError):
            exponent_metadata(7, 0)
        with self.assertRaises(ValueError):
            common_source_kernel_coefficients(Fraction(1, 2), 2, Fraction(1, 16), 8,
                                              multiplier=0, remainder=1)


if __name__ == "__main__":
    unittest.main()
