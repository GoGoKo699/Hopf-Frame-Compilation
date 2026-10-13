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
    dyadic_precision,
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


def _selected_term(length, shift, coefficient, separate_flags=True):
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
    target_word = np.linalg.matrix_power(J, shift)
    for column in range(size):
        flag, source, target, dirty = np.unravel_index(column, shape)
        flag ^= (2 if separate_flags else 1) * (source < shift)
        destination = (source - shift) % length
        for output in range(2):
            row = np.ravel_multi_index((flag, destination, output, dirty), shape)
            permutation[row, column] = target_word[output, target]
    return permutation @ scalar


class HeterogeneousPrecisionSourceTests(unittest.TestCase):
    def test_exact_finite_source_shift_error_including_cap(self):
        for beta in (Fraction(1, 2), Fraction(1, 4), Fraction(3, 8)):
            for length in (2, 4, 7):
                weights = [beta ** j for j in range(length)]
                norm_squared = sum((value * value for value in weights), Fraction())
                for shift in range(length + 1):
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
        # Two aligned four-label components, with a full eight-label feature
        # bank. Feature inputs are ordinary occupied data in this check.
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

        # Even an identity merge must not introduce this cross-boundary
        # splitter: it changes which geometric source the mode requires.
        crossing = plane(0, 8 + 4, .37)
        self.assertGreater(np.linalg.norm(projector @ crossing - crossing @ projector, 2), .3)
        wrong = np.kron(np.eye(length), crossing) @ source - source @ crossing
        self.assertGreater(np.linalg.norm(wrong, 2), .03)
        self.assertGreater(np.linalg.norm(wrong[:, 8 + 4]), .03)

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
            source_shift_residual_squared(Fraction(1, 2), 4, 5)
        for order in ("2", True, 0):
            with self.subTest(order=order), self.assertRaises(ValueError):
                kernel_source_error_bound(Fraction(1, 4), Fraction(1, 2), order, 4,
                                          Fraction(1, 16), 8)


if __name__ == "__main__":
    unittest.main()
