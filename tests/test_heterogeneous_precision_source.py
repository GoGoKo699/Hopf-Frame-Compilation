"""Finite small-ratio source checks, including all-column failure history.

The dense dilations below are diagnostics, not native synthesis routines.
Endpoint checks manipulate exponents and small integers, never 2**(-2**n).
"""
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.heterogeneous_precision_source import (
    bottom_precision_trace,
    cayley_kernel_coefficients,
    cayley_pair,
    clipped_branch_decision,
    clipped_source_parameters,
    common_source_kernel_coefficients,
    common_source_kernel_error_bound,
    common_source_parameters,
    common_source_ratio,
    dyadic_precision,
    exponent_metadata,
    kernel_source_error_bound,
    small_ratio_parameters,
    source_shift_residual_squared,
    splitter_gauge_coefficients,
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


def _embed_modes(count, modes, local, dirty_size=2):
    result = np.eye(count * dirty_size, dtype=complex)
    indices = [dirty_size * mode + dirty for mode in modes for dirty in range(dirty_size)]
    result[np.ix_(indices, indices)] = local
    return result


def _mixed_chart_fixture(zero_left=False):
    """Actual three-port unitaries and an independent scalar tree recursion."""
    size, height, identity = 4, 2, np.eye(2)
    eta = 1 / (1000 * height ** 2)
    ell = np.array([1, -1j]) / np.sqrt(2)
    eigenbasis = np.column_stack((ell, ell.conj()))
    pair_basis = np.kron(_direct_sum((eigenbasis, np.ones((1, 1)))), identity)
    polars = {1: np.array([[0, 1j], [1, 0]]), 2: np.diag([1, 1j])}
    nodes, gauges = {}, {1: 1 + 0j}
    target = np.eye(size, dtype=complex)
    coarse = target.copy()
    coefficients, corrections = np.ones(size, complex), np.zeros(size, complex)
    for vertex, (left, detail, angle) in enumerate(((0, 2, .37), (0, 1, .49), (2, 3, .23)), 1):
        if zero_left and vertex == 2:
            angle = 0.
        c, s = np.cos(angle), np.sin(angle)
        rotation = np.array([[c, -s], [s, c]])
        target[[left, detail]] = rotation @ target[[left, detail]]
        node = dict(left=left, detail=detail, target=(c, s), gauge=gauges[vertex])
        if angle == 0:
            gate, feature = np.eye(2), np.zeros(2)
            loss, accepted, rejected, scalar = 0., 1., np.zeros((2, 2)), np.eye(4)
            local_word = np.eye(6)
            gauges[2 * vertex] = gauges[2 * vertex + 1] = gauges[vertex]
            coefficients[detail], corrections[detail] = 1 / gauges[vertex], 0
        elif vertex != 3:
            tau = eta * min(c, s) / 200
            z = (1 - tau) * np.exp(2j * angle + 1j * tau / 200)
            a, b = (1 + z) / 2, (z - 1) / (2j)
            gate, feature = np.array([[a, -b], [b, a]]), ell
            loss, accepted = 1 - abs(z) ** 2, z.conjugate()
            rejected = np.sqrt(loss) * polars[vertex]
            # Upper-right rejection is scalar, whereas the lower-left one is
            # not. Reversing replay chronology therefore loses the dirty polar.
            scalar = np.block([[accepted * identity, -np.sqrt(loss) * identity],
                               [rejected, accepted.conjugate() * polars[vertex]]])
            local_word = pair_basis @ _embed_modes(3, (0, 2), scalar) @ pair_basis.conj().T
            gauges[2 * vertex] = gauges[vertex] * c / a
            gauges[2 * vertex + 1] = gauges[vertex] * s / b
            t = np.exp(2j * angle)
            coefficients[detail] = (z - np.cos(2 * angle)) / (2j * gauges[vertex] * c * s)
            corrections[detail] = -(z - t) * (z - t.conjugate()) / (
                2 * np.sqrt(2) * gauges[vertex] * c * s)
        else:
            coarse_angle = angle + eta * s / 20
            a, b = np.cos(coarse_angle), np.sin(coarse_angle)
            loss = 1 / (1200 * height)
            carrier = np.sqrt(1 - loss) * np.exp(.61j)
            split = np.sqrt(loss) * np.exp(-.37j)
            scalar = np.kron(np.array([[carrier, -split.conjugate()],
                                      [split, carrier.conjugate()]]), identity)
            real_gate = np.array([[a, -b], [b, a]])
            gate, feature = real_gate @ np.diag([carrier.conjugate(), 1]), np.array([a, b])
            accepted, rejected = carrier, split * identity
            local_word = _embed_modes(3, (0, 2), scalar) @ _embed_modes(
                3, (0, 1), np.kron(real_gate.T, identity))
            gauges[2 * vertex] = gauges[vertex] * c / (carrier.conjugate() * a)
            gauges[2 * vertex + 1] = gauges[vertex] * s / (carrier.conjugate() * b)
            coefficients[detail], corrections[detail] = splitter_gauge_coefficients(
                (c, s), (a, b), carrier, gauges[vertex])
            node.update(carrier=carrier, coarse_pair=(a, b))
        coarse[[left, detail]] = gate @ coarse[[left, detail]]
        node.update(gate=gate, feature=feature, loss=loss, accepted=accepted,
                    rejected=rejected, scalar=scalar, local_word=local_word)
        nodes[vertex] = node
    roots = {size + leaf: np.eye(size)[:, leaf] for leaf in range(size)}
    injections, row_word = np.zeros((size, size), complex), np.eye(2 * size * 2, dtype=complex)
    for vertex in (3, 2, 1):
        node = nodes[vertex]
        left, right = roots[2 * vertex], roots[2 * vertex + 1]
        roots[vertex] = node['gate'][0, 0] * left + node['gate'][1, 0] * right
        injections[:, node['detail']] = node['feature'][0] * left + node['feature'][1] * right
        row_word = _embed_modes(2 * size, (node['left'], node['detail'], size + node['detail']),
                                node['local_word']) @ row_word
    return dict(nodes=nodes, target=target, coarse=coarse, coefficients=coefficients,
                corrections=corrections, diagonal=np.array([gauges[size + j] for j in range(size)]),
                injections=injections, row_word=row_word)


def _mixed_source_words(fixture):
    """Full four-bank words using actual scalar replay, not completed columns."""
    size, identity, kappa, scale = 4, np.eye(2), 8 / 5, 5 / 4
    pword = np.eye(4 * size * 2, dtype=complex)
    detail_nodes = {node['detail']: node for node in fixture['nodes'].values()}
    for label in range(size):
        node = detail_nodes.get(label)
        x = fixture['coefficients'][label] / kappa
        y = (fixture['corrections'][label] / (kappa * node['loss'])
             if node is not None and node['loss'] else 0j)
        loss = node['loss'] if node is not None else 0.
        radial = np.sqrt(abs(x) ** 2 + loss * abs(y) ** 2)
        attenuation = np.array([[radial, -np.sqrt(1 - radial ** 2)],
                                [np.sqrt(1 - radial ** 2), radial]])
        pword = _embed_modes(4 * size, (label, 2 * size + label),
                             np.kron(attenuation, identity)) @ pword
        if y == 0:
            replay = np.kron(np.diag([x / radial, 1]), identity)
        else:
            accepted, scalar = node['accepted'], node['scalar']
            phase = np.exp(2j * np.arcsin(abs(y / radial) / (2 * abs(accepted))))
            first = 1 + (phase - 1) * abs(accepted) ** 2
            second = (phase - 1) * accepted.conjugate()
            phases = np.exp(1j * np.angle([x / (radial * first), y / (radial * second)]))
            replay = (np.kron(np.diag(phases), identity) @ scalar @
                      np.kron(np.diag([phase, 1]), identity) @ scalar.conj().T)
        pword = _embed_modes(4 * size, (label, size + label), replay) @ pword
    qword = np.eye(4 * size * 2, dtype=complex)
    for label, diagonal in enumerate(fixture['diagonal']):
        value = diagonal.conjugate() / scale
        reject = np.sqrt(1 - abs(value) ** 2)
        attenuation = np.array([[value, -reject], [reject, value.conjugate()]])
        qword = _embed_modes(4 * size, (label, 3 * size + label),
                             np.kron(attenuation, identity)) @ qword
    qword = _direct_sum((fixture['row_word'], np.eye(2 * size * 2))) @ qword
    return pword, qword


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

    def test_all_height_clipped_ledger_without_logical_padding(self):
        # Tiny heights use the finite fallback. Rounding a scalar loop order
        # never changes the logical height, including immediately past powers
        # of two. Huge heights exercise exponent-only accounting.
        heights = tuple(range(1, 97)) + ((1 << 20) - 1, (1 << 20) + 1,
                                       (1 << 80) - 3, (1 << 80) + 1)
        for height in heights:
            with self.subTest(height=height):
                parameters = clipped_source_parameters(height, suffix_constant=3)
                order = parameters.order
                logarithm = order.bit_length() - 1
                self.assertEqual(parameters.height, height)
                self.assertLessEqual(height, order)
                self.assertLess(order, 2 * height)
                self.assertEqual(order & (order - 1), 0)
                self.assertEqual(parameters.rare_probability * 8 * order, 1)
                self.assertEqual(parameters.source_length, 4 * order ** 2)
                self.assertEqual(parameters.rescaling_bits, 6 + 2 * logarithm)
                self.assertEqual(parameters.spacing_log2 + 2 * logarithm, height)
                self.assertEqual(parameters.minimum_expensive_exponent_log2 + logarithm,
                                 height + 1)
                self.assertEqual(parameters.suffix_bits,
                                 min(height, 3 * (order + 1).bit_length()))
                self.assertEqual(parameters.clipping_bound_relative_to_endpoint,
                                 Fraction(height, 16 * order))
                self.assertLessEqual(parameters.clipping_bound_relative_to_endpoint,
                                     Fraction(1, 16))
                if height < 14:
                    self.assertTrue(parameters.uses_finite_fallback)
                    self.assertFalse(parameters.stage_precision_at_most_twice_endpoint(20))
                    with self.assertRaises(ValueError):
                        parameters.coefficient_precision_overhead(20)
                    continue
                self.assertFalse(parameters.uses_finite_fallback)
                mantissa_bits = (16000 * order ** 2 - 1).bit_length()
                overhead = parameters.coefficient_precision_overhead(mantissa_bits)
                self.assertEqual(overhead,
                                 order * (mantissa_bits + order ** 2 * parameters.rescaling_bits))
                if height <= 96:
                    # N is small here. The separate huge-height cases above
                    # never construct N, much less a denominator 2**N.
                    endpoint = 1 << height
                    spacing = 1 << parameters.spacing_log2
                    cap = endpoint + parameters.clipping_exponent_offset - 1
                    self.assertEqual(parameters.source_length * spacing, 4 * endpoint)
                    self.assertEqual(cap // spacing, parameters.maximum_multiplier)
                    self.assertTrue(parameters.retained_exponent_within_cap(cap))
                    self.assertFalse(parameters.retained_exponent_within_cap(cap + 1))
                    self.assertGreaterEqual(parameters.source_length *
                                            (spacing - parameters.rescaling_bits), 3 * endpoint)
                    if parameters.stage_precision_at_most_twice_endpoint(mantissa_bits):
                        self.assertLessEqual(height * (order * spacing + overhead), 2 * endpoint)
                else:
                    self.assertTrue(parameters.stage_precision_at_most_twice_endpoint(mantissa_bits))
                    self.assertTrue(parameters.retained_exponent_within_cap(height ** 3))

    def test_certified_clipping_controls_complete_frame_by_layer_maximum(self):
        for height in (2, 3, 4):
            size = 1 << height
            parameters = clipped_source_parameters(height, suffix_constant=3)
            exponent = size + parameters.clipping_exponent_offset
            threshold = Fraction(1, 1 << exponent)
            original, clipped = np.eye(size), np.eye(size)
            layer_bounds, retained_depths = [], []
            for depth in range(height):
                width = size >> depth
                original_layer, clipped_layer = np.eye(size), np.eye(size)
                local_errors, retained = [], []
                for node, left in enumerate(range(0, size, width)):
                    right = left + width // 2
                    if depth == 0:
                        parameter, offset, supplied_bin = Fraction(1, 4), Fraction(), 1
                    elif depth == height - 1:
                        parameter, offset, supplied_bin = threshold / 4, -threshold / 8, None
                    else:
                        parameter = 3 * threshold / 4
                        offset = threshold / 8 if node % 2 == 0 else -threshold / 8
                        supplied_bin = exponent - 1
                    cosine, sine = cayley_pair(parameter)
                    approximation = sine + offset
                    self.assertLessEqual(abs(approximation - sine), threshold / 8)
                    drop = clipped_branch_decision(approximation, threshold)
                    local = np.array([[float(cosine), -float(sine)],
                                      [float(sine), float(cosine)]])
                    original_layer[np.ix_((left, right), (left, right))] = local
                    if drop:
                        self.assertLessEqual(sine, 13 * threshold / 8)
                        local_errors.append(np.linalg.norm(local - np.eye(2), 2))
                    else:
                        self.assertGreater(sine, threshold)
                        self.assertIsNotNone(supplied_bin)
                        self.assertLess(Fraction(1, 1 << (supplied_bin + 2)), sine)
                        self.assertLess(sine, Fraction(1, 1 << supplied_bin))
                        self.assertTrue(parameters.retained_exponent_within_cap(supplied_bin))
                        retained.append(supplied_bin)
                        clipped_layer[np.ix_((left, right), (left, right))] = local
                        local_errors.append(0.)
                maximum = max(local_errors)
                self.assertAlmostEqual(np.linalg.norm(original_layer - clipped_layer, 2),
                                       maximum, delta=ATOL)
                if depth == height - 1:
                    # Several disjoint tiny rotations cost their maximum in
                    # operator norm. Summing nodes would overcharge this layer.
                    self.assertGreater(sum(local_errors), 1.5 * maximum)
                original = original_layer @ original
                clipped = clipped_layer @ clipped
                layer_bounds.append(maximum)
                retained_depths.append(tuple(retained))
            for frame in (original, clipped):
                np.testing.assert_allclose(frame.T @ frame, np.eye(size), atol=ATOL, rtol=0)
            error = np.linalg.norm(original - clipped, 2)
            self.assertGreater(error, float(threshold / 4))
            self.assertLessEqual(error, sum(layer_bounds) + ATOL)
            self.assertLessEqual(sum(layer_bounds), float(4 * height * threshold))
            self.assertLessEqual(4 * height * threshold, Fraction(1, 16 * (1 << size)))
            self.assertEqual(bottom_precision_trace(retained_depths, 1), 0)

        # Certified comparisons can deliberately overlap: the same branch
        # admits either decision, while an exact dyadic branch needs no exact
        # logarithm. These are supplied open-width-four bins.
        threshold = Fraction(1, 1024)
        branch = 3 * threshold / 2
        self.assertTrue(clipped_branch_decision(branch - threshold / 8, threshold))
        self.assertFalse(clipped_branch_decision(branch + threshold / 8, threshold))
        self.assertTrue(clipped_branch_decision(3 * threshold / 2, threshold))
        for power in (1, 4, 9):
            branch, supplied_bin = Fraction(1, 1 << power), power - 1
            self.assertLess(Fraction(1, 1 << (supplied_bin + 2)), branch)
            self.assertLess(branch, Fraction(1, 1 << supplied_bin))

    def test_clipped_bottom_trace_counts_depths_and_empty_layers(self):
        depths = ((60,), (19, 19), (), (4, 4, 4), (7, 1, 7), ())
        parameters = clipped_source_parameters(len(depths), suffix_constant=3)
        self.assertTrue(all(parameters.retained_exponent_within_cap(exponent)
                            for depth in depths for exponent in depth))
        self.assertEqual(bottom_precision_trace(depths, 4), 11)
        self.assertEqual(bottom_precision_trace(depths, 0), 0)
        self.assertEqual(bottom_precision_trace(depths, 6), 90)
        repeated = tuple(depth * 17 for depth in depths)
        self.assertEqual(bottom_precision_trace(repeated, 4), 11)
        self.assertNotEqual(sum(sum(depth) for depth in repeated[-4:]), 11)
        self.assertNotEqual(4 * max(max(depth, default=0) for depth in depths), 11)

    def test_all_height_source_remainders_reach_clipped_cap_and_saturate(self):
        parameters = clipped_source_parameters(14, suffix_constant=3)
        self.assertNotEqual(parameters.height, parameters.order)
        spacing = 1 << parameters.spacing_log2
        beta = common_source_ratio(spacing, parameters.rescaling_bits)
        endpoint = 1 << parameters.height
        cap = endpoint + parameters.clipping_exponent_offset - 1
        order, probability = parameters.order, parameters.rare_probability
        mantissa = Fraction(3, 4)
        for exponent in (2 * endpoint // order + 1, endpoint // 4 + 3, cap):
            multiplier, remainder = exponent_metadata(exponent, spacing)
            self.assertGreater(remainder, 0)
            self.assertLess(remainder, spacing)
            self.assertLessEqual(multiplier, parameters.maximum_multiplier)
            coefficients = common_source_kernel_coefficients(
                mantissa, order, probability, parameters.rescaling_bits,
                multiplier=multiplier, remainder=remainder)
            prefactor = 2 * order * (1 - probability) / probability
            target_parameter = mantissa / (1 << exponent)
            for degree in (1, 3, order):
                coefficient = coefficients[degree - 1]
                shift = multiplier * degree
                self.assertEqual(coefficient * beta ** shift,
                                 prefactor * target_parameter ** degree)
                self.assertLessEqual(abs(coefficient), 1)
                self.assertLessEqual(dyadic_precision(coefficient), degree *
                                     (2 + remainder + multiplier * parameters.rescaling_bits))
                if shift >= parameters.source_length:
                    # S_shift=0, not a cyclic shift: its complete source
                    # residual is the entire desired beta**shift column.
                    self.assertEqual(source_shift_residual_squared(
                        beta, parameters.source_length, shift), beta ** (2 * shift))
            if exponent == cap:
                self.assertEqual(multiplier, parameters.maximum_multiplier)
                self.assertGreater(multiplier * order, parameters.source_length)

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
        # The earlier UNIFORM-ATTENUATION chart, not the mixed local chart
        # tested below: a connected tree, occupied features, and a single
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

        # In this uniform-attenuation chart the identity-merge splitter remains.
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

    def test_mixed_local_charts_full_tree_phase_and_dirty_row_word(self):
        identity = np.eye(2)
        for zero_left in (False, True):
            with self.subTest(zero_left=zero_left):
                fixture = _mixed_chart_fixture(zero_left)
                row_word, nodes = fixture['row_word'], fixture['nodes']
                if not zero_left:
                    first = nodes[1]['rejected'] / np.sqrt(nodes[1]['loss'])
                    second = nodes[2]['rejected'] / np.sqrt(nodes[2]['loss'])
                    self.assertGreater(np.linalg.norm(first @ second - second @ first, 2), 1.)
                coarse, injected = fixture['coarse'], fixture['injections']
                diagonal, coefficients, corrections = (fixture[key] for key in
                                                       ('diagonal', 'coefficients', 'corrections'))
                np.testing.assert_allclose(row_word.conj().T @ row_word, np.eye(16),
                                           atol=ATOL, rtol=0)
                np.testing.assert_allclose(row_word[:8, :8], np.kron(coarse.conj().T, identity),
                                           atol=ATOL, rtol=0)
                np.testing.assert_array_equal(row_word[8:10, :8], np.zeros((2, 8)))
                for node in nodes.values():
                    detail = node['detail']
                    expected = node['rejected'] @ np.kron(injected[:, detail].conj()[None, :], identity)
                    np.testing.assert_allclose(row_word[2 * (4 + detail):2 * (5 + detail), :8],
                                               expected, atol=ATOL, rtol=0)
                    np.testing.assert_allclose(np.eye(2) - node['gate'] @ node['gate'].conj().T,
                                               node['loss'] * np.outer(node['feature'],
                                                                       node['feature'].conj()),
                                               atol=ATOL, rtol=0)
                wanted = np.diag(diagonal) @ (coarse @ np.diag(coefficients) +
                                              injected @ np.diag(corrections))
                np.testing.assert_allclose(wanted, fixture['target'], atol=ATOL, rtol=0)
                self.assertTrue(np.all((abs(diagonal) > .999) & (abs(diagonal) < 1.002)))
                self.assertTrue(np.all((abs(coefficients) > .98) & (abs(coefficients) < 1.02)))
                for vertex, node in nodes.items():
                    if node['loss']:
                        bound = 4.8 / 2 if vertex == 3 else .8
                        self.assertLess(abs(corrections[node['detail']]) / node['loss'], bound)
                if zero_left:
                    zero = nodes[2]
                    np.testing.assert_array_equal(zero['local_word'], np.eye(6))
                    self.assertEqual(corrections[zero['detail']], 0)

                # The actual row target includes dirty action on occupied
                # feature columns. It cannot be replaced by a logical word
                # tensored with identity on the dirty register.
                dirty_z = np.kron(np.eye(8), np.diag([1, -1]))
                self.assertGreater(np.linalg.norm(row_word @ dirty_z - dirty_z @ row_word, 2), .5)
                source = np.kron(_source(Fraction(1, 4), 4)[:, None], np.eye(16))
                np.testing.assert_allclose(np.kron(np.eye(4), row_word) @ source,
                                           source @ row_word, atol=ATOL, rtol=0)

                expensive = nodes[3]
                wrong_c, wrong_a = splitter_gauge_coefficients(
                    expensive['target'], expensive['coarse_pair'], abs(expensive['carrier']),
                    expensive['gauge'])
                bad_c, bad_a = coefficients.copy(), corrections.copy()
                bad_c[expensive['detail']], bad_a[expensive['detail']] = wrong_c, wrong_a
                phase_dropped = np.diag(diagonal) @ (coarse @ np.diag(bad_c) + injected @ np.diag(bad_a))
                self.assertGreater(np.linalg.norm(phase_dropped - fixture['target'], 2), .5)

    def test_mixed_actual_replay_overlap_and_recovery_excludes_dirty_reflection(self):
        for zero_left in (False, True):
            with self.subTest(zero_left=zero_left):
                fixture = _mixed_chart_fixture(zero_left)
                pword, qword = _mixed_source_words(fixture)
                identity = np.eye(2)
                for word in (pword, qword):
                    np.testing.assert_allclose(word.conj().T @ word, np.eye(32), atol=ATOL, rtol=0)
                psource, qsource = pword[:, :8], qword[:, :8]
                np.testing.assert_allclose(psource[:8], np.kron(np.diag(fixture['coefficients'] / 1.6),
                                                               identity), atol=ATOL, rtol=0)
                wrong_source = psource.copy()
                for vertex, node in fixture['nodes'].items():
                    detail, loss = node['detail'], node['loss']
                    rows, columns = slice(2 * (4 + detail), 2 * (5 + detail)), slice(2 * detail, 2 * detail + 2)
                    expected = (fixture['corrections'][detail] * node['rejected'] / (1.6 * loss)
                                if loss else np.zeros((2, 2)))
                    np.testing.assert_allclose(psource[rows, columns], expected, atol=ATOL, rtol=0)
                    if loss and vertex != 3:
                        wrong_source[rows, columns] = (fixture['corrections'][detail] /
                                                       (1.6 * np.sqrt(loss))) * identity
                        scalar, accepted = node['scalar'], node['accepted']
                        phase = np.exp(.83j)
                        signal = np.kron(np.diag([phase, 1]), identity)
                        replayed = (scalar @ signal @ scalar.conj().T)[:, :2]
                        wanted = np.vstack(((1 + (phase - 1) * abs(accepted) ** 2) * identity,
                                            (phase - 1) * accepted.conjugate() * node['rejected']))
                        np.testing.assert_allclose(replayed, wanted, atol=ATOL, rtol=0)
                        wrong_order = (scalar.conj().T @ signal @ scalar)[:, :2]
                        self.assertGreater(np.linalg.norm(wrong_order - wanted, 2), 1e-4)
                expected = np.kron(fixture['target'], identity)
                np.testing.assert_allclose(qsource.conj().T @ psource, expected / 2, atol=ATOL, rtol=0)
                # Replacing dirty rejection by its norm preserves column
                # normalization but loses the required source overlap.
                np.testing.assert_allclose(wrong_source.conj().T @ wrong_source, np.eye(8),
                                           atol=ATOL, rtol=0)
                self.assertGreater(np.linalg.norm(qsource.conj().T @ wrong_source - expected / 2, 2), 1e-8)
                overlap = qword.conj().T @ pword
                reflection = np.diag([-1] * 8 + [1] * 24)
                recovered = (-overlap @ reflection @ overlap.conj().T @ reflection @ overlap)[:, :8]
                full_target = np.vstack((expected, np.zeros((24, 8))))
                np.testing.assert_allclose(recovered, full_target, atol=2 * ATOL, rtol=0)
                # This projector illegally demands dirty=0 as well as clean
                # flags=0. The correct projector includes every dirty column.
                wrong_reflection = np.eye(32)
                wrong_reflection[np.arange(0, 8, 2), np.arange(0, 8, 2)] = -1
                wrong = (-overlap @ wrong_reflection @ overlap.conj().T @
                         wrong_reflection @ overlap)[:, :8]
                self.assertGreater(np.linalg.norm(wrong - full_target, 2), 1.)

    def test_exact_integer_metadata_does_not_overflow_numpy_scalars(self):
        for bits in (10, 64):
            arguments = (Fraction(1, 2), Fraction(3, 4), 2, 4, Fraction(1, 64), bits)
            expected = common_source_kernel_error_bound(*arguments, multiplier=2, remainder=1)
            actual = common_source_kernel_error_bound(
                Fraction(np.int64(1), np.int64(2)), Fraction(np.int64(3), np.int64(4)),
                np.int64(2), np.int64(4), Fraction(np.int64(1), np.int64(64)), np.int64(bits),
                multiplier=np.int64(2), remainder=np.int64(1))
            self.assertEqual(actual, expected)
            self.assertGreater(actual.denominator, 0)
            self.assertEqual(kernel_source_error_bound(*arguments),
                             kernel_source_error_bound(*arguments[:5], np.int64(bits)))

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
        with self.assertRaises(ValueError):
            splitter_gauge_coefficients((1., 0.), (.8, .6), .99)


if __name__ == "__main__":
    unittest.main()
