"""Finite fixtures for scalar gaps, a rare frame cycle, and interface bounds."""
from itertools import combinations, product
from math import comb, cos, gcd, inf, pi, sin, sqrt
import unittest

import numpy as np

from compiler_robust_hopf.frames import complex_magnitude_frame_matrix, real_frame_matrix
from tests.test_four_layer_crossing_rank import _binary_basis, _weighted_phase_fixture
from tests.test_dirty_indicator_depth import (
    _ccz_two_layers, _exact_phase_action, _parity_basis,
)
from tests.test_one_band_walsh_blocks import _affine_spaces
from tests.test_operator_source_compiler import I2, X, Y, Z, _pauli, _word_matrix
from tests.test_outer_axis_normalizer import _tableau
from tests.test_t_depth import _inverse, _word


def _valuation(vector):
    """Exponent of 1-zeta in the integral basis 1,zeta,zeta²,zeta³."""
    if not any(vector):
        return inf
    order = 0
    while sum(vector) % 2 == 0:
        first = (vector[0] - sum(vector[1:])) // 2
        vector = [first, first + vector[1], first + sum(vector[1:3]),
                  first + sum(vector[1:])]
        order += 1
    return order


def _mobius_table(values):
    coefficients = values.copy()
    for bit in range((len(values) - 1).bit_length()):
        step = 1 << bit
        for offset in range(0, len(values), 2 * step):
            coefficients[offset + step:offset + 2 * step] -= coefficients[offset:offset + step]
    return coefficients


class NativeGapRefinementTests(unittest.TestCase):
    def test_complete_hopf_suffix_gives_the_rare_signed_three_cycle(self):
        prefix_bits = 3
        dimension = 1 << (prefix_bits + 2)
        for marked_prefix in (0, 5, 7):
            angles = np.zeros(dimension - 1)
            penultimate_node = (1 << prefix_bits) + marked_prefix
            deepest_node = (1 << (prefix_bits + 1)) + 2 * marked_prefix
            angles[[penultimate_node - 1, deepest_node - 1]] = pi / 2

            expected = np.eye(dimension)
            base = 4 * marked_prefix
            expected[base:base + 4, base:base + 4] = np.eye(4)[:, [2, 0, 1, 3]]
            expected[:, [base + 1, base + 2]] *= -1
            real = real_frame_matrix(angles)
            complex_frame = complex_magnitude_frame_matrix(angles, np.zeros(dimension))
            np.testing.assert_allclose(real, expected, atol=2e-12, rtol=0)
            np.testing.assert_allclose(complex_frame, expected, atol=2e-12, rtol=0)
            np.testing.assert_array_equal(expected @ expected @ expected, np.eye(dimension))
            self.assertAlmostEqual(np.linalg.norm(real - np.eye(dimension), "fro")**2 / dimension,
                                   1.5 / (1 << prefix_bits))

    def test_affine_support_incidence_and_scalar_threshold(self):
        spaces = [set(points) for _, _, points in _affine_spaces(3)]
        for support in spaces:
            near = [other for other in spaces
                    if 4 * len(support & other) ** 2 > len(support) * len(other)]
            for other in near:
                self.assertTrue(support <= other or other <= support)
                self.assertIn(len(other) / len(support), (0.5, 1, 2))
            for first, second in combinations(near, 2):
                if len(first) == len(second) == 2 * len(support):
                    self.assertFalse((first - support) & (second - support))
                if 2 * len(first) == len(support) == 2 * len(second) and first | second != support:
                    self.assertGreaterEqual(4 * len(support - (first | second)), len(support))
        bound = 2 * sin(pi / 24)
        block_loss = (1 / sqrt(2) - 0.5)**2 + (1 / sqrt(2) - sqrt(3) / 2)**2
        self.assertAlmostEqual(bound**2, block_loss)
        self.assertAlmostEqual(bound**2 / 16, (1 - cos(pi / 12)) / 8)
        self.assertLess(8 / 128, bound**2)
        self.assertGreater(8 / 64, bound**2)

    def test_native_sparse_probe_has_two_orthogonal_clifford_terms(self):
        phase = np.diag([1, np.exp(1j*pi/4)])
        band = np.kron(phase, phase.conj())
        hadamard = (X + Z) / sqrt(2)
        cnot = np.eye(4)[[0, 1, 3, 2]]
        core = cnot @ np.kron(hadamard, Z) @ cnot
        image = band @ core @ band.conj().T
        axes = (I2, X, Z, Y)
        paulis = {}
        for x, z in product(range(4), repeat=2):
            factors = {w: axes[((x >> w) & 1) + 2 * ((z >> w) & 1)] for w in range(2)}
            paulis[x, z] = _pauli(2, factors)
        for shift in (0, 1, 2):
            for mask in range(4):
                probe = paulis[shift, mask]
                pulled = band.conj().T @ probe @ band
                coefficients = [np.trace(p @ pulled)/4 for p in paulis.values()]
                terms = [(c, band @ core @ p @ core @ band.conj().T)
                         for c, p in zip(coefficients, paulis.values()) if abs(c) > 1e-10]
                self.assertEqual(len(terms), 1 + shift.bit_count())
                np.testing.assert_allclose(sum(c * m for c, m in terms),
                                           image @ probe @ image, atol=2e-12)
                for _, matrix in terms:
                    for p in paulis.values():
                        output = matrix @ p @ matrix.conj().T
                        self.assertAlmostEqual(max(abs(np.trace(output @ q)) for q in paulis.values()), 4)
                if len(terms) == 2:
                    self.assertAlmostEqual(abs(np.trace(terms[0][1] @ terms[1][1])), 0)

    def test_exponent_digits_after_quadratic_carry(self):
        width = 7
        rng = np.random.default_rng(75103)
        points = range(1 << width)
        for _ in range(12):
            parities = rng.integers(1, 1 << width, size=12)
            quadratics = rng.choice([m for m in points if m.bit_count() <= 2], size=7)
            values = []
            for x in points:
                total = sum((int(mask) & x).bit_count() % 2 for mask in parities)
                for digit in range(3):
                    self.assertEqual((total >> digit) & 1, comb(total, 1 << digit) % 2)
                correction = sum((x & int(mask)) == mask for mask in quadratics) % 2
                values.append((total - 2 * correction) % 8)
            for digit in range(3):
                bits = (np.array(values) >> digit) & 1
                coefficients = _mobius_table(bits) % 2
                self.assertTrue(all(x.bit_count() <= 1 << digit
                                    for x in points if coefficients[x]))

    def test_four_counter_corners_expose_false_row_overlap(self):
        identity = np.eye(4)
        fourier_zero = np.ones(4) / 2
        projector = np.outer(fourier_zero, fourier_zero)
        shift = np.diag([1, 1, 1, -1])
        conjugated = shift.conj().T @ projector @ shift
        flag_a = (np.kron(projector, np.kron(X, I2))
                  + np.kron(identity - projector, identity))
        flag_b = (np.kron(conjugated, np.kron(I2, X))
                  + np.kron(identity - conjugated, identity))
        corner = np.kron(identity, np.diag([1, 1, 1, -1]))
        comparison = (flag_b @ corner @ flag_a @ corner
                      @ flag_b @ corner @ flag_a @ corner)

        overlap = np.linalg.norm(projector @ shift @ projector, 2)
        self.assertAlmostEqual(overlap, 0.5)
        monitored = np.kron(fourier_zero, [1, 0, 0, 0])
        self.assertAlmostEqual(np.vdot(monitored, comparison @ monitored).real,
                               1 - 2 * overlap**2)
        self.assertAlmostEqual(np.linalg.norm(comparison @ monitored - monitored),
                               2 * overlap)

        inactive_clock = np.array([1, -1, 0, 0]) / sqrt(2)
        inactive_columns = np.kron(inactive_clock[:, None], identity)
        np.testing.assert_allclose(comparison @ inactive_columns, inactive_columns,
                                   atol=2e-12, rtol=0)

    def test_cyclotomic_filtration_and_minimum_stratum(self):
        rng = np.random.default_rng(74025)
        roots = np.vstack((np.eye(4, dtype=int), -np.eye(4, dtype=int)))
        for width in range(2, 7):
            size = 1 << width
            marked = np.zeros((size, 4), dtype=int)
            pi_power = roots[0].copy()
            for _ in range(width):
                shifted = np.roll(pi_power, 1)
                shifted[0] *= -1
                pi_power = pi_power - shifted
            marked[-1] = pi_power
            self.assertEqual(_valuation(marked[-1]), width)
            for _ in range(6):
                samples = marked.copy()
                for _ in range(3):
                    phase = _weighted_phase_fixture(width, rng)
                    for x in range(size):
                        exp = sum(c for mask, c in phase.items() if (mask & x) == mask)
                        samples[x] += roots[exp % 8]
                coefficients = _mobius_table(samples)
                self.assertTrue(all(_valuation(v) >= mask.bit_count()
                                    for mask, v in enumerate(coefficients)))
                orders = [_valuation(v) for v in samples]
                lowest = min(orders)
                self.assertEqual(lowest, min(map(_valuation, coefficients)))
                self.assertGreaterEqual(orders.count(lowest) * (1 << lowest), size)

    def test_occupied_quadratic_support_and_outer_crossings(self):
        width = 4
        monomials = [mask for mask in range(1 << width) if mask.bit_count() <= 2]
        for coeffs in range(1, 1 << len(monomials)):
            terms = [m for i, m in enumerate(monomials) if (coeffs >> i) & 1]
            support = {x for x in range(1 << width)
                       if sum(x & m == m for m in terms) % 2}
            reached = {next(iter(support))}
            while True:
                more = {x for x in support if any((x ^ y).bit_count() <= 3 for y in reached)}
                if more == reached:
                    break
                reached = more
            self.assertEqual(reached, support)
        for rows in product(range(8), repeat=3):
            form = [sum((((x & y).bit_count() % 2) << j) for j, y in enumerate(rows))
                    for x in rows]
            self.assertLessEqual(len(_binary_basis(form)), len(_binary_basis(rows)))
        self.assertEqual(sum(6**i for i in range(4)), 259)
        self.assertEqual(sum(126**i for i in range(64)), (126**64 - 1) // 125)


class TwoBandSignedGapTests(unittest.TestCase):
    def test_full_native_normal_form_and_sharp_sparse_fraction(self):
        # This actual two-band CCZ returns its arbitrary helper. The native
        # Clifford suffix supplies a common omega phase, so a forward call
        # cannot silently replace the actual adjoint in the pair signal.
        width = 4
        schedule = _ccz_two_layers(0, 1, 2, 3)
        phase_word = [('S', 2), ('H', 2)] * 3
        c0 = _word_matrix(width, schedule[0][1])
        t1 = _word_matrix(width, schedule[1][1])
        c1 = _word_matrix(width, schedule[2][1] + schedule[3][1])
        t2 = _word_matrix(width, schedule[4][1])
        c2 = _word_matrix(width, schedule[5][1] + phase_word)
        source = _word_matrix(width, _word(schedule) + phase_word)
        self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 2)
        np.testing.assert_allclose(source, c2 @ t2 @ c1 @ t1 @ c0,
                                   atol=3e-13)
        signal = _pauli(width, {2: X})
        j = t1 @ c0
        a = t2.conj().T @ c2.conj().T @ signal @ c2 @ t2
        b = j @ signal @ j.conj().T
        k = c1.conj().T @ a @ c1 @ b
        word = source.conj().T @ signal @ source @ signal
        expected = np.diag([(-1) ** ((z & 1) * ((z >> 1) & 1))
                            for z in range(1 << width)])
        np.testing.assert_allclose(word, expected, atol=3e-13)
        np.testing.assert_allclose(word, j.conj().T @ k @ j, atol=3e-13)
        self.assertGreater(np.linalg.norm(source @ signal @ source @ signal
                                         - expected, 2), 1)
        # Test the spectral statement on every logical/helper column, not
        # only a selected clean sector or a trace of the intended table.
        symplectic, _ = _tableau(k, width)
        difference = (symplectic + np.eye(2 * width, dtype=int)) % 2
        rows = [sum(int(bit) << i for i, bit in enumerate(row))
                for row in difference]
        rank = len(_binary_basis(rows))
        tau = np.trace(k) / (1 << width)
        spectrum = np.linalg.eigvalsh((k + k.conj().T) / 2)
        epsilon = np.count_nonzero(spectrum < 0) / len(spectrum)
        gamma = min(abs(spectrum))
        self.assertEqual(rank, 2)
        self.assertAlmostEqual(abs(tau) ** 2, 2 ** -rank)
        self.assertEqual(epsilon, 0.25)
        self.assertAlmostEqual(gamma, 1)
        self.assertAlmostEqual(gamma, sqrt(epsilon) / (1 - sqrt(epsilon)))
        # The conservative Pauli-balanced rounding threshold really starts
        # at five marked address bits.
        clifford_gap = sqrt(2 - sqrt(2))
        self.assertLess(4 * sqrt(2 ** -5), clifford_gap)
        self.assertGreaterEqual(4 * sqrt(2 ** -4), clifford_gap)

    def test_next_band_pair_source_has_a_real_gap_outside_two_band_bound(self):
        # A one-band Pauli image realizes CNOT(c,y1) H(q); a Clifford bank
        # change turns that CNOT into the full fanout. Both U and its actual
        # inverse are paid below, giving a three-band raw source.
        width, m, q = 5, 3, 4
        basis = _parity_basis([0, 1, q], [5, 6, 7])
        pre = [('S', q)] * 3 + [('H', q)] + basis
        post = (list(reversed(basis)) + [('H', q), ('S', q), ('H', 1)]
                + [('CX', 1, wire) for wire in range(2, m + 1)])
        u_schedule = [('C', pre),
                      ('T', [('TDG', 0), ('TDG', 1), ('T', q)]),
                      ('C', post)]
        e_schedule = [('T', [('T', wire) for wire in range(1, m + 1)])]
        source_schedule = u_schedule + e_schedule + _inverse(u_schedule)
        self.assertEqual(sum(kind == 'T' for kind, _ in source_schedule), 3)
        u = _word_matrix(width, _word(u_schedule))
        e = _word_matrix(width, _word(e_schedule))
        c = _word_matrix(width, [('CX', 0, wire)
                                 for wire in range(1, m + 1)])
        signal = _pauli(width, {q: X})
        expected_image = c @ _pauli(width, {q: (X + Z) / sqrt(2)})
        np.testing.assert_allclose(u @ signal @ u.conj().T,
                                   expected_image, atol=4e-13)
        source = _word_matrix(width, _word(source_schedule))
        k = e.conj().T @ c @ e @ c
        word = source.conj().T @ signal @ source @ signal
        np.testing.assert_allclose(word, u.conj().T @ k @ u, atol=4e-13)
        spectrum = np.linalg.eigvalsh((word + word.conj().T) / 2)
        epsilon = np.count_nonzero(spectrum < 0) / len(spectrum)
        gamma = min(abs(spectrum))
        self.assertEqual(epsilon, 1 / 8)
        self.assertAlmostEqual(gamma, 1 / sqrt(2))
        self.assertGreater(gamma, sqrt(epsilon) / (1 - sqrt(epsilon)))
        self.assertAlmostEqual(np.trace(k).real / len(k),
                               0.5 + 2 ** (-m / 2 - 1))

    def test_gapped_next_band_trace_does_not_bound_effective_spectral_width(self):
        # Count the exact whole-space spectrum without allocating matrices
        # of dimension 2**(m+1). The small emitted word independently checks
        # the diagonal action underlying the binomial multiplicities.
        m = 9
        fanout = [('CX', 0, wire) for wire in range(1, m + 1)]
        word = (fanout + [('T', wire) for wire in range(1, m + 1)]
                + fanout + [('TDG', wire) for wire in range(1, m + 1)])
        for z in range(1 << (m + 1)):
            expected = ((z & 1) * (m - 2 * (z >> 1).bit_count())) % 8
            self.assertEqual(_exact_phase_action(z, word), (z, expected))
        for m in (9, 17):
            residues = [sum(comb(m, weight) for weight in range(m + 1)
                            if weight % 4 == residue) for residue in range(4)]
            high, low = 2 ** (m - 2), 2 ** ((m - 3) // 2)
            self.assertEqual(residues, [high + low, high + low,
                                        high - low, high - low])
            multiplicities = [2 ** m] + residues
            common = gcd(*multiplicities)
            self.assertEqual(common, low)
            effective_bits = m + 1 - (common.bit_length() - 1)
            self.assertEqual(effective_bits, (m + 5) // 2)
            dimension = 2 ** (m + 1)
            eigenvalues = np.exp(1j * pi * (1 - 2 * np.arange(4)) / 4)
            tau = (2 ** m + sum(count * value for count, value
                               in zip(residues, eigenvalues))) / dimension
            self.assertAlmostEqual(tau.real, 0.5 + 2 ** (-m / 2 - 1))
            self.assertAlmostEqual(tau.imag, 0)
            self.assertAlmostEqual(min(abs(eigenvalues.real)), 1 / sqrt(2))
            self.assertEqual((residues[2] + residues[3]) / dimension,
                             0.25 - 2 ** (-(m + 3) / 2))
