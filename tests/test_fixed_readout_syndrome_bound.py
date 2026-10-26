"""Finite witnesses for fixed-readout syndrome covering, not a source compiler."""
import math
import unittest

import numpy as np

from tests.test_minimal_tree_probe import _tree_forward


def _rank(rows):
    pivots = {}
    for row in rows:
        while row:
            bit = row.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = row
                break
            row ^= pivots[bit]
    return len(pivots)


def _labels(p, degree):
    return [sum(1 << subset for subset in range(1 << p)
                if subset & x == subset and subset.bit_count() <= degree)
            for x in range(1 << p)]


def _syndrome(labels, basis, signs=0):
    return sum((((label & basis).bit_count() & 1) ^ ((signs >> x) & 1)) << x
               for x, label in enumerate(labels))


def _covering_radius(code, length):
    return max(min((word ^ candidate).bit_count() for candidate in code)
               for word in range(1 << length))


def _z_matrix(label, width):
    return np.diag([(-1) ** ((label & basis).bit_count() & 1)
                    for basis in range(1 << width)])


def _tree_query_label(x, a, b, p):
    m = 1 << p
    mask = (1 << m) - 1
    source_x = a & mask
    target_z = b >> m
    a ^= _tree_forward(x, source_x, p) << m
    transpose = sum((((_tree_forward(x, 1 << j, p) & target_z).bit_count() & 1) << j)
                    for j in range(m))
    return a, b ^ transpose


def _fixed_fourier_wall(a, b, p):
    m = 1 << p
    def reverse_bank_coordinates(value):
        return sum(((value >> (bank * m + subset)) & 1)
                   << (bank * m + (m - 1 - subset))
                   for bank in range(2) for subset in range(m))
    return reverse_bank_coordinates(b), reverse_bank_coordinates(a)


class FixedReadoutSyndromeBoundTests(unittest.TestCase):
    def test_complete_lookup_labels_realize_every_target_syndrome(self):
        for p in range(1, 7):
            m = 1 << p
            labels = _labels(p, p)
            self.assertEqual(_rank(labels), m)
            if p <= 3:
                self.assertEqual(
                    {_syndrome(labels, basis) for basis in range(1 << m)},
                    set(range(1 << m)),
                )

    def test_low_degree_rank_and_arbitrary_dirty_spectator_do_not_add_syndromes(self):
        for p in range(1, 4):
            m = 1 << p
            signs = sum((x.bit_count() & 1) << x for x in range(m)) ^ (1 << (m - 1))
            for degree in range(p + 1):
                labels = _labels(p, degree)
                expected = sum(math.comb(p, j) for j in range(degree + 1))
                self.assertEqual(_rank(labels), expected)
                # An invertible fixed Pauli-coordinate change mixes an arbitrary
                # spectator into the constant generator; no helper is zeroed.
                mixed = [label ^ ((label & 1) << m) for label in labels]
                code = {_syndrome(mixed, basis, signs)
                        for basis in range(1 << (m + 1))}
                self.assertEqual(_rank(mixed), expected)
                self.assertEqual(len(code), 1 << expected)

    def test_exact_signed_codes_have_the_claimed_covering_deficits(self):
        p, m = 3, 8
        for degree, expected_radius in enumerate((4, 2, 1, 0)):
            labels = _labels(p, degree)
            code = {_syndrome(labels, basis) for basis in range(1 << m)}
            shifted = {word ^ 0b10110100 for word in code}
            self.assertEqual(_covering_radius(code, m), expected_radius)
            self.assertEqual(_covering_radius(shifted, m), expected_radius)
            volume = sum(math.comb(m, j) for j in range(expected_radius + 1))
            self.assertGreaterEqual(len(code) * volume, 1 << m)
        # Rank-four first-degree syndromes cannot cover with radius one.
        self.assertLess((1 << 4) * (1 + 8), 1 << 8)

    def test_commuting_pauli_product_signs_give_an_affine_code(self):
        x = np.array([[0, 1], [1, 0]], dtype=complex)
        y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        z = np.diag([1, -1]).astype(complex)
        p0, p1 = np.kron(x, x), np.kron(z, z)
        family = (p0, p1, np.kron(y, y), -p0)
        np.testing.assert_allclose(p0 @ p1, -family[2], atol=1e-14)
        actual = set()
        for s in (-1, 1):
            for t in (-1, 1):
                projector = (np.eye(4) + s * p0) @ (np.eye(4) + t * p1) / 4
                self.assertAlmostEqual(np.trace(projector).real, 1)
                bits = [round((1 - np.trace(projector @ operator).real) / 2)
                        for operator in family]
                self.assertEqual(bits[2], bits[0] ^ bits[1] ^ 1)
                self.assertEqual(bits[3], bits[0] ^ 1)
                actual.add(sum(bit << i for i, bit in enumerate(bits)))
        self.assertEqual(len(actual), 4)
        self.assertNotIn(0, actual)
        offset = next(iter(actual))
        self.assertEqual(_rank(word ^ offset for word in actual), 2)

    def test_anticommuting_expectation_screen_has_the_displayed_threshold(self):
        x = np.array([[0., 1.], [1., 0.]])
        z = np.diag([1., -1.])
        h = (x + z) / math.sqrt(2)
        _, eigenvectors = np.linalg.eigh(h)
        state = eigenvectors[:, -1]
        expectations = [float(state @ operator @ state) for operator in (x, z)]
        self.assertAlmostEqual(sum(value ** 2 for value in expectations), 1)
        for operator in (x, z):
            self.assertAlmostEqual(np.linalg.norm((operator - np.eye(2)) @ state) ** 2,
                                   2 - math.sqrt(2))
        gamma = math.sqrt(2 - math.sqrt(2))
        self.assertGreater(2 * (1 - (gamma - 0.01) ** 2 / 2) ** 2, 1)

    def test_full_non_gaussian_readout_obeys_the_mismatch_identity(self):
        p, m, signs = 2, 4, 0b1010
        source_labels, target_labels = _labels(p, 1), _labels(p, p)
        h = np.array([[1., 1.], [1., -1.]]) / math.sqrt(2)
        fourier = np.array([[1.]])
        for _ in range(m):
            fourier = np.kron(fourier, h)
        # Fixed T and CCZ phase after a complete carrier Fourier wall.
        phase = np.array([np.exp(1j * math.pi / 4 *
                                ((basis & 1) + 4 * ((basis & 7) == 7)))
                          for basis in range(1 << m)])
        unitary = phase[:, None] * fourier
        source = [(-1) ** ((signs >> x) & 1) * _z_matrix(label, m)
                  for x, label in enumerate(source_labels)]
        target = [_z_matrix(label, m) for label in target_labels]
        images = [unitary @ operator @ unitary.conj().T for operator in source]
        errors = [np.linalg.norm(image - desired, ord=2)
                  for image, desired in zip(images, target)]
        code = {_syndrome(source_labels, basis, signs) for basis in range(1 << m)}
        for basis in range(1 << m):
            desired_word = _syndrome(target_labels, basis)
            state = unitary.conj().T[:, basis]
            expected_distance = 0.
            for x, operator in enumerate(source):
                desired_sign = (-1) ** ((desired_word >> x) & 1)
                mismatch = sum(abs(state[j]) ** 2 for j in range(1 << m)
                               if operator[j, j] != desired_sign)
                residual = np.linalg.norm((operator - desired_sign * np.eye(1 << m)) @ state) ** 2
                self.assertAlmostEqual(residual, 4 * mismatch)
                self.assertLessEqual(residual, errors[x] ** 2 + 1e-12)
                expected_distance += mismatch
            distance = min((desired_word ^ word).bit_count() for word in code)
            self.assertLessEqual(distance, expected_distance + 1e-12)

    def test_actual_tree_xor_queries_have_bounded_label_degree_with_fourier_walls(self):
        p, m = 3, 8
        witnessed_degree_two = False
        for stages in range(4):
            values = []
            for x in range(1 << p):
                a, b = 0, 1 << m
                for _ in range(stages):
                    a, b = _tree_query_label(x, a, b, p)
                    a, b = _fixed_fourier_wall(a, b, p)
                values.append(a | (b << (2 * m)))
            coefficients = values.copy()
            for i in range(p):
                for x in range(1 << p):
                    if x >> i & 1:
                        coefficients[x] ^= coefficients[x ^ (1 << i)]
            degree = max((x.bit_count() for x, coefficient in enumerate(coefficients)
                          if coefficient), default=0)
            self.assertLessEqual(degree, stages)
            self.assertLessEqual(_rank(values), sum(math.comb(p, j)
                                                  for j in range(min(stages, p) + 1)))
            witnessed_degree_two |= degree >= 2
        self.assertTrue(witnessed_degree_two)


if __name__ == '__main__':
    unittest.main()
