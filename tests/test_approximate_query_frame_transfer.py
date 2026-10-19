"""One exact native fixture for robust pair-sector Majorana geometry.

The perturbation mixes arbitrary helper inputs. Its full operator error
is bounded from an exact squared-error identity, without numerical norm
estimates or an inference from exact denominator bounds.
"""
from fractions import Fraction
from itertools import combinations
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, Q2, add, adjoint, eye, native_columns, scale, zeros,
)
from tests.test_majorana_walsh_query import (
    _adjacent_swap, _controlled_pauli, _inverse, _majoranas, _pauli_word,
    _shift, _walsh_schedule, _word,
)
from tests.test_two_layer_galois_return import _conjugate, _product, _tau


def _commutator(left, right):
    # Chronological word for left * right * left^dagger * right^dagger.
    return _inverse(right) + _inverse(left) + right + left


def _overlap(left, right):
    product = _product(left, right)
    return sum((product[j][j] for j in range(len(product))), CZ) / len(product)


def _majorana_action(unitary, width):
    gamma = [_tau(label, width) for label in _majoranas(width)]
    result = zeros(2 * width)
    for column, pauli in enumerate(gamma):
        image = _conjugate(unitary, pauli)
        reconstructed = zeros(1 << width)
        for row, output in enumerate(gamma):
            coefficient = _overlap(output, image)
            assert coefficient.imag == 0
            result[row][column] = coefficient
            reconstructed = add(reconstructed, scale(output, coefficient))
        assert reconstructed == image
    return result


def _projector(column):
    return [[left * right.conj() for right in column] for left in column]


def _counts(word):
    return (sum(gate[0] in ('T', 'TDG') for gate in word),
            sum(gate[0] not in ('T', 'TDG') for gate in word))


class ApproximateQueryFrameTransferTests(unittest.TestCase):
    def test_native_helper_mixing_preserves_robust_pair_planes_below_one_third(self):
        # Four original modes and two modes from one arbitrary helper qubit.
        width, m = 3, 4
        route = [gate for position in (3, 2, 1) for gate in _adjacent_swap(position)]
        first = route + [('T', 0)] + _inverse(route)  # plane {gamma_0,gamma_4}
        second = [('T', 2)]                         # plane {gamma_4,gamma_5}
        perturbation_word = _commutator(first, second)
        for _ in range(2):
            perturbation_word = _commutator(perturbation_word, first)
        walsh_word = _word(_walsh_schedule(2))
        source_word = walsh_word + perturbation_word
        perturbation = native_columns(width, perturbation_word)
        target = native_columns(width, walsh_word)
        source = native_columns(width, source_word)
        self.assertEqual(source, _product(perturbation, target))
        self.assertEqual(native_columns(width, _inverse(source_word)), adjoint(source))
        self.assertEqual(_counts(perturbation_word), (22, 504))
        for orientation in (source_word, _inverse(source_word)):
            self.assertEqual(_counts(orientation), (26, 534))
            self.assertTrue(all(gate[0] in ('H', 'S', 'CX', 'T', 'TDG')
                                for gate in orientation))

        # This scalar identity gives ||E-I||^2 exactly. Conjugating any
        # unitary sector by E changes it by at most tau=2||E-I||.
        error = add(perturbation, scale(eye(1 << width), -1))
        c = Q2(Fraction(23, 16), -1)
        self.assertEqual(_product(adjoint(error), error), scale(eye(1 << width), c))
        tau_squared = 4 * c
        self.assertGreater(c, 0)
        self.assertLess(tau_squared, Fraction(1, 9))
        delta_squared = 9 * c  # (3 tau / 2)^2, strictly below 1/4.
        self.assertLess(delta_squared, Fraction(1, 4))

        actual_modes = _majorana_action(source, width)
        target_modes = _majorana_action(target, width)
        self.assertEqual(_product(adjoint(actual_modes), actual_modes), eye(2 * width))
        expected_target = [[CO * Fraction((-1) ** ((row & column).bit_count()), 2)
                            if row < m and column < m else
                            (CO if row == column and row >= m else CZ)
                            for column in range(2 * width)] for row in range(2 * width)]
        self.assertEqual(target_modes, expected_target)

        labels = _majoranas(2)
        total = 0
        for label in labels:
            total ^= label
        dual = [total ^ label for label in labels]
        controlled = [gate for j in reversed(range(m))
                      for gate in _controlled_pauli(j, dual[j], 2, m)]
        readout_word = (_shift(_inverse(source_word), m) + controlled
                        + _shift(source_word, m))
        for orientation in (readout_word, _inverse(readout_word)):
            self.assertEqual(_counts(orientation), (52, 1086))
            self.assertLess(max(wire for gate in orientation for wire in gate[1:]), m + width)

        pair_planes = {}
        for left, right in combinations(range(m), 2):
            # Actual positive-Pauli representatives retain the pair sector's
            # potentially nonreal scalar; the rightmost factor acts first.
            center_word = _pauli_word(dual[right], 2) + _pauli_word(dual[left], 2)
            center = native_columns(width, center_word)
            center_inverse = native_columns(width, _inverse(center_word))
            self.assertEqual(center_inverse, adjoint(center))
            self.assertEqual(center_inverse, scale(center, -1))
            actual_sector = _conjugate(source, center)
            actual_sector_inverse = _conjugate(source, center_inverse)
            self.assertEqual(actual_sector_inverse, adjoint(actual_sector))
            target_sector = _conjugate(target, center)
            phase = (CO, CI, -CO, -CI)[(left + right) % 4]
            gauged_sector = scale(actual_sector, phase)
            self.assertEqual(scale(actual_sector_inverse, phase.conj()), adjoint(gauged_sector))
            sector_error = add(gauged_sector, scale(target_sector, -phase))
            squared_error = _product(adjoint(sector_error), sector_error)
            # Here the two planes share a line, so the sector's squared
            # error is scalar. This verifies the full operator bound on
            # every occupied and unoccupied helper column.
            sector_squared_norm = squared_error[0][0]
            self.assertEqual(squared_error, scale(eye(1 << width), sector_squared_norm))
            self.assertEqual(sector_squared_norm.imag, 0)
            self.assertGreater(sector_squared_norm.real, 0)
            self.assertLessEqual(sector_squared_norm.real, tau_squared)
            reflection = _majorana_action(gauged_sector, width)
            plane = scale(add(eye(2 * width), scale(reflection, -1)), Fraction(1, 2))
            expected_plane = add(_projector([row[left] for row in actual_modes]),
                                 _projector([row[right] for row in actual_modes]))
            self.assertEqual(plane, expected_plane)
            pair_planes[left, right] = plane

        for index in range(m):
            first_other, second_other = [j for j in range(m) if j != index][:2]
            plane_one = pair_planes[tuple(sorted((index, first_other)))]
            plane_two = pair_planes[tuple(sorted((index, second_other)))]
            plane_three = pair_planes[tuple(sorted((first_other, second_other)))]
            recovered_line = scale(add(add(plane_one, plane_two), scale(plane_three, -1)),
                                   Fraction(1, 2))
            column = [row[index] for row in actual_modes]
            target_column = [row[index] for row in target_modes]
            self.assertEqual(recovered_line, _projector(column))
            inner = sum((left.conj() * right for left, right in zip(column, target_column)), CZ)
            sine_squared = 1 - inner.abs2()
            self.assertGreater(sine_squared, 0)
            self.assertLessEqual(sine_squared, delta_squared)

            probabilities = [entry.abs2() for entry in column]
            reference = [Q2(Fraction(1, m)) if j < m else Q2() for j in range(2 * width)]
            variation = sum((abs(value - reference[j]) for j, value in enumerate(probabilities)), Q2()) / 2
            self.assertLessEqual(variation * variation, sine_squared)
            self.assertEqual(sum(probabilities, Q2()), 1)
            helper_mass = sum(probabilities[m:], Q2())
            self.assertGreater(helper_mass, 0)
            # Coarse-graining helpers splits one original probability 1/4
            # into two positive pieces, leaving the other three unchanged.
            # Thus H(coarse)=2+h_2(4*p_0)/4 >= 2 exactly, without evaluating logs.
            self.assertEqual(probabilities[1:m], [Q2(Fraction(1, m))] * (m - 1))
            self.assertEqual(probabilities[0] + helper_mass, Fraction(1, m))
            self.assertGreater(probabilities[0], 0)
            self.assertEqual(variation, helper_mass)
            support = sum(value != 0 for value in probabilities)
            self.assertGreater(support, m)
            original_support = sum(value != 0 for value in probabilities[:m])
            self.assertLessEqual(inner.abs2(), Fraction(original_support, m) * (1 - helper_mass))
            # The robust flat-overlap bound requires support greater than
            # m/2, so an integer matching-layer depth is still at least p=2.
            self.assertGreater(m * (1 - delta_squared), m // 2)


if __name__ == '__main__':
    unittest.main()
