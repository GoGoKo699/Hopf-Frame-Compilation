"""Mechanism checks for a six-layer Pauli-comparison cube detector.

The low-degree baseline, native moving-probe signs and full-operator
comparison error are checked separately. These examples do not enumerate
returned six-layer queries or prove arbitrary-width bounds numerically.
"""
from __future__ import annotations

import itertools
import math
import random
import unittest

import numpy as np

from tests.test_baseline_word_obstruction import _coordinate_word, _word_matrix
from tests.test_five_layer_native_defect import (
    _anf, _degree, _evaluate, _multiply, _substitute, _xor,
)
from tests.test_one_band_axis_matching import _projective_distance
from tests.test_outer_axis_normalizer import _image, _tableau, _tau
from tests.test_read_only_three_layer_rigidity import H, S, T, X, Z, _cx, _wire


def _product(operators):
    result = np.eye(len(operators[0]), dtype=complex)
    for operator in operators:
        result = result @ operator
    return result


def _cube_addresses(base, directions):
    result = []
    for mask in range(1 << len(directions)):
        address = base
        for bit, direction in enumerate(directions):
            if mask >> bit & 1:
                address ^= direction
        result.append(address)
    return result


def _fixed_action(clifford, labels):
    width = len(labels) // 2
    matrix, _ = _tableau(clifford, width)
    signs = _anf([_image(clifford, probe, width)[1]
                  for probe in range(1 << (2 * width))])
    outgoing = [_xor(*(labels[column] for column in range(2 * width)
                       if matrix[row, column])) for row in range(2 * width)]
    return outgoing, _substitute(signs, labels)


def _one_band_image_action(coordinate, family, labels, orientation):
    labels, sign = _fixed_action(coordinate.conj().T, labels)
    width = len(labels) // 2
    outgoing = list(labels)
    for wire in range(width):
        a, b = labels[wire], labels[width + wire]
        u, w = family[wire], family[width + wire]
        local = _xor(_multiply(_multiply(u, b), _xor({0}, a)),
                     _multiply(w, a),
                     _multiply(u, a) if orientation[wire] else set())
        sign = _xor(sign, local)
        outgoing[width + wire] = _xor(b, _multiply(u, a))
    outgoing, last_sign = _fixed_action(coordinate, outgoing)
    return outgoing, _xor(sign, last_sign)


class SixLayerPauliComparison(unittest.TestCase):
    def test_matched_native_baseline_and_literal_eight_point_cube(self):
        diagonal = _wire(T, 0) @ _wire(T.conj().T, 1)
        q = _wire(Z, 0)
        baseline_labels = []
        comparisons = []
        targets = []
        for address in range(16):
            a, b, c, d = [(address >> bit) & 1 for bit in range(4)]
            left_label = a | (b << 1) | (c << 2) | (d << 3)
            right_label = left_label ^ (a << 2) ^ (b << 3)
            pauli_label = c | (d << 1) | (b << 2) | (a << 3)
            left = diagonal @ _tau(left_label, 2) @ diagonal.conj().T
            right = diagonal @ _tau(right_label, 2) @ diagonal.conj().T
            baseline = left.conj().T @ _tau(pauli_label, 2) @ right
            matrix, signs = _tableau(baseline, 2)
            np.testing.assert_array_equal(matrix, np.eye(4, dtype=int))
            # Pauli generator signs recover the binary Pauli label.
            label = signs[2] | (signs[3] << 1) | (signs[0] << 2) | (signs[1] << 3)
            baseline_labels.append(label)
            target = a * b * c * d
            targets.append(target)
            comparisons.append(baseline @ (q if target else np.eye(4)))
        self.assertEqual(max(_degree(_anf([(label >> bit) & 1 for label in baseline_labels]))
                             for bit in range(4)), 2)
        seen = set()
        for address in range(16):
            cube = _cube_addresses(address, (1, 2, 4))
            derivative = sum(targets[value] for value in cube) % 2
            seen.add(derivative)
            self.assertEqual(np.bitwise_xor.reduce([baseline_labels[value] for value in cube]), 0)
            actual = _product([comparisons[value] for value in cube])
            expected = q if derivative else np.eye(4)
            self.assertLess(_projective_distance(actual, expected), 8e-12)
            phase = np.trace(expected.conj().T @ actual) / 4
            np.testing.assert_allclose(actual, phase * expected, atol=8e-14)
        self.assertEqual(seen, {0, 1})

    def test_two_one_band_factors_on_a_degree_fourteen_moving_probe(self):
        left_coordinate = _cx(0, 1) @ _wire(H, 1) @ _wire(S, 0)
        right_coordinate = _wire(H, 0) @ _cx(1, 0) @ _wire(S, 1)
        family_left = [{1 << bit} for bit in range(4)]
        family_right = [{1 << (4 + bit)} for bit in range(4)]
        middle = [_xor(family_left[bit], family_right[bit]) for bit in range(4)]
        incoming = [{sum(1 << (8 + 14 * bit + offset) for offset in range(14))}
                    for bit in range(4)]
        labels, sign = _one_band_image_action(right_coordinate, family_right,
                                             incoming, (False, True))
        self.assertLessEqual(max(map(_degree, labels)), 15)
        for wire in range(2):
            sign = _xor(sign, _multiply(middle[wire], labels[2 + wire]),
                        _multiply(middle[2 + wire], labels[wire]))
        labels, next_sign = _one_band_image_action(left_coordinate, family_left,
                                                  labels, (True, False))
        sign = _xor(sign, next_sign)
        self.assertLessEqual(max(map(_degree, labels)), 16)
        self.assertLessEqual(_degree(sign), 48)
        # This concrete native example exceeds the naive 3d sign estimate:
        # the variable Pauli-family labels contribute additional degree.
        self.assertGreater(_degree(sign), 3 * 14)

        left_layer = _wire(T.conj().T, 0) @ _wire(T, 1)
        right_layer = _wire(T, 0) @ _wire(T.conj().T, 1)
        # Substitutions include dense assignments, so degree14 incoming
        # monomials are tested both as zero and as one.
        rng = random.Random(63051)
        assignments = [0, (1 << 64) - 1]
        assignments += [((1 << 64) - 1) ^ (1 << bit) for bit in range(64)]
        assignments += [rng.getrandbits(64) for _ in range(32)]
        for assignment in assignments:
            left_label = sum(_evaluate(poly, assignment) << bit for bit, poly in enumerate(family_left))
            right_label = sum(_evaluate(poly, assignment) << bit for bit, poly in enumerate(family_right))
            middle_label = sum(_evaluate(poly, assignment) << bit for bit, poly in enumerate(middle))
            probe = sum(_evaluate(poly, assignment) << bit for bit, poly in enumerate(incoming))
            left = left_coordinate @ left_layer @ _tau(left_label, 2) @ left_layer.conj().T @ left_coordinate.conj().T
            right = right_coordinate @ right_layer @ _tau(right_label, 2) @ right_layer.conj().T @ right_coordinate.conj().T
            actual = left @ _tau(middle_label, 2) @ right
            out_label, out_sign = _image(actual, probe, 2)
            self.assertEqual(out_label, sum(_evaluate(poly, assignment) << bit for bit, poly in enumerate(labels)))
            self.assertEqual(out_sign, _evaluate(sign, assignment))

    def test_two_cube_error_budget_uses_the_actual_source_inverse(self):
        source = T @ H @ T
        q = X
        returned = [q] + [np.eye(2, dtype=complex) for _ in range(7)]
        native = [source.conj().T @ value @ source for value in returned]
        reference = source.conj().T @ q @ source
        np.testing.assert_allclose(_product(native), reference, atol=4e-14)
        eta = 0.001
        cubes = []
        for orientation in (1, -1):
            approximate = []
            for index, value in enumerate(native):
                angle = orientation * eta * (index + 1) / 8
                perturbation = math.cos(angle) * np.eye(2) + 1j * math.sin(angle) * Z
                candidate = perturbation @ value
                self.assertLessEqual(np.linalg.norm(returned[index] - source @ candidate @ source.conj().T,
                                                   ord=2), 2 * eta)
                approximate.append(candidate)
            cube = _product(approximate)
            self.assertLessEqual(np.linalg.norm(cube - reference, ord=2), 16 * eta)
            cubes.append(cube)
        self.assertLessEqual(np.linalg.norm(cubes[0] - cubes[1], ord=2), 32 * eta)
        radius = math.sqrt(2) / 64
        self.assertAlmostEqual(32 * radius, 1 / math.sqrt(2))
        self.assertLess(16 * radius, math.sqrt(2))
        # This test checks the norm telescope. Perturbed cube factors are
        # not assumed Clifford; the next test checks the Clifford detector.

    def test_nonidentity_clifford_target_has_a_signed_or_symplectic_detector(self):
        branches = set()
        for target in (Z, H, T @ X @ T.conj().T):
            matrix, signs = _tableau(target, 1)
            self.assertGreaterEqual(_projective_distance(target, np.eye(2)) + 1e-12, 1 / math.sqrt(2))
            different = np.argwhere(matrix != np.eye(2, dtype=int))
            for bit in (0, 1):
                value = np.exp(0.173j) * (target if bit else np.eye(2))
                actual_matrix, actual_signs = _tableau(value, 1)
                if len(different):
                    branches.add('symplectic')
                    row, column = different[0]
                    self.assertEqual(int(actual_matrix[row, column]) ^ int(row == column), bit)
                else:
                    branches.add('sign')
                    probe = signs.index(1)
                    self.assertEqual(actual_signs[probe], bit)
        self.assertEqual(branches, {'symplectic', 'sign'})

    def test_nonpauli_odd_order_baseline_does_not_cancel_in_the_cube(self):
        work = 2
        coordinate = _word_matrix(_coordinate_word(work), work)
        diagonal = _word_matrix([('T', wire) for wire in range(work)], work)
        all_x = _word_matrix([('X', wire) for wire in range(work)], work)
        a = diagonal @ all_x @ diagonal.conj().T
        n = coordinate @ a @ coordinate.conj().T
        baseline = a @ n
        values = [baseline if address & 1 else np.eye(4) for address in range(8)]
        cube = _product([values[address] for address in _cube_addresses(0, (1, 2, 4))])
        np.testing.assert_allclose(cube, np.linalg.matrix_power(baseline, 4), atol=5e-14)
        matrix, _ = _tableau(cube, work)
        self.assertFalse(np.array_equal(matrix, np.eye(2 * work, dtype=int)))
        # A linear address label does not cancel a non-Pauli group-valued
        # baseline. The eight-positive-factor cube uses the Pauli hypothesis.

    def test_three_derivatives_detect_any_degree_above_fifty_one(self):
        polynomial = {sum(1 << bit for bit in range(52)), 3, 1 << 53}
        for variable in (0, 1, 2):
            polynomial = {term ^ (1 << variable) for term in polynomial if term >> variable & 1}
        self.assertEqual(_degree(polynomial), 49)
        self.assertGreater(_degree(polynomial), 48)


if __name__ == '__main__':
    unittest.main()
