"""Finite operator receipts for the logarithmic-precision construction.

The SELECT and router are specified by their bounded blocks and
permutations. These receipts do not estimate their native synthesis costs.
"""

from bisect import bisect_right
from fractions import Fraction
from itertools import accumulate, product
from math import cos, isqrt, sin
import unittest

import numpy as np

from tests.exact_matrices import CI, CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _ccz, _x


def _tensor(left, right):
    return [[a * b for a in row for b in inner]
            for row in left for inner in right]


def _clifford_menu():
    identity = eye(2)
    x = [[CZ, CO], [CO, CZ]]
    y = [[CZ, -CI], [CI, CZ]]
    z = [[CO, CZ], [CZ, -CO]]
    axes = [(y, identity, identity), (x, x, identity), (x, z, identity),
            (x, y, x), (x, y, z)]
    return [eye(8)] + [scale(_tensor(_tensor(a, b), c), CI) for a, b, c in axes]


def _dyadic_matrix(menu, coefficients, denominator):
    result = zeros(8)
    for matrix, coefficient in zip(menu, coefficients):
        result = add(result, scale(matrix, Fraction(coefficient, denominator)))
    return result


def _four_squares(integer):
    pairs = {}
    for a in range(isqrt(integer) + 1):
        for b in range(a, isqrt(integer - a*a) + 1):
            pairs.setdefault(a*a + b*b, (a, b))
    for subtotal, pair in pairs.items():
        if integer - subtotal in pairs:
            return pair + pairs[integer - subtotal]
    raise AssertionError(integer)


def _numeric(matrix):
    return np.array([[entry.display() for entry in row] for row in matrix])


def _hadamard_blocks(blocks):
    blocks = list(blocks)
    stride = 1
    while stride < len(blocks):
        for start in range(0, len(blocks), 2*stride):
            for index in range(start, start + stride):
                lo, hi = blocks[index], blocks[index + stride]
                blocks[index] = add(lo, hi)
                blocks[index + stride] = add(lo, scale(hi, -CO))
        stride *= 2
    # The fixture uses four selector bits, so the factor is rational.
    return [scale(block, Fraction(1, 4)) for block in blocks]


def _slot_swap(keys):
    labels = np.arange(8**3)
    payloads = [labels//64, labels//8 % 8, labels % 8]
    active = [j + 1 for j, key in enumerate(keys) if key == 0]
    selected = active[:1]
    for j in selected:
        payloads[j] = payloads[j] ^ payloads[0]
    for j in selected:
        payloads[0] = payloads[0] ^ payloads[j]
    for j in selected:
        payloads[j] = payloads[j] ^ payloads[0]
    return 64*payloads[0] + 8*payloads[1] + payloads[2]


def _apply_payload(matrix, state, index):
    ordered = np.moveaxis(state, index, 0)
    return np.moveaxis(np.tensordot(matrix, ordered, axes=1), 0, index)


def _rotation(angle):
    return np.array([[cos(angle), sin(angle)], [-sin(angle), cos(angle)]])


def _cover_rotation(bad_one, bad_two, angle):
    seed_identity = np.eye(len(bad_one))
    good_one = seed_identity - bad_one
    good_two = seed_identity - bad_two
    first = np.kron(bad_one, np.eye(2)) + np.kron(good_one, _rotation(angle))
    second = np.kron(bad_two, np.eye(2)) + np.kron(good_two, _rotation(angle/2))
    mirror = np.kron(bad_one, np.eye(2)) + np.kron(good_one, np.diag([1, -1]))
    return first @ second @ mirror @ second @ mirror.conj().T


class LogarithmicPrecisionDepthTests(unittest.TestCase):
    def test_exact_dyadic_four_square_unitary_and_accuracy(self):
        menu = _clifford_menu()
        identity = eye(8)
        for j, gamma in enumerate(menu[1:], start=1):
            self.assertEqual(adjoint(gamma), scale(gamma, -CO))
            self.assertEqual(multiply(gamma, gamma), scale(identity, -CO))
            for other in menu[1:j]:
                self.assertEqual(add(multiply(gamma, other), multiply(other, gamma)),
                                 zeros(8))

        for accuracy, angle in [(3, .331), (4, -2.24), (4, 2.7)]:
            with self.subTest(accuracy=accuracy, angle=angle):
                denominator = 1 << (2*accuracy + 4)
                a, b = int(denominator*cos(angle)), int(denominator*sin(angle))
                residual = denominator**2 - a*a - b*b
                coefficients = (a, b) + _four_squares(residual)
                self.assertEqual(sum(value*value for value in coefficients), denominator**2)
                spare = 2*denominator - sum(abs(value) for value in coefficients)
                self.assertGreaterEqual(spare, 0)
                self.assertEqual(spare % 2, 0)
                unitary = _dyadic_matrix(menu, coefficients, denominator)
                self.assertEqual(multiply(adjoint(unitary), unitary), identity)
                self.assertEqual(multiply(unitary, adjoint(unitary)), identity)
                target = np.kron(_rotation(angle), np.eye(4))
                self.assertLess(np.linalg.norm(_numeric(unitary) - target, ord=2), 2**-accuracy)

    def test_exact_half_lcu_amplifies_and_returns_all_payload_columns(self):
        menu = _clifford_menu()
        coefficients = (7, -3, 2, 1, 1, 0)
        denominator = 8
        unitary = _dyadic_matrix(menu, coefficients, denominator)
        signed_menu = [signed for matrix in menu for signed in (matrix, scale(matrix, -CO))]
        counts = [max(sign*value, 0) for value in coefficients for sign in (1, -1)]
        padding = (2*denominator - sum(counts))//2
        counts[0] += padding
        counts[1] += padding
        cutoffs = list(accumulate(counts))
        self.assertEqual(cutoffs[-1], 16)
        selected = [signed_menu[bisect_right(cutoffs, label)] for label in range(16)]
        selected_inverse = [adjoint(matrix) for matrix in selected]
        clean_columns = [eye(8)] + [zeros(8) for _ in range(15)]

        def select_block(blocks, inverted):
            ordered = selected_inverse if inverted else selected
            return _hadamard_blocks([multiply(matrix, block) for matrix, block in
                                     zip(ordered, _hadamard_blocks(blocks))])

        self.assertEqual(select_block(clean_columns, False)[0], scale(unitary, Fraction(1, 2)))
        for reverse in (False, True):
            blocks = clean_columns
            for step, inverted in enumerate((reverse, not reverse, reverse)):
                blocks = select_block(blocks, inverted)
                if step < 2:
                    blocks[0] = scale(blocks[0], -CO)
            intended = adjoint(unitary) if reverse else unitary
            self.assertEqual(blocks, [scale(intended, -CO)] + clean_columns[1:])

    def test_priority_swap_and_returned_echo_on_entangled_payload(self):
        rng = np.random.default_rng(313)
        shape = (2, 2, 8, 8, 8, 2)
        incoming = rng.normal(size=shape) + 1j*rng.normal(size=shape)
        incoming /= np.linalg.norm(incoming)
        b = np.kron(np.diag([1, -1]), np.eye(4))
        menu = _clifford_menu()
        approximants = [_numeric(_dyadic_matrix(menu, coeffs, 8))
                        for coeffs in ((7, -3, 2, 1, 1, 0), (7, 3, 1, 2, 0, 1))]
        reflections = [matrix.conj().T @ b @ matrix for matrix in approximants]

        for keys in product(range(2), repeat=2):
            with self.subTest(keys=keys):
                permutation = _slot_swap(keys)
                reverse = np.argsort(permutation)
                labels = np.arange(8**3).reshape(8, 8, 8)
                selected = keys.index(0) + 1 if 0 in keys else None
                expected_permutation = (labels.swapaxes(0, selected)
                                        if selected is not None else labels)
                np.testing.assert_array_equal(permutation, expected_permutation.ravel())

                def swap(state, inverted=False):
                    columns = permutation if inverted else reverse
                    return state.reshape(8**3, 2)[columns].reshape(8, 8, 8, 2)

                def masked_reflection(state):
                    for j, key in enumerate(keys):
                        if key == 0:
                            state = _apply_payload(reflections[j], state, j+1)
                    return state

                actual = _apply_payload(b, incoming[keys], 0)
                # E dagger in chronological order, followed by B and E.
                actual = swap(masked_reflection(swap(masked_reflection(actual), True)))
                actual = _apply_payload(b, actual, 0)
                actual = masked_reflection(swap(masked_reflection(swap(actual, True))))
                target = np.eye(8) if selected is None else (
                    reflections[selected-1] @ b @ reflections[selected-1] @ b)
                expected = _apply_payload(target, incoming[keys], 0)
                np.testing.assert_allclose(actual, expected, atol=2e-14, rtol=0)

    def test_complementary_mirror_and_near_orthogonal_bad_projectors(self):
        basis = np.eye(8)
        first = np.outer(basis[0], basis[0])
        orthogonal = np.outer(basis[1], basis[1])
        flat = np.ones(8)/np.sqrt(8)
        nearly_flat = .013j*basis[0] + np.sqrt(1-.013**2)*basis[1]
        for second_vector in (basis[1], flat, nearly_flat):
            second = np.outer(second_vector, second_vector.conj())
            mu = np.linalg.norm(second @ first, ord=2)
            projected = (np.eye(8) - first) @ second_vector
            projected /= np.linalg.norm(projected)
            repaired = np.outer(projected, projected.conj())
            self.assertAlmostEqual(np.linalg.norm(second - repaired, ord=2), mu, places=13)
            for angle in (.723, -1.31):
                target = np.kron(np.eye(8), _rotation(angle))
                np.testing.assert_allclose(_cover_rotation(first, orthogonal, angle), target,
                                           atol=2e-14, rtol=0)
                np.testing.assert_allclose(_cover_rotation(first, repaired, angle), target,
                                           atol=2e-14, rtol=0)
                error = np.linalg.norm(_cover_rotation(first, second, angle) - target, ord=2)
                self.assertLessEqual(error, 4*mu + 2e-14)

    def test_native_four_corner_hash_phase_and_actual_inverse(self):
        z, mask_lo, mask_hi, enable = range(4)
        bank_u, bank_v = (4, 5), (6, 7)

        def indicator(bank):
            return _x(z) + [('CX', z, bank[0])] + _x(z) + [('CX', z, bank[1])]

        # F is the zero test of the two arbitrary mask bits.
        negate_mask = _x(mask_lo) + _x(mask_hi)
        hash_xor = (negate_mask + [('H', z)] + _ccz((mask_lo, mask_hi, z))
                    + [('H', z)] + _inverse(negate_mask))
        add_u = indicator(bank_u)
        add_v = hash_xor + indicator(bank_v) + _inverse(hash_xor)
        phase = [gate for u, v in zip(bank_u, bank_v) for gate in _ccz((enable, u, v))]
        word = phase + add_u + phase + add_v + phase + _inverse(add_u) + phase + _inverse(add_v)

        expected = eye(1 << 8)
        for column in range(len(expected)):
            hash_is_zero = bool(column & ((1 << mask_lo) | (1 << mask_hi)))
            if hash_is_zero and column >> enable & 1:
                expected[column][column] = -CO
        actual = native_columns(8, word)
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(8, _inverse(word)), adjoint(actual))

    def test_capped_geometric_msb_conditionals(self):
        for width in range(1, 6):
            cap = 1 << width
            for index in range(cap):
                all_ones = True
                probability = Fraction(1)
                for bit in reversed(range(width)):
                    q = Fraction(1, 1 << (1 << bit))
                    one_probability = q if all_ones else q/(1+q)
                    bit_value = index >> bit & 1
                    probability *= one_probability if bit_value else 1-one_probability
                    all_ones &= bool(bit_value)
                self.assertEqual(probability, Fraction(1, 1 << min(index+1, cap-1)))


if __name__ == '__main__':
    unittest.main(verbosity=2)
