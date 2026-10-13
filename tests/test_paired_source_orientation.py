"""Finite all-column checks of the actual paired-source rejection gauge.

Raw native words and Clifford-algebra identities use Q(sqrt(2), i).
Only normalization by a general sqrt(lambda_squared) uses floating point.
These are bounded witnesses, not a scalable orientation-bank emitter.
"""
from fractions import Fraction
from itertools import product
import unittest

import numpy as np

from compiler_robust_hopf.native_residual_rotation import (
    _inverse, _loader, _mask, _scalar,
)
from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns,
    scale, zeros,
)


X = [[CZ, CO], [CO, CZ]]
Z = [[CO, CZ], [CZ, -CO]]


def _product(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = multiply(result, matrix)
    return result


def _kron(*matrices):
    result = [[CO]]
    for matrix in matrices:
        result = [[a * b for a in left_row for b in right_row]
                  for left_row in result for right_row in matrix]
    return result


def _weights(q):
    weights = [Fraction(0)] * (2 * (q + 1))
    weights[0] = Fraction(1, 4)
    for j in range(q):
        tail = Fraction(1, 1 << min(j + 1, q - 1))
        weights[2 * j + 1], weights[2 * j + 2] = tail / 2, tail / 4
    return weights


def _coefficient(weight):
    if not weight:
        return Q2()
    assert weight.numerator == 1
    exponent = weight.denominator.bit_length() - 1
    assert weight.denominator == 1 << exponent
    if exponent % 2:
        return Q2(0, Fraction(1, 1 << ((exponent + 1) // 2)))
    return Q2(Fraction(1, 1 << (exponent // 2)))


def _majorana(core, index):
    qubit, odd = divmod(index, 2)
    result = zeros(1 << core)
    for column in range(1 << core):
        prefix = -1 if (column & ((1 << qubit) - 1)).bit_count() % 2 else 1
        phase = CI * (-1 if column & (1 << qubit) else 1) if odd else CO
        result[column ^ (1 << qubit)][column] = prefix * phase
    return result


def _vector(coefficients, core):
    result = zeros(1 << core)
    for index, coefficient in enumerate(coefficients):
        if coefficient:
            result = add(result, scale(_majorana(core, index), coefficient))
    return result


def _native_raw(q, signs):
    """Retained raw A-dagger B A word, before five-call amplification."""
    signal, target = q + 1, q + 2
    fixed = tuple(int(2 <= index <= 2 * q and index % 2 == 0)
                  for index in range(2 * (q + 1)))
    loader = _loader(q)
    cz = (("H", (target,)), ("CX", (signal, target)), ("H", (target,)))
    cx = (("CX", (signal, target)),)
    left = cz + _scalar(loader, _mask(fixed), signal) + cz
    middle = cx + _scalar(loader, _mask(tuple(int(sign < 0) for sign in signs)), signal) + cx
    word = left + middle + _inverse(left)
    native = native_columns(q + 3, [(name, *wires) for name, wires in word])
    # Native order is target, signal, core. Reorder to signal, target, core
    # so each signal block contains EVERY target and arbitrary-core column.
    order = [((column & ((1 << (q + 1)) - 1))
              | (((column >> (q + 2)) & 1) << (q + 1))
              | (((column >> (q + 1)) & 1) << (q + 2)))
             for column in range(1 << (q + 3))]
    return [[native[row][column] for column in order] for row in order]


def _fixture(q, signs):
    weights = _weights(q)
    coefficients = [_coefficient(weight) for weight in weights]
    fixed_signs = [-1 if 2 <= index <= 2 * q and index % 2 == 0 else 1
                   for index in range(len(weights))]
    source = _vector(coefficients, q + 1)
    fixed = _vector([a * sign for a, sign in zip(coefficients, fixed_signs)], q + 1)
    programmed = _vector([a * sign for a, sign in zip(coefficients, signs)], q + 1)
    s = sum(w * sign for w, sign in zip(weights, signs))
    tau = sum(w * sign * fixed_sign for w, sign, fixed_sign
              in zip(weights, signs, fixed_signs))
    t = s / 2 - tau
    l_vector = add(add(programmed, scale(fixed, 2 * t)), scale(source, -(s + t)))
    accepted = _kron(add(scale(eye(2), s), scale(multiply(X, Z), t)), eye(len(source)))
    rejected = _kron(X, multiply(source, l_vector))
    return source, fixed, programmed, s, t, l_vector, accepted, rejected


def _complex(matrix):
    return np.array([[entry.display() for entry in row] for row in matrix])


class PairedSourceOrientationTests(unittest.TestCase):
    def test_actual_raw_word_has_the_exact_rejection_orientation(self):
        # q=2 is a small algebra fixture, not the q>=5 approximation theorem.
        q, signs = 2, (1, -1, 1, 1, -1, 1)
        source, fixed, programmed, s, t, l_vector, accepted, rejected = _fixture(q, signs)
        identity = eye(len(source))
        self.assertNotEqual(t, 0)
        self.assertEqual(add(multiply(source, l_vector), multiply(l_vector, source)),
                         zeros(len(source)))
        lambda_squared = 1 - s * s - t * t
        self.assertGreater(lambda_squared, 0)
        self.assertEqual(multiply(l_vector, l_vector), scale(identity, lambda_squared))
        df = add(multiply(source, fixed), scale(identity, -Fraction(1, 2)))
        dg = add(multiply(source, programmed), scale(identity, -s))
        self.assertEqual(add(dg, scale(df, 2 * t)), multiply(source, l_vector))

        whole_identity = eye(4 * len(source))
        left = add(scale(whole_identity, Fraction(1, 2)), _kron(X, Z, df))
        middle = add(scale(whole_identity, s), _kron(X, X, dg))
        raw = _product(adjoint(left), middle, left)
        self.assertNotEqual(raw, _product(left, middle, left))
        size = len(accepted)
        expected = [accepted[row] + rejected[row] for row in range(size)]
        expected += [rejected[row] + accepted[row] for row in range(size)]
        self.assertEqual(raw, expected)
        self.assertEqual(adjoint(rejected), scale(rejected, -1))
        self.assertEqual(multiply(adjoint(rejected), rejected), scale(eye(size), lambda_squared))
        # Exact, unnormalized form of the lower-block gauge identity.
        self.assertEqual(_product(adjoint(rejected), accepted, rejected),
                         scale(adjoint(accepted), lambda_squared))
        # Also evaluate the literal retained gate word independently. A q=1
        # witness keeps exact all-column native evaluation small; the q=2
        # algebra above includes two coordinates in each geometric tail.
        native_signs = (1, -1, 1, 1)
        *_, native_accepted, native_rejected = _fixture(1, native_signs)
        native_size = len(native_accepted)
        native_expected = [native_accepted[row] + native_rejected[row]
                           for row in range(native_size)]
        native_expected += [native_rejected[row] + native_accepted[row]
                            for row in range(native_size)]
        self.assertEqual(_native_raw(1, native_signs), native_expected)

    def test_complete_bank_gauge_has_adjoint_lower_block_and_literal_sign(self):
        # General nonzero rejection, then the two zero-rejection mask endpoints.
        for signs in ((1, -1, 1, 1, -1, 1), (1,) * 6, (-1,) * 6):
            with self.subTest(signs=signs):
                source, _, _, s, t, _, accepted_exact, rejected_exact = _fixture(2, signs)
                accepted, rejected = _complex(accepted_exact), _complex(rejected_exact)
                size = len(accepted)
                identity = np.eye(size)
                zero = np.zeros((size, size))
                lambda_squared = 1 - s * s - t * t
                lam = np.sqrt(float(lambda_squared))
                orientation = (rejected / lam if lam else _complex(_kron(X, eye(len(source)))))
                np.testing.assert_allclose(orientation.conj().T @ orientation, identity,
                                           atol=2e-14, rtol=0)
                raw = np.block([[accepted, rejected], [rejected, accepted]])
                gauge = np.block([[identity, zero], [zero, orientation.conj().T]])
                canonical = gauge @ raw @ gauge.conj().T
                expected = np.block([[accepted, -lam * identity],
                                     [lam * identity, accepted.conj().T]])
                np.testing.assert_allclose(canonical, expected, atol=2e-14, rtol=0)
                if lam:
                    self.assertGreater(np.linalg.norm(canonical[size:, size:] - accepted), .1)
                    self.assertGreater(np.linalg.norm(canonical[:size, size:] - lam * identity), 1)
                else:
                    self.assertEqual(rejected_exact, zeros(size))
                    self.assertEqual(accepted_exact, scale(eye(size), s))
                    self.assertIn(s, (-1, 1))

    def test_odd_even_sign_classes_give_at_most_five_orthogonal_directions(self):
        attained_five = False
        for q in (2, 3):
            weights = _weights(q)
            coefficients = [_coefficient(weight) for weight in weights]
            # Include empty classes, both head signs and equal final weights.
            for signs in product((-1, 1), repeat=2 * q + 1):
                groups = [[0]]
                labels = [(1, signs[0])]
                for parity in (1, 0):
                    indices = range(1 if parity else 2, 2 * q + 1, 2)
                    for sign in (-1, 1):
                        group = [index for index in indices if signs[index] == sign]
                        if group:
                            groups.append(group)
                            labels.append((1 if parity else -1, sign))
                self.assertLessEqual(len(groups), 5)
                attained_five |= len(groups) == 5
                vectors = [[coefficient if index in group else Q2()
                            for index, coefficient in enumerate(coefficients)] for group in groups]
                # Their Gram matrix is diagonal with positive entries, so the
                # class-direction rank is exactly the number of nonempty groups.
                for i, vector in enumerate(vectors):
                    for j, other in enumerate(vectors):
                        overlap = sum((a * b for a, b in zip(vector, other)), Q2())
                        if i == j:
                            self.assertEqual(overlap, sum(weights[index] for index in groups[i]))
                            self.assertGreater(overlap, 0)
                        else:
                            self.assertEqual(overlap, 0)
                for index, coefficient in enumerate(coefficients):
                    reconstructed = [sum((vector[index] * factor for vector, factor in zip(vectors, factors)), Q2())
                                     for factors in ([1] * len(groups),
                                                     [label[0] for label in labels],
                                                     [label[1] for label in labels])]
                    fixed_sign = -1 if index and index % 2 == 0 else 1
                    sign = signs[index] if index < len(signs) else 1
                    self.assertEqual(reconstructed, [coefficient, coefficient * fixed_sign, coefficient * sign])
        self.assertTrue(attained_five)


if __name__ == "__main__":
    unittest.main()
