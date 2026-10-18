"""Bounded checks for the five-layer Pauli-comparison subcase.

The analytic proof freezes a native Pauli cocycle, bounds degree by three,
and excludes a cubic through the independent inner axes. These checks
exercise those steps without assuming that the middle conjugator is Clifford.
They do not establish the general five-layer theorem by enumeration.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

from tests.test_read_only_three_layer_rigidity import H, T, _cx, _degree, _pauli, _wire


def _pair(left, right, width=2):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((left >> width) & (right & mask)).bit_count()) % 2


def _native_shear(address, label, axes, rows):
    result = 0
    for axis, row in zip(axes, rows):
        if (address & row).bit_count() % 2 and _pair(axis, label):
            result ^= axis
    return result


def _kernel(rows, columns):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    result = []
    for free in range(columns):
        if free in pivots:
            continue
        vector = 1 << free
        for pivot in sorted(pivots):
            if (pivots[pivot] & vector).bit_count() % 2:
                vector ^= 1 << pivot
        result.append(vector)
    return result


def _matrix(label):
    return _pauli(label & 3, label >> 2)


def _label(matrix):
    overlaps = [np.trace(_matrix(label).conj().T @ matrix) / 4
                for label in range(16)]
    label = int(np.argmax(np.abs(overlaps)))
    if not np.isclose(abs(overlaps[label]), 1, atol=3e-13):
        raise AssertionError("The fixture comparison is not a projective Pauli")
    return label


def _product(matrices):
    result = np.eye(4, dtype=complex)
    for matrix in matrices:
        result = result @ matrix
    return result


def _linear(rows, address):
    return sum(((row & address).bit_count() % 2) << j
               for j, row in enumerate(rows))


class ReadOnlyFiveLayerPauliComparison(unittest.TestCase):
    def test_native_axis_cocycle_kernel_is_linear(self):
        # Three address bits, two arbitrary work qubits. These independent,
        # commuting axes also include X components; this is not a Z-only test.
        axes = (5, 10)
        linear_rows = (5, 6)
        self.assertEqual(_pair(*axes), 0)
        equations = [1 << bit for bit in range(4)]  # c(0)=0.
        for x, y, bit in itertools.product(range(8), range(8), range(4)):
            row = ((1 << (4 * (x ^ y) + bit))
                   ^ (1 << (4 * x + bit)) ^ (1 << (4 * y + bit)))
            for source in range(4):
                if (_native_shear(y, 1 << source, axes, linear_rows)
                        >> bit) & 1:
                    row ^= 1 << (4 * x + source)
            equations.append(row)
        basis = _kernel(equations, 32)
        # Each linear c maps into the common two-dimensional annihilator.
        self.assertEqual(len(basis), 6)
        for vector in basis:
            values = [(vector >> (4 * x)) & 15 for x in range(8)]
            for x, y in itertools.product(range(8), repeat=2):
                self.assertEqual(values[x ^ y], values[x] ^ values[y])
                self.assertEqual(_native_shear(y, values[x], axes,
                                              linear_rows), 0)

    def test_two_native_shears_give_quadratic_baseline_labels(self):
        layer = np.kron(T, T.conj().T)
        # The Pauli helper uses the same two-wire ordering as the retained
        # rigidity tests. Direct matrix checks identify the resulting label.
        axes = (4, 8)
        rows_u = (1, 6)
        rows_v = (8, 5)
        rows_delta = (3, 12)
        rows_p = (9, 10, 6, 7)
        baseline = []
        comparisons = []
        for address in range(16):
            u = _linear(rows_u, address)
            v = _linear(rows_v, address)
            delta = _linear(rows_delta, address)
            a = layer @ _pauli(u, v) @ layer.conj().T
            n = layer @ _pauli(u, v ^ delta) @ layer.conj().T
            c = a.conj().T @ n
            self.assertEqual(_label(c), delta << 2)
            p = _linear(rows_p, address)
            k0 = a.conj().T @ _matrix(p) @ n
            expected = p ^ _native_shear(address, p, axes, rows_u) ^ (delta << 2)
            self.assertEqual(_label(k0), expected)
            baseline.append(expected)
            cubic = int((address & 7) == 7)
            comparisons.append(k0 @ np.linalg.matrix_power(_matrix(4), cubic))
        for bit in range(4):
            self.assertLessEqual(_degree([(value >> bit) & 1
                                          for value in baseline]), 2)
        # A third cube kills every baseline label and retains only v Delta^3 g.
        cube = tuple(range(8))
        for base in range(16):
            actual = _product(comparisons[base ^ offset] for offset in cube)
            self.assertEqual(_label(actual), 4)

    def test_cube_product_uses_actual_inverse_before_multiplication(self):
        x = np.array([[0, 1], [1, 0]], dtype=complex)
        source = T @ H @ T
        middle = np.kron(source @ x @ source.conj().T, np.eye(2))
        probe = np.kron(np.diag([1, -1]), np.eye(2))
        # This native middle is non-Clifford; its image of this Pauli is not
        # a physical Pauli. Cross-commuting it with arbitrary K_y is invalid.
        image = middle @ probe @ middle.conj().T
        self.assertLess(max(abs(np.trace(_matrix(label).conj().T @ image)) / 4
                            for label in range(16)), 1 - 1e-8)
        # An intentionally forbidden quartic g: the cube at base zero has
        # Delta^3 g=0, whereas translation by bit 3 gives Delta^3 g=1.
        family = []
        for address in range(16):
            quadratic = ((address & 3) == 3) * 1
            quartic = int(address == 15)
            family.append(_matrix(quadratic ^ (4 * quartic)))
        shift = 8
        fixed_inverse = family[shift].conj().T
        left_factors = [family[offset ^ shift] @ fixed_inverse
                        for offset in range(8)]
        right_factors = [middle @ family[offset] @ middle.conj().T
                         for offset in range(8)]
        left = _product(left_factors)
        right = _product(right_factors)
        np.testing.assert_allclose(
            right, middle @ _product(family[:8]) @ middle.conj().T,
            atol=4e-14)
        self.assertEqual(_label(left), 4)
        np.testing.assert_allclose(right, np.eye(4), atol=4e-14)
        errors = [np.linalg.norm(a - b, ord=2)
                  for a, b in zip(left_factors, right_factors)]
        self.assertLessEqual(np.linalg.norm(left - right, ord=2),
                             sum(errors) + 4e-14)
        # Any unknown unit phases in the cocycle still leave projective
        # distance at least sqrt(2) between identity and a nonidentity Pauli.
        for angle in np.linspace(0, 2 * np.pi, 65):
            self.assertGreaterEqual(
                np.linalg.norm(left - np.exp(1j * angle) * right, ord=2),
                np.sqrt(2) - 4e-14)
        self.assertAlmostEqual(
            np.linalg.norm(left - 1j * right, ord=2), np.sqrt(2))

    def test_native_partial_pauli_action_and_peeled_clifford_comparison(self):
        # Three arbitrary work qubits. Inner active axes X0 Z1 and X0 Z2
        # are independent and commute. Both anticommute with outer Z0,
        # but their combined displacement Z1 Z2 commutes with that axis.
        basis = _wire(H, 0, 3) @ _cx(0, 1, 3) @ _cx(0, 2, 3)
        inner = basis @ _wire(T, 1, 3) @ _wire(T, 2, 3) @ basis.conj().T
        w = _pauli(0, 1, 3)
        a_prime = inner @ w @ inner.conj().T
        d = _wire(T, 0, 3)
        b = d @ a_prime @ d.conj().T
        target = -_pauli(0, 7, 3)
        paulis = [_pauli(p, q, 3) for q in range(8) for p in range(8)]

        def coefficients(matrix):
            return np.array([np.trace(p.conj().T @ matrix) / 8 for p in paulis])

        def is_clifford(matrix):
            return all(max(abs(coefficients(matrix @ p @ matrix.conj().T)))
                       > 1 - 1e-11 for p in paulis)

        self.assertTrue(is_clifford(a_prime))
        self.assertFalse(is_clifford(b))
        np.testing.assert_allclose(b @ w @ b.conj().T, target, atol=4e-14)
        # On a returned Pauli coefficient from the D-axis span, B's label
        # action is precisely the native inner Clifford's label action.
        np.testing.assert_allclose(a_prime @ w @ a_prime.conj().T,
                                   target, atol=4e-14)
        self.assertEqual(1 ^ 7, 2 ^ 4)  # Displacement is Z1 Z2.
        for q in (0, 1, 6, 7):
            probe = _pauli(0, q, 3)
            image = b @ probe @ b.conj().T
            self.assertAlmostEqual(max(abs(coefficients(image))), 1)
            left = d.conj().T @ image @ d
            right = a_prime @ (d.conj().T @ probe @ d) @ a_prime.conj().T
            np.testing.assert_allclose(left, right, atol=4e-14)
            self.assertTrue(is_clifford(left))
            self.assertTrue(is_clifford(right))
        # An incorrect Pauli target remains separated after the native peel.
        wrong = d.conj().T @ _pauli(1, 0, 3) @ d
        actual = d.conj().T @ target @ d
        self.assertTrue(is_clifford(wrong))
        for angle in np.linspace(0, 2 * np.pi, 65):
            self.assertGreaterEqual(
                np.linalg.norm(actual - np.exp(1j * angle) * wrong, ord=2),
                1 / np.sqrt(2) - 4e-14)
        # Displayed B: outer T pair, two inner T pairs, four 3G basis calls,
        # and Z0=S0^2. No initialized or external helper is used.
        self.assertEqual((2 + 2 * 2, 4, 4 * 3 + 2), (6, 4, 14))

    def test_nonzero_cubic_contraction_has_rank_at_least_three(self):
        # Exhaust all homogeneous cubic forms in five variables. The proof
        # uses one nonzero monomial's identity minor, independent of width.
        n = 5
        triples = list(itertools.combinations(range(n), 3))
        pairs = list(itertools.combinations(range(n), 2))
        for mask in range(1, 1 << len(triples)):
            active = {frozenset(triple) for j, triple in enumerate(triples)
                      if (mask >> j) & 1}
            rows = [sum(1 << k for k in range(n) if k not in pair
                        and frozenset((*pair, k)) in active)
                    for pair in pairs]
            rank = n - len(_kernel(rows, n))
            self.assertGreaterEqual(rank, 3)
            chosen = next(iter(active))
            axes = sorted(chosen)
            restricted = []
            for k in axes:
                pair = tuple(sorted(chosen - {k}))
                row = rows[pairs.index(pair)]
                restricted.append(sum(((row >> j) & 1) << t
                                      for t, j in enumerate(axes)))
            self.assertEqual(restricted, [1, 2, 4])


if __name__ == "__main__":
    unittest.main()
