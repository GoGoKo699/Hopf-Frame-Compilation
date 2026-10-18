"""Literal returned-circuit checks for the six-layer odd-order diagnostic.

This deliberately redundant identity circuit is not a synthesis improvement.
It shows that an actual six-slot comparison can contain projective odd order,
and that repeatedly differentiating its defects need not produce Paulis.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

from tests.test_outer_axis_normalizer import _tableau, _tau
from tests.test_read_only_three_layer_rigidity import H, S, T, X, Z, _cx, _wire


Y = 1j * X @ Z
I2 = np.eye(2, dtype=complex)


def _power(operator, bit):
    return operator if bit else I2


def _product(operators):
    result = np.eye(len(operators[0]), dtype=complex)
    for operator in operators:
        result = result @ operator
    return result


def _native_word(inverse=False):
    """Left-to-right matrix factors; address u=0, v=1, arbitrary work=2."""
    phase = [("S", 0)] * 2
    controlled_y = [("S", 2), ("CX", 0, 2)] + [("S", 2)] * 3
    if inverse:
        return ([("CX", 1, 2), ("Td", 2), ("H", 2), ("Td", 2),
                 ("CX", 0, 2), ("S", 2)] + controlled_y
                + [("Td", 2), ("H", 2), ("T", 2), ("CX", 1, 2)]
                + phase)
    return (phase + [("CX", 1, 2), ("Td", 2), ("H", 2), ("T", 2)]
            + controlled_y + [("S", 2)] * 3
            + [("CX", 0, 2), ("T", 2), ("H", 2), ("T", 2),
               ("CX", 1, 2)])


def _word_matrix(word):
    single = {"H": H, "S": S, "T": T, "Td": T.conj().T}
    operators = [_cx(gate[1], gate[2], 3) if gate[0] == "CX"
                 else _wire(single[gate[0]], gate[1], 3) for gate in word]
    return _product(operators)


def _fibre(u, v):
    return (-1) ** u * _product([
        _power(X, v), T.conj().T, H, T, _power(Y, u), S.conj().T,
        _power(X, u), T, H, T, _power(X, v),
    ])


def _comparison(u, v):
    # These are the actual peeled outer family and target comparison from
    # the six-slot word, including its branch phase lambda=(-1)^u.
    outer = H @ T @ _power(Y, u) @ T.conj().T @ H
    target_comparison = (-1) ** u * T @ _power(X, v) @ T.conj().T
    return outer.conj().T @ target_comparison


class ReturnedSixLayerOddGroup(unittest.TestCase):
    def test_literal_full_columns_controls_inverse_and_charges(self):
        native = _native_word()
        inverse = _native_word(inverse=True)
        for word, expected_g in ((native, 15), (inverse, 13)):
            t_count = sum(gate[0] in ("T", "Td") for gate in word)
            self.assertEqual((t_count, len(word) - t_count), (4, expected_g))
            # All four displayed T gates share the one arbitrary work wire.
            self.assertEqual([gate[1] for gate in word if gate[0] in ("T", "Td")],
                             [2, 2, 2, 2])
        actual = _word_matrix(native)
        actual_inverse = _word_matrix(inverse)
        np.testing.assert_allclose(actual, np.eye(8), atol=4e-14)
        np.testing.assert_allclose(actual_inverse, actual.conj().T, atol=4e-14)
        np.testing.assert_allclose(actual_inverse @ actual, np.eye(8), atol=5e-14)
        for u, v in itertools.product(range(2), repeat=2):
            indices = [u + 2 * v, u + 2 * v + 4]
            np.testing.assert_allclose(actual[np.ix_(indices, indices)],
                                       _fibre(u, v), atol=4e-14)
            np.testing.assert_allclose(_fibre(u, v), I2, atol=4e-14)
        # Omitting the charged address phase breaks literal operator return.
        without_phase = _word_matrix(native[2:])
        np.testing.assert_allclose(without_phase, _wire(Z, 0, 3), atol=4e-14)
        self.assertGreater(np.linalg.norm(without_phase - np.eye(8), ord=2), 1.9)

    def test_actual_comparison_has_projective_order_three(self):
        a = T.conj().T @ X @ T
        b = H @ T @ X @ T.conj().T @ H
        prefix = H @ T @ S.conj().T @ T
        np.testing.assert_allclose(prefix, H, atol=3e-14)
        expected = {(1, 0): [[1, 1], [0, 1]],
                    (0, 1): [[1, 0], [1, 1]],
                    (1, 1): [[0, 1], [1, 1]]}
        for u, v in itertools.product(range(2), repeat=2):
            actual = _comparison(u, v)
            inner = _power(a, u) @ _power(b, v)
            np.testing.assert_allclose(actual, prefix @ inner @ prefix.conj().T,
                                       atol=4e-14)
            if u or v:
                matrix, _ = _tableau(actual, 1)
                np.testing.assert_array_equal(matrix, expected[u, v])
        product = _comparison(1, 1)
        np.testing.assert_allclose(np.linalg.matrix_power(product, 3), -I2,
                                   atol=5e-14)
        self.assertLess(abs(np.trace(product)) / 2, 1 - 1e-8)

    def test_actual_defect_has_literal_order_three(self):
        a = T.conj().T @ X @ T
        b = H @ T @ X @ T.conj().T @ H
        inner = a @ b @ a @ b
        actual = (_comparison(1, 1) @ _comparison(1, 0).conj().T
                  @ _comparison(0, 1).conj().T)
        np.testing.assert_allclose(actual, H @ inner @ H, atol=5e-14)
        np.testing.assert_allclose(np.linalg.matrix_power(inner, 3), I2,
                                   atol=5e-14)
        matrix, _ = _tableau(actual, 1)
        np.testing.assert_array_equal(matrix, [[1, 1], [1, 0]])
        np.testing.assert_array_equal(matrix @ matrix @ matrix % 2,
                                      np.eye(2, dtype=int))
        self.assertFalse(np.array_equal(matrix, np.eye(2, dtype=int)))

    def test_repeated_multiplicative_defects_stabilize_away_from_paulis(self):
        # Fix y=(1,0), restrict x=(0,u), and form the actual ordered defect.
        values = [(_comparison(1, u) @ _comparison(1, 0).conj().T
                   @ _comparison(0, u).conj().T) for u in range(2)]
        c = values[1]
        np.testing.assert_allclose(values[0], I2, atol=4e-14)
        np.testing.assert_allclose(np.linalg.matrix_power(c, 3), I2, atol=6e-14)
        for _ in range(6):
            values = [values[u ^ 1] @ values[u].conj().T for u in range(2)]
            for u, actual in enumerate(values):
                np.testing.assert_allclose(actual, c if u == 0 else c.conj().T,
                                           atol=2e-12)
                self.assertLess(max(abs(np.trace(_tau(label, 1) @ actual)) / 2
                                    for label in range(4)), 1 - 1e-8)


if __name__ == "__main__":
    unittest.main()
