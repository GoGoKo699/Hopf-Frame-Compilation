"""Bounded checks for the grouped programmed precision tail.

The graph completion uses exact integer numerators. The selected operator
check retains every payload and auxiliary lane, and the two-slot displacement
retains every address, selector and arbitrary bank column. The fixtures use
matrix identities, without certifying asymptotic native schedules.
"""
from itertools import product
from math import isqrt
import unittest

import numpy as np


def _four_squares(value):
    """This finite search is only for the fixed small fixture integers."""
    for a in range(isqrt(value), -1, -1):
        remainder_a = value - a * a
        for b in range(isqrt(remainder_a), -1, -1):
            remainder_b = remainder_a - b * b
            for c in range(isqrt(remainder_b), -1, -1):
                remainder_c = remainder_b - c * c
                d = isqrt(remainder_c)
                if d * d == remainder_c:
                    return a, b, c, d
    raise AssertionError(f"No four-square representation for {value}")


def _graph_words(scale, rounded):
    """Return exact dyadic graph numerators over scale squared."""
    two = np.eye(2, dtype=np.int64)
    defect = scale ** 2 * two - rounded.T @ rounded
    off_diagonal = int(defect[0, 1])
    off_sign = 1 if off_diagonal >= 0 else -1
    rows = [(x, off_sign * x) for x in _four_squares(abs(off_diagonal))]
    for column in range(2):
        remainder = int(defect[column, column]) - abs(off_diagonal)
        rows.extend(
            (x, 0) if column == 0 else (0, x)
            for x in _four_squares(remainder)
        )
    graph = np.concatenate((scale * two, rounded, np.array(rows, dtype=np.int64)))
    identity = np.eye(len(graph), dtype=np.int64)
    reflection = graph @ graph.T - scale ** 2 * identity
    swap_reflection = -identity
    swap_reflection[:2, :2] = 0
    swap_reflection[2:4, 2:4] = 0
    swap_reflection[:2, 2:4] = two
    swap_reflection[2:4, :2] = two
    return defect, graph, reflection, reflection @ swap_reflection


def _permutation(images):
    size = len(images)
    return np.eye(size, dtype=complex)[:, images]


def _orientation_echo(selected, main_flip):
    return (main_flip @ selected @ main_flip).conj().T @ selected


def _ideal_graph_word(root):
    ideal = np.eye(16, dtype=complex)
    ideal[:2, :2] = root.T
    ideal[2:4, 2:4] = root
    return ideal


def _selected_residual(word, key_ready, suffix_guard, second_enable):
    """Fixture order: main lane+payload, then slot q+lane+payload."""
    size_slot = len(word)
    slot_columns = np.arange(size_slot)
    selected_q = slot_columns < 4
    centers = slot_columns.copy()
    if key_ready and suffix_guard and second_enable:
        centers[selected_q] ^= 2
    center = _permutation(centers)

    columns = np.arange(4 * size_slot)
    main_lanes, slot_lanes = columns // size_slot, columns % size_slot
    swap_columns = columns.copy()
    if key_ready:
        selected_columns = slot_lanes < 4
        swap_columns[selected_columns] = (
            slot_lanes[selected_columns] * size_slot + main_lanes[selected_columns]
        )
    swap = _permutation(swap_columns)
    main_flip = _permutation(columns ^ (2 * size_slot))
    middle = word.conj().T @ center @ word @ center
    selected = swap.conj().T @ np.kron(np.eye(4), middle) @ swap
    return _orientation_echo(selected, main_flip)


class GroupedProgrammedPrecisionTests(unittest.TestCase):
    def test_exact_dyadic_graph_completion_keeps_all_auxiliary_lanes(self):
        scale = 256
        root = np.array([[4095, -128], [128, 4095]]) / 4097
        contraction = 1 - 16 / scale
        for adjustment in (-2, -1, 1, 2):
            with self.subTest(adjustment=adjustment):
                rounded = np.rint(scale * contraction * root).astype(np.int64)
                rounded[0, 1] += adjustment
                defect, graph, reflection, word = _graph_words(scale, rounded)
                self.assertNotEqual(defect[0, 1], 0)
                identity = np.eye(16, dtype=np.int64)
                np.testing.assert_array_equal(
                    graph.T @ graph, 2 * scale**2 * np.eye(2, dtype=np.int64)
                )
                np.testing.assert_array_equal(reflection @ reflection, scale**4 * identity)
                np.testing.assert_array_equal(word.T @ word, scale**4 * identity)
                ideal = _ideal_graph_word(root)
                delta = word / scale**2 - ideal
                self.assertLess(np.linalg.norm(delta, ord=2), 2.1 * np.sqrt(1 - contraction))

    def test_residual_orientation_returns_slot_and_disables_on_occupied_guards(self):
        root = np.array([[143, -24], [24, 143]]) / 145
        ideal = _ideal_graph_word(root)
        target = np.linalg.matrix_power(root, 4)
        oriented_main = np.zeros((4, 4), dtype=complex)
        oriented_main[:2, :2] = target
        oriented_main[2:, 2:] = target.T
        success_q = np.diag([1] * 4 + [0] * 12)
        idle_q = np.eye(16) - success_q
        expected_active = np.kron(oriented_main, success_q) + np.kron(np.eye(4), idle_q)

        for ready, guard, enable in product((False, True), repeat=3):
            with self.subTest(ready=ready, guard=guard, enable=enable):
                # A malformed occupied key can execute an arbitrary local word.
                slot_word = ideal if ready else _permutation(np.arange(16) ^ 9)
                actual = _selected_residual(slot_word, ready, guard, enable)
                expected = expected_active if ready and guard and enable else np.eye(64)
                np.testing.assert_allclose(actual, expected, atol=2e-15)

        scale = 2**14
        rounded = np.rint((scale - 16) * root).astype(np.int64)
        rounded[0, 1] += 1
        _, _, _, numerator = _graph_words(scale, rounded)
        exact_dyadic = numerator / scale**2
        # The emitted selection can disturb q. Check the actual adjoint both
        # in the inactive cancellation and in the active replacement bound.
        inactive = _selected_residual(exact_dyadic, True, True, False)
        np.testing.assert_allclose(inactive, np.eye(64), atol=2e-15)
        active = _selected_residual(exact_dyadic, True, True, True)
        one_call_error = np.linalg.norm(exact_dyadic - ideal, ord=2)
        echo_error = np.linalg.norm(active - expected_active, ord=2)
        self.assertLess(4 * one_call_error, 0.3)
        self.assertLessEqual(echo_error, 4 * one_call_error + 2e-15)

    def test_coherent_slot_loads_displace_one_shared_arbitrary_bank(self):
        # Low-to-high order: x, two-word bank, alpha_0, alpha_1, two scratch bits.
        columns = np.arange(128, dtype=np.int64)
        alpha_0, alpha_1 = (columns >> 3) & 1, (columns >> 4) & 1

        def fetch(words):
            bit_0 = (words >> (1 + alpha_0)) & 1
            bit_1 = (words >> (1 + alpha_1)) & 1
            return words ^ (bit_0 << 5) ^ (bit_1 << 6)

        for table_word in range(16):
            with self.subTest(table=table_word):
                desired = (table_word >> (2 * (columns & 1))) & 3
                loaded = fetch(columns)
                loaded ^= desired << 1
                loaded = fetch(loaded)
                loaded ^= desired << 1
                expected = columns ^ (((desired >> alpha_0) & 1) << 5)
                expected ^= ((desired >> alpha_1) & 1) << 6
                np.testing.assert_array_equal(loaded, expected)


if __name__ == "__main__":
    unittest.main()
