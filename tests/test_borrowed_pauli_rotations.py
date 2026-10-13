"""Exact finite witnesses for borrowed-signal selected Pauli rotations.

The comparisons retain every row, signal and core column, with literal
Q(sqrt(2), i) phases. They audit an operator identity, not a native emitter,
an arbitrary-angle approximation algorithm or an asymptotic resource bound.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, scale, zeros,
)


HALF_ROOT = Q2(0, Fraction(1, 2))
H = [[CQ2(HALF_ROOT), CQ2(HALF_ROOT)],
     [CQ2(HALF_ROOT), CQ2(-HALF_ROOT)]]
Z = [[CO, CZ], [CZ, -CO]]


def _product(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = multiply(result, matrix)
    return result


def _kron(left, right):
    return [[a * b for a in left_row for b in right_row]
            for left_row in left for right_row in right]


def _direct_sum(*blocks):
    result = zeros(sum(len(block) for block in blocks))
    start = 0
    for block in blocks:
        for row, entries in enumerate(block):
            result[start + row][start:start + len(block)] = entries
        start += len(block)
    return result


def _mask_pauli(width, x_mask, z_mask, phase):
    """Literal i**phase X**x Z**z, with bit zero the last tensor factor."""
    result = zeros(1 << width)
    scalar = (CO, CI, -CO, -CI)[phase % 4]
    for column in range(1 << width):
        sign = -1 if (column & z_mask).bit_count() % 2 else 1
        result[column ^ x_mask][column] = scalar * sign
    return result


def _anticommuting_axis(width, x_mask, z_mask):
    if x_mask:
        return _mask_pauli(width, 0, x_mask & -x_mask, 0)
    if z_mask:
        return _mask_pauli(width, z_mask & -z_mask, 0, 0)
    raise ValueError("The borrowed-signal identity requires a nonidentity axis")


def _rotation(pauli, cosine, sine):
    return add(scale(eye(len(pauli)), cosine), scale(pauli, CI * sine))


def _parts(pauli, complement, cosine, sine):
    # Register order within a row: borrowed signal, then arbitrary core.
    core_identity = eye(len(pauli))
    h_signal = _kron(H, core_identity)
    f = _product(h_signal, _direct_sum(core_identity, pauli), h_signal)
    b = _direct_sum(core_identity, complement)
    r = _kron(_rotation(Z, cosine, sine), core_identity)
    return b, f, r


class BorrowedPauliRotationTests(unittest.TestCase):
    def test_unequal_selected_angles_on_all_borrowed_and_core_columns(self):
        # Four coherent rows include Y, negative axes and a two-Y phase.
        # Angles are pi/4, -pi/4, pi/2 and zero, respectively.
        rows = (
            (1, 1, 1, HALF_ROOT, HALF_ROOT),
            (3, 2, 3, HALF_ROOT, -HALF_ROOT),
            (0, 3, 2, Q2(0), Q2(1)),
            (3, 3, 0, Q2(1), Q2(0)),
        )
        parts, expected_rows = [], []
        for x_mask, z_mask, phase, cosine, sine in rows:
            pauli = _mask_pauli(2, x_mask, z_mask, phase)
            complement = _anticommuting_axis(2, x_mask, z_mask)
            self.assertEqual(adjoint(pauli), pauli)
            self.assertEqual(multiply(pauli, pauli), eye(4))
            self.assertEqual(_product(complement, pauli, complement),
                             scale(pauli, -1))
            parts.append(_parts(pauli, complement, cosine, sine))
            expected_rows.append(_kron(eye(2), _rotation(pauli, cosine, sine)))
        b, f, r = (_direct_sum(*(parts[row][index] for row in range(len(rows))))
                   for index in range(3))
        actual = _product(b, f, r, f, b)
        expected = _direct_sum(*expected_rows)
        # Full operator equality covers coherent row selection and arbitrary
        # signal/core states, including entanglement with an external reference.
        self.assertEqual(actual, expected)

    def test_omitting_complement_leaves_opposite_borrowed_signal_branches(self):
        pauli = _mask_pauli(1, 1, 1, 1)  # Y, including its literal i phase.
        complement = _anticommuting_axis(1, 1, 1)
        b, f, r = _parts(pauli, complement, HALF_ROOT, HALF_ROOT)
        target = _rotation(pauli, HALF_ROOT, HALF_ROOT)
        without_complement = _product(f, r, f)
        self.assertEqual(without_complement, _direct_sum(target, adjoint(target)))
        self.assertNotEqual(without_complement, _kron(eye(2), target))
        self.assertEqual(_product(b, f, r, f, b), _kron(eye(2), target))

    def test_mod_four_mask_phase_cannot_be_discarded(self):
        x = [[CZ, CO], [CO, CZ]]
        y = [[CZ, -CI], [CI, CZ]]
        self.assertEqual(_mask_pauli(2, 3, 2, 3), _kron(scale(y, -1), x))
        self.assertEqual(_mask_pauli(2, 3, 2, 7), _mask_pauli(2, 3, 2, 3))
        for x_mask, z_mask, phase in ((3, 2, 3), (0, 3, 2)):
            with self.subTest(x_mask=x_mask, z_mask=z_mask, phase=phase):
                pauli = _mask_pauli(2, x_mask, z_mask, phase)
                complement = _anticommuting_axis(2, x_mask, z_mask)
                b, f, r = _parts(pauli, complement, HALF_ROOT, HALF_ROOT)
                correct = _product(b, f, r, f, b)
                phase_dropped = _mask_pauli(2, x_mask, z_mask, 0)
                wrong_b, wrong_f, wrong_r = _parts(
                    phase_dropped, complement, HALF_ROOT, HALF_ROOT)
                self.assertNotEqual(_product(wrong_b, wrong_f, wrong_r, wrong_f, wrong_b),
                                    correct)

    def test_actual_inverse_reverses_noncommuting_stages_and_boundary_order(self):
        paulis = (_mask_pauli(2, 3, 2, 1), _mask_pauli(2, 0, 2, 0))
        complements = (_anticommuting_axis(2, 3, 2),
                       _anticommuting_axis(2, 0, 2))
        stages, inverse_stages = [], []
        for pauli, complement in zip(paulis, complements):
            b, f, r = _parts(pauli, complement, HALF_ROOT, HALF_ROOT)
            stage = _product(b, f, r, f, b)
            actual_inverse = _product(b, f, adjoint(r), f, b)
            self.assertEqual(actual_inverse, adjoint(stage))
            # F B is the right boundary; its inverse is B F, not F B.
            self.assertNotEqual(_product(f, b, r, f, b), stage)
            stages.append(stage)
            inverse_stages.append(actual_inverse)
        word = multiply(stages[1], stages[0])
        inverse = multiply(inverse_stages[0], inverse_stages[1])
        self.assertEqual(inverse, adjoint(word))
        self.assertEqual(multiply(inverse, word), eye(8))
        wrong_order = multiply(inverse_stages[1], inverse_stages[0])
        self.assertNotEqual(multiply(wrong_order, word), eye(8))


if __name__ == "__main__":
    unittest.main()
