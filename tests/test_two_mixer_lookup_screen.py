"""Exact finite evidence for the scoped two-Clifford lookup-source screen.

The rank-one exception evaluates its required diagonal truth table directly.
It verifies the complete identity, not a native compiler for that diagonal.
No finite fixture here excludes unrestricted native T-depth.
"""
import unittest

from tests.exact_matrices import (
    CO, CZ, adjoint, eye, multiply, native_columns,
)
from tests.test_common_probe_source import PHASES, _native, _syndrome
from tests.test_joint_lookup_pauli_source import _inverse


def _qbit(bank, p, address):
    selected = address ^ _syndrome(bank, p)
    return (bank.bit_count() ^ (bank >> selected)) & 1


def _fwht(values):
    result = list(values)
    step = 1
    while step < len(result):
        for start in range(0, len(result), 2 * step):
            for offset in range(step):
                i, j = start + offset, start + offset + step
                result[i], result[j] = result[i] + result[j], result[i] - result[j]
        step *= 2
    return result


def _diag(values):
    return [[value if row == column else CZ
             for column in range(len(values))]
            for row, value in enumerate(values)]


def _z(width, mask):
    return _diag([-CO if (row & mask).bit_count() & 1 else CO
                  for row in range(1 << width)])


def _conjugate(unitary, operator):
    return multiply(multiply(unitary, operator), adjoint(unitary))


def _diagonal_pauli_label(matrix):
    """Return the binary Z label, ignoring only the Pauli's scalar."""
    size = len(matrix)
    if any(matrix[row][column] != CZ
           for row in range(size) for column in range(size) if row != column):
        return None
    scalar = matrix[0][0]
    if scalar == CZ:
        return None
    for mask in range(size):
        if all(matrix[row][row] == (
                -scalar if (row & mask).bit_count() & 1 else scalar)
               for row in range(size)):
            return mask
    return None


def _rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def _linear_map(columns, value):
    result = 0
    for index, column in enumerate(columns):
        if (value >> index) & 1:
            result ^= column
    return result


def _rank_one_source(state, p, address, reverse=False):
    """C2 D C1, with D an explicit lookup-dependent sign truth table."""
    m = 1 << p
    gather = [('CX', j, 0) for j in range(1, m)]
    spread = [('CX', 0, j) for j in range(1, m)]
    before, after = gather + [('H', 0)], [('H', 0)] + spread
    if reverse:
        before, after = _inverse(after), _inverse(before)
    state = _native(before, state)
    # L(0,z_rest) has row-zero bit zero and all other bits unchanged.
    state = {bank: (-value if ((bank & 1)
                              and _qbit(bank & ~1, p, address)) else value)
             for bank, value in state.items()}
    return _native(after, state)


class TwoMixerLookupScreenTests(unittest.TestCase):
    def test_complete_p3_fourier_support_and_every_translation(self):
        p, m = 3, 8
        full = (1 << m) - 1
        for address in range(m):
            with self.subTest(address=address):
                truth = [_qbit(bank, p, address) for bank in range(1 << m)]
                transform = _fwht([1 - 2 * bit for bit in truth])
                expected = {}
                for row in range(m):
                    for tag in range(m):
                        h_tag = sum(1 << j for j in range(m)
                                    if (tag & j).bit_count() & 1)
                        mask = full ^ h_tag ^ (1 << row)
                        self.assertNotIn(mask, expected)
                        sign = -1 if (tag & (address ^ row)).bit_count() & 1 else 1
                        expected[mask] = sign * ((1 << m) // m)
                self.assertEqual({mask: value for mask, value in enumerate(transform)
                                  if value}, expected)
                structures = {}
                for shift in range(1 << m):
                    derivative = truth[0] ^ truth[shift]
                    if all((truth[bank] ^ truth[bank ^ shift]) == derivative
                           for bank in range(1 << m)):
                        structures[shift] = derivative
                self.assertEqual(structures, {0: 0, full: 1})
                support = list(expected)
                self.assertEqual(_rank([mask ^ support[0] for mask in support]), m - 1)

    def test_actual_clifford_source_diagonal_has_predicted_affine_support(self):
        width = 3
        c1 = native_columns(width, [('H', 0)])
        c2 = native_columns(width, [('H', 0), ('H', 1), ('CX', 0, 2), ('S', 1)])
        d0 = native_columns(width, [('T', 2)])
        d2 = native_columns(width, [('T', 0), ('T', 1)])
        # A nontrivial diagonal phase, including an exact cubic sign.
        d1 = _diag([PHASES[(bank & 1) + 2 * ((bank >> 2) & 1)]
                    * (-CO if bank == 7 else CO) for bank in range(1 << width)])
        p = _z(width, (1 << width) - 1)
        source = multiply(d2, multiply(c2, multiply(d1, multiply(c1, d0))))
        image = _conjugate(source, p)
        self.assertEqual(multiply(adjoint(source), source), eye(1 << width))

        pauli_base = _conjugate(c2, _conjugate(c1, p))
        diagonal_group, allowed = set(), set()
        for mask in range(1 << width):
            group_element = _conjugate(c2, _z(width, mask))
            label = _diagonal_pauli_label(group_element)
            if label is not None:
                diagonal_group.add(label)
            label = _diagonal_pauli_label(multiply(pauli_base, group_element))
            if label is not None:
                allowed.add(label)
        self.assertEqual(_rank(diagonal_group), 1)
        self.assertEqual(len(allowed), 2)
        u0 = min(allowed)
        self.assertEqual(allowed, {u0 ^ value for value in diagonal_group})
        diagonal = [image[bank][bank] for bank in range(1 << width)]
        spectrum = _fwht(diagonal)
        actual = {mask for mask, value in enumerate(spectrum) if value != CZ}
        self.assertTrue(actual)
        self.assertEqual(actual, allowed)
        for shift in range(1 << width):
            if all(not ((shift & mask).bit_count() & 1) for mask in diagonal_group):
                sign = -1 if (shift & u0).bit_count() & 1 else 1
                self.assertTrue(all(diagonal[bank ^ shift] == sign * diagonal[bank]
                                    for bank in range(1 << width)))

    def test_common_affine_complement_direction_has_even_input_parity(self):
        p, m = 3, 8
        full = (1 << m) - 1
        for special in range(m):
            columns = [full if j == special else 1 << j for j in range(m)]
            self.assertEqual(_rank(columns), m)
            self.assertEqual(columns[special].bit_count() & 1, 0)
            for translation in [0, 1, 0b10010110]:
                truth = [_qbit(_linear_map(columns, value) ^ translation, p, 5)
                         for value in range(1 << m)]
                self.assertTrue(all((truth[value] ^ truth[value ^ (1 << special)]) == 1
                                    for value in range(1 << m)))

    def test_rank_one_exception_has_full_lookup_action_and_actual_inverse(self):
        # This is the genuine exception for the p=3 target, with its missing
        # middle supplied as a complete diagonal truth table, not native gates.
        p, m, address = 3, 8, 5
        for bank in range(1 << m):
            initial = {bank: CO}
            forward = _rank_one_source(initial, p, address)
            returned = _rank_one_source(forward, p, address, reverse=True)
            self.assertEqual(returned, initial)
            state = _rank_one_source(initial, p, address, reverse=True)
            state = {value: (-amplitude if value.bit_count() & 1 else amplitude)
                     for value, amplitude in state.items()}
            state = _rank_one_source(state, p, address)
            self.assertEqual(state, {bank: -CO if _qbit(bank, p, address) else CO})

    def test_rank_one_middle_retains_high_boolean_degree(self):
        p, m, address = 3, 8, 5
        q = [_qbit(bank, p, address) for bank in range(1 << m)]
        g = [_qbit(rest << 1, p, address) for rest in range(1 << (m - 1))]
        middle = [(bank & 1) * _qbit(bank & ~1, p, address)
                  for bank in range(1 << m)]
        for values, degree in [(q, p + 1), (g, p + 1), (middle, p + 2)]:
            coefficients = list(values)
            for bit in range(len(values).bit_length() - 1):
                for mask in range(len(values)):
                    if (mask >> bit) & 1:
                        coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
            self.assertEqual(max(mask.bit_count() for mask, value
                                 in enumerate(coefficients) if value), degree)


if __name__ == '__main__':
    unittest.main()
