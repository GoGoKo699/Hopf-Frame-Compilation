"""Exact bounded checks of the native two-layer one-clean Galois invariant.

All matrices use Q(sqrt(2), i). The fixtures retain actual inverse words,
occupied and unoccupied columns, and physical T-band independence. They
do not prove cubic-query exclusion or infer C3 membership from rationality.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CQ2, Q2, CZ, CO, CI, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_dirty_indicator_depth import _toffoli_two_layers
from tests.test_t_depth import _word


def _sigma(matrix):
    return [[CQ2(Q2(value.real.a, -value.real.b),
                 Q2(value.imag.a, -value.imag.b)) for value in row]
            for row in matrix]


def _product(*factors):
    result = eye(len(factors[0]))
    for factor in factors:
        result = multiply(result, factor)
    return result


def _conjugate(unitary, operator):
    return _product(unitary, operator, adjoint(unitary))


def _inverse(word):
    inverse = {'H': 'H', 'X': 'X', 'CX': 'CX', 'CZ': 'CZ',
               'T': 'TDG', 'TDG': 'T', 'S': 'SDG', 'SDG': 'S'}
    return [(inverse[gate[0]], *gate[1:]) for gate in reversed(word)]


def _tau(label, width):
    """Hermitian tensor Pauli: x in low width bits, z in the next bits."""
    size = 1 << width
    x, z = label & (size - 1), label >> width
    matrix = [[CZ for _ in range(size)] for _ in range(size)]
    phase = (CO, CI, -CO, -CI)[(x & z).bit_count() % 4]
    for column in range(size):
        matrix[column ^ x][column] = phase * (-1 if (z & column).bit_count() % 2 else 1)
    return matrix


def _coefficients(operator, width):
    size = 1 << width
    result = []
    for label in range(1 << (2 * width)):
        product = multiply(_tau(label, width), operator)
        result.append(sum((product[j][j] for j in range(size)), CZ) / size)
    return result


def _clifford_columns(operator, width):
    paulis = [_tau(label, width) for label in range(1 << (2 * width))]
    columns = []
    for bit in range(2 * width):
        image = _conjugate(operator, paulis[1 << bit])
        matches = [label for label, pauli in enumerate(paulis)
                   if image == pauli or image == scale(pauli, -1)]
        if len(matches) != 1:
            raise ValueError('not a Clifford')
        columns.append(matches[0])
    return columns


def _symplectic(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) % 2


def _displacement_values(columns, width):
    values = []
    for label in range(1 << (2 * width)):
        image = 0
        for bit, column in enumerate(columns):
            if label >> bit & 1:
                image ^= column
        values.append(_symplectic(label, image ^ label, width))
    return values


def _controlled_h_word(clean=0, target=1):
    # Chronological A†, CZ, A, with A=S H T H S† and A Z A†=H.
    # Expand S† and CZ into H/S/CX so the elementary ledger is literal.
    sdg = [('S', target)] * 3
    before = sdg + [('H', target), ('TDG', target), ('H', target), ('S', target)]
    cz = [('H', target), ('CX', clean, target), ('H', target)]
    after = sdg + [('H', target), ('T', target), ('H', target), ('S', target)]
    return before + cz + after


class TwoLayerGaloisReturnTests(unittest.TestCase):
    def test_literal_native_automorphism(self):
        h = native_columns(1, [('H', 0)])
        s = native_columns(1, [('S', 0)])
        t = native_columns(1, [('T', 0)])
        z = _tau(2, 1)
        self.assertEqual(_sigma(h), scale(h, -1))
        self.assertEqual(_sigma(s), s)
        self.assertEqual(_sigma(t), multiply(z, t))
        self.assertEqual(_sigma(adjoint(t)), multiply(z, adjoint(t)))
        cx = native_columns(2, [('CX', 0, 1)])
        self.assertEqual(_sigma(cx), cx)

    def test_controlled_h_exact_return_inverse_and_column_scope(self):
        word = _controlled_h_word()
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 2)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word), 15)
        actual = native_columns(2, word)
        expected = [[CZ for _ in range(4)] for _ in range(4)]
        h = native_columns(1, [('H', 0)])
        for target in range(2):
            expected[2 * target][2 * target] = CO
            for output in range(2):
                expected[1 + 2 * output][1 + 2 * target] = h[output][target]
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(2, _inverse(word)), adjoint(actual))
        primitive_inverse = [expanded for gate in _inverse(word)
                             for expanded in ([('S', gate[1])] * 3
                                              if gate[0] == 'SDG' else [gate])]
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in primitive_inverse), 2)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in primitive_inverse), 31)
        self.assertEqual(native_columns(2, primitive_inverse), adjoint(actual))
        z_clean = _tau(4, 2)
        self.assertEqual(_product(adjoint(actual), _sigma(actual)), z_clean)
        # Every input Pauli commuting with clean Z has rational PTM entries.
        for label in range(16):
            if label & 1 == 0:
                coefficients = _coefficients(_conjugate(actual, _tau(label, 2)), 2)
                self.assertTrue(all(value.imag == 0 and value.real.b == 0
                                    for value in coefficients))
        # The clean-flipping column is X_clean tensor H_target, not rational.
        coefficients = _coefficients(_conjugate(actual, _tau(1, 2)), 2)
        expected_coefficients = [CZ] * 16
        expected_coefficients[3] = expected_coefficients[9] = CQ2(Q2(0, Fraction(1, 2)))
        self.assertEqual(coefficients, expected_coefficients)

    def test_dirty_toffoli_and_unused_clean_phase_keep_literal_return(self):
        schedule = _toffoli_two_layers(0, 1, 2, 3)
        word = _word(schedule)
        self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 2)
        plain = native_columns(5, word)
        expected = [[CZ for _ in range(32)] for _ in range(32)]
        for column in range(32):
            row = column ^ 4 if column & 3 == 3 else column
            expected[row][column] = CO
        self.assertEqual(plain, expected)
        self.assertEqual(_product(adjoint(plain), _sigma(plain)), eye(32))
        # The extra T acts only on the unused clean sector and can occupy
        # the first physical band, disjoint from its four existing targets.
        phased_schedule = [(kind, list(gates)) for kind, gates in schedule]
        first = next(j for j, (kind, _) in enumerate(phased_schedule) if kind == 'T')
        phased_schedule[first][1].append(('T', 4))
        phased_word = _word(phased_schedule)
        self.assertEqual(sum(kind == 'T' for kind, _ in phased_schedule), 2)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in phased_word), 9)
        phased = native_columns(5, phased_word)
        self.assertEqual(phased, multiply(native_columns(5, [('T', 4)]), plain))
        self.assertEqual(native_columns(5, _inverse(phased_word)), adjoint(phased))
        self.assertEqual(_product(adjoint(phased), _sigma(phased)), _tau(1 << 9, 5))
        for column in range(16):
            self.assertEqual([row[column] for row in phased], [row[column] for row in plain])

    def test_native_independence_and_dependent_triangle_are_distinct(self):
        independent = native_columns(2, [('S', 0), ('S', 1)])
        independent_columns = _clifford_columns(independent, 2)
        self.assertIn(1, _displacement_values(independent_columns, 2))
        # Phase S on Z0, Z1 and Z0 Z1: the dependent triangle has zero sum.
        triangle = native_columns(2, [('S', 0), ('S', 1),
                                      ('CX', 0, 1), ('S', 1), ('CX', 0, 1)])
        triangle_columns = _clifford_columns(triangle, 2)
        self.assertEqual(_displacement_values(triangle_columns, 2), [0] * 16)
        self.assertNotEqual(triangle_columns, [1, 2, 4, 8])
        expected = native_columns(2, [('Z', 0), ('Z', 1), ('CZ', 0, 1)])
        self.assertEqual(triangle, expected)

    def test_all_signed_pauli_clean_code_extensions_have_zero_displacement_form(self):
        for label in range(4):
            for sign in (-1, 1):
                pauli = scale(_tau(label, 1), sign)
                extension = [[CZ for _ in range(4)] for _ in range(4)]
                for target in range(2):
                    extension[2 * target][2 * target] = CO
                    for output in range(2):
                        extension[1 + 2 * output][1 + 2 * target] = pauli[output][target]
                self.assertEqual(adjoint(extension), extension)
                self.assertEqual(_sigma(extension), extension)
                columns = _clifford_columns(extension, 2)
                self.assertEqual(_displacement_values(columns, 2), [0] * 16)

    def test_rational_native_two_layer_unitary_need_not_be_c3(self):
        # This scope fixture is not a Boolean query: its columns interfere.
        c = [('CX', 0, 1), ('H', 0), ('CX', 0, 1)]
        d = [('T', 0), ('T', 1)]
        word = d + c + d + c
        self.assertEqual(sum(g[0] == 'T' for g in word), 4)
        self.assertEqual(sum(g[0] != 'T' for g in word), 6)
        actual = native_columns(2, word)
        self.assertEqual(_sigma(actual), actual)
        self.assertEqual(native_columns(2, _inverse(word)), adjoint(actual))
        self.assertTrue(any(sum(row[column] != CZ for row in actual) > 1
                            for column in range(4)))
        p = _tau(1, 2)
        image = _conjugate(actual, p)
        labels = [1, 5, 6, 14]  # X_a, Y_a, Z_a X_b, Z_a Y_b.
        expected = [[CZ for _ in range(4)] for _ in range(4)]
        second_expected = [[CZ for _ in range(4)] for _ in range(4)]
        for index, label in enumerate(labels):
            expected = add(expected, scale(_tau(label, 2), Fraction(1, 2)))
            second_expected = add(second_expected, scale(
                _tau(label, 2), Fraction(-1 if index == 0 else 1, 2)))
        self.assertEqual(image, expected)
        self.assertEqual(_conjugate(image, p), second_expected)
        with self.assertRaises(ValueError):
            _clifford_columns(image, 2)
        # An idle one-clean wire supplies exact rational occupied return and G=I.
        extended = native_columns(3, word)
        self.assertEqual(_product(adjoint(extended), _sigma(extended)), eye(8))

    def test_two_clean_bell_code_can_have_nonpauli_galois_relative(self):
        # W tensor W* fixes the Bell state for every unitary W. Here W=THT,
        # so B†(W tensor W*)B fixes |00> with two physical T bands. A third
        # arbitrary data/reference wire is untouched on every column.
        bell = [('H', 0), ('CX', 0, 1)]
        band = [('T', 0), ('TDG', 1)]
        word = bell + band + [('H', 0), ('H', 1)] + band + _inverse(bell)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 4)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word), 6)
        actual = native_columns(3, word)
        self.assertEqual(native_columns(3, _inverse(word)), adjoint(actual))
        identity = eye(8)
        for column in (0, 4):
            self.assertEqual([row[column] for row in actual],
                             [row[column] for row in identity])
        relative = _product(adjoint(actual), _sigma(actual))
        columns = _clifford_columns(relative, 3)
        self.assertNotEqual(columns, [1 << bit for bit in range(6)])
        self.assertIn(1, _displacement_values(columns, 3))
        for column in (0, 4):
            self.assertEqual([row[column] for row in relative],
                             [row[column] for row in identity])
        # Even a column commuting with both clean checks can be irrational:
        # the coefficient of X0 X1 in V Z0 V† is sqrt(2)/4.
        image = _conjugate(actual, _tau(1 << 3, 3))
        product = multiply(_tau(3, 3), image)
        coefficient = sum((product[j][j] for j in range(8)), CZ) / 8
        self.assertEqual(coefficient, CQ2(Q2(0, Fraction(1, 4))))


if __name__ == '__main__':
    unittest.main()
