"""Exact complete-column checks for coupled quarter-phase acceptance.

The six-wire native word fixes one active row, not a clean row bank.
Separate integer arithmetic checks every two-bit arbitrary-bank fiber.
Neither fixture supplies a complete lookup source or a general depth bound.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, CQ2, Q2, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_half_bilinear_gram_source import _private_toffoli
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts
from tests.test_synchronized_address_stage import _private_cs


def _one_row_operators():
    # Logical ordering is address x0+2*x1, then the one active selector.
    a, d, flip, z = (zeros(8) for _ in range(4))
    for selector in (0, 1):
        for address in range(4):
            column = address + 4*selector
            ac = selector * (address & 1)
            dc = selector * (address >> 1)
            a[(address ^ (ac << 1)) + 4*selector][column] = CI if ac else CO
            d[(address ^ dc) + 4*selector][column] = CI if dc else CO
            flip[address + 4*(selector ^ 1)][column] = CO
            z[column][column] = CO if selector == 0 else -CO
    # H_s D H_s, evaluated by its two selector branches without a square root.
    b = zeros(8)
    for row_selector in (0, 1):
        for col_selector in (0, 1):
            sign = 1 if row_selector == col_selector else -1
            for x in range(4):
                for y in range(4):
                    b[x + 4*row_selector][y + 4*col_selector] = (
                        (CO if x == y else CZ) + sign*d[x+4][y+4]
                    ) * Fraction(1, 2)
    return a, b, flip, z


def _gaussian_add(a, b):
    return a[0]+b[0], a[1]+b[1]


def _gaussian_phase(a, exponent):
    # All amplitudes below are pairs of ordinary integers, not float complex.
    return (a, (-a[1], a[0]), (-a[0], -a[1]), (a[1], -a[0]))[exponent % 4]


def _trim(state):
    return {key: value for key, value in state.items() if value != (0, 0)}


def _row_product(state, bank, control, target, inverse=False):
    """Apply each literal anchored row factor, independently of count formulas."""
    result = {}
    for (address, selectors), amplitude in state.items():
        exponent = 0
        labels = range(3, -1, -1) if inverse else range(4)
        for label in labels:
            if (bank >> label & 1) and (selectors >> label & 1):
                mismatch = (address ^ label) >> control & 1
                exponent += (-1 if inverse else 1) * mismatch
                address ^= mismatch << target
        key = address, selectors
        result[key] = _gaussian_add(result.get(key, (0, 0)),
                                    _gaussian_phase(amplitude, exponent))
    return _trim(result)


def _walsh(state):
    # Four selector Hadamards have normalization 1/4. Track this separately.
    result = {}
    for (address, selectors), amplitude in state.items():
        for output in range(16):
            key = address, output
            exponent = 2*((selectors & output).bit_count() % 2)
            result[key] = _gaussian_add(result.get(key, (0, 0)),
                                        _gaussian_phase(amplitude, exponent))
    return _trim(result)


def _flip_all(state):
    return {(x, s ^ 15): amplitude for (x, s), amplitude in state.items()}


def _z_all(state):
    return {(x, s): _gaussian_phase(amplitude, 2*(s.bit_count() % 2))
            for (x, s), amplitude in state.items()}


def _b_prime(state, bank):
    return _z_all(_walsh(_row_product(_walsh(_z_all(state)), bank, 1, 0)))


def _b_inverse(state, bank):
    return _walsh(_row_product(_walsh(state), bank, 1, 0, inverse=True))


def _relative_explicit(address, selectors, bank, crossed=False):
    state = {(address, selectors): (1, 0)}
    if crossed:
        state = _b_prime(state, bank)
    state = _flip_all(_row_product(_flip_all(state), bank, 0, 1))
    if not crossed:
        state = _b_prime(state, bank)
    state = _b_inverse(state, bank)
    return _row_product(state, bank, 0, 1, inverse=True)


def _relative_count_formula(address, selectors, bank):
    def mask(control, value):
        return sum(1 << u for u in range(4)
                   if bank >> u & 1 and (u >> control & 1) != value)
    c_mask = mask(0, address & 1)
    k = c_mask.bit_count() - (c_mask & selectors).bit_count()
    x1 = address ^ ((k & 1) << 1)
    d_mask = mask(1, x1 >> 1)
    n_d = d_mask.bit_count()
    s1 = selectors ^ d_mask
    x2 = x1 ^ (n_d & 1)
    k2 = (mask(0, x2 & 1) & s1).bit_count()
    return x2 ^ ((k2 & 1) << 1), s1, (k+n_d-k2) % 4


class CoupledQuarterPhaseSourceTests(unittest.TestCase):
    def test_native_word_and_actual_reverse_keep_every_helper_column(self):
        # Signal0, selector1, address2,3, arbitrary dirty helpers4,5.
        c_stage = _private_cs(1, 2, 4, 5) + _private_toffoli(1, 2, 3, 5)
        d_stage = _private_cs(1, 3, 4, 5) + _private_toffoli(1, 3, 2, 5)
        fan = [('CX', 0, 1)]
        word = ([('H', 0)] + fan + c_stage + fan + [('H', 1)]
                + fan + d_stage + fan
                + [('H', 1), ('S', 0), ('S', 0), ('H', 0)])
        self.assertEqual(_counts(word), (40, 54))
        self.assertEqual(_counts(_inverse(word)), (40, 54))
        actual = native_columns(6, word)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))
        a, b, flip, z = _one_row_operators()
        a_prime = multiply(multiply(flip, a), flip)
        b_prime = multiply(multiply(z, b), z)
        w0, w1 = multiply(b, a), multiply(b_prime, a_prime)
        for row in range(64):
            for column in range(64):
                expected = CZ
                if row >> 4 == column >> 4:
                    rr = (row >> 2 & 3) + 4*(row >> 1 & 1)
                    cc = (column >> 2 & 3) + 4*(column >> 1 & 1)
                    sign = -1 if (row ^ column ^ 1) & 1 else 1
                    expected = (w0[rr][cc] + sign*w1[rr][cc])*Fraction(1, 2)
                self.assertEqual(actual[row][column], expected)
        k = scale(add(w0, scale(w1, -CO)), Fraction(1, 2))
        q = multiply(adjoint(w0), w1)
        expected_gram = scale(add(scale(eye(8), 2),
                                  scale(add(q, adjoint(q)), -CO)), Fraction(1, 4))
        self.assertEqual(multiply(adjoint(k), k), expected_gram)

    def test_all_two_bit_bank_selector_columns_match_complete_relative_action(self):
        for bank in range(16):
            images = set()
            for selectors in range(16):
                for address in range(4):
                    output, final_selectors, exponent = _relative_count_formula(
                        address, selectors, bank)
                    # Four Walsh transforms contribute exactly 1/4**4.
                    self.assertEqual(_relative_explicit(address, selectors, bank),
                                     {(output, final_selectors):
                                      _gaussian_phase((256, 0), exponent)})
                    images.add((output, final_selectors))
            self.assertEqual(len(images), 64)

    def test_one_row_three_cycles_have_an_unwanted_fixed_vector(self):
        a, b, flip, z = _one_row_operators()
        a_prime = multiply(multiply(flip, a), flip)
        b_prime = multiply(multiply(z, b), z)
        q = multiply(adjoint(multiply(b, a)), multiply(b_prime, a_prime))
        self.assertEqual(multiply(multiply(q, q), q), eye(8))
        # Marked x=0 and an unmarked vector have the same eigenvalue one.
        marked = [[CO if j == 0 else CZ] for j in range(8)]
        unmarked = [[CZ] for _ in range(8)]
        unmarked[1][0], unmarked[6][0], unmarked[3][0] = CO, -CO, -CI
        self.assertEqual(multiply(q, marked), marked)
        self.assertEqual(multiply(q, unmarked), unmarked)
        target = [[CO*Fraction(1, 2) if i == j and i % 4 else CZ
                   for j in range(8)] for i in range(8)]
        self.assertEqual(multiply(target, marked), zeros(8, 1))
        self.assertEqual(multiply(target, unmarked), scale(unmarked, Fraction(1, 2)))
        self.assertEqual(q[6][1], -CO)
        self.assertEqual(q[3][6], CI)
        self.assertEqual(q[1][3], CI)

    def test_wrong_flip_cancels_second_stage_and_crossing_alone_still_fails(self):
        a, b, flip, z = _one_row_operators()
        self.assertEqual(multiply(multiply(flip, b), flip), b)
        self.assertNotEqual(multiply(multiply(z, b), z), b)
        # Crossed relative order A†B†A'B', two active rows, marked x=0.
        # Every surviving selector/address output is orthogonal to the input.
        expected = {(2, 1): (0, 128), (0, 2): (128, 0),
                    (3, 0): (128, 0), (1, 3): (0, 128)}
        actual = _relative_explicit(0, 0, 3, crossed=True)
        self.assertEqual(actual, expected)
        self.assertNotIn((0, 0), actual)
        # Its Gram expectation is 1/2, whereas Y_0=1 requires zero.
        self.assertEqual(Fraction(256-actual.get((0, 0), (0, 0))[0], 512),
                         Fraction(1, 2))

    def test_row_selective_difference_has_complete_local_gram(self):
        phases = (CO, CI, -CO, -CI)
        for bank in range(16):
            for selectors in range(16):
                for label in range(4):
                    for control in (0, 1):
                        mismatching = [u for u in range(4) if (u & 1) != control]
                        counts = [sum((bank >> u & 1)*(s >> u & 1)
                                      for u in mismatching)
                                  for s in (selectors, selectors ^ (1 << label))]
                        k = zeros(2)
                        for target in (0, 1):
                            for sign, count in zip((1, -1), counts):
                                k[target ^ (count & 1)][target] += (
                                    phases[count % 4] * Fraction(sign, 2))
                        value = Fraction((bank >> label & 1)
                                         * int(control != (label & 1)), 2)
                        self.assertEqual(multiply(adjoint(k), k),
                                         scale(eye(2), value))

    def test_four_path_native_source_correct_diagonals_do_not_certify_gram(self):
        c_stage = _private_cs(1, 2, 4, 5) + _private_toffoli(1, 2, 3, 5)
        d_stage = _private_cs(1, 3, 4, 5) + _private_toffoli(1, 3, 2, 5)
        fan = [('CX', 0, 1)]
        # With one active bank row, the paid bank-quarter phases are S_a/S_a†.
        word = ([('H', 0)] + fan + c_stage + fan
                + [('S', 0), ('H', 0), ('H', 1)] + fan + d_stage + fan
                + [('H', 1), ('S', 0), ('S', 0), ('S', 0), ('H', 0)])
        self.assertEqual(_counts(word), (40, 57))
        self.assertEqual(_counts(_inverse(word)), (40, 57))
        actual = native_columns(6, word)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))
        a, b, flip, z = _one_row_operators()
        a_prime = multiply(multiply(flip, a), flip)
        b_prime = multiply(multiply(z, b), z)
        terms = (multiply(b, a), scale(multiply(b, a_prime), CI),
                 scale(multiply(b_prime, a), -CI), multiply(b_prime, a_prime))
        root_eighth = CQ2(Q2(0, Fraction(1, 4)))
        for row in range(64):
            for column in range(64):
                expected = CZ
                if row >> 4 == column >> 4:
                    rr = (row >> 2 & 3) + 4*(row >> 1 & 1)
                    cc = (column >> 2 & 3) + 4*(column >> 1 & 1)
                    signal_in, signal_out = column & 1, row & 1
                    expected = (terms[0][rr][cc]
                                + (-1)**signal_in*terms[1][rr][cc]
                                + (-1)**signal_out*terms[2][rr][cc]
                                - (-1)**(signal_in+signal_out)*terms[3][rr][cc]
                                ) * root_eighth
                self.assertEqual(actual[row][column], expected)
        numerator = add(add(terms[0], terms[1]), add(terms[2], scale(terms[3], -CO)))
        gram = scale(multiply(adjoint(numerator), numerator), Fraction(1, 8))
        for index in range(8):
            self.assertEqual(gram[index][index],
                             CO*Fraction(1, 2) if index % 4 else CZ)
        self.assertEqual(gram[1][3], -CI*Fraction(1, 4))
        self.assertEqual(multiply(gram, gram), scale(gram, Fraction(3, 4)))
        unmarked_kernel = [[CZ] for _ in range(8)]
        unmarked_kernel[1][0] = CO
        unmarked_kernel[3][0], unmarked_kernel[6][0] = -CI, CI
        self.assertEqual(multiply(gram, unmarked_kernel), zeros(8, 1))
        second_kernel = [[CZ] for _ in range(8)]
        second_kernel[2][0], second_kernel[5][0], second_kernel[7][0] = CO, CI, CO
        self.assertEqual(multiply(gram, second_kernel), zeros(8, 1))
        self.assertEqual(multiply(adjoint(unmarked_kernel), unmarked_kernel), [[CO*3]])
        self.assertEqual(multiply(adjoint(second_kernel), second_kernel), [[CO*3]])
        self.assertEqual(multiply(adjoint(unmarked_kernel), second_kernel), [[CZ]])
        for marked_index in (0, 4):
            marked = [[CO if index == marked_index else CZ] for index in range(8)]
            self.assertEqual(multiply(gram, marked), zeros(8, 1))
        p_unmarked = [[CO if i == j and i % 4 else CZ
                       for j in range(8)] for i in range(8)]
        kernel_sum = add(multiply(unmarked_kernel, adjoint(unmarked_kernel)),
                         multiply(second_kernel, adjoint(second_kernel)))
        self.assertEqual(gram, add(scale(p_unmarked, Fraction(3, 4)),
                                   scale(kernel_sum, Fraction(-1, 4))))
        self.assertEqual(sum(gram[j][j] for j in range(8)), CO*3)


if __name__ == '__main__':
    unittest.main()
