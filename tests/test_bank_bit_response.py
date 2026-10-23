"""Exact bank-response and row-selector echo checks in bounded dimensions.

The seven-wire native source fixes two active bank rows and checks every
signal, selector, address and dirty-helper input, plus the actual inverse.
Other fixtures are complete small monomial or Gram identities. They do not
emit a general source or establish an unrestricted native-depth obstruction.
"""
from fractions import Fraction
from itertools import combinations
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, CQ2, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_half_bilinear_gram_source import _private_toffoli
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts, _x
from tests.test_quarter_phase_matching_limits import _compose, _matching_action
from tests.test_synchronized_address_stage import _private_cs


def _identity(p):
    return tuple((x, 0) for x in range(1 << p))


def _row(p, pairs, label):
    return tuple(_matching_action(x, pairs, label) for x in range(1 << p))


def _adjoint_action(word):
    result = [None] * len(word)
    for x, (y, phase) in enumerate(word):
        result[y] = (x, -phase % 4)
    return tuple(result)


def _product(*words):
    result = tuple((x, 0) for x in range(len(words[0])))
    for word in words:
        result = _compose(result, word)
    return result


def _stage(p, pairs, bank, selectors):
    result = _identity(p)
    for label in range(1 << p):
        if bank & selectors & (1 << label):
            result = _compose(_row(p, pairs, label), result)
    return result


def _matrix(word):
    result = zeros(len(word))
    phases = (CO, CI, -CO, -CI)
    for x, (y, phase) in enumerate(word):
        result[y][x] = phases[phase]
    return result


def _subtract(left, right):
    return add(left, scale(right, -1))


def _gram(relative):
    return scale(_subtract(scale(eye(len(relative)), 2),
                           add(relative, adjoint(relative))), Fraction(1, 4))


def _trace(matrix):
    return sum((matrix[i][i] for i in range(len(matrix))), CZ)


def _two_row_relative():
    a, au = (_row(2, ((0, 1),), u) for u in (0, 1))
    b = _row(2, ((1, 0),), 0)
    return (_product(b, a), _product(_adjoint_action(au), b, au, a))


def _disjoint_matchings(p):
    edges = tuple((c, t) for c in range(p) for t in range(p) if c != t)
    for size in range(1, p // 2 + 1):
        for pairs in combinations(edges, size):
            if len({wire for pair in pairs for wire in pair}) == 2 * size:
                yield pairs


def _affine_cube(anchor, directions):
    labels = {anchor}
    for direction in directions:
        labels |= {label ^ direction for label in tuple(labels)}
    return sum(1 << label for label in labels)


class BankBitResponseTests(unittest.TestCase):
    def test_two_stage_echo_on_every_two_bit_bank_and_selector_column(self):
        first, second = ((0, 1),), ((1, 0),)
        for bank in range(16):
            a = {s: _stage(2, first, bank, s) for s in range(16)}
            b = {s: _stage(2, second, bank, s) for s in range(16)}
            w = {s: _product(b[s], a[s]) for s in range(16)}
            for selectors in range(16):
                for v in range(4):
                    actual = _product(_adjoint_action(w[selectors]),
                                      w[selectors ^ (1 << v)])
                    if bank >> v & 1:
                        background = _stage(2, first, bank & ~(1 << v), selectors)
                        expected = _product(
                            _adjoint_action(background), _row(2, second, v),
                            background, _row(2, first, v))
                        if selectors >> v & 1:
                            expected = _adjoint_action(expected)
                    else:
                        expected = _identity(2)
                    with self.subTest(bank=bank, selectors=selectors, row=v):
                        # Each tuple includes all address columns and literal phases.
                        self.assertEqual(actual, expected)

    def test_two_row_coherent_defect_and_sharp_selector_dependence(self):
        actions = _two_row_relative()
        self.assertEqual(actions[0], ((0, 0), (2, 2), (3, 1), (1, 1)))
        self.assertEqual(actions[1], ((3, 2), (0, 1), (2, 0), (1, 1)))
        grams = [_gram(_matrix(word)) for word in actions]
        for gram in grams:
            self.assertEqual(gram, adjoint(gram))
            self.assertEqual(multiply(gram, gram), scale(gram, Fraction(3, 4)))
            self.assertEqual(_trace(gram), CQ2(Fraction(3, 2)))
        difference = _subtract(grams[1], grams[0])
        expected = scale([
            [2*CO, -CI, CZ, CO], [CI, CZ, -CO, CZ],
            [CZ, -CO, -2*CO, -CI], [CO, CZ, CI, CZ],
        ], Fraction(1, 4))
        self.assertEqual(difference, expected)
        self.assertEqual(multiply(multiply(difference, difference), difference),
                         scale(difference, Fraction(1, 2)))
        self.assertEqual(_trace(difference), CZ)
        self.assertEqual(_trace(multiply(difference, difference)), CO)
        # For s_u=1, this is also the genuine mixed bank response at zero:
        # the row-v echo vanishes when v is inactive, even if u is active.
        face = []
        for bank in range(4):
            branches = [_product(_stage(2, ((1, 0),), bank, selectors),
                                 _stage(2, ((0, 1),), bank, selectors))
                        for selectors in (2, 3)]
            face.append(_gram(_matrix(_product(_adjoint_action(branches[0]),
                                                branches[1]))))
        self.assertEqual(face[0], zeros(4))
        self.assertEqual(face[2], zeros(4))
        self.assertEqual(add(_subtract(_subtract(face[3], face[2]), face[1]), face[0]),
                         difference)
        # Hermiticity and these polynomial/trace identities give eigenvalues
        # -1/sqrt(2),0,0,1/sqrt(2): the common-operator floor is 1/(2sqrt(2)).
        target = zeros(4)
        target[2][2] = target[3][3] = CQ2(Fraction(1, 2))
        marked = [[CO], [CI], [CZ], [CZ]]
        self.assertEqual(multiply(grams[1], marked), scale(marked, Fraction(3, 4)))
        self.assertEqual(multiply(target, marked), zeros(4, 1))
        self.assertEqual(grams[0][0][0], CZ)
        self.assertEqual(grams[1][0][0], CQ2(Fraction(1, 2)))
        # 0<=G_j<=3I/4 and 0<=target<=I/2 supply the matching upper bound 3/4.

    def test_same_bank_marked_and_unmarked_inputs_have_equal_relative_eigenvalue(self):
        q0, q1 = map(_matrix, _two_row_relative())
        marked = [[CO], [CZ], [CZ], [CZ]]
        unmarked = [[CZ], [CZ], [CO], [CZ]]
        self.assertEqual(multiply(q0, marked), marked)
        self.assertEqual(multiply(q1, unmarked), unmarked)
        target = zeros(4)
        target[2][2] = target[3][3] = CQ2(Fraction(1, 2))
        self.assertEqual(multiply(target, marked), zeros(4, 1))
        self.assertEqual(multiply(target, unmarked), scale(unmarked, Fraction(1, 2)))
        # Signal-only processing of Q and Q† cannot distinguish these +1 inputs.

    def test_fixed_row_echo_cancels_every_other_two_bit_bank_row(self):
        for pairs in (((0, 1),), ((1, 0),)):
            for bank in range(16):
                for selectors in range(16):
                    stage = _stage(2, pairs, bank, selectors)
                    for v in range(4):
                        flipped = _stage(2, pairs, bank, selectors ^ (1 << v))
                        actual = _product(_adjoint_action(stage), flipped)
                        expected = (_row(2, pairs, v) if bank >> v & 1
                                    else _identity(2))
                        if selectors >> v & 1:
                            expected = _adjoint_action(expected)
                        self.assertEqual(actual, expected)
                        control = pairs[0][0]
                        expected_gram = zeros(4)
                        for x in range(4):
                            expected_gram[x][x] = CQ2(Fraction(
                                (bank >> v & 1) * (((x ^ v) >> control) & 1), 2))
                        self.assertEqual(_gram(_matrix(actual)), expected_gram)

    def test_affine_cube_global_echo_is_identity_only_on_zero_selector_fiber(self):
        fixtures = ((3, _affine_cube(0, (1, 2, 4))),
                    (4, _affine_cube(8, (1, 2, 4))),
                    (4, _affine_cube(2, (3, 5, 8))))
        for p, bank in fixtures:
            self.assertEqual(bank.bit_count(), 8)
            matchings = list(_disjoint_matchings(p))
            self.assertEqual(len(matchings), 6 if p == 3 else 24)
            for pairs in matchings:
                zero = _stage(p, pairs, bank, 0)
                flipped = _stage(p, pairs, bank, (1 << (1 << p)) - 1)
                self.assertEqual(_product(_adjoint_action(zero), flipped), _identity(p))
                for x in range(1 << p):
                    for control, _ in pairs:
                        count = sum((bank >> u & 1) * (((x ^ u) >> control) & 1)
                                    for u in range(1 << p))
                        self.assertIn(count, (0, 4, 8))
        # The completed echo preserves selectors; the identity is not promised
        # on arbitrary selector inputs. One selected row gives a literal minus.
        bank, selectors, pairs = 255, 1, ((0, 1),)
        actual = _product(_adjoint_action(_stage(3, pairs, bank, selectors)),
                          _stage(3, pairs, bank, selectors ^ 255))
        self.assertEqual(actual[1], (1, 2))
        self.assertNotEqual(actual, _identity(3))

    def test_response_errors_accumulate_even_for_positive_integrable_grams(self):
        m, epsilon = 4, Fraction(1, 8)
        target = lambda bank: [Fraction(1 - (bank >> x & 1), 2) for x in range(m)]
        gram = lambda bank: [value + epsilon*bank.bit_count() for value in target(bank)]
        self.assertEqual(gram(0), [Fraction(1, 2)]*m)
        for bank in range(1 << m):
            values = gram(bank)
            self.assertTrue(all(0 <= value <= 1 for value in values))
            self.assertEqual(max(abs(a-b) for a, b in zip(values, target(bank))),
                             epsilon*bank.bit_count())
            absent = [v for v in range(m) if not bank >> v & 1]
            for v in absent:
                actual = [b-a for a, b in zip(values, gram(bank | (1 << v)))]
                self.assertEqual(actual, [epsilon-Fraction(x == v, 2) for x in range(m)])
            for u, v in combinations(absent, 2):
                mixed = [d-b-c+a for a, b, c, d in zip(
                    values, gram(bank | (1 << u)), gram(bank | (1 << v)),
                    gram(bank | (1 << u) | (1 << v)))]
                self.assertEqual(mixed, [0]*m)
        self.assertEqual(gram(15), [Fraction(1, 2)]*m)

    def test_native_two_row_source_and_actual_inverse_on_all_seven_wire_inputs(self):
        # Signal0, selectors sv1/su2, address3/4, arbitrary helpers5/6.
        def row_word(selector, control, target, label_control):
            mask = _x(control) if label_control else []
            return (mask + _private_cs(selector, control, 5, 6)
                    + _private_toffoli(selector, control, target, 6) + mask)

        first = row_word(1, 3, 4, 0) + row_word(2, 3, 4, 1)
        second = row_word(1, 4, 3, 0) + row_word(2, 4, 3, 0)
        word = ([('H', 0), ('CX', 0, 1)] + first + second
                + [('CX', 0, 1), ('S', 0), ('S', 0), ('H', 0)])
        first_t, first_g = _counts(first)
        second_t, second_g = _counts(second)
        self.assertEqual(_counts(word), (first_t+second_t, first_g+second_g+6))
        self.assertEqual(_counts(word), (80, 102))
        self.assertEqual(_counts(_inverse(word)), (80, 102))
        actual = native_columns(7, word)
        self.assertEqual(native_columns(7, _inverse(word)), adjoint(actual))
        branches = []
        for flip in (0, 1):
            branch = zeros(16)
            for selectors in range(4):
                action = _product(_stage(2, ((1, 0),), 3, selectors ^ flip),
                                  _stage(2, ((0, 1),), 3, selectors ^ flip))
                matrix = _matrix(action)
                for x in range(4):
                    for y in range(4):
                        branch[x+4*selectors][y+4*selectors] = matrix[x][y]
            branches.append(branch)
        for row in range(128):
            for column in range(128):
                expected = CZ
                if row >> 5 == column >> 5:
                    rr = (row >> 3 & 3) + 4*(row >> 1 & 3)
                    cc = (column >> 3 & 3) + 4*(column >> 1 & 3)
                    sign = -1 if (row ^ column ^ 1) & 1 else 1
                    expected = (branches[0][rr][cc] + sign*branches[1][rr][cc])*Fraction(1, 2)
                self.assertEqual(actual[row][column], expected)


if __name__ == '__main__':
    unittest.main()
