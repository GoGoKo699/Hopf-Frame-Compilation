"""Exact fixtures for a scoped co-moving quarter-phase lookup trial.

The logical fixtures retain address, bank, selectors and a borrowed parity
bit.  Phases are integer powers of i.  The small native fixture checks the
new shared-control privatization, not a full arbitrary-size stage emitter.
Neither fixture supplies a complete lookup source or a general depth bound.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, CQ2, Q2, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_majorana_walsh_query import _inverse


def _linear(x, pairs):
    for control, target in pairs:
        x ^= ((x >> control) & 1) << target
    return x


class _Matching:
    def __init__(self, p, pairs):
        self.p, self.pairs, self.m = p, pairs, 1 << p
        self.mask = (1 << self.m) - 1
        self.controls = tuple(control for control, _ in pairs)
        self.images = tuple(_linear(u, pairs) for u in range(self.m))
        self.swaps = tuple((u, v) for u, v in enumerate(self.images) if u < v)

    def moments(self, x, bank, selectors):
        distances = tuple(sum(((x ^ u) >> c) & 1 for c in self.controls)
                          for u in range(self.m))
        n = sum(distances[u] for u in range(self.m)
                if (bank >> u) & (selectors >> u) & 1)
        total = sum(distances[u] for u in range(self.m) if (bank >> u) & 1)
        return (bank & selectors).bit_count() % 2, n, total, bank.bit_count()

    def relabel(self, bank):
        return sum(((bank >> u) & 1) << self.images[u]
                   for u in range(self.m))

    def parity_query(self, state):
        x, bank, selectors, borrowed, phase = state
        return x, bank, selectors, borrowed ^ ((bank & selectors).bit_count() % 2), phase

    def controlled_joint_literal(self, state):
        x, bank, selectors, borrowed, phase = state
        if borrowed:
            # Physical address CNOTs and disjoint row swaps on BOTH banks.
            x = _linear(x, self.pairs)
            for u, v in self.swaps:
                bank ^= (((bank >> u) ^ (bank >> v)) & 1) * ((1 << u) | (1 << v))
                selectors ^= (((selectors >> u) ^ (selectors >> v)) & 1) * (
                    (1 << u) | (1 << v))
        return x, bank, selectors, borrowed, phase

    def quarter_phases_literal(self, state, inverse=False):
        x, bank, selectors, borrowed, phase = state
        sign = -1 if inverse else 1
        for control in self.controls:
            for u in range(self.m):
                phase += sign * ((bank >> u) & 1) * ((selectors >> u) & 1) * (
                    ((x ^ u) >> control) & 1)
        return x, bank, selectors, borrowed, phase % 4

    def stage_literal(self, state, inverse=False):
        # Forward chronology A,B,A^dagger,B,D; A and B are literal
        # full-space involutions.  Reverse the actual factor sequence.
        if inverse:
            state = self.quarter_phases_literal(state, inverse=True)
            for operation in (self.controlled_joint_literal, self.parity_query,
                              self.controlled_joint_literal, self.parity_query):
                state = operation(state)
            return state
        for operation in (self.parity_query, self.controlled_joint_literal,
                          self.parity_query, self.controlled_joint_literal):
            state = operation(state)
        return self.quarter_phases_literal(state)

    def stage_formula(self, state):
        x, bank, selectors, borrowed, phase = state
        parity, n, _, _ = self.moments(x, bank, selectors)
        if parity:
            x, bank, selectors = self.images[x], self.relabel(bank), self.relabel(selectors)
        return x, bank, selectors, borrowed, (phase + n) % 4

    def flip(self, state):
        x, bank, selectors, borrowed, phase = state
        return x, bank, selectors ^ self.mask, borrowed, phase

    def echo_literal(self, state):
        state = self.flip(state)
        state = self.stage_literal(state)
        state = self.flip(state)
        return self.stage_literal(state, inverse=True)

    def echo_formula(self, state):
        x, bank, selectors, borrowed, phase = state
        _, n, total, weight = self.moments(x, bank, selectors)
        if weight % 2:
            x, bank, selectors = self.images[x], self.relabel(bank), self.relabel(selectors)
        return x, bank, selectors, borrowed, (phase + total - 2 * n) % 4


def _ccz_7t(a, b, c):
    # Seven parity phases, with a six-CNOT walk restoring all three wires.
    return [('T', a), ('T', b), ('T', c),
            ('CX', a, b), ('TDG', b), ('CX', b, c), ('T', c),
            ('CX', a, c), ('TDG', c), ('CX', b, c), ('TDG', c),
            ('CX', a, c), ('CX', a, b)]


def _shared_control_toffoli(control, u, v, dirty):
    return ([('H', v), ('CX', control, dirty)] + _ccz_7t(dirty, u, v)
            + [('CX', control, dirty)] + _ccz_7t(dirty, u, v) + [('H', v)])


class CoMovingLookupStageTests(unittest.TestCase):
    def assert_complete_action(self, matching, x, bank, selectors, borrowed):
        state = x, bank, selectors, borrowed, 0
        actual = matching.stage_literal(state)
        self.assertEqual(actual, matching.stage_formula(state))
        self.assertEqual(matching.stage_literal(actual, inverse=True), state)
        reverse = matching.stage_literal(state, inverse=True)
        self.assertEqual(matching.stage_literal(reverse), state)
        self.assertEqual((bank >> x) & 1, (actual[1] >> actual[0]) & 1)
        self.assertEqual(matching.echo_literal(state), matching.echo_formula(state))

    def test_every_two_bit_basis_column_and_actual_inverse(self):
        for pairs in (((0, 1),), ((1, 0),)):
            matching = _Matching(2, pairs)
            for x, bank, selectors, borrowed in product(range(4), range(16), range(16), range(2)):
                self.assert_complete_action(matching, x, bank, selectors, borrowed)

    def test_three_bit_backgrounds_and_four_bit_changed_partitions(self):
        for p, families in (
                (3, (((0, 1),), ((1, 2),), ((2, 0),))),
                (4, (((0, 1), (2, 3)), ((1, 0), (3, 2)), ((0, 2), (1, 3))))):
            m, full = 1 << p, (1 << (1 << p)) - 1
            even_rows = sum(1 << u for u in range(m) if u % 2 == 0)
            cube = (1 << 8) - 1
            banks = {0, full, 1, 1 << (m - 1), 3, 1 | (1 << (m - 1)),
                     1 | (1 << (m // 2)), even_rows, full ^ even_rows,
                     cube, full ^ cube}
            selectors = {0, full, even_rows, full ^ even_rows}
            selectors.update(1 << u for u in range(m))
            selectors.update(full ^ (1 << u) for u in range(m))
            for pairs in families:
                matching = _Matching(p, pairs)
                self.assertEqual(len(matching.swaps), (m - m // (1 << len(pairs))) // 2)
                for x, bank, selector, borrowed in product(range(m), sorted(banks),
                                                          sorted(selectors), range(2)):
                    self.assert_complete_action(matching, x, bank, selector, borrowed)

    def test_even_bank_echo_square_and_half_gram_floor(self):
        # Twice a diagonal Gram entry is 1-Re(E); twice the target is 1-Y_x.
        real_i_power = (1, 0, -1, 0)
        for p, pairs in ((2, ((0, 1),)), (3, ((2, 0),)),
                         (4, ((0, 2), (1, 3)))):
            matching = _Matching(p, pairs)
            full = matching.mask
            banks = (range(1, full) if p <= 3 else
                     (3, 5, 17, full ^ 3, 255, full ^ 255))
            for bank in banks:
                if bank.bit_count() % 2:
                    continue
                self.assertTrue(0 < bank < full)
                bank_parity = sum(((u >> c) & 1)
                                  for u in range(matching.m) if (bank >> u) & 1
                                  for c in matching.controls) % 2
                error_twice = 0
                for x, selectors in product(range(matching.m), (0, bank, full)):
                    state = x, bank, selectors, 1, 0
                    echoed = matching.echo_literal(state)
                    self.assertEqual(echoed[:4], state[:4])
                    self.assertEqual(2 * echoed[4] % 4, 2 * bank_parity)
                    gram_twice = 1 - real_i_power[echoed[4]]
                    target_twice = 1 - ((bank >> x) & 1)
                    error_twice = max(error_twice, abs(gram_twice - target_twice))
                self.assertGreaterEqual(error_twice, 1)

    def test_two_rows_have_both_imaginary_eigenvalues_at_every_address(self):
        for p, pairs in ((2, ((0, 1),)), (3, ((2, 1),)),
                         (4, ((0, 2), (1, 3)))):
            matching = _Matching(p, pairs)
            control = pairs[0][0]
            for u in (0, matching.m - 1):
                v = u ^ (1 << control)
                bank = (1 << u) | (1 << v)
                marked, unmarked = set(), set()
                for x in range(matching.m):
                    exponents = set()
                    for first, second in product(range(2), repeat=2):
                        selectors = (first << u) | (second << v)
                        result = matching.echo_literal((x, bank, selectors, 0, 0))
                        self.assertEqual(result[:4], (x, bank, selectors, 0))
                        exponents.add(result[4])
                    self.assertEqual(exponents, {1, 3})
                    (marked if (bank >> x) & 1 else unmarked).update(exponents)
                self.assertEqual(marked, {1, 3})
                self.assertEqual(unmarked, {1, 3})

    def test_four_path_transplant_has_marked_acceptance_three_quarters(self):
        # Address0 and bank e0 stay fixed; every one of the16 arbitrary
        # selector columns remains.  The two matchings permute different
        # inactive selector rows even on this marked address sector.
        stages = []
        for pairs in (((0, 1),), ((1, 0),)):
            matching = _Matching(2, pairs)
            matrix = zeros(16)
            for column in range(16):
                x, bank, output, borrowed, phase = matching.stage_literal((0, 1, column, 0, 0))
                self.assertEqual((x, bank, borrowed, phase), (0, 1, 0, 0))
                matrix[output][column] = CO
            stages.append(matrix)
        a, d = stages
        walsh = [[CO * Fraction((-1)**((row & column).bit_count()), 4)
                  for column in range(16)] for row in range(16)]
        b = multiply(multiply(walsh, d), walsh)
        flip, all_z, p_swap, q_swap, x0, z0 = (zeros(16) for _ in range(6))
        for column in range(16):
            flip[column ^ 15][column] = CO
            all_z[column][column] = CO * (-1)**column.bit_count()
            p_swap[column ^ ((((column >> 1) ^ (column >> 3)) & 1) * 10)][column] = CO
            q_swap[column ^ ((((column >> 2) ^ (column >> 3)) & 1) * 12)][column] = CO
            x0[column ^ 1][column] = CO
            z0[column][column] = CO * (-1)**(column & 1)
        p_plus = scale(add(eye(16), p_swap), Fraction(1, 2))
        p_minus = scale(add(eye(16), scale(p_swap, -CO)), Fraction(1, 2))
        q_plus = scale(add(eye(16), q_swap), Fraction(1, 2))
        q_minus = scale(add(eye(16), scale(q_swap, -CO)), Fraction(1, 2))
        self.assertEqual(a, add(p_plus, multiply(p_minus, z0)))
        self.assertEqual(b, add(q_plus, multiply(q_minus, x0)))
        a_prime = multiply(multiply(flip, a), flip)
        b_prime = multiply(multiply(all_z, b), all_z)
        self.assertEqual(a_prime, add(p_plus, scale(multiply(p_minus, z0), -CO)))
        self.assertEqual(b_prime, add(q_plus, scale(multiply(q_minus, x0), -CO)))
        numerator = add(
            add(multiply(b, a), scale(multiply(b, a_prime), CI)),
            add(scale(multiply(b_prime, a), -CI), scale(multiply(b_prime, a_prime), -CO)))
        root_half = CQ2(Q2(0, Fraction(1, 2)))
        accepted = scale(numerator, CQ2(Q2(0, Fraction(1, 4))))
        simplified = add(
            scale(multiply(multiply(q_plus, p_minus), z0), root_half * (CO - CI)),
            scale(multiply(multiply(q_minus, p_plus), x0), root_half * (CO + CI)))
        self.assertEqual(accepted, simplified)
        gram = scale(multiply(adjoint(numerator), numerator), Fraction(1, 8))
        self.assertEqual(gram, add(multiply(multiply(p_minus, q_plus), p_minus),
                                   multiply(multiply(p_plus, q_minus), p_plus)))
        for selector0 in (0, 1):
            phi = zeros(16, 1)
            phi[2 | selector0][0], phi[8 | selector0][0] = CO, -CO
            self.assertEqual(multiply(p_swap, phi), scale(phi, -CO))
            self.assertEqual(multiply(adjoint(phi), phi), [[CO * 2]])
            self.assertEqual(multiply(multiply(adjoint(phi), q_swap), phi), [[CO]])
            # This unnormalized vector has norm squared2, so acceptance is3/4.
            self.assertEqual(multiply(multiply(adjoint(phi), gram), phi), [[CO * Fraction(3, 2)]])

    def test_native_private_toffoli_and_fredkin_all_columns_and_reverses(self):
        # Shared control0, pair1/2, arbitrary private dirty3.  No T gate
        # touches the shared control, so disjoint pairs can reuse its wire.
        toffoli = _shared_control_toffoli(0, 1, 2, 3)
        fredkin = [('CX', 1, 2)] + _shared_control_toffoli(0, 2, 1, 3) + [('CX', 1, 2)]
        for name, word, cliffords in (('toffoli', toffoli, 16), ('fredkin', fredkin, 18)):
            for direction in (word, _inverse(word)):
                self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in direction), 14)
                self.assertEqual(sum(gate[0] not in ('T', 'TDG') for gate in direction), cliffords)
                self.assertTrue(all(gate[1] != 0 for gate in direction if gate[0] in ('T', 'TDG')))
            actual = native_columns(4, word)
            self.assertEqual(native_columns(4, _inverse(word)), adjoint(actual))
            for column in range(16):
                expected = column
                if column & 1:
                    if name == 'toffoli':
                        expected ^= ((column >> 1) & 1) << 2
                    else:
                        expected ^= (((column >> 1) ^ (column >> 2)) & 1) * 6
                for row in range(16):
                    self.assertEqual(actual[row][column], CO if row == expected else CZ)


if __name__ == '__main__':
    unittest.main()
