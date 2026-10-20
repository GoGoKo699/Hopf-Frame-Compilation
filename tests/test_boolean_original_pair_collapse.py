"""Small exact guards for Boolean collapse of original physical pairs.

These fixed algebra fixtures protect literal phases, scalar-sector
restriction, actual-list independence and complete clean return. They
do not enumerate native circuits or prove the dimension-uniform theorem.
"""
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_all_active_singleton_closure import _quarter_turn
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_minimal_physical_pair_completion import _product, _t_phase
from tests.test_original_middle_pair_insertions import _ghz, _pair_coordinates
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _symplectic, _tau


def _projector(pauli, sign):
    return scale(add(eye(len(pauli)), scale(pauli, sign)), CO / 2)


def _signed_t(pauli, orientation):
    result = _t_phase(pauli)
    return result if orientation == 1 else adjoint(result)


class BooleanOriginalPairCollapseTests(unittest.TestCase):
    def test_actual_paired_boolean_identity_retains_unused_clean_action(self):
        # A three-wire, three-native-band source with exactly one pair
        # and no middle singleton. Its occupied action is literal I;
        # the unused clean block is deliberately not replaced by I.
        width, size = 3, 8
        k_word = _ghz(2)
        ea_word = [('T', 0), ('T', 1), ('TDG', 2)]
        ub_word = _pair_coordinates(2) + [('H', 2), ('CX', 0, 2), ('H', 2)]
        uc_word = ([('H', 2), ('S', 2), ('S', 2), ('H', 2)]
                   + [('CX', 0, 2), ('CX', 1, 2)])
        eb_word = _inverse(ub_word) + [('T', 0), ('T', 1)] + ub_word
        ec_word = _inverse(uc_word) + [('T', bit) for bit in range(3)] + uc_word
        word = k_word + ea_word + eb_word + _inverse(ec_word) + _inverse(k_word)
        k, ea, eb, ec, ub, uc = [native_columns(width, source) for source in
                                 (k_word, ea_word, eb_word, ec_word, ub_word, uc_word)]
        clean, pa = _tau(8, width), _tau(3, width)
        b, d = _tau(35, width), _tau(59, width)
        gamma = scale(_tau(24, width), -1)
        silent_a, silent_c = _tau(32, width), scale(_tau(56, width), -1)
        self.assertEqual(_conjugate(k, clean), pa)
        self.assertEqual(_conjugate(ub, _tau(8, width)), b)
        self.assertEqual(_conjugate(ub, _tau(16, width)), d)
        self.assertEqual(_conjugate(uc, _tau(32, width)), silent_c)
        self.assertEqual(multiply(b, d), gamma)
        self.assertEqual(_rank([8, 16, 32]), 3)
        self.assertEqual(_rank([8, 16, 56]), 3)
        self.assertEqual(_rank([35, 59]), 2)
        j = _conjugate(ea, pa)
        self.assertEqual(_conjugate(j, b), d)
        self.assertEqual(_conjugate(j, d), b)
        self.assertEqual(_conjugate(eb, j), _conjugate(ec, pa))

        actual = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(size))
        self.assertEqual(_conjugate(actual, clean), clean)
        self.assertEqual([len(band) for band in _native_t_bands(word)], [3, 2, 3])
        for column in range(size):
            phase = CI * (-1 if column & 4 else 1) if column & 3 == 3 else CO
            self.assertEqual([actual[row][column] for row in range(size)],
                             [phase if row == column else CZ for row in range(size)])

        # An apparent canceling negative singleton U=Z0 X2 commutes
        # the active axes and the pair, but fails on the ACTUAL silent
        # input axis. It therefore cannot be inserted into that band.
        proposed_u = 12
        self.assertEqual(_symplectic(proposed_u, 3, width), 1)
        for label in (8, 16, 35, 59):
            self.assertEqual(_symplectic(proposed_u, label, width), 0)
        self.assertEqual(_symplectic(proposed_u, 32, width), 1)
        self.assertEqual(_conjugate(_tau(proposed_u, width), silent_a),
                         scale(silent_a, -1))

    def test_silent_list_remains_independent_mod_clean_and_two_pair_centers(self):
        # Four actual active Z axes and two actual silent Z axes. The
        # clean Pauli is X0 X1 X2 X3; each pair difference is active-only.
        width = 6
        active = [1 << (width + bit) for bit in range(4)]
        silent = [1 << (width + bit) for bit in (4, 5)]
        clean = 15
        differences = [active[0] ^ active[1], active[2] ^ active[3]]
        constraints = [clean] + differences
        self.assertEqual(_rank(active + silent), 6)
        self.assertEqual(_rank(constraints), 3)
        self.assertEqual(_rank(silent + constraints), 5)
        self.assertEqual(_symplectic(active[0], clean, width), 1)
        for label in active + silent + differences:
            self.assertEqual(_symplectic(active[0], label, width), 0)
        # Thus the entire silent span intersects span(clean,r1,r2)
        # trivially, not merely each individual silent axis.
        self.assertEqual(_rank(silent + constraints) - _rank(constraints),
                         len(silent))

    def test_dependent_quotient_can_hide_nonempty_silent_coupling(self):
        # Deliberately violate the needed active/silent separation:
        # Gamma=Z1 lies in span(q1=Z0,q2=Z0 Z1). The physical list is
        # independent, but its Gamma=-1 quotient is not. On that sector
        # the two T rotations cancel up to their literal scalar.
        width = 2
        q1, q2 = _tau(1 << width, width), _tau(3 << width, width)
        gamma, h = _tau(2 << width, width), _tau(1, width)
        band = _product(_t_phase(q1), _t_phase(q2))
        dressed = _conjugate(adjoint(band), h)
        minus, plus = _projector(gamma, -1), _projector(gamma, 1)
        self.assertEqual(_rank([4, 12]), 2)
        self.assertEqual(_rank([4, 12, 8]) - _rank([8]), 1)
        for axis in (q1, q2):
            self.assertEqual(_conjugate(h, axis), scale(axis, -1))
        self.assertEqual(multiply(dressed, minus), multiply(h, minus))
        self.assertEqual(multiply(multiply(h, dressed), minus), minus)
        anticommutator = add(multiply(dressed, h), multiply(h, dressed))
        self.assertEqual(multiply(anticommutator, plus), scale(plus, 0))
        self.assertNotEqual(dressed, h)
        self.assertEqual(multiply(adjoint(band), band), eye(4))

    def test_independent_silent_axis_has_nonboolean_scalar_sector_derivative(self):
        # With no quotient relation, one nonempty active silent axis
        # keeps its odd S power. Its square is not identity even on the
        # entire Gamma=-1 block; clean return alone cannot remove it.
        width = 2
        q, gamma, h = (_tau(label, width) for label in (4, 8, 1))
        band = _t_phase(q)
        dressed = _conjugate(adjoint(band), h)
        derivative = multiply(h, dressed)
        minus = _projector(gamma, -1)
        self.assertEqual(_rank([4, 8]), 2)
        self.assertEqual(multiply(derivative, derivative), scale(q, -CI))
        self.assertNotEqual(multiply(multiply(derivative, derivative), minus),
                            minus)
        for column in (2, 3):
            self.assertNotIn(derivative[column][column], (CO, -CO))
        self.assertEqual(multiply(adjoint(band), band), eye(4))

    def test_two_literal_pairs_keep_every_signed_parity_branch_and_adjoint(self):
        # Two fixed commuting independent pairs. This is the literal
        # pair-factor identity, not a claim that this diagonal fixture
        # itself is an occupied Boolean original-pair source.
        width, size = 4, 16
        identity = eye(size)
        omega = native_columns(1, [('T', 0)])[1][1]
        b1, d1, b2 = (_tau(1 << (width + bit), width)
                       for bit in range(3))
        d2 = scale(_tau(1 << (width + 3), width), -1)
        pairs = [(b1, d1), (b2, d2)]
        t = _tau(3, width)  # Active on pair one only; preserves both parities.
        for settings in (((1, 1), (-1, -1)), ((-1, 1), (1, -1))):
            with self.subTest(settings=settings):
                factors, gammas, prefactors, branch_phases = [], [], [], []
                for (b, d), (swap_sign, epsilon) in zip(pairs, settings):
                    gamma = scale(multiply(b, d), swap_sign)
                    phase = omega if epsilon == 1 else omega.conj()
                    prefactor = CO if swap_sign == 1 else phase.conj()
                    factor = _product(_signed_t(b, epsilon),
                                      _signed_t(d, swap_sign * epsilon))
                    expected = scale(add(
                        multiply(_projector(gamma, 1), _quarter_turn(b, epsilon)),
                        scale(_projector(gamma, -1), phase)), prefactor)
                    self.assertEqual(factor, expected)
                    factors.append(factor)
                    gammas.append(gamma)
                    prefactors.append(prefactor)
                    branch_phases.append(phase)
                actual = _product(*factors)
                actual_inverse = _product(adjoint(factors[1]), adjoint(factors[0]))
                self.assertEqual(actual_inverse, adjoint(actual))
                self.assertEqual(multiply(actual_inverse, actual), identity)
                dressed = _conjugate(actual, t)
                for first_sign, second_sign in ((-1, -1), (-1, 1),
                                                (1, -1), (1, 1)):
                    sector = _product(_projector(gammas[0], first_sign),
                                      _projector(gammas[1], second_sign))
                    expected = identity
                    scalar = CO
                    for (b, _), (_, epsilon), sign, prefactor, phase in zip(
                            pairs, settings, (first_sign, second_sign),
                            prefactors, branch_phases):
                        scalar *= prefactor
                        if sign == -1:
                            scalar *= phase
                        else:
                            expected = multiply(expected, _quarter_turn(b, epsilon))
                    self.assertEqual(multiply(actual, sector),
                                     scale(multiply(expected, sector), scalar))
                    if first_sign == -1:
                        self.assertEqual(multiply(dressed, sector), multiply(t, sector))
                    else:
                        anticommutator = add(multiply(dressed, t), multiply(t, dressed))
                        self.assertEqual(multiply(anticommutator, sector),
                                         scale(sector, 0))

    def test_rank_two_five_negative_source_returns_clean_but_is_not_boolean(self):
        # The retained five-wire fixture: active ranks (2,5), one
        # original pair and three negative singletons. It satisfies
        # complete clean return, but even one occupied column is not
        # monomial. Therefore it does not contradict Boolean collapse.
        width = 5
        ka_word = _ghz(2)
        kc_word = ka_word + [gate for bit in (2, 3, 4)
                             for gate in [('H', bit), ('CX', 0, bit), ('H', bit)]]
        ub_word = _pair_coordinates(2) + [('CX', 0, bit) for bit in (2, 3, 4)]
        uc_word = ([('CX', bit, 0) for bit in (2, 3, 4)]
                   + [('CX', 0, bit) for bit in (2, 3, 4)]
                   + [('CX', bit, 0) for bit in (2, 3, 4)]
                   + [('CX', 0, bit) for bit in (2, 3, 4)])
        ea_word = [('T', 0), ('T', 1)]
        eb_word = _inverse(ub_word) + [('T', bit) for bit in range(width)] + ub_word
        ec_word = _inverse(uc_word) + [('T', bit) for bit in range(width)] + uc_word
        word = ka_word + ea_word + eb_word + _inverse(ec_word) + _inverse(kc_word)
        actual = native_columns(width, word)
        clean, pa = _tau(1 << width, width), _tau(3, width)
        j = _conjugate(native_columns(width, ea_word), pa)
        b_labels = [31, 31 | (3 << width), 5 << width, 9 << width, 17 << width]
        c_labels = [29 << width, 2 << width, 25 << width, 21 << width, 13 << width]
        pc_label = 3 | (28 << width)
        ub, uc = native_columns(width, ub_word), native_columns(width, uc_word)
        self.assertEqual(_conjugate(native_columns(width, ka_word), clean), pa)
        self.assertEqual(_conjugate(native_columns(width, kc_word), clean),
                         _tau(pc_label, width))
        for bit, (b_label, c_label) in enumerate(zip(b_labels, c_labels)):
            coordinate = _tau(1 << (width + bit), width)
            self.assertEqual(_conjugate(ub, coordinate), _tau(b_label, width))
            self.assertEqual(_conjugate(uc, coordinate), _tau(c_label, width))
            self.assertEqual(_symplectic(c_label, pc_label, width), 1)
        self.assertEqual(_rank([1 << width, 2 << width]), 2)
        self.assertEqual(_rank(b_labels), 5)
        self.assertEqual(_rank(c_labels), 5)
        self.assertEqual(_conjugate(j, _tau(b_labels[0], width)),
                         _tau(b_labels[1], width))
        for label in b_labels[2:]:
            self.assertEqual(_conjugate(j, _tau(label, width)),
                             scale(_tau(label, width), -1))
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(32))
        self.assertEqual(_conjugate(actual, clean), clean)
        self.assertEqual([len(band) for band in _native_t_bands(word)], [2, 5, 5])
        # Fixed occupied input c=0,a=1 has two nonzero output amplitudes.
        self.assertEqual(actual[2][2], (CO + CI) / 2)
        self.assertEqual(actual[30][2], (CO + CI) / 2)
        self.assertTrue(all(actual[row][2] == CZ for row in range(32)
                            if row not in (2, 30)))


if __name__ == '__main__':
    unittest.main()
