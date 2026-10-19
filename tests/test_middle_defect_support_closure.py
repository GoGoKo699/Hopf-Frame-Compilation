"""Finite checks for the middle-defect support argument.

The proof supplies the arbitrary-width conclusions. These checks separate
actual signed native permutation defects from a nonpermutation orthogonal
defect, and retain the Clifford square-root phase before taking a Boolean
cubic symbol. No center fiber is substituted for a full native word.
The independent native-list skeleton is rejected by its exact phase
defect; its one-core fiber is a falsification check, not a synthesis.
"""
from itertools import combinations
import unittest

from tests.exact_matrices import CZ, CO, adjoint, native_columns
from tests.test_bounded_relation_native_matching import (
    _native_phase_word, _phase, _rank,
)
from tests.test_two_layer_galois_return import (
    _clifford_columns, _conjugate, _inverse, _product, _symplectic, _tau,
)


def _wedge_linear(two_form, linear):
    """Exterior product over F2, with monomials stored as bit masks."""
    result = set()
    for pair in two_form:
        for bit in range(linear.bit_length()):
            axis = 1 << bit
            if linear & axis and not pair & axis:
                triple = pair | axis
                result.symmetric_difference_update((triple,))
    return result


def _polar_rank(two_form, width):
    columns = [0] * width
    for pair in two_form:
        first, second = [bit for bit in range(width) if pair >> bit & 1]
        columns[first] ^= 1 << second
        columns[second] ^= 1 << first
    return _rank(columns)


def _wedge_vectors(left, right, width):
    return {(1 << first) | (1 << second)
            for first, second in combinations(range(width), 2)
            if (((left >> first) & 1) * ((right >> second) & 1)
                ^ ((left >> second) & 1) * ((right >> first) & 1))}


class MiddleDefectSupportClosureTests(unittest.TestCase):
    def test_native_clean_gram_skeleton_still_fails_boolean_matching(self):
        # Nine diagonal coordinates a,c,p,z1,...,z6, plus spectator d.
        # This is a physical native-list/clean-return skeleton, NOT a
        # Boolean source. Fresh d repairs its logical Galois products;
        # a cubic coefficient and an exact one-core fiber still fail.
        a, c, p, d = 1, 2, 4, 512
        z = [1 << bit for bit in range(3, 9)]
        w = sum(z)
        h, v = a ^ c, p ^ z[0] ^ z[1]
        active_a = [a] + [a ^ axis for axis in z]
        active_c = [c] + [c ^ p ^ axis for axis in z]
        negative = [a ^ w, c ^ p ^ w, c, c ^ p]
        positive = [h ^ z[0], h ^ z[2], h ^ z[3],
                    h ^ p ^ z[1] ^ z[2] ^ z[3]]

        def gram(axes):
            return tuple(sum(((axis >> first) & 1) * ((axis >> second) & 1)
                             for axis in axes) % 2
                         for first in range(9) for second in range(9))

        self.assertEqual(tuple(left ^ right for left, right
                               in zip(gram(active_a), gram(active_c))), gram(negative))
        self.assertEqual(_rank(active_a), 7)
        self.assertEqual(_rank(active_c), 7)
        self.assertEqual(_rank(negative + positive), 8)
        self.assertEqual(_rank(positive + [v]), 4)
        for axes in (active_a, active_c, negative):
            self.assertEqual(_rank(axes + [v]), _rank(axes) + 1)
        # v is nevertheless clean-even and belongs to R=R_A+R_C.
        even_a = z
        even_c = [p ^ axis for axis in z]
        self.assertEqual(_rank(even_a + even_c + [v]), 7)

        width = 10
        zp = lambda axis: axis << width
        silent_a = [c, c ^ zp(v), zp(d)]
        # Original C axes, before removing C0=Z_c Z_p S_w CZ_(a,w).
        silent_c = [a ^ zp(w), a ^ zp(w ^ v), zp(d ^ v)]
        original_a = [zp(axis) for axis in active_a] + silent_a
        original_c = [zp(axis) for axis in active_c] + silent_c
        original_b = [zp(axis) for axis in negative + positive]
        for axes, rank in ((original_a, 10), (original_c, 10), (original_b, 8)):
            self.assertEqual(_rank(axes), rank)
            for left in axes:
                for right in axes:
                    self.assertEqual(_symplectic(left, right, width), 0)
        clean_a, clean_c = a ^ c, a ^ c ^ zp(w)
        for axis in silent_a:
            self.assertEqual(_symplectic(axis, clean_a, width), 0)
        for axis in silent_c:
            self.assertEqual(_symplectic(axis, clean_c, width), 0)
        for axis in negative:
            self.assertEqual((axis & (a ^ c)).bit_count() % 2, 1)
        for axis in positive:
            self.assertEqual((axis & (a ^ c)).bit_count() % 2, 0)

        # Exact phase identity on every assignment of all nine intrinsic
        # coordinates. No Hilbert-space matrix or initialized center is used.
        parity = lambda axis, state: (axis & state).bit_count() % 2
        for state in range(512):
            native_phase = (sum(parity(axis, state) for axis in active_a + negative)
                            - sum(parity(axis, state) for axis in active_c)) % 8
            bits = [parity(axis, state) for axis in z]
            quadratic = sum(bits[i] * bits[j] for i, j in combinations(range(6), 2)) % 2
            clifford_phase = (4 * parity(c ^ p, state) + 2 * parity(w, state)
                               + 4 * parity(a, state) * parity(w, state)) % 8
            self.assertEqual(native_phase,
                             (clifford_phase + 4 * parity(a ^ c ^ p, state)
                              * quadratic) % 8)
            flipped = state ^ a ^ c
            flipped_phase = (sum(parity(axis, flipped) for axis in active_a + negative)
                             - sum(parity(axis, flipped) for axis in active_c)) % 8
            self.assertEqual((flipped_phase - native_phase) % 8,
                             4 * (1 ^ parity(w, state)))

        # J_Zt's actual native displacement activates the two X_c axes.
        # It preserves the positive span but permutes none of its axes.
        def displacement(label):
            image = label
            for axis in original_a:
                if (_symplectic(axis, zp(h), width)
                        and _symplectic(axis, label, width)):
                    image ^= axis
            return image

        self.assertEqual(displacement(zp(h)), zp(h ^ v))
        for axis in positive:
            self.assertEqual(displacement(zp(axis)), zp(axis ^ v))
            self.assertNotIn(axis ^ v, positive)
        # Fresh d gives literal logical products A=v+d, C=d, B=v.
        total_a = total_c = 0
        for axis in silent_a:
            total_a ^= axis
        for axis in silent_c:
            total_c ^= axis
        self.assertEqual(total_a, zp(v ^ d))
        self.assertEqual(total_c, zp(d))
        total_positive = 0
        for axis in positive:
            total_positive ^= axis
        self.assertEqual(total_positive, v)
        self.assertEqual(total_a, total_c ^ zp(total_positive))

        # The middle defect and its shifted axes contain no z5 or z6.
        # The controlled-Pauli cubic has a nonzero p*z5*z6 coefficient.
        beta = {left | right for left, right in combinations(z, 2)}
        self.assertIn(p | z[4] | z[5], _wedge_linear(beta, v))
        for axis in positive + [axis ^ v for axis in positive]:
            self.assertEqual(axis & (z[4] | z[5]), 0)

        # Direct exact phase witness at p=1,z=0,d=0: E=omega I,
        # B=omega S_Z, and the two added safe rotations give omega^-1.
        # The core output is exactly S, whose two phases are not Boolean.
        tx = native_columns(1, [('H', 0), ('T', 0), ('H', 0)])
        omega = native_columns(1, [('T', 0)])[1][1]
        tminus_x = [[omega * value for value in row] for row in adjoint(tx)]
        e = _product(tx, tminus_x)
        self.assertEqual(e, [[omega if row == column else CZ for column in range(2)]
                             for row in range(2)])
        s_gate = native_columns(1, [('S', 0)])
        b_fiber = [[omega * value for value in row] for row in s_gate]
        occupied = _product(adjoint(e), b_fiber, e)
        occupied = [[omega.conj() * value for value in row] for row in occupied]
        self.assertEqual(occupied, s_gate)
        self.assertNotIn(occupied[1][1], (CO, -CO))

    def test_small_exterior_support_and_two_control_directions(self):
        # All 1,024 alternating forms in dimension five, with every
        # coordinate overlap dimension and every possible linear control.
        # This is a bounded falsification fixture, not a proof by search.
        width = 5
        pairs = [sum(1 << bit for bit in pair)
                 for pair in combinations(range(width), 2)]
        mixed_cases = 0
        for encoding in range(1 << len(pairs)):
            beta = {pair for index, pair in enumerate(pairs)
                    if encoding >> index & 1}
            rank = _polar_rank(beta, width)
            if rank < 4:
                continue
            for dimension in range(width + 1):
                center = (1 << dimension) - 1
                allowed = []
                for linear in range(1, 1 << width):
                    cubic = _wedge_linear(beta, linear)
                    if any(term & ~center for term in cubic):
                        continue
                    allowed.append(linear)
                    self.assertEqual(linear & ~center, 0)
                    # beta is in Lambda^2(C) + linear wedge R.
                    transverse = {term for term in beta if term & ~center}
                    self.assertTrue(any(
                        {term for term in _wedge_vectors(linear, other, width)
                         if term & ~center} == transverse
                        for other in range(1 << width)
                    ))
                    self.assertLessEqual(rank, dimension + 2)
                if _rank(allowed) >= 2:
                    self.assertTrue(all(term & ~center == 0 for term in beta))
                if allowed and any(term & ~center for term in beta):
                    mixed_cases += 1
                    self.assertEqual(_rank(allowed), 1)
        self.assertGreater(mixed_cases, 0)

    def test_literal_signed_permutation_defects_are_clifford(self):
        width = 3
        axes = (1, 2, 4)
        # J = X0 X1 SWAP01 is an actual involution, sending the first
        # two native axes to each other's negatives and fixing the third.
        j_word = [('CX', 0, 1), ('CX', 1, 0), ('CX', 0, 1)]
        for bit in (0, 1):
            j_word += [('H', bit), ('S', bit), ('S', bit), ('H', bit)]
        j = native_columns(width, j_word)
        self.assertEqual(j, adjoint(j))
        self.assertEqual(_conjugate(j, _tau(8, width)),
                         [[-value for value in row] for row in _tau(16, width)])
        self.assertEqual(_conjugate(j, _tau(16, width)),
                         [[-value for value in row] for row in _tau(8, width)])
        omega = native_columns(1, [('T', 0)])[1][1]
        roots = [CO]
        for _ in range(7):
            roots.append(roots[-1] * omega)
        for signs, daggers in ((0, 0), (1, 2), (5, 6)):
            b_word = _native_phase_word(axes, signs, daggers)
            b = native_columns(width, b_word)
            defect_word = _inverse(j_word) + _inverse(b_word) + j_word + b_word
            defect = native_columns(width, defect_word)
            self.assertEqual(defect, _product(b, _conjugate(j, adjoint(b))))
            phase = _phase(axes, width, signs, daggers)
            image = [((state & 4) | ((state & 1) << 1)
                      | ((state & 2) >> 1)) ^ 3 for state in range(8)]
            expected = [[roots[(phase[column] - phase[image[column]]) % 8]
                         if row == column else CZ for column in range(8)]
                        for row in range(8)]
            self.assertEqual(defect, expected)
            _clifford_columns(defect, width)
            self.assertEqual(native_columns(width, _inverse(defect_word)),
                             adjoint(defect))

    def test_independent_orthogonal_defect_needs_boolean_normalization(self):
        width = 4
        # O = I + 11^T has O^2=I and O^T O=I. It takes each basis
        # axis to the other three, so none is a native permutation orbit.
        j_word = [('CX', 3, 2), ('CX', 2, 3), ('CX', 2, 1),
                  ('CX', 2, 0), ('CX', 3, 2), ('CX', 1, 0),
                  ('CX', 0, 3), ('CX', 0, 2), ('CX', 0, 1),
                  ('CX', 1, 0)]
        j = native_columns(width, j_word)
        self.assertEqual(j, adjoint(j))
        for bit in range(width):
            self.assertEqual(_conjugate(j, _tau(1 << (width + bit), width)),
                             _tau((15 ^ (1 << bit)) << width, width))
        b_word = [('T', bit) for bit in range(width)]
        defect_word = _inverse(j_word) + _inverse(b_word) + j_word + b_word
        defect = native_columns(width, defect_word)
        parity = [('CX', bit, 3) for bit in range(3)]
        square_root_word = parity + [('S', 3)] + list(reversed(parity))
        square_root = native_columns(width, square_root_word)
        self.assertEqual(_product(defect, defect),
                         _product(square_root, square_root))
        boolean_word = defect_word + _inverse(square_root_word)
        boolean = native_columns(width, boolean_word)
        expected = [[CZ for _ in range(16)] for _ in range(16)]
        for state in range(16):
            bits = [(state >> bit) & 1 for bit in range(width)]
            h = (sum(bits) + sum(bits[i] * bits[k] * bits[l]
                                for i, k, l in combinations(range(width), 3))) % 2
            expected[state][state] = -CO if h else CO
        self.assertEqual(boolean, expected)
        self.assertEqual(native_columns(width, _inverse(boolean_word)), adjoint(boolean))
        with self.assertRaises(ValueError):
            _clifford_columns(boolean, width)


if __name__ == '__main__':
    unittest.main()
