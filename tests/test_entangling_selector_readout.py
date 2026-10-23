"""Exact controlled selector-zeta and partial-Fourier component audits.

All control/selector columns, literal phases and actual inverses are kept in
small native fixtures. Larger binary checks certify complete linear maps.
These components do not supply a complete lookup source or a depth bound.
"""
import unittest
from fractions import Fraction

from tests.exact_matrices import CO, CZ, CQ2, Q2, adjoint, native_columns
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts
from tests.test_t_depth import _ccz_batch
from tests.test_two_layer_diagonal_rigidity import _reduce_adjacent_s


def _zeta_word(p, offset=0, transpose=False):
    """Classical oriented Boolean-zeta transform, in chronological order."""
    word = [('CX', offset + u, offset + (u | (1 << j)))
            for j in range(p) for u in range(1 << p) if not u & (1 << j)]
    return [('CX', b, a) for _, a, b in word] if transpose else word


def _conjugator_word(p, offset=0, transpose=False):
    """P=Q diag(I,M_{p-1}) Q, with P=P^{-1}."""
    q = 1 << (p - 1)
    qword = [('CX', offset + q + i, offset + i) for i in range(q)]
    word = qword + _zeta_word(p - 1, offset + q) + qword
    return [('CX', b, a) for _, a, b in word] if transpose else word


def _controlled_zeta_schedule(p, transpose=False, partial_h=False):
    """Wire 0 is an arbitrary control, wires 1..2^p are arbitrary selectors."""
    q = 1 << (p - 1)
    pairs = [(1 + i, 1 + q + i) for i in range(q)]
    if transpose:
        pairs = [(b, a) for a, b in pairs]
    basis = _conjugator_word(p, 1, transpose)
    walls = [('H', b) for _, b in pairs]
    partial = [('H', 1+i) for i in range(q)] if partial_h else []
    batch = _ccz_batch(0, pairs)
    # The retained emitter displays T_control^q as S^(q//2); reduce only
    # adjacent S powers modulo four before charging the common bound.
    batch = [(kind, _reduce_adjacent_s(word)) for kind, word in batch]
    return ([('C', partial + basis + walls)] + batch
            + [('C', walls + _inverse(basis) + partial)])


def _flatten(schedule):
    return [gate for _, word in schedule for gate in word]


def _binary_image(word, state):
    for name, control, target in word:
        assert name == 'CX'
        state ^= ((state >> control) & 1) << target
    return state


def _zeta_image(state, p, transpose=False):
    m = 1 << p
    result = 0
    for v in range(m):
        bit = 0
        for u in range(m):
            if (v & u) == (v if transpose else u):
                bit ^= (state >> u) & 1
        result |= bit << v
    return result


def _binary_rank(columns):
    basis = {}
    for column in columns:
        while column:
            pivot = column.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = column
                break
            column ^= basis[pivot]
    return len(basis)


class SelectorZetaNativeTests(unittest.TestCase):
    def test_three_bit_selected_native_columns_and_actual_inverse(self):
        # Nine physical wires: selected columns keep this native check small.
        # All columns remain covered at p=1,2, and the binary identities below
        # are complete linear-map checks, not sampling of selector strings.
        p, m = 3, 8
        selectors = (0, 1, 3, 0x55, 0xAA, 0xFF)
        columns = [(s << 1) | control for s in selectors for control in (0, 1)]
        for transpose in (False, True):
            word = _flatten(_controlled_zeta_schedule(p, transpose))
            expected = [[CZ for _ in columns] for _ in range(1 << (m+1))]
            for j, column in enumerate(columns):
                image = ((_zeta_image(column >> 1, p, transpose) << 1) | 1
                         if column & 1 else column)
                expected[image][j] = CO
            for orientation in (word, _inverse(word)):
                self.assertEqual(native_columns(m+1, orientation, columns), expected)

    def test_partial_h_zeta_dense_literal_kernel_and_actual_inverse(self):
        # C=H_lower M H_lower; its exact complete kernel is dense and real.
        # This breaks the selector-zero premise, unlike global H M H.
        for p in (1, 2):
            m, q = 1 << p, 1 << (p-1)
            halfmask = (1 << q)-1
            factor = CQ2(Q2(Fraction(1, 1 << q)))
            word = _flatten(_controlled_zeta_schedule(p, partial_h=True))
            expected = [[CZ for _ in range(1 << (m+1))]
                        for _ in range(1 << (m+1))]
            for column in range(1 << (m+1)):
                if not column & 1:
                    expected[column][column] = CO
                    continue
                a, b = (column >> 1) & halfmask, column >> (q+1)
                ab = _zeta_image(b, p-1)
                for row in range(1, 1 << (m+1), 2):
                    v, w = (row >> 1) & halfmask, row >> (q+1)
                    aw = _zeta_image(w, p-1)
                    parity = ((a & aw).bit_count() + (a & b).bit_count()
                              + (v & w).bit_count() + (v & ab).bit_count()) % 2
                    expected[row][column] = -factor if parity else factor
            actual = native_columns(m+1, word)
            self.assertEqual(actual, expected)
            self.assertEqual(native_columns(m+1, _inverse(word)), adjoint(expected))
            self.assertTrue(all(expected[row][1] != CZ for row in range(1, 1 << (m+1), 2)))
            for orientation in (word, _inverse(word)):
                t_count, clifford_count = _counts(orientation)
                self.assertEqual(t_count, 3*m+q%2)
                self.assertLessEqual(clifford_count, (p+17)*m//2+3)

    def test_linear_conjugacy_involution_rank_and_literal_zeta(self):
        for p in range(1, 8):
            m, q = 1 << p, 1 << (p - 1)
            for transpose in (False, True):
                zeta = _zeta_word(p, transpose=transpose)
                basis = _conjugator_word(p, transpose=transpose)
                jword = [('CX', i, q+i) for i in range(q)]
                if transpose:
                    jword = [('CX', b, a) for _, a, b in jword]
                self.assertEqual(len(zeta), p*m//2)
                self.assertEqual(len(basis), m+(p-1)*m//4)
                columns = []
                for wire in range(m):
                    state = 1 << wire
                    image = _binary_image(zeta, state)
                    self.assertEqual(image, _zeta_image(state, p, transpose))
                    self.assertEqual(_binary_image(zeta, image), state)
                    self.assertEqual(_binary_image(basis+basis, state), state)
                    self.assertEqual(_binary_image(basis+jword+_inverse(basis), state), image)
                    columns.append(image ^ state)
                self.assertEqual(_binary_rank(columns), q)

    def test_native_forward_and_actual_inverse_all_columns(self):
        for p in (1, 2):
            m = 1 << p
            for transpose in (False, True):
                word = _flatten(_controlled_zeta_schedule(p, transpose))
                expected = [[CZ for _ in range(1 << (m+1))]
                            for _ in range(1 << (m+1))]
                for column in range(1 << (m+1)):
                    image = ((_zeta_image(column >> 1, p, transpose) << 1) | 1
                             if column & 1 else column)
                    expected[image][column] = CO
                actual = native_columns(m+1, word)
                self.assertEqual(actual, expected)
                self.assertEqual(native_columns(m+1, _inverse(word)), adjoint(expected))

    def test_positive_s_ledgers_and_disjoint_t_targets(self):
        for p in range(1, 10):
            m, q = 1 << p, 1 << (p-1)
            for transpose in (False, True):
                schedule = _controlled_zeta_schedule(p, transpose)
                self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 4)
                for kind, gates in schedule:
                    self.assertTrue(all(gate[0] in ('CX', 'H', 'S', 'T', 'TDG')
                                        for gate in gates))
                    if kind == 'T':
                        wires = [gate[1] for gate in gates]
                        self.assertEqual(len(wires), len(set(wires)))
                        if p >= 2:
                            self.assertNotIn(0, wires)
                for word in (_flatten(schedule), _inverse(_flatten(schedule))):
                    t_count, clifford_count = _counts(word)
                    self.assertEqual(t_count, 3*m + q % 2)
                    self.assertLessEqual(clifford_count, (p+15)*m//2 + 3)


if __name__ == '__main__':
    unittest.main()
