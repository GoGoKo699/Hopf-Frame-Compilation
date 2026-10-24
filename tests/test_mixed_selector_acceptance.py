"""Exact common-fixed-vector obstruction for three cell-CCZ dressings.

On x=0,Y=e_0 every literal matching stage fixes H0={s_0=0}.
For a cell product J, its origin-cell factor is a conjugate of CCZ x I,
so rank(J_cell-I)=2.  Consequently rank(P1 J P0)<=D/8.  Three such
off-diagonal blocks leave a common kernel of dimension at least D/8
inside H0, fixed by all six branches R_j and J_j R_j J_j.

This is a specified source-grammar obstruction, not a native-depth lower
bound.  Arbitrary complex path coefficients relax the coefficients of a
native one-signal circuit; the relaxation can prove failure, not synthesis.
The exact complete p=2,3 fixtures retain literal stage phases and dirty
wire return.  A local product witness also explains the scalable example.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product
import unittest

import numpy as np

from tests.exact_matrices import CO, CQ2, CZ, Q2, native_columns
from tests.test_co_moving_lookup_stage import _Matching
from tests.test_entangling_selector_readout import _binary_image, _zeta_word
from tests.test_t_depth import _ccz_batch, _word


# (omitted corner, Hadamard mask, X-conjugation mask, Z-conjugation mask).
# The same four-row-cell frame is repeated across the selector register.
_FRAMES = ((1, 1, 1, 1), (2, 15, 2, 1), (3, 3, 4, 5))


@lru_cache(None)
def _cell_frame(omitted, hmask, xmask=0, zmask=0):
    """Exact numerator / 2^exponent of (XZ) H P H (XZ)^dagger."""
    triple = 15 ^ (1 << omitted)
    hmask &= triple  # Hadamards on the omitted wire cancel.
    submasks = [s for s in range(16) if not s & ~hmask]
    exponent = hmask.bit_count()
    matrix = np.zeros((16, 16), dtype=np.int64)
    for column, delta in product(range(16), submasks):
        outside = (column ^ xmask) & ~hmask
        total = 0
        for inner in submasks:
            state = outside | inner
            parity = ((state & triple) == triple)
            parity += (state & delta).bit_count() + (zmask & delta).bit_count()
            total += (-1) ** parity
        matrix[column ^ delta, column] = total
    while exponent and np.all(matrix % 2 == 0):
        matrix //= 2
        exponent -= 1
    return matrix, exponent


def _frame(p, spec):
    local, local_exponent = _cell_frame(*spec)
    matrix = np.ones((1, 1), dtype=np.int64)
    for _ in range((1 << p) // 4):
        matrix = np.kron(local, matrix)
    return matrix, local_exponent * ((1 << p) // 4)


def _exact_rank(matrix):
    """Rational elimination for these small exact audit maps."""
    rows = [[Fraction(int(value)) for value in row] for row in matrix]
    rank = 0
    for column in range(matrix.shape[1]):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [value / scale for value in rows[rank]]
        for i in range(rank + 1, len(rows)):
            factor = rows[i][column]
            if factor:
                rows[i] = [a - factor*b for a, b in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def _native_cell_word(spec):
    omitted, hmask, xmask, zmask = spec
    xs = [gate for wire in range(4) if (xmask >> wire) & 1
          for gate in [('H', wire), ('S', wire), ('S', wire), ('H', wire)]]
    zs = [('S', wire) for wire in range(4) if (zmask >> wire) & 1 for _ in range(2)]
    hs = [('H', wire) for wire in range(4) if (hmask >> wire) & 1]
    a, b, c = [wire for wire in range(4) if wire != omitted]
    return xs + zs + hs + _word(_ccz_batch(a, [(b, c)])) + hs + zs + xs


def _literal_origin_stage(matching, inverse=False):
    dimension = 1 << matching.m
    matrix = np.zeros((dimension, dimension), dtype=np.int64)
    for selector, borrowed, phase in product(range(dimension), (0, 1), (0, 3)):
        state = (0, 1, selector, borrowed, phase)
        actual = matching.stage_literal(state, inverse=inverse)
        assert actual[:2] == (0, 1) and actual[3:] == (borrowed, phase)
        assert matching.stage_literal(actual, inverse=not inverse) == state
        matrix[actual[2], selector] = 1
    return matrix


class MixedSelectorAcceptanceTests(unittest.TestCase):
    def test_local_rank_two_defect_and_off_diagonal_rank_bound(self):
        # All omitted corners and partial Fourier masks, with several affine
        # Pauli shifts.  Rank invariance proves the statement beyond fixtures.
        for omitted, hmask, (xmask, zmask) in product(
                range(4), range(16), ((0, 0), (5, 10), (15, 15))):
            matrix, exponent = _cell_frame(omitted, hmask, xmask, zmask)
            identity = np.eye(16, dtype=np.int64) * (1 << exponent)
            self.assertTrue(np.array_equal(matrix, matrix.T))
            self.assertTrue(np.array_equal(matrix @ matrix, identity @ identity))
            self.assertEqual(_exact_rank(matrix - identity), 2)
            self.assertLessEqual(_exact_rank(matrix[1::2, 0::2]), 2)

    def test_cell_native_words_preserve_literal_phases_on_every_column(self):
        for spec in _FRAMES:
            matrix, exponent = _cell_frame(*spec)
            expected = [[CQ2(Q2(Fraction(int(value), 1 << exponent)))
                         for value in row] for row in matrix]
            self.assertEqual(native_columns(4, _native_cell_word(spec)), expected)

    def test_two_and_three_bit_factorization_pays_the_whole_off_diagonal_rank(self):
        for p in (2, 3):
            dimension = 1 << (1 << p)
            ranks = []
            for spec in _FRAMES:
                local, exponent = _cell_frame(*spec)
                matrix, _ = _frame(p, spec)
                rest = (np.ones((1, 1), dtype=np.int64) if p == 2 else local)
                # Other-cell factors are invertible; the Kronecker rank is
                # rank(B_origin) times their full dimension D/16.
                self.assertTrue(np.array_equal(matrix[1::2, 0::2],
                                               np.kron(rest, local[1::2, 0::2])))
                self.assertTrue(np.array_equal(local @ local,
                                               np.eye(16, dtype=np.int64)*(1 << (2*exponent))))
                ranks.append(_exact_rank(local[1::2, 0::2]) * (dimension // 16))
            self.assertTrue(all(0 < rank <= dimension // 8 for rank in ranks))
            self.assertGreaterEqual(dimension // 2 - sum(ranks), dimension // 8)
            self.assertEqual(dimension // 2 - 4*(dimension // 8), 0)

    def test_tiny_product_witness_and_actual_six_branch_fixed_vector(self):
        local_witness = np.zeros(16, dtype=np.int64)
        local_witness[[0, 2]] = 1  # |+> on selector wire 1, zero elsewhere.
        for spec in _FRAMES:
            local, exponent = _cell_frame(*spec)
            self.assertTrue(np.array_equal(local @ local_witness,
                                           local_witness * (1 << exponent)))
            self.assertTrue(np.any(local[1::2, 0::2]))  # J does not preserve H0.
        first, exponent = _cell_frame(*_FRAMES[0])
        self.assertFalse(np.array_equal(first @ np.ones(16, dtype=np.int64),
                                        np.ones(16, dtype=np.int64)*(1 << exponent)))
        for p in (2, 3):
            pairs = (((0, 1),), ((1, 0),), ((0, 1),)) if p == 2 else (
                ((0, 1),), ((1, 2),), ((2, 0),))
            witness = np.zeros(1 << (1 << p), dtype=np.int64)
            witness[[0, 2]] = 1
            for spec, pair in zip(_FRAMES, pairs):
                frame, exponent = _frame(p, spec)
                moved = frame @ witness
                self.assertFalse(np.any(moved[1::2]))
                for inverse in (False, True):
                    stage = _literal_origin_stage(_Matching(p, pair), inverse)
                    self.assertTrue(np.array_equal(stage @ witness, witness))
                    self.assertTrue(np.array_equal(frame @ stage @ moved,
                                                   witness * (1 << (2*exponent))))

    def test_empty_bank_all_branches_identity_and_relaxed_coefficient_floor(self):
        for p in (2, 3):
            pairs = (((0, 1),), ((1, 0),)) if p == 2 else (
                ((0, 1),), ((1, 2),), ((2, 0),))
            for pair in pairs:
                matching = _Matching(p, pair)
                for x, selector, borrowed in product(range(1 << p),
                        range(1 << matching.m), (0, 1)):
                    state = (x, 0, selector, borrowed, 3)
                    self.assertEqual(matching.stage_literal(state), state)
                    self.assertEqual(matching.stage_literal(state, inverse=True), state)
        # J^2=I above then gives every empty-bank branch/word exactly I.
        # On both the empty fiber and the normalized marked witness, K has
        # scalar sum(c).  Thus max(q,|q-1/2|)>=1/4 for q=|sum(c)|^2.
        # These arbitrary rational complex coefficients need not be realizable
        # by one-signal native interlayers; using the superset is intentional.
        for coefficients in ([(Fraction(1, 16), Fraction(0))]*8,
                             [(Fraction((-1)**k, 8), Fraction(1, 16)) for k in range(8)],
                             [(Fraction(0), Fraction(0))]*7 + [(Fraction(1, 4), Fraction(1, 2))]):
            q = sum(real for real, _ in coefficients)**2 + sum(imag for _, imag in coefficients)**2
            self.assertGreaterEqual(max(q, abs(q-Fraction(1, 2))), Fraction(1, 4))
        self.assertEqual(max(Fraction(1, 4), abs(Fraction(1, 4)-Fraction(1, 2))),
                         Fraction(1, 4))

    def test_transpose_zeta_wrapped_candidate_has_a_different_exact_witness(self):
        # N=M_3^T is the actual chronological CX word; use its reversed word
        # for N^dagger even though its complete binary map is involutive.
        word = _zeta_word(3, transpose=True)
        inverse_word = list(reversed(word))
        images = np.array([_binary_image(word, selector) for selector in range(256)])
        self.assertEqual([_binary_image(inverse_word, int(image)) for image in images],
                         list(range(256)))
        expected = [[CO if row == images[column] else CZ for column in range(256)]
                    for row in range(256)]
        self.assertEqual(native_columns(8, word), expected)
        self.assertEqual(native_columns(8, inverse_word), list(map(list, zip(*expected))))
        preimages = np.argsort(images)
        witness = np.zeros(256, dtype=np.int64)
        witness[[68, 70]] = 1  # s_1 is |+>, s_2=s_6=1, all other wires zero.
        for spec, pairs in zip(_FRAMES, (((0, 1),), ((1, 2),), ((2, 0),))):
            frame, exponent = _frame(3, spec)
            # J'=N [(XZ) H P H (XZ)^dagger] N^dagger: the Pauli words
            # remain conjugations, not independently multiplied interlayers.
            wrapped = frame[np.ix_(preimages, preimages)]
            self.assertEqual(_exact_rank(wrapped[1::2, 0::2]), 48)
            self.assertGreater(48, 256 // 8)  # Escapes the cell-local budget.
            for inverse in (False, True):
                stage = _literal_origin_stage(_Matching(3, pairs), inverse)
                self.assertTrue(np.array_equal(stage @ witness, witness))
                self.assertTrue(np.array_equal(wrapped @ stage @ wrapped @ witness,
                                               witness * (1 << (2*exponent))))


if __name__ == '__main__':
    unittest.main()
