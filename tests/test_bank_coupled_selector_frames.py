"""Exact Bell-cell witness for the specified bank-coupled selector frames.

The arbitrary second selector bank t belongs to the complete Gram domain;
it need not return after a source block.  Choosing t0=0 below supplies a
legal counterexample input, never a promised initialization.  Implementation
helpers of each literal stage do return.  The crossing-code witness survives
the asymmetric branch V^dagger B_Y R_C V, despite failure of an older witness.
A separate full-column Gram certificate rejects the specified bank-H extension.
"""
from itertools import product
from fractions import Fraction
import unittest

import numpy as np

from tests.exact_matrices import CO, CQ2, CZ, Q2, native_columns
from tests.test_co_moving_lookup_stage import _Matching
from tests.test_entangling_selector_readout import _binary_image, _binary_rank, _zeta_word
from tests.test_mixed_selector_acceptance import (
    _FRAMES, _cell_frame, _frame, _literal_origin_stage, _native_cell_word,
)
from tests.test_t_depth import _fredkin_batch, _inverse, _word


def _code_columns(p):
    mask = sum(1 << u for u in range(1 << p) if u % 4 == 1)
    return [s for s in range(1 << (1 << p)) if not s & ~mask]


def _wrapped_frame(p, spec):
    frame, exponent = _frame(p, spec)
    word = _zeta_word(p, transpose=True)
    images = np.array([_binary_image(word, s) for s in range(len(frame))])
    return frame[np.ix_(np.argsort(images), np.argsort(images))], exponent


def _swap_images(p):
    """Origin bank e0 activates only SWAP(s0,t0); t0 is effective wire m."""
    m = 1 << p
    return np.array([v ^ ((((v >> m) ^ v) & 1) * ((1 << m) | 1))
                     for v in range(1 << (m+1))])


def _on_selectors(matrix, vector):
    # Both t0 blocks are retained; every other t wire is an identity factor.
    return (vector.reshape(2, len(matrix)) @ matrix.T).reshape(-1)


def _asymmetric_branch(vector, frame, stage, swap, inverse=False):
    vector = _on_selectors(frame, vector)
    if inverse:
        vector = vector[swap]
        vector = _on_selectors(stage.T, vector)
    else:
        vector = _on_selectors(stage, vector)
        vector = vector[swap]
    return _on_selectors(frame.T, vector)


def _h0_numerator(values):
    """H on selector0, with its single factor 1/sqrt(2) kept separate."""
    result = values.copy()
    result[0::2] = values[0::2] + values[1::2]
    result[1::2] = values[0::2] - values[1::2]
    return result


def _bank_h_branch(vector, frame, stage, swap, inverse=False):
    vector = _on_selectors(frame, vector)
    if inverse:
        vector = _h0_numerator(vector)[swap]
        vector = _on_selectors(stage.T, vector)
    else:
        vector = _h0_numerator(_on_selectors(stage, vector)[swap])
    return _on_selectors(frame.T, vector)


def _remove_twos(matrix, exponent):
    while exponent and np.all(matrix % 2 == 0):
        matrix //= 2
        exponent -= 1
    return matrix, exponent


def _complete_bank_h_gram():
    """All 512 columns and eight paths; exact G=(real+sqrt(2)radical)/den.

    Each path is W_b=A_3,b3 A_2,b2 A_1,b1, in lexicographic bit order.
    Branch1 chronology is V,R,B,H0,V^dagger on the origin singleton.
    Shared prefixes, permutation rows and Hadamard butterflies avoid dense
    512-by-512 matrix multiplication. All t0 inputs are included.
    """
    swap = _swap_images(3)
    paths = [(np.eye(512, dtype=np.int64), 0, 0)]
    exponent_budget = 0
    for spec, pair in zip(_FRAMES, (((0, 1),), ((1, 2),), ((2, 0),))):
        frame, exponent = _wrapped_frame(3, spec)
        assert np.array_equal(frame, frame.T)
        exponent_budget += 2*exponent
        stage = _literal_origin_stage(_Matching(3, pair))
        preimages = np.argmax(stage, axis=1)
        rows = np.concatenate((preimages, preimages+256))
        def apply_frame(matrix):
            return np.vstack((frame @ matrix[:256], frame @ matrix[256:]))
        following = []
        for matrix, scale, count in paths:
            following.append((matrix[rows], scale, count))
            moved = apply_frame(matrix)[rows][swap]
            moved = apply_frame(_h0_numerator(moved))
            moved, new_scale = _remove_twos(moved, scale+2*exponent)
            following.append((moved, new_scale, count+1))
        paths = following
    # A prefix's squared Frobenius norm is at most this integer. It also
    # bounds absolute dot-product accumulation, keeping int64 exact.
    assert 512*(1 << (2*exponent_budget+3)) < (1 << 63)
    words = []
    for matrix, exponent, count in paths:
        matrix, exponent = _remove_twos(matrix, exponent+(count+1)//2)
        words.append((matrix, exponent, count % 2))
    scale = max(exponent for _, exponent, _ in words)
    # Also bound the common-denominator dot products and the factor of two
    # used for the real Gram entries when both paths carry sqrt(2).
    assert 1024*(1 << (2*scale)) < (1 << 63)
    flat = np.stack([(matrix*(1 << (scale-exponent))).reshape(-1)
                     for matrix, exponent, _ in words])
    dots = flat @ flat.T
    real, radical = np.zeros((8, 8), dtype=np.int64), np.zeros((8, 8), dtype=np.int64)
    for i, j in product(range(8), repeat=2):
        left, right = words[i][2], words[j][2]
        if left != right:
            radical[i, j] = dots[i, j]
        else:
            real[i, j] = dots[i, j]*(2 if left else 1)
    return real, radical, 512*(1 << (2*scale))


class BankCoupledSelectorFrameTests(unittest.TestCase):
    def test_bank_hadamard_normalization_and_actual_branch_inverse(self):
        root_half = CQ2(Q2(Fraction(0), Fraction(1, 2)))
        self.assertEqual(native_columns(1, [('H', 0)]),
                         [[root_half, root_half], [root_half, -root_half]])
        self.assertEqual(native_columns(1, [('H', 0), ('H', 0)]), [[CO, CZ], [CZ, CO]])
        swap = _swap_images(3)
        for spec, pair in zip(_FRAMES, (((0, 1),), ((1, 2),), ((2, 0),))):
            frame, exponent = _wrapped_frame(3, spec)
            stage = _literal_origin_stage(_Matching(3, pair))
            for selector, t0 in product((0, 1, 34, 70, 255), (0, 1)):
                vector = np.zeros(512, dtype=np.int64)
                vector[selector+256*t0] = 1
                output = _bank_h_branch(vector, frame, stage, swap)
                returned = _bank_h_branch(output, frame, stage, swap, inverse=True)
                self.assertTrue(np.array_equal(returned, vector*(1 << (4*exponent+1))))

    def test_every_local_native_frame_fixes_the_bell_cell_with_literal_phase(self):
        bell = np.zeros(16, dtype=np.int64)
        bell[[0, 3]] = 1  # |Phi+> on corners0,1; corners2,3 are zero.
        for spec in _FRAMES:
            frame, exponent = _cell_frame(*spec)
            self.assertTrue(np.array_equal(frame @ bell, bell * (1 << exponent)))
            columns = native_columns(4, _native_cell_word(spec), [0, 3])
            self.assertEqual([row[0] + row[1] for row in columns],
                             [CO if row in (0, 3) else CZ for row in range(16)])

    def test_native_bank_swap_and_actual_reverse_on_all_control_dirty_columns(self):
        schedule = _fredkin_batch(0, [(1, 2)])
        expected = [[CZ for _ in range(8)] for _ in range(8)]
        for column in range(8):
            image = column
            if column & 1:
                image ^= (((column >> 1) ^ (column >> 2)) & 1) * 6
            expected[image][column] = CO
        self.assertEqual(native_columns(3, _word(schedule)), expected)
        self.assertEqual(native_columns(3, _word(_inverse(schedule))), expected)

    def test_transpose_zeta_sends_uniform_corner_code_to_product_bell_cells(self):
        for p in (2, 3):
            m = 1 << p
            columns = _code_columns(p)
            word = _zeta_word(p, transpose=True)
            images = [_binary_image(word, column) for column in columns]
            bell_columns = {sum(3 << (4*cell) for cell in range(m//4)
                                if (selection >> cell) & 1)
                            for selection in range(1 << (m//4))}
            self.assertEqual(set(images), bell_columns)
            native = native_columns(m, word, columns)
            self.assertEqual([sum(row, CZ) for row in native],
                             [CO if row in bell_columns else CZ for row in range(1 << m)])
        # Sparse p=4 generator checks expose the all-size linear argument:
        # higher-coordinate zeta is an invertible basis change within the
        # corner1 code, while low-coordinate zeta produces Bell pairs.
        p, cells = 4, 4
        word = _zeta_word(p, transpose=True)
        images = [_binary_image(word, 1 << (4*cell+1)) for cell in range(cells)]
        expected = [sum(3 << (4*lower) for lower in range(cells)
                        if lower & cell == lower) for cell in range(cells)]
        self.assertEqual(images, expected)
        self.assertEqual(_binary_rank(images), cells)

    def test_wrapped_frames_fix_the_new_code_and_preserve_actual_stage_adjoints(self):
        for p in (2, 3):
            dimension = 1 << (1 << p)
            psi = np.zeros(dimension, dtype=np.int64)
            psi[_code_columns(p)] = 1
            swap = _swap_images(p)
            self.assertTrue(np.array_equal(swap[swap], np.arange(2*dimension)))
            pairs = (((0, 1),), ((1, 0),), ((0, 1),)) if p == 2 else (
                ((0, 1),), ((1, 2),), ((2, 0),))
            for spec, pair in zip(_FRAMES, pairs):
                frame, exponent = _wrapped_frame(p, spec)
                self.assertTrue(np.array_equal(frame @ psi, psi * (1 << exponent)))
                stage = _literal_origin_stage(_Matching(p, pair))
                reverse = _literal_origin_stage(_Matching(p, pair), inverse=True)
                self.assertTrue(np.array_equal(reverse, stage.T))
                self.assertTrue(np.array_equal(stage @ psi, psi))
                witness = np.concatenate((psi, np.zeros(dimension, dtype=np.int64)))
                self.assertTrue(np.array_equal(witness[swap], witness))
                for inverse in (False, True):
                    self.assertTrue(np.array_equal(
                        _asymmetric_branch(witness, frame, stage, swap, inverse),
                        witness * (1 << (2*exponent))))
                # Both arbitrary t0 inputs are part of the effective operator
                # domain. Exact component adjoints establish the full-domain
                # inverse; these probes additionally distinguish B/R order.
                for selector, t0 in product((0, 1, dimension-1), (0, 1)):
                    vector = np.zeros(2*dimension, dtype=np.int64)
                    vector[selector+t0*dimension] = 1
                    output = _asymmetric_branch(vector, frame, stage, swap)
                    returned = _asymmetric_branch(output, frame, stage, swap, inverse=True)
                    self.assertTrue(np.array_equal(returned, vector*(1 << (4*exponent))))

    def test_inside_frame_bank_swap_breaks_old_witness_but_not_new_code(self):
        p, dimension = 3, 256
        old = np.zeros(512, dtype=np.int64)
        old[[68, 70]] = 1
        new = np.zeros(512, dtype=np.int64)
        new[_code_columns(p)] = 1
        swap = _swap_images(p)
        broken = []
        for spec, pairs in zip(_FRAMES, (((0, 1),), ((1, 2),), ((2, 0),))):
            frame, exponent = _wrapped_frame(p, spec)
            stage = _literal_origin_stage(_Matching(p, pairs))
            broken.append(not np.array_equal(
                _asymmetric_branch(old, frame, stage, swap), old*(1 << (2*exponent))))
            self.assertTrue(np.array_equal(
                _asymmetric_branch(new, frame, stage, swap), new*(1 << (2*exponent))))
        self.assertTrue(any(broken))

    def test_bank_h_complete_marked_gram_has_a_rational_diagonal_dominance_floor(self):
        real, radical, denominator = _complete_bank_h_gram()
        expected_real = np.array([
            [131072,0,0,17920,0,16128,20096,0],
            [0,131072,35584,0,20352,0,0,20096],
            [0,35584,131072,0,26304,0,0,11136],
            [17920,0,0,131072,0,26304,10336,0],
            [0,20352,26304,0,131072,0,0,17920],
            [16128,0,0,26304,0,131072,35584,0],
            [20096,0,0,10336,0,35584,131072,0],
            [0,20096,11136,0,17920,0,0,131072]], dtype=np.int64)
        expected_radical = np.array([
            [0,2048,0,0,2048,0,0,2624],
            [2048,0,0,0,0,2048,1952,0],
            [0,0,0,2048,0,1184,2048,0],
            [0,0,2048,0,704,0,0,2048],
            [2048,0,0,704,0,2048,0,0],
            [0,2048,1184,0,2048,0,0,0],
            [0,1952,2048,0,0,0,0,2048],
            [2624,0,0,2048,0,0,2048,0]], dtype=np.int64)
        self.assertEqual(denominator, 131072)
        self.assertTrue(np.array_equal(real, expected_real))
        self.assertTrue(np.array_equal(radical, expected_radical))
        self.assertTrue(np.array_equal(real, real.T) and np.array_equal(radical, radical.T))
        self.assertTrue(np.all(real >= 0) and np.all(radical >= 0))
        self.assertTrue(np.all(np.diag(real) == denominator))
        self.assertFalse(np.any(np.diag(radical)))
        # sqrt(2)<3/2 gives exact rational upper bounds for every off-diagonal
        # row sum. Their maximum is <2/3, so symmetric diagonal dominance
        # proves G >= I/3. This is a Gram certificate, not a common-vector test.
        upper_twice = 2*(real.sum(axis=1)-np.diag(real))+3*radical.sum(axis=1)
        self.assertFalse(np.any(upper_twice % 2))
        upper = upper_twice // 2
        self.assertEqual(upper.tolist(), [64224,85104,80944,61760,71776,85936,75088,59232])
        self.assertTrue(np.all(3*upper < 2*denominator))
        # Empty-bank branches are I (both bank-controlled walls are inactive).
        # For arbitrary complex path coefficients, a=|sum(c)|^2<=8||c||^2.
        # Marked error is >=a/24 and empty error >=|a-1/2|, hence >=1/50.
        # The stated dyadic floor is weaker; coefficient realizability is not
        # assumed and t remains part of the complete arbitrary Gram domain.
        self.assertEqual(Fraction(1, 3)/(2*(8+Fraction(1, 3))), Fraction(1, 50))
        self.assertGreater(Fraction(1, 50), Fraction(1, 64))


if __name__ == '__main__':
    unittest.main()
