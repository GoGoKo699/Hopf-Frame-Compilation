"""Exact Bell-cell witness for the specified bank-coupled selector frames.

The arbitrary second selector bank t belongs to the complete Gram domain;
it need not return after a source block.  Choosing t0=0 below supplies a
legal counterexample input, never a promised initialization.  Implementation
helpers of each literal stage do return.  The crossing-code witness survives
the asymmetric branch V^dagger B_Y R_C V, despite failure of an older witness.
A separate full-column Gram certificate rejects the specified bank-H extension.
The full-bank fixtures check cell bounds, complete small tensor-swap paths,
and cancellation of the redundant signal angle. The all-size trace bound
and its centered-recovery consequence are analytic, not a numerical search.
The single-hole fixture checks the signed zeta restriction and all reduced
path traces, including the exact wrong-sign mean of the native polar word.
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


def _literal_full_bank_origin_stage(matching, inverse=False):
    """All selector columns, including the actual quarter phases and inverse."""
    dimension = 1 << matching.m
    matrix = np.zeros((dimension, dimension), dtype=complex)
    for selectors, borrowed in product(range(dimension), (0, 1)):
        state = (0, matching.mask, selectors, borrowed, 0)
        output = matching.stage_literal(state, inverse=inverse)
        assert output[:2] == state[:2] and output[3] == borrowed
        assert matching.stage_literal(output, inverse=not inverse) == state
        matrix[output[2], selectors] = 1j ** output[4]
    return matrix


def _full_bank_small_components():
    """One-cell calibration of the all-size tensor identity, on all s,t inputs.

    These three p=2 stages are a bounded algebra check, not the p>=3
    all-coordinate matching family. Native frame phases are checked in the
    retained cell fixtures; the stages use their literal complete action.
    """
    dimension = 16
    frames, stages, reverses = [], [], []
    for spec, pair in zip(_FRAMES, (((0, 1),), ((1, 0),), ((0, 1),))):
        frame, exponent = _wrapped_frame(2, spec)
        frames.append(frame / (1 << exponent))
        matching = _Matching(2, pair)
        stages.append(_literal_full_bank_origin_stage(matching))
        reverses.append(_literal_full_bank_origin_stage(matching, inverse=True))
    hadamard = np.array([[(-1) ** ((row & col).bit_count()) / 4
                         for col in range(dimension)] for row in range(dimension)])
    swap = np.array([value // dimension + dimension * (value % dimension)
                     for value in range(dimension ** 2)])

    def on_s(operator, columns):
        return np.einsum('ab,tbc->tac', operator,
                         columns.reshape(dimension, dimension, -1)).reshape(
                             dimension ** 2, -1)

    def branch(j, enabled, columns, inverse=False):
        if not enabled:
            return on_s(reverses[j] if inverse else stages[j], columns)
        columns = on_s(frames[j], columns)
        if inverse:
            columns = on_s(hadamard, columns)[swap]
            columns = on_s(reverses[j], columns)
        else:
            columns = on_s(stages[j], columns)[swap]
            columns = on_s(hadamard, columns)
        return on_s(frames[j].conj().T, columns)

    return dimension, frames, stages, reverses, hadamard, swap, branch


class BankCoupledSelectorFrameTests(unittest.TestCase):
    def test_full_bank_cell_projectors_and_exact_tensor_statistics(self):
        cells = [_cell_frame(*spec) for spec in _FRAMES]
        identity = np.eye(16, dtype=np.int64)

        def pauli(wire, z=False):
            result = np.zeros((16, 16), dtype=np.int64)
            for column in range(16):
                row = column if z else column ^ (1 << wire)
                result[row, column] = (-1) ** ((column >> wire) & 1) if z else 1
            return result

        x0, x1, x3 = [pauli(wire) for wire in (0, 1, 3)]
        z2, z3 = [pauli(wire, z=True) for wire in (2, 3)]
        formulas = [
            4 * identity - (identity+x0) @ (identity-z2) @ (identity-z3),
            4 * identity - (identity+x0) @ (identity-x1) @ (identity-x3),
            4 * identity - (identity+x0) @ (identity-x1) @ (identity+z2),
        ]
        for (matrix, exponent), numerator in zip(cells, formulas):
            self.assertTrue(np.array_equal(4 * matrix, (1 << exponent) * numerator))

        def word(indices):
            matrix, exponent = identity.copy(), 0
            for index in indices:
                matrix = matrix @ cells[index][0]
                exponent += cells[index][1]
                matrix, exponent = _remove_twos(matrix, exponent)
            return matrix, 1 << exponent

        for i, j in product(range(3), repeat=2):
            if i == j:
                continue
            matrix, denominator = word((i, j))
            absolute = np.abs(matrix)
            for axis in (0, 1):
                self.assertEqual(Fraction(int(absolute.max(axis=axis).sum()),
                                          16 * denominator), Fraction(3, 4))
                expected = Fraction(2) if {i, j} == {0, 2} else Fraction(5, 2)
                self.assertEqual(Fraction(int(absolute.sum(axis=axis).max()),
                                          denominator), expected)

        j212, denominator = word((1, 0, 1))
        formula = 4 * identity - (identity+x0) @ (identity-z2) @ (identity-x1 @ z3)
        self.assertTrue(np.array_equal(4 * j212, denominator * formula))
        for matrix, scale in ((j212, denominator),
                              (cells[2][0], 1 << cells[2][1])):
            for axis in (0, 1):
                sums = np.abs(matrix).sum(axis=axis)
                self.assertEqual(sorted(sums.tolist()), [scale] * 8 + [2 * scale] * 8)
                self.assertEqual(Fraction(int(sums @ sums), 16 * scale ** 2),
                                 Fraction(5, 2))
        matrix, denominator = word((1, 0, 1, 2))
        for axis in (0, 1):
            self.assertTrue(np.all(2 * np.abs(matrix).max(axis=axis) == denominator))
            self.assertTrue(np.all(np.abs(matrix).sum(axis=axis) == 2 * denominator))

    def test_full_bank_tensor_paths_actual_reverses_and_signal_cancellation(self):
        dimension, frames, stages, reverses, hadamard, swap, branch = (
            _full_bank_small_components())
        size = dimension ** 2
        identity = np.eye(size, dtype=complex)
        for stage, reverse in zip(stages, reverses):
            self.assertTrue(np.array_equal(reverse, stage.conj().T))
        paths, factors = {}, {}
        for bits in product((0, 1), repeat=3):
            actual = identity.copy()
            left = right = np.eye(dimension, dtype=complex)
            for j, enabled in enumerate(bits):
                actual = branch(j, enabled, actual)
                if enabled:
                    left, right = (frames[j] @ hadamard @ right,
                                   stages[j] @ frames[j] @ left)
                else:
                    left = stages[j] @ left
            expected = np.kron(right, left)
            if sum(bits) % 2:
                expected = expected[:, swap]
            # Entries here are Gaussian dyadics with exactly representable
            # numerators; no floating-point tolerance hides a path phase.
            self.assertTrue(np.array_equal(actual, expected))
            returned = actual
            for j in range(2, -1, -1):
                returned = branch(j, bits[j], returned, inverse=True)
            self.assertTrue(np.array_equal(returned, identity))
            paths[bits], factors[bits] = actual, (left, right)

        for left_bits, right_bits in product(paths, repeat=2):
            a = factors[left_bits][0].conj().T @ factors[right_bits][0]
            b = factors[left_bits][1].conj().T @ factors[right_bits][1]
            actual = np.vdot(paths[left_bits], paths[right_bits]) / size
            if sum(left_bits) % 2 == sum(right_bits) % 2:
                expected = np.trace(a) * np.trace(b) / size
                if left_bits != right_bits:
                    self.assertLessEqual(abs(actual), Fraction(3, 4))
            else:
                expected = np.trace(a @ b) / size
                self.assertLessEqual(abs(actual), Fraction(1, dimension))
            self.assertEqual(actual, expected)

        signal_h = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)

        def rz(angle):
            return np.diag([np.exp(-0.5j * angle), np.exp(0.5j * angle)])

        def source(alpha, beta):
            signal = [signal_h, rz(alpha) @ signal_h, rz(beta) @ signal_h]
            signal.append((signal[2] @ signal[1] @ signal[0]).conj().T)
            columns = np.eye(2 * size, dtype=complex).reshape(2, size, 2 * size)
            columns = np.einsum('ab,bcd->acd', signal[0], columns)
            for j in range(3):
                columns = np.stack([branch(j, enabled, columns[enabled])
                                    for enabled in (0, 1)])
                columns = np.einsum('ab,bcd->acd', signal[j+1], columns)
            return columns.reshape(2 * size, 2 * size)

        # Fixed phase choices test the complete signal identity; this is not
        # an angle optimization or a claim of successful decision bias.
        for alpha in (0, np.pi / 2):
            actual = source(alpha, np.pi / 4)
            np.testing.assert_allclose(actual, source(alpha, 0), atol=2e-13, rtol=0)
            coefficients = {bits: np.exp(1j * alpha * bits[1])
                             * (-1) ** (bits[1] * (bits[0] + bits[2])) / 4
                             for bits in paths}
            expected = sum(coefficients[bits] * paths[bits] for bits in paths)
            np.testing.assert_allclose(actual[:size, :size], expected,
                                       atol=2e-13, rtol=0)
            self.assertAlmostEqual(abs(sum(coefficients.values()) - 1), 0, places=13)
            self.assertAlmostEqual(sum(abs(c) ** 2 for c in coefficients.values()),
                                   0.5, places=13)

        # A fixed, exact native polar-angle change escapes only the balanced
        # path-norm argument. This checks its calibration, not a signed gap.
        forward = [('S', 0)] * 3 + [('H', 0), ('T', 0), ('H', 0), ('S', 0)]
        reverse = [('S', 0)] * 3 + [('H', 0), ('TDG', 0), ('H', 0), ('S', 0)]
        w = native_columns(1, forward)
        last = native_columns(1, reverse * 3)
        coefficients = [last[0][b3] * w[b3][b2] * w[b2][b1] * w[b1][0]
                        for b1, b2, b3 in product((0, 1), repeat=3)]
        self.assertEqual(sum(coefficients, CZ), CO)
        self.assertEqual(sum((c.conj() * c for c in coefficients), CZ),
                         CQ2(Q2(Fraction(3, 8))))
        self.assertEqual(native_columns(1, forward * 3 + reverse * 3),
                         [[CO, CZ], [CZ, CO]])

    def test_single_hole_signed_restriction_and_exact_polar_mean(self):
        # x=0, Y=1-e0 and X_s0=-1 form an invariant unmarked sector.
        # t0 remains an arbitrary spectator; all other s,t inputs enter
        # the normalized trace through the complete tensor-swap factors.
        p, dimension, bank = 3, 128, 254
        signed_zeta = np.zeros((dimension, dimension), dtype=np.int64)
        zeta = _zeta_word(p, transpose=True)
        for tail in range(dimension):
            image = _binary_image(zeta, tail << 1)
            self.assertEqual(image & 1, tail.bit_count() % 2)
            signed_zeta[image >> 1, tail] = (-1) ** tail.bit_count()
        self.assertTrue(np.array_equal(
            signed_zeta @ signed_zeta, np.eye(dimension, dtype=np.int64)))

        frames, stages, exponents = [], [], []
        pairs = (((0, 1),), ((1, 2),), ((2, 0),))
        for spec, pair in zip(_FRAMES, pairs):
            physical, exponent = _wrapped_frame(p, spec)
            self.assertTrue(np.array_equal(physical[0::2, 0::2],
                                           physical[1::2, 1::2]))
            self.assertTrue(np.array_equal(physical[0::2, 1::2],
                                           physical[1::2, 0::2]))
            reduced = physical[0::2, 0::2] - physical[0::2, 1::2]
            cell, cell_exponent = _cell_frame(*spec)
            expected = signed_zeta @ np.kron(
                cell, np.eye(8, dtype=np.int64)) @ signed_zeta.T
            self.assertTrue(np.array_equal(
                reduced, expected * (1 << (exponent - cell_exponent))))
            frames.append(expected)
            exponents.append(cell_exponent)

            matching = _Matching(p, pair)
            stage = np.zeros((dimension, dimension), dtype=complex)
            for tail, s0, borrowed in product(range(dimension), (0, 1), (0, 1)):
                state = (0, bank, (tail << 1) | s0, borrowed, 0)
                output = matching.stage_literal(state)
                self.assertEqual(output[:2], (0, bank))
                self.assertEqual(output[2] & 1, s0)
                self.assertEqual(output[3], borrowed)
                self.assertEqual(matching.stage_literal(output, inverse=True), state)
                entry = 1j ** output[4]
                previous = stage[output[2] >> 1, tail]
                self.assertTrue(previous == 0 or previous == entry)
                stage[output[2] >> 1, tail] = entry
            self.assertTrue(np.array_equal(stage.conj().T @ stage,
                                           np.eye(dimension)))
            stages.append(stage)

        # H^(m-1) has numerator H_num, scalar sqrt(2)/16. Keep every
        # radical and denominator outside the Gaussian-integer matrices.
        hadamard = np.array([[(-1) ** ((row & col).bit_count())
                              for col in range(dimension)]
                             for row in range(dimension)], dtype=np.int64)
        paths = {}
        for bits in product((0, 1), repeat=3):
            left = right = np.eye(dimension, dtype=complex)
            for j, enabled in enumerate(bits):
                if enabled:
                    left, right = (frames[j] @ hadamard @ right,
                                   stages[j] @ frames[j] @ left)
                else:
                    left = stages[j] @ left
            count = sum(bits)
            exponent = (2 * sum(e * b for e, b in zip(exponents, bits))
                        + 4 * count - count // 2)
            paths[bits] = left, right, exponent, count % 2

        # Unitarity bounds all absolute dot-product accumulations below
        # 2^50, so these Gaussian-integer complex128 operations are exact.
        maximum = max(item[2] for item in paths.values())
        self.assertLessEqual(4 * dimension ** 2 * (1 << (2 * maximum)), 1 << 50)
        gram = []
        for left, right, exponent, parity in paths.values():
            row = []
            for other_left, other_right, other_exponent, other_parity in paths.values():
                a = left.conj().T @ other_left
                b = right.conj().T @ other_right
                numerator = (np.trace(a) * np.trace(b) if parity == other_parity
                             else np.trace(a @ b))
                self.assertEqual(numerator.real, round(numerator.real))
                self.assertEqual(numerator.imag, round(numerator.imag))
                factor = 2 if parity and other_parity else 1
                denominator = dimension ** 2 * (1 << (exponent + other_exponent))
                real = Fraction(factor * int(numerator.real), denominator)
                imag = Fraction(factor * int(numerator.imag), denominator)
                row.append(CQ2(Q2(0, real), Q2(0, imag))
                           if parity != other_parity else CQ2(real, imag))
            gram.append(row)
        for i, j in product(range(8), repeat=2):
            self.assertEqual(gram[i][j], gram[j][i].conj())
            if i == j:
                self.assertEqual(gram[i][j], CO)

        forward = [('S', 0)] * 3 + [('H', 0), ('T', 0), ('H', 0), ('S', 0)]
        reverse = [('S', 0)] * 3 + [('H', 0), ('TDG', 0), ('H', 0), ('S', 0)]
        wall, last = native_columns(1, forward), native_columns(1, reverse * 3)
        coefficients = [last[0][b3] * wall[b3][b2] * wall[b2][b1] * wall[b1][0]
                        for b1, b2, b3 in paths]
        mean = 2 * sum((coefficients[i].conj() * gram[i][j] * coefficients[j]
                        for i, j in product(range(8), repeat=2)), CZ) - CO
        self.assertEqual(mean, CQ2(Q2(Fraction(-260847, 1048576),
                                      Fraction(-5139, 4194304))))
        self.assertLess(mean.real, 0)
        # This exact p=3 mean is a finite corroboration. The all-size
        # wrong-sign mean and arbitrary-wall query bounds are analytic.

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
