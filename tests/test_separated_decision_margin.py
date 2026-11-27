"""Exact operator audits for the two-address-block margin screen.

Proof home: research/depth/PRESCRIBED_SCALAR_TABLES.md,
#two-address-block-margin. The bounded fixtures check full compressed
operators, clean-code trace normalization, and a literal native escape
with three address blocks (four when dirty work must return). They do
not infer an asymptotic or unrestricted native-depth lower bound.
"""

from collections import Counter
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns, scale,
    zeros,
)
from tests.test_operator_source_compiler import _adjoint


INV_SQRT_TWO = CQ2(Q2(0, Fraction(1, 2)))


def _majoranas(width):
    """Jordan--Wigner X/Y modes, with least-significant work wire first."""
    size = 1 << width
    result = []
    for wire in range(width):
        for use_y in (False, True):
            matrix = zeros(size)
            for basis in range(size):
                sign = -1 if (basis & ((1 << wire) - 1)).bit_count() % 2 else 1
                phase = CI * (-1 if (basis >> wire) & 1 else 1) if use_y else CO
                matrix[basis ^ (1 << wire)][basis] = sign * phase
            result.append(matrix)
    return result


def _signal_rotation(core, opposite=False):
    """The complete block unitary [[I,-R],[R,I]]/sqrt(2), or its adjoint."""
    size = len(core)
    identity = eye(size)
    sign = -1 if opposite else 1
    blocks = [
        [identity[row] + [-sign * value for value in core[row]]
         for row in range(size)],
        [[sign * value for value in core[row]] + identity[row]
         for row in range(size)],
    ]
    return scale(blocks[0] + blocks[1], INV_SQRT_TWO)


def _private_ccz_batch(signal, address, bank, helpers):
    """Eight legal T bands; Gray phases have only private dirty targets."""
    bands = []
    word = []
    for step in range(8):
        band = [('T' if step % 2 == 0 else 'TDG', helper)
                for helper in helpers]
        bands.append(band)
        word += band
        for u, w, helper in zip(address, bank, helpers):
            # 000,100,110,010,011,111,101,001,000 on (u,w,c).
            control = (u, w, u, signal, u, w, u, signal)[step]
            word.append(('CX', control, helper))
    return word, bands


def _walsh_echo(half_width, complete):
    """Chronological native word; matrix products in the proof act rightmost first."""
    signal = 0
    u = list(range(1, 1 + half_width))
    v = list(range(1 + half_width, 1 + 2 * half_width))
    bank = list(range(1 + 2 * half_width, 1 + 3 * half_width))
    helpers = list(range(1 + 3 * half_width, 1 + 4 * half_width))
    phase, bands = _private_ccz_batch(signal, u, bank, helpers)
    shift = [('CX', control, target) for control, target in zip(v, bank)]
    # A_u=H_c P_u H_c and B_v=shift. The middle H_c pair commutes
    # through B_v and cancels, so this chronological word is A_u B_v A_u.
    raw = [('H', signal)] + phase + shift + phase + [('H', signal)]
    # Prepending chronological B_v gives the matrix A_u B_v A_u B_v.
    return (shift + raw if complete else raw), bands + bands, helpers


def _trace(matrix):
    return sum((row[index] for index, row in enumerate(matrix)), CZ)


class SeparatedDecisionMarginTests(unittest.TestCase):
    def test_majorana_construction_saturates_full_compressed_margin(self):
        for work_width in (1, 2):
            modes = _majoranas(work_width)
            m, work_size = len(modes), 1 << work_width
            identity = eye(work_size)
            inverse_root_m = INV_SQRT_TWO if m == 2 else CQ2(Fraction(1, 2))
            signal = zeros(2 * work_size)
            for basis in range(2 * work_size):
                signal[basis][basis] = CO if basis < work_size else -CO

            # Check the actual finite carrier algebra, including the Y phases.
            for left, gamma in enumerate(modes):
                self.assertEqual(adjoint(gamma), gamma)
                for right, other in enumerate(modes):
                    self.assertEqual(
                        add(multiply(gamma, other), multiply(other, gamma)),
                        scale(identity, 2 if left == right else 0))

            row_rotations = []
            for u in range(m):
                core = zeros(work_size)
                for v, gamma in enumerate(modes):
                    sign = -1 if (u & v).bit_count() % 2 else 1
                    core = add(core, scale(gamma, sign * inverse_root_m))
                rotation = _signal_rotation(core)
                self.assertEqual(multiply(adjoint(rotation), rotation), eye(2 * work_size))
                row_rotations.append(rotation)

            for v, gamma in enumerate(modes):
                column_rotation = _signal_rotation(gamma, opposite=True)
                self.assertEqual(multiply(adjoint(column_rotation), column_rotation),
                                 eye(2 * work_size))
                for u, row_rotation in enumerate(row_rotations):
                    with self.subTest(modes=m, row=u, column=v):
                        source = multiply(column_rotation, row_rotation)
                        observable = multiply(adjoint(source), multiply(signal, source))
                        compressed = [row[:work_size] for row in observable[:work_size]]
                        sign = -1 if (u & v).bit_count() % 2 else 1
                        self.assertEqual(compressed, scale(identity, sign * inverse_root_m))

    def test_native_three_block_decision_and_four_block_return_on_all_columns(self):
        for half_width in (1, 2):
            width = 1 + 4 * half_width
            size = 1 << width
            mask = (1 << half_width) - 1
            for complete in (False, True):
                with self.subTest(half_width=half_width, complete=complete):
                    word, bands, helpers = _walsh_echo(half_width, complete)
                    counts = Counter(gate[0] for gate in word)
                    self.assertEqual(counts['T'] + counts['TDG'], 16 * half_width)
                    self.assertEqual(counts['CX'], (18 if complete else 17) * half_width)
                    self.assertEqual(counts['H'], 2)
                    self.assertEqual(len(bands), 16)
                    self.assertEqual(sum(map(len, bands)), counts['T'] + counts['TDG'])
                    for band in bands:
                        self.assertEqual(len({gate[1] for gate in band}), half_width)
                        self.assertEqual({gate[1] for gate in band}, set(helpers))
                    self.assertTrue(all(gate[1] in helpers for gate in word
                                        if gate[0] in ('T', 'TDG')))
                    self.assertTrue(all(gate[1] == 0 for gate in word if gate[0] == 'H'))

                    expected = zeros(size)
                    for basis in range(size):
                        u = (basis >> 1) & mask
                        v = (basis >> (1 + half_width)) & mask
                        output = basis ^ ((u & v).bit_count() % 2)
                        if not complete:
                            output ^= v << (1 + 2 * half_width)
                        expected[output][basis] = CO
                    # All address, signal, dirty-bank and phase-helper columns
                    # are evaluated. The literal inverse must have the same
                    # involutory permutation, without a residual scalar phase.
                    self.assertEqual(native_columns(width, word), expected)
                    self.assertEqual(native_columns(width, _adjoint(word)), expected)

    def test_clean_code_trace_centering_with_noncommuting_native_boundaries(self):
        width, size = 3, 8
        preparation = [('H', 0), ('T', 0), ('CX', 1, 0)]
        a = native_columns(width, preparation)
        for readout in ('Z', 'X'):
            readout_word = [('S', 0), ('S', 0)]
            finish = [('TDG', 0), ('H', 0)]
            if readout == 'X':
                readout_word = [('H', 0)] + readout_word + [('H', 0)]
                finish = [('TDG', 0)]
            observable = native_columns(width, readout_word)
            b = native_columns(width, finish)
            source = native_columns(width, preparation + finish)
            pulled_back = multiply(adjoint(source), multiply(observable, source))
            output_factor = multiply(adjoint(b), multiply(observable, b))
            self.assertEqual(_trace(output_factor), CZ)

            for clean in (1, 2):
                with self.subTest(readout=readout, clean=clean):
                    code = list(range(0, size, 1 << clean))
                    compressed = [[pulled_back[row][column] for column in code]
                                  for row in code]
                    expected = zeros(len(code))
                    for index, basis in enumerate(code):
                        expected[index][index] = CZ if basis & 2 else CO
                    self.assertEqual(compressed, expected)
                    # One clean bit leaves the second wire arbitrary and gives
                    # bias 1/2; initializing it too changes the bias to one.
                    actual_bias = _trace(compressed) / len(code)
                    self.assertEqual(actual_bias, CQ2(Fraction(1, 2) if clean == 1 else 1))

                    centered = zeros(size)
                    for basis in range(size):
                        centered[basis][basis] = CQ2((1 << clean) - 1 if basis in code else -1)
                    input_factor = multiply(a, multiply(centered, adjoint(a)))
                    self.assertEqual(_trace(multiply(input_factor, input_factor)) / size,
                                     CQ2((1 << clean) - 1))
                    self.assertEqual(_trace(multiply(input_factor, output_factor)) / size,
                                     actual_bias)


if __name__ == '__main__':
    unittest.main()
