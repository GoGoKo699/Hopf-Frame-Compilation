"""Exact marked-sector screens for specified entangling selector sources.

The common-module identities have all-size proofs in the research note.
The two-bit Gram certificate only screens the displayed four-branch source
at p=2; it does not exclude a family patched at small sizes, other accepted
words, or general native T-depth. All amplitudes below are exact fractions.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.test_co_moving_lookup_stage import _Matching


_ONE = Fraction(1)


def _accumulate(vector, key, value):
    vector[key] = vector.get(key, Fraction(0)) + value
    if not vector[key]:
        del vector[key]


def _hadamard(vector, rows):
    """Even-sized tensor Hadamard, retaining every coherent output."""
    rows = tuple(rows)
    if len(rows) % 2:
        raise ValueError("this rational fixture requires an even row count")
    mask = sum(1 << row for row in rows)
    denominator = 1 << (len(rows) // 2)
    outputs = tuple(sum(((word >> j) & 1) << row for j, row in enumerate(rows))
                    for word in range(1 << len(rows)))
    result = {}
    for selectors, amplitude in vector.items():
        for local in outputs:
            target = (selectors & ~mask) | local
            sign = (-1) ** ((selectors & local).bit_count())
            _accumulate(result, target, amplitude * Fraction(sign, denominator))
    return result


def _zeta(vector, p, reverse=False):
    result = {}
    for selectors, amplitude in vector.items():
        for coordinate in range(p):
            for row in range(1 << p):
                if (row >> coordinate) & 1:
                    continue
                other = row ^ (1 << coordinate)
                control, target = (other, row) if reverse else (row, other)
                selectors ^= ((selectors >> control) & 1) << target
        _accumulate(result, selectors, amplitude)
    return result


def _partial_fourier(vector, p, coordinate):
    """C_j=H_{u_j=0} M H_{u_j=0}; this is an involution."""
    rows = tuple(row for row in range(1 << p) if not ((row >> coordinate) & 1))
    return _hadamard(_zeta(_hadamard(vector, rows), p), rows)


def _fixed_origin_stage(vector, matching, bank=1, inverse=False):
    """Use the literal existing stage, at x=0 and Y=e0 or the empty bank."""
    if bank not in (0, 1):
        raise ValueError("this fixture only promises the two fixed bank inputs")
    result = {}
    for selectors, amplitude in vector.items():
        x, output_bank, output, borrowed, phase = matching.stage_literal(
            (0, bank, selectors, 1, 0), inverse=inverse)
        if (x, output_bank, borrowed, phase) != (0, bank, 1, 0):
            raise AssertionError("the fixed-origin permutation premise failed")
        _accumulate(result, output, amplitude)
    return result


def _graph_phase(selectors, p, coordinate):
    return sum(((selectors >> row) & 1) *
               ((selectors >> (row ^ (1 << coordinate))) & 1)
               for row in range(1 << p) if not ((row >> coordinate) & 1)) % 2


def _graph_state(p, coordinate):
    m = 1 << p
    return {s: Fraction((-1) ** _graph_phase(s, p, coordinate), 1 << (m // 2))
            for s in range(1 << m)}


def _four_branch(selectors, a, b, full_h, bank=1):
    """The corrected, moved-direction choice: C0 around R_C, C1 around R_D."""
    first, second = _Matching(2, ((0, 1),)), _Matching(2, ((1, 0),))
    vector = {selectors: _ONE}
    if a:
        vector = _partial_fourier(vector, 2, 0)
    vector = _fixed_origin_stage(vector, first, bank)
    if a:
        vector = _partial_fourier(vector, 2, 0)
    if b:
        vector = _partial_fourier(vector, 2, 1)
    if full_h:
        vector = _hadamard(vector, range(4))
    vector = _fixed_origin_stage(vector, second, bank)
    if full_h:
        vector = _hadamard(vector, range(4))
    if b:
        vector = _partial_fourier(vector, 2, 1)
    return vector


class EntanglingMarkedAcceptanceTests(unittest.TestCase):
    def test_common_zero_plus_module_and_actual_inverse(self):
        for p, families in ((2, (((0, 1),), ((1, 0),))),
                            (3, (((0, 1),), ((1, 2),), ((2, 0),)))):
            m = 1 << p
            zero = {0: _ONE}
            plus = {s: Fraction(1, 1 << (m // 2)) for s in range(1 << m)}
            for vector, reverse in product((zero, plus), (False, True)):
                self.assertEqual(_zeta(vector, p, reverse), vector)
            self.assertEqual(_hadamard(zero, range(m)), plus)
            self.assertEqual(_hadamard(plus, range(m)), zero)
            for pairs, vector, inverse, bank in product(families, (zero, plus),
                                                       (False, True), (0, 1)):
                self.assertEqual(_fixed_origin_stage(vector, _Matching(p, pairs), bank,
                                                    inverse=inverse), vector)

    def test_partial_fourier_graph_state_and_moved_direction_response(self):
        for p in (2, 3):
            matching = _Matching(p, ((0, 1),))
            for coordinate in (0, 1):
                graph = _graph_state(p, coordinate)
                self.assertEqual(_partial_fourier({0: _ONE}, p, coordinate), graph)
                self.assertEqual(_partial_fourier(graph, p, coordinate), {0: _ONE})
                # A tensor of two-vertex graph states is fixed by full H.
                self.assertEqual(_hadamard(graph, range(1 << p)), graph)
                actual = _fixed_origin_stage(graph, matching)
                expected = {}
                for selectors, amplitude in graph.items():
                    relabeled = matching.relabel(selectors)
                    contrast = (_graph_phase(selectors, p, coordinate)
                                + _graph_phase(relabeled, p, coordinate)) % 2
                    expected[selectors] = amplitude * (-1) ** ((selectors & 1) * contrast)
                self.assertEqual(actual, expected)
                if coordinate == 1:
                    # C e_j=e_j: a target direction is not a genuine escape.
                    self.assertEqual(actual, graph)
                else:
                    # The control direction moves, and this exact state changes.
                    self.assertNotEqual(actual, graph)

    def test_corrected_two_bit_marked_gram_and_empty_bank_screen(self):
        half = Fraction(1, 2)
        expected_grams = (
            [[16, 9, 9, 7], [9, 16, 6, 9], [9, 6, 16, 9], [7, 9, 9, 16]],
            [[16, 9, 10, 13*half], [9, 16, 6, 10],
             [10, 6, 16, 9], [13*half, 10, 9, 16]],
        )
        for full_h, expected in enumerate(expected_grams):
            gram = [[Fraction(0) for _ in range(4)] for _ in range(4)]
            for selectors in range(16):
                outputs = [_four_branch(selectors, a, b, full_h)
                           for b, a in product(range(2), repeat=2)]
                for output in outputs:
                    self.assertEqual(sum(v*v for v in output.values()), 1)
                for i, j in product(range(4), repeat=2):
                    gram[i][j] += sum(value * outputs[j].get(key, 0)
                                      for key, value in outputs[i].items())
                for b, a in product(range(2), repeat=2):
                    self.assertEqual(_four_branch(selectors, a, b, full_h, bank=0),
                                     {selectors: _ONE})
            self.assertEqual(gram, expected)
            # G-8*11^T is Hermitian, strictly diagonally dominant and positive.
            # Hence average marked acceptance >= |sum c|^2/2 for complex c.
            for i in range(4):
                self.assertEqual(gram[i][i] - 8, 8)
                self.assertLessEqual(sum(abs(gram[i][j]-8) for j in range(4) if j != i), 5)
            # Empty-bank energy is |sum c|^2, with target1/2.  If the complete
            # Gram error is epsilon, epsilon >= (1/2)(1/2-epsilon), hence1/6.
            self.assertEqual(Fraction(1, 4) / (1 + half), Fraction(1, 6))


if __name__ == '__main__':
    unittest.main()
