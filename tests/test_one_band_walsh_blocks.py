"""Exact affine-support and native sign checks for one-band Walsh blocks.

Finite enumeration checks the linear-algebra ingredients, not an all-size
copy compiler or an unrestricted native-depth lower bound.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, CQ2, Q2, adjoint, multiply, native_columns,
)


def _span(generators):
    result = {0}
    for generator in generators:
        result |= {value ^ generator for value in tuple(result)}
    return frozenset(result)


def _orthogonal(space, dimension):
    return frozenset(value for value in range(1 << dimension)
                     if all((value & vector).bit_count() % 2 == 0 for vector in space))


def _subspaces(dimension):
    spaces = {frozenset({0})}
    pending = [frozenset({0})]
    while pending:
        space = pending.pop()
        for vector in range(1 << dimension):
            if vector not in space:
                extended = space | frozenset(vector ^ item for item in space)
                if extended not in spaces:
                    spaces.add(extended)
                    pending.append(extended)
    return sorted(spaces, key=lambda space: (len(space), tuple(sorted(space))))


def _affine_spaces(dimension):
    result = []
    for linear in _subspaces(dimension):
        seen = set()
        for offset in range(1 << dimension):
            affine = frozenset(offset ^ value for value in linear)
            if affine not in seen:
                seen.add(affine)
                result.append((offset, linear, tuple(sorted(affine))))
    return result


def _rank(matrix):
    if not matrix:
        return 0
    work = [[Fraction(value) for value in row] for row in matrix]
    pivot = 0
    for column in range(len(work[0])):
        found = next((row for row in range(pivot, len(work)) if work[row][column]), None)
        if found is None:
            continue
        work[pivot], work[found] = work[found], work[pivot]
        value = work[pivot][column]
        work[pivot] = [entry / value for entry in work[pivot]]
        for row in range(pivot + 1, len(work)):
            value = work[row][column]
            if value:
                work[row] = [entry - value * lead
                             for entry, lead in zip(work[row], work[pivot])]
        pivot += 1
        if pivot == len(work):
            break
    return pivot


def _walsh(rows, columns):
    return [[-1 if (row & column).bit_count() & 1 else 1 for column in columns]
            for row in rows]


def _matmul(left, right):
    return [[sum(a * b for a, b in zip(row, column)) for column in zip(*right)]
            for row in left]


def _pauli(width, x, z):
    size = 1 << width
    phase = (CO, CI, -CO, -CI)[(x & z).bit_count() % 4]
    return [[phase * (-1 if (z & column).bit_count() & 1 else 1)
             if row == column ^ x else CZ for column in range(size)]
            for row in range(size)]


def _coefficient(pauli, operator):
    size = len(pauli)
    return sum((pauli[row][column].conj() * operator[row][column]
                for row in range(size) for column in range(size)), CZ) / size


class OneBandWalshBlockTests(unittest.TestCase):
    def _check_affine_pair(self, dimension, source, target):
        source_offset, u, s = source
        target_offset, v, t = target
        v_perp = _orthogonal(v, dimension)
        intersection = u & v_perp
        h = len(intersection).bit_length() - 1
        dim_u = len(u).bit_length() - 1
        dim_v = len(v).bit_length() - 1
        compatible = v_perp <= u
        # Kernel of the rows outside T is exactly the simultaneous
        # time/Fourier support intersection; all entries here are integers.
        outside = sorted(set(range(1 << dimension)) - set(t))
        kernel_dimension = len(s) - _rank(_walsh(outside, s))
        wanted = (1 << (dim_u + dim_v - dimension)) if compatible else 0
        self.assertEqual(kernel_dimension, wanted)

        block = _walsh(t, s)
        gram = _matmul([list(column) for column in zip(*block)], block)
        multiplier = 1 << (dim_v + h)
        self.assertEqual(_matmul(gram, gram),
                         [[multiplier * entry for entry in row] for row in gram])
        rank = 1 << (dim_u - h)
        self.assertEqual(_rank(block), rank)
        self.assertEqual(sum(gram[j][j] for j in range(len(s))), multiplier * rank)
        eigenvalue = Fraction(multiplier, 1 << dimension)
        self.assertEqual(eigenvalue == 1, compatible)
        if not compatible:
            self.assertLessEqual(eigenvalue, Fraction(1, 2))

        # Character-offset signs in the Gram matrix are literal, even
        # when neither affine subspace passes through zero.
        for i, a in enumerate(s):
            for j, b in enumerate(s):
                expected = (len(v) * (-1 if (target_offset & (a ^ b)).bit_count() & 1 else 1)
                            if (a ^ b) in v_perp else 0)
                self.assertEqual(gram[i][j], expected)
        if compatible:
            # A disjoint-coset basis gives all admissible functions, with
            # f(z+h)=(-1)^(target_offset.h) f(z) on V^perp.
            unseen = set(s)
            witnesses = []
            while unseen:
                representative = min(unseen)
                coset = {representative ^ value for value in v_perp}
                unseen -= coset
                values = {z: (-1 if (target_offset & (z ^ representative)).bit_count() & 1 else 1)
                          for z in coset}
                witnesses.append(values)
                for frequency in outside:
                    self.assertEqual(sum(value * (-1 if (frequency & z).bit_count() & 1 else 1)
                                         for z, value in values.items()), 0)
            self.assertEqual(len(witnesses), wanted)

    def test_every_affine_support_pair_through_dimension_three(self):
        counts = []
        for dimension in range(4):
            spaces = _affine_spaces(dimension)
            counts.append(len(spaces))
            for source in spaces:
                for target in spaces:
                    self._check_affine_pair(dimension, source, target)
        self.assertEqual(counts, [1, 3, 11, 51])

    def test_selected_four_dimensional_translated_and_oblique_spaces(self):
        linear_spaces = [_span(vectors) for vectors in
                         ((), (1,), (3,), (1, 2), (3, 5), (3, 12),
                          (1, 2, 4), (3, 5, 9), (1, 2, 4, 8))]
        for u, v in product(linear_spaces, repeat=2):
            for source_offset, target_offset in ((0, 0), (7, 11), (15, 6)):
                source = source_offset, u, tuple(sorted(source_offset ^ z for z in u))
                target = target_offset, v, tuple(sorted(target_offset ^ z for z in v))
                self._check_affine_pair(4, source, target)

    def test_literal_native_t_and_inverse_t_signed_walsh_blocks(self):
        # All eight orientations of three simultaneous physical T gates.
        # A coefficient block has one coordinate for each active input X.
        width = 3
        for negative_mask in range(8):
            word = [('TDG' if negative_mask >> wire & 1 else 'T', wire)
                    for wire in range(width)]
            unitary = native_columns(width, word)
            for x in range(8):
                active = [wire for wire in range(width) if x >> wire & 1]
                inactive = [wire for wire in range(width) if not x >> wire & 1]
                count = len(active)
                normalization = CO
                for _ in range(count):
                    normalization *= CQ2(Q2(0, Fraction(1, 2)))
                # Both inactive Z patterns are sampled, including a
                # nonzero spectator phase when an inactive position exists.
                spectators = (0, sum(1 << wire for wire in inactive))
                for spectator in set(spectators):
                    for input_bits in range(1 << count):
                        z_in = spectator | sum(((input_bits >> i) & 1) << wire
                                               for i, wire in enumerate(active))
                        incoming = _pauli(width, x, z_in)
                        image = multiply(multiply(unitary, incoming), adjoint(unitary))
                        reconstructed = [[CZ for _ in range(8)] for _ in range(8)]
                        for output_bits in range(1 << count):
                            z_out = spectator | sum(((output_bits >> i) & 1) << wire
                                                    for i, wire in enumerate(active))
                            exponent = (input_bits & output_bits).bit_count()
                            for i, wire in enumerate(active):
                                exponent += ((output_bits if negative_mask >> wire & 1 else input_bits) >> i) & 1
                            coefficient = normalization * (-1 if exponent & 1 else 1)
                            outgoing = _pauli(width, x, z_out)
                            self.assertEqual(_coefficient(outgoing, image), coefficient)
                            for row in range(8):
                                for column in range(8):
                                    reconstructed[row][column] += coefficient * outgoing[row][column]
                        self.assertEqual(image, reconstructed)


if __name__ == '__main__':
    unittest.main()
