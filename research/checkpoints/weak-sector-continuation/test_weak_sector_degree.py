"""Scratch exact fixtures for the weak-probe sector seed-degree theorem.

No frozen repository files are modified. The shared tableau helper is a
bounded exact evaluator, not a native all-size compiler.
"""
from functools import reduce
from itertools import product
import unittest

from tests.exact_matrices import (CQ2, Q2, CO, CZ, CI, adjoint, eye,
                                 multiply, native_columns, scale, zeros)
from tests.field_sector_tableaux import (BitslicedTableau, active_derivative,
                                        helper_shift, query_word, sector_tableau)


def vector_degree(values):
    coefficients = list(values)
    bits = (len(values) - 1).bit_length()
    assert len(values) == 1 << bits
    for bit in range(bits):
        for mask in range(1 << bits):
            if mask >> bit & 1:
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return max((mask.bit_count() for mask, value in enumerate(coefficients)
                if value), default=-1)


def translate_mask(mask, address, size):
    return sum(((mask >> j) & 1) << (j ^ address) for j in range(size))


def label_action(columns, label):
    result = 0
    while label:
        bit = label & -label
        result ^= columns[bit.bit_length() - 1]
        label ^= bit
    return result


def tableau_columns(tableau, character):
    return [x | (z << tableau.n) for x, z, sign in tableau.at(character)]


def inverse_z_label(tableau, character, wire):
    # M^-1=J M^T J: the pullback of Z_wire is the corresponding X row,
    # transposed with its two column halves exchanged.
    n = tableau.n
    return (sum(((tableau.x[wire][n+i] >> character) & 1) << i
                for i in range(n))
            | (sum(((tableau.x[wire][i] >> character) & 1) << i
                   for i in range(n)) << n))


def permutation_matrix(width, permutation):
    result = zeros(1 << width)
    for column in range(1 << width):
        result[permutation(column)][column] = CO
    return result


def wire_translation_matrix(width, address):
    return permutation_matrix(width,
                              lambda word: translate_mask(word, address, width))


def z_matrix(width, wire):
    result = zeros(1 << width)
    for word in range(1 << width):
        result[word][word] = -CO if word >> wire & 1 else CO
    return result


def tensor(left, right):
    rows, columns = len(right), len(right[0])
    return [[left[i // rows][j // columns] * right[i % rows][j % columns]
             for j in range(len(left[0]) * columns)]
            for i in range(len(left) * rows)]


def block_diagonal(blocks):
    block_size = len(blocks[0])
    result = zeros(block_size * len(blocks))
    for index, block in enumerate(blocks):
        for i, row in enumerate(block):
            for j, value in enumerate(row):
                result[index * block_size+i][index * block_size+j] = value
    return result


class WeakSectorDegreeTests(unittest.TestCase):
    def test_split_translation_family_attains_every_degree_tradeoff(self):
        # Q translates the first d address coordinates. A helper-controlled
        # wire translation handles the other p-d coordinates. This saturates
        # both degree and helper-rank bounds and does not use a square-zero
        # query when d>=2.
        for p in range(1, 6):
            m = 1 << p
            for d in range(p+1):
                h = p-d
                seeds = [1 << (r << d) for r in range(1 << h)]
                self.assertEqual(vector_degree(seeds), h)
                query_entry = [int((x & ((1 << d)-1)) == 0)
                               for x in range(m)]
                self.assertEqual(vector_degree(query_entry), d)
                for r in range(1 << h):
                    total = 0
                    for x in range(m):
                        query_address = x & ((1 << d)-1)
                        shift = x >> d
                        moved_seed = translate_mask(seeds[r ^ shift],
                                                    query_address, m)
                        pullback = 1 << (x ^ (r << d))
                        self.assertEqual(moved_seed, pullback)
                        total ^= moved_seed
                    self.assertEqual(total, (1 << m)-1)
        # Two independent translation deviations have a nonzero product.
        first = translate_mask(1, 2, 4) ^ 1
        second = translate_mask(first, 1, 4) ^ first
        self.assertEqual(second, 15)

    def test_augmented_seed_cocycle_and_top_active_coefficient(self):
        for p in (2, 3, 4):
            table = sector_tableau(p)
            n = table.n
            seed_bits = table.x[0][n:] + table.x[0][:n]
            top = [active_derivative(table, value, p) for value in seed_bits]
            self.assertEqual(reduce(int.__or__, top, 0), table.ones)
            seeds = [inverse_z_label(table, r, 0) for r in range(table.sectors)]
            self.assertEqual(vector_degree(seeds), p-1)
            query_columns = [tableau_columns(
                BitslicedTableau(n, 0).apply(query_word(p, x)), 0)
                for x in range(1 << p)]
            top_address = 1 << (p-1)
            for r in range(table.sectors):
                total = 0
                for x in range(1 << p):
                    shifted = r ^ helper_shift(p, x)
                    actual = label_action(query_columns[x], seeds[shifted])
                    wanted = inverse_z_label(table, r, x)
                    self.assertEqual(actual, wanted)
                    total ^= actual
                self.assertNotEqual(total, 0)
                coefficient = reduce(int.__xor__,
                                     (seeds[r ^ helper_shift(p, a)]
                                      for a in range(top_address)), 0)
                derivative = (label_action(query_columns[top_address], coefficient)
                              ^ coefficient)
                self.assertEqual(derivative, total)

    def test_exact_nonclifford_sectors_with_projective_helper_query(self):
        # p=2, Q_x translates the low bit and one helper supplies the high bit.
        # E_r=T_0 W_(2r) is not Clifford, but every Z pullback is a Pauli.
        # The helper query is X^a Z^b in the physical helper basis, hence
        # Z^a X^b in its X basis. The projective sign is retained literally.
        n = 4
        diagonal = native_columns(n, [('T', 0)])
        blocks = [multiply(diagonal, wire_translation_matrix(n, 2*r))
                  for r in (0, 1)]
        decoder = block_diagonal(blocks)
        inverse = adjoint(decoder)
        seed = multiply(multiply(inverse, tensor(eye(2), z_matrix(n, 0))), decoder)
        for r in (0, 1):
            self.assertEqual(multiply(multiply(adjoint(blocks[r]), z_matrix(n, 0)),
                                      blocks[r]), z_matrix(n, 2*r))
        helper_queries = []
        for address in range(4):
            a, b = address & 1, address >> 1
            helper = zeros(2)
            for r in (0, 1):
                helper[r ^ b][r] = -CO if a and (r ^ b) else CO
            helper_queries.append(helper)
            source = tensor(helper, wire_translation_matrix(n, a))
            image = multiply(multiply(multiply(multiply(decoder, source), seed),
                                      adjoint(source)), inverse)
            self.assertEqual(image, tensor(eye(2), z_matrix(n, address)))
        self.assertEqual(multiply(helper_queries[1], helper_queries[2]),
                         scale(multiply(helper_queries[2], helper_queries[1]), -CO))
        x0 = permutation_matrix(n, lambda word: word ^ 1)
        moved_x = multiply(multiply(blocks[0], x0), adjoint(blocks[0]))
        pauli_entries = {CZ, CO, -CO, CI, -CI}
        self.assertTrue(any(value not in pauli_entries for row in moved_x for value in row))

    def test_nonclifford_pauli_seed_does_not_allow_uniform_pauli_rounding(self):
        # E rotates core wire1 around X and leaves the seed Z0 exactly Pauli.
        # With the wire-translation query, one output is a non-Pauli observable
        # at distance sqrt(2-sqrt(2)) from Z1, below the Clifford snap radius.
        n = 4
        decoder = native_columns(n, [('H', 1), ('T', 1), ('H', 1)])
        inverse = adjoint(decoder)
        self.assertEqual(multiply(multiply(inverse, z_matrix(n, 0)), decoder),
                         z_matrix(n, 0))
        image = multiply(multiply(decoder, z_matrix(n, 1)), inverse)
        wanted = z_matrix(n, 1)
        difference = [[a-b for a, b in zip(row, other)]
                      for row, other in zip(image, wanted)]
        self.assertNotEqual(image, wanted)
        self.assertEqual(multiply(difference, difference),
                         scale(eye(1 << n), CQ2(Q2(2, -1))))
        self.assertLess(Q2(2, -1), 2)
        self.assertTrue(all(sum(value != CZ for value in row) == 2 for row in image))

    def test_anticommutation_and_parity_mismatch_bounds(self):
        x = [[CZ, CO], [CO, CZ]]
        z = [[CO, CZ], [CZ, -CO]]
        total = [[a+b for a, b in zip(row, other)] for row, other in zip(x, z)]
        self.assertEqual(multiply(total, total), scale(eye(2), CQ2(2)))
        # A commuting Pauli relation has product sigma. Every target character
        # with opposite product is at squared Hamming distance at least4.
        for m in (2, 4, 8):
            for sigma in (-1, 1):
                minimum = None
                for prefix in product((-1, 1), repeat=m-1):
                    a = prefix + (sigma * reduce(int.__mul__, prefix, 1),)
                    for b in product((-1, 1), repeat=m):
                        if reduce(int.__mul__, b, 1) != -sigma:
                            continue
                        distance = sum((u-v)**2 for u, v in zip(a, b))
                        self.assertGreaterEqual(distance, 4)
                        minimum = distance if minimum is None else min(minimum, distance)
                        # Standard Q0=I makes A0=B0 exactly; the same mismatch
                        # can then be charged entirely to the other m-1 terms.
                        if a[0] == b[0]:
                            self.assertEqual(distance, sum((u-v)**2 for u, v in zip(a[1:], b[1:])))
                self.assertEqual(minimum, 4)

    def test_reversible_helper_tree_raises_degree_and_returns_helpers(self):
        # Two forward Toffoli layers compute into arbitrary, not clean,
        # internal helpers. One affine core control reads the root; two actual
        # inverse helper layers restore every input. The final core condition
        # has degree4 in the original helper coordinates.
        conditions = []
        for r in range(128):
            bits = [(r >> j) & 1 for j in range(7)]
            original = list(bits)
            bits[4] ^= bits[0] & bits[1]
            bits[5] ^= bits[2] & bits[3]
            bits[6] ^= bits[4] & bits[5]
            conditions.append(bits[6])
            bits[6] ^= bits[4] & bits[5]
            bits[5] ^= bits[2] & bits[3]
            bits[4] ^= bits[0] & bits[1]
            self.assertEqual(bits, original)
        self.assertEqual(vector_degree(conditions), 4)


if __name__ == '__main__':
    unittest.main()
