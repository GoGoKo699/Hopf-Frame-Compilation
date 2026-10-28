"""Exact fixtures for the Clifford-dual seed fingerprint screens.

The all-size proofs are in NATIVE_COPY_INTERACTION.md, Sections 7--8.
These native and signed-Pauli fixtures do not supply a copy or an emitter.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations
import unittest

from tests.exact_matrices import CO, CZ, Q2, adjoint, multiply, native_columns, scale
from tests.test_clifford_orbit_query_rank import _one_qubit_cliffords
from tests.test_nonclifford_seed_spread import _pauli, _trace, _conjugate
from tests.test_field_shear_intertwiner import _pencil_columns


def q_formula(width, active, label):
    answer = F(1)
    for bit in range(active):
        local = ((label >> bit) & 1) | (((label >> (width + bit)) & 1) << 1)
        answer *= (F(1), F(-1, 2), F(0), F(-1, 2))[local]
    return answer


def symp(width, left, right):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) % 2


def linear(columns, value):
    answer = 0
    for bit, image in enumerate(columns):
        if value >> bit & 1:
            answer ^= image
    return answer


def sp4():
    basis = [1, 2, 4, 8]
    gram = [[symp(2, a, b) for b in basis] for a in basis]
    for columns in permutations(range(1, 16), 4):
        if all(symp(2, columns[i], columns[j]) == gram[i][j]
               for i in range(4) for j in range(i)):
            yield columns


def field_image(p, address, x, z, phase):
    k = 1 << (p - 1)
    mask = (1 << k) - 1
    xu, zv = x & mask, (z >> k) & mask
    cols = _pencil_columns(p, address)
    x ^= linear(cols, xu) << k
    z ^= sum(((column & zv).bit_count() & 1) << i
             for i, column in enumerate(cols))
    return x, z, phase


def normal_product(factors):
    x = z = phase = 0
    for xx, zz, pp in factors:
        phase += pp + 2 * (z & xx).bit_count()
        x ^= xx
        z ^= zz
    return x, z, phase % 4


class SeedFingerprintTests(unittest.TestCase):
    def test_literal_native_fingerprint_with_arbitrary_helper(self):
        width = 2
        r = native_columns(width, [('TDG', 0), ('H', 0), ('T', 0)])
        self.assertEqual(adjoint(r), r)
        for label in range(1 << (2 * width)):
            p = _pauli(width, label)
            actual = _trace(multiply(multiply(multiply(r, p), r), p)) / 4
            self.assertEqual(actual, CO * Q2(q_formula(width, 1, label)))
        # Every helper input is retained in these complete matrices.
        self.assertEqual({label for label in range(16)
                          if q_formula(2, 1, label) == 1}, {0, 2, 8, 10})
        for active in range(1, 8):
            spectrum = {F(0)} | {F(-1, 2) ** k for k in range(active + 1)}
            self.assertGreaterEqual(min(abs(v - 1) for v in spectrum if v != 1),
                                    F(3, 4))
            self.assertGreaterEqual(min(abs(v + F(1, 2)) for v in spectrum
                                        if v != F(-1, 2)), F(3, 8))

    def test_exhaustive_two_wire_symplectic_fingerprint_normalizers(self):
        group = list(sp4())
        self.assertEqual(len(group), 720)
        for active, expected in ((1, 12), (2, 8)):
            q = [q_formula(2, active, label) for label in range(16)]
            special = {F(1), F(-1, 2)}
            preserving = [cols for cols in group if all(
                q[linear(cols, label)] == value
                for label, value in enumerate(q) if value in special)]
            self.assertEqual(len(preserving), expected)
            helper = {label for label in range(16) if q[label] == 1}
            active_space = {label for label in range(16)
                            if all(symp(2, label, h) == 0 for h in helper)}
            active_xy = {1 << i for i in range(active)} | {
                (1 << i) | (1 << (2 + i)) for i in range(active)}
            for columns in preserving:
                self.assertEqual({linear(columns, h) for h in helper}, helper)
                self.assertEqual({linear(columns, v) for v in active_space}, active_space)
                self.assertEqual({linear(columns, v) for v in active_xy}, active_xy)
                self.assertTrue(all(q[linear(columns, v)] == q[v] for v in range(16)))

    def test_all_literal_one_cell_cliffords_and_product_overlap_gap(self):
        a = native_columns(1, [('TDG', 0), ('H', 0), ('T', 0)])
        cliffords = list(_one_qubit_cliffords().values())
        self.assertEqual(len(cliffords), 24)
        normalizer_overlaps = []
        exact_signs = []
        for u in cliffords:
            image = _conjugate(u, a)
            fingerprint_ok = True
            for label in (1, 3):  # active X and Y
                p = _conjugate(u, _pauli(1, label))
                q = _trace(multiply(multiply(multiply(a, p), a), p)) / 2
                fingerprint_ok &= q == CO * Q2(F(-1, 2))
            if not fingerprint_ok:
                continue
            overlap = _trace(multiply(a, image)) / 2
            self.assertEqual(overlap.imag, 0)
            normalizer_overlaps.append(overlap.real)
            for sign in (-1, 1):
                if image == scale(a, sign):
                    exact_signs.append(sign)
                    # The unique maximal coefficient's Pauli has the same sign.
                    self.assertEqual(_conjugate(u, _pauli(1, 2)), scale(_pauli(1, 2), sign))
        self.assertEqual(len(normalizer_overlaps), 8)
        self.assertEqual(Counter(normalizer_overlaps), Counter({
            Q2(1): 1, Q2(-1): 1, Q2(F(1, 2)): 2,
            Q2(F(-1, 2)): 2, Q2(0): 2}))
        self.assertEqual(sorted(exact_signs), [-1, 1])
        for count in range(1, 6):
            values = {F(1)}
            for _ in range(count):
                values = {left * right for left in values
                          for right in (F(-1), F(-1, 2), F(0), F(1, 2), F(1))}
            self.assertTrue(all(abs(value) <= F(1, 2) for value in values if abs(value) != 1))

    def test_affine_field_four_pauli_word_keeps_literal_scalar_phase(self):
        totals = Counter()
        for p in (2, 3, 4):
            width = 1 << p
            mask = (1 << width) - 1
            labels = (range(1 << (2 * width)) if p == 2 else
                      [0, mask, mask << width, mask | (mask << width)] +
                      [1 << bit for bit in range(2 * width)] +
                      [(1 << bit) | (1 << ((bit + 3) % (2 * width)))
                       for bit in range(2 * width)])
            for label in labels:
                x, z = label & mask, label >> width
                phase = (x & z).bit_count() % 4
                for origin in range(1 << p):
                    addresses = [origin, origin ^ 1, origin ^ 2, origin ^ 3]
                    word = [field_image(p, address, x, z, phase) for address in addresses]
                    xx, zz, pp = normal_product(word)
                    self.assertEqual((xx, zz), (0, 0))
                    totals[pp] += 1
                    # Exactly one address equals origin: scalar is negated, never fixed.
                    requested_phase = (pp + 2 * sum(a == origin for a in addresses)) % 4
                    self.assertNotEqual(requested_phase, pp)
        self.assertGreater(sum(totals.values()), 2000)
        self.assertEqual(set(totals), {0})
        # A Clifford label family of affine degree can still have a -I scalar
        # cube; literal phases must not be silently replaced by +I.
        word = [(1, 0, 0), (0, 1, 0), (1, 0, 0), (0, 1, 0)]
        self.assertEqual(normal_product(word), (0, 0, 2))


if __name__ == '__main__':
    unittest.main(verbosity=2)
