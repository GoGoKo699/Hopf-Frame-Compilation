"""Exact receipts for the complete-probe affine-support count constraint.

The robust count proof is analytic. A one-band fixture explicitly prevents
interpreting affine support dimension itself as a growing native depth cost.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CO
from tests.test_native_indicator_sharing import _affine_rank, _propagate
from tests.test_three_layer_native_mechanisms import _native_t_bands


def _signs(p, dirty=0):
    m = 1 << p
    return [1 - 2 * ((basis >> (p + (basis & (m - 1)))) & 1)
            for basis in range(1 << (p + m + dirty))]


def _walsh(values):
    values = list(values)
    half = 1
    while half < len(values):
        for start in range(0, len(values), 2 * half):
            for offset in range(half):
                a, b = values[start + offset], values[start + offset + half]
                values[start + offset], values[start + offset + half] = a + b, a - b
        half *= 2
    return [Fraction(value, len(values)) for value in values]


def _symplectic(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) & 1


class CompleteProbeSupportCountTests(unittest.TestCase):
    def test_all_signed_translation_symmetries_and_dirty_directions(self):
        p, m, dirty = 2, 4, 1
        width, signs = p + m + dirty, _signs(p, dirty)
        actual = {}
        for translation in range(1 << width):
            for scalar in (-1, 1):
                distance = max(abs(signs[basis ^ translation] - scalar * signs[basis])
                               for basis in range(1 << width))
                address = translation & (m - 1)
                bank = (translation >> p) & ((1 << m) - 1)
                expected = address == 0 and bank in (0, (1 << m) - 1)
                expected = expected and scalar == (1 if bank == 0 else -1)
                self.assertEqual(distance, 0 if expected else 2)
                if expected:
                    actual[translation] = scalar
        self.assertEqual(len(actual), 1 << (dirty + 1))
        labels = [x | (z << width) for x in actual for z in range(1 << width)]
        self.assertEqual(_affine_rank(labels), width + dirty + 1)
        self.assertEqual(2 * width - _affine_rank(labels), m + p - 1)
        # Z masks, including Y representatives, change no conjugation of
        # a diagonal sign; both Pauli directions on work remain available.
        self.assertIn(1 << (p + m), actual)

    def test_exact_lookup_pauli_expansion_and_isotropic_difference_space(self):
        for p, dirty in ((2, 0), (2, 1), (3, 0)):
            m, width = 1 << p, p + (1 << p) + dirty
            coefficients = _walsh(_signs(p, dirty))
            expected = {}
            for u in range(m):
                for t in range(m):
                    z = t | (1 << (p + u))
                    expected[z] = Fraction((-1) ** ((u & t).bit_count() & 1), m)
            self.assertEqual({z: value for z, value in enumerate(coefficients) if value},
                             expected)
            support = [z << width for z in expected]
            self.assertEqual(_affine_rank(support), m + p - 1)
            differences = [label ^ support[0] for label in support]
            self.assertTrue(all(not _symplectic(a, b, width)
                                for a in differences for b in differences))
            for label in differences:
                z = label >> width
                self.assertEqual((z >> p).bit_count() & 1, 0)
                self.assertEqual(z >> (p + m), 0)
            # Independent ANF extraction confirms the nonquadratic target
            # premise used by the existing Clifford-distance threshold.
            anf = [(1 - value) // 2 for value in _signs(p, dirty)]
            for bit in range(width):
                for mask in range(len(anf)):
                    if mask >> bit & 1:
                        anf[mask] ^= anf[mask ^ (1 << bit)]
            self.assertEqual(max(mask.bit_count() for mask, value in enumerate(anf) if value),
                             p + 1)

    def test_one_native_band_has_large_support_without_lookup(self):
        p, m, width = 2, 4, 6
        rank = m + p - 1
        # Start with the specified root Z. A paid Clifford makes its X
        # product, followed by a single physical band of rank T gates.
        word = [('H', 0)] + [('CX', 0, j) for j in range(1, rank)]
        word += [('T', j) for j in range(rank)]
        poly, active = {1 << width: CO}, 0
        for gate in word:
            before = poly
            poly = _propagate(poly, width, gate)
            if gate[0] in ('T', 'TDG'):
                active += poly != before
            self.assertLessEqual(_affine_rank(poly), active)
        self.assertEqual([len(band) for band in _native_t_bands(word)], [rank])
        self.assertEqual(active, rank)
        self.assertEqual(_affine_rank(poly), rank)
        self.assertEqual(len(poly), 1 << rank)
        self.assertEqual(len(set(poly.values())), 1)
        # The support rank alone matches the lower-bound target, but the
        # actual source has an anticommuting pair of terms in its support.
        first, second = next(iter(poly)), list(poly)[1]
        self.assertEqual(_symplectic(first, second, width), 1)
        self.assertTrue(all(label & ((1 << width) - 1) for label in poly))

    def test_independent_noncommuting_axes_grow_support_exactly_when_active(self):
        width = 3
        # These five independent axes are X0, Z0, X0 X1, Z1, Z2.
        # The first four act on the evolving probe; the final Z2 is silent.
        axes = [1, 1 << width, 3, 1 << (width + 1), 1 << (width + 2)]
        self.assertEqual(_affine_rank([0, *axes]), len(axes))
        self.assertEqual(_symplectic(axes[0], axes[1], width), 1)
        modules = [
            [('H', 0), ('T', 0), ('H', 0)],
            [('T', 0)],
            [('H', 0), ('H', 1), ('CX', 0, 1), ('T', 1),
             ('CX', 0, 1), ('H', 1), ('H', 0)],
            [('T', 1)], [('T', 2)],
        ]
        poly, active = {1 << width: CO}, []
        for module in modules:
            before_module = poly
            for gate in module:
                before, old_rank = poly, _affine_rank(poly)
                poly = _propagate(poly, width, gate)
                if gate[0] == 'T':
                    changed = poly != before
                    active.append(changed)
                    self.assertEqual(_affine_rank(poly), old_rank + changed)
                else:
                    self.assertEqual(_affine_rank(poly), old_rank)
            if module == modules[-1]:
                self.assertEqual(poly, before_module)
        self.assertEqual(active, [True, True, True, True, False])
        self.assertEqual(_affine_rank(poly), 4)


if __name__ == '__main__':
    unittest.main()
