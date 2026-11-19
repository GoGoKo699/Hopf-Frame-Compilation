"""Finite checks of the two-wide-layer and layer-profile counting bounds."""

from fractions import Fraction
from itertools import combinations
from random import Random
import unittest

from tests.test_clifford_orbit_query_rank import _rank


def _pair(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width))
            ^ ((right & mask) & (left >> width))).bit_count() % 2


def _gram_rank(axes, width):
    return _rank(sum(_pair(a, b, width) << i for i, b in enumerate(axes))
                 for a in axes)


def _profiles(total, maximum):
    if total == 0:
        yield ()
    for k in range(min(total, maximum), 0, -1):
        for rest in _profiles(total - k, k):
            yield (k, *rest)


class NativeLayerCountDistributionTests(unittest.TestCase):
    def test_isotropic_layers_and_augmented_rank(self):
        rng = Random(773)
        for width in range(1, 8):
            for _ in range(50):
                layers = []
                for _ in range(rng.randint(1, 8)):
                    axes = [1 << i for i in range(rng.randint(1, width))]
                    for _ in range(4 * width):
                        pivot = rng.randrange(1, 1 << (2 * width))
                        axes = [a ^ pivot if _pair(a, pivot, width) else a for a in axes]
                    layers.append(axes)
                axes = [a for layer in layers for a in layer]
                t = len(axes)
                v = _rank(axes)
                twice_r = _gram_rank(axes, width)
                largest = sorted(map(len, layers), reverse=True)[:2]
                s = t - sum(largest)
                self.assertLessEqual(twice_r, 2 * (v - largest[0]))
                self.assertLessEqual(twice_r, s + max(0, 2 * v - t))
                boundary = [rng.randrange(1 << (2 * width)) for _ in range(6)]
                self.assertLessEqual(_gram_rank(axes + boundary, width),
                                     twice_r + 2 * len(boundary))

    def test_commuting_layer_encoding_repeated_signed_affine_axes(self):
        rng = Random(931)
        for width in range(1, 7):
            for _ in range(15):
                parities = [(rng.randrange(1 << width), rng.choice((-1, 1)), rng.randrange(2))
                            for _ in range(30)]
                coefs = [0] * (1 << width)
                for mask, sign, offset in parities:
                    coefs[0] += sign * offset
                    bits = [1 << j for j in range(width) if mask & (1 << j)]
                    for degree in range(1, 4):
                        for subset in combinations(bits, degree):
                            coefs[sum(subset)] += sign * (-1)**offset * (-2)**(degree - 1)
                for x in range(1 << width):
                    direct = sum(sign * (((x & mask).bit_count() % 2) ^ offset)
                                 for mask, sign, offset in parities)
                    model = sum(a for mask, a in enumerate(coefs) if x & mask == mask)
                    self.assertEqual(direct % 8, model % 8)

    def test_radial_profile_envelopes_at_all_breakpoints(self):
        for t in range(1, 22):
            for parts in _profiles(t, t):
                upper = min(Fraction(t, 4), Fraction(t, 2) - parts[0])
                if upper < 0:
                    continue
                s = t - sum(parts[:2])
                profile = Fraction(3 * t*t, 8) - sum(Fraction(k*k, 2) for k in parts)
                two = (Fraction(s*t, 2) + Fraction(s*s, 4) if 2*s <= t else
                       Fraction(3*t*t, 8) - Fraction((t-s)**2, 4))
                for r in {Fraction(0), upper, *(Fraction(k) for k in parts if k < upper)}:
                    penalty = sum((k*k - max(0, k-r)**2)/Fraction(2) for k in parts)
                    self.assertLessEqual(2*r*t - 2*r*r - penalty, min(profile, two))


if __name__ == '__main__':
    unittest.main()
