"""Bounded mechanism checks for the exact overlap-two source theorem.

These fixtures support the arbitrary-width analytic proof. They neither
search for optimal depth nor infer C3 membership from returned-query action.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

from tests.exact_matrices import adjoint, native_columns
from tests.test_native_band_source_structure import (
    _band, _flatten, _inverse, _is_c3, _ledger, _pair, _rank,
)


def _transvection(label, width):
    size = 2 * width
    vector = np.array([(label >> bit) & 1 for bit in range(size)], dtype=np.uint8)
    functional = np.array([_pair(label, 1 << bit, width)
                           for bit in range(size)], dtype=np.uint8)
    return np.eye(size, dtype=np.uint8) ^ np.outer(vector, functional)


class TwoOverlapSourceStructure(unittest.TestCase):
    def test_activation_relation_types_and_triple_orbit_obstruction(self):
        planes = {frozenset((r, s, r ^ s))
                  for r in range(1, 32) for s in range(r + 1, 32)}
        seen = set()
        for relation in planes:
            zero = {i for i in range(5) if 1 << i in relation}
            classes = []
            for i in range(5):
                if i in zero:
                    continue
                for group in classes:
                    if (1 << i) ^ (1 << group[0]) in relation:
                        group.append(i)
                        break
                else:
                    classes.append([i])
            repeated = tuple(sorted(len(group) for group in classes if len(group) > 1))
            kind = (len(zero), repeated)
            seen.add(kind)
            self.assertLessEqual(len(zero) + sum(m - 1 for m in repeated), 2)
        self.assertEqual(len(planes), 155)
        self.assertEqual(seen, {(0, ()), (0, (2,)), (0, (3,)), (0, (2, 2)),
                               (1, ()), (1, (2,)), (2, ())})
        # Two orbit differences in a repeated triple violate C isotropy.
        pairs = [0b011, 0b101, 0b110]
        for left, right in itertools.combinations(pairs, 2):
            self.assertEqual((left & right).bit_count() & 1, 1)

    def test_full_space_involution_obstruction_including_radical_couplings(self):
        count = 0
        for width in (2, 3):
            # Target is qubit 0; h=Z1 is a nonzero radical. A further
            # spectator is included at width3. H may anticommute with h.
            a, c, h = 1 << width, 1, 1 << (width + 1)
            b = a ^ h
            m = (_transvection(c, width) @ _transvection(b, width)
                 @ _transvection(a, width)) % 2
            identity = np.eye(2 * width, dtype=np.uint8)
            fixed = [np.array([(v >> bit) & 1 for bit in range(2 * width)],
                              dtype=np.uint8) for v in (a, c, h)]
            rest = list(range(1, width))
            entries = [(i, j) for i in rest for j in rest if i <= j]
            for values in itertools.product((0, 1), repeat=len(entries)):
                diagonal = identity.copy()
                for (i, j), value in zip(entries, values):
                    diagonal[width + i, j] = diagonal[width + j, i] = value
                self.assertTrue(all(np.array_equal(diagonal @ v % 2, v) for v in fixed))
                for label in range(1 << (2 * width)):
                    if _pair(label, a, width):
                        continue
                    action = m @ _transvection(label, width) @ diagonal % 2
                    self.assertFalse(np.array_equal(action @ action % 2, identity))
                    count += 1
        self.assertEqual(count, 272)

    def test_two_controlled_cores_and_cross_control_signs(self):
        width = 4
        b, c = [32, 48, 128, 192], [2, 18, 8, 72]
        for alpha, beta in itertools.product((0, 1), repeat=2):
            a = [2 ^ (64 * alpha), 18 ^ (64 * alpha),
                 8 ^ (16 * beta), 72 ^ (16 * beta)]
            self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 2)
            self.assertEqual(_rank(a), len(a))
            self.assertTrue(all(_pair(p, q, width) == 0
                                for p, q in itertools.combinations(a, 2)))
            source = _band(c, width) @ _band(b, width) @ _band(a, width)
            # These explicit cross-control choices violate the full C3
            # premise; the relation budget alone is not a certificate.
            self.assertEqual(_is_c3(source, width), alpha == beta == 0)

    def test_nonduplicate_zero_and_mixed_cores_fail_c3_for_all_orientations(self):
        fixtures = [([48, 56, 16], [32, 40], [4, 20]),
                    ([16, 8], [24, 48, 56], [2, 10])]
        count = 0
        for a, b, c in fixtures:
            self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 2)
            for labels in (a, b, c):
                self.assertEqual(_rank(labels), len(labels))
                self.assertTrue(all(_pair(p, q, 3) == 0
                                    for p, q in itertools.combinations(labels, 2)))
            for signs in itertools.product((-1, 1), repeat=len(a) + len(b) + len(c)):
                cut1, cut2 = len(a), len(a) + len(b)
                source = (_band(c, 3, signs[cut2:]) @ _band(b, 3, signs[cut1:cut2])
                          @ _band(a, 3, signs[:cut1]))
                self.assertFalse(_is_c3(source, 3))
                count += 1
        self.assertEqual(count, 256)

    def test_two_zero_duplicates_keep_all_columns_inverse_and_fair_baseline(self):
        # C axes X0 Z2 and X1 Z2 are independent and commute. Both zero B
        # directions have nontrivial native C spread, but A/B cancel exactly.
        pre = [('H', 0), ('H', 1), ('CX', 2, 0), ('CX', 2, 1)]
        post = list(reversed(pre))
        outer = [('C', pre), ('T', [('T', 0), ('T', 1)]), ('C', post)]
        source = [('T', [('T', 0), ('T', 1)]),
                  ('T', [('TDG', 0), ('TDG', 1)])] + outer
        self.assertEqual(_ledger(source), (6, 8, 3))
        self.assertEqual(_ledger(_inverse(source)), (6, 8, 3))
        self.assertEqual(_ledger(outer), (2, 8, 1))
        exact = native_columns(3, _flatten(source))
        self.assertEqual(exact, native_columns(3, _flatten(outer)))
        self.assertEqual(native_columns(3, _flatten(_inverse(source))), adjoint(exact))
        self.assertEqual(native_columns(3, _flatten(_inverse(outer))), adjoint(exact))
        self.assertTrue(_is_c3(_band([33, 34], 3), 3))

    def test_mixed_duplicate_retains_noncommuting_pair_source(self):
        # A tensor factor supplies the one-overlap noncommuting source; a
        # separate zero B axis cancels while its outer C axis couples dirty q3.
        a, b, c = [2, 18, 64], [32, 48, 64], [2, 18, 132]
        self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 2)
        source = _band(c, 4) @ _band(b, 4, [1, 1, -1]) @ _band(a, 4)
        simplified = _band(c, 4) @ _band(b[:2], 4) @ _band(a[:2], 4)
        np.testing.assert_allclose(source, simplified, atol=3e-14)
        np.testing.assert_allclose(source.conj().T, simplified.conj().T, atol=3e-14)
        self.assertTrue(_is_c3(source, 4))
        self.assertFalse(_is_c3(_band(c[:2], 4) @ _band(b[:2], 4), 4))
        first = [('C', [('H', 1), ('CX', 1, 0)]),
                 ('T', [('T', 0), ('T', 1), ('T', 2)]),
                 ('C', [('CX', 1, 0), ('H', 1)])]
        middle = [('C', [('CX', 1, 0)]),
                  ('T', [('T', 0), ('T', 1), ('TDG', 2)]),
                  ('C', [('CX', 1, 0)])]
        pre = [('H', 1), ('H', 2), ('CX', 1, 0), ('CX', 3, 2)]
        outer = [('C', pre), ('T', [('T', 0), ('T', 1), ('T', 2)]),
                 ('C', list(reversed(pre)))]
        displayed = first + middle + outer
        self.assertEqual(_ledger(displayed), (9, 14, 3))
        self.assertEqual(_ledger(_inverse(displayed)), (9, 14, 3))
        exact = native_columns(4, _flatten(displayed))
        np.testing.assert_allclose([[v.display() for v in row] for row in exact],
                                   source, atol=4e-14)
        self.assertEqual(native_columns(4, _flatten(_inverse(displayed))), adjoint(exact))

        # The exact one-overlap fair baseline shares its T slot with the
        # disjoint dirty outer phase. Both forward and inverse are charged.
        c0 = [('H', 1), ('CX', 0, 1), ('H', 1), ('S', 1), ('S', 1)]
        phase_pre = [('S', 1)]*3 + [('H', 1), ('CX', 1, 0),
                                  ('H', 2), ('CX', 3, 2)]
        phase_post = [('CX', 3, 2), ('H', 2), ('CX', 1, 0), ('H', 1), ('S', 1)]
        phases = [('C', phase_pre), ('T', [('T', 0), ('T', 1), ('T', 2)]),
                  ('C', phase_post)]
        baseline = [('C', c0)] + phases + [('C', [('S', 0)])]
        inverse_phases = [('C', phase_pre),
                          ('T', [('TDG', 0), ('TDG', 1), ('TDG', 2)]),
                          ('C', phase_post)]
        baseline_inverse = [('C', [('S', 0)]*3)] + inverse_phases + [('C', c0)]
        self.assertEqual(_ledger(baseline), (3, 18, 1))
        self.assertEqual(_ledger(baseline_inverse), (3, 20, 1))
        self.assertEqual(native_columns(4, _flatten(baseline)), exact)
        self.assertEqual(native_columns(4, _flatten(baseline_inverse)), adjoint(exact))


if __name__ == '__main__':
    unittest.main()
