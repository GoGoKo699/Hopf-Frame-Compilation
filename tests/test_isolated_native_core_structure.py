"""Bounded checks for isolated native cores at arbitrary control count.

The arbitrary-width theorem is analytic. These checks retain native scalar
phases and exercise signs involving more than two controls; they do not
search for optimal depth or infer C3 membership from workspace return.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

from tests.exact_matrices import CO, adjoint, native_columns
from tests.test_native_band_source_structure import (
    _band, _flatten, _inverse, _is_c3, _ledger, _native_phase, _pair, _rank,
)
from tests.test_outer_axis_normalizer import _image, _tau


def _anf_support(values):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for mask in range(len(values)):
            if mask >> bit & 1:
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return {mask for mask, value in enumerate(coefficients) if value}


def _core_schedules(index):
    control, target = 2 * index, 2 * index + 1
    b = [('C', [('CX', target, control)]),
         ('T', [('T', control), ('T', target)]),
         ('C', [('CX', target, control)])]
    c = [('C', [('H', target)])] + b + [('C', [('H', target)])]
    original = c + b + c
    middle = [('H', target), ('CX', control, target), ('H', target),
              ('S', target), ('S', target)]
    pre = [('S', target)] * 3 + [('H', target), ('CX', target, control)]
    post = [('CX', target, control), ('H', target), ('S', target)]
    phase = [('C', pre), ('T', [('T', control), ('T', target)]), ('C', post)]
    inverse_phase = [('C', pre), ('T', [('TDG', control), ('TDG', target)]),
                     ('C', post)]
    forward = [('C', middle)] + phase + [('C', [('S', control)])]
    inverse = [('C', [('S', control)] * 3)] + inverse_phase + [('C', middle)]
    return original, forward, inverse


def _parallel_schedules(schedules):
    """Align disjoint cores' T slots, retaining every chronological gate."""
    layers = []
    for schedule in schedules:
        groups = [[]]
        for kind, gates in schedule:
            if kind == 'T':
                groups.extend([gates, []])
            else:
                groups[-1].extend(gates)
        layers.append(groups)
    return [('C' if j % 2 == 0 else 'T',
             [gate for row in layers for gate in row[j]])
            for j in range(len(layers[0]))]


class IsolatedNativeCoreStructure(unittest.TestCase):
    def test_literal_eighth_root_scalar_identity(self):
        omega = native_columns(1, [('T', 0)])[1][1]
        powers = [CO]
        for _ in range(7):
            powers.append(powers[-1] * omega)
        count = 0
        for e, f, ell, u, eigenvalue in itertools.product(
                (-1, 1), (-1, 1), (0, 1), (0, 1), (-1, 1)):
            activation = (1 ^ u) if e == f else u
            left_exp = e * int(((-1) ** ell) * eigenvalue == -1)
            left_exp += f * int(((-1) ** (ell ^ u)) * eigenvalue == -1)
            right_exp = e * ell + f * (ell ^ u)
            right_exp += 2 * e * activation * int(eigenvalue == -1)
            right_exp += 4 * activation * ell * int(eigenvalue == -1)
            self.assertEqual(powers[left_exp % 8], powers[right_exp % 8])
            count += 1
        self.assertEqual(count, 32)

    def test_four_controls_all_affine_parities_and_axis_signs(self):
        count = 0
        for axis, e, f, mask, constant in itertools.product(
                range(1, 4), (-1, 1), (-1, 1), range(16), (0, 1)):
            tables = {probe: [] for probe in range(1, 4)}
            labels = {probe: [set(), set()] for probe in range(1, 4)}
            pauli = _tau(axis, 1)
            for word in range(16):
                own_control = word & 1
                ell = ((word & mask).bit_count() & 1) ^ constant
                paired = (_native_phase(((-1) ** ell) * pauli, e)
                          @ _native_phase(((-1) ** (ell ^ own_control)) * pauli, f))
                activation = (1 ^ own_control) if e == f else own_control
                gamma = np.exp(1j * np.pi * (e * ell + f * (ell ^ own_control)) / 4)
                quarter = _native_phase(pauli, 2 * e)
                expected = (gamma * (quarter if activation else np.eye(2))
                            @ (pauli if activation * ell else np.eye(2)))
                np.testing.assert_allclose(paired, expected, atol=2e-14)
                for probe in tables:
                    label, sign = _image(paired, probe, 1)
                    tables[probe].append(sign)
                    labels[probe][own_control].add(label)
                count += 1
            for probe, values in tables.items():
                self.assertTrue(all(len(group) == 1 for group in labels[probe]))
                support = _anf_support(values)
                self.assertTrue(all(term.bit_count() <= 2 for term in support))
                self.assertTrue(all(term & 1 for term in support if term.bit_count() == 2))
        self.assertEqual(count, 6144)

    def test_three_core_literal_source_inverse_and_fair_baseline(self):
        width = 6
        # Conjugating by two literal CZ words decorates the first target
        # with the parity of BOTH other controls. No wire is initialized.
        a = [1282, 1346, 8, 264, 32, 1056]
        b = [128, 192, 512, 768, 2048, 3072]
        c = a
        for axes in (a, b, c):
            self.assertEqual(_rank(axes), len(axes))
            self.assertTrue(all(_pair(p, q, width) == 0
                                for p, q in itertools.combinations(axes, 2)))
        self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 3)
        original, baseline, baseline_inverse = [
            _parallel_schedules([_core_schedules(i)[j] for i in range(3)])
            for j in range(3)]
        wrapper = [('H', 1), ('CX', 2, 1), ('H', 1),
                   ('H', 1), ('CX', 4, 1), ('H', 1)]

        def wrap(schedule):
            return [('C', wrapper)] + schedule + [('C', list(reversed(wrapper)))]

        original, baseline, baseline_inverse = map(
            wrap, (original, baseline, baseline_inverse))
        self.assertEqual(_ledger(original), (18, 42, 3))
        self.assertEqual(_ledger(_inverse(original)), (18, 42, 3))
        self.assertEqual(_ledger(baseline), (6, 54, 1))
        self.assertEqual(_ledger(baseline_inverse), (6, 60, 1))
        source = native_columns(width, _flatten(original))
        self.assertEqual(native_columns(width, _flatten(baseline)), source)
        inverse = native_columns(width, _flatten(_inverse(original)))
        self.assertEqual(inverse, adjoint(source))
        self.assertEqual(native_columns(width, _flatten(baseline_inverse)), inverse)

    def test_single_orbit_quotient_derivative_rank(self):
        for width in (2, 4, 6, 8, 10):
            columns = [1 << j for j in range(width - 1)]
            columns.append((1 << (width - 1)) - 1)
            self.assertEqual(_rank(columns), width - 1)
            # With this one-relation basis, the derivative is I + 11^T.
            matrix = (np.eye(width - 1, dtype=np.uint8)
                      ^ np.ones((width - 1, width - 1), dtype=np.uint8))
            labels = [sum(int(matrix[row, col]) << row
                          for row in range(width - 1))
                      for col in range(width - 1)]
            self.assertEqual(_rank(labels), width - 2)

    def test_weight_four_orbit_activation_does_not_imply_c3(self):
        width = 4
        b = [16, 32, 64, 128]
        control = 240
        c = [15, 255]
        even = [z << width for z in range(16) if z.bit_count() % 2 == 0]
        count = 0
        for offsets in itertools.combinations(even, 4):
            a = [15 ^ d for d in offsets]
            if _rank(a) != 4 or a[0] ^ a[1] ^ a[2] ^ a[3] != control:
                continue
            for axes in (a, b, c):
                self.assertEqual(_rank(axes), len(axes))
                self.assertTrue(all(_pair(p, q, width) == 0
                                    for p, q in itertools.combinations(axes, 2)))
            self.assertEqual(_rank(a) + _rank(b) - _rank(a + b), 3)
            self.assertFalse(_is_c3(_band(c, width) @ _band(b, width)
                                    @ _band(a, width), width))
            count += 1
        self.assertEqual(count, 8)

    def test_coupled_rows_localize_all_involutory_corrections(self):
        rows = [0b001111, 0b111100]
        b_width = 6
        reduced = []
        for row in rows:
            for pivot, old in reduced:
                if row >> pivot & 1:
                    row ^= old
            reduced.append(((row & -row).bit_length() - 1, row))
        free = [j for j in range(b_width)
                if j not in [pivot for pivot, _ in reduced]]
        width = len(free)

        def project(word):
            for pivot, row in reduced:
                if word >> pivot & 1:
                    word ^= row
            return sum(((word >> p) & 1) << j for j, p in enumerate(free))

        beta = [project(1 << j) << width for j in range(b_width)]
        c = [sum(((row >> p) & 1) << j for j, p in enumerate(free))
             for row in rows]
        self.assertEqual(_rank(beta), b_width - len(rows))
        self.assertTrue(all(_pair(c[j], beta[i], width) == (rows[j] >> i & 1)
                            for i in range(b_width) for j in range(len(rows))))
        identity = tuple(1 << j for j in range(2 * width))

        def image(matrix, value):
            result = 0
            for j, column in enumerate(matrix):
                if value >> j & 1:
                    result ^= column
            return result

        def product(left, right):
            return tuple(image(left, column) for column in right)

        def transvection(axis):
            return tuple(e ^ (axis if _pair(axis, e, width) else 0)
                         for e in identity)

        active_count = correction_count = uniform_count = even_count = 0
        for u in range(1 << b_width):
            active = [j for j, row in enumerate(rows) if (u & row).bit_count() % 2]
            p = q = identity
            for j in active:
                p = product(transvection(c[j]), p)
            for i, axis in enumerate(beta):
                if u >> i & 1:
                    q = product(transvection(axis), q)
            base = product(p, q)
            ws = []
            for j in active:
                w = 0
                for i, axis in enumerate(beta):
                    if (u & rows[j]) >> i & 1:
                        w ^= axis
                ws.append(w)
            square = product(base, base)
            self.assertEqual(_rank([a ^ b for a, b in zip(square, identity)]),
                             2 * _rank(ws))
            if not active:
                continue
            active_count += 1
            c_sum = 0
            for j in active:
                c_sum ^= c[j]
            for h in range(1 << (2 * width)):
                corrected = product(base, transvection(h))
                if product(corrected, corrected) != identity:
                    continue
                correction_count += 1
                even_count += len(active) % 2 == 0
                self.assertEqual(_rank(ws), 1)
                self.assertNotEqual(h, 0)
                self.assertIn(h, (ws[0], c_sum, ws[0] ^ c_sum))
                # Native A isotropy requires one bit across the whole class.
                pairings = {_pair(h, axis, width) for axis in beta}
                if len(pairings) == 1 and h != ws[0]:
                    uniform_count += 1
        self.assertEqual(active_count, 48)
        self.assertEqual(correction_count, 104)
        self.assertGreater(even_count, 0)
        self.assertEqual(uniform_count, 0)

    def test_zero_axis_incidence_with_zero_diagonal_gram_pairing(self):
        width = 3
        # H=(1,1,0) gives one paired nonzero class and a zero B axis.
        # K=(1,1,1) also makes that zero axis touch both C members.
        beta = [1 << width, 1 << width, 1 << (width + 1)]
        c = 0b011
        identity = tuple(1 << j for j in range(2 * width))

        def image(matrix, value):
            result = 0
            for j, column in enumerate(matrix):
                if value >> j & 1:
                    result ^= column
            return result

        def product(left, right):
            return tuple(image(left, column) for column in right)

        def transvection(axis):
            return tuple(e ^ (axis if _pair(axis, e, width) else 0)
                         for e in identity)

        spectators = (identity, transvection(1 << (width + 2)))
        active_count = correction_count = zero_diagonal_count = compatible_count = 0
        for spectator in spectators:
            self.assertTrue(all(image(spectator, axis) == axis
                                for axis in beta + [c]))
            for u in range(8):
                if (u & 0b011).bit_count() % 2 == 0:
                    continue
                active_count += 1
                q, w = identity, 0
                for j, axis in enumerate(beta):
                    if u >> j & 1:
                        q = product(transvection(axis), q)
                        w ^= axis
                self.assertNotEqual(w, 0)
                base = product(transvection(c), q)
                inverse_base = product(q, transvection(c))
                self.assertEqual(_rank([a ^ b for a, b in zip(base, inverse_base)]), 2)
                diagonal_pairing = _pair(c, w, width)
                for h in range(1 << (2 * width)):
                    corrected = product(product(base, transvection(h)), spectator)
                    if product(corrected, corrected) != identity:
                        continue
                    correction_count += 1
                    zero_diagonal_count += diagonal_pairing == 0
                    self.assertNotEqual(h, 0)
                    self.assertIn(h, (w, c, w ^ c))
                    # h=w has a nonzero class character and is excluded
                    # separately. Other native A images must fix zero B.
                    if h != w and _pair(h, beta[-1], width) == 0:
                        compatible_count += 1
        self.assertEqual(active_count, 8)
        self.assertEqual(correction_count, 20)
        self.assertEqual(zero_diagonal_count, 8)
        self.assertEqual(compatible_count, 0)


if __name__ == '__main__':
    unittest.main()
