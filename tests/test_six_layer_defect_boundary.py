"""Exact bounded evidence for the next native defect boundary.

The six-layer normalization is an identity for arbitrary native fibres.
The odd-order example concerns its inner Clifford families only: it is
not a returned six-layer query or a compiler, and supplies no lower bound.
All matrix checks retain literal phases and every work column.
"""
from __future__ import annotations

import itertools
import unittest

from tests.exact_matrices import (
    CI, CO, adjoint, eye, multiply, native_columns, scale,
)


def _product(*operators):
    result = eye(len(operators[0]))
    for operator in operators:
        result = multiply(result, operator)
    return result


def _conjugate(unitary, operator):
    return _product(unitary, operator, adjoint(unitary))


def _x_word(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


def _pauli(p, q, width):
    # Literal X^p Z^q, including the phase when X and Z overlap.
    word = []
    for wire in range(width):
        if q >> wire & 1:
            word.extend([('S', wire)] * 2)
        if p >> wire & 1:
            word.extend(_x_word(wire))
    return native_columns(width, word)


def _inverse_word(word):
    inverse = []
    for name, *wires in reversed(word):
        if name == 'S':
            inverse.extend([('S', *wires)] * 3)
        else:
            inverse.append(({'T': 'TDG', 'TDG': 'T'}.get(name, name), *wires))
    return inverse


def _inner_word(u_controls, v_controls, work):
    # Chronological native word for (T† X^u T)(H T X^v T† H).
    return ([('H', work), ('TDG', work)]
            + [('CX', control, work) for control in v_controls]
            + [('T', work), ('H', work), ('T', work)]
            + [('CX', control, work) for control in u_controls]
            + [('TDG', work)])


class SixLayerDefectBoundary(unittest.TestCase):
    def test_literal_six_layer_comparison_and_three_band_pullback(self):
        width = 2
        constants = [
            [('H', 0), ('CX', 0, 1), ('S', 1)],
            [('S', 0), ('H', 1)],
            [('H', 0), ('CX', 1, 0)],
            [('CX', 0, 1), ('H', 1), ('S', 0)],
            [('H', 0), ('S', 1), ('CX', 1, 0)],
            [('S', 0), ('CX', 0, 1), ('H', 1)],
            [('H', 0), ('CX', 1, 0), ('S', 1)],
        ]
        c = [native_columns(width, word) for word in constants]
        d = [eye(4)] + [native_columns(width, [
            ('T' if j & 1 else 'TDG', 0),
            ('TDG' if j % 3 else 'T', 1)]) for j in range(1, 7)]
        rows = [((1, 2), (3, 1)), ((3, 2), (1, 2)),
                ((1, 3), (2, 1)), ((2, 3), (1, 3)),
                ((3, 1), (2, 3)), ((2, 1), (3, 2)),
                ((3, 2), (2, 1))]
        phases = (CO, CI, -CO, -CI)

        def p(j, address):
            labels = [sum(((row & address).bit_count() & 1) << wire
                          for wire, row in enumerate(part)) for part in rows[j]]
            return _pauli(*labels, width)

        def phase(j, address):
            u, v = address & 1, address >> 1
            return phases[((j + 1) * u + j * v + 2 * u * v) % 4]

        def fibre(address):
            value = scale(_product(p(0, address), c[0]), phase(0, address))
            for j in range(1, 7):
                value = scale(_product(p(j, address), c[j], d[j], value),
                              phase(j, address))
            return value

        zero = fibre(0)
        self.assertNotEqual(zero, eye(4))
        self.assertEqual(_product(zero, adjoint(zero)), eye(4))
        w = _product(c[6], d[6])
        prefix_a = _product(c[5], d[5])
        prefix_b = _product(prefix_a, c[4], d[4])
        prefix = _product(prefix_b, c[3], d[3])
        s = _product(c[3], d[3])
        e = _product(c[2], d[2])
        pauligroup = [scale(_pauli(p0, q0, width), scalar)
                     for p0, q0, scalar in itertools.product(range(4), range(4), phases)]
        generators = [_pauli(1 << i, 0, width) for i in range(width)]
        generators += [_pauli(0, 1 << i, width) for i in range(width)]
        for address in range(4):
            a = _conjugate(prefix_a, p(4, address))
            b = _conjugate(prefix_b, p(3, address))
            j = _conjugate(prefix, p(2, address))
            o = _conjugate(_product(prefix, e), p(1, address))
            r = _conjugate(adjoint(c[0]), p(0, address))
            scalar = CO
            for band in range(7):
                scalar *= phase(band, address)
            actual = fibre(address)
            normalized = scale(_product(
                p(6, address), w, p(5, address), a, b, j, o,
                adjoint(w), zero, r), scalar)
            self.assertEqual(actual, normalized)
            # This uses the actual zero fibre's adjoint, not a target value.
            extracted = scale(_product(
                adjoint(a), adjoint(p(5, address)), adjoint(w),
                adjoint(p(6, address)), actual, adjoint(r),
                adjoint(zero), w), scalar.conj())
            left = _conjugate(adjoint(s), p(3, address))
            right = _conjugate(e, p(1, address))
            inner = _product(left, p(2, address), right)
            self.assertEqual(extracted, _product(b, j, o))
            self.assertEqual(extracted, _conjugate(prefix, inner))
            for factor in (left, p(2, address), right):
                self.assertEqual(_product(factor, adjoint(factor)), eye(4))
                for probe in generators:
                    self.assertIn(_conjugate(factor, probe), pauligroup)

    def test_first_inner_defect_can_have_symplectic_order_three(self):
        h = native_columns(1, [('H', 0)])
        t = native_columns(1, [('T', 0)])
        x = _pauli(1, 0, 1)
        a = _conjugate(adjoint(t), x)
        y = _conjugate(_product(h, t), x)
        values = [_product(a if address & 1 else eye(2),
                           y if address & 2 else eye(2)) for address in range(4)]
        # x=e2 and y=e1 in Z(x,y)=L(x+y)L(y)^-1 L(x)^-1.
        defect = _product(values[3], adjoint(values[1]), adjoint(values[2]))
        self.assertEqual(defect, _product(a, y, a, y))
        hermitian = [eye(2), x, _pauli(0, 1, 1), scale(_pauli(1, 1, 1), CI)]
        permutation = []
        for probe in hermitian:
            image = _conjugate(defect, probe)
            matches = [label for label, pauli in enumerate(hermitian)
                       if image in (pauli, scale(pauli, -CO))]
            self.assertEqual(len(matches), 1)
            permutation.append(matches[0])
        self.assertNotEqual(permutation, list(range(4)))
        self.assertEqual([permutation[permutation[permutation[j]]]
                          for j in range(4)], list(range(4)))
        self.assertFalse(any(defect == scale(pauli, scalar)
                             for pauli, scalar in itertools.product(
                                 hermitian, (CO, CI, -CO, -CI))))
        # Fair fixed-branch baseline: the literal defect is Clifford-only.
        # This 12G word is an upper bound, not a new optimal synthesis claim.
        baseline = [('S', 0), ('H', 0), ('S', 0), ('S', 0), ('S', 0),
                    ('H', 0), ('S', 0), ('S', 0), ('H', 0), ('S', 0),
                    ('S', 0), ('H', 0)]
        self.assertEqual(native_columns(1, baseline), defect)
        self.assertEqual((len(baseline), sum(g[0] in ('T', 'TDG') for g in baseline)),
                         (12, 0))

    def test_coherent_native_family_actual_inverse_and_gate_ledgers(self):
        # Two address controls, one arbitrary work wire, no clean/helper wire.
        word = _inner_word([0], [1], 2)
        actual = native_columns(3, word)
        inverse_word = _inverse_word(word)
        self.assertEqual(native_columns(3, inverse_word), adjoint(actual))
        self.assertEqual(_product(actual, adjoint(actual)), eye(8))
        self.assertEqual((sum(g[0] in ('T', 'TDG') for g in word),
                          sum(g[0] not in ('T', 'TDG') for g in word)), (4, 4))
        self.assertEqual(len(word), len(inverse_word))
        self.assertTrue(all(g[2] == 2 for g in word if g[0] == 'CX'))
        h = native_columns(1, [('H', 0)])
        t = native_columns(1, [('T', 0)])
        a = _conjugate(adjoint(t), _pauli(1, 0, 1))
        y = _conjugate(_product(h, t), _pauli(1, 0, 1))
        fibres = [_product(a if address & 1 else eye(2),
                           y if address & 2 else eye(2)) for address in range(4)]
        for address, source, target in itertools.product(range(4), range(2), range(2)):
            self.assertEqual(actual[address | (target << 2)][address | (source << 2)],
                             fibres[address][target][source])
        # Four read-only address wires hold x=(0,1), y=(2,3).
        # Actual inverse order implements L(x+y)L(y)^-1 L(x)^-1.
        defect_word = (_inverse_word(_inner_word([0], [1], 4))
                       + _inverse_word(_inner_word([2], [3], 4))
                       + _inner_word([0, 2], [1, 3], 4))
        self.assertEqual((sum(g[0] in ('T', 'TDG') for g in defect_word),
                          sum(g[0] not in ('T', 'TDG') for g in defect_word)), (12, 14))
        # All T gates use the same work wire, giving the emitted upper
        # depths 4 and 12; fixed branches may simplify, as checked above.
        self.assertTrue(all(g[1] == 4 for g in defect_word if g[0] in ('T', 'TDG')))
        defect_native = native_columns(5, defect_word)
        for xx, yy in itertools.product(range(4), repeat=2):
            address = xx | (yy << 2)
            expected = _product(fibres[xx ^ yy], adjoint(fibres[yy]), adjoint(fibres[xx]))
            for source, target in itertools.product(range(2), repeat=2):
                self.assertEqual(defect_native[address | (target << 4)][address | (source << 4)],
                                 expected[target][source])


if __name__ == '__main__':
    unittest.main()
