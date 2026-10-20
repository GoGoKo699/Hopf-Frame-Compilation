"""Small literal native receipt for the fermionic moving-bank identity.

The six-wire case checks every occupation sector and the actual inverse.
It supports the analytic identity, not the asymptotic depth obstructions.
"""
from __future__ import annotations

import unittest

from tests.exact_matrices import CO, CZ, native_columns


def _inverse(word):
    result = []
    for name, *wires in reversed(word):
        if name == 'S':
            result.extend([('S', *wires)] * 3)
        else:
            result.append(({'T': 'TDG', 'TDG': 'T'}.get(name, name),
                           *wires))
    # Reduce only adjacent powers of the same S; cross no other gate.
    reduced = []
    for gate in result:
        reduced.append(gate)
        if gate[0] == 'S' and reduced[-4:] == [gate] * 4:
            del reduced[-4:]
    return reduced


def _cz(a, b):
    return [('H', b), ('CX', a, b), ('H', b)]


def _fswap(a, b):
    return [('CX', a, b), ('CX', b, a), ('CX', a, b)] + _cz(a, b)


def _walsh_layers(b, offset):
    m = 1 << b
    result = []
    for bit in range(b):
        pairs = [(u, u ^ (1 << bit)) for u in range(m)
                 if not u >> bit & 1]
        order = list(range(m))
        desired = [u for pair in pairs for u in pair]
        route = []
        for position, wanted in enumerate(desired):
            cursor = order.index(wanted)
            while cursor > position:
                route.extend(_fswap(offset + cursor - 1, offset + cursor))
                order[cursor - 1], order[cursor] = (
                    order[cursor], order[cursor - 1])
                cursor -= 1
        before, layer, after = list(route), [], []
        for i in range(0, m, 2):
            a, c = offset + i, offset + i + 1
            align = [('CX', a, c), ('H', a), ('CX', a, c), ('S', c)]
            before.extend([('S', c), ('S', c)] + _inverse(align))
            layer.extend([('T', a), ('TDG', c)])
            after.extend(align)
        after.extend(_inverse(route))
        result.extend([('C', before), ('T', layer), ('C', after)])
    return result


class BankMotionWalshIndicator(unittest.TestCase):
    def test_six_wire_complete_indicator_and_actual_inverse(self):
        b, m = 2, 4
        schedule = _walsh_layers(b, b)
        layers = [word for kind, word in schedule if kind == 'T']
        self.assertEqual(len(layers), b)
        self.assertEqual(sum(map(len, layers)), b * m)
        for layer in layers:
            self.assertEqual(len({gate[1] for gate in layer}), m)
        transform = [gate for _, word in schedule for gate in word]
        inverse = _inverse(transform)
        self.assertEqual(sum(g[0] not in {'T', 'TDG'} for g in transform), 60)
        self.assertEqual(sum(g[0] not in {'T', 'TDG'} for g in inverse), 60)
        diagonal = [gate for i in range(b) for u in range(m)
                    if u >> i & 1 for gate in _cz(i, b + u)]
        routing = transform + diagonal + inverse
        bank_h = [('H', b + u) for u in range(m)]
        word = bank_h + routing + [('S', b), ('S', b)] + _inverse(routing) + bank_h
        self.assertEqual(sum(g[0] in {'T', 'TDG'} for g in word), 32)
        self.assertEqual(sum(g[0] not in {'T', 'TDG'} for g in word), 258)
        for actual_word in (word, _inverse(word)):
            actual = native_columns(b + m, actual_word)
            for column in range(1 << (b + m)):
                address = column & (m - 1)
                row = column ^ (1 << (b + address))
                expected = [CO if i == row else CZ
                            for i in range(1 << (b + m))]
                self.assertEqual([actual[i][column]
                                  for i in range(1 << (b + m))], expected)


if __name__ == '__main__':
    unittest.main()
