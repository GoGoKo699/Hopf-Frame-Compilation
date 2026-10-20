"""Exact, bounded native fixture for a whole-bank special-angle call.

This is a detuning component, not a selective indicator. Integer phase
exponents verify every basis column, including literal phases and the
actual inverse, without a large matrix or an asymptotic simulation claim.
"""
from __future__ import annotations

import unittest

try:
    from .test_t_depth import _ccz_batch
except ImportError:
    from test_t_depth import _ccz_batch


def _positive_s(word):
    return [native for gate in word
            for native in ([('S', gate[1])] * 3 if gate[0] == 'SDG' else [gate])]


def _inverse(word):
    result = []
    for gate in reversed(word):
        if gate[0] == 'S':
            result += [('S', gate[1])] * 3
        elif gate[0] in ('T', 'TDG'):
            result.append(('TDG' if gate[0] == 'T' else 'T', gate[1]))
        else:
            assert gate[0] == 'CX'
            result.append(gate)
    return result


def _cz(left, right):
    # 2a+2b-2(a XOR b) = 4ab modulo eight; positive S only.
    return ([('S', left), ('S', right), ('CX', left, right)]
            + [('S', right)] * 3 + [('CX', left, right)])


def _fixture():
    # b=2, m=4: x=(0,1), Y=(2,3,4,5), signal=(6,7,8,9).
    x, bank, signal = (0, 1), tuple(range(2, 6)), tuple(range(6, 10))
    signs = []
    for u, s in enumerate(signal):
        for j, address in enumerate(x):
            if (u >> j) & 1:
                signs += ([('CX', address, s)] + [('S', s)] * 3
                          + [('CX', address, s)])
    # T_x^m = T_x^4 = Z_x. The exact scalar is omega^8 = 1.
    signs += [('S', address) for address in x for _ in range(2)]
    local = [('CX', s, y) for s, y in zip(signal, bank)]
    schedule = [('C', signs + local),
                ('T', [('TDG', wire) for wire in bank + signal]),
                ('C', list(reversed(local)))]
    px, ps = x[0], signal[0]
    extract_x = [('CX', x[1], px)]
    extract_s = [('CX', s, ps) for s in signal[1:]]
    schedule += [('C', extract_x + extract_s),
                 ('T', [('TDG', px), ('TDG', ps)]),
                 ('C', [('CX', px, ps)]),
                 ('T', [('T', ps)]),
                 ('C', [('CX', px, ps)])]
    # In parity coordinates Q_X = px*x1 XOR x1.
    schedule += _ccz_batch(ps, [(px, x[1])])
    schedule += [('C', _cz(ps, x[1]) + list(reversed(extract_s)))]
    # For original signals a,b,c,d, set a'=a+c+d, b'=b+c+d.
    # Then Q_S = a'b' XOR cd XOR c XOR d.
    a, b, c, d = signal
    quadratic = [('CX', c, a), ('CX', d, a),
                 ('CX', c, b), ('CX', d, b)]
    schedule += [('C', quadratic)]
    schedule += _ccz_batch(px, [(a, b), (c, d)])
    schedule += [('C', _cz(px, c) + _cz(px, d)
                  + list(reversed(quadratic)) + list(reversed(extract_x)))]
    return [(kind, _positive_s(gates)) for kind, gates in schedule]


def _basis_action(word, basis):
    exponent = 0
    for gate in word:
        name, wire = gate[:2]
        if name == 'CX':
            basis ^= ((basis >> wire) & 1) << gate[2]
        else:
            exponent += {'S': 2, 'T': 1, 'TDG': -1}[name] * ((basis >> wire) & 1)
    return basis, exponent % 8


class JointDetuningResynthesisTests(unittest.TestCase):
    def test_eleven_native_bands_literal_all_input_word_and_actual_inverse(self):
        schedule = _fixture()
        self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 11)
        self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'), 30)
        for kind, gates in schedule:
            allowed = ('T', 'TDG') if kind == 'T' else ('CX', 'S')
            self.assertTrue(all(gate[0] in allowed for gate in gates))
            if kind == 'T':
                targets = [gate[1] for gate in gates]
                self.assertEqual(len(targets), len(set(targets)))
            self.assertTrue(all(0 <= q < 10 for gate in gates for q in gate[1:]))
        word = [gate for _, gates in schedule for gate in gates]
        inverse = _inverse(word)
        for basis in range(1 << 10):
            address = basis & 3
            expected = sum(((address ^ u).bit_count() + 1 - ((basis >> (2 + u)) & 1))
                           * (1 - 2 * ((basis >> (6 + u)) & 1))
                           for u in range(4)) % 8
            self.assertEqual(_basis_action(word, basis), (basis, expected))
            self.assertEqual(_basis_action(inverse, basis), (basis, -expected % 8))


if __name__ == '__main__':
    unittest.main()
