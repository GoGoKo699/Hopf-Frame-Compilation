"""Exact composition receipt for one shared HRT carrier and all bank outputs.

The source F and its actual inverse are represented by their retained exact
permutations, whose native realization is checked in earlier fixtures. This
is not a new native F receipt. The zeta/copy word is expanded into actual
CNOTs, and every phase column of the 13-wire composition is checked.
"""
import unittest

from tests.test_balanced_algebra_permutation import _hrt_f, _hrt_k


def _inverse(word):
    return [({'F': 'Finv', 'Finv': 'F'}.get(gate[0], gate[0]), *gate[1:])
            for gate in reversed(word)]


def _words(b):
    m, r = 1 << b, b + 1
    bank = list(range(b, b + m))
    carrier = {subset: b + m + subset - 1
               for subset in range(1, 1 << r)}
    zeta = [('CX', bank[subset ^ (1 << bit)], bank[subset])
            for bit in range(b) for subset in range(m)
            if subset & (1 << bit)]
    copy = [('CX', bank[subset], carrier[((1 << r) - 1) ^ subset])
            for subset in range(1, m)]
    shift = zeta + copy + _inverse(zeta)
    jword = [('Finv',)] + [('CX', bit, carrier[1 << bit])
                          for bit in range(b)] + [('F',)]
    cword = jword + [('Z', carrier[(1 << r) - 1])] + _inverse(jword)
    cword += [('Z', carrier[(1 << r) - 1])]
    phase = cword + shift + _inverse(cword) + _inverse(shift)
    # This Z is the phase-form equivalent of the final XOR X correction.
    phase += [('Z', bank[0])]
    return phase, shift, bank, carrier


def _column(word, column, carrier, tables):
    state, phase = column, 0
    carrier_mask = sum(1 << wire for wire in carrier.values())
    for gate in word:
        if gate[0] in tables:
            value = sum(((state >> wire) & 1) << subset
                        for subset, wire in carrier.items())
            value = tables[gate[0]][value >> 1]
            state &= ~carrier_mask
            state |= sum(((value >> subset) & 1) << wire
                         for subset, wire in carrier.items())
        elif gate[0] == 'CX':
            state ^= ((state >> gate[1]) & 1) << gate[2]
        elif gate[0] == 'Z':
            phase ^= (state >> gate[1]) & 1
        else:
            raise ValueError(gate)
    return state, phase


class FullBankCharacterTransportTests(unittest.TestCase):
    def test_shared_carrier_all_phase_columns_and_actual_reverse(self):
        b, r, m = 2, 3, 4
        word, shift, bank, carrier = _words(b)
        width = b + m + (1 << r) - 1
        self.assertEqual(width, 13)
        self.assertEqual(len(shift), b * m + m - 1)
        self.assertEqual(len(shift), 11)
        self.assertEqual(sum(g[0] in ('F', 'Finv') for g in word), 8)
        self.assertEqual(sum(g[0] == 'CX' for g in word),
                         4 * b + 2 * len(shift))
        tables = {'F': [_hrt_f(raw << 1, r) for raw in range(1 << 7)],
                  'Finv': [_hrt_k(raw << 1, r) for raw in range(1 << 7)]}
        inverse = _inverse(word)
        for column in range(1 << width):
            address = column & (m - 1)
            expected_phase = (column >> bank[address]) & 1
            self.assertEqual(_column(word, column, carrier, tables),
                             (column, expected_phase))
            self.assertEqual(_column(inverse, column, carrier, tables),
                             (column, expected_phase))
            # The native CNOT shear itself returns the arbitrary bank and
            # toggles precisely the complementary coefficient coordinates.
            expected_state = column
            for subset in range(1, m):
                bit = sum((column >> bank[u]) & 1 for u in range(m)
                          if u & subset == u) % 2
                expected_state ^= bit << carrier[((1 << r) - 1) ^ subset]
            self.assertEqual(_column(shift, column, carrier, tables),
                             (expected_state, 0))


if __name__ == '__main__':
    unittest.main()
