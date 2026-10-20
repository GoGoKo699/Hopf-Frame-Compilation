"""All-column composition witness with output aliases and quantum helper fibers.

F and its inverse use the retained exact coefficient permutations, not a
new native source emitter. The source's input-dependent helper action is
an explicit native controlled-H/CZ word. It deliberately does not return
the helper; the completed character and its actual reverse must do so.
"""
import unittest

from tests.exact_matrices import CO, CZ, Gate, apply_word
from tests.test_balanced_algebra_permutation import _hrt_f, _hrt_k


def _layout(b):
    m, full = 1 << b, (1 << (b + 1)) - 1
    bank = list(range(b, b + m))
    carrier = {full ^ subset: bank[subset] for subset in range(1, m)}
    extra = b + m
    for subset in range(1, full + 1):
        if subset not in carrier:
            carrier[subset] = extra
            extra += 1
    return bank, carrier, extra


def _fiber_word(carrier, helper):
    # Controlled-H = A† CZ A, A=S H T H S†, in chronological order.
    before = [('SDG', helper), ('H', helper), ('TDG', helper),
              ('H', helper), ('S', helper)]
    after = [('SDG', helper), ('H', helper), ('T', helper),
             ('H', helper), ('S', helper)]
    cz = [('H', helper), ('CX', carrier[1], helper), ('H', helper)]
    phase = [('H', helper), ('CX', carrier[max(carrier)], helper),
             ('H', helper)]
    return before + cz + after + phase


def _apply(word, state, carrier, tables):
    mask = sum(1 << wire for wire in carrier.values())
    for name, *wires in word:
        if name in tables:
            result = {}
            for basis, amplitude in state.items():
                c = sum(((basis >> wire) & 1) << subset
                        for subset, wire in carrier.items())
                value = tables[name][c >> 1]
                image = (basis & ~mask) | sum(
                    ((value >> subset) & 1) << wire
                    for subset, wire in carrier.items())
                result[image] = amplitude
            state = result
        else:
            native = ([Gate('S', tuple(wires))] * 2 if name == 'Z'
                      else [Gate(name, tuple(wires))])
            state = apply_word(native, state)
    return {basis: amplitude for basis, amplitude in state.items()
            if amplitude != CZ}


def _reverse(word):
    inverses = {'F': 'Finv', 'Finv': 'F', 'T': 'TDG', 'TDG': 'T',
                'S': 'SDG', 'SDG': 'S'}
    return [(inverses.get(name, name), *wires)
            for name, *wires in reversed(word)]


class OutputAliasedCharacterTests(unittest.TestCase):
    def test_all_columns_with_nonreturning_quantum_source_fiber(self):
        b, r, m = 2, 3, 4
        bank, carrier, helper = _layout(b)
        width = helper + 1
        self.assertEqual(width, 11)  # b + m outputs + m carrier + one fiber
        self.assertEqual(len(set(carrier.values()) & set(bank)), m - 1)
        self.assertNotIn(bank[0], carrier.values())
        zeta = [('CX', bank[s ^ (1 << bit)], bank[s])
                for bit in range(b) for s in range(m) if s & (1 << bit)]
        source = _fiber_word(carrier, helper) + [('F',)]
        jword = _reverse(source)
        jword += [('CX', bit, carrier[1 << bit]) for bit in range(b)]
        jword += source
        character = jword + [('Z', carrier[(1 << r) - 1])]
        character += _reverse(jword) + [('Z', carrier[(1 << r) - 1])]
        phase = zeta + character + _reverse(zeta) + [('Z', bank[0])]
        self.assertEqual(sum(g[0] in ('F', 'Finv') for g in phase), 4)
        self.assertEqual(2 * len(zeta) + 2 * b + 2 * m + 6, 26)
        tables = {'F': [_hrt_f(raw << 1, r) for raw in range(1 << 7)],
                  'Finv': [_hrt_k(raw << 1, r) for raw in range(1 << 7)]}
        witness = _apply(source, {1 << carrier[1]: CO}, carrier, tables)
        self.assertEqual(len(witness), 2)  # The source leaves coherent work.
        inverse = _reverse(phase)
        for column in range(1 << width):
            address = column & (m - 1)
            sign = -1 if (column >> bank[address]) & 1 else 1
            expected = {column: sign * CO}
            self.assertEqual(_apply(phase, {column: CO}, carrier, tables), expected)
            self.assertEqual(_apply(inverse, {column: CO}, carrier, tables), expected)


if __name__ == '__main__':
    unittest.main()
