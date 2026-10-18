"""Exact native fixtures for the read-only three/four-layer frontier.

The tests use the actual CNOT basis words, their reverses, physical T
layers and endpoint Hadamards. Complete matrices retain literal phases
on every work input; initialized-sector correctness is not substituted
for full return of a borrowed helper.
"""
from itertools import product
import unittest

from tests.exact_matrices import CI, CO, CZ, native_columns


def _phase_layer(compute, phases):
    # Every compute gate here is a CNOT and hence is its own inverse.
    return [('C', compute), ('T', phases), ('C', list(reversed(compute)))]


def _zero_clean_four_layers(a=0, b=1, t=2, d=3):
    schedule = [('C', [('H', t)])]
    for s1, s2 in product((0, 1), repeat=2):
        compute = ([('CX', a, d)] if s1 else [])
        compute += ([('CX', b, d)] if s2 else [])
        compute += [('CX', d, t)]
        positive = (s1 + s2) % 2 == 0
        phases = [('T' if positive else 'TDG', d),
                  ('TDG' if positive else 'T', t)]
        schedule += _phase_layer(compute, phases)
    return schedule + [('C', [('H', t)])]


def _one_clean_three_layers(a=0, b=1, t=2, d=3, c=4):
    # The active parities correspond respectively to
    # s=000,001,100; 010,101,110; 011,111 in the eight-term cube.
    first = [('CX', d, t), ('CX', d, c), ('CX', a, c)]
    second = [('CX', a, c), ('CX', c, t), ('CX', d, t),
              ('CX', d, c), ('CX', b, d), ('CX', b, c)]
    third = [('CX', t, d), ('CX', b, d), ('CX', d, c),
             ('CX', a, c)]
    schedule = [('C', [('H', t)])]
    schedule += _phase_layer(first, [('T', d), ('TDG', t), ('TDG', c)])
    schedule += _phase_layer(second, [('TDG', d), ('T', t), ('T', c)])
    schedule += _phase_layer(third, [('T', d), ('TDG', c)])
    return schedule + [('C', [('H', t)])]


def _flatten(schedule):
    return [gate for _, gates in schedule for gate in gates]


def _one_clean_relative_two_layers(a=0, b=1, t=2, c=3):
    before = [('H', t), ('CX', t, c)]
    after = list(reversed(before))
    cz = [('H', t), ('CX', b, t), ('H', t)]
    return [
        ('C', [('CX', a, c)] + before),
        ('T', [('T', t), ('TDG', c)]),  # actual D inverse
        ('C', after + cz + before),
        ('T', [('TDG', t), ('T', c)]),
        ('C', after + cz + [('CX', a, c)]),
    ]


class ReadOnlyQueryLayerFixtureTests(unittest.TestCase):
    def test_zero_clean_four_layers_are_literal_toffoli_on_every_input(self):
        actual = native_columns(4, _flatten(_zero_clean_four_layers()))
        for column in range(16):
            row = column ^ (4 if (column & 1) and (column & 2) else 0)
            self.assertEqual([actual[i][column] for i in range(16)],
                             [CO if i == row else CZ for i in range(16)])

    def test_one_clean_three_layers_have_the_claimed_full_space_extension(self):
        actual = native_columns(5, _flatten(_one_clean_three_layers()))
        for column in range(32):
            a = column & 1
            b = (column >> 1) & 1
            c = (column >> 4) & 1
            row = column ^ (4 if b and (a ^ c) else 0)
            self.assertEqual([actual[i][column] for i in range(32)],
                             [CO if i == row else CZ for i in range(32)])

    def test_native_layer_and_elementary_gate_ledgers(self):
        for schedule, depth, cliffords, width in (
                (_zero_clean_four_layers(), 4, 18, 4),
                (_one_clean_three_layers(), 3, 28, 5)):
            layers = [gates for kind, gates in schedule if kind == 'T']
            self.assertEqual(len(layers), depth)
            self.assertEqual(sum(map(len, layers)), 8)
            self.assertEqual(sum(len(gates) for kind, gates in schedule
                                 if kind == 'C'), cliffords)
            for layer in layers:
                targets = [gate[1] for gate in layer]
                self.assertEqual(len(targets), len(set(targets)))
                self.assertTrue(all(2 <= target < width for target in targets))
            for gate in _flatten(schedule):
                if gate[0] == 'CX':
                    self.assertGreaterEqual(gate[2], 2)
                elif gate[0] == 'H':
                    self.assertEqual(gate[1], 2)

    def test_two_layer_candidate_retains_its_nonclifford_address_phase(self):
        schedule = _one_clean_relative_two_layers()
        actual = native_columns(4, _flatten(schedule))
        for column in range(16):
            a = column & 1
            b = (column >> 1) & 1
            c = (column >> 3) & 1
            active = b * (a ^ c)
            row = column ^ (4 if active else 0)
            phase = -CI if active else CO
            self.assertEqual([actual[i][column] for i in range(16)],
                             [phase if i == row else CZ for i in range(16)])
        layers = [gates for kind, gates in schedule if kind == 'T']
        self.assertEqual(len(layers), 2)
        self.assertEqual(sum(map(len, layers)), 4)
        self.assertEqual(sum(len(gates) for kind, gates in schedule
                             if kind == 'C'), 16)
        # On clean c=0, address a=b=1 witnesses literal failure, not leakage.
        self.assertEqual(actual[7][3], -CI)
        self.assertNotEqual(actual[7][3], CO)


if __name__ == '__main__':
    unittest.main()
