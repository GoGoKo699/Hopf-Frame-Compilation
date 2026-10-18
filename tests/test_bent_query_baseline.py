"""Exact paid baseline for the six-address cubic bent XOR query.

Five serialized dirty-helper Toffolis give a complete native word, using
only the seven arbitrary logical wires. The explicit resource ledger is
a finite upper bound, not an optimality claim or a two-layer synthesis.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_dirty_indicator_depth import _exact_phase_action
from tests.test_t_depth import _inverse, _word


def _compact_ccz(a, b, c, dirty):
    """Eight T gates and sixteen CNOTs in two physical T bands."""
    assert len({a, b, c, dirty}) == 4
    first = [('CX', dirty, a), ('CX', dirty, b), ('CX', dirty, c)]
    second = [('CX', a, b), ('CX', b, dirty), ('CX', dirty, c),
              ('CX', c, a), ('CX', a, b)]
    # First parity basis: (a+d,b+d,c+d,d).
    # Second: (b+c+d,a+c+d,a+b+c+d,a+b+d).
    return [
        ('C', first),
        ('T', [('TDG', a), ('TDG', b), ('TDG', c), ('T', dirty)]),
        ('C', list(reversed(first))),
        ('C', second),
        ('T', [('T', a), ('T', b), ('TDG', c), ('T', dirty)]),
        ('C', list(reversed(second))),
    ]


def _compact_toffoli(a, b, target, dirty):
    return ([('C', [('H', target)])]
            + _compact_ccz(a, b, target, dirty)
            + [('C', [('H', target)])])


def _bent_baseline():
    # Logical order: x1,x2,x3,y1,y2,y3,z; no extra physical wires.
    # Chronological A,B1,A,B2,B3. Every helper is an arbitrary address
    # bit and is returned literally before it is used by another gate.
    triples_and_helpers = [
        (1, 2, 3, 4),
        (0, 3, 6, 4),
        (1, 2, 3, 4),
        (1, 4, 6, 3),
        (2, 5, 6, 3),
    ]
    return [stage for args in triples_and_helpers
            for stage in _compact_toffoli(*args)]


def _bent_value(address):
    x, y = address & 7, (address >> 3) & 7
    return ((x & y).bit_count() % 2) ^ int(x == 7)


class BentQueryBaselineTests(unittest.TestCase):
    def test_compact_toffoli_exact_parity_phase_and_dirty_return(self):
        phase_schedule = _compact_ccz(0, 1, 2, 3)
        for schedule in (phase_schedule, _inverse(phase_schedule)):
            for basis in range(16):
                expected_phase = 4 if basis & 7 == 7 else 0
                self.assertEqual(_exact_phase_action(basis, _word(schedule)),
                                 (basis, expected_phase))

        schedule = _compact_toffoli(0, 1, 2, 3)
        actual = native_columns(4, _word(schedule))
        expected = [[CZ for _ in range(16)] for _ in range(16)]
        for column in range(16):
            row = column ^ 4 if column & 3 == 3 else column
            expected[row][column] = CO
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(4, _word(_inverse(schedule))),
                         adjoint(actual))

    def test_full_bent_query_all_columns_literal_phase_and_actual_inverse(self):
        schedule = _bent_baseline()
        actual = native_columns(7, _word(schedule))
        expected = [[CZ for _ in range(128)] for _ in range(128)]
        for column in range(128):
            row = column ^ (_bent_value(column & 63) << 6)
            expected[row][column] = CO
        # Full matrix equality includes both output values, every borrowed
        # address value and off-diagonal coherences with external references.
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(7, _word(_inverse(schedule))),
                         adjoint(expected))

    def test_paid_native_resource_ledger_and_simultaneous_allocation(self):
        primitive = _word(_compact_toffoli(0, 1, 2, 3))
        self.assertEqual(Counter(gate[0] for gate in primitive),
                         Counter({'T': 4, 'TDG': 4, 'CX': 16, 'H': 2}))
        schedule = _bent_baseline()
        for actual in (schedule, _inverse(schedule)):
            word = _word(actual)
            counts = Counter(gate[0] for gate in word)
            self.assertEqual(counts,
                             Counter({'T': 20, 'TDG': 20, 'CX': 80, 'H': 10}))
            bands = [gates for kind, gates in actual if kind == 'T']
            self.assertEqual(len(bands), 10)
            for band in bands:
                self.assertEqual(len(band), 4)
                self.assertEqual(len({gate[1] for gate in band}), 4)
            self.assertEqual({wire for gate in word for wire in gate[1:]},
                             set(range(7)))


if __name__ == '__main__':
    unittest.main()
