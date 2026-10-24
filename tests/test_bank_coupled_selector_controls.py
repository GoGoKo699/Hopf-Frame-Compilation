"""Exact native bank-controlled selector swaps and Hadamards.

The data bank, both selector banks, signal and borrowed helpers are
arbitrary. Small native fixtures retain all literal phases and actual
inverses; larger fixtures check physical T layers and complete counts.
These are paid components, not complete addressed acceptance sources.
"""
from itertools import product
import unittest

from tests.exact_matrices import CO, adjoint, native_columns, zeros
from tests.test_co_moving_selector_readout import _basis_change, _expected_hadamards
from tests.test_fixed_row_transport_readout import _positive_s_schedule
from tests.test_majorana_walsh_query import _inverse as _inverse_word
from tests.test_nonlinear_selector_control import _nonlinear_selector_control
from tests.test_phase_signal_query_limits import _counts
from tests.test_t_depth import _ccz_batch, _fredkin_batch, _inverse, _word


def _parallel_schedules(schedules):
    """Combine corresponding layers of disjoint native gate schedules."""
    kinds = tuple(kind for kind, _ in schedules[0])
    assert all(tuple(kind for kind, _ in schedule) == kinds
               for schedule in schedules)
    return [(kind, [gate for schedule in schedules for gate in schedule[i][1]])
            for i, kind in enumerate(kinds)]


def _bank_swap_schedule(q, controlled=False):
    """Bank Y, selectors s,t, and optionally a signal plus q dirty helpers."""
    assert q >= 1
    offset = int(controlled)
    bank = list(range(offset, offset+q))
    first = list(range(offset+q, offset+2*q))
    second = list(range(offset+2*q, offset+3*q))
    if not controlled:
        schedule = _parallel_schedules([
            _fredkin_batch(y, [(s, t)]) for y, s, t in zip(bank, first, second)])
        return 3*q, _positive_s_schedule(schedule)
    width, phase = _nonlinear_selector_control(q)
    wires = {0: 0}
    for row, (y, s, t) in enumerate(zip(bank, first, second)):
        # Map the component's disjoint triples to (Y,t,s); its last q
        # arbitrary borrowed wires already have the required positions.
        wires.update({1+3*row: y, 2+3*row: t, 3+3*row: s,
                      1+3*q+row: 1+3*q+row})
    phase = [(kind, [(gate[0], *(wires[wire] for wire in gate[1:]))
                     for gate in gates]) for kind, gates in phase]
    outer = [('CX', s, t) for s, t in zip(first, second)]
    walls = [('H', s) for s in first]
    # Only the CCZ middle needs signal control. For signal zero, these
    # unconditional Fredkin Clifford walls cancel exactly.
    schedule = ([('C', outer + walls)] + phase
                + [('C', walls + list(reversed(outer)))])
    return width, _positive_s_schedule(schedule)


def _basis_wall(selectors, inverse=False):
    """One native T band for the disjoint V or actual V-adjoint words."""
    before, band, after = [], [], []
    for selector in selectors:
        word = _basis_change(selector)
        if inverse:
            word = _inverse_word(word)
        positions = [i for i, gate in enumerate(word) if gate[0] in ('T', 'TDG')]
        assert len(positions) == 1
        pivot = positions[0]
        before.extend(word[:pivot])
        band.append(word[pivot])
        after.extend(word[pivot+1:])
    return [('C', before), ('T', band), ('C', after)]


def _bank_hadamard_schedule(q, controlled=False):
    """V Z V^dagger=H, with one distinct bank control for every selector."""
    assert q >= 1
    offset = int(controlled)
    bank = list(range(offset, offset+q))
    selectors = list(range(offset+q, offset+2*q))
    if controlled:
        middle = _ccz_batch(0, list(zip(bank, selectors)))
    else:
        middle = [('C', [gate for y, s in zip(bank, selectors)
                         for gate in [('H', s), ('CX', y, s), ('H', s)]])]
    schedule = (_basis_wall(selectors, inverse=True) + middle
                + _basis_wall(selectors))
    return offset+2*q, _positive_s_schedule(schedule)


def _swap_columns(q, controlled=False, columns=None):
    width, offset = (1+4*q if controlled else 3*q), int(controlled)
    columns = list(range(1 << width)) if columns is None else list(columns)
    expected = zeros(1 << width, len(columns))
    for index, column in enumerate(columns):
        output = column
        if not controlled or column & 1:
            for row in range(q):
                y, s, t = offset+row, offset+q+row, offset+2*q+row
                if (column >> y) & 1 and ((column >> s) ^ (column >> t)) & 1:
                    output ^= (1 << s) | (1 << t)
        expected[output][index] = CO
    return expected


class BankCoupledSelectorControlTests(unittest.TestCase):
    def test_disjoint_bank_swaps_all_one_and_two_row_columns(self):
        for q in (1, 2):
            width, schedule = _bank_swap_schedule(q)
            expected = _swap_columns(q)
            actual = native_columns(width, _word(schedule))
            self.assertEqual(actual, expected)
            reverse = _positive_s_schedule(_inverse(schedule))
            self.assertEqual(native_columns(width, _word(reverse)), adjoint(actual))

    def test_controlled_swap_all_one_row_columns_and_actual_inverse(self):
        width, schedule = _bank_swap_schedule(1, controlled=True)
        expected = _swap_columns(1, controlled=True)
        actual = native_columns(width, _word(schedule))
        self.assertEqual(actual, expected)
        reverse = _positive_s_schedule(_inverse(schedule))
        self.assertEqual(native_columns(width, _word(reverse)), adjoint(actual))

    def test_two_row_controlled_swap_returns_every_borrowed_assignment(self):
        # Every signal, bank and borrowed assignment participates. The selected
        # principal columns include equal, opposite and differently marked rows.
        columns = [signal | (bank << 1) | (first << 3) | (second << 5) | (dirty << 7)
                   for signal, bank, dirty, (first, second) in product(
                       (0, 1), range(4), range(4), ((0, 0), (3, 0), (1, 2), (2, 1)))]
        width, schedule = _bank_swap_schedule(2, controlled=True)
        expected = _swap_columns(2, controlled=True, columns=columns)
        for orientation in (schedule, _positive_s_schedule(_inverse(schedule))):
            self.assertEqual(native_columns(width, _word(orientation), columns), expected)

    def test_bank_hadamards_all_small_columns_and_actual_inverse(self):
        for q, controlled in product((1, 2), (False, True)):
            width, schedule = _bank_hadamard_schedule(q, controlled)
            offset = int(controlled)
            expected = _expected_hadamards(width, range(1 << width), lambda column: [
                (offset+q+row, bool(((column >> (offset+row)) & 1)
                                   and (not controlled or column & 1)))
                for row in range(q)])
            actual = native_columns(width, _word(schedule))
            self.assertEqual(actual, expected)
            reverse = _positive_s_schedule(_inverse(schedule))
            self.assertEqual(native_columns(width, _word(reverse)), adjoint(actual))

    def test_native_layers_and_common_forward_inverse_resource_bounds(self):
        for q in (1, 2, 3, 4, 8, 16, 32, 64):
            cases = (
                (_bank_swap_schedule, False, 4, 7*q, 14*q),
                (_bank_swap_schedule, True, 16, 26*q+2*(q % 2), 48*q+6),
                (_bank_hadamard_schedule, False, 2, 2*q, 15*q),
                (_bank_hadamard_schedule, True, 6, 8*q+q % 2, 22*q+3),
            )
            for emitter, controlled, depth, expected_t, max_clifford in cases:
                width, schedule = emitter(q, controlled)
                for orientation in (schedule, _positive_s_schedule(_inverse(schedule))):
                    self.assertEqual(sum(kind == 'T' for kind, _ in orientation), depth)
                    for kind, gates in orientation:
                        self.assertTrue(all(0 <= wire < width
                                            for gate in gates for wire in gate[1:]))
                        if kind == 'T':
                            self.assertTrue(gates)
                            self.assertTrue(all(gate[0] in ('T', 'TDG') for gate in gates))
                            targets = [gate[1] for gate in gates]
                            self.assertEqual(len(targets), len(set(targets)))
                            if controlled and q % 2 == 0:
                                self.assertNotIn(0, targets)
                        else:
                            self.assertTrue(all(gate[0] in ('H', 'S', 'CX')
                                                for gate in gates))
                    t_count, clifford_count = _counts(_word(orientation))
                    self.assertEqual(t_count, expected_t)
                    self.assertLessEqual(clifford_count, max_clifford)


if __name__ == '__main__':
    unittest.main()
