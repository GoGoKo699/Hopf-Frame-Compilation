"""Exact full-space checks for bank-controlled selector/dirty swaps.

The second row register t is an arbitrary principal input, not initialized
workspace.  It may move within a source or its accepted fiber; the complete
Gram target must retain its identity factor.  Only the stage implementation's
borrowed parity wire is required to return between these components.

Exhaustive two-bit columns and structured larger columns corroborate the
analytic formulas, including literal powers of i and actual inverse order.
They do not establish a complete accepted lookup source or a depth bound.
"""
from itertools import product
import unittest

from tests.test_co_moving_lookup_stage import _Matching


def _bank_swap_literal(m, state):
    """A chronological product of disjoint Fredkins controlled by current Y."""
    x, bank, selectors, work, borrowed, phase = state
    for row in range(m):
        if (bank >> row) & 1:
            difference = ((selectors >> row) ^ (work >> row)) & 1
            selectors ^= difference << row
            work ^= difference << row
    return x, bank, selectors, work, borrowed, phase


def _stage_literal(matching, state, inverse=False):
    x, bank, selectors, work, borrowed, phase = state
    x, bank, selectors, borrowed, phase = matching.stage_literal(
        (x, bank, selectors, borrowed, phase), inverse=inverse)
    return x, bank, selectors, work, borrowed, phase


def _phase_moment(matching, x, bank, register):
    """The integer n_C, independently of the native phase-factor chronology."""
    return sum(sum(((x ^ row) >> control) & 1
                   for control in matching.controls)
               for row in range(matching.m) if (bank & register) >> row & 1)


def _paired_literal(matching, state, inverse=False):
    state = _bank_swap_literal(matching.m, state)
    state = _stage_literal(matching, state, inverse)
    return _bank_swap_literal(matching.m, state)


def _paired_formula(matching, state, inverse=False):
    x, bank, selectors, work, borrowed, phase = state
    active = (bank & work).bit_count() % 2
    parked = ((matching.mask ^ bank) & work) ^ (bank & selectors)
    moment = _phase_moment(matching, x, bank, work)
    if active:
        address = matching.images[x]
        moved_bank = matching.relabel(bank)
        moved_selectors = matching.relabel(selectors)
        moved_work = matching.relabel(work)
    else:
        address, moved_bank = x, bank
        moved_selectors, moved_work = selectors, work
    inactive = matching.mask ^ moved_bank
    output_selectors = (inactive & moved_selectors) ^ (moved_bank & parked)
    output_work = (inactive & parked) ^ (moved_bank & moved_work)
    sign = -1 if inverse else 1
    return (address, moved_bank, output_selectors, output_work, borrowed,
            (phase + sign * moment) % 4)


def _unpaired_literal(first, second, state, inverse=False):
    if inverse:
        first, second = second, first
    state = _stage_literal(first, state, inverse)
    state = _bank_swap_literal(first.m, state)
    return _stage_literal(second, state, inverse)


def _unpaired_formula(first, second, state, inverse=False):
    if inverse:
        first, second = second, first
    x, bank, selectors, work, borrowed, phase = state
    h = (bank & selectors).bit_count() % 2
    moment = _phase_moment(first, x, bank, selectors)
    if h:
        address, moved_bank = first.images[x], first.relabel(bank)
        moved_selectors = first.relabel(selectors)
    else:
        address, moved_bank, moved_selectors = x, bank, selectors
    k = (moved_bank & work).bit_count() % 2
    moment += _phase_moment(second, address, moved_bank, work)
    parked = ((first.mask ^ moved_bank) & work) ^ (moved_bank & moved_selectors)
    if k:
        address = second.images[address]
        moved_bank = second.relabel(moved_bank)
        moved_selectors = second.relabel(moved_selectors)
        work = second.relabel(work)
    output_selectors = ((first.mask ^ moved_bank) & moved_selectors) ^ (moved_bank & work)
    sign = -1 if inverse else 1
    return (address, moved_bank, output_selectors, parked, borrowed,
            (phase + sign * moment) % 4)


class BankCoupledSelectorDynamicsTests(unittest.TestCase):
    def test_paired_every_two_bit_column_and_actual_inverse(self):
        for pairs in (((0, 1),), ((1, 0),)):
            matching = _Matching(2, pairs)
            fixed_rows = sum(1 << row for row, image in enumerate(matching.images)
                             if row == image)
            for x, bank, selectors, work, borrowed in product(
                    range(4), range(16), range(16), range(16), (0, 1)):
                phase = (x + bank + selectors + work) % 4
                state = x, bank, selectors, work, borrowed, phase
                actual = _paired_literal(matching, state)
                self.assertEqual(actual, _paired_formula(matching, state))
                self.assertEqual(_paired_literal(matching, actual, inverse=True), state)
                self.assertEqual(_paired_literal(matching, state, inverse=True),
                                 _paired_formula(matching, state, inverse=True))
                self.assertEqual(actual[4], borrowed)
                self.assertEqual((actual[3] ^ work) & fixed_rows, 0)
                # The paired permutation leaves its trigger and integer phase
                # invariant; these are needed for the same-permutation adjoint.
                self.assertEqual((actual[1] & actual[3]).bit_count() % 2,
                                 (bank & work).bit_count() % 2)
                self.assertEqual(_phase_moment(matching, actual[0], actual[1], actual[3]),
                                 _phase_moment(matching, x, bank, work))

    def test_unpaired_every_two_bit_column_and_reversed_stage_order(self):
        first, second = _Matching(2, ((0, 1),)), _Matching(2, ((1, 0),))
        for first, second in ((first, second), (second, first)):
            for x, bank, selectors, work, borrowed in product(
                    range(4), range(16), range(16), range(16), (0, 1)):
                phase = (x + bank + selectors + work) % 4
                state = x, bank, selectors, work, borrowed, phase
                actual = _unpaired_literal(first, second, state)
                self.assertEqual(actual, _unpaired_formula(first, second, state))
                self.assertEqual(_unpaired_literal(first, second, actual, inverse=True), state)
                self.assertEqual(_unpaired_literal(first, second, state, inverse=True),
                                 _unpaired_formula(first, second, state, inverse=True))
                self.assertEqual(actual[4], borrowed)

    def test_structured_larger_banks_retain_literal_phases_and_dirty_motion(self):
        families = ((3, ((0, 1),), ((1, 2),)),
                    (3, ((2, 0),), ((0, 1),)),
                    (4, ((0, 1), (2, 3)), ((0, 2), (1, 3))))
        phases_seen, moved_columns = set(), 0
        for p, first_pairs, second_pairs in families:
            first, second = _Matching(p, first_pairs), _Matching(p, second_pairs)
            m, full = first.m, first.mask
            even = sum(1 << row for row in range(m) if row % 2 == 0)
            banks = (0, full, 1, 2, 3, 1 << (m - 1), even, full ^ even)
            registers = (0, full, 1, 2, even, full ^ even)
            for x, bank, selectors, work, borrowed, phase in product(
                    (0, 1, m // 2, m - 1), banks, registers, registers, (0, 1), (0, 3)):
                state = x, bank, selectors, work, borrowed, phase
                for literal, formula, matchings in (
                        (_paired_literal, _paired_formula, (first,)),
                        (_unpaired_literal, _unpaired_formula, (first, second))):
                    actual = literal(*matchings, state)
                    self.assertEqual(actual, formula(*matchings, state))
                    self.assertEqual(literal(*matchings, actual, inverse=True), state)
                    self.assertEqual(literal(*matchings, state, inverse=True),
                                     formula(*matchings, state, inverse=True))
                    phases_seen.add((actual[-1] - phase) % 4)
                    moved_columns += actual[3] != work
        self.assertEqual(phases_seen, {0, 1, 2, 3})
        self.assertGreater(moved_columns, 0)

    def test_nonreturn_phase_and_actual_inverse_are_distinct_obligations(self):
        first, second = _Matching(2, ((0, 1),)), _Matching(2, ((1, 0),))
        # Even trigger gives identity permutation but a nontrivial i phase.
        state = (0, 3, 5, 3, 1, 0)
        self.assertEqual(_paired_literal(first, state), (0, 3, 5, 3, 1, 1))
        self.assertEqual(_paired_literal(first, state, inverse=True), (0, 3, 5, 3, 1, 3))
        # A moved occupied row transports t; t is not an implementation helper.
        state = (0, 2, 0, 2, 1, 0)
        actual = _paired_literal(first, state)
        self.assertEqual(actual, (0, 8, 0, 8, 1, 1))
        self.assertEqual(_paired_literal(first, actual, inverse=True), state)
        # Merely negating phases without reversing the stage order is not U†.
        state = (0, 2, 2, 0, 1, 0)
        actual = _unpaired_literal(first, second, state)
        self.assertEqual(_unpaired_literal(first, second, actual, inverse=True), state)
        wrong = _stage_literal(first, actual, inverse=True)
        wrong = _bank_swap_literal(first.m, wrong)
        wrong = _stage_literal(second, wrong, inverse=True)
        self.assertNotEqual(wrong, state)

    def test_origin_joint_zero_code_and_empty_bank_are_pointwise_fixed(self):
        for p, pairs in ((2, ((0, 1),)), (2, ((1, 0),)),
                         (3, ((0, 1),)), (3, ((1, 2),)), (3, ((2, 0),))):
            matching = _Matching(p, pairs)
            registers = range(1 << matching.m)
            for selectors, work, borrowed in product(
                    range(0, 1 << matching.m, 2), range(0, 1 << matching.m, 2), (0, 1)):
                state = (0, 1, selectors, work, borrowed, 3)
                self.assertEqual(_bank_swap_literal(matching.m, state), state)
                self.assertEqual(_stage_literal(matching, state), state)
                self.assertEqual(_stage_literal(matching, state, inverse=True), state)
            # Empty-bank identity holds at every p=2 address and row-register
            # column; larger empty cases include sparse and dense registers.
            selected = (registers if p == 2 else
                        (0, 1, 2, matching.mask, matching.mask ^ 1, 0x55, 0xAA))
            for x, selectors, work, borrowed in product(
                    range(matching.m), selected, selected, (0, 1)):
                state = (x, 0, selectors, work, borrowed, 1)
                self.assertEqual(_bank_swap_literal(matching.m, state), state)
                self.assertEqual(_stage_literal(matching, state), state)
                self.assertEqual(_stage_literal(matching, state, inverse=True), state)


if __name__ == '__main__':
    unittest.main()
