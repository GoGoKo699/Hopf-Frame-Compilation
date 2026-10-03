"""Exact bounded audits of cached prefix selectors and consumed suffix enables.

Predicate words are literal X/CX/CCX circuits. The active fixture uses exact
rational rotations and retains an entangled source bit across noncommuting
stages; it is a selector-interface test, not a native amplified frame emitter.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

try:
    from .test_grouped_program_prefetch import _basis_word
    from .test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis,
    )
except ImportError:
    from test_grouped_program_prefetch import _basis_word
    from test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis,
    )


def _prefix_step(parent, target_bit, children, copies):
    """Private controls make the simultaneous Toffolis genuinely disjoint."""
    copy_word, toffolis = [], []
    for row, parent_wire in enumerate(parent):
        for value in (0, 1):
            child = 2 * row + value
            a, b = copies[2 * child:2 * child + 2]
            copy_word += [('CX', parent_wire, a), ('CX', target_bit, b)]
            if value == 0:
                copy_word.append(('X', b))
            toffolis.append(('CCX', a, b, children[child]))
    return copy_word + toffolis + _adjoint(copy_word), toffolis


class _CachedPredicates:
    def __init__(self, g):
        self.g = g
        self.h = 0
        self.local = list(range(1, 1 + g))
        self.source = 1 + g
        start = self.source + 1
        self.u = list(range(start, start + g))
        start += g
        self.prefix = []
        for ell in range(g):
            self.prefix.append(list(range(start, start + (1 << ell))))
            start += 1 << ell
        self.copies = list(range(start, start + (1 << g)))
        self.width = start + (1 << g)
        self.steps, self.rounds = [], []
        for ell in range(g - 1):
            step, rounds = _prefix_step(
                self.prefix[ell], self.local[ell],
                self.prefix[ell + 1], self.copies,
            )
            self.steps.append(step)
            self.rounds.append(rounds)

    def suffix_gate(self, ell):
        if ell == self.g - 1:
            return [('X', self.u[ell])]
        return [('X', self.local[ell + 1]),
                ('CCX', self.local[ell + 1], self.u[ell + 1], self.u[ell]),
                ('X', self.local[ell + 1])]

    def prepare(self):
        return ([gate for ell in reversed(range(self.g))
                 for gate in self.suffix_gate(ell)]
                + [('X', self.prefix[0][0])])

    def finish_prefix(self, wrong_order=False):
        indices = range(self.g - 1) if wrong_order else reversed(range(self.g - 1))
        return ([gate for ell in indices for gate in _adjoint(self.steps[ell])]
                + [('X', self.prefix[0][0])])

    def inactive_word(self):
        word = self.prepare()
        for ell in range(self.g):
            # The intervening h-guarded stage is exactly I for h=0.
            word += self.suffix_gate(ell)
            if ell < self.g - 1:
                word += self.steps[ell]
        return word + self.finish_prefix()


def _add(state, basis, coefficient):
    state[basis] = state.get(basis, Fraction(0)) + coefficient
    if not state[basis]:
        del state[basis]


def _exact_word(state, word):
    labels = list(state)
    outputs, signs = _basis_word(word, labels)
    return {int(output): int(sign) * state[label]
            for label, output, sign in zip(labels, outputs, signs)}


def _selected(basis, fixture, ell, row, cached):
    if not ((basis >> fixture.h) & 1):
        return False
    if cached:
        return bool(((basis >> fixture.u[ell]) & 1)
                    and ((basis >> fixture.prefix[ell][row]) & 1))
    later_zero = all(not ((basis >> q) & 1) for q in fixture.local[ell + 1:])
    # Tree children are ordered by append, with the earliest bit most significant.
    prefix = 0
    for q in fixture.local[:ell]:
        prefix = 2 * prefix + ((basis >> q) & 1)
    return later_zero and prefix == row


def _exact_stage(state, fixture, ell, cached):
    """A literal controlled rational rotation followed by source leakage.

    Each selected row rotates t_ell, then toggles the source conditional on
    its new value. No source reset is performed. All other logical bits
    and every cached predicate are preserved by the completed stage.
    """
    target = fixture.local[ell]
    for row in range(1 << ell):
        c, s = ((Fraction(3, 5), Fraction(4, 5)) if (row + ell) % 2 == 0
                else (Fraction(5, 13), Fraction(12, 13)))
        result = {}
        for basis, amplitude in state.items():
            if not _selected(basis, fixture, ell, row, cached):
                _add(result, basis, amplitude)
                continue
            bit = (basis >> target) & 1
            for image, factor in ((basis, c), (basis ^ (1 << target), s * (-1) ** bit)):
                if (image >> target) & 1:
                    image ^= 1 << fixture.source
                _add(result, image, amplitude * factor)
        state = result
    return state


def _cached_group(state, fixture, delayed_suffix=False, wrong_prefix_order=False):
    state = _exact_word(state, fixture.prepare())
    for ell in range(fixture.g):
        state = _exact_stage(state, fixture, ell, cached=True)
        if not delayed_suffix:
            state = _exact_word(state, fixture.suffix_gate(ell))
        if ell < fixture.g - 1:
            state = _exact_word(state, fixture.steps[ell])
    if delayed_suffix:
        for ell in range(fixture.g):
            state = _exact_word(state, fixture.suffix_gate(ell))
    return _exact_word(state, fixture.finish_prefix(wrong_prefix_order))


def _direct_group(state, fixture):
    for ell in range(fixture.g):
        state = _exact_stage(state, fixture, ell, cached=False)
    return state


class GroupedSelectorReuseTests(unittest.TestCase):
    def test_active_two_and_three_stage_operator_with_arbitrary_source(self):
        for g in (1, 2, 3):
            fixture = _CachedPredicates(g)
            # All local/source input basis vectors: equality extends to arbitrary
            # input superpositions and references, with exact Fraction arithmetic.
            for data in range(1 << (g + 1)):
                state = {(data << 1) | 1: Fraction(1)}
                actual = _cached_group(state, fixture)
                self.assertEqual(actual, _direct_group(state, fixture))
                self.assertTrue(all(basis >> (fixture.source + 1) == 0 for basis in actual))
                self.assertEqual(sum(value * value for value in actual.values()), 1)

    def test_fixture_leaks_source_and_uses_noncommuting_stages(self):
        fixture = _CachedPredicates(2)
        state = {1: Fraction(1)}
        first = _exact_stage(state, fixture, 0, cached=False)
        self.assertTrue(any((basis >> fixture.source) & 1 for basis in first))
        # The two target/source amplitudes certify an entangled, unreset source.
        self.assertEqual(len(first), 2)
        self.assertNotEqual(
            _exact_stage(first, fixture, 1, cached=False),
            _exact_stage(_exact_stage(state, fixture, 1, cached=False),
                         fixture, 0, cached=False),
        )

    def test_inactive_identity_on_all_arbitrary_work_inputs(self):
        fixture = _CachedPredicates(2)
        inputs = range(0, 1 << fixture.width, 2)
        images, phases = _basis_word(fixture.inactive_word(), inputs)
        np.testing.assert_array_equal(images, list(inputs))
        np.testing.assert_array_equal(phases, 1)

    def test_cleanup_order_negative_controls(self):
        fixture = _CachedPredicates(3)
        state = {1: Fraction(1)}
        expected = _direct_group(state, fixture)
        self.assertNotEqual(_cached_group(state, fixture, delayed_suffix=True), expected)
        self.assertNotEqual(_cached_group(state, fixture, wrong_prefix_order=True), expected)

    def test_literal_private_copies_native_phase_and_parallel_support(self):
        word, round_gates = _prefix_step([0], 1, [2, 3], [4, 5, 6, 7])
        support = [q for gate in round_gates for q in gate[1:]]
        self.assertEqual(len(support), len(set(support)))
        native = _expand_toffolis(word)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in native), 14)
        initial = np.zeros((256, 4), dtype=complex)
        initial[np.arange(4), np.arange(4)] = 1
        images, phases = _basis_word(word, range(4))
        expected = np.zeros_like(initial)
        expected[np.asarray(images, dtype=int), np.arange(4)] = phases
        actual = _apply_native_word(8, native, initial)
        np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=0)
        np.testing.assert_allclose(_apply_native_word(8, _adjoint(native), actual),
                                   initial, atol=2e-12, rtol=0)
        for basis in images:
            self.assertEqual(int(basis) >> 4, 0)

    def test_exact_count_and_reservation_ledger(self):
        for m in range(2, 30):
            for g in range(1, m + 1):
                reservation = 11 * m * (1 << g) + 5 * m + (1 << (g + 1)) + g - 1
                self.assertLessEqual(reservation, 15 * m * (1 << g))
                self.assertLessEqual(reservation, 16 * m * (1 << g))
        for g in range(1, 8):
            fixture = _CachedPredicates(g)
            toffolis = sum(gate[0] == 'CCX' for gate in fixture.inactive_word())
            self.assertEqual(toffolis, 2 * ((1 << g) + g - 3))


if __name__ == '__main__':
    unittest.main()
