"""Bounded exact audits of consume-before-change window predicates.

Small block predicates and flag words are literal X/CX/CCX circuits.
Completed group bodies are reduced rational unitaries retaining an arbitrary
source bit; they are not native group, predicate-depth, or QROM emitters.
Two arbitrary dirty helpers are explicit. Active cache inputs are zero;
inactive H=h=0 inputs have arbitrary cache, source, logical, and dirty bits.
"""
from __future__ import annotations

from fractions import Fraction
import unittest

import numpy as np

try:
    from .test_grouped_program_prefetch import _basis_word
    from .test_grouped_selector_reuse import _add, _exact_word
    from .test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis,
    )
except ImportError:
    from test_grouped_program_prefetch import _basis_word
    from test_grouped_selector_reuse import _add, _exact_word
    from test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis,
    )


def _read(basis, wires):
    return sum(((basis >> wire) & 1) << j for j, wire in enumerate(wires))


def _zero_word(controls, target):
    """Only one- and two-bit blocks occur in this finite fixture."""
    assert 1 <= len(controls) <= 2
    negative = [('X', q) for q in controls]
    gate = ('CX' if len(controls) == 1 else 'CCX', *controls, target)
    return negative + [gate] + _adjoint(negative)


def _flag_word(H, v, u, h, dirty):
    left = ('CCX', H, v, dirty)
    right = ('CCX', dirty, u, h)
    return [left, right, left, right]


class _Window:
    def __init__(self, sizes):
        self.sizes = sizes
        self.H, self.h = 0, 1
        self.logical = list(range(2, 2 + sum(sizes)))
        self.source, tail = 2 + sum(sizes), 3 + sum(sizes)
        self.tail = [tail]
        start = tail + 1
        self.e = list(range(start, start + len(sizes)))
        start += len(sizes)
        self.u = list(range(start, start + len(sizes)))
        self.v = start + len(sizes)
        self.dirty = [self.v + 1, self.v + 2]
        self.width = self.v + 3
        self.blocks, start = [], 0
        for size in sizes:
            self.blocks.append(self.logical[start:start + size])
            start += size

    def chain(self, i):
        if i == len(self.sizes) - 1:
            return [('X', self.u[i])]
        return [('CCX', self.e[i + 1], self.u[i + 1], self.u[i])]

    def block_zero(self, i):
        return _zero_word(self.blocks[i], self.e[i])

    def flag(self, i, missing_H=False):
        if missing_H:
            return [('CCX', self.v, self.u[i], self.h)]
        return _flag_word(self.H, self.v, self.u[i], self.h, self.dirty[0])

    def prepare(self):
        return (_zero_word(self.tail, self.v)
                + [gate for i in range(len(self.sizes)) for gate in self.block_zero(i)]
                + [gate for i in reversed(range(len(self.sizes))) for gate in self.chain(i)])

    def stage(self, state, position, end, direct=False):
        """Row-dependent rational rotation followed by retained source leakage."""
        target = self.logical[position]
        result = {}
        for basis, amplitude in state.items():
            if direct:
                enabled = (bool((basis >> self.H) & 1)
                           and _read(basis, self.tail + self.logical[position + 1:]) == 0)
            else:
                enabled = (bool((basis >> self.h) & 1)
                           and _read(basis, self.logical[position + 1:end]) == 0)
            if not enabled:
                _add(result, basis, amplitude)
                continue
            row = _read(basis, self.logical[:position])
            c, s = ((Fraction(3, 5), Fraction(4, 5)) if (row + position) % 2 == 0
                    else (Fraction(5, 13), Fraction(12, 13)))
            bit = (basis >> target) & 1
            for image, factor in ((basis, c), (basis ^ (1 << target), s * (-1) ** bit)):
                if (image >> target) & 1:
                    image ^= 1 << self.source
                _add(result, image, amplitude * factor)
        return result

    def body(self, state, i, transient=True):
        start, end = sum(self.sizes[:i]), sum(self.sizes[:i + 1])
        # A future logical wire is temporarily borrowed and restored inside
        # the completed body. No cached predicate is retested mid-body.
        future = self.logical[end] if end < len(self.logical) else self.tail[0]
        scratch = [('CX', self.h, future)] if transient else []
        state = _exact_word(state, scratch)
        for position in range(start, end):
            state = self.stage(state, position, end)
        return _exact_word(state, _adjoint(scratch))

    def run(self, state, *, late_e=False, late_chain=False, missing_H=False):
        state = _exact_word(state, self.prepare())
        for i in range(len(self.sizes)):
            if not late_e:
                state = _exact_word(state, _adjoint(self.block_zero(i)))
            flag = self.flag(i, missing_H)
            state = _exact_word(state, flag)
            state = self.body(state, i)
            state = _exact_word(state, _adjoint(flag))
            if late_e:
                state = _exact_word(state, _adjoint(self.block_zero(i)))
            if not late_chain:
                state = _exact_word(state, _adjoint(self.chain(i)))
        if late_chain:
            for i in range(len(self.sizes)):
                state = _exact_word(state, _adjoint(self.chain(i)))
        return _exact_word(state, _adjoint(_zero_word(self.tail, self.v)))

    def direct(self, state):
        # Independent full logical suffix predicates: no e/u/v cache or h.
        for position in range(len(self.logical)):
            state = self.stage(state, position, len(self.logical), direct=True)
        return state

    def inputs(self, H=1, include_cache=False):
        wires = self.logical + [self.source] + self.tail + self.dirty
        if include_cache:
            wires += self.e + self.u + [self.v]
        return [H | sum(((data >> j) & 1) << q for j, q in enumerate(wires))
                for data in range(1 << len(wires))]


class WindowedGroupPredicateTests(unittest.TestCase):
    def test_exact_active_unequal_windows_on_every_data_and_source_column(self):
        for sizes in ((1,), (2,), (1, 2), (2, 1), (1, 2, 1)):
            f = _Window(sizes)
            cache = sum(1 << q for q in f.e + f.u + [f.v, f.h])
            for basis in f.inputs():
                initial = {basis: Fraction(1)}
                actual = f.run(initial)
                self.assertEqual(actual, f.direct(initial))
                self.assertTrue(all(not (image & cache) for image in actual))
                self.assertEqual(sum(value * value for value in actual.values()), 1)
        # Exact equality on all columns extends to coherent reference inputs.

    def test_reduced_bodies_retain_source_and_restore_transient_future_work(self):
        f = _Window((1, 2))
        initial = {1: Fraction(1)}
        first = f.stage(initial, 0, len(f.logical), direct=True)
        self.assertEqual(len(first), 2)
        self.assertTrue(any((basis >> f.source) & 1 for basis in first))
        self.assertNotEqual(
            f.stage(first, 1, len(f.logical), direct=True),
            f.stage(f.stage(initial, 1, len(f.logical), direct=True),
                    0, len(f.logical), direct=True),
        )
        enabled = {1 | (1 << f.h): Fraction(1)}
        borrowed = _exact_word(enabled, [('CX', f.h, f.blocks[1][0])])
        self.assertTrue(all((basis >> f.blocks[1][0]) & 1 for basis in borrowed))
        self.assertEqual(f.body(enabled, 0), f.body(enabled, 0, transient=False))

    def test_inactive_identity_for_every_arbitrary_cache_and_dirty_input(self):
        for sizes in ((1,), (1, 2), (1, 2, 1)):
            f = _Window(sizes)
            initial = f.inputs(H=0, include_cache=True)
            labels, phases = _basis_word(f.prepare(), initial)
            for i in range(len(sizes)):
                labels, signs = _basis_word(_adjoint(f.block_zero(i)) + f.flag(i), labels)
                phases *= signs
                # Each complete reduced body is literally identity here.
                self.assertTrue(all(not ((int(basis) >> f.h) & 1) for basis in labels))
                labels, signs = _basis_word(_adjoint(f.flag(i)) + _adjoint(f.chain(i)), labels)
                phases *= signs
            labels, signs = _basis_word(_adjoint(_zero_word(f.tail, f.v)), labels)
            np.testing.assert_array_equal(labels, initial)
            np.testing.assert_array_equal(phases * signs, 1)

    def test_native_constant_arity_flag_has_literal_phase_and_dirty_return(self):
        # H,v,u,h,d plus the second untouched returned helper.
        word = _flag_word(0, 1, 2, 3, 4)
        native = _expand_toffolis(word)
        self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in native), 28)
        initial = np.eye(64, dtype=complex)
        expected = np.zeros_like(initial)
        for basis in range(64):
            image = basis ^ ((((basis >> 0) & 1) & ((basis >> 1) & 1)
                              & ((basis >> 2) & 1)) << 3)
            expected[image, basis] = 1
        actual = _apply_native_word(6, native, initial)
        np.testing.assert_allclose(actual, expected, atol=3e-12, rtol=0)
        np.testing.assert_allclose(_apply_native_word(6, _adjoint(native), actual),
                                   initial, atol=3e-12, rtol=0)

    def test_late_cache_erasure_and_missing_original_H_are_detected(self):
        f = _Window((1, 2, 1))
        initial = {1: Fraction(1)}
        expected = f.direct(initial)
        self.assertNotEqual(f.run(initial, late_e=True), expected)
        self.assertNotEqual(f.run(initial, late_chain=True), expected)
        inactive = {0: Fraction(1)}
        self.assertEqual(f.run(inactive), inactive)
        self.assertNotEqual(f.run(inactive, missing_H=True), inactive)


if __name__ == '__main__':
    unittest.main()
