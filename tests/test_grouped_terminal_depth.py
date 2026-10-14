"""Bounded native audits of whole-word coarse prefetch in terminal groups.

The fixed-alphabet interpreter uses literal Clifford+T words, the existing
dirty multiple-control emitter, and a one-bit-address XOR query. Only a
small primitive has a square matrix; whole-word tests propagate columns.
These are not a shallow table emitter or an asymptotic depth proof. The
integer schedule uses illustrative constants and ceiling logarithms; it
checks rounding and balance, not a universal numerical crossover.
"""
from __future__ import annotations

import math
import unittest

import numpy as np

from compiler_robust_hopf.native_ordered_banks import CompactEmitter

try:
    from .test_operator_source_compiler import (
        H, _adjoint, _apply_native_word, _expand_toffolis, _word_matrix,
    )
except ImportError:
    from test_operator_source_compiler import (
        H, _adjoint, _apply_native_word, _expand_toffolis, _word_matrix,
    )


ATOL = 2e-10
ALPHABET = (np.eye(2), H, _word_matrix(1, [('T', 0)]),
            _word_matrix(1, [('TDG', 0)]))


def _mcx(controls, target, borrow):
    """Reuse the exact dirty-control emitter; no allocated bank is touched."""
    width = max([*controls, target, borrow]) + 1
    emitter = CompactEmitter(width)
    emitter.mcx(controls, target, borrow)
    word = [(name.upper(), *wires) for name, *wires in emitter.gates]
    assert all(max(gate[1:]) < width for gate in word)
    return word


def _controlled_h(controls, target, borrow):
    # V^dagger, C_F(Z), V, with V Z V^dagger=H, including literal phase.
    before = [('SDG', target), ('H', target), ('TDG', target),
              ('H', target), ('S', target)]
    central = [('H', target)] + _mcx(controls, target, borrow) + [('H', target)]
    return before + central + _adjoint(before)


def _phase_s(controls, dirty, borrow):
    toggle = _mcx(controls, dirty, borrow)
    return (toggle + [('S', dirty)] + _adjoint(toggle) + [('SDG', dirty)]
            + [('H', dirty)] + toggle + [('H', dirty)])


def _phase_t(controls, z, w, borrow, *, omit_correction=False):
    toggle = _mcx(controls, z, borrow)
    word = toggle + [('T', z)] + _adjoint(toggle) + [('TDG', z)]
    if not omit_correction:
        word += _phase_s([*controls, z], w, borrow)
    return word


def _letter(target, activity, predicate, code, dirty, *, missing_h=False,
            omit_phase_correction=False):
    """Two-bit codes 0:I, 1:H, 2:T, 3:T^dagger; all work is returned."""
    z, w, borrow = dirty
    common = ([activity] if not missing_h else []) + [predicate]
    low, high = code
    # Only code 01 enables H. The negative high-bit control is restored.
    word = ([('X', high)] + _controlled_h([*common, low, high], target, borrow)
            + [('X', high)])
    # For high=1, T followed by low-controlled S^dagger gives T or T^dagger.
    word += _phase_t([*common, high, target], z, w, borrow,
                     omit_correction=omit_phase_correction)
    word += _adjoint(_phase_s([*common, high, low, target], z, borrow))
    return word


def _load_word(program, activity, address, tables):
    word = []
    for row, packed in enumerate(tables):
        negative = [('X', address)] if row == 0 else []
        word += negative
        word += [('CCX', activity, address, wire) for bit, wire in enumerate(program)
                 if (packed >> bit) & 1]
        word += _adjoint(negative)
    return _expand_toffolis(word)


def _whole_word(*, missing_h=False, mutable_address=False):
    # One target, prefix, local suffix, external H, and reusable predicate.
    target, prefix, suffix, activity, predicate = range(5)
    program, dirty = list(range(5, 9)), [9, 10, 11]
    programs = ((1, 2), (3, 1))  # H,T and T^dagger,H, in chronological order.
    tables = [sum(code << (2 * slot) for slot, code in enumerate(row))
              for row in programs]
    address = target if mutable_address else prefix
    query = _load_word(program, activity, address, tables)
    condition = [('X', suffix), ('CX', suffix, predicate), ('X', suffix)]
    interpreter = []
    for slot in range(2):
        interpreter += _letter(target, activity, predicate, program[2 * slot:2 * slot + 2],
                               dirty, missing_h=missing_h)
    word = condition + query + interpreter + _adjoint(query) + _adjoint(condition)
    return dict(word=word, programs=programs, width=12, target=target,
                prefix=prefix, suffix=suffix, activity=activity,
                predicate=predicate, program=program, dirty=dirty)


def _active_columns(fixture):
    """Every target/prefix/suffix input, entangled with a returned dirty bank."""
    result = np.zeros((1 << fixture['width'], 8), dtype=complex)
    for column in range(8):
        logical = column | (1 << fixture['activity'])
        for dirty in range(8):
            # All dirty amplitudes are nonzero; phases also depend on the
            # logical column, so this includes coherent dirty references.
            result[logical | (dirty << fixture['dirty'][0]), column] = (
                np.exp(1j * (dirty + 1) * (column + 1) / 7) / np.sqrt(8))
    return result


def _ideal_whole_word(fixture, columns):
    result = columns.copy()
    matrices = []
    for row in fixture['programs']:
        matrix = np.eye(2, dtype=complex)
        for code in row:
            matrix = ALPHABET[code] @ matrix
        matrices.append(matrix)
    for low in range(0, len(result), 2):
        if ((low >> fixture['activity']) & 1) and not ((low >> fixture['suffix']) & 1):
            address = (low >> fixture['prefix']) & 1
            result[[low, low + 1]] = matrices[address] @ result[[low, low + 1]]
    return result


def _terminal_partition(n, precision, total, a0=2, a1=4):
    """An integer upper schedule: ell=ceil(log2(n+2)), s=ceil(sqrt(L+ell))."""
    assert 0 <= total <= n
    ell = (n + 1).bit_length()
    p = precision + ell
    size = math.isqrt(p - 1) + 1
    reserved = min(total, a0 * ell + a1 * size)
    remaining = reserved
    backward = []
    while remaining < total:
        group = min(size, total - remaining)
        backward.append((remaining, group))
        remaining += group
    return ell, p, size, reserved, list(reversed(backward))


class GroupedTerminalDepthTests(unittest.TestCase):
    def test_native_alphabet_literal_phase_dirty_return_and_actual_inverse(self):
        # target,H,predicate,code0,code1,z,w,borrow: every input is retained.
        word = _letter(0, 1, 2, [3, 4], [5, 6, 7])
        identity = np.eye(256, dtype=complex)
        expected = identity.copy()
        for low in range(0, 256, 2):
            if (low & 6) == 6:
                expected[np.ix_([low, low + 1], [low, low + 1])] = ALPHABET[(low >> 3) & 3]
        actual = _apply_native_word(8, word, identity)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_apply_native_word(8, _adjoint(word), actual),
                                   identity, atol=ATOL, rtol=0)
        # The first dirty T echo alone has the wrong phase when z=1.
        witness = 1 | 2 | 4 | (2 << 3) | (1 << 5)
        wrong = _letter(0, 1, 2, [3, 4], [5, 6, 7], omit_phase_correction=True)
        got = _apply_native_word(8, wrong, identity[:, [witness]])
        self.assertGreater(np.linalg.norm(got - expected[:, [witness]]), 1)

    def test_native_prefetch_preserves_address_and_returns_zero_program(self):
        fixture = _whole_word()
        columns = _active_columns(fixture)
        expected = _ideal_whole_word(fixture, columns)
        actual = _apply_native_word(fixture['width'], fixture['word'], columns)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(
            _apply_native_word(fixture['width'], _adjoint(fixture['word']), actual),
            columns, atol=ATOL, rtol=0)
        # A target in the query address changes under H. Actual unload then
        # reads the wrong row and cannot return the conditionally clean word.
        wrong_fixture = _whole_word(mutable_address=True)
        wrong = _apply_native_word(fixture['width'], wrong_fixture['word'], columns)
        program_mask = sum(1 << bit for bit in fixture['program'])
        occupied = (np.arange(len(wrong)) & program_mask) != 0
        self.assertGreater(np.linalg.norm(wrong[occupied]), .5)

    def test_inactive_occupied_program_predicate_and_dirty_superpositions(self):
        fixture = _whole_word()
        rng = np.random.default_rng(1813)
        columns = rng.normal(size=(1 << fixture['width'], 2)).astype(complex)
        columns += 1j * rng.normal(size=columns.shape)
        active = (np.arange(len(columns)) >> fixture['activity']) & 1
        columns[active != 0] = 0
        columns /= np.linalg.norm(columns, axis=0)
        # Program, predicate, local suffix and dirty helpers are all occupied
        # coherently. Native alphabet equality above separately covers every
        # inactive primitive column, not only these two complete-word vectors.
        actual = _apply_native_word(fixture['width'], fixture['word'], columns)
        np.testing.assert_allclose(actual, columns, atol=ATOL, rtol=0)
        wrong_fixture = _whole_word(missing_h=True)
        wrong = _apply_native_word(fixture['width'], wrong_fixture['word'], columns)
        self.assertGreater(np.linalg.norm(wrong - columns), .5)

    def test_integer_partition_rounding_execution_order_and_square_root_balance(self):
        # Constants here only test integer bookkeeping; they do not assert
        # the theorem's native workspace constants or finite crossover.
        a0, a1 = 2, 4
        nonempty, shortened = 0, 0
        for n in (1, 2, 7, 8, 15, 16, 257, 1024, 10**6, 1 << 4096):
            for precision in (6, 7, 8, 9, 15, 16, 17, 63, 64, 65):
                ell = (n + 1).bit_length()
                p = precision + ell
                size = math.isqrt(p - 1) + 1
                cutoff = a0 * ell + a1 * size
                totals = {0, 1, cutoff - 1, cutoff, cutoff + 1,
                          cutoff + size - 1, cutoff + size, cutoff + size + 1,
                          8 * (precision + (n + precision + 1).bit_length())}
                for total in sorted({min(n, value) for value in totals}):
                    with self.subTest(n_bits=n.bit_length(), precision=precision, total=total):
                        got_ell, got_p, got_size, reserved, groups = _terminal_partition(
                            n, precision, total, a0, a1)
                        self.assertEqual((got_ell, got_p, got_size), (ell, p, size))
                        self.assertLess(1 << (ell - 1), n + 2)
                        self.assertLessEqual(n + 2, 1 << ell)
                        self.assertLess((size - 1) ** 2, p)
                        self.assertLessEqual(p, size ** 2)
                        self.assertEqual(reserved + sum(g for _, g in groups), total)
                        current = total
                        for position, (r, g) in enumerate(groups):
                            self.assertEqual(r + g, current)
                            self.assertTrue(1 <= g <= size)
                            self.assertGreaterEqual(r, reserved)
                            if position:
                                self.assertEqual(g, size)
                            current = r
                        self.assertEqual(current, reserved)
                        self.assertEqual(len(groups), (total - reserved + size - 1) // size)
                        # One serial source allowance p per group, coarse
                        # allowance g^2, and the layerwise reserved-tail Rp.
                        depth_proxy = sum(p + g * g for _, g in groups) + reserved * p
                        self.assertLessEqual(depth_proxy, 2 * total * size + p + reserved * p)
                        if groups:
                            self.assertEqual(reserved, cutoff)
                            nonempty += 1
                            shortened += groups[0][1] < size
        self.assertGreater(nonempty, 100)
        self.assertGreater(shortened, 100)


if __name__ == '__main__':
    unittest.main()
