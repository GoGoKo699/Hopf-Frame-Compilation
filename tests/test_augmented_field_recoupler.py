"""Exact augmented field recoupling, with arbitrary moving helpers.

Full computational-basis columns check the algebraic primitive word at p=2,3.
A separately emitted p=2 native word and an overlapping shared-target phase
batch verify the native composition point. No all-size depth improvement is
inferred from the bounded matrices; the proof owns the uniform ledger.
"""
from functools import lru_cache
import unittest

import numpy as np

from tests.test_dirty_indicator_depth import _exact_phase_action
from tests.test_field_shear_intertwiner import _translation
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_projective_field_recoupler import (
    fourier_word, linear_word, ring_matrix, equal_ring, adjoint_ring)
from tests.test_recursive_field_shear_decoder import (
    _tower_product as mul, _tower_trace as trace, _address_field, _tower_shear)


@lru_cache(None)
def _addition_value(level, u, v, which):
    ell = 1 << (level - 2)
    low = (1 << ell) - 1
    u0, u1, v0, v1 = u & low, u >> ell, v & low, v >> ell
    if which == 'before':
        return (mul(level, u, v) >> ell) ^ mul(level - 1, u0, u1)
    assert which == 'after'
    return mul(level - 1, u0, v1) ^ mul(level - 1, u1, v0)


@lru_cache(None)
def _fourier_signs(level, multiplier, source):
    return tuple(-1 if trace(level, mul(level, multiplier, mul(level, source, target)))
                 else 1 for target in range(1 << (1 << (level - 1))))


def _algebraic_word(level, top, offset=0):
    if level == 1:
        return [('CX', offset + 1, offset)]
    ell, degree = 1 << (level - 2), 1 << (level - 1)
    helper = (1 << top) + ell - 1
    recoupler = [('F', level, offset + degree, 1),
                 ('ADD', level, offset, helper, 'before'),
                 ('F', level, offset + degree, 1 << ell),
                 ('ADD', level, offset, helper, 'after')]
    recoupler += [gate for j in range(degree)
                  for gate in [('CX', offset + j, offset + degree + j),
                               ('CX', offset + degree + j, offset + j)]]
    recoupler += [('F', level - 1, offset + ell, 1),
                  ('F', level - 1, offset + degree + ell, 1)]
    return (recoupler + _algebraic_word(level - 1, top, offset)
            + _algebraic_word(level - 1, top, offset + degree))


def _apply_algebraic(operations, column):
    state = {column: 1}
    normalization = 0
    for kind, *args in operations:
        if kind == 'CX':
            control, target = args
            state = {row ^ (((row >> control) & 1) << target): value
                     for row, value in state.items()}
        elif kind == 'ADD':
            level, offset, helper, which = args
            ell, degree = 1 << (level - 2), 1 << (level - 1)
            mask = (1 << degree) - 1
            changed = {}
            for row, value in state.items():
                u, v = (row >> offset) & mask, (row >> (offset + degree)) & mask
                delta = _addition_value(level, u, v, which)
                target = row ^ (delta << helper)
                assert target not in changed  # The whole add is a permutation.
                changed[target] = value
            state = changed
        else:
            assert kind == 'F'
            level, offset, multiplier = args
            degree = 1 << (level - 1)
            mask = ((1 << degree) - 1) << offset
            normalization += degree
            changed = {}
            for row, value in state.items():
                source = (row & mask) >> offset
                for target, sign in enumerate(_fourier_signs(level, multiplier, source)):
                    out = (row & ~mask) | (target << offset)
                    changed[out] = changed.get(out, 0) + sign * value
            state = {row: value for row, value in changed.items() if value}
    return state, normalization


def _helper_phase(top, address, basis):
    exponent = 0
    for level in range(2, top + 1):
        ell = 1 << (level - 2)
        helper = (1 << top) + ell - 1
        value = (basis >> helper) & ((1 << ell) - 1)
        coefficient = _address_field(level - 1, address & ((1 << (level - 1)) - 1))
        exponent ^= trace(level - 1, mul(level - 1, coefficient, value))
    return -1 if exponent else 1


def _parity_phase_batch(forms, helpers):
    """Returned eight-phase masks for arbitrary overlapping linear forms."""
    assert len(forms) == len(helpers) and len(set(helpers)) == len(helpers)
    assert all(not (mask >> h) & 1
               for row in forms for mask in row for h in helpers)
    emitted, previous = [], 0
    for step in range(8):
        vertex = step ^ (step >> 1)
        difference = previous ^ vertex
        if difference:
            index = difference.bit_length() - 1
            emitted += [('CX', bit, h)
                        for row, h in zip(forms, helpers)
                        for bit in range(row[index].bit_length())
                        if (row[index] >> bit) & 1]
        emitted += [('TDG' if vertex.bit_count() & 1 else 'T', h)
                    for h in helpers]
        previous = vertex
    emitted += [('CX', bit, h)
                for row, h in zip(forms, helpers)
                for bit in range(row[2].bit_length())
                if (row[2] >> bit) & 1]
    return emitted


def _p2_forms(offset, helper, which):
    u0, u1, v0, v1 = [1 << (offset + j) for j in range(4)]
    h = 1 << helper
    if which == 'before':
        return [(u0, v1 ^ u1, h), (u1, v0 ^ v1, h)]
    assert which == 'after'
    return [(u0, v1, h), (u1, v0, h)]


def _native_p2_word():
    """32 T / 16 bands; four core, one moving, two returned private wires."""
    emitted = []
    for kind, *args in _algebraic_word(2, 2):
        if kind == 'CX':
            emitted.append((kind, *args))
        elif kind == 'ADD':
            level, offset, helper, which = args
            assert level == 2
            emitted += ([('H', helper)]
                        + _parity_phase_batch(_p2_forms(offset, helper, which), [5, 6])
                        + [('H', helper)])
        else:
            level, offset, multiplier = args
            degree = 1 << (level - 1)
            if multiplier != 1:
                emitted += linear_word(
                    [mul(level, multiplier, 1 << j) for j in range(degree)],
                    range(offset, offset + degree))
            emitted += fourier_word(level, offset)
    return emitted


@lru_cache(None)
def _complete_algebraic_columns(top):
    width = (1 << top) + (1 << (top - 1)) - 1
    operations = _algebraic_word(top, top)
    columns, normalization = [], None
    for column in range(1 << width):
        state, exponent = _apply_algebraic(operations, column)
        assert sum(value * value for value in state.values()) == 1 << exponent
        assert normalization is None or normalization == exponent
        normalization = exponent
        columns.append(state)
    return columns, normalization


class AugmentedFieldRecouplerTests(unittest.TestCase):
    def test_relative_trace_two_product_formula_and_zero_cell(self):
        for level in (2, 3):
            ell = 1 << (level - 2)
            low = (1 << ell) - 1
            size = 1 << (2 * ell)
            for u in range(size):
                for v in range(size):
                    u0, u1, v0, v1 = u & low, u >> ell, v & low, v >> ell
                    self.assertEqual(_addition_value(level, u, v, 'before'),
                                     mul(level - 1, u0, v1 ^ u1)
                                     ^ mul(level - 1, u1, v0 ^ v1))
                    self.assertEqual(_addition_value(level, u, v, 'after'),
                                     _addition_value(level, u, v ^ u, 'after'))
            for which in ('before', 'after'):
                self.assertTrue(all(_addition_value(level, 0, v, which) == 0
                                    for v in range(size)))

    def test_native_shared_target_batch_and_arbitrary_private_return(self):
        # Two disjoint four-wire cells share one target and four private masks.
        forms = _p2_forms(0, 8, 'before') + _p2_forms(4, 8, 'before')
        emitted = _parity_phase_batch(forms, [9, 10, 11, 12])
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in emitted), 32)
        self.assertEqual(len(_bands(emitted)), 8)
        for native in (emitted, _inverse(emitted)):
            for basis in range(1 << 13):
                phase = 4 * (sum(all((basis & mask).bit_count() & 1 for mask in row)
                                 for row in forms) & 1)
                self.assertEqual(_exact_phase_action(basis, native), (basis, phase))

    def test_complete_p2_native_word_and_actual_inverse(self):
        emitted = _native_p2_word()
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in emitted), 32)
        self.assertEqual(len(_bands(emitted)), 16)
        columns, exponent = _complete_algebraic_columns(2)
        self.assertEqual(exponent % 2, 0)
        expected = np.zeros((4, 128, 128), dtype=np.int64)
        for private in range(4):
            for column, entries in enumerate(columns):
                for row, value in entries.items():
                    expected[0, row | (private << 5), column | (private << 5)] = value
        actual = ring_matrix(emitted, 7)
        self.assertTrue(equal_ring(actual, (expected, exponent // 2)))
        self.assertTrue(equal_ring(ring_matrix(_inverse(emitted), 7), adjoint_ring(actual)))

    def test_all_computational_helper_columns_and_reversed_word(self):
        for top in (2, 3):
            columns, exponent = _complete_algebraic_columns(top)
            expected_inverse = [{} for _ in columns]
            for column, entries in enumerate(columns):
                for row, value in entries.items():
                    expected_inverse[row][column] = value
            reversed_word = list(reversed(_algebraic_word(top, top)))
            for column in range(len(columns)):
                inverse, inverse_exponent = _apply_algebraic(reversed_word, column)
                self.assertEqual(inverse_exponent, exponent)
                self.assertEqual(inverse, expected_inverse[column], (top, column))

    def test_augmented_intertwining_every_address_and_full_column(self):
        for top in (2, 3):
            columns, _ = _complete_algebraic_columns(top)
            core_mask = (1 << (1 << top)) - 1
            for address in range(1 << top):
                for column in range(len(columns)):
                    source = ((column & ~core_mask)
                              | _tower_shear(top, address, column & core_mask))
                    sign = _helper_phase(top, address, column)
                    left = {row: sign * value for row, value in columns[source].items()}
                    right = {((row & ~core_mask)
                              | _translation(top, address, row & core_mask)):
                             _helper_phase(top, address, row) * value
                             for row, value in columns[column].items()}
                    self.assertEqual(left, right, (top, address, column))

    def test_helper_query_is_essential_and_echo_cancels_it(self):
        columns, _ = _complete_algebraic_columns(2)
        without_query_fails = False
        for address in range(4):
            for column in range(32):
                source = (column & 16) | _tower_shear(2, address, column & 15)
                moved = {(row & 16) | _translation(2, address, row & 15): value
                         for row, value in columns[column].items()}
                without_query_fails |= columns[source] != moved
                # The target augmented query is a signed involution; its
                # actual Z0 echo has exactly Z_address and no helper phase.
                middle = (column & 16) | _translation(2, address, column & 15)
                phase = (_helper_phase(2, address, column)
                         * (-1 if middle & 1 else 1)
                         * _helper_phase(2, address, middle))
                self.assertEqual(phase, -1 if (column >> address) & 1 else 1)
        self.assertTrue(without_query_fails)


if __name__ == '__main__':
    unittest.main()
