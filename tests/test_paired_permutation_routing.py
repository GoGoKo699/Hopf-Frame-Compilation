"""Exact bounded checks for the generic paired permutation wrapper.

The native phase-grid fixture reuses R91's gate evaluator and checks every
occupied input. The paired-array tests include noninvolutive permutations,
inactive/padded selectors, actual inverse words, and a negative selector-
mutation example. The returned-address checks allow internal selector
changes and pack ragged addresses before a full dirty-output echo. Exact
full-basis identities imply reference compatibility; these fixtures do not
provide a scalable emitter or a frame-depth theorem.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

from tests.test_parallel_bilinear_word_query import (
    _bilinear_phase, _inverse, _monomial,
)


def _permute(value, permutation):
    """Output site j receives input site permutation[j]."""
    return sum(((value >> source) & 1) << j
               for j, source in enumerate(permutation))


def _inverse_permutation(permutation):
    return tuple(permutation.index(j) for j in range(len(permutation)))


def _shear_word(a, e, permutation, reverse=False):
    """Three completed shears plus swap, or the actual reversed word."""
    inverse = _inverse_permutation(permutation)
    operations = ('up', 'down', 'up', 'swap')
    if reverse:
        operations = operations[::-1]
    for operation in operations:
        if operation == 'up':
            e ^= _permute(a, permutation)
        elif operation == 'down':
            a ^= _permute(e, inverse)
        else:
            a, e = e, a
    return a, e


class PairedPermutationRoutingTests(unittest.TestCase):
    def test_native_sparse_grid_echo_on_every_occupied_input(self):
        # A three-cycle is not its inverse. Its alternative shares rows
        # nontrivially, so the fixture does not reduce to XOR translation.
        permutations = ((1, 2, 0), (2, 1, 0))
        selector = 0
        y = [1, 2]
        a = [3, 4, 5]
        z = [6, 7, 8]
        dirty = [[9 + 3 * u + j for j in range(3)] for u in range(2)]
        matrices = [[[int(source == permutations[u][j])
                      for source in range(3)]
                     for u in range(2)] for j in range(3)]
        b = _bilinear_phase(y, a, z, dirty, matrices)
        # Exact two-entry dirty indicator: Y_0 ^= not r, Y_1 ^= r.
        indicator = [('X', y[0]), ('CX', selector, y[0]),
                     ('CX', selector, y[1])]
        word = b + indicator + _inverse(b) + _inverse(indicator)
        images, phases = _monomial(word, 15)
        inputs = np.arange(1 << 15, dtype=np.int64)
        np.testing.assert_array_equal(images, inputs)
        expected = np.zeros_like(inputs)
        for j in range(3):
            selected = np.where(
                inputs & 1,
                (inputs >> a[permutations[1][j]]) & 1,
                (inputs >> a[permutations[0][j]]) & 1,
            )
            expected ^= selected & ((inputs >> z[j]) & 1)
        np.testing.assert_array_equal(phases, 4 * expected)
        # Hadamard conjugation on all z targets gives the full XOR shear,
        # including every output column, with no residual dirty phase.

    def test_all_three_site_permutations_pair_with_their_inverse(self):
        for permutation in itertools.permutations(range(3)):
            inverse = _inverse_permutation(permutation)
            for a, e in itertools.product(range(8), repeat=2):
                actual = _shear_word(a, e, permutation)
                self.assertEqual(actual, (_permute(a, permutation),
                                          _permute(e, inverse)))
                self.assertEqual(_shear_word(*actual, permutation, reverse=True),
                                 (a, e))

    def test_body_wrapper_returns_every_dirty_array_and_literal_phase(self):
        permutations = tuple(itertools.permutations(range(3)))
        identity = (0, 1, 2)
        # The body is CNOT A_0 -> A_1 followed by Z on A_2.
        # All selector/activity/dirty basis inputs, including padding.
        for label, active, a, e in itertools.product(range(8), range(2),
                                                    range(8), range(8)):
            permutation = (permutations[label]
                           if active and label < len(permutations)
                           else identity)
            routed_a, routed_e = _shear_word(a, e, permutation)
            if active:
                body_a = routed_a ^ ((routed_a & 1) << 1)
                phase = -1 if (body_a >> 2) & 1 else 1
            else:
                body_a, phase = routed_a, 1
            actual_a, actual_e = _shear_word(body_a, routed_e, permutation,
                                            reverse=True)
            expected_a = _permute(body_a, _inverse_permutation(permutation))
            self.assertEqual((actual_a, actual_e), (expected_a, e))
            self.assertEqual(phase, (-1 if active and ((body_a >> 2) & 1) else 1))
            if not active:
                self.assertEqual((actual_a, actual_e, phase), (a, e, 1))

    def test_mutating_selector_can_prevent_auxiliary_return(self):
        identity, cycle = (0, 1, 2), (1, 2, 0)
        a, e = 0, 1
        routed = _shear_word(a, e, cycle)
        # An invalid body changes the selector from cycle to identity.
        # The actual control-dependent inverse then uses the new value.
        after = _shear_word(*routed, identity, reverse=True)
        self.assertNotEqual(after[1], e)

    def test_boundary_return_suffices_after_native_selector_and_address_changes(self):
        # Wires: selector=0, packed A=1,2,3, private query output=4,
        # final output=5. This actual Clifford word temporarily changes both
        # the selector and A_0, but completes the XOR c += r + A_0 + A_1.
        query = [('CX', 2, 1), ('CX', 1, 0), ('CX', 0, 4),
                 ('CX', 1, 0), ('CX', 2, 1)]
        echo = [('CX', 4, 5)]
        body = query + echo + _inverse(query) + echo
        inputs = np.arange(1 << 6, dtype=np.int64)
        images, phases = _monomial(body, 6)
        target = (inputs & 1) ^ ((inputs >> 1) & 1) ^ ((inputs >> 2) & 1)
        np.testing.assert_array_equal(images, inputs ^ (target << 5))
        np.testing.assert_array_equal(phases, np.zeros_like(inputs))
        during, _ = _monomial(query[:2], 6)
        self.assertTrue(np.any((during ^ inputs) & 1))
        self.assertTrue(np.any((during ^ inputs) & 2))

        permutations = ((0, 1, 2), (1, 2, 0))
        for label, a, e, private, output in itertools.product(
                range(2), range(8), range(8), range(2), range(2)):
            permutation = permutations[label]
            routed_a, routed_e = _shear_word(a, e, permutation)
            body_input = label | (routed_a << 1) | (private << 4) | (output << 5)
            body_output = int(images[body_input])
            actual_label = body_output & 1
            actual_a = (body_output >> 1) & 7
            # Use the actual selector at the inverse boundary, as a circuit
            # does; intermediate selector values need not equal this one.
            restored_a, restored_e = _shear_word(
                actual_a, routed_e, permutations[actual_label], reverse=True)
            selected = label ^ ((a >> permutation[0]) & 1) ^ ((a >> permutation[1]) & 1)
            self.assertEqual((actual_label, restored_a, restored_e), (label, a, e))
            self.assertEqual((body_output >> 4) & 1, private)
            self.assertEqual((body_output >> 5) & 1, output ^ selected)

    def test_ragged_packing_full_dirty_output_echo_and_all_selector_rows(self):
        # Three physical address sites are needed in the union, but at most
        # two in any row. The common P bit is not part of the routed array.
        addresses = ((2,), (0, 1), (1,), ())
        # Explicit total two-output tables, first indexed by selector, then P.
        # Label 2 represents a malformed legacy row with a nonzero extension;
        # label 3 is an explicitly defined invalid zero row. Neither can be
        # dropped as an initialized-input case. H=0 separately specifies zero.
        tables = (((0, 3), (1, 2)),
                  ((2, 1, 3, 0), (3, 0, 2, 1)),
                  ((1, 2), (3, 0)),
                  ((0,), (0,)))
        self.assertEqual(len(set().union(*map(set, addresses))), 3)
        self.assertEqual(max(map(len, addresses)), 2)

        for label, active, prefix, a, e, private, output in itertools.product(
                range(4), range(2), range(2), range(8), range(8), range(4), range(4)):
            address = addresses[label] if active else ()
            permutation = address + tuple(j for j in range(3) if j not in address)
            routed_a, routed_e = _shear_word(a, e, permutation)
            packed_row = routed_a & ((1 << len(address)) - 1)
            loaded = tables[label][prefix][packed_row] if active else 0

            # Complete query/CNOT/actual-unquery/CNOT, retaining the unknown
            # two-bit output of the query. The premise supplies this exact
            # table XOR; this fixture checks its surrounding physical packing.
            actual_private = private ^ loaded
            actual_output = output ^ actual_private
            actual_private ^= loaded
            actual_output ^= actual_private
            restored_a, restored_e = _shear_word(
                routed_a, routed_e, permutation, reverse=True)

            original_row = sum(((a >> wire) & 1) << j
                               for j, wire in enumerate(addresses[label]))
            expected = tables[label][prefix][original_row] if active else 0
            self.assertEqual((restored_a, restored_e, actual_private), (a, e, private))
            self.assertEqual(actual_output, output ^ expected)
            if not active or label == 3:
                self.assertEqual((routed_a, routed_e, actual_output), (a, e, output))


if __name__ == '__main__':
    unittest.main()
