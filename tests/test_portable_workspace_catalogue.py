"""Bounded operator checks for carrying workspace through a catalogue.

Native H/S/T/CNOT words retain literal phase and use their actual inverses.
Complete operators have at most eight wires. The small leakage and frame
fixtures check operator interfaces, not a native approximation emitter.
The separate R91 tests cover dirty-program echoes and the paired router;
these checks do not establish asymptotic counts, depth, or allocation.
"""
from __future__ import annotations

import unittest

import numpy as np

from compiler_robust_hopf.frames import hopf_ry, real_frame_matrix

try:
    from .test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis, _word_matrix,
    )
except ImportError:
    from test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis, _word_matrix,
    )


ATOL = 3e-12
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)


def _shift_word(word, offset):
    return [(gate[0], *(wire + offset for wire in gate[1:]))
            for gate in word]


def _swap_word(pairs, control=None):
    word = []
    for left, right in pairs:
        middle = ('CX', right, left) if control is None else (
            'CCX', control, right, left)
        word += [('CX', left, right), middle, ('CX', left, right)]
    return _expand_toffolis(word)


def _transfer_word(array_word, selected_swap):
    # Chronological order for K = S A_all S A_all^dagger.
    return _adjoint(array_word) + selected_swap + array_word + selected_swap


def _commutator_word(a, b):
    # Chronological order for A B A^dagger B^dagger.
    return _adjoint(b) + _adjoint(a) + b + a


def _native_operator(width, word):
    return _apply_native_word(width, word, np.eye(1 << width, dtype=complex))


def _direct_sum(left, right):
    zero = np.zeros((len(left), len(right)), dtype=complex)
    return np.block([[left, zero], [zero.T, right]])


class PortableWorkspaceCatalogueTests(unittest.TestCase):
    def test_complete_four_field_payload_and_actual_inverse(self):
        """All c,X,d,D and slot columns, including occupied flags, return."""
        local_word = [
            ('H', 0), ('T', 0), ('CX', 0, 1), ('TDG', 1),
            ('CX', 1, 2), ('H', 2), ('T', 2), ('CX', 2, 3),
            ('S', 3), ('CX', 3, 0), ('T', 0),
        ]
        a = _native_operator(4, local_word)
        b = _word_matrix(4, [('X', 0)])
        array_word = _shift_word(local_word, 4)
        swap = _swap_word([(wire, wire + 4) for wire in range(4)])
        transfer = _transfer_word(array_word, swap)
        actual_k = _native_operator(8, transfer)
        # Logical fields are the four low bits, the slot the four high bits.
        np.testing.assert_allclose(actual_k, np.kron(a.conj().T, a),
                                   atol=ATOL, rtol=0)
        actual = _native_operator(8, _commutator_word(transfer, [('X', 0)]))
        expected = np.kron(np.eye(16), a @ b @ a.conj().T @ b)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_native_operator(8, _adjoint(transfer)),
                                   actual_k.conj().T, atol=ATOL, rtol=0)
        # This particular native A is not self-adjoint: forwarding the
        # same word a second time cannot stand in for its actual inverse.
        wrong = _native_operator(8, [('X', 0)] + transfer + [('X', 0)] + transfer)
        self.assertGreater(np.linalg.norm(wrong - expected), .1)

    def test_two_noncommuting_entries_and_coherent_selection(self):
        """A smaller complete payload checks selection on every prefix column."""
        words = ([('H', 0), ('T', 0), ('CX', 0, 1), ('S', 1)],
                 [('H', 1), ('TDG', 1), ('CX', 1, 0), ('T', 0)])
        entries = [_native_operator(2, word) for word in words]
        self.assertGreater(np.linalg.norm(entries[0] @ entries[1]
                                          - entries[1] @ entries[0]), .1)
        # Logical c,d occupy 0,1; slots occupy 2,3 and 4,5; P is bit 6.
        # A controlled slot exchange followed by a fixed swap and its
        # inverse is the selected swap. R91 separately checks the native
        # dirty-program/paired-router realization of this interface.
        shift = _swap_word([(2, 4), (3, 5)], control=6)
        selected_swap = shift + _swap_word([(0, 2), (1, 3)]) + _adjoint(shift)
        array_word = _shift_word(words[0], 2) + _shift_word(words[1], 4)
        transfer = _transfer_word(array_word, selected_swap)
        actual = _native_operator(7, _commutator_word(transfer, [('X', 0)]))
        b = _word_matrix(2, [('X', 0)])
        commutators = [a @ b @ a.conj().T @ b for a in entries]
        expected = _direct_sum(*(np.kron(np.eye(16), u) for u in commutators))
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(commutators[0] - commutators[1]), .1)
        # A complete-operator identity also applies to coherent P and all
        # array/reference entanglement; the unknown slot flags are never
        # restricted to zero in this comparison.

    def test_controlled_hopf_completion_and_initialized_leakage_bound(self):
        """Check all marker columns, branch direction, and a leaking inverse."""
        alpha, beta, gamma = .29, -.43, .67
        frame = real_frame_matrix([alpha, beta, gamma]).astype(complex)
        controlled = real_frame_matrix([0, 0, alpha, 0, 0, beta, gamma])
        expected = _direct_sum(np.eye(4), frame)
        np.testing.assert_allclose(controlled, expected, atol=ATOL, rtol=0)
        flip = np.kron(X, np.eye(4))
        commutator = controlled @ flip @ controlled.conj().T @ flip
        np.testing.assert_allclose(commutator,
                                   _direct_sum(frame.conj().T, frame),
                                   atol=ATOL, rtol=0)

        # Basis order is d,c,X0,X1,D. Clean input constrains only d.
        ideal = np.kron(controlled, I2)
        local_flip = np.kron(flip, I2)
        occupied = _word_matrix(4, [('H', 0), ('CX', 0, 3), ('T', 3), ('S', 1)])
        extension = _direct_sum(ideal, occupied)
        generator = np.kron(X, np.kron(X, np.kron(Z, np.kron(X, Z))))
        angle = .017
        leakage = np.cos(angle) * np.eye(32) + 1j * np.sin(angle) * generator
        actual = leakage @ extension
        embedding = np.eye(32, dtype=complex)[:, :16]
        forward_error = np.linalg.norm(actual @ embedding - embedding @ ideal, ord=2)
        inverse_error = np.linalg.norm(actual.conj().T @ embedding
                                       - embedding @ ideal.conj().T, ord=2)
        self.assertAlmostEqual(forward_error, 2 * np.sin(angle / 2), delta=ATOL)
        self.assertGreater(np.linalg.norm(actual[16:, :16]), .01)
        self.assertGreater(np.linalg.norm(actual - np.kron(I2, ideal)), .1)
        self.assertAlmostEqual(inverse_error, forward_error, delta=ATOL)
        b = np.kron(I2, local_flip)
        result = actual @ b @ actual.conj().T @ b @ embedding
        target = embedding @ (ideal @ local_flip @ ideal.conj().T @ local_flip)
        self.assertLessEqual(np.linalg.norm(result - target, ord=2),
                             2 * forward_error + ATOL)

    def test_partial_payload_does_not_return_the_complete_array(self):
        """Refute exact full-space return at one partial-routing interface."""
        # X is a spectator. A is controlled-d SWAP(c,D), exactly identity
        # on all d=0 inputs, including arbitrary c,D and their references.
        # Logical c,d,D are bits 0,1,2; their slot copies are bits 3,4,5.
        local = _swap_word([(0, 2)], control=1)
        a = _native_operator(3, local)
        clean_columns = [value for value in range(8) if not ((value >> 1) & 1)]
        np.testing.assert_allclose(a[:, clean_columns], np.eye(8)[:, clean_columns],
                                   atol=ATOL, rtol=0)
        array_word = _shift_word(local, 3)
        full_swap = _swap_word([(0, 3), (1, 4), (2, 5)])
        partial_swap = _swap_word([(0, 3), (1, 4)])
        full = _native_operator(6, _commutator_word(
            _transfer_word(array_word, full_swap), [('X', 0)]))
        partial = _native_operator(6, _commutator_word(
            _transfer_word(array_word, partial_swap), [('X', 0)]))
        # Independent exact formulas: the partial route flips the slot D'
        # rather than logical D when logical d=1. It leaves a dirty-array
        # residue even though the local compiler's initialized action is I.
        expected = _word_matrix(6, [('CX', 1, 0), ('CX', 1, 2)])
        wrong_target = _word_matrix(6, [('CX', 1, 0), ('CX', 1, 5)])
        np.testing.assert_allclose(full, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(partial, wrong_target, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(partial - full), .1)
        self.assertEqual(int(np.argmax(abs(partial[:, 2]))), 2 ^ 1 ^ 32)
        # This is not a lower bound on all partial-routing constructions:
        # the two operators agree on this example's initialized d=0 sector.
        initialized = [value for value in range(64) if not ((value >> 1) & 1)]
        np.testing.assert_allclose(partial[:, initialized], full[:, initialized],
                                   atol=ATOL, rtol=0)

    def test_four_mode_group_reads_the_current_selector(self):
        """A retained initial selector fails after an earlier group changes it."""
        alpha, beta, gamma = .31, -.47, .73
        root = np.eye(4, dtype=complex)
        root[np.ix_([0, 2], [0, 2])] = hopf_ry(alpha)
        children = _direct_sum(hopf_ry(beta), hopf_ry(gamma))
        expected = real_frame_matrix([alpha, beta, gamma])
        np.testing.assert_allclose(children @ root, expected, atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(children @ root - root @ children), .1)
        # Save the initial high logical bit, execute the root, then use
        # that frozen copy for the children and try to uncompute the copy.
        cache = _word_matrix(3, [('CX', 1, 2)])
        frozen_children = _direct_sum(np.kron(I2, hopf_ry(beta)),
                                      np.kron(I2, hopf_ry(gamma)))
        embedding = np.eye(8, dtype=complex)[:, :4]
        frozen = cache.conj().T @ frozen_children @ np.kron(I2, root) @ cache
        actual = frozen @ embedding
        target = embedding @ expected
        self.assertGreater(np.linalg.norm(actual - target), .1)
        self.assertGreater(np.linalg.norm(actual[4:]), .1)

    def test_guard_conjugation_preserves_original_flag_without_factoring_it(self):
        """Guard return is exact despite payload leakage; occupied G can matter."""
        # Bits are c,d,G,Y,P,H. Only the logical c,d payload changes in M.
        # C_G toggles G by H*[Y=0], excluding every changed payload bit.
        width, guard, lower, prefix, saved = 6, 2, 3, 4, 5
        compute = _expand_toffolis([
            ('X', lower), ('CCX', saved, lower, guard), ('X', lower),
        ])
        c = _native_operator(width, compute)
        local = [
            _native_operator(2, [('H', 0), ('CX', 0, 1), ('T', 1)]),
            _native_operator(2, [('H', 1), ('CX', 1, 0), ('TDG', 0)]),
        ]
        m = np.zeros((64, 64), dtype=complex)
        target = np.zeros_like(m)
        for fields in range(16):
            g, y, p, h = ((fields >> bit) & 1 for bit in range(4))
            block = slice(4 * fields, 4 * fields + 4)
            # The reserved inactive code is exactly I on every c,d input.
            m[block, block] = local[p] if h and g else np.eye(4)
            active_after_compute = h and (g ^ (h and not y))
            target[block, block] = local[p] if active_after_compute else np.eye(4)
        actual = c.conj().T @ m @ c
        np.testing.assert_allclose(actual, target, atol=ATOL, rtol=0)
        indices = np.arange(64)
        g0 = indices[((indices >> guard) & 1) == 0]
        g1 = indices[((indices >> guard) & 1) == 1]
        np.testing.assert_allclose(actual[np.ix_(g0, g1)], 0, atol=ATOL, rtol=0)
        np.testing.assert_allclose(actual[np.ix_(g1, g0)], 0, atol=ATOL, rtol=0)

        # H=0 is full identity even when c,d,G are occupied or coherent.
        inactive = indices[((indices >> saved) & 1) == 0]
        np.testing.assert_allclose(actual[:, inactive], np.eye(64)[:, inactive],
                                   atol=ATOL, rtol=0)
        # For H=1,G=0 the output G returns to zero on every input column,
        # although this explicit local unitary leaks from both c=0,d=0.
        initialized = g0[((g0 >> saved) & 1) == 1]
        np.testing.assert_allclose(actual[np.ix_(g1, initialized)], 0,
                                   atol=ATOL, rtol=0)
        witness = 1 << saved  # H=1, P=Y=G=c=d=0.
        output = actual[:, witness]
        self.assertGreater(np.linalg.norm(output[(indices & 1) != 0]), .5)
        self.assertGreater(np.linalg.norm(output[(indices & 2) != 0]), .5)

        # Diagonal in original G does not mean an operator tensor I_G:
        # at H=1,Y=0,P=0, G=0 executes local[0], whereas G=1 executes I.
        block0 = np.arange(witness, witness + 4)
        block1 = block0 + (1 << guard)
        np.testing.assert_allclose(actual[np.ix_(block0, block0)], local[0],
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(actual[np.ix_(block1, block1)], np.eye(4),
                                   atol=ATOL, rtol=0)
        self.assertGreater(np.linalg.norm(local[0] - np.eye(4)), .1)

        # A concrete forbidden predicate reads d, which M changes. Its
        # actual uncomputation then leaves a nonzero G=1 component.
        bad_compute = _expand_toffolis([
            ('X', 1), ('CCX', saved, 1, guard), ('X', 1),
        ])
        bad_c = _native_operator(width, bad_compute)
        bad_output = (bad_c.conj().T @ m @ bad_c)[:, witness]
        self.assertGreater(np.linalg.norm(bad_output[g1]), .5)


if __name__ == '__main__':
    unittest.main()
