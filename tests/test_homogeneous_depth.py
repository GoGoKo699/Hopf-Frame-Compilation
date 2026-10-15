"""Bounded algebra audits for a depth-homogeneous full-frame word.

The reference frame is built from suffix-zero-controlled local layers and,
independently, from its literal marker columns. The five-stage tag word is
checked for general complex SU(2) data. A separate physical selected-source
fixture uses actual one-qubit Clifford+T words V, with U=(V Z V^dagger) Z;
this restricted factorization is not asserted for every SU(2) matrix.

The query fixture is its exact reversible truth-table action on arbitrary
dirty masks. These bounded checks supply neither a scalable query emitter
nor a proof of the asymptotic depth or workspace bounds. No scientific
receipt is changed.
"""
from __future__ import annotations

from fractions import Fraction as F
import unittest

import numpy as np

try:
    from .test_operator_source_compiler import _word_matrix
    from .test_subtree_projector_decoder import _geometry, REGULAR, SINGULAR
except ImportError:
    from test_operator_source_compiler import _word_matrix
    from test_subtree_projector_decoder import _geometry, REGULAR, SINGULAR


ATOL = 3e-11
Z = np.diag([1, -1]).astype(complex)
NATIVE_WORDS = (
    [('H', 0), ('T', 0), ('H', 0)],
    [('T', 0), ('H', 0), ('S', 0), ('T', 0)],
    [('H', 0), ('TDG', 0), ('H', 0), ('T', 0)],
)


def _last_one(word, width):
    """One-indexed last occupied position, with the first qubit on the left."""
    return width - ((word & -word).bit_length() - 1) if word else 0


def _tensor(factors):
    result = np.ones((1, 1), dtype=complex)
    for factor in factors:
        result = np.kron(result, factor)
    return result


def _su2(theta, alpha=0.0, beta=0.0):
    a = np.exp(1j * alpha) * np.cos(theta)
    b = np.exp(1j * beta) * np.sin(theta)
    return np.array([[a, -b.conjugate()], [b, a.conjugate()]])


def _layer_frame(unitaries):
    """Independent complete frame from its chronological guarded layers."""
    width = len(unitaries)
    result = np.eye(1 << width, dtype=complex)
    for position, unitary in enumerate(unitaries):
        bit = 1 << (width - position - 1)
        for low in range(len(result)):
            if low & bit or low & (bit - 1):
                continue
            pair = np.stack((result[low].copy(), result[low | bit].copy()))
            result[low], result[low | bit] = unitary @ pair
    return result


def _column_frame(unitaries):
    """Literal prefix, marker, and descendant product columns."""
    width = len(unitaries)
    result = np.zeros((1 << width, 1 << width), dtype=complex)
    for word in range(len(result)):
        tag = _last_one(word, width)
        factors = []
        for position, unitary in enumerate(unitaries, 1):
            bit = (word >> (width - position)) & 1
            vector = np.eye(2, dtype=complex)[:, bit]
            factors.append(vector if position < tag else unitary @ vector)
        column = np.ones(1, dtype=complex)
        for factor in factors:
            column = np.kron(column, factor)
        result[:, word] = column
    return result


def _data(matrix, state):
    """Apply to the first axis of (X,T,dirty,batch)."""
    return np.tensordot(matrix, state, axes=(1, 0))


def _tag_xor(state, *, active):
    width = (state.shape[0] - 1).bit_length()
    result = np.empty_like(state)
    for word in range(state.shape[0]):
        tag = _last_one(word, width) if active else 0
        for old in range(state.shape[1]):
            result[word, old ^ tag] = state[word, old]
    return result


def _selected_reference(state, unitaries, *, active):
    result = np.empty_like(state)
    for tag in range(state.shape[1]):
        matrix = _tensor([unitary if active and tag <= position else np.eye(2)
                          for position, unitary in enumerate(unitaries, 1)])
        result[:, tag] = _data(matrix, state[:, tag])
    return result


def _mask_query(state, *, active):
    """Exact XOR query; its address is the preserved tag, never the data."""
    width = (state.shape[0] - 1).bit_length()
    result = np.empty_like(state)
    for tag in range(state.shape[1]):
        mask = sum(int(active and tag <= position + 1) << position
                   for position in range(width))
        for dirty in range(state.shape[2]):
            result[:, tag, dirty ^ mask] = state[:, tag, dirty]
    return result


def _dirty_cz_bank(state):
    width = (state.shape[0] - 1).bit_length()
    result = state.copy()
    for word in range(state.shape[0]):
        for dirty in range(state.shape[2]):
            parity = sum(((word >> (width - position - 1)) & 1)
                         * ((dirty >> position) & 1) for position in range(width))
            result[word, :, dirty] *= (-1) ** parity
    return result


def _physical_selected_source(state, conjugators, *, active):
    """Selected Z, then selected V Z V^dagger, with dirty-query echoes.

    Both echoes execute on disabled branches, including the actual V and
    its adjoint; they are not replaced with identity in the fixture.
    """
    bank = _tensor(conjugators)
    for conjugated in (False, True):
        for _ in range(2):
            state = _mask_query(state, active=active)
            if conjugated:
                state = _data(bank.conj().T, state)
            state = _dirty_cz_bank(state)
            if conjugated:
                state = _data(bank, state)
    return state


def _five_stage_word(state, unitaries, *, active, conjugators=None):
    state = _tag_xor(state, active=active)
    if conjugators is None:
        state = _selected_reference(state, unitaries, active=active)
    else:
        state = _physical_selected_source(state, conjugators, active=active)
    bank = _tensor(unitaries)
    # These two banks are unconditional actual inverses. On H=0 they
    # cancel even for occupied tags and arbitrary dirty/reference inputs.
    state = _data(bank.conj().T, state)
    state = _tag_xor(state, active=active)
    return _data(bank, state)


def _initialized_tag_columns(width, dirty_width=0):
    size, labels, dirty_size = 1 << width, 1 << width.bit_length(), 1 << dirty_width
    result = np.zeros((size, labels, dirty_size, size * dirty_size), dtype=complex)
    for word in range(size):
        for dirty in range(dirty_size):
            result[word, 0, dirty, word * dirty_size + dirty] = 1
    return result


def _all_columns(width):
    size, labels = 1 << width, 1 << width.bit_length()
    dimension = size * labels * size
    return np.eye(dimension, dtype=complex).reshape(size, labels, size, dimension)


def _sign(value):
    return -1 if value < 0 else 1  # In particular sign(0)=+1.


def _tree_signs(width, pairs):
    leaf = []
    marker = [1] * (1 << width)
    for word in range(1 << width):
        sign = 1
        for depth in range(width):
            prefix = word >> (width - depth)
            pair = pairs[(1 << depth) - 1 + prefix]
            branch = (word >> (width - depth - 1)) & 1
            sign *= _sign(pair[branch])
        leaf.append(sign)
    for depth in range(width):
        for prefix in range(1 << depth):
            incoming = 1
            for ancestor in range(depth):
                node = (1 << ancestor) - 1 + (prefix >> (depth - ancestor))
                branch = (prefix >> (depth - ancestor - 1)) & 1
                incoming *= _sign(pairs[node][branch])
            cosine, sine = pairs[(1 << depth) - 1 + prefix]
            word = (2 * prefix + 1) << (width - depth - 1)
            marker[word] = incoming * _sign(cosine) * _sign(sine)
    return np.asarray(leaf), np.asarray(marker)


class HomogeneousDepthTests(unittest.TestCase):
    def test_complex_frame_layers_match_literal_marker_columns(self):
        cases = (
            [_su2(theta) for theta in (0.37, -1.2, 2.3)],
            [_su2(theta) for theta in (np.pi / 2, np.pi, 0.0)],
            [_su2(*parameters) for parameters in
             ((0.37, 0.7, -0.2), (-1.2, -0.3, 0.9), (2.3, 0.5, -0.8))],
        )
        self.assertGreater(abs(cases[-1][0][0, 0].imag), 0.2)
        for width in (1, 2, 3):
            for case in cases:
                with self.subTest(width=width, complex_diagonal=case[0][0, 0].imag):
                    frame = _layer_frame(case[:width])
                    np.testing.assert_allclose(frame, _column_frame(case[:width]), atol=ATOL)
                    np.testing.assert_allclose(frame.conj().T @ frame,
                                               np.eye(1 << width), atol=ATOL)
                    initial = _initialized_tag_columns(width)
                    actual = _five_stage_word(initial, case[:width], active=True)
                    # Equality of all initialized columns checks literal
                    # phases and tag return, also after adjoining a reference.
                    np.testing.assert_allclose(actual, _data(frame, initial), atol=ATOL)

    def test_selected_native_echo_on_all_dirty_and_tag_columns(self):
        for width in (1, 2, 3):
            conjugators = [_word_matrix(1, word) for word in NATIVE_WORDS[:width]]
            unitaries = [v @ Z @ v.conj().T @ Z for v in conjugators]
            self.assertGreater(max(np.max(np.abs(u.imag)) for u in unitaries), 0.1)
            initial = _all_columns(width)
            for active in (False, True):
                with self.subTest(width=width, active=active):
                    actual = _physical_selected_source(initial, conjugators, active=active)
                    expected = _selected_reference(initial, unitaries, active=active)
                    np.testing.assert_allclose(actual, expected, atol=ATOL)
            # width=2 includes the padded tag 3, on which every selected
            # involution is disabled although the dirty masks are occupied.

    def test_full_native_word_returns_dirty_masks_and_inactive_occupied_tags(self):
        for width in (1, 2, 3):
            conjugators = [_word_matrix(1, word) for word in NATIVE_WORDS[:width]]
            unitaries = [v @ Z @ v.conj().T @ Z for v in conjugators]
            initial = _initialized_tag_columns(width, dirty_width=width)
            actual = _five_stage_word(initial, unitaries, active=True,
                                      conjugators=conjugators)
            expected = _data(_layer_frame(unitaries), initial)
            np.testing.assert_allclose(actual, expected, atol=ATOL)
            occupied = _all_columns(width)
            inactive = _five_stage_word(occupied, unitaries, active=False,
                                        conjugators=conjugators)
            np.testing.assert_allclose(inactive, occupied, atol=ATOL)

    def test_signed_tree_factorization_including_zero_edges(self):
        # The last pattern starts with a zero incoming edge whose descendants
        # still have literal nonzero marker columns. No ratio or sign of an
        # incoming amplitude is used to define their orientation.
        patterns = (REGULAR, SINGULAR,
                    ((F(0), F(-1)), (F(-3, 5), F(4, 5)),
                     (F(-1), F(0)), (F(0), F(1))))
        for width in (1, 2, 3):
            for pattern in patterns:
                with self.subTest(width=width, pattern=pattern):
                    _, frame, _, _ = _geometry(width, pattern)
                    absolute = tuple((abs(cosine), abs(sine)) for cosine, sine in pattern)
                    _, magnitude_frame, _, _ = _geometry(width, absolute)
                    pairs = tuple(pattern[i % len(pattern)] for i in range((1 << width) - 1))
                    leaf, marker = _tree_signs(width, pairs)
                    np.testing.assert_array_equal(
                        frame, leaf[:, None] * magnitude_frame * marker[None, :])
                    np.testing.assert_array_equal(frame.T @ frame, np.eye(1 << width, dtype=int))
                    self.assertEqual(marker[0], 1)


if __name__ == '__main__':
    unittest.main()
