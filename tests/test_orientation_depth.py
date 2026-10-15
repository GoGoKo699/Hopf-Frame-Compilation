"""Bounded full-column audits of affine and deterministic orientation words.

The prescribed frame is reconstructed from path products, independently of
the packed-tag circuits. Actual complex one-qubit Clifford+T words exercise
the affine construction; signed deterministic trees exercise the Fourier
source and last-deviation eraser. Query fixtures implement their reversible
truth-table actions, including dirty phase echoes. They are not scalable
emitters, and these checks imply no asymptotic resource bound.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_homogeneous_depth import (
        NATIVE_WORDS, Z, _data, _five_stage_word, _last_one, _tag_xor,
    )
    from .test_operator_source_compiler import _expand_toffolis, _word_matrix
except ImportError:
    from test_homogeneous_depth import (
        NATIVE_WORDS, Z, _data, _five_stage_word, _last_one, _tag_xor,
    )
    from test_operator_source_compiler import _expand_toffolis, _word_matrix


ATOL = 4e-11
X = np.array([[0, 1], [1, 0]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)


def _bits(word, width):
    return np.array([(word >> (width - position - 1)) & 1
                     for position in range(width)], dtype=int)


def _word(bits):
    return sum(int(bit) << (len(bits) - position - 1) for position, bit in enumerate(bits))


def _inverse_binary(matrix):
    width = len(matrix)
    result = np.column_stack((matrix.copy(), np.eye(width, dtype=int)))
    for column in range(width):
        pivot = next(row for row in range(column, width) if result[row, column])
        result[[column, pivot]] = result[[pivot, column]]
        for row in range(width):
            if row != column and result[row, column]:
                result[row] ^= result[column]
    return result[:, width:]


def _affine_permutation(matrix, offset):
    width = len(matrix)
    result = np.zeros((1 << width, 1 << width))
    for old in range(len(result)):
        new = _word((matrix @ _bits(old, width) + offset) % 2)
        result[new, old] = 1
    return result


def _frame_columns(nodes, width):
    """Literal prefix, marker column, and descendant path amplitudes."""
    size = 1 << width
    result = np.zeros((size, size), dtype=complex)
    for column in range(size):
        tag = _last_one(column, width)
        prefix_length = max(0, tag - 1)
        prefix = column >> (width - prefix_length)
        for leaf in range(size):
            if leaf >> (width - prefix_length) != prefix:
                continue
            amplitude = 1
            for depth in range(prefix_length, width):
                node = (1 << depth) - 1 + (leaf >> (width - depth))
                branch = (leaf >> (width - depth - 1)) & 1
                source_bit = int(tag > 0 and depth == tag - 1)
                amplitude *= nodes[node][branch, source_bit]
            result[leaf, column] = amplitude
    return result


def _initial(width, *, occupied=False):
    size, labels, dirty = 1 << width, 1 << width.bit_length(), 1 << (width + 1)
    if occupied:
        dimension = size * labels * dirty
        return np.eye(dimension, dtype=complex).reshape(size, labels, dirty, dimension)
    result = np.zeros((size, labels, dirty, size * dirty), dtype=complex)
    for word in range(size):
        for old in range(dirty):
            result[word, 0, old, word * dirty + old] = 1
    return result


def _mask_query(state, masks):
    result = np.empty_like(state)
    for tag, mask in enumerate(masks):
        for dirty in range(state.shape[2]):
            result[:, tag, dirty ^ mask] = state[:, tag, dirty]
    return result


def _controlled_bank(state, matrix):
    width = (state.shape[0] - 1).bit_length()
    result = state.copy()
    for position in range(width):
        bit = 1 << (width - position - 1)
        for low in range(state.shape[0]):
            if low & bit:
                continue
            for dirty in range(state.shape[2]):
                if not (dirty >> position) & 1:
                    continue
                pair = np.stack((result[low, :, dirty].copy(),
                                 result[low | bit, :, dirty].copy()))
                transformed = np.tensordot(matrix, pair, axes=(1, 0))
                result[low, :, dirty], result[low | bit, :, dirty] = transformed
    return result


def _selected_bank(state, matrix, mask, *, active=True):
    masks = [mask(tag) if active else 0 for tag in range(state.shape[1])]
    for _ in range(2):
        state = _mask_query(state, masks)
        state = _controlled_bank(state, matrix)
    return state


def _tail_mask(tag, width):
    return sum(int(position + 1 >= tag) << position for position in range(width))


def _phase_echo(state, phase, *, active):
    """Q_f, Z_dirty, Q_f^dagger, Z_dirty, retaining arbitrary dirty input."""
    width = (state.shape[0] - 1).bit_length()
    phase_bit = 1 << width
    for _ in range(2):
        result = np.empty_like(state)
        for word in range(state.shape[0]):
            for tag in range(state.shape[1]):
                bit = phase(word, tag) if active else 0
                for dirty in range(state.shape[2]):
                    result[word, tag, dirty ^ (phase_bit * bit)] = state[word, tag, dirty]
        state = result
        for dirty in range(state.shape[2]):
            if dirty & phase_bit:
                state[:, :, dirty] *= -1
    return state


def _affine_tail(state, inverse, offset, *, active):
    width = len(inverse)
    forward = _inverse_binary(inverse)
    translation = forward @ offset % 2

    def diagonal(word, tag):
        z = _bits(word, width)
        if tag == 0:
            return int(z @ offset % 2)
        if tag > width:
            return 0
        prefix = tag - 1
        phase = 0
        for target in range(prefix, width):
            coefficients = inverse[target, :prefix] @ forward[:prefix, :prefix] % 2
            constant = (inverse[target, :prefix] @ translation[:prefix]
                        + inverse[target, prefix] + offset[target] + int(target == prefix)) % 2
            phase ^= int(z[target] * ((coefficients @ z[:prefix] + constant) % 2))
        return phase

    # The tail Hadamards are deliberately unconditional with respect to E.
    # When E=0 their actual complete dirty-bank words cancel on all tags.
    state = _selected_bank(state, H, lambda tag: _tail_mask(tag, width))
    state = _phase_echo(state, diagonal, active=active)
    return _selected_bank(state, H, lambda tag: _tail_mask(tag, width))


def _affine_q(state, inverse, offset, *, active):
    width = len(inverse)
    forward = _inverse_binary(inverse)
    translation = forward @ offset % 2
    state = _tag_xor(state, active=active)
    if active:
        state = _data(_affine_permutation(inverse, offset), state)
    state = _affine_tail(state, inverse, offset, active=active)

    def sign(word, tag):
        if not 1 <= tag <= width:
            return 0
        prefix = tag - 1
        return int((forward[prefix, :prefix] @ _bits(word, width)[:prefix]
                    + translation[prefix]) % 2)

    state = _phase_echo(state, sign, active=active)
    return _tag_xor(state, active=active)


def _q_reference(inverse, offset):
    width = len(inverse)
    result = np.zeros((1 << width, 1 << width))
    result[0, 0] = 1
    for word in range(1, 1 << width):
        tag = _last_one(word, width)
        prefix = tag - 1
        old = _bits(word, width)[:prefix]
        new = (inverse[:prefix, :prefix] @ old + offset[:prefix]) % 2
        target = (_word(new) << (width - prefix)) | (1 << (width - tag))
        orientation = (inverse[prefix, :prefix] @ old + offset[prefix]) % 2
        result[target, word] = (-1) ** orientation
    return result


def _oriented_nodes(unitaries, inverse, offset):
    nodes = []
    for depth, unitary in enumerate(unitaries):
        for prefix in range(1 << depth):
            orientation = (inverse[depth, :depth] @ _bits(prefix, depth) + offset[depth]) % 2
            nodes.append(X @ unitary @ Z if orientation else unitary)
    return nodes


def _corner_table(nodes, width, tag, prefix):
    """Classical continuation and literal sign, without constructing a frame."""
    leaf, sign = prefix, 1
    start = max(0, tag - 1)
    for depth in range(start, width):
        rotation = nodes[(1 << depth) - 1 + leaf]
        favored = int(rotation[1, 0] != 0)
        scalar = int(rotation[favored, 0].real)
        marker = int(tag > 0 and depth == tag - 1)
        branch = favored ^ marker
        sign *= scalar * ((-1) ** (favored * marker))
        leaf = 2 * leaf + branch
    return leaf, sign


def _last_deviation(nodes, width, word):
    result, prefix = 0, 0
    for depth in range(width):
        rotation = nodes[(1 << depth) - 1 + prefix]
        favored = int(rotation[1, 0] != 0)
        bit = (word >> (width - depth - 1)) & 1
        if bit != favored:
            result = depth + 1
        prefix = 2 * prefix + bit
    return result


def _corner_word(state, nodes, *, active):
    width = (state.shape[0] - 1).bit_length()
    state = _tag_xor(state, active=active)
    state = _selected_bank(state, X,
                           lambda tag: (1 << (tag - 1)) if 1 <= tag <= width else 0,
                           active=active)
    state = _selected_bank(state, H, lambda tag: _tail_mask(tag, width), active=active)

    def source_phase(word, tag):
        if tag > width:
            return 0
        prefix_length = max(0, tag - 1)
        suffix = width - prefix_length
        prefix = word >> suffix
        output, sign = _corner_table(nodes, width, tag, prefix)
        mask = (1 << suffix) - 1
        return int(sign == -1) ^ ((output & word & mask).bit_count() & 1)

    state = _phase_echo(state, source_phase, active=active)
    state = _selected_bank(state, H, lambda tag: _tail_mask(tag, width), active=active)
    result = np.empty_like(state)
    for word in range(state.shape[0]):
        depth = _last_deviation(nodes, width, word) if active else 0
        for tag in range(state.shape[1]):
            result[word, tag ^ depth] = state[word, tag]
    return result


class OrientationDepthTests(unittest.TestCase):
    def test_affine_tail_marker_sign_root_and_outside_prefix(self):
        inverse = np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1]])
        outside = np.array([[1, 0], [0, 1], [1, 1]])
        initial = _initial(3)
        for external in range(4):
            offset = (np.array([1, 0, 1]) + outside @ _bits(external, 2)) % 2
            actual = _affine_q(initial, inverse, offset, active=True)
            expected = _data(_q_reference(inverse, offset), initial)
            np.testing.assert_allclose(actual, expected, atol=ATOL)
        occupied = _initial(2, occupied=True)
        # h=2 has a padded tag 3; disabled phases must cancel there too.
        np.testing.assert_allclose(
            _affine_q(occupied, inverse[:2, :2], np.array([1, 1]), active=False),
            occupied, atol=ATOL)

    def test_complete_affine_stage_with_actual_complex_native_core(self):
        width = 3
        inverse = np.array([[1, 0, 0], [1, 1, 0], [0, 1, 1]])
        forward = _inverse_binary(inverse)
        conjugators = [_word_matrix(1, word) for word in NATIVE_WORDS]
        unitaries = [v @ Z @ v.conj().T @ Z for v in conjugators]
        self.assertGreater(max(np.max(np.abs(u.imag)) for u in unitaries), 0.1)
        for offset in (np.array([1, 0, 1]), np.array([0, 1, 1])):
            initial = _initial(width)
            state = _affine_q(initial, inverse, offset, active=True)
            state = _five_stage_word(state, unitaries, active=True, conjugators=conjugators)
            state = _data(_affine_permutation(forward, forward @ offset % 2), state)
            reference = _frame_columns(_oriented_nodes(unitaries, inverse, offset), width)
            np.testing.assert_allclose(state, _data(reference, initial), atol=ATOL)
            np.testing.assert_allclose(reference.conj().T @ reference, np.eye(8), atol=ATOL)
        occupied = _initial(2, occupied=True)
        disabled = _affine_q(occupied, inverse[:2, :2], np.array([1, 0]), active=False)
        disabled = _five_stage_word(disabled, unitaries[:2], active=False,
                                    conjugators=conjugators[:2])
        np.testing.assert_allclose(disabled, occupied, atol=ATOL)

    def test_unipotent_involutions_and_literal_controlled_basis_words(self):
        jordan = np.array([[1, 0, 0], [1, 1, 0], [0, 1, 1]])
        substitution = np.array([[1, 0, 0], [0, 1, 0], [0, 1, 1]])
        second = substitution @ jordan % 2
        identity = np.eye(3, dtype=int)
        np.testing.assert_array_equal(substitution @ substitution % 2, identity)
        np.testing.assert_array_equal(substitution @ jordan @ substitution % 2,
                                      _inverse_binary(jordan))
        np.testing.assert_array_equal(second @ second % 2, identity)
        bases = (np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]]),
                 np.array([[0, 1, 0], [1, 0, 1], [1, 0, 0]]))
        canonical = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
        controlled = []
        native_pair = _word_matrix(4, _expand_toffolis([('CCX', 3, 1, 2)]))
        for involution, basis in zip((substitution, second), bases):
            np.testing.assert_array_equal(basis @ canonical @ _inverse_binary(basis) % 2,
                                          involution)
            conjugator = np.kron(np.eye(2), _affine_permutation(basis, np.zeros(3, dtype=int)))
            controlled.append(conjugator @ native_pair @ conjugator.T)
        actual = controlled[0] @ controlled[1]
        expected = np.zeros((16, 16))
        expected[:8, :8] = np.eye(8)
        expected[8:, 8:] = _affine_permutation(jordan, np.zeros(3, dtype=int))
        np.testing.assert_allclose(actual, expected, atol=ATOL)
        translation = _word_matrix(4, [('CX', 3, 2), ('CX', 3, 0)])
        expected[8:, 8:] = _affine_permutation(jordan, np.array([1, 0, 1]))
        np.testing.assert_allclose(translation @ actual, expected, atol=ATOL)

    def test_signed_deterministic_source_and_last_deviation_return(self):
        corners = ((1, 0), (0, 1), (-1, 0), (0, -1))
        patterns = (corners,
                    ((0, -1), (-1, 0), (0, 1), (1, 0), (0, 1), (0, -1), (0, 1)))
        for width in (1, 2, 3):
            for pattern in patterns:
                pairs = [pattern[index % len(pattern)] for index in range((1 << width) - 1)]
                nodes = [np.array([[c, -s], [s, c]], dtype=complex) for c, s in pairs]
                reference = _frame_columns(nodes, width)
                np.testing.assert_array_equal(reference.conj().T @ reference, np.eye(1 << width))
                initial = _initial(width)
                actual = _corner_word(initial, nodes, active=True)
                np.testing.assert_allclose(actual, _data(reference, initial), atol=ATOL)
                for column in range(1 << width):
                    leaf = int(np.flatnonzero(reference[:, column])[0])
                    self.assertEqual(_last_deviation(nodes, width, leaf), _last_one(column, width))
        occupied = _initial(2, occupied=True)
        nodes = [np.array([[c, -s], [s, c]], dtype=complex) for c, s in corners[:3]]
        np.testing.assert_allclose(_corner_word(occupied, nodes, active=False), occupied, atol=ATOL)

    def test_simultaneous_dirty_triangular_update_has_changed_address_counterexample(self):
        # x1,x2 occupy bits 1,0; dirty d1,d2 occupy bits 3,2.
        query = [('X', 3), ('CX', 1, 2)]
        update = [('CX', 3, 1), ('CX', 2, 0)]
        actual = _word_matrix(4, query + update + list(reversed(query)) + update)
        self.assertEqual(int(np.argmax(np.abs(actual[:, 0]))), 7)  # x=11, d=01.
        self.assertEqual(actual[7, 0], 1)
        self.assertEqual(actual[2, 0], 0)  # Intended x=10, returned d=00.


if __name__ == '__main__':
    unittest.main()
