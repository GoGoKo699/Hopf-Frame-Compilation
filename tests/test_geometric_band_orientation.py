"""Small coherent audits of geometric-band orientation recovery.

Frames are reconstructed from literal path products, independently of the
recovery word. These tests check its ideal algebra and initialized-tail
premises, not a scalable native emitter or asymptotic resource theorem.
"""
from __future__ import annotations

import math
import unittest

import numpy as np

from tests.test_orientation_depth import _frame_columns
from tests.test_subtree_projector_decoder import _search as _spectral_search
from tests.test_subtree_projector_decoder import _xor_tag as _xor_spectral_tag


ATOL = 4e-13


def _boundaries(height, reserve):
    boundaries = [height]
    remaining, capacity = height, reserve
    while remaining:
        width = min(remaining, capacity // 2)
        remaining -= width
        capacity += width
        boundaries.append(remaining)
    return tuple(reversed(boundaries))


def _last_one(word, width):
    return width - ((word & -word).bit_length() - 1) if word else 0


def _orientation(depth, prefix):
    if depth == 0:
        return 0
    if depth == 1:
        return prefix & 1
    if depth == 2:
        return int(prefix == 3)  # r_3(x_1,x_2)=x_1*x_2, genuinely nonlinear.
    return ((prefix >> (depth - 1)) & 1) ^ int(prefix & 3 == 3)


def _rotation(c, s):
    return np.array([[c, -s], [s, c]], dtype=complex)


def _nodes(pairs, signed=False):
    result, signs = [], []
    pattern = ((1, -1), (-1, 1), (-1, -1), (1, 1))
    for depth, (c, s) in enumerate(pairs):
        for prefix in range(1 << depth):
            a, b = (s, c) if _orientation(depth, prefix) else (c, s)
            sign = pattern[len(result) % len(pattern)] if signed else (1, 1)
            result.append(_rotation(sign[0] * a, sign[1] * b))
            signs.append(sign)
    return result, signs


def _apply_level(state, width, depth, matrices, inverse=False):
    result = state.copy()
    suffix_width = width - depth - 1
    for prefix, matrix in enumerate(matrices):
        matrix = matrix.conj().T if inverse else matrix
        for suffix in range(1 << suffix_width):
            low = (prefix << (suffix_width + 1)) | suffix
            high = low | (1 << suffix_width)
            pair = np.stack((state[low], state[high]))
            result[low], result[high] = np.tensordot(matrix, pair, axes=(1, 0))
    return result


def _descendant(state, nodes, width, boundary, inverse=False):
    depths = range(boundary, width)
    for depth in reversed(depths) if inverse else depths:
        start = (1 << depth) - 1
        state = _apply_level(state, width, depth,
                             nodes[start:start + (1 << depth)], inverse)
    return state


def _local_map(word, prefix, start, stop, inverse=False):
    result, physical_prefix = 0, prefix
    for depth in range(start, stop):
        bit = (word >> (stop - depth - 1)) & 1
        output = bit ^ _orientation(depth, physical_prefix)
        result = 2 * result + output
        physical_prefix = 2 * physical_prefix + (bit if inverse else output)
    return result


def _replace_block(word, value, width, start, stop):
    shift, mask = width - stop, (1 << (stop - start)) - 1
    return (word & ~(mask << shift)) | (value << shift)


def _local_direct(state, pairs, width, start, stop, inverse=False):
    def permutation(value):
        result = np.empty_like(value)
        mask = (1 << (stop - start)) - 1
        for word in range(1 << width):
            block = (word >> (width - stop)) & mask
            target = _local_map(block, word >> (width - start), start, stop, inverse)
            result[_replace_block(word, target, width, start, stop)] = value[word]
        return result

    if inverse:
        state = permutation(state)
    for depth in range(start, stop):
        state = _apply_level(state, width, depth,
                             [_rotation(*pairs[depth])] * (1 << depth), inverse)
    return state if inverse else permutation(state)


def _borrowed_basis(registers, width, start, stop, donors, inverse=False, active=True):
    """Literal save/query/swap/inverse-query/restore on physical X,Z,Y."""
    registers = list(registers)
    mask = (1 << (stop - start)) - 1
    steps = ("save", "forward", "swap", "backward", "save")
    for step in reversed(steps) if inverse else steps:
        if not active:
            continue
        if step == "save":
            for position, (register, bit) in enumerate(donors):
                difference = ((registers[register] >> bit) ^ (registers[2] >> position)) & 1
                registers[register] ^= difference << bit
                registers[2] ^= difference << position
        else:
            block = (registers[0] >> (width - stop)) & mask
            if step == "swap":
                registers[0] = _replace_block(registers[0], registers[2] & mask,
                                              width, start, stop)
                registers[2] = (registers[2] & ~mask) | block
            else:
                # The physical tail may now contain dirty Y. Neither table
                # reads it: both use precisely the current block and prefix.
                registers[2] ^= _local_map(
                    block, registers[0] >> (width - start), start, stop,
                    inverse=step == "backward")
    return tuple(registers)


def _borrowed_action(state, width, start, stop, donors, inverse=False, active=True):
    result = np.empty_like(state)
    for registers in np.ndindex(state.shape[:3]):
        target = _borrowed_basis(registers, width, start, stop, donors, inverse, active)
        result[target] = state[registers]
    return result


def _selected_uniform(state, pairs, width, start, stop, inverse=False, active=True):
    """Actual V Q V-dagger Q word, with nontrivial scalar phases in V."""
    def phase(value):
        result = value.copy()
        mask = ((1 << (stop - start)) - 1) << (width - stop)
        for word in range(1 << width):
            if active and (word & mask).bit_count() % 2:
                result[word] *= -1
        return result

    def half_angles(value, dagger):
        for depth, (c, s) in enumerate(pairs):
            angle = math.atan2(s, c) / 2
            matrix = np.exp(1j * (.13 + .07 * depth)) * _rotation(math.cos(angle), math.sin(angle))
            value = _apply_level(value, width, depth, [matrix] * (1 << depth), dagger)
        return value

    if inverse:
        return phase(half_angles(phase(half_angles(state, True)), False))
    return half_angles(phase(half_angles(phase(state), True)), False)


def _donors(width, stop, count):
    available = [(0, width - position - 1) for position in range(stop, width)]
    available += [(1, position) for position in range(count)]
    return available[:count]


def _bin(depth, boundaries):
    return next((index for index, stop in enumerate(boundaries[1:]) if depth <= stop), 0)


def _fine_recovery(state, pairs, nodes, boundaries, active=True):
    """Axes are X, Zfree, dirty Y, fine T, coarse R, reference."""
    width = len(pairs)
    computed = np.empty_like(state)
    for tag in range(state.shape[3]):
        offset = _bin(tag, boundaries) if active else 0
        for coarse in range(state.shape[4]):
            computed[:, :, :, tag, coarse ^ offset] = state[:, :, :, tag, coarse]
    for coarse, (start, stop) in enumerate(zip(boundaries, boundaries[1:])):
        branch = computed[:, :, :, :, coarse].copy()
        donors = _donors(width, stop, stop - start)
        if active:
            branch = _descendant(branch, nodes, width, stop, inverse=True)
        branch = _borrowed_action(branch, width, start, stop, donors, inverse=True, active=active)
        branch = _selected_uniform(branch, pairs, width, start, stop, inverse=True, active=active)
        tagged = np.empty_like(branch)
        for word in range(1 << width):
            block = (word >> (width - stop)) & ((1 << (stop - start)) - 1)
            depth = start + _last_one(block, stop - start) if active else 0
            for old in range(state.shape[3]):
                tagged[word, :, :, old ^ depth] = branch[word, :, :, old]
        branch = _selected_uniform(tagged, pairs, width, start, stop, active=active)
        branch = _borrowed_action(branch, width, start, stop, donors, active=active)
        if active:
            branch = _descendant(branch, nodes, width, stop)
        computed[:, :, :, :, coarse] = branch
    return computed


class GeometricBandOrientationTests(unittest.TestCase):
    PAIRS = (((4 / 5, 3 / 5), (3 / 5, 4 / 5), (5 / 13, 12 / 13), (8 / 17, 15 / 17)),
             ((1, 0), (3 / 5, 4 / 5), (0, 1), (5 / 13, 12 / 13)))

    def test_integer_boundaries_have_slack_and_logarithmically_many_bins(self):
        for height in (0, 1, 2, 3, 4, 5, 7, 16, 65, 513, 10000):
            for reserve in (2, 3, 4, 6, 7, 32, 1030):
                boundaries = _boundaries(height, reserve)
                self.assertEqual((boundaries[0], boundaries[-1]), (0, height))
                for start, stop in zip(boundaries, boundaries[1:]):
                    self.assertGreater(stop, start)
                    self.assertLessEqual(2 * (stop - start), height - stop + reserve)
                # Every nonfinal backward step multiplies capacity by >=4/3.
                bound = 1 + math.ceil(math.log1p(height / reserve) / math.log(4 / 3))
                self.assertLessEqual(len(boundaries) - 1, bound)

    def test_all_truncated_flags_and_local_bands_include_nonlinear_orientation(self):
        width = 4
        self.assertEqual([_orientation(2, prefix) for prefix in range(4)], [0, 0, 0, 1])
        for pairs in self.PAIRS:
            nodes, _ = _nodes(pairs)
            frame = _frame_columns(nodes, width)
            np.testing.assert_allclose(frame.conj().T @ frame, np.eye(1 << width), atol=ATOL, rtol=0)
            for stop in range(1, width + 1):
                exposed = _descendant(frame, nodes, width, stop, inverse=True)
                truncated = _frame_columns(nodes[:(1 << stop) - 1], stop)
                for column in range(1 << width):
                    depth = _last_one(column, width)
                    if depth > stop:
                        continue
                    expected = np.zeros(1 << width, dtype=complex)
                    expected[::1 << (width - stop)] = truncated[:, column >> (width - stop)]
                    np.testing.assert_allclose(exposed[:, column], expected, atol=ATOL, rtol=0)
                for start in range(stop):
                    transported = _local_direct(exposed, pairs, width, start, stop, inverse=True)
                    for column in range(1 << width):
                        depth = _last_one(column, width)
                        if not (depth <= stop and (start == 0 or depth > start)):
                            continue
                        allowed = []
                        for word in range(1 << width):
                            local = (word >> (width - stop)) & ((1 << (stop - start)) - 1)
                            allowed.append((word & ((1 << (width - stop)) - 1)) == 0
                                           and start + _last_one(local, stop - start) == depth
                                           and word >> (width - start) == column >> (width - start))
                        np.testing.assert_allclose(transported[np.logical_not(allowed), column],
                                                   0, atol=ATOL, rtol=0)

    def test_physical_dirty_buffer_word_returns_every_dirty_basis_input(self):
        width, size = 4, 8
        for start, stop in ((0, 1), (1, 4), (0, 3)):
            donors = _donors(width, stop, stop - start)
            for registers in np.ndindex(1 << width, size, size):
                target = _borrowed_basis(registers, width, start, stop, donors)
                self.assertEqual(target[2], registers[2])
                self.assertEqual(_borrowed_basis(target, width, start, stop, donors, inverse=True), registers)
                self.assertEqual(_borrowed_basis(registers, width, start, stop, donors, active=False), registers)
                if any((registers[register] >> bit) & 1 for register, bit in donors):
                    continue
                word, protected, dirty = registers
                block = (word >> (width - stop)) & ((1 << (stop - start)) - 1)
                output = _local_map(block, word >> (width - start), start, stop)
                self.assertEqual(target, (_replace_block(word, output, width, start, stop), protected, dirty))

    def test_fine_clear_reload_preserves_cross_bin_coherence_and_literal_phases(self):
        width, dirty_size = 4, 8
        boundaries = _boundaries(width, 6)
        self.assertEqual(boundaries, (0, 1, 4))
        coefficients = np.array([column + 1 + 1j * (-1) ** column for column in range(1 << width)])
        coefficients /= np.linalg.norm(coefficients)
        for pairs in self.PAIRS:
            positive_nodes, _ = _nodes(pairs)
            positive = _frame_columns(positive_nodes, width)
            for signed in (False, True):
                nodes, signs = _nodes(pairs, signed)
                frame = _frame_columns(nodes, width)
                leaf_phase = np.exp(1j * np.arange(1 << width) ** 2 / 19) if signed else np.ones(1 << width)
                wrapper = leaf_phase.astype(complex).copy()
                marker_sign = np.ones(1 << width)
                for word in range(1 << width):
                    for depth in range(width):
                        prefix = word >> (width - depth)
                        wrapper[word] *= signs[(1 << depth) - 1 + prefix][(word >> (width - depth - 1)) & 1]
                    tag = _last_one(word, width)
                    if tag:
                        for depth in range(tag - 1):
                            prefix = word >> (width - depth)
                            marker_sign[word] *= signs[(1 << depth) - 1 + prefix][(word >> (width - depth - 1)) & 1]
                        node = (1 << (tag - 1)) - 1 + (word >> (width - tag + 1))
                        marker_sign[word] *= signs[node][0] * signs[node][1]
                frame = leaf_phase[:, None] * frame
                np.testing.assert_allclose(frame, wrapper[:, None] * positive * marker_sign,
                                           atol=ATOL, rtol=0)
                shape = (1 << width, dirty_size, dirty_size, 8, 2, dirty_size)
                source, expected = np.zeros(shape, complex), np.zeros(shape, complex)
                for column, coefficient in enumerate(coefficients):
                    depth = _last_one(column, width)
                    for dirty in range(dirty_size):
                        amplitude = coefficient * frame[:, column] / math.sqrt(dirty_size)
                        source[:, 0, dirty, depth, 0, dirty] += amplitude
                        expected[:, 0, dirty, 0, _bin(depth, boundaries), dirty] += amplitude
                unwrapped = source * wrapper.conj().reshape((-1,) + (1,) * 5)
                actual = _fine_recovery(unwrapped, pairs, positive_nodes, boundaries)
                actual *= wrapper.reshape((-1,) + (1,) * 5)
                np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)

    def test_coarse_search_xor_reverse_clears_five_bins_with_padded_thresholds(self):
        pairs = self.PAIRS[0] + ((20 / 29, 21 / 29),)
        width, bins, labels = len(pairs), 5, 8
        boundaries = tuple(range(width + 1))
        self.assertEqual([_bin(depth, boundaries) for depth in range(width + 1)],
                         [0, 0, 1, 2, 3, 4])
        nodes, _ = _nodes(pairs)
        frame = _frame_columns(nodes, width)
        size = 1 << width
        # R78's search prepares a binary rank threshold m, not a physical
        # depth. Its controlled blocks include the complement-bit sign.
        reflections = np.empty((2, labels, size, size), dtype=complex)
        for rank_threshold in range(labels):
            endpoint = boundaries[min(rank_threshold + 1, bins)]
            loader = _descendant(np.eye(size, dtype=complex), nodes, width, endpoint)
            scaling = loader[:, ::1 << (width - endpoint)]
            reflections[0, rank_threshold] = np.eye(size)
            reflections[1, rank_threshold] = 2 * scaling @ scaling.conj().T - np.eye(size)

        # Axes match the retained physical search: X, source W, fresh Y,
        # threshold D, old coarse R, and independent input-column batch.
        shape = (size, 1, labels, labels, labels, size)
        source, extracted, expected = (np.zeros(shape, complex) for _ in range(3))
        for column in range(size):
            rank = _bin(_last_one(column, width), boundaries)
            source[:, 0, 0, 0, rank, column] = frame[:, column]
            extracted[:, 0, rank, 0, rank, column] = frame[:, column]
            expected[:, 0, 0, 0, 0, column] = frame[:, column]
        searched = _spectral_search(source, reflections)
        np.testing.assert_allclose(searched, extracted, atol=ATOL, rtol=0)
        actual = _spectral_search(_xor_spectral_tag(searched, active=True),
                                  reflections, inverse=True)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)

        # Rank four reaches the virtual threshold m=5. Its endpoint must
        # saturate at h; treating padding as an empty projector changes a bit.
        wrong_padding = reflections.copy()
        wrong_padding[1, bins:] = -np.eye(size)
        self.assertGreater(np.linalg.norm(_spectral_search(source, wrong_padding) - extracted), 1)

    def test_inactive_fine_word_is_identity_with_occupied_protected_records(self):
        pairs = self.PAIRS[0]
        nodes, _ = _nodes(pairs)
        shape = (16, 8, 8, 8, 2, 1)
        indices = np.arange(math.prod(shape)).reshape(shape)
        state = (indices % 13 - 6 + 1j * (indices % 17 - 8)).astype(complex)
        state /= np.linalg.norm(state)
        actual = _fine_recovery(state, pairs, nodes, (0, 1, 4), active=False)
        np.testing.assert_allclose(actual, state, atol=ATOL, rtol=0)


if __name__ == "__main__":
    unittest.main()
