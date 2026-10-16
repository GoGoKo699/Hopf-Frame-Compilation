"""Finite exact partitions and full-frame flags for compact nonuniform coins.

Dense unitary completions are bounded algebraic fixtures, not native emitters.
"""
from __future__ import annotations

from itertools import product
import unittest

import numpy as np

from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import real_frame_matrix, real_tree_data


ATOL = 8e-11


def _compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _compositions(total - first, parts - 1):
            yield (first,) + rest


def _fixed_atoms(outcomes, denominator):
    remaining, weights = denominator, []
    while remaining:
        quotient = remaining // outcomes
        atom = 1 if not quotient else 1 << (quotient.bit_length() - 1)
        weights.append(atom)
        remaining -= atom
    return tuple(weights)


def _partition(weights, capacities):
    remaining = list(capacities)
    assignment = []
    for atom in weights:
        label = max(range(len(remaining)), key=lambda j: remaining[j])
        if atom > remaining[label]:
            raise AssertionError("fixed atom does not fit a remaining capacity")
        assignment.append(label)
        remaining[label] -= atom
    if any(remaining):
        raise AssertionError("partition left an unfilled capacity")
    return np.asarray(assignment, dtype=int)


def _two_outcome_columns(capacities):
    count, denominator = len(capacities), sum(capacities)
    residual = [count * value for value in capacities]
    active = set(range(count))
    acceptance = [denominator] * count
    parent = list(range(count))
    removed = []
    while True:
        small = sorted(j for j in active if residual[j] < denominator)
        if not small:
            break
        large = sorted(j for j in active if residual[j] >= denominator)
        j, target = small[0], large[0]
        acceptance[j], parent[j] = residual[j], target
        residual[target] -= denominator - residual[j]
        active.remove(j)
        removed.append(j)
    if any(residual[j] != denominator for j in active):
        raise AssertionError("surviving column has the wrong residual")
    return acceptance, parent, removed + sorted(active)


def _completion(target, phase):
    target = np.asarray(target, dtype=complex)
    anchor = np.zeros(target.size, dtype=complex)
    anchor[0] = 1
    direction = anchor - target
    if np.linalg.norm(direction) < 1e-13:
        word = np.eye(target.size, dtype=complex)
    else:
        direction /= np.linalg.norm(direction)
        word = np.eye(target.size) - 2 * np.outer(direction, direction.conj())
    # The first column stays literal, while all other input columns change.
    return word * np.exp(1j * phase * np.arange(target.size))[None, :]


def _block_fixture(numerators, blocks=(1, 2), tau=2, phase=0.0):
    q, radix = sum(blocks), 1 << tau
    boundaries = np.cumsum((0,) + blocks)
    weights, partitions, dimensions, transforms = [], [], [], []
    for block, height in enumerate(blocks):
        denominator = radix**height
        atoms = _fixed_atoms(1 << height, denominator)
        dimension = 1 << (len(atoms) - 1).bit_length()
        probabilities = np.zeros(dimension)
        probabilities[:len(atoms)] = np.asarray(atoms) / denominator
        weights.append(probabilities)
        dimensions.append(dimension)
        gamma = np.sqrt(probabilities)
        transforms.append(_completion(gamma, phase * (block + 1)).conj().T)
        by_prefix = {}
        start = int(boundaries[block])
        for prefix in range(1 << start):
            capacities = []
            for suffix in range(1 << height):
                node, mass = (1 << start) + prefix, 1
                for position in range(height):
                    bit = suffix >> (height - position - 1) & 1
                    numerator = numerators[node - 1]
                    mass *= radix - numerator if bit else numerator
                    node = 2 * node + bit
                capacities.append(mass)
            assignment = np.zeros(dimension, dtype=int)
            assignment[:len(atoms)] = _partition(atoms, capacities)
            by_prefix[prefix] = assignment
        partitions.append(by_prefix)

    coin_tuples = list(product(*(range(size) for size in dimensions)))
    decoded, total_weights = [], []
    for coins in coin_tuples:
        prefix, mass = 0, 1.0
        for block, (coin, height) in enumerate(zip(coins, blocks)):
            child = partitions[block][prefix][coin]
            prefix = (prefix << height) + int(child)
            mass *= weights[block][coin]
        decoded.append(prefix)
        total_weights.append(mass)
    decoded = np.asarray(decoded)
    total_weights = np.asarray(total_weights)
    leaf_mass = np.bincount(decoded, weights=total_weights, minlength=1 << q)
    embedding = np.zeros((len(coin_tuples), 1 << q))
    for row, leaf in enumerate(decoded):
        embedding[row, leaf] = np.sqrt(total_weights[row] / leaf_mass[leaf])
    transform = np.ones((1, 1), dtype=complex)
    for word in transforms:
        transform = np.kron(transform, word)
    last = np.asarray([max([i + 1 for i, value in enumerate(coins) if value] or [0])
                       for coins in coin_tuples])
    depths = np.zeros(1 << q, dtype=int)
    bands = np.zeros(1 << q, dtype=int)
    for node in range(1, 1 << q):
        depth = node.bit_length()
        column = marker_label(node, q)
        depths[column] = depth
        bands[column] = int(np.searchsorted(boundaries[1:], depth, side="left")) + 1
    angles = np.arccos(np.sqrt(np.asarray(numerators) / radix))
    return dict(q=q, boundaries=boundaries, decoded=decoded, weights=total_weights,
                embedding=embedding, transform=transform, last=last, bands=bands,
                angles=angles, frame=real_frame_matrix(angles), leaf_mass=leaf_mass)


class CompactCoinFlagTests(unittest.TestCase):
    def test_one_fixed_atom_list_handles_every_small_capacity_vector(self):
        for outcomes, denominator in ((2, 8), (3, 8), (4, 8), (4, 16)):
            weights = _fixed_atoms(outcomes, denominator)
            self.assertEqual(sum(weights), denominator)
            self.assertTrue(all(atom & (atom - 1) == 0 for atom in weights))
            for capacities in _compositions(denominator, outcomes):
                assignment = _partition(weights, capacities)
                actual = [sum(atom for atom, label in zip(weights, assignment) if label == x)
                          for x in range(outcomes)]
                self.assertEqual(actual, list(capacities))

    def test_two_outcome_columns_have_exact_mass_and_no_directed_cycle(self):
        for outcomes, denominator in ((2, 8), (3, 8), (4, 8)):
            binary_atoms = _fixed_atoms(2, denominator)
            for capacities in _compositions(denominator, outcomes):
                acceptance, parent, order = _two_outcome_columns(capacities)
                position = {label: rank for rank, label in enumerate(order)}
                mass = [0] * outcomes
                for column in range(outcomes):
                    assignment = _partition(binary_atoms,
                                            (acceptance[column], denominator - acceptance[column]))
                    for atom, rejected in zip(binary_atoms, assignment):
                        mass[parent[column] if rejected else column] += atom
                    if acceptance[column] != denominator:
                        self.assertLess(position[column], position[parent[column]])
                self.assertEqual(mass, [outcomes * value for value in capacities])

    def test_heavy_fiber_cannot_fit_a_short_label_preserving_rank(self):
        weights = _fixed_atoms(4, 16)
        self.assertLessEqual(max(weights), 4)
        assignment = _partition(weights, (13, 1, 1, 1))
        self.assertGreater(np.count_nonzero(assignment == 0), 2)
        # A one-bit conditional rank has only two basis states, regardless
        # of the following classical permutation.
        self.assertGreater(np.count_nonzero(assignment == 0), 1 << 1)

    def test_equal_weight_swap_preserves_root_but_not_coarse_flag(self):
        plus = np.ones(2) / np.sqrt(2)
        minus = np.array((1.0, -1.0)) / np.sqrt(2)
        hadamard = np.array(((1.0, 1.0), (1.0, -1.0))) / np.sqrt(2)
        transform = np.kron(np.kron(hadamard, hadamard),
                            np.kron(hadamard, hadamard))
        root = np.kron(np.kron(plus, plus), np.kron(plus, plus))
        top_marker = np.kron(np.kron(minus, plus), np.kron(plus, plus))
        # Coordinate order is (j1,z1,j2,z2). Every atom has weight 1/16,
        # and each block-weight tuple is (1/4,1/4), before and after swap.
        permutation = np.arange(16).reshape((2,) * 4).transpose(2, 1, 0, 3).ravel()
        np.testing.assert_allclose(root[permutation], root, atol=ATOL)
        original = transform @ top_marker
        swapped = transform @ top_marker[permutation]
        expected_original, expected_swapped = np.zeros(16), np.zeros(16)
        expected_original[8], expected_swapped[2] = 1, 1
        np.testing.assert_allclose(original, expected_original, atol=ATOL)
        np.testing.assert_allclose(swapped, expected_swapped, atol=ATOL)
        original_blocks = divmod(int(np.argmax(np.abs(original))), 4)
        swapped_blocks = divmod(int(np.argmax(np.abs(swapped))), 4)
        self.assertEqual(max(i + 1 for i, block in enumerate(original_blocks) if block), 1)
        self.assertEqual(max(i + 1 for i, block in enumerate(swapped_blocks) if block), 2)

    def test_partition_and_every_coarse_boundary_scaling_vector(self):
        for numerators in ((1, 2, 3, 1, 3, 2, 1), (3, 1, 1, 3, 1, 3, 3)):
            data = _block_fixture(numerators)
            embedding, q = data["embedding"], data["q"]
            np.testing.assert_allclose(embedding.T @ embedding, np.eye(1 << q), atol=ATOL)
            tree = real_tree_data(data["angles"])
            for depth in data["boundaries"]:
                for prefix in range(1 << int(depth)):
                    node = (1 << int(depth)) + prefix
                    mask = data["decoded"] >> (q - int(depth)) == prefix
                    expected = np.sqrt(data["weights"] * mask)
                    expected /= np.linalg.norm(expected)
                    np.testing.assert_allclose(embedding @ tree.subtree[node], expected, atol=ATOL)

    def test_full_coarse_flag_is_independent_of_coin_completion(self):
        for phase in (0.0, 0.47):
            for numerators in ((1, 2, 3, 1, 3, 2, 1), (3, 1, 1, 3, 1, 3, 3)):
                data = _block_fixture(numerators, phase=phase)
                transformed = data["transform"] @ data["embedding"] @ data["frame"]
                expected_root = np.zeros(len(transformed))
                expected_root[0] = 1
                np.testing.assert_allclose(transformed[:, 0], expected_root, atol=ATOL)
                for column, band in enumerate(data["bands"]):
                    np.testing.assert_allclose(transformed[data["last"] != band, column], 0,
                                               atol=ATOL)
                    self.assertAlmostEqual(np.linalg.norm(transformed[:, column]), 1, delta=ATOL)

    def test_full_tag_intertwining_and_return_include_references(self):
        data = _block_fixture((1, 2, 3, 1, 3, 2, 1), phase=0.61)
        combined = data["transform"] @ data["embedding"]
        tag_size, leaves = 4, 1 << data["q"]
        # The final axis represents all logical columns and may be entangled
        # with any external reference; it is not measured or sampled.
        source = np.zeros((leaves, tag_size, leaves), dtype=complex)
        expected = np.zeros_like(source)
        for column, band in enumerate(data["bands"]):
            source[:, band, column] = data["frame"][:, column]
            expected[:, 0, column] = data["frame"][:, column]
        coin = np.einsum("zx,xtc->ztc", combined, source)
        erased = np.empty_like(coin)
        for z, label in enumerate(data["last"]):
            for tag in range(tag_size):
                erased[z, tag ^ int(label)] = coin[z, tag]
        output = np.einsum("xz,ztc->xtc", combined.conj().T, erased)
        np.testing.assert_allclose(output, expected, atol=ATOL)
        reference = np.exp(0.31j * np.arange(leaves * 3)).reshape(leaves, 3)
        reference /= np.linalg.norm(reference)
        np.testing.assert_allclose(output @ reference, expected @ reference, atol=ATOL)


if __name__ == "__main__":
    unittest.main()
