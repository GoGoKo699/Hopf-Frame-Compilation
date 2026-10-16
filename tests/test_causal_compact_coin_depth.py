"""Bounded causal compact-coin words, not native/asymptotic emitters.

The finite matching catalogues below are explicit unions of the shifts used
by these fixtures. They check the batched word and its causal dependencies;
the polynomial catalogue existence proof is a separate mathematical lemma.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product
import random
import unittest

import numpy as np

from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import real_frame_matrix
from tests.test_compact_coin_flags import (
    ATOL, _completion, _fixed_atoms, _partition, _two_outcome_columns,
)


def _parity(value):
    return value.bit_count() & 1


def _shear(seed, high, low_width, high_width):
    mask = (1 << high_width) - 1
    return sum(_parity((seed >> (row * high_width) & mask) & high) << row
               for row in range(low_width))


def _capacities(numerators, height, radix=4):
    result = []
    for leaf in range(1 << height):
        node, mass = 1, 1
        for position in range(height):
            bit = leaf >> (height - position - 1) & 1
            numerator = numerators[node - 1]
            mass *= radix - numerator if bit else numerator
            node = 2 * node + bit
        result.append(mass)
    return result


class _LocalCoins:
    def __init__(self, height, families):
        self.height, self.outcomes = height, 1 << height
        self.denominator = 4**height
        self.atoms = _fixed_atoms(2, self.denominator)
        self.gamma_size = 1 << (len(self.atoms) - 1).bit_length()
        self.core_size = self.outcomes * self.gamma_size
        self.capacities = {u: _capacities(nums, height) for u, nums in families.items()}
        self.fibers, self.metadata, self.code_inverse = {}, {}, {}
        self.rejected_collision_candidates = 0
        for prefix, capacities in self.capacities.items():
            acceptance, parent, _ = _two_outcome_columns(capacities)
            fibers = [[] for _ in capacities]
            for column in range(self.outcomes):
                rejected = _partition(self.atoms,
                                      (acceptance[column], self.denominator - acceptance[column]))
                for coin, reject in enumerate(rejected):
                    fibers[parent[column] if reject else column].append((column, coin))
            self.fibers[prefix] = fibers
            groups = defaultdict(list)
            for child, fiber in enumerate(fibers):
                groups[(len(fiber) - 1).bit_length()].append(child)
            for group, children in groups.items():
                low_width = min(height, (4 * len(children) - 1).bit_length())
                remaining, stage = set(children), 0
                while remaining:
                    if low_width == height:
                        seed, assigned = 0, sorted(remaining)
                    else:
                        for seed in range(1 << (low_width * (height - low_width))):
                            hashes = {x: (x & ((1 << low_width) - 1)) ^
                                      _shear(seed, x >> low_width, low_width, height - low_width)
                                      for x in remaining}
                            counts = Counter(hashes.values())
                            assigned = sorted(x for x in remaining if counts[hashes[x]] == 1)
                            if len(assigned) * 2 >= len(remaining):
                                break
                            self.rejected_collision_candidates += 1
                        else:
                            raise AssertionError("no finite staged hash isolates half the set")
                    for child in assigned:
                        self.metadata[prefix, child] = (group, low_width, stage, seed)
                    remaining.difference_update(assigned)
                    stage += 1
        maxima = [max(meta[j] for meta in self.metadata.values()) for j in range(4)]
        self.field_widths = [max(1, value.bit_length()) for value in maxima]
        self.payload_width = max(self.core_size.bit_length() - 1, height,
                                 max(group + low for group, low, _, _ in self.metadata.values()))
        self.width = self.payload_width + sum(self.field_widths)
        self.permutations, self.inverses, self.membership = {}, {}, {}
        for prefix in families:
            classes = defaultdict(list)
            for child, fiber in enumerate(self.fibers[prefix]):
                meta = self.metadata[prefix, child]
                group, low_width, stage, seed = meta
                low = (child & ((1 << low_width) - 1)) ^ \
                    _shear(seed, child >> low_width, low_width, height - low_width)
                key = (prefix, meta, low)
                if key in self.code_inverse:
                    raise AssertionError("two children have the same staged code")
                self.code_inverse[key] = child
                for rank, (_, coin) in enumerate(fiber):
                    classes[coin].append(self.atom_code(prefix, child, rank))
            permutation, membership = {}, {}
            for coin, sources in classes.items():
                targets = [column * self.gamma_size + coin for column in range(self.outcomes)]
                random.Random(101 * prefix + 17 * coin + height).shuffle(targets)
                for source, target in zip(sorted(sources), targets):
                    permutation[source] = target
                    membership[target] = self.decode_custom(prefix, source)
            # Complete the finite partial bijection into a total permutation:
            # close every path from source-only to target-only coordinates.
            source_set, target_set = set(permutation), set(permutation.values())
            for start in source_set - target_set:
                endpoint = permutation[start]
                while endpoint in source_set:
                    endpoint = permutation[endpoint]
                permutation[endpoint] = start
            if len(set(permutation.values())) != len(permutation):
                raise AssertionError("local permutation completion collided")
            self.permutations[prefix] = permutation
            self.inverses[prefix] = {value: key for key, value in permutation.items()}
            self.membership[prefix] = membership
        # The same catalogue is used by every prefix family at this block.
        shifts = {0}
        for permutation in self.permutations.values():
            shifts.update(source ^ target for source, target in permutation.items())
        self.catalogue = tuple(sorted(shifts))
        self.catalogue_index = {shift: index for index, shift in enumerate(self.catalogue)}
        gamma = np.zeros(self.gamma_size)
        gamma[:len(self.atoms)] = np.sqrt(np.asarray(self.atoms) / self.denominator)
        self.gamma_transform = _completion(gamma, 0.23 + height / 19).conj().T
        self.hadamard = np.asarray([[(-1.0)**_parity(row & column)
                                     for column in range(self.outcomes)]
                                    for row in range(self.outcomes)]) / np.sqrt(self.outcomes)

    def pack_metadata(self, meta):
        value, shift = 0, self.payload_width
        for part, width in zip(meta, self.field_widths):
            value |= part << shift
            shift += width
        return value

    def unpack_metadata(self, code):
        result, shift = [], self.payload_width
        for width in self.field_widths:
            result.append(code >> shift & ((1 << width) - 1))
            shift += width
        return tuple(result)

    def atom_code(self, prefix, child, rank):
        meta = self.metadata[prefix, child]
        _, low_width, _, seed = meta
        low = (child & ((1 << low_width) - 1)) ^ \
            _shear(seed, child >> low_width, low_width, self.height - low_width)
        return self.pack_metadata(meta) | low | (rank << low_width)

    def decode_custom(self, prefix, code):
        meta = self.unpack_metadata(code)
        low_width = meta[1]
        low = code & ((1 << min(low_width, self.width)) - 1)
        return self.code_inverse.get((prefix, meta, low), 0)

    def forward(self, prefix, code):
        return self.permutations[prefix].get(code, code)

    def inverse(self, prefix, code):
        return self.inverses[prefix].get(code, code)


class _CausalFixture:
    def __init__(self):
        self.blocks = (1, 5)
        later = {prefix: [random.Random(seed + node * 99991).randrange(1, 4)
                          for node in range(31)] for prefix, seed in ((0, 5), (1, 19))}
        # Seed 5 has a cardinality group containing both 1 and 17, with
        # two other points: the first hash leaves a genuine two-point collision.
        later[0] = [3, 2, 3, 2, 3, 3, 3, 3, 1, 2, 1, 3, 1, 1, 1,
                    2, 2, 1, 2, 3, 1, 3, 1, 1, 3, 1, 2, 2, 1, 2, 1]
        self.locals = (_LocalCoins(1, {0: [1]}), _LocalCoins(5, later))
        numerators = np.zeros(63, dtype=int)
        numerators[0] = 1
        for prefix, values in later.items():
            for local_node, numerator in enumerate(values, start=1):
                depth = local_node.bit_length() - 1
                global_node = (1 << (depth + 1)) + (prefix << depth) + local_node - (1 << depth)
                numerators[global_node - 1] = numerator
        self.frame = real_frame_matrix(np.arccos(np.sqrt(numerators / 4)))
        self.bands = np.zeros(64, dtype=int)
        for node in range(1, 64):
            self.bands[marker_label(node, 6)] = 1 if node == 1 else 2
        first, second = self.locals
        self.embedding = np.zeros((first.core_size * second.core_size, 64))
        self.custom_atoms = []
        for leaf in range(64):
            first_child, second_child = divmod(leaf, 32)
            for rank1, rank2 in product(range(len(first.fibers[0][first_child])),
                                        range(len(second.fibers[first_child][second_child]))):
                custom = (first.atom_code(0, first_child, rank1),
                          second.atom_code(first_child, second_child, rank2))
                target, clean = self.match_word(custom)
                if any(clean):
                    raise AssertionError("global matcher did not return its selector")
                coin1 = first.fibers[0][first_child][rank1][1]
                coin2 = second.fibers[first_child][second_child][rank2][1]
                amplitude = np.sqrt(first.atoms[coin1] /
                                    (first.outcomes * first.capacities[0][first_child]))
                amplitude *= np.sqrt(second.atoms[coin2] /
                                     (second.outcomes * second.capacities[first_child][second_child]))
                self.embedding[target[0] * second.core_size + target[1], leaf] = amplitude
                self.custom_atoms.append(custom)
        self.last = np.array([2 if second_coin else (1 if first_coin else 0)
                              for first_coin in range(first.core_size)
                              for second_coin in range(second.core_size)])

    def selectors(self, custom):
        first, second = self.locals
        prefix = first.decode_custom(0, custom[0])
        outputs = (first.forward(0, custom[0]), second.forward(prefix, custom[1]))
        return tuple(local.catalogue_index[old ^ new]
                     for local, old, new in zip(self.locals, custom, outputs))

    def inverse_selectors(self, output):
        first, second = self.locals
        old_first = first.inverse(0, output[0])
        prefix = first.decode_custom(0, old_first)
        old_second = second.inverse(prefix, output[1])
        return self.selectors((old_first, old_second))

    def match_word(self, custom):
        selector = self.selectors(custom)
        output = tuple(value ^ local.catalogue[index]
                       for value, local, index in zip(custom, self.locals, selector))
        returned = tuple(left ^ right for left, right in
                         zip(selector, self.inverse_selectors(output)))
        return output, returned

    def transform(self, states, inverse=False):
        shape = (self.locals[0].outcomes, self.locals[0].gamma_size,
                 self.locals[1].outcomes, self.locals[1].gamma_size) + states.shape[1:]
        result = states.reshape(shape)
        words = [self.locals[0].hadamard, self.locals[0].gamma_transform,
                 self.locals[1].hadamard, self.locals[1].gamma_transform]
        for axis, word in enumerate(words):
            if inverse:
                word = word.conj().T
            result = np.moveaxis(np.tensordot(word, result, axes=(1, axis)), 0, axis)
        return result.reshape(states.shape)


class CausalCompactCoinDepthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = _CausalFixture()

    def test_batched_compressor_clears_high_bits_and_actual_inverse_returns_metadata(self):
        first, second = self.fixture.locals
        self.assertTrue(any(meta[2] > 0 for meta in second.metadata.values()))
        self.assertGreater(second.rejected_collision_candidates, 0)
        self.assertTrue(any(meta[1] < second.height for meta in second.metadata.values()))
        self.assertNotEqual(second.metadata[0, 1], second.metadata[1, 1])
        for leaf in range(64):
            children = divmod(leaf, 32)
            prefixes = (0, children[0])
            metadata = [local.metadata[prefix, child]
                        for local, prefix, child in zip(self.fixture.locals, prefixes, children)]
            payloads, codes = [], []
            for local, child, meta in zip(self.fixture.locals, children, metadata):
                _, low_width, _, seed = meta
                high, low = child >> low_width, child & ((1 << low_width) - 1)
                hashed = low ^ _shear(seed, high, low_width, local.height - low_width)
                payloads.append((high << low_width) | hashed)
                codes.append(local.pack_metadata(meta) | hashed)
            decoded_first = first.decode_custom(0, codes[0])
            decoded = (decoded_first, second.decode_custom(decoded_first, codes[1]))
            self.assertEqual(decoded, children)
            for index, (local, meta) in enumerate(zip(self.fixture.locals, metadata)):
                _, low_width, _, seed = meta
                high_mask = decoded[index] >> low_width << low_width
                payloads[index] ^= high_mask
                self.assertEqual(payloads[index] >> low_width, 0)
                self.assertEqual(local.pack_metadata(meta) | payloads[index], codes[index])
                payloads[index] ^= high_mask
                payloads[index] ^= _shear(seed, payloads[index] >> low_width,
                                           low_width, local.height - low_width)
                self.assertEqual(payloads[index], children[index])
                self.assertEqual(local.pack_metadata(meta) ^
                                 local.pack_metadata(local.metadata[prefixes[index], payloads[index]]), 0)

    def test_global_short_shift_word_and_inverse_selectors(self):
        fixture = self.fixture
        outputs = set()
        for custom in fixture.custom_atoms:
            output, returned = fixture.match_word(custom)
            self.assertEqual(returned, (0, 0))
            outputs.add(output)
            first = fixture.locals[0].inverse(0, output[0])
            prefix = fixture.locals[0].decode_custom(0, first)
            second = fixture.locals[1].inverse(prefix, output[1])
            self.assertEqual((first, second), custom)
        self.assertEqual(len(outputs), len(fixture.custom_atoms))
        # This fixture chose complete local permutations, so its word is
        # additionally checked on occupied, off-image custom coordinates.
        rng = random.Random(1307)
        for _ in range(150):
            custom = tuple(rng.randrange(1 << local.width) for local in fixture.locals)
            output, returned = fixture.match_word(custom)
            self.assertEqual(returned, (0, 0))
            self.assertEqual(fixture.selectors(custom), fixture.inverse_selectors(output))

    def test_every_marker_band_and_full_tag_word_with_reference(self):
        fixture = self.fixture
        embedding, frame = fixture.embedding, fixture.frame
        np.testing.assert_allclose(embedding.T @ embedding, np.eye(64), atol=ATOL)
        transformed = fixture.transform(embedding @ frame)
        for column, band in enumerate(fixture.bands):
            np.testing.assert_allclose(transformed[fixture.last != band, column], 0, atol=ATOL)
            self.assertAlmostEqual(np.linalg.norm(transformed[:, column]), 1, delta=ATOL)
        rng = np.random.default_rng(3301)
        coefficients = rng.normal(size=(64, 4, 3)) + 1j * rng.normal(size=(64, 4, 3))
        coefficients /= np.linalg.norm(coefficients)
        logical = np.einsum('xc,ctr->xtr', frame, coefficients)
        expected_coefficients = np.empty_like(coefficients)
        coarse = np.maximum(fixture.bands - 1, 0)
        for column, label in enumerate(coarse):
            for tag in range(4):
                expected_coefficients[column, tag ^ int(label)] = coefficients[column, tag]
        expected = np.einsum('xc,ctr->xtr', frame, expected_coefficients)
        coin = fixture.transform(np.einsum('zx,xtr->ztr', embedding, logical))
        toggled = np.empty_like(coin)
        for row, band in enumerate(fixture.last):
            for tag in range(4):
                toggled[row, tag ^ max(int(band) - 1, 0)] = coin[row, tag]
        returned_coin = fixture.transform(toggled, inverse=True)
        output = np.einsum('xz,ztr->xtr', embedding.T, returned_coin)
        np.testing.assert_allclose(output, expected, atol=ATOL)
        np.testing.assert_allclose(returned_coin,
                                   np.einsum('zx,xtr->ztr', embedding, expected), atol=ATOL)

    def test_reference_source_work_must_be_part_of_nonzero_block(self):
        plus, minus = np.ones(2) / np.sqrt(2), np.array((1., -1.)) / np.sqrt(2)
        zero, one = np.eye(2)
        completion = np.column_stack((np.kron(plus, zero), np.kron(minus, zero),
                                      np.kron(plus, one), np.kron(minus, one)))
        np.testing.assert_allclose(completion.conj().T @ completion, np.eye(4), atol=ATOL)
        transformed = completion.conj().T @ np.kron(minus, zero)
        np.testing.assert_allclose(transformed, np.array((0., 1., 0., 0.)), atol=ATOL)
        coin_zero_probability = np.linalg.norm(transformed.reshape(2, 2)[0])**2
        self.assertAlmostEqual(coin_zero_probability, 1, delta=ATOL)
        self.assertAlmostEqual(abs(transformed[0])**2, 0, delta=ATOL)


if __name__ == '__main__':
    unittest.main()
