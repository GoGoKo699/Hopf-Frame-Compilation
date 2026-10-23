"""Exact finite witnesses for scoped quarter-phase matching limits.

Monomial phases are represented by powers of i, without floating point.
These are algebra fixtures, not full arbitrary-bank native emitters or
an asymptotic verification of the inherited expander family.
"""
from collections import deque
from itertools import permutations, product
import unittest


def _local_rotation(control, target):
    return tuple((x ^ (((x >> control) & 1) << target),
                  (x >> control) & 1) for x in range(4))


def _compose(left, right):
    return tuple((left[x][0], (phase + left[x][1]) % 4)
                 for x, phase in right)


def _local_group():
    identity = tuple((x, 0) for x in range(4))
    generators = (_local_rotation(0, 1), _local_rotation(1, 0))
    seen = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for generator in generators:
            value = _compose(current, generator)
            if value not in seen:
                seen.add(value)
                queue.append(value)
    return tuple(sorted(seen))


def _matching_action(x, pairs, anchor=0):
    shifted = x ^ anchor
    output, phase = shifted, 0
    for control, target in pairs:
        bit = (shifted >> control) & 1
        output ^= bit << target
        phase += bit
    return output ^ anchor, phase % 4


def _apply_monomial(vector, action):
    # Coefficients are Gaussian integers stored as (real, imaginary).
    phases = ((1, 0), (0, 1), (-1, 0), (0, -1))
    out = {}
    for x, (real, imag) in vector.items():
        y, phase = action(x)
        a, b = phases[phase]
        old_real, old_imag = out.get(y, (0, 0))
        value = (old_real + a * real - b * imag,
                 old_imag + a * imag + b * real)
        if value == (0, 0):
            out.pop(y, None)
        else:
            out[y] = value
    return out


def _cellwise_action(x, words):
    output = phase = 0
    for j, word in enumerate(words):
        y, exponent = word[(x >> (2 * j)) & 3]
        output |= y << (2 * j)
        phase += exponent
    return output, phase % 4


class QuarterPhaseMatchingLimitTests(unittest.TestCase):
    def test_exact_local_group_is_real_orientation_preserving_after_gauge(self):
        group = _local_group()
        self.assertEqual(len(group), 24)
        gauge = (0, 1, 1, 2)  # S tensor S.
        real_forms = set()
        for word in group:
            self.assertEqual(word[0], (0, 0))
            gauged = tuple((y, (phase + gauge[x] - gauge[y]) % 4)
                          for x, (y, phase) in enumerate(word))
            self.assertTrue(all(phase in (0, 2) for _, phase in gauged))
            perm = [gauged[x][0] for x in (1, 2, 3)]
            inversions = sum(perm[i] > perm[j]
                             for i in range(3) for j in range(i + 1, 3))
            signs = sum(gauged[x][1] // 2 for x in (1, 2, 3))
            self.assertEqual((inversions + signs) % 2, 0)
            real_forms.add(gauged)
        self.assertEqual(len(real_forms), 24)

    def test_all_local_group_elements_fix_the_two_cell_singlet(self):
        chi = {5: (-1, 0), 10: (-1, 0), 15: (1, 0)}
        for word in _local_group():
            self.assertEqual(_apply_monomial(
                chi, lambda x: _cellwise_action(x, (word, word))), chi)

    def test_duplicate_signatures_survive_heterogeneous_address_layers(self):
        # 24^2+1 cells force a repeated two-layer signature.  The sparse
        # tensor witness uses only three big-integer computational labels.
        group = _local_group()
        signatures = list(product(group, repeat=2))
        signatures.append(signatures[17])
        first, second = 17, len(signatures) - 1
        chi = {((a << (2 * first)) | (a << (2 * second))): (sign, 0)
               for a, sign in ((1, -1), (2, -1), (3, 1))}
        for layer in range(2):
            words = tuple(signature[layer] for signature in signatures)
            self.assertEqual(_apply_monomial(
                chi, lambda x: _cellwise_action(x, words)), chi)

    def test_all_fixed_bipartition_pairings_fix_even_size_witness(self):
        for p in (4, 8):
            left, right = tuple(range(0, p, 2)), tuple(range(1, p, 2))
            a = sum(1 << bit for bit in left)
            b = sum(1 << bit for bit in right)
            sign = 1 if (len(left) // 2) % 2 == 0 else -1
            base = {a: (1, 0), b: (1, 0), a | b: (sign, 0)}
            for order in permutations(right):
                forward = tuple(zip(left, order))
                reverse = tuple((target, control) for control, target in forward)
                for pairs in (forward, reverse):
                    for anchor in (0, 1, (1 << p) - 1):
                        witness = {x ^ anchor: coefficient
                                   for x, coefficient in base.items()}
                        self.assertNotIn(anchor, witness)
                        self.assertEqual(_matching_action(anchor, pairs, anchor),
                                         (anchor, 0))
                        self.assertEqual(_apply_monomial(
                            witness, lambda x: _matching_action(x, pairs, anchor)),
                            witness)
                        self.assertEqual(_apply_monomial(
                            witness, lambda x: (
                                _matching_action(x, pairs, anchor)[0],
                                -_matching_action(x, pairs, anchor)[1] % 4)),
                            witness)

    def test_changed_control_partition_breaks_the_same_witness(self):
        witness = {5: (1, 0), 10: (1, 0), 15: (-1, 0)}
        actual = _apply_monomial(
            witness, lambda x: _matching_action(x, ((0, 2), (1, 3))))
        self.assertEqual(actual, {1: (0, 1), 2: (0, 1), 3: (1, 0)})
        self.assertTrue(set(actual).isdisjoint(witness))

    def test_matching_square_is_the_nonempty_control_mask(self):
        for p, pairs in ((2, ((0, 1),)), (4, ((0, 2), (1, 3))),
                         (8, ((0, 3), (1, 5), (2, 6), (4, 7)))):
            control_mask = sum(1 << control for control, _ in pairs)
            for x in range(1 << p):
                y, first = _matching_action(x, pairs)
                z, second = _matching_action(y, pairs)
                self.assertEqual(z, x)
                self.assertEqual((first + second) % 4,
                                 2 * ((x & control_mask).bit_count() % 2))


if __name__ == '__main__':
    unittest.main()
