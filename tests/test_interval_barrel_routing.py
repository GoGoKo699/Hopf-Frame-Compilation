"""Bounded logical permutation checks for the contiguous-interval router.

Queries are represented by their exact classical XOR interfaces.  The
emitted words explicitly retain dirty controls and both Fredkin banks of
each echo.  Unique wire tokens test arbitrary routed/padded values without
assuming zeros.  These fixtures are not a native Clifford+T emitter,
phase-synthesis check, or proof of asymptotic resource bounds.
"""
from __future__ import annotations

import itertools
import unittest


def _shift_matchings(width, step):
    """Two fixed involutions whose chronological action moves i to i-step."""
    unseen = set(range(width))
    matchings = [[], []]
    while unseen:
        first = min(unseen)
        cycle = []
        site = first
        while site not in cycle:
            cycle.append(site)
            unseen.remove(site)
            site = (site - step) % width
        for index, physical_site in enumerate(cycle):
            for matching, reflected in enumerate((-index, 1 - index)):
                partner = reflected % len(cycle)
                if index < partner:
                    matchings[matching].append((physical_site, cycle[partner]))
    return matchings


def _matching_echo(pairs, predicates):
    """Each query target is a private dirty bit, reused after its return."""
    queries = [("xor", control, predicate)
               for control, predicate in enumerate(predicates)]
    swaps = [("swap", control, left, right)
             for control, (left, right) in enumerate(pairs)]
    return queries + swaps + list(reversed(queries)) + swaps


def _rotation_word(width, start_name):
    word = []
    for bit in range(width.bit_length() - 1):
        for matching in _shift_matchings(width, 1 << bit):
            word += _matching_echo(
                matching, [("bit", start_name, bit)] * len(matching))
    return word


def _interval_word(width, bank_offset, bank_capacity, start_name, length_name):
    rotate = _rotation_word(width, start_name)
    pairs = [(site, bank_offset + site) for site in range(bank_capacity)]
    predicates = [("prefix", length_name, site) for site in range(bank_capacity)]
    return rotate + _matching_echo(pairs, predicates) + list(reversed(rotate))


def _run(word, tokens, dirty, active, parameters):
    """All primitive gates are involutions; reversing word is the actual inverse."""
    tokens, dirty = list(tokens), list(dirty)
    for gate in word:
        if gate[0] == "xor":
            _, control, (kind, name, index) = gate
            value = parameters[name]
            predicate = ((value >> index) & 1) if kind == "bit" else int(index < value)
            dirty[control] ^= active & predicate
        else:
            _, control, left, right = gate
            if dirty[control]:
                tokens[left], tokens[right] = tokens[right], tokens[left]
    return tokens, dirty


def _direct_exchange(tokens, start, length, bank_offset):
    """Independent specification: exchange corresponding slice positions."""
    expected = list(tokens)
    for index in range(length):
        expected[start + index] = tokens[bank_offset + index]
        expected[bank_offset + index] = tokens[start + index]
    return expected


class IntervalBarrelRoutingTests(unittest.TestCase):
    def test_variable_intervals_restore_arbitrary_padding_and_dirty_controls(self):
        for physical_length in (3, 5):
            width = 1 << (physical_length - 1).bit_length()
            capacity = 3
            word = _interval_word(width, width, capacity, "start", "length")
            tokens = list(range(width + capacity))
            control_count = max(width // 2, capacity)
            for start in range(physical_length + 1):
                for length in range(min(capacity, physical_length - start) + 1):
                    parameters = {"start": start, "length": length}
                    for dirty in itertools.product((0, 1), repeat=control_count):
                        with self.subTest(L=physical_length, start=start,
                                          length=length, dirty=dirty):
                            actual, returned = _run(word, tokens, dirty, 1, parameters)
                            self.assertEqual(actual, _direct_exchange(
                                tokens, start, length, width))
                            self.assertEqual(actual[physical_length:width],
                                             tokens[physical_length:width])
                            self.assertEqual(actual[width + length:],
                                             tokens[width + length:])
                            self.assertEqual(returned, list(dirty))
                            restored, returned = _run(
                                reversed(word), actual, returned, 1, parameters)
                            self.assertEqual(restored, tokens)
                            self.assertEqual(returned, list(dirty))

    def test_inactive_sector_is_identity_for_arbitrary_parameter_words(self):
        width, capacity = 8, 3
        word = _interval_word(width, width, capacity, "start", "length")
        tokens = list(range(width + capacity))
        # Include starts outside the unpadded array and lengths above capacity.
        for start, length in itertools.product(range(width), repeat=2):
            parameters = {"start": start, "length": length}
            for dirty in itertools.product((0, 1), repeat=width // 2):
                actual, returned = _run(word, tokens, dirty, 0, parameters)
                self.assertEqual(actual, tokens)
                self.assertEqual(returned, list(dirty))
                restored, returned = _run(
                    reversed(word), actual, returned, 0, parameters)
                self.assertEqual(restored, tokens)
                self.assertEqual(returned, list(dirty))

    def test_two_disjoint_intervals_preserve_parked_bank_contents(self):
        physical_length, width = 6, 8
        local_capacity, coin_capacity = 3, 4
        local_offset, coin_offset = width, width + local_capacity
        local_word = _interval_word(
            width, local_offset, local_capacity, "a", "q")
        coin_word = _interval_word(
            width, coin_offset, coin_capacity, "b", "r")
        layout = local_word + coin_word
        tokens = list(range(width + local_capacity + coin_capacity))
        for a in range(physical_length):
            for q in range(1, min(local_capacity, physical_length - a) + 1):
                b = a + q
                for r in range(min(coin_capacity, physical_length - b) + 1):
                    parameters = {"a": a, "b": b, "q": q, "r": r}
                    for dirty in itertools.product((0, 1), repeat=coin_capacity):
                        routed, returned = _run(layout, tokens, dirty, 1, parameters)
                        expected = _direct_exchange(tokens, a, q, local_offset)
                        expected = _direct_exchange(expected, b, r, coin_offset)
                        self.assertEqual(routed, expected)
                        self.assertEqual(returned, list(dirty))
                        self.assertEqual(routed[:a], tokens[:a])
                        self.assertEqual(routed[physical_length:width],
                                         tokens[physical_length:width])
                        restored, returned = _run(
                            reversed(layout), routed, returned, 1, parameters)
                        self.assertEqual(restored, tokens)
                        self.assertEqual(returned, list(dirty))

                        # A nontrivial body acts only on active external slots.
                        # Reverse each active bank, then unroute the ACTUAL word.
                        body = list(routed)
                        body[local_offset:local_offset + q] = reversed(
                            body[local_offset:local_offset + q])
                        body[coin_offset:coin_offset + r] = reversed(
                            body[coin_offset:coin_offset + r])
                        actual, returned = _run(
                            reversed(layout), body, dirty, 1, parameters)
                        expected = list(tokens)
                        expected[a:b] = reversed(tokens[a:b])
                        expected[b:b + r] = reversed(tokens[b:b + r])
                        self.assertEqual(actual, expected)
                        self.assertEqual(actual[local_offset:], tokens[local_offset:])
                        self.assertEqual(returned, list(dirty))


if __name__ == "__main__":
    unittest.main()
