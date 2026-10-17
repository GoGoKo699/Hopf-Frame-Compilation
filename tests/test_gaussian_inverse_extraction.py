"""Exact finite interfaces for the Gaussian inverse-extraction proof.

The full graph work is arbitrary in the exhaustive two-bit fixtures. Integer
phases retain the scalar sign; reversed edge words are actual inverses. These
tests do not establish a general query-depth lower bound or scalable emitter.
"""
from __future__ import annotations

import itertools
import unittest


def _words(length):
    return list(itertools.product((0, 1), repeat=length))


def _graph(address_bits, table):
    """Return topologically ordered edges (source, target, literal-or-None)."""
    p = address_bits // 2
    q = address_bits - p
    left, right = {}, {}
    cursor = 0
    for depth in range(p + 1):
        for prefix in _words(depth):
            left[prefix] = cursor
            cursor += 1
    for depth in range(q, -1, -1):
        for prefix in _words(depth):
            right[prefix] = cursor
            cursor += 1
    edges = []
    for depth in range(1, p + 1):
        for prefix in _words(depth):
            edges.append((left[prefix[:-1]], left[prefix],
                          (depth - 1, prefix[-1])))
    for i, u in enumerate(_words(p)):
        for j, v in enumerate(_words(q)):
            if table[(i << q) | j]:
                edges.append((left[u], right[v], None))
    for depth in range(q, 0, -1):
        for prefix in _words(depth):
            edges.append((right[prefix], right[prefix[:-1]],
                          (p + depth - 1, prefix[-1])))
    return cursor, edges, left[()], right[()]


def _address(word, bits):
    return tuple((word >> (bits - 1 - i)) & 1 for i in range(bits))


def _weight(literal, x):
    return literal is None or x[literal[0]] == literal[1]


def _edge_word(z, edges, x):
    for u, v, literal in edges:
        if _weight(literal, x) and ((z >> u) & 1):
            z ^= 1 << v
    return z


def _n_apply(z, edges, x):
    out = 0
    for u, v, literal in edges:
        if _weight(literal, x) and ((z >> u) & 1):
            out ^= 1 << v
    return out


def _inverse_apply(z, width, edges, x):
    """Independent nilpotent Neumann series, including every path length."""
    out, term = z, z
    for _ in range(1, width):
        term = _n_apply(term, edges, x)
        out ^= term
    if _n_apply(term, edges, x):
        raise AssertionError("graph was not nilpotent at its claimed order")
    return out


def _displacement(z, edges, x, source):
    # Chronological C_A, X_b, C_A^dagger. Both use the literal edge word.
    z = _edge_word(z, reversed(edges), x)
    z ^= 1 << source
    return _edge_word(z, edges, x)


def _anf_degree(values):
    coeff = list(values)
    bits = (len(coeff) - 1).bit_length()
    for i in range(bits):
        for mask in range(len(coeff)):
            if mask & (1 << i):
                coeff[mask] ^= coeff[mask ^ (1 << i)]
    return max((mask.bit_count() for mask, value in enumerate(coeff) if value),
               default=-1)


class GaussianInverseExtractionTests(unittest.TestCase):
    def test_all_two_bit_tables_full_dirty_work_and_literal_phase(self):
        for table_word in range(16):
            table = [(table_word >> x) & 1 for x in range(4)]
            width, edges, source, sink = _graph(2, table)
            self.assertEqual(width, 6)
            for address in range(4):
                x = _address(address, 2)
                shift = _inverse_apply(1 << source, width, edges, x)
                self.assertEqual((shift >> sink) & 1, table[address])
                for z in range(1 << width):
                    # The topological edge word is A^-1; its actual reverse is A.
                    forward = _edge_word(z, edges, x)
                    self.assertEqual(forward, _inverse_apply(z, width, edges, x))
                    reverse = _edge_word(z, reversed(edges), x)
                    self.assertEqual(reverse, z ^ _n_apply(z, edges, x))
                    self.assertEqual(_edge_word(forward, reversed(edges), x), z)

                    # Every work basis word returns, with exactly the table sign.
                    first = _displacement(z, edges, x, source)
                    self.assertEqual(first, z ^ shift)
                    second = _displacement(first, edges, x, source)
                    phase = ((first >> sink) ^ (second >> sink)) & 1
                    self.assertEqual(second, z)
                    self.assertEqual(phase, table[address])
                    for control in (0, 1):
                        self.assertEqual(control * phase,
                                         control * table[address])

    def test_cheap_shear_extracts_direct_entry_not_inverse_entry(self):
        width, edges, source, sink = _graph(2, [0, 0, 0, 1])
        x = (1, 1)
        b = 1 << source
        ab = b ^ _n_apply(b, edges, x)
        inverse_b = _inverse_apply(b, width, edges, x)
        self.assertEqual((ab >> sink) & 1, 0)
        self.assertEqual((inverse_b >> sink) & 1, 1)

        def shear(y, z):
            return y ^ z ^ _n_apply(z, edges, x), z

        def conjugate(y, z):
            y, z = shear(y, z)
            z ^= b
            return shear(y, z)

        for y in range(1 << width):
            for z in range(1 << width):
                y1, z1 = conjugate(y, z)
                self.assertEqual((y1, z1), (y ^ ab, z ^ b))
                y2, z2 = conjugate(y1, z1)
                self.assertEqual((y2, z2), (y, z))
                self.assertEqual(((y1 >> sink) ^ (y2 >> sink)) & 1, 0)

    def test_inverse_entry_degree_and_serial_resource_envelope(self):
        for address_bits in range(1, 8):
            size = 1 << address_bits
            table = [0] * size
            table[-1] = 1
            width, edges, source, sink = _graph(address_bits, table)
            p = address_bits // 2
            q = address_bits - p
            self.assertEqual(width, 2 * ((1 << p) + (1 << q)) - 2)
            self.assertEqual(sum(lit is not None for _, _, lit in edges),
                             2 * ((1 << p) + (1 << q)) - 4)
            self.assertEqual(sum(lit is None for _, _, lit in edges), 1)
            inverse_entry, direct_entry = [], []
            for address in range(size):
                x = _address(address, address_bits)
                v = _inverse_apply(1 << source, width, edges, x)
                inverse_entry.append((v >> sink) & 1)
                direct = (1 << source) ^ _n_apply(1 << source, edges, x)
                direct_entry.append((direct >> sink) & 1)
            self.assertEqual(inverse_entry, table)
            self.assertEqual(_anf_degree(inverse_entry), address_bits)
            self.assertLessEqual(_anf_degree(direct_entry), 1)


if __name__ == "__main__":
    unittest.main()
