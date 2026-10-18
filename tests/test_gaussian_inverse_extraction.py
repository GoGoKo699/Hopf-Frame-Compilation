"""Exact finite interfaces for the Gaussian inverse-extraction proof.

The full graph work is arbitrary in the exhaustive two-bit fixtures. Integer
phases retain the scalar sign; reversed edge words are actual inverses. These
tests do not establish a general query-depth lower bound or scalable emitter.
The field fixtures check the dirty-offset interpolation interface, not a
native implementation of a finite-field multiplier.
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


def _gf4_mul(a, b):
    """F_2[t]/(t^2+t+1), with bit zero the coordinate of one."""
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return out


def _gf4_pow(a, exponent):
    out = 1  # In particular, the interpolation endpoint uses 0**0 = 1.
    for _ in range(exponent):
        out = _gf4_mul(out, a)
    return out


def _upper_pairs(order):
    return [(i, j) for i in range(order) for j in range(i + 1, order)]


def _field_upper_inverse(order, coefficients):
    """Independent unit-upper back substitution over F4."""
    upper = [[0] * order for _ in range(order)]
    for (i, j), value in zip(_upper_pairs(order), coefficients):
        upper[i][j] = value
    inverse = [[int(i == j) for j in range(order)] for i in range(order)]
    for i in range(order - 1, -1, -1):
        for j in range(i + 1, order):
            for k in range(i + 1, j + 1):
                inverse[i][j] ^= _gf4_mul(upper[i][k], inverse[k][j])
    return tuple(inverse[i][j] for i, j in _upper_pairs(order))


def _binary_path_inverse(order, coefficients):
    """Binary expected output by direct enumeration of increasing paths."""
    upper = dict(zip(_upper_pairs(order), coefficients))
    out = 0
    for bit, (i, j) in enumerate(_upper_pairs(order)):
        value = 0
        for mask in range(1 << (j - i - 1)):
            vertices = [i] + [k for k in range(i + 1, j)
                              if mask & (1 << (k - i - 1))] + [j]
            product = 1
            for source, target in zip(vertices, vertices[1:]):
                product &= upper[source, target]
            value ^= product
        out |= value << bit
    return out


def _graded_readout(order, coefficients, mask, output):
    """Load/query/unload with one common arbitrary F4 coefficient mask."""
    pairs = _upper_pairs(order)
    for point in range(4):
        load = tuple(_gf4_mul(_gf4_pow(point, j - i), value)
                     for (i, j), value in zip(pairs, coefficients))
        mask = tuple(value ^ offset for value, offset in zip(mask, load))
        inverse = _field_upper_inverse(order, mask)
        for bit, ((i, j), value) in enumerate(zip(pairs, inverse)):
            weight = _gf4_pow(point, 3 - (j - i))
            output ^= (_gf4_mul(weight, value) & 1) << bit
        mask = tuple(value ^ offset for value, offset in zip(mask, load))
    return output, mask


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

    def test_dirty_translation_echo_retains_mask_terms(self):
        for u, v, w, x, y in itertools.product((0, 1), repeat=5):
            before = _field_upper_inverse(3, (u, w, v))[1]
            after = _field_upper_inverse(3, (u ^ x, w, v ^ y))[1]
            self.assertEqual(before ^ after, (x & v) ^ (y & u) ^ (x & y))
        # At x=y=1 the desired scalar is -1, but this occupied mask gives +1.
        before = _field_upper_inverse(3, (1, 0, 0))[1]
        after = _field_upper_inverse(3, (0, 0, 1))[1]
        self.assertEqual(before ^ after, 0)
        self.assertEqual(abs((1 - 2 * (before ^ after)) - (-1)), 2)

    def test_graded_field_readout_all_small_dirty_columns(self):
        for coefficients in itertools.product((0, 1), repeat=3):
            expected = _binary_path_inverse(3, coefficients)
            for mask in itertools.product(range(4), repeat=3):
                for output in range(8):
                    actual, returned = _graded_readout(3, coefficients, mask,
                                                      output)
                    self.assertEqual(actual, output ^ expected)
                    self.assertEqual(returned, mask)

        # Order four reaches k=|F4|-1, including the t=0, exponent-zero case.
        masks = [(value,) * 6 for value in range(4)]
        masks += [tuple(int(i == j) for i in range(6)) for j in range(6)]
        masks += [(0, 1, 2, 3, 1, 2), (3, 0, 2, 1, 3, 1)]
        for coefficients in itertools.product((0, 1), repeat=6):
            expected = _binary_path_inverse(4, coefficients)
            for mask in masks:
                for output in (0, 21, 63):
                    actual, returned = _graded_readout(4, coefficients, mask,
                                                      output)
                    self.assertEqual(actual, output ^ expected)
                    self.assertEqual(returned, mask)

    def test_no_short_fixed_field_linear_coefficient_extractor(self):
        # pi(weight * value) exhausts all four binary-linear F4 readouts.
        readouts = {tuple(_gf4_mul(weight, value) & 1 for value in range(4))
                    for weight in range(4)}
        self.assertEqual(len(readouts), 4)
        for degree in (2, 3):
            def valid(points, weights):
                for power in range(degree + 1):
                    # Lower coefficients are arbitrary F4; the leading
                    # coefficient of the binary-U family is only binary.
                    for coefficient in ((1, 2) if power < degree else (1,)):
                        actual = 0
                        for point, weight in zip(points, weights):
                            value = _gf4_mul(coefficient,
                                            _gf4_pow(point, power))
                            actual ^= _gf4_mul(weight, value) & 1
                        expected = int(power == degree) * (coefficient & 1)
                        if actual != expected:
                            return False
                return True

            for points in itertools.combinations(range(4), degree):
                self.assertFalse(any(valid(points, weights) for weights in
                                     itertools.product(range(4), repeat=degree)))
            points = tuple(range(4))
            weights = tuple(_gf4_pow(point, 3 - degree) for point in points)
            self.assertTrue(valid(points, weights))


if __name__ == "__main__":
    unittest.main()
