"""Exact fixtures for high-rank occupied Clifford extension failures.

The analytic family has one initialized qubit and arbitrary data. These
checks retain the complete logical family, both clean sectors, literal
native phases, the physical axis basis, and the actual inverse ledger.
They establish no hard target: the occupied transition is identity.
"""
from collections import Counter
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, Q2, add, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clean_code_clifford_rigidity import (
    _high_rank_words, _pauli_label, _rank,
)
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _conjugate, _tau


def _high_rank_extension_words(width):
    """Width k>=2: wire zero clean; the other k-1 wires arbitrary."""
    if width < 2:
        raise ValueError('the family needs at least one data wire')
    # K=fanout(c -> data) H_c maps Z_c to X on every physical wire.
    k = [('H', 0)] + [('CX', 0, j) for j in range(1, width)]
    e = [('T', j) for j in range(width)]
    # D Z_1 D†=P and D Z_c D†=Z_c P, with P=product(data X).
    d = ([('CX', 1, 0), ('H', 1)]
         + [('CX', 1, j) for j in range(2, width)])
    e0 = _inverse(d) + [('TDG', 0), ('T', 1)] + d
    ua = k + e
    ub = e0 + ua
    # The two adjacent E bands are scheduled separately before cancellation.
    displayed = ua + _inverse(ub)
    return k, e, d, e0, ua, ub, displayed, _inverse(e0)


def _clifford_labels(matrix, width):
    labels = []
    for bit in range(2 * width):
        label = _pauli_label(_conjugate(matrix, _tau(1 << bit, width)), width)
        if label is None:
            raise ValueError('not a Clifford')
        labels.append(label)
    return labels


class OccupiedLogicalTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.examples = []
        for width in (2, 3, 4):
            words = _high_rank_extension_words(width)
            k, e, d, e0, ua, ub = [native_columns(width, word)
                                    for word in words[:6]]
            clean = _tau(1 << width, width)
            s = _conjugate(ua, clean)
            cls.examples.append((width, words, k, e, d, e0, ua, ub, clean, s))

    def test_native_axes_and_single_clean_image_have_arbitrary_active_rank(self):
        for width, _, k, e, d, e0, ua, ub, clean, s in self.examples:
            with self.subTest(width=width):
                p_label = (1 << width) - 2
                p = _tau(p_label, width)
                self.assertEqual(_conjugate(d, _tau(1 << (width + 1), width)), p)
                self.assertEqual(_conjugate(d, clean), multiply(clean, p))
                self.assertEqual(_rank([p_label, p_label | (1 << width)]), 2)
                self.assertEqual(_conjugate(k, clean), _tau((1 << width) - 1, width))
                self.assertEqual(_conjugate(ub, clean), s)
                self.assertEqual(adjoint(s), s)
                self.assertEqual(multiply(s, s), eye(1 << width))
                columns = _clifford_labels(s, width)
                self.assertEqual(_rank([label ^ (1 << bit)
                                        for bit, label in enumerate(columns)]), width)
                # The code has dimension 2^(k-1), not k clean conditions.
                self.assertEqual(sum((s[j][j] for j in range(1 << width)), CZ), CZ)

    def test_every_logical_generator_is_clifford_only_on_the_occupied_code(self):
        for width, _, k, e, _, e0, ua, ub, clean, s in self.examples:
            twice_code = add(eye(1 << width), s)
            y_clean = _tau(1 | (1 << width), width)
            expected_probe = _conjugate(e, y_clean)
            self.assertEqual(expected_probe,
                             scale(add(y_clean, scale(_tau(1, width), -1)),
                                   Q2(0, 1) / 2))
            self.assertIsNone(_pauli_label(expected_probe, width))
            for data in range(1, width):
                with self.subTest(width=width, data=data):
                    z = _tau(1 << (width + data), width)
                    c, a = _conjugate(ua, z), _conjugate(ub, z)
                    self.assertEqual(c, multiply(clean, z))
                    self.assertEqual(multiply(a, s), multiply(s, a))
                    self.assertEqual(multiply(a, twice_code), multiply(c, twice_code))
                    j = _conjugate(multiply(k, e0), z)
                    _clifford_labels(j, width)
                    self.assertEqual(_conjugate(e, j), a)
                    self.assertEqual(_conjugate(j, clean), y_clean)
                    self.assertEqual(_conjugate(a, clean), expected_probe)
                    with self.assertRaises(ValueError):
                        _clifford_labels(a, width)

    def test_complete_commuting_logical_family_and_all_products(self):
        for width, _, _, _, _, _, ua, ub, _, s in self.examples:
            twice_code = add(eye(1 << width), s)
            generators = [_conjugate(ub, _tau(1 << (width + data), width))
                          for data in range(1, width)]
            for left in generators:
                for right in generators:
                    self.assertEqual(multiply(left, right), multiply(right, left))
            for subset in range(1 << (width - 1)):
                with self.subTest(width=width, subset=subset):
                    z = _tau((subset << 1) << width, width)
                    a, c = _conjugate(ub, z), _conjugate(ua, z)
                    product = eye(1 << width)
                    for bit, generator in enumerate(generators):
                        if subset >> bit & 1:
                            product = multiply(product, generator)
                    self.assertEqual(product, a)
                    self.assertEqual(multiply(a, twice_code), multiply(c, twice_code))
                    # Even products commute E0 and agree on both sectors;
                    # every odd product retains the non-Clifford extension.
                    if subset.bit_count() % 2 == 0:
                        self.assertEqual(a, c)
                    else:
                        self.assertNotEqual(a, c)

    def test_full_transition_literal_unused_block_and_actual_inverse(self):
        inverse_sqrt_two = Q2(0, 1) / 2
        for width, words, _, _, _, _, ua, ub, _, _ in self.examples:
            with self.subTest(width=width):
                displayed, baseline = words[6:]
                actual = native_columns(width, displayed)
                self.assertEqual(actual, multiply(adjoint(ub), ua))
                self.assertEqual(actual, native_columns(width, baseline))
                self.assertEqual(native_columns(width, _inverse(displayed)), adjoint(actual))
                self.assertEqual(native_columns(width, _inverse(baseline)), adjoint(actual))
                p_mask = (1 << width) - 2
                for column in range(1 << width):
                    for row in range(1 << width):
                        if not column & 1:
                            expected = CO if row == column else CZ
                        elif row == column:
                            expected = CO * inverse_sqrt_two
                        elif row == (column ^ p_mask):
                            expected = CI * inverse_sqrt_two
                        else:
                            expected = CZ
                        self.assertEqual(actual[row][column], expected)
                self.assertNotEqual(actual, eye(1 << width))
                self.assertNotEqual(actual, adjoint(actual))

    def test_no_shared_pauli_clean_flip_even_with_all_bad_logical_generators(self):
        for width, _, _, _, _, _, ua, ub, _, s in self.examples:
            # Exhaust the physical Pauli group through k=3, and every
            # possible Z-type common preimage at k=4. Analytically a
            # non-Z Pauli expands into distinct terms through E.
            labels = (range(1 << (2 * width)) if width <= 3
                      else [z << width for z in range(1 << width)])
            common = []
            for label in labels:
                pauli = _tau(label, width)
                first = _pauli_label(_conjugate(adjoint(ua), pauli), width)
                z = label >> width
                is_z_type = label & ((1 << width) - 1) == 0
                self.assertEqual(first is not None, is_z_type)
                if first is None:
                    continue
                second = _pauli_label(_conjugate(adjoint(ub), pauli), width)
                expected = not z & 1 and (z >> 1).bit_count() % 2 == 0
                self.assertEqual(second is not None, expected)
                if expected:
                    common.append(label)
                    self.assertEqual(_conjugate(pauli, s), s)
            self.assertEqual(len(common), 1 << (width - 2))

    def test_scalable_ledgers_include_coordinates_and_both_inverse_words(self):
        for width in (2, 3, 4, 9):
            words = _high_rank_extension_words(width)
            k, e, d, e0, ua, ub, displayed, baseline = words
            self.assertEqual(len(k), width)
            self.assertEqual(len(d), width)
            self.assertEqual([len(band) for band in _native_t_bands(ua)], [width])
            self.assertEqual([len(band) for band in _native_t_bands(ub)], [2, width])
            for word, t_count, cliffords in (
                (displayed, 2 * width + 2, 4 * width),
                (_inverse(displayed), 2 * width + 2, 4 * width),
                (baseline, 2, 2 * width),
                (_inverse(baseline), 2, 2 * width),
            ):
                counts = Counter(gate[0] for gate in word)
                self.assertEqual(counts['T'] + counts['TDG'], t_count)
                self.assertEqual(counts['H'] + counts['CX'], cliffords)
                self.assertEqual(set(counts) - {'T', 'TDG', 'H', 'CX'}, set())
                self.assertEqual({wire for gate in word for wire in gate[1:]},
                                 set(range(width)))
            for word in (baseline, _inverse(baseline)):
                bands = _native_t_bands(word)
                self.assertEqual([len(band) for band in bands], [2])
                self.assertEqual({gate[1] for gate in bands[0]}, {0, 1})

    def test_actual_paired_differences_commute_all_bands_and_transport_globally(self):
        for pairs in (1, 2):
            width = 2 * pairs
            ua_word, ub_word, displayed, _ = _high_rank_words(pairs)
            # In this existing native family, K_A=K_C=K and E_A=E_C=E.
            k = native_columns(width, ua_word[:-width])
            e = native_columns(width, [('T', bit) for bit in range(width)])
            ua, ub = native_columns(width, ua_word), native_columns(width, ub_word)
            v = native_columns(width, displayed)
            clean = _tau(1 << width, width)
            j_clean = _conjugate(ua, clean)
            self.assertEqual(_conjugate(ub, clean), j_clean)
            endpoint_axes = [_tau(1 << (width + bit), width) for bit in range(width)]
            middle_axes = []
            differences = []
            for pair in range(pairs):
                c, d = 2 * pair, 2 * pair + 1
                coordinates = native_columns(width, [('CX', c, d), ('H', c), ('CX', c, d)])
                b = _conjugate(coordinates, endpoint_axes[c])
                other = _conjugate(coordinates, endpoint_axes[d])
                middle_axes.extend([b, other])
                # Actual axes are XX and -YY; their literal product is ZZ.
                r = multiply(b, other)
                expected_r = multiply(endpoint_axes[c], endpoint_axes[d])
                self.assertEqual(r, expected_r)
                self.assertEqual(_conjugate(j_clean, b), scale(other, -1))
                self.assertEqual(_conjugate(j_clean, other), scale(b, -1))
                differences.append(r)
            self.assertEqual(_rank([_pauli_label(r, width) for r in differences]), pairs)
            for pair, r in enumerate(differences):
                for axis in endpoint_axes + middle_axes:
                    self.assertEqual(multiply(r, axis), multiply(axis, r))
                self.assertEqual(_conjugate(e, r), r)
                self.assertEqual(_conjugate(j_clean, r), r)
                p_a = _conjugate(adjoint(k), r)
                p_c = _conjugate(adjoint(k), r)
                self.assertEqual(p_a, endpoint_axes[2 * pair + 1])
                self.assertEqual(_conjugate(v, p_a), p_c)
                self.assertEqual(multiply(p_a, clean), multiply(clean, p_a))
            # Identities hold on every physical column, including clean one.
            self.assertEqual(v, multiply(adjoint(ub), ua))

    def test_noncoordinate_affine_derivatives_remove_only_the_quadratic_part(self):
        # A small independent algebra fixture, not an asserted native query.
        # x=M y is an actual invertible CNOT coordinate change; U=M<e0,e1>.
        width = 5
        cnot_coordinates = [(0, 1), (1, 2), (1, 4)]

        def transform(value, inverse=False):
            gates = reversed(cnot_coordinates) if inverse else cnot_coordinates
            for control, target in gates:
                if value >> control & 1:
                    value ^= 1 << target
            return value

        def q(y):
            bits = [(y >> bit) & 1 for bit in range(width)]
            return (bits[0] * bits[2]) ^ (bits[1] * bits[3]) ^ (bits[0] * bits[1]) ^ bits[0]

        def h(y):
            return ((y >> 2) & 1) * ((y >> 3) & 1) * ((y >> 4) & 1)

        def g(x):
            y = transform(x, inverse=True)
            return q(y) ^ h(y)

        directions = [transform(1), transform(2)]
        self.assertEqual(_rank(directions), 2)
        self.assertTrue(all(direction.bit_count() > 1 for direction in directions))
        for x in range(1 << width):
            self.assertEqual(transform(transform(x), inverse=True), x)
        for subset in range(4):
            direction = (directions[0] if subset & 1 else 0) ^ (directions[1] if subset & 2 else 0)
            values = [g(x) ^ g(x ^ direction) for x in range(1 << width)]
            constant = values[0]
            linear = [values[1 << bit] ^ constant for bit in range(width)]
            for x, value in enumerate(values):
                affine = constant
                for bit in range(width):
                    affine ^= linear[bit] & (x >> bit & 1)
                self.assertEqual(value, affine)
        # Recover the ANF in the common coordinates: every degree >=3
        # coefficient avoids the two U coordinates, and q has no scalar.
        coefficients = [g(transform(y)) for y in range(1 << width)]
        for bit in range(width):
            for mask in range(1 << width):
                if mask >> bit & 1:
                    coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
        self.assertEqual({mask for mask, value in enumerate(coefficients) if value},
                         {1, 3, 5, 10, 28})
        for y in range(1 << width):
            self.assertEqual(g(transform(y)) ^ q(y), h(y))

    def test_clean_controlled_alignment_and_dirty_only_diagonalization(self):
        # Independent normalization fixture: clean c, active a, dirty d,e.
        # It supplies actual Clifford words, not a difficult Boolean query.
        width = 4
        clean, active, dirty1, dirty2 = range(width)
        z = [_tau(1 << (width + bit), width) for bit in range(width)]
        x = [_tau(1 << bit, width) for bit in range(width)]
        p_a = [multiply(z[active], x[dirty1]), x[dirty2]]
        q0 = multiply(z[dirty1], z[dirty2])
        # D=controlled-(Z_d Z_e), expanded only into H and CNOT.
        d_word = [gate for bit in (dirty1, dirty2)
                  for gate in [('H', bit), ('CX', clean, bit), ('H', bit)]]
        d = native_columns(width, d_word)
        g = z[active]
        v = multiply(d, g)
        self.assertEqual(multiply(d, d), eye(1 << width))
        normalized = multiply(d, v)
        self.assertEqual(normalized, g)
        dirty_word = [('H', dirty1), ('H', dirty2)]
        dirty = native_columns(width, dirty_word)
        for index, input_pauli in enumerate(p_a):
            output_pauli = multiply(input_pauli, z[clean])
            self.assertEqual(_conjugate(v, input_pauli), output_pauli)
            self.assertEqual(_conjugate(q0, input_pauli), scale(input_pauli, -1))
            self.assertEqual(_conjugate(d, output_pauli), input_pauli)
            self.assertEqual(_conjugate(normalized, input_pauli), input_pauli)
            diagonal = _conjugate(dirty, input_pauli)
            expected = multiply(z[active], z[dirty1]) if index == 0 else z[dirty2]
            self.assertEqual(diagonal, expected)
        self.assertEqual(_conjugate(dirty, normalized), g)
        for column in range(0, 1 << width, 2):
            for row in range(1 << width):
                self.assertEqual(d[row][column], CO if row == column else CZ)
                self.assertEqual(v[row][column], g[row][column])


if __name__ == '__main__':
    unittest.main()
