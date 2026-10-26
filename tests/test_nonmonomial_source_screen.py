"""Exact receipts for the Hadamard-pair source's cubic-sign screen.

The arbitrary-Clifford/helper theorem is analytic. These fixtures check its
complete first image, support geometry, literal phases and weighted escape;
no asymptotic source or unrestricted depth conclusion is inferred.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import CI, CO, CZ, native_columns
from tests.test_common_probe_source import _native
from tests.test_joint_address_helper_source import _integer_word
from tests.test_joint_lookup_pauli_source import _inverse


X1, X2, X3, X4, U, V, ROOT = range(7)


def _bit(value, wire):
    return (value >> wire) & 1


def _symplectic(left, right):
    x, z = left
    other_x, other_z = right
    return ((x & other_z).bit_count() + (z & other_x).bit_count()) & 1


def _first_image():
    """Independent exact expansion of equation (13), with all CZ signs."""
    zr = 1 << ROOT
    result = {(0, zr | (1 << U) | (1 << V)): Fraction(1, 2)}
    even = set(result)
    odd = set()
    for s, t in product(range(2), repeat=2):
        sign = (-1) ** (s * t)
        first = (1 << U, zr | (1 << V) | (s << X1) | (t << X2))
        second = (1 << V, zr | (1 << U) | (s << X3) | (t << X4))
        result[first] = result[second] = Fraction(sign, 4)
        odd.update((first, second))
    for bits in product(range(2), repeat=4):
        label = ((1 << U) | (1 << V),
                 zr | sum(bit << wire for wire, bit in enumerate(bits)))
        result[label] = Fraction((-1) ** (bits[0] * bits[1] + bits[2] * bits[3]), 8)
        even.add(label)
    return result, even, odd


def _walsh(values):
    out = list(values)
    half = 1
    while half < len(out):
        for start in range(0, len(out), 2 * half):
            for offset in range(half):
                i, j = start + offset, start + offset + half
                out[i], out[j] = out[i] + out[j], out[i] - out[j]
        half *= 2
    return [Fraction(value, len(out)) for value in out]


def _phase(value, terms):
    return sum(all(_bit(value, wire) for wire in term) for term in terms) & 1


class NonmonomialSourceScreenTests(unittest.TestCase):
    def test_first_image_on_every_column_including_arbitrary_helpers(self):
        expansion, _, _ = _first_image()
        d2 = [('CCZ', U, X1, X2), ('CCZ', V, X3, X4)]
        word = d2 + [('Z', ROOT), ('H', U), ('H', V)] + d2
        for basis in range(1 << 7):
            actual, denominator = _integer_word(word, basis)
            actual = {row: Fraction(value, denominator) for row, value in actual.items()}
            expected = {}
            for (x, z), coefficient in expansion.items():
                self.assertEqual(x & z, 0)  # These particular terms contain no Y.
                row = basis ^ x
                value = coefficient * (-1) ** ((basis & z).bit_count() & 1)
                expected[row] = expected.get(row, 0) + value
            expected = {row: value for row, value in expected.items() if value}
            self.assertEqual(actual, expected)

    def test_two_support_classes_and_outlying_singleton(self):
        expansion, even, odd = _first_image()
        self.assertEqual((len(expansion), len(even), len(odd)), (25, 17, 8))
        for labels in (even, odd):
            self.assertTrue(all(_symplectic(p, q) == 0 for p in labels for q in labels))
        self.assertTrue(all(_symplectic(p, q) == 1 for p in even for q in odd))
        singleton = (0, (1 << ROOT) | (1 << U) | (1 << V))
        flat = even - {singleton}
        origin = next(iter(flat))
        differences = {(x ^ origin[0], z ^ origin[1]) for x, z in flat}
        self.assertEqual(len(differences), 16)
        self.assertTrue(all((x ^ xx, z ^ zz) in differences
                            for x, z in differences for xx, zz in differences))
        self.assertNotIn(singleton, flat)
        self.assertEqual(sum(abs(value) for value in expansion.values()), Fraction(9, 2))

    def test_cubic_derivatives_have_flat_commuting_affine_walsh_support(self):
        width = 5
        polynomials = (
            [(0, 1, 2)],
            [(0, 1, 2), (2, 3, 4)],
            [(0,), (1, 3), (0, 1, 4), (0, 2, 3), (1, 2, 4)],
        )
        observed_sizes = set()
        for terms in polynomials:
            for displacement in range(1 << width):
                signs = [(-1) ** (_phase(z, terms) ^ _phase(z ^ displacement, terms))
                         for z in range(1 << width)]
                coefficients = _walsh(signs)
                support = {s for s, value in enumerate(coefficients) if value}
                observed_sizes.add(len(support))
                self.assertIn(len(support), (1, 4, 16))
                self.assertEqual(len({abs(coefficients[s]) for s in support}), 1)
                self.assertEqual(next(abs(coefficients[s]) for s in support) ** 2,
                                 Fraction(1, len(support)))
                origin = next(iter(support))
                directions = {s ^ origin for s in support}
                self.assertTrue(all(s ^ t in directions for s in directions for t in directions))
                self.assertTrue(all((s & displacement).bit_count() % 2 == 0 for s in support))
                # Permit overlapping X/Z masks, hence literal Hermitian Y factors.
                for z_mask in (0, displacement, 19):
                    labels = [(displacement, z_mask ^ s) for s in support]
                    self.assertTrue(all(_symplectic(p, q) == 0 for p in labels for q in labels))
        self.assertEqual(observed_sizes, {1, 4, 16})

    def test_target_coefficient_uses_all_carrier_inputs(self):
        trace = 0
        accepted_addresses = 0
        # Seven independent bits x1,x2,x3,x4,c1,c2,d; Z_r cancels the root sign.
        for value in range(128):
            a = _bit(value, 0) * _bit(value, 1)
            b = _bit(value, 2) * _bit(value, 3)
            c1, c2, d = (_bit(value, wire) for wire in (4, 5, 6))
            trace += (-1) ** (a * c1 + b * c2 + a * b * d)
        for address in range(16):
            a = _bit(address, 0) * _bit(address, 1)
            b = _bit(address, 2) * _bit(address, 3)
            accepted_addresses += a == b == 0
        self.assertEqual(Fraction(trace, 128), Fraction(9, 16))
        self.assertEqual(accepted_addresses, 9)
        self.assertEqual(Fraction(trace, 128) - Fraction(1, 2), Fraction(1, 16))

    def test_weighted_outer_phases_escape_the_half_coefficient_screen(self):
        # Literal pullback C†D†(Xu Xv)DC equals Hu Hv, including all phases.
        clifford = [('H', 0), ('S', 0), ('H', 0),
                    ('H', 1), ('S', 1), ('H', 1)]
        diagonal = [('T', 0), ('T', 1)]
        x = [('H', 0), ('S', 0), ('S', 0), ('H', 0),
             ('H', 1), ('S', 1), ('S', 1), ('H', 1)]
        pullback = clifford + diagonal + x + _inverse(diagonal) + _inverse(clifford)
        self.assertEqual(native_columns(2, pullback), native_columns(2, [('H', 0), ('H', 1)]))
        expansion, _, _ = _first_image()
        coefficient = Fraction(0)
        for x_mask, z_mask in (
            (0, (1 << U) | (1 << V)),
            (1 << U, 1 << V),
            (1 << V, 1 << U),
            ((1 << U) | (1 << V), 0),
        ):
            coefficient += expansion[x_mask, z_mask | (1 << ROOT)] / 2
        self.assertEqual(coefficient, Fraction(9, 16))
        self.assertGreater(coefficient, Fraction(1, 2))
        # Independently evaluate the full weighted word on every active column.
        outer = [('H', U), ('S', U), ('H', U),
                 ('H', V), ('S', V), ('H', V), ('T', U), ('T', V)]
        k_word = [('S', ROOT), ('S', ROOT), ('H', U), ('H', V)]

        def first_phase(state):
            return {row: amplitude * (-1) ** (
                _bit(row, U) * _bit(row, X1) * _bit(row, X2)
                + _bit(row, V) * _bit(row, X3) * _bit(row, X4))
                    for row, amplitude in state.items()}

        traces = [CZ, CZ]
        for basis in range(128):
            state = _native(_inverse(outer), {basis: CO})
            state = first_phase(state)
            state = _native(k_word, state)
            state = first_phase(state)
            actual = _native(outer, state)
            a = _bit(basis, X1) * _bit(basis, X2)
            b = _bit(basis, X3) * _bit(basis, X4)
            phase = CO * (-1) ** _bit(basis, ROOT)
            for product_bit, wire in ((a, U), (b, V)):
                if product_bit:
                    phase *= -CI * (-1) ** _bit(basis, wire)
            output = basis ^ (1 << U) ^ (1 << V)
            self.assertEqual(actual, {output: phase})
            # Literal Hermitian Paulis: Zr Xu Xv and Zr Yu Xv.
            traces[0] += phase * (-1) ** _bit(output, ROOT)
            traces[1] += phase * CI * (-1) ** (_bit(output, ROOT) + _bit(output, U))
        self.assertEqual(traces[0] * Fraction(1, 128), CO * Fraction(9, 16))
        self.assertEqual(traces[1] * Fraction(1, 128), CO * Fraction(-3, 16))
        self.assertEqual(_symplectic(
            ((1 << U) | (1 << V), 1 << ROOT),
            ((1 << U) | (1 << V), (1 << ROOT) | (1 << U))), 1)


if __name__ == '__main__':
    unittest.main()
