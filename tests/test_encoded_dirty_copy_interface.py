"""Complete encoded-copy echo and the fixed-Clifford companion screen.

The native p=1 fixture is a finite contract example, not a growing-size
compiler.  Its noncanonical copy is not an involution, so its actual
inverse and arbitrary helper inputs are checked explicitly.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CI, CO, CZ, adjoint, eye, multiply, native_columns
from tests.test_bent_query_baseline import _compact_toffoli
from tests.test_field_shear_intertwiner import (
    _native_x, _pencil_columns, _permutation_matrix,
)
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_t_depth import _word


def _p1_encoded_words():
    # Address x; carrier u,v; physical outputs y0,y1; arbitrary h and mask.
    query = _word(_compact_toffoli(0, 1, 2, 6))
    seed = [('T', 1)] + _native_x(1) + [('TDG', 1)]
    source = _inverse(query) + seed + query
    canonical = [('CX', 1, 3), ('CX', 2, 3), ('CX', 2, 4)]
    # W=controlled_h(R0) commutes with both R0 and R1=R0 X_v.
    # C=canonical W therefore has the copy action but is not canonical.
    gauge = [('T', 1), ('CX', 5, 1), ('TDG', 1)]
    copy = gauge + canonical
    echo = _inverse(source) + _inverse(copy) + source + copy
    return query, seed, source, copy, echo


def _cnot_label(label, source, target):
    x, z = label
    return (x ^ (((x >> source) & 1) << target),
            z ^ (((z >> target) & 1) << source))


def _label_word(label, word):
    for name, source, target in word:
        assert name == 'CX'
        label = _cnot_label(label, source, target)
    return label


def _p2_fixed_query(address):
    return [('CX', source, 2 + target)
            for source, column in enumerate(_pencil_columns(2, address))
            for target in range(2) if column >> target & 1]


class EncodedDirtyCopyInterfaceTests(unittest.TestCase):
    def test_noncanonical_native_copy_relation_and_actual_inverse(self):
        _, seed, _, copy, _ = _p1_encoded_words()
        # Restrict only the unused address and private query mask.  Every
        # carrier, output and gauge-helper input remains a full matrix.
        remap = lambda word: [(gate[0], *(wire - 1 for wire in gate[1:]))
                              for gate in word]
        r0 = native_columns(5, remap(seed))
        r1 = multiply(r0, native_columns(5, [('H', 1), ('S', 1),
                                            ('S', 1), ('H', 1)]))
        actual = native_columns(5, remap(copy))
        reverse = native_columns(5, remap(_inverse(copy)))
        self.assertEqual(reverse, adjoint(actual))
        self.assertNotEqual(actual, reverse)
        for address, reflection in enumerate((r0, r1)):
            target_x = native_columns(5, _native_x(2 + address))
            self.assertEqual(multiply(multiply(actual, reflection), reverse),
                             multiply(reflection, target_x))
        # Squaring the noncanonical C leaves a helper-controlled output X.
        # Reusing C in place of C^-1 would therefore fail on arbitrary h.
        square = native_columns(5, [('CX', 4, 2)])
        self.assertEqual(multiply(actual, actual), square)
        self.assertNotEqual(square, eye(32))
        self.assertEqual(multiply(r0, r0), eye(32))
        # R0 is not a Pauli: its off-diagonal entries have eighth-root phase.
        one_seed = native_columns(1, [('T', 0)] + _native_x(0) + [('TDG', 0)])
        self.assertNotIn(one_seed[1][0], (CO, -CO, CI, -CI))

    def test_native_indicator_every_column_return_and_ledger(self):
        query, seed, source, copy, echo = _p1_encoded_words()
        expected = _permutation_matrix([
            basis ^ (1 << (3 + (basis & 1))) for basis in range(128)])
        self.assertEqual(native_columns(7, echo), expected)
        self.assertEqual(native_columns(7, _inverse(echo)), adjoint(expected))
        # The target changes only the selected physical output.  Both
        # carrier wires and both arbitrary private inputs return literally.
        for basis in range(128):
            row = next(index for index, values in enumerate(expected)
                       if values[basis] != CZ)
            self.assertEqual((row ^ basis) & 0b1100111, 0)
        t_count = lambda word: sum(gate[0] in ('T', 'TDG') for gate in word)
        self.assertEqual(t_count(query), 8)
        self.assertEqual(t_count(seed), 2)
        self.assertEqual(t_count(copy), 2)
        self.assertEqual(t_count(echo), 40)
        for native in (echo, _inverse(echo)):
            self.assertEqual(t_count(native),
                             4 * t_count(query) + 2 * t_count(seed)
                             + 2 * t_count(copy))
            self.assertLessEqual(len(_bands(native)),
                                 4 * len(_bands(query)) + 2 * len(_bands(seed))
                                 + 2 * len(_bands(copy)))
            self.assertEqual({wire for gate in native for wire in gate[1:]},
                             set(range(7)))
            for band in _bands(native):
                self.assertEqual(len({gate[1] for gate in band}), len(band))
        self.assertEqual(len(source), 2 * len(query) + len(seed))

    def test_clifford_copy_support_with_genuinely_nonclifford_seed(self):
        # R=THT^-1=(X+Y)/2+Z/sqrt(2), on source u0.  Its conjugate
        # of X is neither diagonal nor antidiagonal, hence not a Pauli.
        seed = native_columns(1, [('TDG', 0), ('H', 0), ('T', 0)])
        x = native_columns(1, _native_x(0))
        image = multiply(multiply(seed, x), seed)
        self.assertTrue(all(value != CZ for row in image for value in row))
        labels = [(1, 0), (1, 1), (0, 1)]
        weights = [Fraction(1, 4), Fraction(1, 4), Fraction(1, 2)]
        # A fixed Clifford copy makes the X/Y output label affine in x.
        # It matches e_0,e_1,e_2 and misses e_3; the Z term always misses.
        copy = [('CX', 0, 4), ('CX', 2, 4), ('CX', 2, 5),
                ('CX', 3, 4), ('CX', 3, 6)]
        matched = []
        output_labels = []
        for address in range(4):
            query = _p2_fixed_query(address)
            outputs = []
            for label in labels:
                transformed = _label_word(_label_word(label, query), copy)
                outputs.append((transformed[0] >> 4, transformed[1] >> 4))
            output_labels.append(outputs)
            matched.append(sum(weight for weight, output in zip(weights, outputs)
                               if output == (1 << address, 0)))
        self.assertEqual([row[0] for row in output_labels],
                         [(1, 0), (2, 0), (4, 0), (7, 0)])
        self.assertEqual(matched, [Fraction(1, 2)] * 3 + [Fraction(0)])
        self.assertEqual(sum(matched) / 4, Fraction(3, 8))
        self.assertLessEqual(sum(matched) / 4, Fraction(3, 4))
        # At address three, output X signs +,+,+,- make actual and
        # target carrier reflections opposites. Pull both back through
        # the same Q_3 to compare the seed with its negative.
        actual_sector = seed
        target_sector = [[-entry for entry in row] for row in seed]
        difference = [[a - b for a, b in zip(row, other)]
                      for row, other in zip(actual_sector, target_sector)]
        self.assertEqual(multiply(adjoint(difference), difference),
                         [[CO * 4, CZ], [CZ, CO * 4]])

    def test_pauli_table_companions_have_only_affine_signs(self):
        # For R0=T^-1 X_u0 T, every Z-type Pauli is a symmetry;
        # its sign is precisely its Z_u0 label.  Field conjugation makes
        # these table signs affine, including all eight affine p=2 masks.
        tables = set()
        for z_mask in range(16):
            signs = tuple(_label_word((0, z_mask),
                                      list(reversed(_p2_fixed_query(address))))[1] & 1
                          for address in range(4))
            self.assertEqual(signs[0] ^ signs[1] ^ signs[2] ^ signs[3], 0)
            tables.add(signs)
        self.assertEqual(len(tables), 8)
        self.assertNotIn((0, 0, 0, 1), tables)


if __name__ == '__main__':
    unittest.main()
