"""Native copy components and exact screens for one field-query orbit.

The signed-Pauli pilot checks every coefficient of a prospective seed;
it is not a matrix simulation or a general one-band impossibility claim.
"""
import unittest

from tests.exact_matrices import adjoint, eye, multiply, native_columns, scale
from tests.test_field_shear_intertwiner import _native_x, _pencil_columns
from tests.test_joint_lookup_pauli_source import _bands, _inverse, _parity_basis
from tests.test_weak_pauli_orbit_compression import _binary_rank


def _one_band_copy(carrier=0, output=1):
    # C_H copies X and Z (and hence H), preserving Y.
    copy_h = ([('H', output)] + [('S', carrier)] * 3
              + [('CX', output, carrier), ('S', carrier), ('H', output)])
    # D sends Z_carrier to itself and Z_output to Z_carrier X_output.
    diagonalizer = [('CX', carrier, output), ('H', output)]
    return (copy_h + _inverse(diagonalizer)
            + [('T', carrier), ('TDG', output)] + diagonalizer)


def _pauli_word(width, x, z):
    word = []
    for wire in range(width):
        if (x >> wire) & 1:
            if (z >> wire) & 1:
                word += [('S', wire)] * 3 + _native_x(wire) + [('S', wire)]
            else:
                word += _native_x(wire)
        elif (z >> wire) & 1:
            word += [('S', wire)] * 2
    return word


def _field_label(address, x, z, active_blocks=1):
    """Hermitian Pauli action on one or two disjoint GF(4) shear blocks."""
    before = (x & z).bit_count()
    columns = ((0, 0), (1, 2), (2, 3), (3, 1))[address]
    for block in range(active_blocks):
        offset = 4 * block
        for source, column in enumerate(columns):
            if (x >> (offset + source)) & 1:
                x ^= column << (offset + 2)
            if ((z >> (offset + 2)) & column).bit_count() & 1:
                z ^= 1 << (offset + source)
    sign = ((before - (x & z).bit_count()) % 4) // 2
    return x, z, sign


def _magic_label(mask, x, z):
    """Adjoint action of B_mask = T^w X^mask T^-w."""
    return x, z ^ (x & mask), (z & ~x & mask).bit_count() & 1


def _coded_pilot(active_blocks=1):
    width, size = 8, 1 << 16
    columns = (255, 178, 212, 232)
    parents, parity, bad = list(range(size)), [0] * size, [False] * size
    sizes = [1] * size

    def find(label):
        if parents[label] != label:
            ancestor, bit = find(parents[label])
            parity[label] ^= bit
            parents[label] = ancestor
        return parents[label], parity[label]

    constraints = 0
    for address in range(4):
        for index, mask in enumerate(columns):
            for label in range(size):
                x, z, first = _field_label(address, label & 255, label >> width,
                                                  active_blocks)
                x, z, middle = _magic_label(mask, x, z)
                x, z, last = _field_label(address, x, z, active_blocks)
                other = x | (z << width)
                sign = first ^ middle ^ last ^ (address == index)
                left, lp = find(label)
                right, rp = find(other)
                if left == right:
                    bad[left] |= (lp ^ rp) != sign
                else:
                    if sizes[left] < sizes[right]:
                        left, right = right, left
                    parents[right] = left
                    parity[right] = lp ^ rp ^ sign
                    sizes[left] += sizes[right]
                    bad[left] |= bad[right]
                constraints += 1
    roots = {find(label)[0] for label in range(size)}
    return constraints, len(roots), sum(not bad[root] for root in roots)


class NativeCopyInteractionTests(unittest.TestCase):
    def test_literal_nonclifford_copy_in_one_physical_band(self):
        word = _one_band_copy()
        actual = native_columns(2, word)
        reverse = native_columns(2, _inverse(word))
        self.assertEqual(reverse, adjoint(actual))
        self.assertEqual(multiply(actual, reverse), eye(4))
        seed_word = [('TDG', 0), ('H', 0), ('T', 0)]
        seed = native_columns(2, seed_word)
        target_x = native_columns(2, _native_x(1))
        self.assertEqual(multiply(multiply(actual, seed), reverse),
                         multiply(seed, target_x))
        for axis in (_native_x(0), [('Z', 0)]):
            generator = native_columns(2, [('TDG', 0)] + axis + [('T', 0)])
            self.assertEqual(multiply(multiply(actual, generator), reverse),
                             multiply(generator, target_x))
        untouched = native_columns(2, [('TDG', 0)]
                                   + _pauli_word(2, 1, 1) + [('T', 0)])
        self.assertEqual(multiply(multiply(actual, untouched), reverse), untouched)
        self.assertEqual(len(_bands(word)), 1)
        self.assertEqual(len(_bands(word)[0]), 2)
        self.assertEqual({gate[1] for gate in _bands(word)[0]}, {0, 1})

    def test_coded_two_band_copy_and_dual_readout_all_columns(self):
        # Columns (011,101) extend to an invertible carrier basis.
        columns = (3, 5, 1)
        rows = [sum(((column >> bit) & 1) << index
                    for index, column in enumerate(columns)) for bit in range(3)]
        basis = [('H', bit) for bit in range(3)] + _parity_basis(rows)
        encoding = basis + [('T', bit) for bit in range(3)]
        copy = (_inverse(encoding) + [('CX', 0, 3), ('CX', 1, 4)] + encoding)
        actual = native_columns(5, copy)
        reverse = native_columns(5, _inverse(copy))
        self.assertEqual(reverse, adjoint(actual))
        self.assertEqual(len(_bands(copy)), 2)
        self.assertEqual([len(band) for band in _bands(copy)], [3, 3])
        for mask in range(4):
            carrier_mask = ((columns[0] if mask & 1 else 0)
                            ^ (columns[1] if mask & 2 else 0))
            z_output = [('Z', 3 + bit) for bit in range(2) if (mask >> bit) & 1]
            expected = native_columns(5, [('TDG', bit) for bit in range(3)]
                                      + _pauli_word(3, carrier_mask, 0)
                                      + [('T', bit) for bit in range(3)] + z_output)
            self.assertEqual(multiply(multiply(reverse, native_columns(5, z_output)), actual),
                             expected)
        for bit in range(2):
            reflection = native_columns(5, _inverse(encoding) + _native_x(bit) + encoding)
            target = native_columns(5, _native_x(3 + bit))
            self.assertEqual(multiply(multiply(actual, reflection), reverse),
                             multiply(reflection, target))

    def test_signed_label_action_against_literal_native_words(self):
        # Audit signs before using the sparse coefficient calculation.
        for address in range(4):
            query = [('CX', source, 2 + target)
                     for source, column in enumerate(_pencil_columns(2, address))
                     for target in range(2) if (column >> target) & 1]
            q = native_columns(4, query)
            for x, z in ((0, 0), (1, 4), (5, 5), (15, 15), (11, 14)):
                tx, tz, sign = _field_label(address, x, z)
                pauli = native_columns(4, _pauli_word(4, x, z))
                expected = native_columns(4, _pauli_word(4, tx, tz))
                self.assertEqual(multiply(multiply(q, pauli), adjoint(q)),
                                 scale(expected, -1 if sign else 1))
        magic = native_columns(1, [('TDG', 0)] + _native_x(0) + [('T', 0)])
        for x, z in ((0, 0), (1, 0), (0, 1), (1, 1)):
            tx, tz, sign = _magic_label(1, x, z)
            pauli = native_columns(1, _pauli_word(1, x, z))
            expected = native_columns(1, _pauli_word(1, tx, tz))
            self.assertEqual(multiply(multiply(magic, pauli), magic),
                             scale(expected, -1 if sign else 1))

    def test_coded_pilot_has_no_nonzero_seed_coefficients(self):
        self.assertEqual(_coded_pilot(), (1048576, 2760, 0))
        # Activate the second field core instead of leaving identity helpers.
        # This finite case lies outside the single-core analytic obstruction.
        self.assertEqual(_coded_pilot(active_blocks=2), (1048576, 1097, 0))

    def test_odd_masks_have_full_translation_rank_and_sparse_core_witnesses(self):
        for p in (1, 2, 3):
            m = 1 << p
            for mask in range(1 << m):
                if mask.bit_count() & 1:
                    shifts = [sum(((mask >> (address ^ shift)) & 1) << address
                                  for address in range(m)) for shift in range(m)]
                    self.assertEqual(_binary_rank(shifts), m)
        # The rank-(m-1) core map kills precisely the all-one mask.
        m = 8
        columns = [3 << bit for bit in range(m - 1)] + [129]
        self.assertEqual(_binary_rank(columns), m - 1)
        kernel, sparse_odd = [], []
        for mask in range(1 << m):
            image = 0
            for bit, column in enumerate(columns):
                if (mask >> bit) & 1:
                    image ^= column
            if image == 0:
                kernel.append(mask)
            if mask.bit_count() & 1 and image.bit_count() <= 2:
                sparse_odd.append(mask)
        self.assertEqual(kernel, [0, (1 << m) - 1])
        self.assertTrue(sparse_odd)


if __name__ == '__main__':
    unittest.main()
