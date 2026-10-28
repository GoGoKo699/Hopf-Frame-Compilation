"""Native evidence for the overlapping triple-wall one-band copy attempt.

The all-size zero-kernel argument and fixed-error floor are in
ONE_BAND_COPY_FEASIBILITY.md. These checks retain native phases, arbitrary
helper inputs, the actual inverse and the exact gap certificates. They
do not exclude general one-band copies.
"""
from itertools import product
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, adjoint, multiply, native_columns, scale,
)
from tests.test_field_shear_intertwiner import _pencil_columns
from tests.test_joint_lookup_pauli_source import _inverse


def _x_word(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


def _local_k(hadamard=False):
    word = [('CX', 2, 1), ('CX', 1, 2), ('CX', 1, 0), ('CX', 0, 2)]
    return word + ([('H', 0)] if hadamard else [])


def _local_sector(y, following, hadamard=False):
    k = _local_k(hadamard)
    beta = (y, following, y ^ following)
    return (_inverse(k)
            + [gate for wire, bit in enumerate(beta) if bit for gate in _x_word(wire)]
            + [('T', wire) for wire in range(3)]
            + [gate for wire in range(3) if y for gate in _x_word(wire)] + k)


def _pauli(width, x, z):
    phase = (CO, CI, -CO, -CI)[(x & z).bit_count() % 4]
    return [[phase * (-1 if (z & column).bit_count() & 1 else 1)
             if row == column ^ x else CZ for column in range(1 << width)]
            for row in range(1 << width)]


def _clifford_label(word, x, z, sign=0):
    """Literal Hermitian-Pauli adjoint for the H,S,CNOT native alphabet."""
    for name, *wires in word:
        wire = wires[0]
        a, b = x >> wire & 1, z >> wire & 1
        if name == 'H':
            sign ^= a & b
            if a != b:
                x ^= 1 << wire
                z ^= 1 << wire
        elif name in ('S', 'SDG'):
            sign ^= a & (b if name == 'S' else 1 ^ b)
            if a:
                z ^= 1 << wire
        elif name == 'CX':
            target = wires[1]
            c, d = x >> target & 1, z >> target & 1
            sign ^= a & d & (c ^ b ^ 1)
            if a:
                x ^= 1 << target
            if d:
                z ^= 1 << wire
        else:
            raise ValueError(name)
    return x, z, sign


def _native_walls(p):
    m = 1 << p
    k, alpha, beta = [], [], []
    for j in range(m):
        wires = (j, m + 2 * j, m + 2 * j + 1)
        k.extend((name, *(wires[wire] for wire in local))
                 for name, *local in _local_k())
        for local in range(3):
            alpha.append(('CX', 3 * m + j, wires[local]))
        beta.extend([('CX', 3 * m + j, wires[0]),
                     ('CX', 3 * m + (j + 1) % m, wires[1]),
                     ('CX', 3 * m + j, wires[2]),
                     ('CX', 3 * m + (j + 1) % m, wires[2])])
    k += [('H', wire) for wire in range(m // 2)]
    output_h = [('H', 3 * m + j) for j in range(m)]
    # Lists are chronological; C=F D G executes G, D, F.
    return alpha + k + output_h, output_h + _inverse(k) + beta


def _rank(values):
    pivots = {}
    for value in values:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def _coefficient_probe(m, t, w, z):
    previous = (t - 1) % m
    wt, wp = w >> t & 1, w >> previous & 1
    x = (wt << (2 * t + 1)) | (wp << (2 * previous))
    zmask = (1 << (2 * t)) | (1 << (2 * previous + 1))
    sign = (z >> t & 1) ^ wt ^ wp
    return x, zmask, sign


class JointWallCopyFamilyTests(unittest.TestCase):
    def test_literal_native_sector_defect_on_every_surviving_local_pauli(self):
        # Compute C_y C_0^dagger from the actual T/TDG words. There is no
        # use of the claimed Clifford defect simplification in this check.
        for hadamard in (False, True):
            coordinate = native_columns(3, [('H', 0)] if hadamard else [])
            zero = _local_sector(0, 0, hadamard)
            for y, following in product(range(2), repeat=2):
                word = _inverse(zero) + _local_sector(y, following, hadamard)
                defect = native_columns(3, word)
                defect_inverse = native_columns(3, _inverse(word))
                self.assertEqual(defect_inverse, adjoint(defect))
                for w, z, r1, r2 in product(range(2), repeat=4):
                    incoming = _pauli(3, w, z | (r1 << 1) | (r2 << 2))
                    incoming = multiply(multiply(coordinate, incoming), adjoint(coordinate))
                    actual = multiply(multiply(defect, incoming), defect_inverse)
                    out1, out2 = r1 ^ (w & following), r2 ^ (w & y)
                    exponent = (y * z + y * r1 + following * r2
                                + w * (y + following + y * following)) & 1
                    expected = _pauli(3, w, z | (out1 << 1) | (out2 << 2))
                    expected = multiply(multiply(coordinate, expected), adjoint(coordinate))
                    self.assertEqual(actual, scale(expected, -1 if exponent else 1))

    def test_native_walls_axes_independence_and_paid_ledger(self):
        for p in (1, 2, 3):
            m, width = 1 << p, 4 << p
            f, g = _native_walls(p)
            self.assertEqual(len(f) + len(g), 18 * m)
            native = g + [('T', wire) for wire in range(3 * m)] + f
            self.assertEqual(sum(gate[0] == 'T' for gate in native), 3 * m)
            self.assertEqual(len(native), 21 * m)
            for wall, use_alpha in ((f, True), (_inverse(g), False)):
                labels = []
                for j in range(m):
                    wires = (j, m + 2 * j, m + 2 * j + 1)
                    helpers = (1 << wires[1], 1 << wires[2],
                               (1 << wires[1]) | (1 << wires[2]))
                    output = (1 << j, 1 << ((j + 1) % m),
                              (1 << j) | (1 << ((j + 1) % m)))
                    for local, wire in enumerate(wires):
                        x, z = (1 << j, 0) if j < m // 2 else (0, 1 << j)
                        z |= helpers[local]
                        x |= (1 << j if use_alpha else output[local]) << (3 * m)
                        self.assertEqual(_clifford_label(wall, 0, 1 << wire), (x, z, 0))
                        labels.append(x | (z << width))
                self.assertEqual(_rank(labels), 3 * m)
                for left, right in product(labels, repeat=2):
                    mask = (1 << width) - 1
                    self.assertEqual(((left & mask & (right >> width)).bit_count()
                                      + (right & mask & (left >> width)).bit_count()) % 2, 0)

    def test_tensor_prefilter_keeps_every_core_label_but_no_helper_x(self):
        m = 4
        # z labels do not enter these symplectic pairings. Exhaust every
        # core X label and all helper X labels in the fixed a coordinates.
        for core_x in range(1 << m):
            for helper_x in range(1 << (2 * m)):
                tensor = 0
                for j in range(m):
                    a = core_x >> j & 1
                    h1, h2 = helper_x >> (2 * j) & 1, helper_x >> (2 * j + 1) & 1
                    activations = (a ^ h1, a ^ h2, a ^ h1 ^ h2)
                    beta = (1 << j, 1 << ((j + 1) % m),
                            (1 << j) | (1 << ((j + 1) % m)))
                    for active, row in zip(activations, beta):
                        if active:
                            tensor ^= row << (m * j)
                self.assertEqual(tensor == 0, helper_x == 0)

    def test_bad_local_walsh_blocks_have_exact_half_squared_norm(self):
        # In l coordinates, beta=(e0,e1,e0+e1), alpha=(e0,e0,e0).
        # Fix a native-D X block and every spectator/output coordinate.
        # Each bad triple has one or two active axes; its beta restriction
        # is injective, so each nonempty source slice is a singleton.
        beta = (1, 2, 3)
        checked = 0
        for active in range(1, 7):
            dimension = active.bit_count()
            for spectator in range(8):
                if spectator & active:
                    continue
                block = [q for q in range(8) if q & (7 ^ active) == spectator]
                for incoming_output, target_output in product(range(4), repeat=2):
                    source = [q for q in block
                              if ((beta[0] if q & 1 else 0)
                                  ^ (beta[1] if q & 2 else 0)
                                  ^ (beta[2] if q & 4 else 0)) == incoming_output]
                    target = [q for q in block
                              if q.bit_count() % 2 == (incoming_output ^ target_output)]
                    if not source or not target:
                        continue
                    self.assertEqual(len(source), 1)
                    # Native Walsh entries, including either orientation's
                    # row/column signs, have squared magnitude 2^-dimension.
                    self.assertEqual(2 * len(target), 1 << dimension)
                    checked += 1
        self.assertEqual(checked, 96)

    def test_four_address_literal_identity_contradicts_each_lookup_character(self):
        p, m, k = 2, 4, 2
        base = [('H', wire) for wire in range(k)]
        for w, z in product(range(1 << m), repeat=2):
            orbit = []
            for address in range(1 << p):
                query = [('CX', source, k + target)
                         for source, column in enumerate(_pencil_columns(p, address))
                         for target in range(k) if column >> target & 1]
                wx, zx, _ = _clifford_label(base + query + base, w, z)
                orbit.append((wx, zx))
            for target in range(m):
                labels = [_coefficient_probe(m, target, *u) for u in orbit]
                x, zmask, sign = 0, 0, 0
                for a, b, phase in labels:
                    # All these individual Paulis have disjoint X and Z
                    # supports, so Hermitian and normal phases agree.
                    self.assertEqual(a & b, 0)
                    sign ^= phase ^ ((zmask & a).bit_count() & 1)
                    x ^= a
                    zmask ^= b
                self.assertEqual((x, zmask, sign), (0, 0, 0))
                self.assertEqual(sum(address == target for address in range(4)) % 2, 1)

        # Product +I restricts simultaneous eigenvalues to even parity.
        # A one-hot target has odd parity: its squared residual sum is
        # at least four, giving ||good||^2 <= 4*eta^2 from four errors 2*eta.
        for target in range(4):
            wanted = [-1 if address == target else 1 for address in range(4)]
            residuals = []
            for signs in product((-1, 1), repeat=4):
                if signs[0] * signs[1] * signs[2] * signs[3] == 1:
                    residuals.append(sum((a - b) ** 2 for a, b in zip(signs, wanted)))
            self.assertEqual(min(residuals), 4)


if __name__ == '__main__':
    unittest.main()
