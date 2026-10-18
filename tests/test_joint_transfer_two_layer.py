"""Exact bounded evidence for the unrestricted two-layer joint-row bound.

The proof uses two entries of one orthonormal physical Pauli-transfer row,
not two independently chosen alphabet values.  Exact finite enumerations
below have explicit tail bounds; no floating comparison certifies a gap.
Sparse permutation actions audit the AND target on every occupied input.
Returned native encoders also test the stronger shared-axis single alphabet;
the cubic-bent scalar boundary is explicitly not a circuit realization.
"""
from fractions import Fraction
from itertools import product
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, Q2, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_dirty_indicator_depth import _parity_basis, _toffoli_two_layers
from tests.test_t_depth import _word
from tests.test_two_layer_galois_return import (
    _conjugate, _controlled_h_word, _inverse,
)


ZERO, ONE = Q2(), Q2(1)


def _alphabet(last_exponent):
    result = [ZERO]
    for exponent in range(last_exponent + 1):
        value = (Q2(Fraction(1, 2 ** (exponent // 2)))
                 if exponent % 2 == 0
                 else Q2(0, Fraction(1, 2 ** ((exponent + 1) // 2))))
        result.extend((value, -value))
    return result


def _row_admissible_sums(last_exponent):
    alphabet = _alphabet(last_exponent)
    return {left + right for left, right in product(alphabet, repeat=2)
            if left * left + right * right <= ONE}


def _pauli(width, x_mask, z_mask):
    """Hermitian i^(x.z) X^x Z^z with an exact signed permutation matrix."""
    size = 1 << width
    matrix = [[CZ for _ in range(size)] for _ in range(size)]
    phase = (CO, CI, -CO, -CI)[(x_mask & z_mask).bit_count() % 4]
    for column in range(size):
        matrix[column ^ x_mask][column] = (
            -phase if (column & z_mask).bit_count() % 2 else phase)
    return matrix


def _transfer_row(unitary, output_pauli):
    size = len(unitary)
    width = size.bit_length() - 1
    image = multiply(multiply(adjoint(unitary), output_pauli), unitary)
    row = {}
    for x_mask, z_mask in product(range(size), repeat=2):
        phase = (CO, CI, -CO, -CI)[(x_mask & z_mask).bit_count() % 4]
        coefficient = sum((
            (-phase if (basis & z_mask).bit_count() % 2 else phase)
            * image[basis][basis ^ x_mask]
            for basis in range(size)), CZ) / size
        if coefficient.imag != ZERO:
            raise AssertionError("Hermitian Pauli-transfer entry is not real")
        row[x_mask, z_mask] = coefficient.real
    return row, image


def _and_query_images(address_bits, dirty_bits=0):
    """Full permutation, with an arbitrary output and untouched dirty wires."""
    target = 1 << address_bits
    mask = target - 1
    return [basis ^ target if (basis & mask) == mask else basis
            for basis in range(1 << (address_bits + 1 + dirty_bits))]


def _binary_span(labels):
    span = {0}
    for label in labels:
        span |= {value ^ label for value in span}
    return span


def _normalized_trace(matrix):
    return sum((matrix[j][j] for j in range(len(matrix))), CZ) / len(matrix)


def _integer_walsh(values):
    spectrum = [(-1 if value else 1) for value in values]
    stride = 1
    while stride < len(spectrum):
        for start in range(0, len(spectrum), 2 * stride):
            for offset in range(stride):
                left, right = spectrum[start + offset], spectrum[start + offset + stride]
                spectrum[start + offset] = left + right
                spectrum[start + offset + stride] = left - right
        stride *= 2
    return spectrum


class JointTransferTwoLayerTests(unittest.TestCase):
    def test_native_row_norm_and_initialized_projector_sum(self):
        # Qubit 0 is logical and qubit 1 is the sole initialized input.
        # These words need not return the clean qubit; the row constraint and
        # initialized-projector identity apply to every physical unitary.
        words = [
            [('H', 0), ('T', 0), ('CX', 0, 1), ('H', 1),
             ('TDG', 1), ('H', 0)],
            [('H', 1), ('T', 1), ('CX', 1, 0), ('S', 0),
             ('H', 0), ('T', 0), ('CX', 0, 1)],
        ]
        alphabet = set(_alphabet(4))  # q=2 gives j <= 2q.
        for word in words:
            unitary = native_columns(2, word)
            for output_label in ((0, 1), (1, 0), (1, 1)):
                row, image = _transfer_row(unitary, _pauli(2, *output_label))
                self.assertTrue(set(row.values()) <= alphabet)
                self.assertEqual(sum((value * value for value in row.values()), ZERO), ONE)
                for logical_x, logical_z in product(range(2), repeat=2):
                    c_i = row[logical_x, logical_z]
                    c_z = row[logical_x, logical_z | 2]
                    self.assertLessEqual(c_i * c_i + c_z * c_z, ONE)
                    compressed = [[image[i][j] for j in (0, 1)] for i in (0, 1)]
                    logical = _pauli(1, logical_x, logical_z)
                    product_matrix = multiply(logical, compressed)
                    initialized = sum((product_matrix[i][i] for i in range(2)), CZ) / 2
                    self.assertEqual(initialized.imag, ZERO)
                    self.assertEqual(initialized.real, c_i + c_z)

    def test_raw_alphabet_false_positive_and_exact_rational_gap(self):
        for arity in range(4, 17):
            small = Q2(Fraction(1, 2 ** (arity - 1)))
            target = ONE - small
            # A+A alone accepts the target; a common physical row cannot.
            self.assertEqual(ONE + (-small), target)
            self.assertGreater(ONE + small * small, ONE)

        admissible = _row_admissible_sums(12)
        self.assertFalse(any(value.b == 0 and Q2(Fraction(3, 4)) < value < ONE
                             for value in admissible))
        self.assertIn(Q2(Fraction(3, 4)), admissible)
        self.assertIn(ONE, admissible)
        # Any omitted entry has magnitude < 1/64. If its partner is +/-1,
        # row norm excludes it. Otherwise the partner has magnitude <=1/sqrt2,
        # so their sum cannot reach 3/4. This closes the infinite tail.
        self.assertLess(Q2(0, Fraction(1, 2)) + Fraction(1, 64), Fraction(3, 4))

    def test_closest_joint_alphabet_points_and_all_arity_gap_formulas(self):
        admissible = _row_admissible_sums(6)
        targets = {
            4: (Q2(0, Fraction(5, 8)), Q2(Fraction(-7, 16), Fraction(5, 16))),
            5: (Q2(Fraction(1, 4), Fraction(1, 2)),
                Q2(Fraction(-11, 32), Fraction(1, 4))),
            6: (Q2(Fraction(1, 4), Fraction(1, 2)),
                Q2(Fraction(23, 64), Fraction(-1, 4))),
        }
        tail_maximum = Q2(0, Fraction(1, 2)) + Fraction(1, 8)
        for arity in range(4, 65):
            target = ONE - Fraction(1, 2 ** (arity - 1))
            if arity in targets:
                nearest, operator_floor = targets[arity]
            else:
                nearest, operator_floor = ONE, Q2(Fraction(1, 2 ** arity))
            distance = min(abs(target - value) for value in admissible)
            self.assertEqual(distance, abs(target - nearest))
            self.assertEqual(distance / 2, operator_floor)
            self.assertGreater(operator_floor, ZERO)
            # Beyond j=6, a nonzero omitted entry is <1/8. The +/-1
            # partner is forbidden, so every omitted sum is <=tail_maximum.
            self.assertLess(tail_maximum, target - distance)

        # For all a>=7, b>=63/64. This exact inequality makes 1 nearer
        # than the largest admissible point below 1; points above 1 cannot
        # be nearer than 1. Thus the test has an analytic unbounded-arity tail.
        below_one = max(value for value in admissible if value < ONE)
        self.assertEqual(below_one, Q2(Fraction(1, 4), Fraction(1, 2)))
        self.assertGreater(Q2(Fraction(63, 64)), (ONE + below_one) / 2)
        self.assertGreater(targets[4][1], Fraction(1, 256))

    def test_sparse_full_and_query_transfer_coefficients(self):
        for arity in range(1, 7):
            target_bit = 1 << arity
            mask = target_bit - 1
            images = _and_query_images(arity, dirty_bits=1)
            self.assertEqual(sorted(images), list(range(len(images))))
            self.assertTrue(all(images[images[basis]] == basis for basis in range(len(images))))
            for basis, image in enumerate(images):
                self.assertEqual((basis ^ image) & ~target_bit, 0)
            # Row P=Z_target and column Q=Z_target Z^frequency is W_f.
            # This trace uses every output and dirty input, not |0> columns.
            for frequency in range(1 << arity):
                trace = sum((-1 if (((image >> arity) ^ (basis >> arity)
                                      ^ (basis & frequency).bit_count()) & 1) else 1)
                            for basis, image in enumerate(images))
                coefficient = Fraction(trace, len(images))
                expected = (Fraction(1) if frequency == 0 else Fraction(0))
                expected -= Fraction(2, 1 << arity) * (
                    -1 if (mask & frequency).bit_count() % 2 else 1)
                self.assertEqual(coefficient, expected)
            self.assertEqual(Fraction(sum(1 if ((basis ^ image) & target_bit) == 0 else -1
                                          for basis, image in enumerate(images)), len(images)),
                             Fraction(1) - Fraction(1, 1 << (arity - 1)))

    def test_walsh_affine_distance_corollary_and_sparse_boundary(self):
        # Exhaust every three-input Boolean table, including both constants
        # and all affine functions. This audits the Walsh/distance conversion;
        # passing this necessary condition does not certify a native circuit.
        for table in range(1 << 8):
            values = [(table >> x) & 1 for x in range(8)]
            distances = []
            coefficients = []
            for frequency in range(8):
                linear = [(x & frequency).bit_count() % 2 for x in range(8)]
                distance = sum(f != ell for f, ell in zip(values, linear))
                complement_distance = sum(f == ell for f, ell in zip(values, linear))
                coefficient = Fraction(sum((-1) ** (f ^ ell)
                                           for f, ell in zip(values, linear)), 8)
                self.assertEqual(coefficient, 1 - Fraction(2 * distance, 8))
                self.assertEqual(-coefficient, 1 - Fraction(2 * complement_distance, 8))
                distances.extend((distance, complement_distance))
                coefficients.append(coefficient)
            is_affine = min(distances) == 0
            self.assertEqual(is_affine, max(abs(w) for w in coefficients) == 1)
            if is_affine:
                self.assertEqual(sum(w != 0 for w in coefficients), 1)
            else:
                self.assertGreaterEqual(Fraction(min(distances), 8), Fraction(1, 8))
                self.assertLessEqual(max(abs(w) for w in coefficients), Fraction(3, 4))

        # One marked input is exactly at the permitted arithmetic boundary
        # for a=3, and strictly in the forbidden gap for every a>=4.
        # Complementing the table checks the negative gap and dense predicates.
        for arity in range(3, 17):
            mean = Fraction(1) - Fraction(2, 1 << arity)
            distance_to_constant = Fraction(1, 1 << arity)
            self.assertEqual(distance_to_constant, (1 - mean) / 2)
            if arity == 3:
                self.assertEqual(mean, Fraction(3, 4))
                self.assertEqual(distance_to_constant, Fraction(1, 8))
            else:
                self.assertTrue(Fraction(3, 4) < mean < 1)
                self.assertTrue(-1 < -mean < Fraction(-3, 4))
                self.assertLess(distance_to_constant, Fraction(1, 8))

    def test_retained_two_layer_dirty_toffoli_remains_allowed(self):
        schedule = _toffoli_two_layers(0, 1, 2, 3)
        self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 2)
        self.assertEqual(sum(len(gates) for kind, gates in schedule if kind == 'T'), 8)
        unitary = native_columns(4, _word(schedule))
        images = _and_query_images(2, dirty_bits=1)
        self.assertTrue(all(unitary[row][column] == (CO if row == images[column] else CZ)
                            for row in range(16) for column in range(16)))
        # The arbitrary helper may also be designated as the one clean bit.
        # The target-Z row has c_I=1/2 and c_Z=0, consistent with both bounds.
        c_i = Fraction(sum(1 if ((basis ^ image) & 4) == 0 else -1
                           for basis, image in enumerate(images)), 16)
        c_z = Fraction(sum((-1 if basis & 8 else 1)
                           * (1 if ((basis ^ image) & 4) == 0 else -1)
                           for basis, image in enumerate(images)), 16)
        self.assertEqual((c_i, c_z), (Fraction(1, 2), Fraction(0)))
        self.assertLessEqual(c_i * c_i + c_z * c_z, 1)

    def test_nonpauli_shared_midpoint_cross_encoder_axes_and_trace_cancellation(self):
        # Clean wire 0 is encoded into X0 X1 X2 X3, then all four native
        # T gates act. The second encoder uses a genuinely different active
        # basis, not just a permutation: Zi maps to minus product_{j!=i} Zj.
        width, size = 4, 16
        clifford = [('H', 0)] + [('CX', 0, j) for j in range(1, width)]
        band = [('T', j) for j in range(width)]
        axis_change = _parity_basis(list(range(width)), [15 ^ (1 << j) for j in range(width)])
        # x -> (I+J_4)x + 1111; expand X as H S S H in the native alphabet.
        axis_change += [gate for j in range(width)
                        for gate in [('H', j), ('S', j), ('S', j), ('H', j)]]
        u_a = native_columns(width, clifford + band)
        u_b = native_columns(width, clifford + band + axis_change)
        w = native_columns(width, axis_change)
        self.assertEqual(w, adjoint(w))
        clean_z = _pauli(width, 0, 1)
        shared = _conjugate(u_a, clean_z)
        self.assertEqual(shared, _conjugate(u_b, clean_z))
        self.assertEqual(_conjugate(shared, clean_z), scale(clean_z, -1))

        first_axes, second_axes, displacements = [], [], []
        for j in range(width):
            bit = 1 << j
            # These eight exact generator equations identify F_S. It sends
            # Xj to Yj and Zj to -Zj, so im(F_S-I)=span{Z0,...,Z3}.
            self.assertEqual(_conjugate(shared, _pauli(width, bit, 0)),
                             _pauli(width, bit, bit))
            self.assertEqual(_conjugate(shared, _pauli(width, 0, bit)),
                             scale(_pauli(width, 0, bit), -1))
            self.assertEqual(_conjugate(w, _pauli(width, 0, bit)),
                             scale(_pauli(width, 0, 15 ^ bit), -1))
            first_axes.append(bit << width)
            second_axes.append((15 ^ bit) << width)
            displacements.extend((bit << width, 0))
        self.assertNotEqual(set(first_axes), set(second_axes))
        self.assertEqual(_binary_span(first_axes), _binary_span(second_axes))
        self.assertEqual(_binary_span(first_axes), _binary_span(displacements))

        # W is its actual inverse, as checked on its full matrix above.
        word = clifford + band + axis_change + _inverse(band) + _inverse(clifford)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 8)
        self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in word), 38)
        actual = native_columns(width, word)
        self.assertEqual(actual, multiply(adjoint(u_b), u_a))
        self.assertEqual(actual, adjoint(actual))
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        self.assertTrue(all(actual[row][column] == CZ for row in range(size)
                            for column in range(size) if (row ^ column) & 1))
        occupied = [[actual[2 * row][2 * column] for column in range(8)] for row in range(8)]
        expected = eye(8)
        expected[0][0] = -CO
        self.assertEqual(occupied, expected)
        # This is a cubic phase on three logical wires. A target-H conversion
        # retains a control-dependent sign; it is not a three-address cubic XOR.
        # The full literal matrix above is the actual returned native fixture.
        alphabet = set(_alphabet(6))
        for output in ((1, 0), (0, 1), (3, 0), (1, 1)):
            row, _ = _transfer_row(occupied, _pauli(3, *output))
            self.assertTrue(set(row.values()) <= alphabet)

        surviving_terms = []
        for logical_x, logical_z in ((1, 0), (0, 1), (3, 0)):
            lifted = _pauli(width, logical_x << 1, logical_z << 1)
            a = _conjugate(u_a, lifted)
            b = _conjugate(u_b, lifted)
            g = multiply(b, a)
            gs = multiply(g, shared)
            trace_g, trace_gs = _normalized_trace(g), _normalized_trace(gs)
            self.assertTrue(trace_g == CZ or trace_gs == CZ)
            epsilon = 1 if _conjugate(clean_z, g) == g else -1
            self.assertEqual(_conjugate(clean_z, g), scale(g, epsilon))
            self.assertEqual(_conjugate(clean_z, gs), scale(gs, -epsilon))
            self.assertEqual(trace_g if epsilon == -1 else trace_gs, CZ)
            logical = _pauli(3, logical_x, logical_z)
            direct = _normalized_trace(multiply(logical, _conjugate(occupied, logical)))
            self.assertEqual(direct, trace_g + trace_gs)
            surviving_terms.append((trace_g != CZ, trace_gs != CZ))
        self.assertIn((True, False), surviving_terms)
        self.assertIn((False, True), surviving_terms)

    def test_pauli_midpoint_branch_allows_irrational_occupied_transfer(self):
        # Controlled-H has a Pauli clean midpoint; a T on a separate logical
        # wire can share its first physical band. Clean-zero action is I x T,
        # so the theorem cannot assume a rational occupied target.
        word = _controlled_h_word(clean=0, target=1)
        split = next(j for j, gate in enumerate(word) if gate[0] == 'TDG') + 1
        first, last = word[:split] + [('T', 2)], word[split:]
        u_a = native_columns(3, first)
        u_b = native_columns(3, _inverse(last))
        clean_z = _pauli(3, 0, 1)
        self.assertEqual(_conjugate(u_a, clean_z), clean_z)
        self.assertEqual(_conjugate(u_b, clean_z), clean_z)
        actual = native_columns(3, first + last)
        self.assertEqual(actual, multiply(adjoint(u_b), u_a))
        occupied = [[actual[2 * row][2 * column] for column in range(4)] for row in range(4)]
        self.assertEqual(occupied, native_columns(2, [('T', 1)]))
        row, _ = _transfer_row(occupied, _pauli(2, 2, 0))
        self.assertTrue(set(row.values()) <= set(_alphabet(4)))
        self.assertIn(Q2(0, Fraction(1, 2)), row.values())

    def test_single_alphabet_and_gap_and_r75_return_threshold(self):
        alphabet = _alphabet(2)
        tau_squared = Q2(Fraction(1, 16), Fraction(-1, 32))
        for arity in range(3, 65):
            bias = ONE - Fraction(1, 2 ** (arity - 1))
            nearest = Q2(0, Fraction(1, 2)) if arity == 3 else ONE
            distance = min(abs(bias - candidate) for candidate in alphabet)
            self.assertEqual(distance, abs(bias - nearest))
            floor = (Q2(Fraction(3, 8), Fraction(-1, 4)) if arity == 3
                     else Q2(Fraction(1, 2 ** arity)))
            self.assertEqual(distance / 2, floor)
            self.assertLess(floor * floor, tau_squared)
            # Every omitted alphabet point has magnitude <1/2. This exact
            # inequality certifies the infinite alphabet tail uniformly.
            self.assertGreater(bias - distance, Fraction(1, 2))
        self.assertNotIn(Q2(Fraction(3, 4)), alphabet)
        self.assertGreater(Q2(Fraction(3, 8), Fraction(-1, 4)), Fraction(1, 64))
        # For every a>=4, b>=7/8 is above the midpoint of 1 and 1/sqrt2.
        self.assertGreater(Q2(Fraction(7, 8)), (ONE + Q2(0, Fraction(1, 2))) / 2)

    def test_cubic_bent_full_scalar_transfer_alphabet_is_not_a_realization(self):
        # f(x,y)=x.y+x0*x1*x2 on six addresses. No native circuit is claimed.
        # This tests ALL possible nonzero PTM-entry magnitudes of its XOR
        # query via the reduced Walsh formula; the remaining entries vanish
        # by their unmatched address-shift or target-Z labels.
        size = 64
        values = [(((z & 7) & (z >> 3)).bit_count() & 1) ^ int((z & 7) == 7)
                  for z in range(size)]
        anf = list(values)
        for bit in range(6):
            for mask in range(size):
                if mask & (1 << bit):
                    anf[mask] ^= anf[mask ^ (1 << bit)]
        self.assertEqual({mask for mask, coefficient in enumerate(anf) if coefficient},
                         {7, 9, 18, 36})
        spectrum = _integer_walsh(values)
        self.assertEqual(set(spectrum), {-8, 8})
        magnitudes = set()
        for shift in range(size):
            derivative = [values[x] ^ values[x ^ shift] for x in range(size)]
            derivative_spectrum = _integer_walsh(derivative)
            if shift == 0:
                self.assertEqual(derivative_spectrum, [size] + [0] * (size - 1))
            elif shift & 7:
                self.assertEqual(sum(v != 0 for v in derivative_spectrum), 4)
                self.assertLessEqual(set(derivative_spectrum), {0, -32, 32})
            else:
                self.assertEqual(sum(v != 0 for v in derivative_spectrum), 1)
                self.assertEqual(max(abs(v) for v in derivative_spectrum), size)
            if shift:
                self.assertEqual(derivative_spectrum[0], 0)
            for target_z in (0, 1):
                for frequency in range(size):
                    left = spectrum[frequency] if target_z else (size if frequency == 0 else 0)
                    right = ((-1 if (shift & frequency).bit_count() % 2 else 1)
                             * spectrum[frequency] if target_z else derivative_spectrum[frequency])
                    for output_shift_difference in (0, 1):
                        numerator = left + (-1 if output_shift_difference else 1) * right
                        magnitudes.add(abs(Fraction(numerator, 2 * size)))
        self.assertEqual(magnitudes, {Fraction(0), Fraction(1), Fraction(1, 2),
                                      Fraction(1, 4), Fraction(1, 8)})

        # Independently validate complete representative PTM rows against
        # the literal 128-input XOR permutation, retaining Hermitian Pauli
        # phases. This checks the reduced formula, not merely its spectrum.
        query = [basis ^ (64 if values[basis & 63] else 0) for basis in range(128)]
        for shift, input_target_shift, input_frequency, target_z in (
                (0, 0, 0, 0), (1, 1, 5, 0), (8, 0, 17, 0),
                (7, 1, 31, 1), (63, 0, 42, 1)):
            input_x = shift | (input_target_shift << 6)
            input_z = input_frequency | (target_z << 6)
            derivative_spectrum = _integer_walsh([values[x] ^ values[x ^ shift]
                                                   for x in range(size)])
            for frequency, output_target_shift in product(range(size), range(2)):
                output_x = shift | (output_target_shift << 6)
                output_z = (input_frequency ^ frequency) | (target_z << 6)
                real_sum = 0
                for basis in range(128):
                    before_pauli = query[basis]
                    image = query[before_pauli ^ input_x]
                    if image ^ output_x == basis:
                        parity = ((input_z & before_pauli).bit_count()
                                  + (output_z & image).bit_count()) % 2
                        real_sum += -1 if parity else 1
                phase = ((input_x & input_z).bit_count()
                         + (output_x & output_z).bit_count()) % 4
                if phase % 2:
                    self.assertEqual(real_sum, 0)
                actual_magnitude = abs(Fraction(real_sum, 128))
                left = spectrum[frequency] if target_z else (size if frequency == 0 else 0)
                right = ((-1 if (shift & frequency).bit_count() % 2 else 1)
                         * spectrum[frequency] if target_z else derivative_spectrum[frequency])
                predicted = abs(Fraction(left + (-1 if input_target_shift ^ output_target_shift
                                                   else 1) * right, 2 * size))
                self.assertEqual(actual_magnitude, predicted)


if __name__ == '__main__':
    unittest.main(verbosity=2)
