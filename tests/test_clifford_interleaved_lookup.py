"""Exact fixtures for Clifford-interleaved lookup and completed cube shears.

The analytic results have separate hypotheses: two cubic blocks with arbitrary
Clifford frames, pulse-only propagation, and completed mask-frame words.
These fixtures do not establish a growing native T-depth lower bound.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CI, CO, CZ, adjoint, eye, multiply
from tests.test_bipartite_shear_mixing import _edge_word, _lookup_sign
from tests.test_common_probe_source import OMEGA, ROOT_HALF, _native
from tests.test_coordinate_weighted_pulses import _even_edge, _norm_squared
from tests.test_dirty_cubic_phase_compiler import _dirty_cubic_phase
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_weak_pauli_orbit_compression import _binary_rank


def _edge_specs(p, parity, shift=0):
    return [(subset, bit, 3 * subset + 5 * bit + 1 + 2 * shift, bool(bit % 2))
            for subset in range(1 << p) if subset.bit_count() % 2 == parity
            for bit in range(p) if not (subset >> bit) & 1]


def _tail_layer(p, parity, edges):
    """Paid quarter-angle pulse with independent tail weights and edge types."""
    size = 1 << p
    linear, quadratic, cubic, even = [0] * (p + size), [], [], []
    for subset, bit, weight, signed in edges:
        assert subset.bit_count() % 2 == parity and not (subset >> bit) & 1
        a, b = p + subset, p + (subset | (1 << bit))
        half, odd = divmod(weight % 8, 2)
        even += _even_edge(p, subset, bit, signed, half)
        if odd and signed:
            linear[b] -= 1
            quadratic += [(bit, b, 1), (a, b, 1)]
            cubic.append((bit, a, b))
        elif odd:
            linear[b] += 1
            quadratic.append((a, b, -1))
    phase, width = _dirty_cubic_phase(linear, quadratic, cubic)
    wall = [('H', p + subset) for subset in range(size)
            if subset.bit_count() % 2 == parity]
    return even + wall + phase + list(reversed(wall)), width


def _separate_edges(p, edges):
    return [gate for subset, bit, weight, signed in edges
            for gate in _edge_word(p, subset, bit, signed) * (weight % 8)]


def _completed_half(p, parity, relative_phase=False):
    """An exact completed CCX bank, optionally retaining its i^(x*b) phase."""
    size = 1 << p
    quadratic, cubic = [], []
    for subset in range(size):
        if subset.bit_count() % 2 != parity:
            continue
        for bit in range(p):
            if not (subset >> bit) & 1:
                a, b = p + subset, p + (subset | (1 << bit))
                cubic.append((bit, a, b))
                if relative_phase:
                    quadratic.append((bit, b, 1))
    phase, width = _dirty_cubic_phase([0] * (p + size), quadratic, cubic)
    wall = [('H', p + subset) for subset in range(size)
            if subset.bit_count() % 2 == parity]
    return wall + phase + list(reversed(wall)), width


def _unit_layer(p, parity, signed):
    edges = [(subset, bit, 1, signed)
             for subset in range(1 << p) if subset.bit_count() % 2 == parity
             for bit in range(p) if not (subset >> bit) & 1]
    return _tail_layer(p, parity, edges)


def _physical_coordinate_frame(p):
    return [gate for bit in range(p)
            for subset in range(1 << p) if not (subset >> bit) & 1
            for gate in _even_edge(p, subset, bit, False, 1)]


def _anf(values):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for subset in range(len(values)):
            if (subset >> bit) & 1:
                coefficients[subset] ^= coefficients[subset ^ (1 << bit)]
    return coefficients


def _coefficient_sign(address, bank):
    return sum((bank >> subset) & 1 for subset in range(1 << 3)
               if subset & address == subset) & 1


def _identity(p):
    return tuple(1 << row for row in range(1 << p))


def _compose(left, right):
    """Binary matrix product, storing each row as an integer mask."""
    result = []
    for row in left:
        value = 0
        while row:
            bit = row & -row
            value ^= right[bit.bit_length() - 1]
            row ^= bit
        result.append(value)
    return tuple(result)


def _inverse_rows(rows):
    size = len(rows)
    work = [row | (1 << (size + i)) for i, row in enumerate(rows)]
    for bit in range(size):
        pivot = next(i for i in range(bit, size) if (work[i] >> bit) & 1)
        work[bit], work[pivot] = work[pivot], work[bit]
        for i in range(size):
            if i != bit and (work[i] >> bit) & 1:
                work[i] ^= work[bit]
    assert [row & ((1 << size) - 1) for row in work] == [1 << i for i in range(size)]
    return tuple(row >> size for row in work)


def _shear(p, parity, mask):
    rows = list(_identity(p))
    for subset in range(1 << p):
        if subset.bit_count() % 2 == parity:
            for bit in range(p):
                if (mask >> bit) & 1 and not (subset >> bit) & 1:
                    rows[subset] ^= 1 << (subset | (1 << bit))
    return tuple(rows)


def _wedge(p, u, v):
    return tuple(sum(((((u >> i) & 1) * ((v >> j) & 1))
                      ^ (((v >> i) & 1) * ((u >> j) & 1))) << j
                     for j in range(p)) for i in range(p))


def _xor_rows(*matrices):
    return tuple(_xor(entries) for entries in zip(*matrices))


def _xor(values):
    result = 0
    for value in values:
        result ^= value
    return result


def _central(p, matrix):
    result = _identity(p)
    for i in range(p):
        for j in range(i + 1, p):
            if (matrix[i] >> j) & 1:
                pair = (1 << i) | (1 << j)
                factor = tuple((1 << subset) ^ (0 if subset & pair else 1 << (subset | pair))
                               for subset in range(1 << p))
                result = _compose(result, factor)
    return result


def _pfaffian(matrix, subset):
    if not subset:
        return 1
    if subset.bit_count() % 2:
        return 0
    first = subset & -subset
    i, rest = first.bit_length() - 1, subset ^ first
    return _xor(((matrix[i] >> j) & 1) * _pfaffian(matrix, rest ^ (1 << j))
                for j in range(len(matrix)) if (rest >> j) & 1)


def _frame_normal_form(p, factors):
    a = b = 0
    central = (0,) * p
    actual = _identity(p)
    for parity, mask in factors:
        actual = _compose(actual, _shear(p, parity, mask))
        if parity == 0:
            central = _xor_rows(central, _wedge(p, mask, b))
            a ^= mask
        else:
            b ^= mask
    return actual, central, a, b


def _clique(p, address):
    return tuple((address ^ (1 << i)) if (address >> i) & 1 else 0 for i in range(p))


class CliffordInterleavedLookupTests(unittest.TestCase):
    def test_literal_clifford_columns_exhibit_both_rounding_thresholds(self):
        delta_squared = 2 * CO - 2 * ROOT_HALF
        hadamard_column = _native([('H', 0)], {0: CO})
        difference = {0: hadamard_column[0] - CO, 1: hadamard_column[1]}
        self.assertEqual(_norm_squared(difference), delta_squared)

        x = [('H', 0), ('S', 0), ('S', 0), ('H', 0)]
        # T X T^dagger = omega X S^dagger is Hermitian and Clifford,
        # but its literal nonzero entries are eighth rather than fourth roots.
        clifford = [('TDG', 0)] + x + [('T', 0)]
        for basis, phase in ((0, OMEGA), (1, OMEGA.conj())):
            actual = _native(clifford, {basis: CO})
            self.assertEqual(actual, {basis ^ 1: phase})
            self.assertEqual(_native(clifford, actual), {basis: CO})
            self.assertEqual(_native(_inverse(clifford), {basis: CO}), actual)
            self.assertEqual(_norm_squared({basis ^ 1: phase - CO}), delta_squared)
            # The literal diagonal conjugation cancels the eighth-root phase.
            rounded = [('T', 0)] + clifford + [('TDG', 0)]
            self.assertEqual(_native(rounded, {basis: CO}), {basis ^ 1: CO})

    def test_affine_final_wall_preserves_degree_and_partial_hadamard_support(self):
        p, helper = 3, 11
        wall = [('CX', 10, 5), ('CX', 4, 10), ('CX', 0, 10),
                ('H', 6), ('S', 6), ('S', 6), ('H', 6),
                ('S', 10), ('CX', helper, 10)]

        def transformed_bank(basis, helper_value=None):
            address, bank = basis & 7, (basis >> p) & 255
            work = ((basis >> helper) & 1) if helper_value is None else helper_value
            bank ^= ((bank >> 7) & 1) << 2
            bank ^= (((bank >> 1) & 1) ^ (address & 1) ^ work) << 7
            return bank ^ (1 << 3)

        for address in range(8):
            for bank in (0, 130, 209):
                basis = address | (bank << p) | ((bank & 1) << helper)
                state = _native([('H', helper)] + wall, {basis: CO})
                state = _lookup_sign(state, p)
                actual = _native(_inverse(wall) + [('H', helper)], state)
                phase = (-1) ** _coefficient_sign(address, transformed_bank(basis, 0))
                output = basis ^ ((address == 7) << helper)
                self.assertEqual(actual, {output: phase * CO})
        coefficients = _anf([_coefficient_sign(basis & 7, transformed_bank(basis))
                             for basis in range(1 << 12)])
        self.assertEqual(max(subset.bit_count() for subset, value in enumerate(coefficients) if value), 4)
        flip_coefficients = _anf([int(address == 7) for address in range(8)])
        self.assertEqual([subset for subset, value in enumerate(flip_coefficients) if value], [7])

    def test_tail_weighted_mixed_edge_compiler_pays_phases_return_and_ledgers(self):
        for p in (3, 5):
            size, edges_count = 1 << p, p * (1 << p) // 4
            for parity in (0, 1):
                edges = _edge_specs(p, parity)
                word, width = _tail_layer(p, parity, edges)
                self.assertLessEqual(width - p - size, 5 * edges_count)
                for native in (word, _inverse(word)):
                    counts = Counter(gate[0] for gate in native)
                    self.assertLessEqual(counts['T'] + counts['TDG'], size + 32 * edges_count)
                    self.assertLessEqual(len(native), 6 * size + 71 * edges_count)
                    self.assertLessEqual(len(_bands(native)), 13)
                    for band in _bands(native):
                        self.assertEqual(len({gate[1] for gate in band}), len(band))
                if p == 3:
                    pool = ((1 << (width - p - size)) - 1) << (p + size)
                    reference = 1 << width
                    state = {5 | (37 << p): ROOT_HALF,
                             2 | (213 << p) | pool | reference: CI * ROOT_HALF}
                    actual = _native(word, state)
                    self.assertEqual(actual, _native(_separate_edges(p, edges), state))
                    self.assertEqual(_native(_inverse(word), actual), state)

    def test_tail_filtration_has_complete_error_two_before_enough_alternations(self):
        p, top = 3, 7 | (1 << 10)
        witness = {top: ROOT_HALF, top | (1 << 13): CI * ROOT_HALF}
        for parities, expected_threshold in (((0, 1), 1), ((1, 0, 0, 0, 0), 2)):
            backward, threshold = [], p
            for index, parity in enumerate(parities):
                backward += _separate_edges(p, _edge_specs(p, parity, index))
                if parity == (threshold - 1) % 2:
                    threshold -= 1
            self.assertEqual(threshold, expected_threshold)
            pulled = _native(backward, witness)
            self.assertTrue(all(not (basis >> p) & 1 for basis in pulled))
            echo = backward + [('S', p), ('S', p)] + _inverse(backward)
            self.assertEqual(_native(echo, witness), witness)
            self.assertEqual(_lookup_sign(witness, p),
                             {basis: -value for basis, value in witness.items()})
            self.assertEqual(_native(backward + _inverse(backward), witness), witness)

    def test_a_zero_t_clifford_chain_crosses_the_filtration_but_is_not_lookup(self):
        p, backward = 3, []
        for rank in range(p, 0, -1):
            backward += _even_edge(p, (1 << (rank - 1)) - 1, rank - 1, False, 1)
        self.assertFalse(_bands(backward))
        top = 7 | (1 << 10)
        propagated_bank = sum(1 << subset for subset in (0, 1, 3, 7))
        self.assertEqual(_native(backward, {top: CO}), {7 | (propagated_bank << p): -CI})
        echo = backward + [('S', p), ('S', p)] + _inverse(backward)
        self.assertEqual(_native(echo, {top: CO}), {top: -CO})
        # This address-independent Clifford passes the propagation witness
        # only at selected addresses; at address zero it still has gap two.
        zero_address = 1 << 10
        self.assertEqual(_native(echo, {zero_address: CO}), {zero_address: -CO})
        self.assertEqual(_lookup_sign({zero_address: CO}, p), {zero_address: CO})

    def test_completed_frame_normal_form_retains_pfaffians_and_actual_inverse_row(self):
        for p in (2, 3, 5, 6):
            factors = []
            for i in range(0, p - 1, 2):
                factors += [(0, 1 << i), (1, 1 << (i + 1))] * 2
            factors += [(0, (1 << p) - 1), (1, (1 << p) // 3)]
            actual, matrix, a, b = _frame_normal_form(p, factors)
            central = _central(p, matrix)
            normal = _compose(_compose(central, _shear(p, 0, a)), _shear(p, 1, b))
            self.assertEqual(actual, normal)
            self.assertEqual(central[0], sum(_pfaffian(matrix, subset) << subset
                                             for subset in range(1 << p)))
            inverse = _inverse_rows(actual)
            one_plus_n = _compose(_shear(p, 0, a), _shear(p, 1, a))
            self.assertEqual(inverse[0], _compose(central, one_plus_n)[0])
            self.assertEqual(_compose(actual, inverse), _identity(p))
            v = (1 << p) - 2
            for parity, other in ((0, b), (1, a)):
                conjugate = _compose(_compose(actual, _shear(p, parity, v)), inverse)
                expected = _compose(_central(p, _wedge(p, v, other)), _shear(p, parity, v))
                self.assertEqual(conjugate, expected)
            if p >= 5:
                self.assertEqual((central[0] >> 15) & 1, 1)

    def test_lookup_derivative_fixes_a_and_k_but_leaves_the_right_mask_free(self):
        for p in (3, 5):
            for address in (0, 1, (1 << p) - 1, (1 << p) - 3):
                target = sum(1 << subset for subset in range(1 << p) if subset & address == subset)
                for b in (0, 1, (1 << p) // 3):
                    matrix = _clique(p, address)
                    source = _compose(_compose(_central(p, matrix), _shear(p, 0, address)),
                                      _shear(p, 1, b))
                    self.assertEqual(_inverse_rows(source)[0], target)
                    if p == 3 and address == 7 and b == 1:
                        self.assertNotEqual(source[0], target)
                wrong = list(_clique(p, address))
                wrong[0] ^= 2
                wrong[1] ^= 1
                source = _compose(_central(p, wrong), _shear(p, 0, address))
                # One pair coefficient already distinguishes the wrong K.
                self.assertEqual((_inverse_rows(source)[0] ^ target) & (1 << 3), 1 << 3)

    def test_two_completed_shears_with_free_frames_miss_a_rank_six_target(self):
        p, v1, v2, c1, c2 = 6, 63, 21, 11, 49
        frame1, frame2 = _shear(p, 1, c1), _shear(p, 0, c2)
        first = _compose(_compose(frame2, _shear(p, 1, v2)), frame2)
        second = _compose(_compose(frame1, _shear(p, 0, v1)), frame1)
        source = _compose(first, second)
        matrix = _xor_rows(_wedge(p, c2, v2), _wedge(p, v1, c1), _wedge(p, v1, v2))
        normal = _compose(_compose(_central(p, matrix), _shear(p, 0, v1)), _shear(p, 1, v2))
        self.assertEqual(source, normal)
        self.assertLessEqual(_binary_rank(matrix), 4)
        self.assertEqual(_binary_rank(_clique(p, v1)), 6)
        pair = next((1 << i) | (1 << j) for i in range(p) for j in range(i + 1, p)
                    if not (matrix[i] >> j) & 1)
        # At the all-one address, one occupied pair bit must have opposite
        # signs under the actual derivative and the complete lookup target.
        self.assertEqual((_inverse_rows(source)[0] >> pair) & 1, 0)
        self.assertEqual((_central(p, _clique(p, v1))[0] >> pair) & 1, 1)

    def test_two_cubic_halves_pay_j_x_with_one_returned_helper_pool(self):
        for p in (2, 3, 5):
            size, edges = 1 << p, p * (1 << p) // 4
            pieces = [_completed_half(p, parity) for parity in (0, 1)]
            word = pieces[0][0] + pieces[1][0]
            width = max(piece[1] for piece in pieces)
            self.assertEqual(width - p - size, edges)
            for native in (word, _inverse(word)):
                counts = Counter(gate[0] for gate in native)
                self.assertEqual(counts['T'] + counts['TDG'], 4 * p * size)
                self.assertLessEqual(len(native), 8 * p * size + 2 * size)
                self.assertEqual(len(_bands(native)), 16)
                for band in _bands(native):
                    self.assertEqual(len({gate[1] for gate in band}), len(band))
            if p <= 3:
                pool = ((1 << edges) - 1) << (p + size)
                for address in (0, 1, (1 << p) - 1):
                    for bank in (0, 1 << (size - 1), (1 << size) - 1):
                        output = bank
                        for subset in range(size):
                            for bit in range(p):
                                if (address >> bit) & 1 and not (subset >> bit) & 1:
                                    output ^= ((bank >> (subset | (1 << bit))) & 1) << subset
                        basis = address | (bank << p) | pool
                        expected = address | (output << p) | pool
                        actual = _native(word, {basis: CO})
                        self.assertEqual(actual, {expected: CO})
                        self.assertEqual(_native(_inverse(word), actual), {basis: CO})

    def test_physical_frame_trial_and_its_calibrated_fusion_keep_their_gaps(self):
        # The analytic zero-address relative block is a tensor product of
        # two SU(2) matrices with trace one. Its squared defect has spectrum
        # {0,0,3,3}; this checks the sqrt(3) block gap without adjoining sqrt(3).
        b_plus = [[ROOT_HALF, (CO + CI) / 2], [(-CO + CI) / 2, ROOT_HALF]]
        b_minus = [[ROOT_HALF, (-CO + CI) / 2], [(CO + CI) / 2, ROOT_HALF]]
        left, right = multiply(b_plus, adjoint(b_minus)), multiply(adjoint(b_plus), b_minus)
        for factor in (left, right):
            self.assertEqual(factor[0][0] + factor[1][1], CO)
            self.assertEqual(multiply(factor, adjoint(factor)), eye(2))
        relative = [[left[i // 2][j // 2] * right[i % 2][j % 2] for j in range(4)]
                    for i in range(4)]
        defect = [[relative[i][j] - (CO if i == j else CZ) for j in range(4)]
                  for i in range(4)]
        squared = multiply(adjoint(defect), defect)
        self.assertEqual(multiply(squared, squared), [[3 * value for value in row] for row in squared])
        self.assertEqual(sum((squared[i][i] for i in range(4)), CZ), 6 * CO)
        for p, pair in ((2, (0, 1)), (3, (0, 2))):
            size, edges = 1 << p, p * (1 << p) // 4
            even, even_width = _unit_layer(p, 0, False)
            even_x, signed_width = _unit_layer(p, 0, True)
            odd, odd_width = _unit_layer(p, 1, False)
            odd_x, _ = _unit_layer(p, 1, True)
            middle, middle_width = _completed_half(p, 1, relative_phase=True)
            frame = _physical_coordinate_frame(p)
            trial = even + _inverse(frame) + odd_x + frame + even_x
            repaired = even + _inverse(frame) + middle + frame + even_x
            width = max(even_width, signed_width, odd_width, middle_width)
            self.assertLessEqual(width - p - size, 5 * edges)
            for native in (repaired, _inverse(repaired)):
                counts = Counter(gate[0] for gate in native)
                self.assertLessEqual(counts['T'] + counts['TDG'], 2 * size + 84 * edges)
                self.assertLessEqual(len(_bands(native)), 38)
                self.assertLessEqual(len(native), 13 * size + 194 * edges)
            pool = ((1 << (width - p - size + 1)) - 1) << (p + size)
            i, j = pair
            top = 1 << ((1 << i) | (1 << j))
            zero_address = (top << p) | pool
            witness = {zero_address: CO}
            self.assertEqual(_native(repaired, witness), witness)
            actual = _native(_inverse(trial) + [('S', p), ('S', p)] + trial, witness)
            difference = {basis: actual.get(basis, CZ) - witness.get(basis, CZ)
                          for basis in actual.keys() | witness.keys()}
            self.assertEqual(_norm_squared(difference), 2 * CO)
            marked = (1 << i) | (1 << j) | ((top | (1 << (1 << i))) << p) | pool
            state = {marked: ROOT_HALF, marked ^ pool: CI * ROOT_HALF}
            self.assertEqual(_native(middle, state), _native(odd + odd_x, state))
            echo = _inverse(repaired) + [('S', p), ('S', p)] + repaired
            self.assertEqual(_native(echo, state), {basis: -value for basis, value in state.items()})
            self.assertEqual(_lookup_sign(state, p), state)
            self.assertEqual(_native(repaired + _inverse(repaired), state), state)

    def test_nonzero_address_commutators_and_mixed_frame_normalization(self):
        # Every nonzero address translation leaves an independent cubic
        # coefficient containing a_123; even a bank translation cannot cancel it.
        for shift in range(1, 8):
            for bank_shift in (0, 129):
                values = [_coefficient_sign(basis & 7, basis >> 3)
                          ^ _coefficient_sign((basis & 7) ^ shift, (basis >> 3) ^ bank_shift)
                          for basis in range(1 << 11)]
                coefficients = _anf(values)
                self.assertEqual(max(mask.bit_count() for mask, value in enumerate(coefficients) if value), 3)
                for bit in range(3):
                    if (shift >> bit) & 1:
                        monomial = (1 << 10) | (7 ^ (1 << bit))
                        self.assertEqual(coefficients[monomial], 1)
        # A Y factor and a helper translation retain the same complete sign
        # commutator. Their native scalar phases cancel against the actual inverse.
        x_gate = lambda bit: [('H', bit), ('S', bit), ('S', bit), ('H', bit)]
        pauli = [('S', 0)] * 3 + x_gate(0) + [('S', 0)]
        pauli += x_gate(2) + x_gate(3) + x_gate(10) + x_gate(11)
        state = {5 | (129 << 3): ROOT_HALF,
                 2 | (57 << 3) | (1 << 11) | (1 << 12): CI * ROOT_HALF}
        actual = _native(_inverse(pauli), state)
        actual = _lookup_sign(actual, 3)
        actual = _lookup_sign(_native(pauli, actual), 3)
        expected = {}
        for basis, value in state.items():
            address, bank = basis & 7, (basis >> 3) & 255
            exponent = _coefficient_sign(address, bank) ^ _coefficient_sign(address ^ 5, bank ^ 129)
            expected[basis] = (-1) ** exponent * value
        self.assertEqual(actual, expected)

        # Three wires x,bank,helper: C mixes address and helper through a
        # signed affine wall B, yet C*B is an address-preserving Clifford.
        b_wall = [('CX', 2, 0)] + x_gate(0)
        prime = [('H', 1), ('CX', 0, 1), ('H', 2), ('CX', 0, 2), ('H', 2)]
        clifford = _inverse(b_wall) + prime
        for basis in range(8):
            self.assertEqual(_native(b_wall + clifford, {basis: CO}),
                             _native(prime, {basis: CO}))
            pulled_z = clifford + [('S', 0), ('S', 0)] + _inverse(clifford)
            sign = (-1) ** (1 + (basis & 1) + ((basis >> 2) & 1))
            self.assertEqual(_native(pulled_z, {basis: CO}), {basis: sign * CO})
            normalized = b_wall + pulled_z + _inverse(b_wall)
            self.assertEqual(_native(normalized, {basis: CO}), {basis: (-1) ** (basis & 1) * CO})
            bank_image = prime + [('S', 1), ('S', 1)] + _inverse(prime)
            self.assertEqual(_native(bank_image, {basis: CO}), {basis ^ 2: (-1) ** (basis & 1) * CO})


if __name__ == '__main__':
    unittest.main()
