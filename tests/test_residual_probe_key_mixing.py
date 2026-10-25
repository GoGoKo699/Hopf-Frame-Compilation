"""Exact checks for residual-probe extraction and later-key-mixing screens.

These fixtures audit complete small operators and algebraic witnesses. They
are not a scalable implementation of the still-missing residual phase.
"""
from collections import Counter
from itertools import combinations, product
from random import Random
import unittest

from tests.test_complete_observable_bases import (
    _address_key, _address_unit, _controlled_h, _hadamard_numerator, _hrt_unit, _pack,
    _permutation, _quadratic_inverse, _registers, _ring_product, _upper_zeta,
)
from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_common_probe_source import _native
from tests.test_joint_lookup_pauli_source import _inverse


def _mobius(values, modulus):
    coefficients = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for subset in range(len(values)):
            if subset & (1 << bit):
                coefficients[subset] = (coefficients[subset]
                                        - coefficients[subset ^ (1 << bit)]) % modulus
    return coefficients


def _degree(values, modulus):
    return max((subset.bit_count() for subset, value in enumerate(_mobius(values, modulus))
                if value), default=-1)


def _integer_sqrt_two_sign(a, b):
    """Sign of a+b*sqrt(2), using integer comparisons only."""
    if not b:
        return (a > 0) - (a < 0)
    if not a:
        return (b > 0) - (b < 0)
    if (a > 0) == (b > 0):
        return 1 if a > 0 else -1
    comparison = a * a - 2 * b * b
    return ((comparison > 0) - (comparison < 0)) * (1 if a > 0 else -1)


def _one_bank_word(address, key, borrowed, output, p):
    k = lambda bank: _quadratic_inverse(bank, p)
    u = _hrt_unit(key, p) ^ 1
    v = _hrt_unit(key ^ _address_key(address, p), p) ^ 1
    output ^= k(borrowed) ^ borrowed             # R,L
    borrowed ^= u
    output ^= 1                                # U, including its unit toggle
    output ^= k(borrowed)                       # R
    borrowed ^= v
    output ^= 1                                # V
    output ^= k(borrowed) ^ borrowed             # R,L
    borrowed ^= u
    output ^= 1                                # actual U inverse, ideal action
    output ^= k(borrowed)                       # R
    borrowed ^= v
    output ^= 1                                # actual V inverse, ideal action
    output ^= 1                                # final literal X on unit output
    return key, borrowed, output


def _apply_aliased_word(state, word):
    # p=1: physical x,a,B,Y_empty,Y_1,extra_0,extra_1. The Xi key
    # is z=(a,Y_empty,B), its carrier y=(Y_1,extra_0,extra_1).
    z_wires, y_wires = (1, 3, 2), (4, 5, 6)

    def extract(basis, wires):
        return sum(((basis >> wire) & 1) << bit for bit, wire in enumerate(wires))

    def xor_into(basis, value, wires):
        return basis ^ sum(((value >> bit) & 1) << wire for bit, wire in enumerate(wires))

    for operation in word:
        name, *args = operation
        if name == 'H':
            state = _hadamard_numerator(state, args[0])
        elif name == 'X':
            state = _permutation(state, lambda basis: basis ^ (1 << args[0]))
        elif name == 'CX':
            state = _permutation(state, lambda basis: basis ^ (
                ((basis >> args[0]) & 1) << args[1]))
        elif name == 'SWAP':
            state = _permutation(state, lambda basis: basis ^ (
                (((basis >> args[0]) ^ (basis >> args[1])) & 1)
                * ((1 << args[0]) | (1 << args[1]))))
        elif name == 'K':
            state = _permutation(state, lambda basis: xor_into(
                basis, _quadratic_inverse(extract(basis, y_wires) << 1, 2) >> 1, z_wires))
        elif name == 'XI':
            updated = {}
            for basis, amplitude in state.items():
                key = extract(basis, z_wires) << 1
                carrier = extract(basis, y_wires) << 1
                sign = (_hrt_unit(key ^ _quadratic_inverse(carrier, 2), 2) >> 3) & 1
                updated[basis] = -amplitude if sign else amplitude
            state = updated
        else:
            raise AssertionError(name)
    return state


class ResidualProbeKeyMixingTests(unittest.TestCase):
    def test_one_bank_polarization_retains_all_old_values(self):
        for p in (1, 2, 3):
            m = 1 << p
            values = range(0, 1 << m, 2)
            inverse = {value: _quadratic_inverse(value, p) for value in values}
            for old, u, v in product(values, repeat=3):
                four_corners = inverse[old] ^ inverse[old ^ u]
                four_corners ^= inverse[old ^ u ^ v] ^ inverse[old ^ v]
                self.assertEqual(four_corners, _ring_product(u, v, p))
                self.assertEqual(four_corners ^ old ^ (old ^ u ^ v),
                                 _ring_product(u, v, p) ^ u ^ v)

    def test_one_bank_complete_word_and_constant_toggle(self):
        for p in (1, 2):
            m = 1 << p
            for address, key, borrowed, output in product(
                    range(m), range(0, 1 << m, 2), range(0, 1 << m, 2), range(1 << m)):
                actual = _one_bank_word(address, key, borrowed, output, p)
                self.assertEqual(actual[:2], (key, borrowed))
                self.assertEqual(actual[2], output ^ _address_unit(address, p))
                self.assertEqual(_upper_zeta(actual[2] ^ output, p), 1 << address)
                self.assertEqual((actual[2] ^ output) & 1, 1)

    def test_literal_phase_query_with_output_carrier_alias_and_actual_reverse(self):
        # This is the physical complement/Hadamard/Q_k/Xi/Q_k query,
        # including a completely arbitrary old constant output and carrier.
        query = [('SWAP', 3, 2), ('H', 3), ('H', 2), ('K',), ('XI',),
                 ('K',), ('H', 2), ('H', 3), ('SWAP', 3, 2)]
        shifted = [('CX', 0, 1)] + query + [('CX', 0, 1)]
        r = [('CX', 2, 4)]
        zeta = [('CX', 4, 3)]
        complete = (zeta + r + r + query + r + shifted + r + r
                    + list(reversed(query)) + r + list(reversed(shifted))
                    + [('X', 3)] + zeta)
        for basis in range(128):
            expected = basis ^ (1 << 3) ^ (((basis >> 1) & 1) << 2)
            self.assertEqual(_apply_aliased_word({basis: 1}, query), {expected: 4})
            self.assertEqual(_apply_aliased_word({basis: 1}, list(reversed(query))),
                             {expected: 4})
            target = basis ^ (1 << (3 + (basis & 1)))
            self.assertEqual(_apply_aliased_word({basis: 1}, complete), {target: 256})
            self.assertEqual(_apply_aliased_word({basis: 1}, list(reversed(complete))),
                             {target: 256})

    def test_weighted_cubic_phase_rounding_population_and_exact_radius(self):
        rows = [((value >> 0) & 1, (value >> 1) & 1, (value >> 2) & 1)
                for value in range(8)]
        monomials = [(a, b, c, a * b, a * c, b * c, a * b * c) for a, b, c in rows]
        largest_population = 0
        nonconstant = 0
        for raw in product(range(8), range(8), range(8),
                           range(4), range(4), range(4), range(2)):
            coefficients = raw[:3] + tuple(2 * value for value in raw[3:6]) + (4 * raw[6],)
            values = [sum(a * b for a, b in zip(coefficients, row)) % 8 for row in monomials]
            counts = Counter(values)
            if len(counts) == 1:
                continue
            nonconstant += 1
            largest_population = max(largest_population, max(counts.values()))
            # 128 |mean(omega**values)|**2 = q+c sqrt(2), exactly.
            a, b = counts[0] - counts[4], counts[1] - counts[3] - counts[5] + counts[7]
            c, d = counts[2] - counts[6], counts[1] + counts[3] - counts[5] - counts[7]
            rational = 2 * a * a + 2 * c * c + b * b + d * d
            irrational = 2 * a * b + 2 * c * d
            # |mean| < 31/32 makes every optimally phased column distance >1/4.
            self.assertLess(_integer_sqrt_two_sign(8 * rational - 961, 8 * irrational), 0)
        self.assertEqual(nonconstant, 65535)
        self.assertEqual(largest_population, 7)

    def test_rounding_equal_and_separated_columns_with_a_moving_prefix(self):
        def ccz(state):
            return {basis: -amplitude if basis & 7 == 7 else amplitude
                    for basis, amplitude in state.items()}

        for basis in range(8):
            # D X_0 D = X_0 CZ_(1,2), including its literal signs.
            shifted = _permutation(ccz({basis: CO}), lambda value: value ^ 1)
            phase = -CO if basis & 6 == 6 else CO
            self.assertEqual(ccz(shifted), {basis ^ 1: phase})
        differences = []
        for basis in (6, 7):
            observed = ccz(_native([('H', 0)], ccz({basis: CO})))
            target = _native([('H', 0)], {basis: CO})
            differences.append({row: observed.get(row, CZ) - target.get(row, CZ)
                                for row in (6, 7)})
        for i, j in product(range(2), repeat=2):
            inner = sum((differences[i][row].conj() * differences[j][row]
                         for row in (6, 7)), CZ)
            self.assertEqual(inner, 2 * CO if i == j else CZ)
        # U=CH H_key moves the key, while its chosen probe image still holds it.
        ch = _controlled_h(0, 1)
        prefix = [('H', 0)] + ch
        matrix = native_columns(2, prefix)
        self.assertEqual(native_columns(2, _inverse(prefix)), adjoint(matrix))
        self.assertTrue(any(matrix[row][0] != CZ for row in (1, 3)))
        echo = _inverse(prefix) + [('S', 1), ('S', 1)] + prefix
        held_echo = _inverse(ch) + [('S', 1), ('S', 1)] + ch
        self.assertEqual(native_columns(2, echo), native_columns(2, held_echo))
        for basis in range(4):
            self.assertTrue(all((row & 1) == (basis & 1)
                                for row in _native(echo, {basis: CO})))

    def test_endpoint_modulo_eight_difference_and_degree_six_carry(self):
        rng = Random(251250)
        p, bank_bits = 7, 3
        masks = [mask for mask in range(1 << p) if mask.bit_count() <= 3]
        address_pairs = tuple(combinations(range(p), 2))
        for _ in range(20):
            labels = [rng.sample(masks, 10) for _ in range(bank_bits)]
            constants = [rng.randrange(8) for _ in range(bank_bits)]
            linear = [[rng.randrange(4) for _ in range(p)] for _ in range(bank_bits)]
            quadratic = [[rng.randrange(2) for _ in address_pairs] for _ in range(bank_bits)]
            pairs = {(j, k): (rng.randrange(4), [rng.randrange(2) for _ in range(p)])
                     for j, k in combinations(range(bank_bits), 2)}
            cubic = rng.randrange(2)
            differences = [[] for _ in range(bank_bits)]
            for key in range(1 << p):
                x = [(key >> bit) & 1 for bit in range(p)]
                a = [sum((key & mask) == mask for mask in terms) & 1 for terms in labels]
                ell = [(constants[j] + 2 * sum(v * bit for v, bit in zip(linear[j], x))
                        + 4 * sum(v * x[i] * x[k] for v, (i, k) in zip(quadratic[j], address_pairs))) % 8
                       for j in range(bank_bits)]
                q = {pair: (constant + 2 * sum(v * bit for v, bit in zip(slope, x))) % 4
                     for pair, (constant, slope) in pairs.items()}

                def phase(value):
                    z = [(value >> bit) & 1 for bit in range(bank_bits)]
                    return (sum(v * bit for v, bit in zip(ell, z))
                            + 2 * sum(v * z[j] * z[k] for (j, k), v in q.items())
                            + 4 * cubic * z[0] * z[1] * z[2]) % 8

                translation = sum(bit << j for j, bit in enumerate(a))
                for j in range(bank_bits):
                    delta = (phase(translation ^ (1 << j)) - phase(translation)
                             - phase(1 << j) + phase(0)) % 8
                    others = [k for k in range(bank_bits) if k != j]
                    qj = lambda k: q[tuple(sorted((j, k)))]
                    formula = (-ell[j] * a[j] + sum(qj(k) * a[k] for k in others)
                               - 2 * sum(qj(k) * a[j] * a[k] for k in others)
                               + 2 * cubic * a[others[0]] * a[others[1]]) % 4
                    self.assertEqual(delta % 2, 0)
                    self.assertEqual(delta // 2, formula)
                    differences[j].append(formula)
            self.assertTrue(all(_degree(values, 4) <= 6 for values in differences))
        # The carry bound is needed: two cubic Boolean monomials give degree
        # six in the integer lift, and a cubic endpoint gives 2*a_1*a_2.
        lifted = [int(((key & 7) == 7) ^ ((key & 56) == 56)) for key in range(64)]
        endpoint = [2 * int((key & 7) == 7) * int((key & 56) == 56) for key in range(64)]
        self.assertEqual((_degree(lifted, 2), _degree(lifted, 4)), (3, 6))
        self.assertEqual(_degree(endpoint, 4), 6)

    def test_residual_probe_contains_degree_seven_carrier_label(self):
        # r=8, y_{8} label: coefficient e_1...e_7 of f_8(z), restricted
        # to its seven singleton key bits. Multiply the formal factors
        # sparsely; no initialized physical bank is used by the theorem.
        values = []
        for key in range(128):
            support = {0}
            for bit in range(7):
                if (key >> bit) & 1:
                    factor = 1 << bit
                    shifted = {subset | factor for subset in support if not subset & factor}
                    support.symmetric_difference_update(shifted)
            values.append(int(127 in support))
        self.assertEqual(values, [int(key == 127) for key in range(128)])
        self.assertEqual(_degree(values, 2), 7)
        self.assertEqual(_mobius(values, 2)[127], 1)

    def test_group_bell_completed_word_and_full_kernel(self):
        r, width = 2, 3
        size = 1 << width
        inverse = [_quadratic_inverse(value << 1, r) >> 1 for value in range(size)]

        def group(a, b):
            return a ^ b ^ (_ring_product(a << 1, b << 1, r) >> 1)

        def query(state, control):
            return _permutation(state, lambda basis: _pack(
                _registers(basis, width)[0], _registers(basis, width)[1],
                _registers(basis, width)[2] ^ inverse[_registers(basis, width)[control]], width))

        def multiply_group(state):
            return _permutation(state, lambda basis: _pack(
                _registers(basis, width)[0], group(*_registers(basis, width)[:2]),
                _registers(basis, width)[2], width))

        for basis in range(size ** 3):
            state = {basis: 1}
            for wire in range(width, 2 * width):
                state = _hadamard_numerator(state, wire)
            reduced = multiply_group(query(state, 0))
            literal = multiply_group(query(state, 1))
            for wire in range(width):
                reduced = _hadamard_numerator(reduced, wire)
                literal = _hadamard_numerator(literal, wire)
            literal = query(literal, 1)
            c, a, z = _registers(basis, width)
            kernel = {_pack(t, u, z ^ inverse[c], width):
                      -1 if ((a & group(c, u)).bit_count() + (t & c).bit_count()) & 1 else 1
                      for t, u in product(range(size), repeat=2)}
            self.assertEqual(literal, reduced)
            self.assertEqual(literal, kernel)

    def test_parallel_choi_carrier_readout_loses_both_input_keys(self):
        # Two disjoint width-one copies, plus an untouched arbitrary output.
        # Each program is the complete one-query Choi normal form, and its
        # reverse is emitted explicitly. The two query tables are different.
        def circuit(state, c, a, b, table, reverse=False):
            def cycle(value, backwards=False):
                cv, av, bv = [(value >> wire) & 1 for wire in (c, a, b)]
                target = (av, bv, cv) if backwards else (bv, cv, av)
                return value ^ sum((old ^ new) << wire for wire, old, new
                                   in zip((c, a, b), (cv, av, bv), target))

            operations = ['HA', 'CX', 'Q', 'HB', 'CYCLE']
            if reverse:
                operations.reverse()
            for operation in operations:
                if operation == 'HA':
                    state = _hadamard_numerator(state, a)
                elif operation == 'HB':
                    state = _hadamard_numerator(state, b)
                elif operation == 'CX':
                    state = _permutation(state, lambda value: value ^ (((value >> b) & 1) << c))
                elif operation == 'Q':
                    state = _permutation(state, lambda value: value ^ (table[(value >> a) & 1] << c))
                else:
                    state = _permutation(state, lambda value: cycle(value, reverse))
            return state

        def readout(state, left, right, xor):
            if xor:
                return _permutation(state, lambda value: value ^ (
                    (((value >> left) & 1) * ((value >> right) & 1)) << 6))
            return {value: -amplitude if (value >> left) & (value >> right) & 1 else amplitude
                    for value, amplitude in state.items()}

        for xor, basis in product((False, True), range(128)):
            actual = circuit({basis: 1}, 0, 1, 2, (0, 1))
            actual = circuit(actual, 3, 4, 5, (1, 0))
            actual = readout(actual, 2, 5, xor)
            actual = circuit(actual, 3, 4, 5, (1, 0), reverse=True)
            actual = circuit(actual, 0, 1, 2, (0, 1), reverse=True)
            expected = _hadamard_numerator({basis: 1}, 1)
            expected = _hadamard_numerator(expected, 4)
            expected = readout(expected, 1, 4, xor)
            expected = _hadamard_numerator(expected, 4)
            expected = _hadamard_numerator(expected, 1)
            # Actual denominator16 versus expected denominator4.
            self.assertEqual(actual, {value: 4 * amplitude for value, amplitude in expected.items()})
            self.assertTrue(all((value & 0b001001) == (basis & 0b001001) for value in actual))


if __name__ == '__main__':
    unittest.main()
