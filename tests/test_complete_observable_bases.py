"""Exact fixtures for complete observable sources and their basis interfaces.

The finite checks retain arbitrary input columns, phases and actual inverses.
They do not supply an asymptotic compiler or an unrestricted depth lower bound.
"""
from collections import Counter
from functools import lru_cache
from itertools import product
from random import Random
import unittest

from tests.exact_matrices import CI, CO, CZ, adjoint, eye, multiply, native_columns
from tests.test_common_probe_source import _native
from tests.test_joint_lookup_pauli_source import _bands, _inverse


def _controlled_h(control, target):
    # R = S H T H S-dagger, so R Z R-dagger is literally H.
    r = [('S', target)] * 3 + [('H', target), ('T', target),
                              ('H', target), ('S', target)]
    cz = [('H', target), ('CX', control, target), ('H', target)]
    return _inverse(r) + cz + r


def _controlled_s(control, target):
    return [('T', control), ('T', target), ('CX', control, target),
            ('TDG', target), ('CX', control, target)]


def _three_clifford_fibers():
    return (_controlled_h(0, 3) + _controlled_s(1, 3)
            + _controlled_h(2, 3))


def _anf(values):
    values = list(values)
    for bit in range((len(values) - 1).bit_length()):
        for subset in range(len(values)):
            if subset & (1 << bit):
                values[subset] ^= values[subset ^ (1 << bit)]
    return values


def _hadamard_numerator(state, wire):
    """Unnormalized H; the caller retains its common power-of-two scale."""
    result = {}
    for basis, amplitude in state.items():
        zero = basis & ~(1 << wire)
        result[zero] = result.get(zero, 0) + amplitude
        one = zero | (1 << wire)
        result[one] = result.get(one, 0) + (
            -amplitude if basis & (1 << wire) else amplitude)
    return {basis: amplitude for basis, amplitude in result.items() if amplitude}


def _permutation(state, transform):
    result = {}
    for basis, amplitude in state.items():
        target = transform(basis)
        result[target] = result.get(target, 0) + amplitude
    return {basis: amplitude for basis, amplitude in result.items() if amplitude}


def _registers(basis, width):
    mask = (1 << width) - 1
    return basis & mask, (basis >> width) & mask, basis >> (2 * width)


def _pack(c, a, b, width):
    return c | (a << width) | (b << (2 * width))


def _choi_columns(width, table, basis, reduced=False):
    """Apply the two literal register words, with denominator 2**width."""
    state = {basis: 1}
    for wire in range(width, 2 * width):
        state = _hadamard_numerator(state, wire)
    if reduced:
        # H_A; CX(B,C); Q_k(A;C); H_B; (C,A,B) -> (B,C,A).
        state = _permutation(state, lambda n: _pack(
            _registers(n, width)[0] ^ _registers(n, width)[2]
            ^ table[_registers(n, width)[1]],
            _registers(n, width)[1], _registers(n, width)[2], width))
        for wire in range(2 * width, 3 * width):
            state = _hadamard_numerator(state, wire)
        return _permutation(state, lambda n: _pack(
            _registers(n, width)[2], _registers(n, width)[0],
            _registers(n, width)[1], width))
    # Q_k(A;B), SWAP(A,B), then the Bell CX(C,A), H_C.
    state = _permutation(state, lambda n: _pack(
        _registers(n, width)[0],
        _registers(n, width)[2] ^ table[_registers(n, width)[1]]
        ^ _registers(n, width)[0], _registers(n, width)[1], width))
    for wire in range(width):
        state = _hadamard_numerator(state, wire)
    # The coherent Z_k correction and the ordinary CZ(C,A).
    result = {}
    for output, amplitude in state.items():
        t, u, y = _registers(output, width)
        parity = ((t & table[y]).bit_count() + (t & u).bit_count()) & 1
        result[output] = -amplitude if parity else amplitude
    return result


@lru_cache(maxsize=None)
def _ring_product(left, right, p):
    """Square-zero multiplication, retaining the arbitrary unit coefficients."""
    answer = 0
    for a in range(1 << p):
        if (left >> a) & 1:
            for b in range(1 << p):
                if not (a & b) and (right >> b) & 1:
                    answer ^= 1 << (a | b)
    return answer


@lru_cache(maxsize=None)
def _hrt_unit(key, p):
    result = 1
    for subset in range(1, 1 << p):
        if (key >> subset) & 1:
            result = _ring_product(result, 1 | (1 << subset), p)
    return result


def _quadratic_inverse(key, p):
    result = key
    for a in range(1, 1 << p):
        for b in range(a + 1, 1 << p):
            if not (a & b) and (key >> a) & 1 and (key >> b) & 1:
                result ^= 1 << (a | b)
    return result


def _pairs(p):
    return tuple((1 << i) | (1 << j)
                 for i in range(p) for j in range(i + 1, p))


@lru_cache(maxsize=None)
def _pair_unit(mask, p):
    result = 1
    for bit, pair in enumerate(_pairs(p)):
        if (mask >> bit) & 1:
            result = _ring_product(result, 1 | (1 << pair), p)
    return result


def _affine_pair_action(d, c, mask, p):
    g = _pair_unit(mask, p)
    return (_ring_product(g, 1 ^ d, p) ^ 1,
            _ring_product(g, c, p))


def _partial_pair_basis(d, c, p):
    s = sum(((d >> pair) & 1) << bit for bit, pair in enumerate(_pairs(p)))
    residual_d, residual_c = _affine_pair_action(d, c, s, p)
    return s, residual_d, residual_c


def _address_unit(address, p):
    # Independent coefficient definition, not the HRT product routine.
    return sum(1 << subset for subset in range(1 << p)
               if subset & address == subset)


def _address_key(address, p):
    return sum(((address >> bit) & 1) << (1 << bit) for bit in range(p))


def _upper_zeta(bank, p):
    for bit in range(p):
        for subset in range(1 << p):
            if not subset & (1 << bit):
                bank ^= ((bank >> (subset | (1 << bit))) & 1) << subset
    return bank


class CompleteObservableBasesTests(unittest.TestCase):
    def test_native_three_block_columns_actual_inverse_and_ledger(self):
        word = _three_clifford_fibers()
        actual = native_columns(4, word)
        self.assertEqual(native_columns(4, _inverse(word)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(16))
        for address in range(8):
            fiber = ([('H', 3)] if address & 1 else [])
            fiber += [('S', 3)] if address & 2 else []
            fiber += [('H', 3)] if address & 4 else []
            for bank in (0, 1):
                basis = address | (bank << 3)
                expected = _native(fiber, {basis: CO})
                self.assertEqual({row: actual[row][basis] for row in range(16)
                                  if actual[row][basis] != CZ}, expected)
        for native in (word, _inverse(word)):
            counts = Counter(gate[0] for gate in native)
            self.assertEqual(counts['T'] + counts['TDG'], 7)
            self.assertLessEqual(len(native), 39)
            self.assertLessEqual(len(_bands(native)), 6)
            self.assertTrue(all(len({gate[1] for gate in band}) == len(band)
                                for band in _bands(native)))

    def test_three_block_pauli_labels_reach_degree_three_with_literal_sign(self):
        word = _three_clifford_fibers()
        echo = _inverse(word) + [('S', 3), ('S', 3)] + word
        x_labels, z_labels, signs = [], [], []
        for address in range(8):
            x0, x1, x2 = [(address >> bit) & 1 for bit in range(3)]
            sign = x0 * x1 * x2
            a = x0 ^ x2 ^ sign
            b = 1 ^ x0 ^ (x0 * x1) ^ x2 ^ sign
            x_labels.append(a)
            z_labels.append(b)
            signs.append(sign)
            for bank in (0, 1):
                basis = address | (bank << 3)
                phase = (CI if a * b else CO) * (-CO if sign ^ (b * bank) else CO)
                self.assertEqual(_native(echo, {basis: CO}),
                                 {address | ((bank ^ a) << 3): phase})
        for values in (x_labels, z_labels, signs):
            self.assertEqual(_anf(values)[7], 1)

    def test_zero_block_square_root_two_lookup_boundary(self):
        # X on the empty coefficient anticommutes with every lookup fiber.
        # Two applications of X-Q give precisely 2 I on each column.
        for p in range(1, 5):
            m = 1 << p
            banks = range(1 << m) if p <= 2 else (0, 1, 17, (1 << m) - 1)
            for address, bank in product(range(m), banks):
                label = _address_unit(address, p)
                q = -1 if (bank & label).bit_count() & 1 else 1
                q_flipped = -1 if ((bank ^ 1) & label).bit_count() & 1 else 1
                self.assertEqual(q_flipped, -q)
                # The off-diagonal entries of (X-Q)^2 cancel literally.
                self.assertEqual(-(q + q_flipped), 0)
                self.assertEqual(1 + q * q, 2)

    def test_dirty_choi_normal_form_all_boolean_maps_and_every_column(self):
        for width in (1, 2):
            size = 1 << width
            for table in product(range(size), repeat=size):
                for basis in range(1 << (3 * width)):
                    literal = _choi_columns(width, table, basis)
                    reduced = _choi_columns(width, table, basis, reduced=True)
                    c, a, z = _registers(basis, width)
                    kernel = {_pack(t, c ^ z ^ table[y], y, width):
                              -1 if ((a & y).bit_count() + (t & z).bit_count()) & 1
                              else 1 for t, y in product(range(size), repeat=2)}
                    self.assertEqual(literal, reduced)
                    self.assertEqual(literal, kernel)
                    self.assertEqual(sum(value * value for value in literal.values()), size * size)

    def test_coherent_phase_correction_is_the_full_query_hadamard_dual(self):
        width, size = 2, 4
        for table in product(range(size), repeat=size):
            for t, y in product(range(size), repeat=2):
                state = {t | (y << width): 1}
                for wire in range(width):
                    state = _hadamard_numerator(state, wire)
                state = _permutation(state, lambda n: n ^ table[n >> width])
                for wire in range(width):
                    state = _hadamard_numerator(state, wire)
                expected = -size if (t & table[y]).bit_count() & 1 else size
                self.assertEqual(state, {t | (y << width): expected})

    def test_partial_pair_gauge_and_inverse_on_every_three_label_unit(self):
        p, m = 3, 8
        pair_bits = sum(1 << pair for pair in _pairs(p))
        for d, c in product(range(0, 1 << m, 2), range(1 << m)):
            s, rd, rc = _partial_pair_basis(d, c, p)
            self.assertEqual(rd & pair_bits, 0)
            self.assertEqual(_affine_pair_action(rd, rc, s, p), (d, c))
        # All pair characters, with each arbitrary coefficient direction and
        # dense fibers, check S A_K = (X^K tensor I) S independently.
        words = (0, 2, 6, 42, 126, 254)
        for d, c, mask in product(words, (0, 1, 85, 170, 255), range(8)):
            s, rd, rc = _partial_pair_basis(d, c, p)
            moved_d, moved_c = _affine_pair_action(d, c, mask, p)
            self.assertEqual(_partial_pair_basis(moved_d, moved_c, p),
                             (s ^ mask, rd, rc))

    def test_affine_pair_double_commutator_keeps_all_dirty_inputs(self):
        p, m = 4, 16
        pairs = _pairs(p)
        first = 1 << pairs.index(3)
        second = 1 << pairs.index(12)
        rng = Random(250)
        columns = [(0, 0), ((1 << m) - 2, (1 << m) - 1)]
        columns += [(1 << bit, 1 << (15 - bit)) for bit in range(1, m)]
        columns += [(rng.getrandbits(m) & ~1, rng.getrandbits(m)) for _ in range(40)]

        def apply(state, token):
            d, c, private = state
            if token == 'R':
                return d, c ^ 1, private
            nd, nc = _affine_pair_action(d, c, first if token == 'I' else second, p)
            return nd, nc, private

        # [R,A_i] is A_i,R,A_i,R chronologically; the outer commutator
        # uses its actual reversed word rather than an involution guess.
        inner = ['I', 'R', 'I', 'R']
        double = ['J'] + list(reversed(inner)) + ['J'] + inner
        for (d, c), private in product(columns, range(8)):
            state = d, c, private
            for token in double:
                state = apply(state, token)
            self.assertEqual(state, (d, c ^ (1 << 15), private))
        # Three disjoint pair increments retain a nonzero witness as well.
        vector = 1
        for pair in (3, 12, 48):
            vector = _ring_product(1 << pair, vector, 6)
        self.assertEqual(vector, 1 << 63)

    def test_hrt_probe_residual_and_quadratic_untwist_on_full_small_keys(self):
        for p in (2, 3):
            m = 1 << p
            values = range(0, 1 << m, 2)
            for key, carrier in product(values, repeat=2):
                inverse = _quadratic_inverse(carrier, p)
                self.assertEqual(_hrt_unit(inverse, p) ^ 1, carrier)
                shift = _hrt_unit(key, p) ^ 1
                affine = carrier ^ shift ^ _ring_product(carrier, shift, p)
                self.assertEqual(affine, _hrt_unit(key ^ inverse, p) ^ 1)
                untwisted = _hrt_unit(key ^ inverse ^ inverse, p)
                self.assertEqual(untwisted, _hrt_unit(key, p))
        # This correction changes the chosen Z probe already on y=0.
        self.assertEqual((_hrt_unit((1 << 1) | (1 << 2), 2) >> 3) & 1, 1)

    def test_scalar_top_phase_hadamard_dual_retains_unital_coefficient(self):
        for p in (1, 2):
            m = 1 << p
            for key, bank in product(range(0, 1 << m, 2), range(1 << m)):
                state = {bank: 1}
                for wire in range(m):
                    state = _hadamard_numerator(state, wire)
                phased = {}
                for dual_bank, value in state.items():
                    # Complement the physical b labels around S_(p+1).
                    reversed_bank = sum(((dual_bank >> subset) & 1)
                                        << (m - 1 - subset) for subset in range(m))
                    full_key = key | (reversed_bank << m)
                    sign = (_hrt_unit(full_key, p + 1) >> ((1 << (p + 1)) - 1)) & 1
                    phased[dual_bank] = -value if sign else value
                for wire in range(m):
                    phased = _hadamard_numerator(phased, wire)
                self.assertEqual(phased, {bank ^ _hrt_unit(key, p): 1 << m})
                self.assertEqual(_hrt_unit(key, p) & 1, 1)

    def test_cubic_bent_dual_walsh_signs_and_dense_operator_kernel(self):
        for p in (1, 2, 3):
            n = (1 << p) - 1
            size = 1 << n
            mask = size - 1
            forward = [(_hrt_unit(a << 1, p) ^ 1) >> 1 for a in range(size)]
            inverse = [_quadratic_inverse(v << 1, p) >> 1 for v in range(size)]
            signs = [-1 if ((basis >> n) & forward[basis & mask]).bit_count() & 1
                     else 1 for basis in range(size * size)]
            # Integer FWHT: dividing this vector by 2**n is precisely the
            # normalized Walsh transform on 2n bits, not an operator entry.
            spectrum = list(signs)
            step = 1
            while step < len(spectrum):
                for start in range(0, len(spectrum), 2 * step):
                    for offset in range(step):
                        left, right = spectrum[start + offset], spectrum[start + offset + step]
                        spectrum[start + offset] = left + right
                        spectrum[start + offset + step] = left - right
                step *= 2
            dual_signs = [-1 if ((basis & mask) & inverse[basis >> n]).bit_count() & 1
                          else 1 for basis in range(size * size)]
            self.assertEqual(spectrum, [size * sign for sign in dual_signs])
            if p <= 2:
                # Literal H D_(g*) H on every column has denominator 2**(2n).
                # Its nonzero off-diagonal entries distinguish the dense
                # Fourier kernel from the required diagonal D_g operation.
                for column in range(size * size):
                    state = {column: 1}
                    for wire in range(2 * n):
                        state = _hadamard_numerator(state, wire)
                    state = {row: value * dual_signs[row] for row, value in state.items()}
                    for wire in range(2 * n):
                        state = _hadamard_numerator(state, wire)
                    expected = {row: size * signs[row ^ column] for row in range(size * size)}
                    self.assertEqual(state, expected)
                    self.assertNotEqual(state[column ^ 1], 0)

    def test_four_corner_source_returns_every_key_and_bank(self):
        for p in (1, 2):
            m = 1 << p
            values = range(1 << m)
            for address, key, b, c, output in product(
                    range(m), range(0, 1 << m, 2), values, values, values):
                old = key, b, c, output
                u = _hrt_unit(key, p)
                v = _hrt_unit(key ^ _address_key(address, p), p)
                output ^= _ring_product(b, c, p)
                b ^= u
                output ^= _ring_product(b, c, p)
                c ^= v
                output ^= _ring_product(b, c, p)
                b ^= u
                output ^= _ring_product(b, c, p)
                c ^= v
                self.assertEqual((key, b, c), old[:3])
                gain = output ^ old[3]
                self.assertEqual(gain, _address_unit(address, p))
                self.assertEqual(_upper_zeta(gain, p), 1 << address)


if __name__ == '__main__':
    unittest.main()
