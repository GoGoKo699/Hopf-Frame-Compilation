"""Exact guards for the normalized XOR-carry phase source.

These component identities and p=1,2 fixtures do not establish a scalable
lookup source or an asymptotic exclusion of other interference words.
"""
from fractions import Fraction
import unittest

import numpy as np

from tests.exact_matrices import (
    CO, CQ2, CZ, Q2, adjoint, multiply, native_columns, scale,
)
from tests.test_gram_four_call_readout import _partial_gram_source
from tests.test_half_bilinear_gram_source import (
    _binary_image, _private_toffoli, _walsh_integer,
)
from tests.test_majorana_walsh_query import _inverse
from tests.test_nilpotent_pencil_observable_source import _field_matrix, _power
from tests.test_phase_signal_query_limits import _counts


def _p1_mixer():
    word = []
    for z, t in ((0, 2), (1, 3)):
        word.extend([('CX', z, t), ('CX', t, z), ('CX', z, t)])
    return word + [('CX', 2, 0), ('CX', 3, 0), ('CX', 3, 1)] + [
        ('H', wire) for wire in range(4)]


def _commutator(diagonal, mixer):
    crossed = mixer + diagonal + _inverse(mixer)
    return crossed + diagonal + _inverse(crossed) + _inverse(diagonal)


def _xor_carry_expectation(p, transpose=False):
    """Exact Q(sqrt(2)) receipt at x=1,z=e_p,t=0, with four H walls.

    Values are numerator polynomials in omega, modulo omega^4=-1.
    Both complex embeddings have norm at most 2^(4n), so each integer
    coefficient has magnitude at most 2^(4n)<=2^32. Products and sums
    below use Python integers, not fixed-width squared numerators.
    """
    m, n = 1 << p, p * (1 << p)
    size, dimension = 1 << n, 1 << (2 * n)
    basis = ((3, 2) if p == 1 else
             _field_matrix([[_power(u, j) for u in range(m)] for j in range(m)]))
    images = np.array([_binary_image(basis, z) for z in range(size)],
                      dtype=np.int64)
    labels = np.arange(dimension, dtype=np.int64)
    z, t = labels & (size - 1), labels >> n
    permutation = (images[z] ^ t) | (z << n)
    weights = np.array([value.bit_count() for value in range(size)],
                       dtype=np.int64)
    image = (z << p) & (size - 1) if transpose else z >> p
    exponent = weights[t ^ image] - weights[t]
    values = np.zeros((4, dimension), dtype=np.int64)
    values[0, 1 << p] = 1

    def phase(sign):
        nonlocal values
        new = np.zeros_like(values)
        powers = (sign * exponent) % 8
        for power in range(8):
            mask = powers == power
            for degree in range(4):
                target = (degree + power) % 8
                new[target % 4, mask] = ((-1 if target >= 4 else 1)
                                        * values[degree, mask])
        values = new

    def mixer(inverse=False):
        nonlocal values
        if not inverse:
            new = np.zeros_like(values)
            new[:, permutation] = values
            values = new
        for row in values:
            _walsh_integer(row)
        if inverse:
            values = values[:, permutation]

    mixer()
    phase(1)
    mixer(True)
    phase(1)
    mixer()
    phase(-1)
    mixer(True)
    phase(-1)
    norm_a = norm_b = probe_a = probe_b = 0
    for label, column in enumerate(values.T):
        a, b, c, d = map(int, column)
        rational = a*a + b*b + c*c + d*d
        radical = a*b + b*c + c*d - a*d
        sign = -1 if (label ^ (label >> n)) & 1 else 1
        norm_a += rational
        norm_b += radical
        probe_a += sign * rational
        probe_b += sign * radical
    denominator = 1 << (8 * n)
    return tuple(Fraction(value, denominator) for value in
                 (norm_a, norm_b, probe_a, probe_b))


class TWallXorGramSourceTests(unittest.TestCase):
    def test_literal_phase_includes_every_address_and_dirty_column(self):
        # x,z0,z1,t0,t1,h.  The whole native query is a returned Toffoli.
        query = _private_toffoli(0, 2, 3, 5)
        word = ([('TDG', 3), ('TDG', 4)] + query
                + [('T', 3), ('T', 4)] + query)
        self.assertEqual(_counts(word), (20, 20))
        actual = native_columns(6, word)
        omega = CQ2(Q2(0, Fraction(1, 2)), Q2(0, Fraction(1, 2)))
        powers = [CO]
        for _ in range(7):
            powers.append(powers[-1] * omega)
        for column in range(64):
            toggle = (column & 1) & ((column >> 2) & 1)
            t0 = (column >> 3) & 1
            expected = powers[((t0 ^ toggle) - t0) % 8]
            for row in range(64):
                self.assertEqual(actual[row][column], expected if row == column else CZ)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))

    def test_complete_one_clean_source_has_exact_half_projector_gram(self):
        # Fixed address x=1: both arbitrary carrier banks remain in the test.
        diagonal = [('TDG', 2), ('TDG', 3), ('CX', 1, 2),
                    ('T', 2), ('T', 3), ('CX', 1, 2)]
        transport = _commutator(diagonal, _p1_mixer())
        base = [tuple([gate[0]] + [{0: 0, 1: 2, 2: 4}[wire]
                                  for wire in gate[1:]])
                for gate in _partial_gram_source()]
        actual = native_columns(5, transport + base + _inverse(transport))
        accepted = [row[:16] for row in actual[:16]]
        gram = multiply(adjoint(accepted), accepted)
        v = native_columns(4, transport)
        projector = [[(CO if not ((row ^ (row >> 2)) & 1) else CZ)
                      if row == column else CZ for column in range(16)]
                     for row in range(16)]
        expected = scale(multiply(multiply(adjoint(v), projector), v),
                         CQ2(Fraction(1, 2)))
        self.assertEqual(gram, expected)
        self.assertEqual(multiply(gram, gram), scale(gram, CQ2(Fraction(1, 2))))
        self.assertEqual(gram[2][2], CQ2(Fraction(1, 2)))

    def test_relative_kernel_echo_has_no_unpaid_scalar_or_helper_phase(self):
        # At p=1, P deletes z0 and Q_P is the fixed CNOT z1 -> t1.
        # K_x=(I+B_x)P adds z1 to t1 and x*z1 to t0.
        query = _private_toffoli(0, 2, 3, 5)
        fixed = [('CX', 2, 4)]
        wall = [('T', 3), ('T', 4)]

        def phase(q):
            return _inverse(wall) + q + wall + _inverse(q)

        relative = _inverse(phase(fixed)) + phase(query + fixed)
        conjugated = fixed + phase(query) + fixed
        self.assertEqual(native_columns(6, relative),
                         native_columns(6, conjugated))

    def test_actual_address_hadamard_variant_still_fails_p1(self):
        mixer = [tuple([gate[0]] + [wire + 1 for wire in gate[1:]])
                 for gate in _p1_mixer()] + [('H', 0)]
        query = _private_toffoli(0, 2, 3, 5)
        diagonal = ([('TDG', 3), ('TDG', 4)] + query
                    + [('T', 3), ('T', 4)] + query)
        # x=1,z=e1,t=0; check both arbitrary helper values.
        actual = native_columns(6, _commutator(diagonal, mixer),
                                columns=[5, 5 | 32])
        for column, helper in enumerate((0, 1)):
            expectation = CZ
            for row in range(64):
                if row >> 5 != helper:
                    self.assertEqual(actual[row][column], CZ)
                sign = -1 if ((row >> 1) ^ (row >> 3)) & 1 else 1
                expectation += actual[row][column].conj() * actual[row][column] * sign
            self.assertEqual(expectation, CO)

    def test_exact_cyclotomic_p1_and_p2_sector_receipts(self):
        self.assertEqual(_xor_carry_expectation(1), (1, 0, 1, 0))
        self.assertEqual(_xor_carry_expectation(2),
                         (1, 0, Fraction(57, 128), Fraction(33, 256)))
        self.assertEqual(_xor_carry_expectation(2, transpose=True),
                         (1, 0, Fraction(29, 128), Fraction(21, 128)))


if __name__ == '__main__':
    unittest.main()
