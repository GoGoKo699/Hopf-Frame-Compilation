"""Exact guards for flat inverse kernels and their normalized forward Gram.

These checks retain arbitrary dirty target columns and literal native
phases.  They are finite fixtures, not a scalable inverse-column source.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, CQ2, Q2, adjoint, eye, multiply, native_columns,
)
from tests.test_majorana_walsh_query import _inverse as word_inverse
from tests.test_nilpotent_pencil_observable_source import (
    _field_matrix, _inverse, _mul, _power, _product,
)


def parity(value):
    return value.bit_count() & 1


def hadamard_sum(values):
    result = list(values)
    stride = 1
    while stride < len(result):
        for start in range(0, len(result), 2 * stride):
            for offset in range(stride):
                left, right = start + offset, start + offset + stride
                a, b = result[left], result[right]
                result[left], result[right] = a + b, a - b
        stride *= 2
    return result


def anf_degree(values):
    coefficients = list(values)
    bit = 1
    while bit < len(coefficients):
        for index in range(len(coefficients)):
            if index & bit:
                coefficients[index] ^= coefficients[index ^ bit]
        bit <<= 1
    return max((i.bit_count() for i, value in enumerate(coefficients) if value),
               default=0)


def apply_binary(rows, vector):
    return sum(parity(row & vector) << i for i, row in enumerate(rows))


def ry_t(signal, inverse=False, cube=False):
    # L=SH takes Z to Y.  In chronological order L^dagger,T,L.
    middle = [('S', signal)] if cube else []
    return ([('S', signal)] * 3 + [('H', signal)] + middle
            + [('TDG' if inverse else 'T', signal), ('H', signal), ('S', signal)])


class BentGramSourceTests(unittest.TestCase):
    def test_exact_mm_dual_really_contains_higher_degree_inverse(self):
        # A quadratic triangular permutation with inverse coordinate degree3.
        # pi_3=b3+b0*b1 and pi_4=b4+b2*b3; its inverse fourth
        # coordinate has the cubic b0*b1*b2 term.
        n, size = 5, 32
        permutation = [b ^ ((((b >> 0) & 1) & ((b >> 1) & 1)) << 3)
                         ^ ((((b >> 2) & 1) & ((b >> 3) & 1)) << 4)
                       for b in range(size)]
        self.assertEqual(len(set(permutation)), size)
        inverse = [permutation.index(u) for u in range(size)]
        self.assertEqual(anf_degree([(b >> 4) & 1 for b in inverse]), 3)
        self.assertEqual(anf_degree([(b >> 4) & 1 for b in permutation]), 2)
        # Index a+2^n*b, Fourier label u+2^n*v.  Integer Walsh
        # sums verify every dual coefficient with no rounded amplitudes.
        phase = [(-1) ** parity(a & permutation[b])
                 for b in range(size) for a in range(size)]
        transform = hadamard_sum(phase)
        for u in range(size):
            for v in range(size):
                self.assertEqual(transform[u + size * v],
                                 size * (-1) ** parity(v & inverse[u]))
        # The normalized Fourier-sandwich kernel has denominator d,
        # not one: there are d^2 entries of squared modulus 1/d^2.
        self.assertEqual(sum(Fraction(1, size * size)
                             for _ in range(size * size)), 1)

    def test_literal_normalized_echo_all_columns_and_actual_inverse(self):
        # B0,B1 are arbitrary data; A1 is arbitrary dirty target;
        # A0=s is the sole clean signal.  The query changes A1 too,
        # so cancellation is checked rather than assumed.
        signal, width = 3, 4
        query = [('CX', 0, signal), ('CX', 1, signal), ('CX', 1, 2)]
        source = (word_inverse(query) + ry_t(signal, inverse=True)
                  + query + ry_t(signal, cube=True))
        actual = native_columns(width, source)
        self.assertEqual(native_columns(width, word_inverse(source)), adjoint(actual))
        self.assertEqual(multiply(adjoint(actual), actual), eye(1 << width))
        half = CQ2(Q2(Fraction(1, 2)))
        alpha = CQ2(Q2(Fraction(1, 2)), Q2(Fraction(1, 2)))
        root_half = Q2(0, Fraction(1, 2))
        omega_inverse = CQ2(root_half, -root_half)
        # R^2=S_Y on the q=0 branch; omega^-1 Y on q=1.
        r_squared = native_columns(1, ry_t(0) + ry_t(0))
        y = [[CZ, CQ2(Q2(0), Q2(-1))],
             [CQ2(Q2(0), Q2(1)), CZ]]
        for column in range(16):
            q = parity(column & 3)
            for row in range(16):
                if row & 7 != column & 7:
                    expected = CZ
                elif q:
                    expected = omega_inverse * y[row >> 3][column >> 3]
                else:
                    expected = r_squared[row >> 3][column >> 3]
                self.assertEqual(actual[row][column], expected)
        block = [row[:8] for row in actual[:8]]
        gram = multiply(adjoint(block), block)
        for row in range(8):
            for column in range(8):
                accepted = row == column and parity(column & 3) == 0
                self.assertEqual(block[row][column], alpha if accepted else CZ)
                self.assertEqual(gram[row][column], half if accepted else CZ)
        self.assertEqual(sum(gate[0] in {'T', 'TDG'} for gate in source), 2)

    def test_gf4_pencil_forward_column_differs_on_complete_bank(self):
        p, m, n = 2, 4, 8
        basis = _field_matrix([[_power(u, j) for u in range(m)] for j in range(m)])
        basis_inverse = _inverse(basis)
        mismatches = 0
        for x in range(m):
            standard = [[int(i == j) ^ (_mul(x, 1) if i == j + 1 else 0)
                         for j in range(m)] for i in range(m)]
            pencil = _product(_product(basis_inverse, _field_matrix(standard)), basis)
            forward = apply_binary(pencil, 1)
            inverse = apply_binary(_inverse(pencil), 1)
            self.assertEqual(inverse, 1 << (p * x))
            if x == 0:
                self.assertEqual(forward, inverse)
            else:
                self.assertNotEqual(forward, inverse)
                mismatches += 1
                # C xor D may be any value, including private carrier
                # coordinates.  A differing character gives Gram gap1/2.
                witness = next(1 << bit for bit in range(n)
                               if ((forward ^ inverse) >> bit) & 1)
                trial = Fraction(1 + (-1) ** parity(forward & witness), 4)
                target = Fraction(1 + (-1) ** parity(inverse & witness), 4)
                self.assertEqual(abs(trial - target), Fraction(1, 2))
        self.assertEqual(mismatches, 3)

    def test_one_clean_clifford_compression_singular_alternatives(self):
        # Data wire0, signal wire1, with signal0 embedding columns0,1.
        examples = [
            ([('CX', 0, 1)], [CO, CZ]),  # Pauli half-projector
            ([('H', 1)], [CQ2(Q2(Fraction(1, 2)))] * 2),
            ([('S', 1), ('S', 1)], [CO, CO]),
            ([('H', 1), ('S', 1), ('S', 1), ('H', 1)], [CZ, CZ]),
        ]
        for word, diagonal in examples:
            matrix = native_columns(2, word)
            self.assertEqual(matrix, adjoint(matrix))
            self.assertEqual(multiply(matrix, matrix), eye(4))
            block = [row[:2] for row in matrix[:2]]
            gram = multiply(adjoint(block), block)
            self.assertEqual(gram, [[diagonal[r] if r == c else CZ
                                     for c in range(2)] for r in range(2)])

    def test_lookup_projector_pauli_orbit_gap_is_at_least_root_half(self):
        # p=1, address bit0, C0,C1 bits1,2, D0,D1 bits3,4.
        width, size = 5, 32
        q = [((z >> (1 + (z & 1))) ^ (z >> (3 + (z & 1)))) & 1
             for z in range(size)]
        self.assertEqual(anf_degree(q), 2)
        for x_mask in range(size):
            for z_mask in range(size):
                if not x_mask and not z_mask:
                    continue
                if not x_mask:
                    for sign in (0, 1):
                        pauli_projector = [sign ^ parity(z_mask & z) for z in range(size)]
                        self.assertNotEqual(pauli_projector, q)
                else:
                    # For each two-element orbit, Q entries equal ->
                    # norm^2=1; unequal -> norm^2=1/2.  Pauli phases
                    # affect no singular value of this 2x2 block.
                    squared = max(Fraction(1) if q[z] == q[z ^ x_mask]
                                  else Fraction(1, 2) for z in range(size))
                    self.assertGreaterEqual(squared, Fraction(1, 2))


if __name__ == '__main__':
    unittest.main()
