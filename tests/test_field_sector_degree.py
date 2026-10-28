"""Independent finite checks of the sector-degree barrier and its scope.

The independent group-algebra and escape checks are supplemented by exact
active-coset and local-affinity identities in the retained sector tableaux.
"""
import unittest
from fractions import Fraction

from tests.exact_matrices import CQ2, CO, CZ, adjoint, eye, multiply as exact_multiply, scale
from tests.field_sector_tableaux import (BitslicedTableau, recoupler_word,
                                        sector_tableau, top_derivative_identity)


def identity(n):
    return [1 << i for i in range(n)]


def add(a, b):
    return [x ^ y for x, y in zip(a, b)]


def multiply(a, b):
    out = []
    for row in a:
        value = 0
        while row:
            low = row & -row
            value ^= b[low.bit_length() - 1]
            row ^= low
        out.append(value)
    return out


def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            bit = row.bit_length() - 1
            if bit not in pivots:
                pivots[bit] = row
                break
            row ^= pivots[bit]
    return len(pivots)


def translation(p, x):
    return [1 << (row ^ x) for row in range(1 << p)]


def span(basis):
    values = [0]
    for vector in basis:
        values += [value ^ vector for value in values]
    return values


def complement(p, vector):
    basis = [vector]
    rest = []
    for bit in range(p):
        value = 1 << bit
        if rank(basis + [value]) > len(basis):
            basis.append(value)
            rest.append(value)
    assert len(basis) == p
    return rest


def anf_degree(values):
    coefficients = list(values)
    variables = (len(values) - 1).bit_length()
    assert len(values) == 1 << variables
    for bit in range(variables):
        for mask in range(1 << variables):
            if mask >> bit & 1:
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return max((mask.bit_count() for mask, value in enumerate(coefficients)
                if value), default=-1)


class FieldSectorDegreeTests(unittest.TestCase):
    def test_regular_group_norm_and_independent_deviations(self):
        for p in range(1, 7):
            n = 1 << p
            unit = identity(n)
            total = [0] * n
            for x in range(n):
                total = add(total, translation(p, x))
            self.assertEqual(total, [(1 << n) - 1] * n)
            self.assertEqual(rank(total), 1)
            if p >= 2:
                product = multiply(add(translation(p, 1), unit),
                                   add(translation(p, 2), unit))
                self.assertEqual(rank(product), n // 4)
                self.assertTrue(any(product))

    def test_every_nonzero_kernel_direction_and_complement(self):
        for p in range(1, 6):
            n = 1 << p
            unit = identity(n)
            wanted = [(1 << n) - 1] * n
            for b in range(1, n):
                directions = complement(p, b)
                coefficient = [0] * n
                deviation = add(translation(p, b), unit)
                for t in span(directions):
                    coefficient = add(coefficient,
                                      multiply(translation(p, t), deviation))
                self.assertEqual(coefficient, wanted, (p, b))
                # Multiplication by a nontrivial invertible core change cannot
                # annihilate the top-degree coefficient. Elementary row/column
                # operations suffice for a direct independent finite check.
                change = identity(n)
                if n > 1:
                    change[0] ^= change[1]
                self.assertTrue(any(multiply(coefficient, change)))

    def test_full_rank_helper_case_needs_degree_p(self):
        for p in range(1, 7):
            n = 1 << p
            # Every entry of the regular translation family is one address
            # indicator, of degree p in the p active helper-shift variables.
            for output in (0, n - 1):
                for source in (0, n - 1):
                    self.assertEqual(anf_degree([
                        int((source ^ x) == output) for x in range(n)]), p)

    def test_square_zero_source_assumption_is_essential(self):
        # If the source itself is the regular bank, E=I and Lambda=0 is a
        # constant Clifford-sector solution. Its deviations violate N_b N_c=0.
        p = 3
        n = 1 << p
        unit = identity(n)
        source = {x: translation(p, x) for x in range(n)}
        for x in range(n):
            self.assertEqual(multiply(unit, source[x]),
                             multiply(translation(p, x), unit))
        self.assertTrue(any(multiply(add(source[1], unit),
                                     add(source[2], unit))))

    def test_nonlinear_helper_motion_invalidates_affine_stage_count(self):
        # Reversible helper Toffoli: r2 -> r2 + r0*r1. It is its own inverse.
        # Between its two occurrences, a single core Clifford shear is
        # controlled AFFINELY by the current r2. In the original helper frame
        # the completed block has a quadratic coefficient, despite containing
        # only one affine core stage. The helper motion must be charged and is
        # outside the fixed-helper-frame stage hypothesis.
        values = []
        permutation = []
        for r in range(8):
            shifted = r ^ (((r & 1) * ((r >> 1) & 1)) << 2)
            permutation.append(shifted)
            values.append((shifted >> 2) & 1)
        self.assertEqual(sorted(permutation), list(range(8)))
        self.assertTrue(all(permutation[permutation[r]] == r for r in range(8)))
        self.assertEqual(anf_degree([(r >> 2) & 1 for r in range(8)]), 1)
        self.assertEqual(anf_degree(values), 2)

    def test_extra_target_spectators_cannot_cancel_the_bank_block(self):
        p = 3
        n = 1 << p
        dimension = n + 2
        unit = identity(dimension)
        queries = {}
        for x in range(n):
            # A moving residual 2-wire swap, on an arbitrary address character.
            residual_bit = (x & 3).bit_count() & 1
            queries[x] = (translation(p, x)
                          + [1 << (n + (j ^ residual_bit)) for j in range(2)])
        b = 4
        deviation = add(queries[b], unit)
        total = [0] * dimension
        for t in span([1, 2]):
            total = add(total, multiply(queries[t], deviation))
        self.assertEqual(total[:n], [(1 << n) - 1] * n)
        self.assertTrue(any(total))


    def test_active_coset_top_coefficient_and_local_affinity(self):
        for p in (2, 3, 4):
            self.assertTrue(top_derivative_identity(p))
            table = sector_tableau(p)
            self.assertEqual(max(table.degree_profile()['histograms']['symplectic']), p - 1)
            ell = 1 << (p - 2)
            offset = ell - 1
            local_word = [((*gate[:-1], gate[-1] >> offset)
                           if gate[0] in ('ZC', 'CZC') else gate)
                          for gate in recoupler_word(p)]
            local = BitslicedTableau(1 << p, ell).apply(local_word)
            self.assertEqual(max(local.degree_profile()['histograms']['symplectic']), 1)

    def test_signed_pauli_gap_does_not_fix_scalar_sector_phase(self):
        # Orthogonal Pauli directions attain the sqrt(2) gap exactly.
        x = [[CZ, CO], [CO, CZ]]
        z = [[CO, CZ], [CZ, -CO]]
        difference = [[a - b for a, b in zip(row, other)]
                      for row, other in zip(x, z)]
        self.assertEqual(exact_multiply(difference, difference), scale(eye(2), CQ2(2)))
        self.assertEqual(exact_multiply(scale(x, CQ2(2)), scale(x, CQ2(2))),
                         scale(eye(2), CQ2(4)))
        # A nontrivial rational unit scalar is closer than 1/sqrt(2), but its
        # Pauli action is exactly identical. Robust rounding therefore gives
        # the symplectic cocycle; it cannot set sector scalars to +1.
        phase = CQ2(Fraction(15, 17), Fraction(8, 17))
        self.assertEqual(phase.abs2(), 1)
        delta = phase - 1
        self.assertEqual(delta.abs2(), Fraction(4, 17))
        self.assertLess(delta.abs2(), Fraction(1, 2))
        unitary = scale(eye(2), phase)
        self.assertNotEqual(unitary, eye(2))
        self.assertEqual(exact_multiply(exact_multiply(unitary, x), adjoint(unitary)), x)
        self.assertEqual(exact_multiply(exact_multiply(unitary, z), adjoint(unitary)), z)


if __name__ == '__main__':
    unittest.main()
