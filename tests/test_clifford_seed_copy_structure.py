"""Exact witnesses for the nested-cube Clifford-seed copy obstruction.

The proof is in CLIFFORD_SEED_COPY_STRUCTURE.md.  Sparse Pauli arithmetic
below evaluates a complete operator on 128 carrier wires, without treating
its exponentially large computational basis as initialized or sampled.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, Q2, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_clifford_orbit_query_rank import _compose, _degree, _pack


def _pauli_product(left, right):
    """Multiply normal X^x Z^z expansions, retaining literal coefficients."""
    answer = {}
    for (x, z), value in left.items():
        for (a, b), other in right.items():
            label = x ^ a, z ^ b
            sign = -1 if (z & a).bit_count() & 1 else 1
            answer[label] = answer.get(label, CZ) + value * other * sign
    return {label: value for label, value in answer.items() if value != CZ}


class CliffordSeedCopyStructureTests(unittest.TestCase):
    def test_entry_degree_two_can_require_group_degree_four(self):
        # Multiplication by the two variables in F2[n1,n2]/(n1^2,n2^2).
        identity = (1, 2, 4, 8)
        n1, n2 = (2, 0, 8, 0), (4, 8, 0, 0)
        self.assertEqual(_compose(n1, n1), (0,) * 4)
        self.assertEqual(_compose(n2, n2), (0,) * 4)
        self.assertEqual(_compose(n1, n2), _compose(n2, n1))
        values = []
        for x in range(32):
            a, b = (x & 3) == 3, (x & 12) == 12
            values.append(tuple(i ^ (u if a else 0) ^ (v if b else 0)
                                for i, u, v in zip(identity, n1, n2)))
        self.assertEqual(_degree([_pack(value) for value in values]), 2)
        for a in values:
            self.assertEqual(_compose(a, a), identity)
            for b in values:
                self.assertEqual(_compose(a, b), _compose(b, a))
        product = identity
        for value in values[:16]:
            product = _compose(product, value)
        expected = tuple(a ^ b for a, b in zip(identity, _compose(n1, n2)))
        self.assertEqual(product, expected)
        self.assertNotEqual(product, identity)
        # A fifth cube direction cancels the genuine degree-four factor.
        for value in values[16:]:
            product = _compose(product, value)
        self.assertEqual(product, identity)

    def test_four_translated_paulis_can_have_literal_negative_scalar(self):
        # Q_(a,b)=H0^a H1^b is an additive, degree-one Clifford family.
        # Its symplectic hypothesis alone must not silently set phases to +1.
        pauli = native_columns(2, [('H', 0), ('Z', 0), ('H', 0),
                                   ('H', 1), ('Z', 1), ('H', 1)])
        product = eye(4)
        for address in range(4):
            query = native_columns(2, [('H', bit) for bit in range(2)
                                       if address >> bit & 1])
            translated = multiply(multiply(query, pauli), adjoint(query))
            product = multiply(product, translated)
        self.assertEqual(product, scale(eye(4), -1))

    def test_complete_128_address_field_word_has_nonzero_output_parity(self):
        p, k = 7, 64
        root = CO * Q2(0, Fraction(1, 2))
        identity = {(0, 0): CO}
        # R=H on u0.  For every degree-64 polynomial field basis,
        # M_x e0=x for seven-bit x: no reduction is encountered.
        # Thus R_x=(X_u0 X_v^x + Z_u0)/sqrt(2), as a full operator.
        orbit = [{(1 | (x << k), 0): root, (0, 1): root}
                 for x in range(1 << p)]
        axes = [(1 | (x << k), 1) for x in range(1 << p)]
        for x, z in axes:
            for a, b in axes:
                self.assertEqual(((x & b).bit_count() + (z & a).bit_count()) % 2, 0)
        # These pairwise commuting transvection axes certify the premise,
        # although the actual Clifford reflections do not all commute.
        self.assertNotEqual(_pauli_product(orbit[0], orbit[1]),
                            _pauli_product(orbit[1], orbit[0]))
        product, output_parity, count = identity, 0, 0
        for translate in (0, 32, 64, 96):
            block = identity
            for low in range(32):
                address = translate ^ low
                self.assertEqual(_pauli_product(orbit[address], orbit[address]), identity)
                block = _pauli_product(block, orbit[address])
                output_parity ^= 1 << address
                count += 1
            # This example is a scalar special case of the Pauli conclusion.
            self.assertEqual(block, identity)
            product = _pauli_product(product, block)
        self.assertEqual(product, identity)
        self.assertEqual(count, 128)
        self.assertEqual(output_parity, (1 << 128) - 1)
        self.assertEqual(Fraction(2, count), Fraction(1, 64))


if __name__ == '__main__':
    unittest.main()
