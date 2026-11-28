"""Exact native baselines for signed MOD3-controlled two-bit Clifford banks.

Every address, block input and optional external signal is arbitrary; the
construction uses no clean or borrowed helpers. Small exact fixtures retain
literal phases and actual positive-S inverses. The algebraic-degree screen
concerns reversible-classical multiplicative depth only, not native T-depth
when Hadamards and interference are permitted. These are components, not a
complete accepted source or a new asymptotic depth upper bound.
The four-label eigenbasis fixture separately checks the two-angle identity
in floating-point arithmetic.
"""
from itertools import product
import unittest

import numpy as np

from tests.exact_matrices import CO, adjoint, native_columns, zeros
from tests.test_fixed_row_transport_readout import _positive_s_schedule
from tests.test_phase_signal_query_limits import _counts
from tests.test_t_depth import _ccz_batch, _fredkin_batch, _inverse, _word


def _ccx_batch(control, pairs):
    walls = [('H', target) for _, target in pairs]
    return ([('C', walls)] + _ccz_batch(control, pairs)
            + [('C', walls)])


def _mod3_native_batch(p, signs, controlled=False):
    """U(x)=tensor_j A**sum_i signs[j][i]*x_i, A(s,t)=(t,s xor t)."""
    assert p >= 1 and signs
    assert all(len(row) == p and set(row) <= {-1, 1} for row in signs)
    offset, q = int(controlled), len(signs)
    pairs = [(offset+p+2*j, offset+p+2*j+1) for j in range(q)]
    word = []
    for i in range(p):
        # A is CX(s;t), CX(t;s) chronologically; reversing gives A^-1.
        oriented = [(s, t) if signs[j][i] == 1 else (t, s)
                    for j, (s, t) in enumerate(pairs)]
        word += _ccx_batch(offset+i, oriented)
        word += _ccx_batch(offset+i, [(t, s) for s, t in oriented])
    if controlled:
        reflection = _fredkin_batch(0, pairs)
        # R U R=U^-1 and U^3=I. In chronology U^-1,R^a,U,R^a,
        # signal one gives R U R U^-1=U^-2=U; zero gives identity.
        word = _inverse(word) + reflection + word + reflection
    return offset+p+2*q, _positive_s_schedule(word)


def _act_pair(s, t, exponent):
    for _ in range(exponent % 3):
        s, t = t, s ^ t
    return s, t


def _mod3_columns(p, signs, controlled=False, inverse=False, columns=None):
    offset, q = int(controlled), len(signs)
    width = offset+p+2*q
    columns = list(range(1 << width)) if columns is None else list(columns)
    expected = zeros(1 << width, len(columns))
    for index, column in enumerate(columns):
        output = column
        if not controlled or column & 1:
            for j, row in enumerate(signs):
                exponent = sum(sign*((column >> (offset+i)) & 1)
                               for i, sign in enumerate(row))
                if inverse:
                    exponent = -exponent
                left, right = offset+p+2*j, offset+p+2*j+1
                s, t = _act_pair((column >> left) & 1,
                                 (column >> right) & 1, exponent)
                output = ((output & ~((1 << left) | (1 << right)))
                          | (s << left) | (t << right))
        expected[output][index] = CO
    return expected


def _anf(values):
    coefficients = list(values)
    for bit in range((len(values)-1).bit_length()):
        for mask in range(len(values)):
            if mask & (1 << bit):
                coefficients[mask] ^= coefficients[mask ^ (1 << bit)]
    return coefficients


class Mod3NativeBatchTests(unittest.TestCase):
    def test_two_angle_eigenbasis_and_signed_cycle_powers(self):
        hadamard_pair = np.eye(4, dtype=complex)
        hadamard_pair[1:3, 1:3] = np.array([[1, 1], [1, -1]]) / np.sqrt(2)

        theta = np.arccos(1 / np.sqrt(3))
        private_phase = np.kron(np.eye(2), np.diag([1, np.exp(1j * theta)]))
        private_phase *= np.exp(0.37j)
        parity = np.eye(4)[:, [1, 0, 2, 3]]  # Flip the low bit when the high bit is zero.
        signed_phase = parity @ private_phase.conj().T @ parity @ private_phase
        y_basis = np.kron(np.eye(2), np.array([[1, 1], [1j, -1j]]) / np.sqrt(2))
        rotation = y_basis @ signed_phase @ y_basis.conj().T
        basis = hadamard_pair @ rotation @ np.diag([1, -1, 1j, 1]) @ hadamard_pair

        omega = (-1 + 1j * np.sqrt(3)) / 2
        expected = np.eye(4, dtype=complex)
        expected[:3, :3] = np.array([[1, 1, 1], [1, omega, omega.conjugate()],
                                     [1, omega.conjugate(), omega]]) / np.sqrt(3)
        np.testing.assert_allclose(basis, expected, atol=1e-14, rtol=0)

        cycle = np.eye(4)[:, [1, 2, 0, 3]]
        public_phase = np.exp(-0.23j) * np.array([1, omega])
        for mismatch in (0, 1):
            row_phases = []
            for high, low in product((0, 1), repeat=2):
                exponent = (mismatch ^ high) - (mismatch ^ low) - high + low
                self.assertEqual(exponent % 3, mismatch * (high - low) % 3)
                row_phases.append(public_phase[mismatch ^ high]
                                  * public_phase[mismatch ^ low].conjugate()
                                  * public_phase[high].conjugate() * public_phase[low])
            np.testing.assert_allclose(basis @ np.diag(row_phases) @ basis.conj().T,
                                       np.linalg.matrix_power(cycle, mismatch),
                                       atol=1e-14, rtol=0)

        for power in (-1, 0, 1, 2):
            phases = np.diag([1, omega**(-power), omega**power, 1])
            np.testing.assert_allclose(basis @ phases @ basis.conj().T,
                                       np.linalg.matrix_power(cycle, power),
                                       atol=1e-14, rtol=0)

    def test_one_block_all_columns_signs_and_actual_inverse(self):
        for p in (1, 2):
            for row, controlled in product(product((-1, 1), repeat=p),
                                           (False, True)):
                signs = (row,)
                with self.subTest(p=p, row=row, controlled=controlled):
                    width, schedule = _mod3_native_batch(p, signs, controlled)
                    actual = native_columns(width, _word(schedule))
                    self.assertEqual(actual, _mod3_columns(p, signs, controlled))
                    reverse = _positive_s_schedule(_inverse(schedule))
                    inverse = native_columns(width, _word(reverse))
                    self.assertEqual(inverse, adjoint(actual))
                    self.assertEqual(inverse, _mod3_columns(
                        p, signs, controlled, inverse=True))

    def test_two_signed_blocks_shared_address_and_signal(self):
        p, signs = 2, ((1, 1), (1, -1))
        for controlled in (False, True):
            width, schedule = _mod3_native_batch(p, signs, controlled)
            offset = int(controlled)
            # All address and signal values, including the zero and all-one
            # block inputs; the other patterns distinguish both orientations.
            columns = [((address | (blocks << p)) << offset) | signal
                       for address, blocks, signal in product(
                           range(4), (0, 1, 6, 9, 15), range(1 << offset))]
            for inverse, orientation in (
                    (False, schedule),
                    (True, _positive_s_schedule(_inverse(schedule)))):
                self.assertEqual(native_columns(width, _word(orientation), columns),
                                 _mod3_columns(p, signs, controlled, inverse, columns))

    def test_native_layers_complete_counts_and_no_helpers(self):
        cases = [(p, 1 << (p-1)) for p in range(1, 9)]
        cases += [(p, q) for p, q in product((1, 2, 5), (1, 2, 3, 4, 16, 64))]
        for p, q in cases:
            signs = tuple(tuple(1 if i == 0 or not ((j >> (i-1)) & 1) else -1
                                for i in range(p)) for j in range(q))
            chi, k, back = q % 2, (q//2) % 4, (-(q//2)) % 4
            for controlled in (False, True):
                width, schedule = _mod3_native_batch(p, signs, controlled)
                self.assertEqual(width, p+2*q+int(controlled))
                for reverse, orientation in (
                        (False, schedule),
                        (True, _positive_s_schedule(_inverse(schedule)))):
                    residue = back if reverse else k
                    if controlled:
                        expected_t = (2*p+1)*(12*q+2*chi)
                        expected_g = 48*p*q+28*q+2*p*(k+back)+2*residue
                        common_g = 48*p*q+28*q+8*p+6
                        expected_depth = 16*p+8
                    else:
                        expected_t = p*(12*q+2*chi)
                        expected_g = 24*p*q+2*p*residue
                        common_g = 24*p*q+6*p
                        expected_depth = 8*p
                    self.assertEqual(_counts(_word(orientation)),
                                     (expected_t, expected_g))
                    self.assertLessEqual(expected_g, common_g)
                    self.assertEqual(sum(kind == 'T' for kind, _ in orientation),
                                     expected_depth)
                    for kind, gates in orientation:
                        self.assertTrue(all(0 <= wire < width
                                            for gate in gates for wire in gate[1:]))
                        if kind == 'T':
                            self.assertTrue(gates)
                            self.assertTrue(all(gate[0] in ('T', 'TDG')
                                                for gate in gates))
                            targets = [gate[1] for gate in gates]
                            self.assertEqual(len(targets), len(set(targets)))
                        else:
                            self.assertTrue(all(gate[0] in ('H', 'S', 'CX')
                                                for gate in gates))

    def test_exact_top_degree_for_every_small_signed_mod3_block(self):
        powers = (((1, 0), (0, 1)), ((0, 1), (1, 1)), ((1, 1), (1, 0)))
        for p in range(1, 7):
            for row in product((-1, 1), repeat=p):
                # Over F2, A^sigma+I=A^-sigma. Thus the coefficient of
                # x_1...x_p in the matrix polynomial is A^-sum(sigma).
                top = powers[(-sum(row)) % 3]
                for controlled in (False, True):
                    width = p+2+int(controlled)
                    values = [[], []]
                    for column in range(1 << width):
                        exponent = sum(sign*((column >> i) & 1)
                                       for i, sign in enumerate(row))
                        if controlled and not ((column >> (p+2)) & 1):
                            exponent = 0
                        output = _act_pair((column >> p) & 1,
                                           (column >> (p+1)) & 1, exponent)
                        for coordinate in range(2):
                            values[coordinate].append(output[coordinate])
                    for coordinate in range(2):
                        coefficients = _anf(values[coordinate])
                        self.assertEqual(max(mask.bit_count() for mask, value
                                             in enumerate(coefficients) if value),
                                         p+1+int(controlled))
                        for work in range(2):
                            mask = ((1 << p)-1) | (1 << (p+work))
                            if controlled:
                                mask |= 1 << (p+2)
                            self.assertEqual(coefficients[mask], top[coordinate][work])


if __name__ == '__main__':
    unittest.main()
