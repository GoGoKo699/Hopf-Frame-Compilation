"""Exact finite checks for specified probe-asymmetric source repairs.

These fixtures reject three explicit p=3 identities; they do not establish
all-size source exclusions. Proof home: COUPLED_FINAL_ENCODER.md, Section 10.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CQ2, Q2
from tests.test_common_probe_source import (
    CO, CI, _native, _diagonal, _d_phase, _m_map, _syndrome,
)
from tests.test_joint_lookup_pauli_source import _inverse


def _relative_zero_column(repair):
    p, m, address = 3, 8, 0
    t = ([('H', j) for j in range(m)]
         + [('T', j) for j in range(m)]
         + [('H', j) for j in range(m)])
    upper = sum(1 << j for j in range(m) if j & 1)
    phases = [CO, CI, -CO, -CI]
    state = {0: CO}
    # E_+=T†D†TD; E_-=V†T†DVT D†.
    # F=E_+ Z_pi E_-†. The middle D Z_pi D is R_s.
    for op in ['v', 't', 'dd', 'vd', 'td', 'r', 't', 'dd', 'td']:
        if op == 't':
            state = _native(t, state)
        elif op == 'td':
            state = _native(_inverse(t), state)
        elif op in ['v', 'vd']:
            sign = 1 if op == 'v' else -1
            if repair == 'parity':
                state = _diagonal(state, lambda y:
                                  -CO if y.bit_count() & 1 else CO)
            elif repair == 'm_map':
                state = {_m_map(y, p): z for y, z in state.items()}
            elif repair == 'paired_shear':
                state = {
                    y ^ ((y & upper) >> 1):
                    z * phases[(sign * (y & upper).bit_count()) % 4]
                    for y, z in state.items()
                }
            else:
                raise ValueError(repair)
        elif op == 'dd':
            state = _diagonal(state, lambda y:
                              _d_phase(y, p, address).conj())
        else:
            state = _diagonal(
                state,
                lambda y: (-CO if y.bit_count() & 1 else CO)
                * _d_phase(y, p, address) * _d_phase(y, p, address),
            )
    return state[0]


def _rank(vectors):
    pivots = {}
    for value in vectors:
        while value:
            bit = value.bit_length() - 1
            if bit in pivots:
                value ^= pivots[bit]
            else:
                pivots[bit] = value
                break
    return len(pivots)


def _syndrome_involution_basis(p):
    """Columns u_i,v_i,k_j exposing the actual M permutation."""
    m = 1 << p
    u = [(1 << m) - 1]
    u += [sum(1 << j for j in range(m) if (j >> i) & 1)
          for i in range(p)]
    v = [1] + [1 ^ (1 << (1 << i)) for i in range(p)]
    extra = []
    span = list(u)
    for j in range(m):
        k = 1 << j
        for i in range(p):
            if (j >> i) & 1:
                k ^= 1 << (1 << i)
        if not (j.bit_count() & 1):
            k ^= 1
        if _rank(span + [k]) > len(span):
            span.append(k)
            extra.append(k)
    return u, v, extra


class ProbeAsymmetricTrialTests(unittest.TestCase):
    def test_three_relative_fibers_have_exact_nonreal_diagonal(self):
        expected = {
            'parity': CQ2(
                Q2(Fraction(-7919, 65536), Fraction(747, 8192)),
                Q2(Fraction(-21, 2048), Fraction(-21, 4096))),
            'm_map': CQ2(
                Q2(Fraction(-8115, 65536), Fraction(-1743, 16384)),
                Q2(Fraction(423, 8192), Fraction(-37, 16384))),
            'paired_shear': CQ2(
                Q2(Fraction(-13159, 65536), Fraction(-4199, 65536)),
                Q2(Fraction(343, 32768), Fraction(-1213, 65536))),
        }
        for repair, want in expected.items():
            with self.subTest(repair=repair):
                got = _relative_zero_column(repair)
                self.assertEqual(got, want)
                self.assertNotEqual(got.imag, Q2())

    def test_controlled_m_has_disjoint_cnot_coordinates(self):
        for p in range(3, 7):
            with self.subTest(p=p):
                m, q = 1 << p, p + 1
                u, v, extra = _syndrome_involution_basis(p)
                self.assertEqual(len(u), q)
                self.assertEqual(len(extra), m - 2 * q)
                self.assertEqual(_rank(u + v + extra), m)
                for left, right in zip(u, v):
                    self.assertEqual(_m_map(left, p), left)
                    self.assertEqual(_m_map(right, p), right ^ left)
                for vector in extra:
                    self.assertEqual(_m_map(vector, p), vector)
                differences = [_m_map(1 << j, p) ^ (1 << j)
                               for j in range(m)]
                self.assertEqual(_rank(differences), q)
                for vector in differences:
                    self.assertEqual(_m_map(vector, p), vector)

    def test_fixed_syndrome_target_retains_top_degree(self):
        # Exact restrictions to row 0 and the basis rows challenge the
        # coefficient used in the all-p analytic degree argument.
        for p in range(2, 7):
            for address in range(1 << p):
                truth = []
                for variables in range(1 << (p + 1)):
                    bank = variables & 1
                    for i in range(p):
                        if (variables >> (i + 1)) & 1:
                            bank |= 1 << (1 << i)
                    row = address ^ _syndrome(bank, p)
                    truth.append((bank.bit_count() ^ (bank >> row)) & 1)
                for bit in range(p + 1):
                    for mask in range(len(truth)):
                        if (mask >> bit) & 1:
                            truth[mask] ^= truth[mask ^ (1 << bit)]
                self.assertEqual(truth[-1], 1, (p, address))


if __name__ == '__main__':
    unittest.main()
