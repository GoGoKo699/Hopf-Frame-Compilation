"""Exact full-space guards for four-path acceptance polynomials.

The native example is a three-qubit diagnostic, not an asymptotic lookup
source. No comparison quotients by a scalar phase or drops rejected outputs.
"""
from fractions import Fraction
from itertools import combinations, product
import unittest

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns, scale,
)


HALF = CQ2(Fraction(1, 2))
QUARTER = CQ2(Fraction(1, 4))
INVROOT = CQ2(Q2(0, Fraction(1, 2)))
ROOT = CQ2(Q2(0, 1))
OMEGA = INVROOT * (CO + CI)


def _sum(*matrices):
    result = scale(matrices[0], CZ)
    for matrix in matrices:
        result = add(result, matrix)
    return result


def _commutator(a, b):
    return add(multiply(a, b), scale(multiply(b, a), -CO))


def _reverse(word):
    result = []
    for name, *qubits in reversed(word):
        if name == 'S':
            result.extend([('S', *qubits)] * 3)
        else:
            result.append(({'T': 'TDG', 'TDG': 'T'}.get(name, name), *qubits))
    return result


def _ccx_bands(a, b, d):
    # Three independent parities in the middle T band, with literal uncompute.
    change = [('CX', b, a), ('CX', a, d), ('CX', d, b)]
    return [('C', [('H', d)]),
            ('T', [('T', a), ('T', b), ('T', d)]),
            ('C', change), ('T', [('TDG', a), ('TDG', b), ('T', d)]),
            ('C', _reverse(change) + [('CX', b, d)]),
            ('T', [('TDG', d)]), ('C', [('CX', b, d), ('H', d)])]


def _source_bands():
    # Data is the low bit; a,b are physical signal bits, in that order.
    d, a, b = 0, 1, 2
    return ([('C', [('H', a), ('H', b), ('CX', a, d),
                    ('H', d), ('CX', b, d), ('H', d)])]
            + _ccx_bands(a, b, d)
            + [('C', [('S', a)]), ('T', [('T', b)]),
               ('C', [('H', a), ('H', b)])])


def _source():
    return [gate for _, gates in _source_bands() for gate in gates]


def _accepted(source, size):
    return [row[:size] for row in source[:size]]


def _paths(a, b):
    return [eye(len(a)), scale(a, OMEGA), scale(b, CI),
            scale(multiply(multiply(b, a), b), CI * OMEGA)]


def _det2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


class FourPathAcceptedSourceTests(unittest.TestCase):
    def test_calibrated_gram_full_noncommuting_identity(self):
        a = native_columns(1, [('Z', 0)])
        x = native_columns(1, [('H', 0), ('Z', 0), ('H', 0)])
        choices = [eye(2), scale(eye(2), -CO), a, x,
                   scale(add(a, x), INVROOT)]
        for b in choices:
            with self.subTest(b=b):
                self.assertEqual(multiply(b, b), eye(2))
                bab = multiply(multiply(b, a), b)
                k = scale(add(scale(add(eye(2), a), OMEGA),
                              scale(add(b, bab), OMEGA.conj())), QUARTER)
                expected = scale(_sum(scale(eye(2), 4), scale(a, 2),
                                      multiply(a, b), multiply(b, a),
                                      scale(add(_commutator(b, a),
                                                _commutator(bab, a)), CI)),
                                 CQ2(Fraction(1, 16)))
                g = multiply(adjoint(k), k)
                self.assertEqual(g, expected)
                flatness = multiply(g, add(g, scale(eye(2), -HALF)))
                self.assertEqual(flatness == scale(eye(2), CZ), b == eye(2))

    def test_calibrated_commuting_scalar_edges(self):
        for a, b in product((-1, 1), repeat=2):
            k = (OMEGA * (1 + a) + OMEGA.conj() * (b + a)) * QUARTER
            g = k.conj() * k
            self.assertEqual(g * (g - HALF) == CZ, b == 1)
            if b == 1:
                self.assertEqual(g, CQ2(Fraction(1 + a, 4)))

    def test_robust_gap_uses_an_exact_positive_interval_certificate(self):
        # Coefficients in increasing powers of t. For 0<=t<=1/4 the
        # right factorization is at least 3/16, even if target Q mixes blocks.
        f = Fraction
        lower_square = [f(49, 64), f(7, 8), f(1, 4)]
        upper_square = [f(1, 2), f(1), f(1)]
        difference = [a - b for a, b in zip(lower_square, upper_square)]
        factored = [f(3, 16) + f(1, 4) * f(5, 16),
                    f(1, 4) * f(3, 4) - f(5, 16), -f(3, 4)]
        self.assertEqual(difference, factored)
        # A=-Z,B=-I has G=I/4; arbitrary non-diagonal targets still have
        # squared Gram error I/16, with no block-preservation assumption.
        q = native_columns(1, [('H', 0)])
        error = scale(q, -QUARTER)
        self.assertEqual(multiply(error, error), scale(eye(2), CQ2(f(1, 16))))

    def test_four_paths_have_full_half_gram_and_moving_polar_part(self):
        # Conjugate an anticommuting pair through an arbitrary-work coupling.
        encoder = native_columns(2, [('H', 1), ('T', 1), ('CX', 1, 0)])
        a0 = native_columns(2, [('Z', 0)])
        b0 = native_columns(2, [('H', 0), ('Z', 0), ('H', 0)])
        a = multiply(multiply(encoder, a0), adjoint(encoder))
        b = multiply(multiply(encoder, b0), adjoint(encoder))
        self.assertEqual(add(multiply(a, b), multiply(b, a)), scale(eye(4), CZ))
        k = scale(_sum(*_paths(a, b)), QUARTER)
        q = scale(add(a, scale(multiply(a, b), CI)), INVROOT)
        qp = scale(add(a, scale(multiply(a, b), -CI)), INVROOT)
        self.assertEqual(multiply(q, q), eye(4))
        self.assertEqual(multiply(qp, qp), eye(4))
        self.assertNotEqual(q, qp)
        self.assertEqual(multiply(adjoint(k), k), scale(add(eye(4), q), QUARTER))
        self.assertEqual(multiply(k, adjoint(k)), scale(add(eye(4), qp), QUARTER))
        polar = scale(add(eye(4), scale(b, CI)), INVROOT)
        self.assertEqual(k, scale(multiply(polar, add(eye(4), q)), INVROOT * HALF))
        self.assertEqual(multiply(k, add(eye(4), scale(q, -CO))), scale(eye(4), CZ))
        # Exact trace certifies half the whole private/data space is accepted.
        self.assertEqual(sum(q[j][j] for j in range(4)), CZ)

    def test_wrong_sector_cancellation_is_not_pairwise(self):
        a = native_columns(1, [('Z', 0)])
        b = native_columns(1, [('H', 0), ('Z', 0), ('H', 0)])
        paths = _paths(a, b)
        q = scale(add(a, scale(multiply(a, b), CI)), INVROOT)
        wrong = scale(add(eye(2), scale(q, -CO)), HALF)
        k = scale(_sum(*paths), QUARTER)
        self.assertEqual(multiply(k, wrong), scale(eye(2), CZ))
        for j, l in combinations(range(4), 2):
            self.assertNotEqual(_det2(add(paths[j], paths[l])), CZ)
            self.assertNotEqual(multiply(adjoint(paths[j]), paths[l]), q)

    def test_independent_reflection_split_and_trace_saturation(self):
        ident = eye(4)
        zero = scale(ident, CZ)
        a = native_columns(2, [('Z', 0)])
        x = native_columns(2, [('H', 0), ('Z', 0), ('H', 0)])
        mixing_b = native_columns(2, [('H', 1)])
        block = [[CO if i == j and i < 2 else CZ for j in range(4)]
                 for i in range(4)]
        z2 = native_columns(1, [('Z', 0)])
        x2 = native_columns(1, [('H', 0), ('Z', 0), ('H', 0)])
        noncommuting_d = [[(z2 if i < 2 else x2)[i % 2][j % 2]
                          if i // 2 == j // 2 else CZ for j in range(4)]
                         for i in range(4)]
        same_sign_d = multiply(a, add(scale(block, 2), scale(ident, -CO)))
        # The A,D blocks are invariant, but B creates superpositions across
        # them. Thus the rejection fixture does not assume B preserves them.
        self.assertEqual(_commutator(a, block), zero)
        self.assertEqual(_commutator(noncommuting_d, block), zero)
        self.assertNotEqual(_commutator(a, noncommuting_d), zero)
        self.assertNotEqual(_commutator(mixing_b, block), zero)

        def trace(matrix):
            return sum(matrix[j][j] for j in range(4))

        cases = [(noncommuting_d, mixing_b, False),
                 (same_sign_d, mixing_b, False),
                 (scale(a, -CO), multiply(mixing_b, x), True)]
        for d, b, flat in cases:
            with self.subTest(flat=flat, d=d):
                for reflection in (a, b, d):
                    self.assertEqual(adjoint(reflection), reflection)
                    self.assertEqual(multiply(reflection, reflection), ident)
                s = _sum(ident, scale(a, OMEGA), scale(b, CI),
                         scale(d, CI * OMEGA))
                k = scale(s, QUARTER)
                u = scale(add(ident, scale(b, CI)), INVROOT)
                r = scale(add(a, scale(d, CI)), OMEGA * INVROOT)
                v = multiply(adjoint(u), r)
                t = add(ident, v)
                self.assertEqual(s, scale(add(u, r), ROOT))
                self.assertEqual(k, scale(multiply(u, t), INVROOT * HALF))
                self.assertEqual(multiply(adjoint(u), u), ident)
                self.assertEqual(multiply(adjoint(r), r),
                                 add(ident, scale(_commutator(a, d), CI * HALF)))
                self.assertEqual(trace(multiply(adjoint(r), r)), CQ2(4))
                self.assertEqual(trace(multiply(adjoint(v), v)), CQ2(4))
                self.assertEqual(trace(multiply(adjoint(t), t)),
                                 trace(t) + trace(t).conj())
                g = multiply(adjoint(k), k)
                self.assertEqual(multiply(g, add(g, scale(ident, -HALF))) == zero,
                                 flat)
                if flat:
                    # Extract the projector from the independently summed
                    # paths; exact saturation identifies T with its positive
                    # square root, rather than merely matching singular values.
                    projector = scale(g, 2)
                    self.assertEqual(adjoint(projector), projector)
                    self.assertEqual(multiply(projector, projector), projector)
                    self.assertEqual(trace(projector), CQ2(2))
                    self.assertEqual(multiply(adjoint(t), t), scale(projector, 4))
                    self.assertEqual(trace(t), 2 * trace(projector))
                    self.assertEqual(t, scale(projector, 2))
                    polar = scale(multiply(t, projector), HALF)
                    self.assertEqual(multiply(adjoint(polar), polar), projector)
                    defect = add(polar, scale(projector, -CO))
                    self.assertEqual(trace(multiply(adjoint(defect), defect)), CZ)
                    self.assertEqual(adjoint(v), v)
                    self.assertEqual(multiply(v, v), ident)
                    self.assertEqual(d, scale(a, -CO))
                    self.assertEqual(add(multiply(a, b), multiply(b, a)), zero)

        # The remaining same-sign sector cannot satisfy equation (30), even
        # with this B mixing sectors: its compressed commutator has trace zero
        # against F, whereas the required right side has positive trace four.
        e = multiply(a, add(ident, scale(block, -CO)))
        f = multiply(a, block)
        r = scale(add(a, scale(same_sign_d, CI)), OMEGA * INVROOT)
        self.assertEqual(r, add(e, scale(f, CI)))
        self.assertEqual(multiply(e, f), zero)
        self.assertEqual(multiply(f, f), block)
        lhs = add(add(multiply(e, mixing_b), multiply(mixing_b, e)),
                  scale(_commutator(mixing_b, f), CI))
        compressed = multiply(multiply(block, lhs), block)
        self.assertEqual(trace(multiply(f, compressed)), CZ)
        self.assertEqual(trace(multiply(f, scale(f, 2))), CQ2(4))
        self.assertNotEqual(lhs, scale(f, 2))

    def test_literal_three_qubit_source_block_and_actual_reverse(self):
        word = _source()
        bands = _source_bands()
        self.assertEqual(sum(kind == 'T' for kind, _ in bands), 4)
        self.assertEqual(sum(kind == 'T' for kind, _ in _ccx_bands(1, 2, 0)), 3)
        for kind, gates in bands:
            self.assertTrue(all((gate[0] in ('T', 'TDG')) == (kind == 'T')
                                for gate in gates))
            if kind == 'T':
                self.assertEqual(len({gate[1] for gate in gates}), len(gates))
        self.assertEqual(len(word), 27)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in word), 8)
        self.assertTrue(all(g[0] in ('H', 'S', 'T', 'TDG', 'CX') for g in word))
        source = native_columns(3, word)
        inverse_word = _reverse(word)
        self.assertEqual(len(inverse_word), 29)
        self.assertEqual(sum(g[0] in ('T', 'TDG') for g in inverse_word), 8)
        reverse = native_columns(3, inverse_word)
        self.assertEqual(reverse, adjoint(source))
        self.assertEqual(multiply(reverse, source), eye(8))
        a = native_columns(1, [('Z', 0)])
        b = native_columns(1, [('H', 0), ('Z', 0), ('H', 0)])
        expected = scale(_sum(eye(2), scale(a, ROOT), scale(b, CI)), QUARTER)
        self.assertEqual(_accepted(source, 2), expected)
        # Spectator is arbitrary, not a clean input: compare every joint column.
        extended = native_columns(4, word)
        self.assertEqual(extended,
                         [[source[i & 7][j & 7] if i >> 3 == j >> 3 else CZ
                           for j in range(16)] for i in range(16)])


if __name__ == '__main__':
    unittest.main()
