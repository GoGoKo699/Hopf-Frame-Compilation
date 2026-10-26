"""Exact full-space guards for actual-signal compression and one biased pilot.

The width-independent error bounds are analytic. These fixtures retain
literal source phases, arbitrary dirty inputs and rejected amplitudes.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_clifford_involution_source import _phase_core, _word
from tests.test_co_moving_lookup_stage import _ccz_7t
from tests.test_joint_lookup_pauli_source import _bands, _inverse
from tests.test_two_path_precision_bypass import _minus_identity


def _sub(left, right):
    return add(left, scale(right, -1))


def _product(*matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = multiply(result, matrix)
    return result


def _block(matrix, rows, columns):
    return [[matrix[row][column] for column in columns] for row in rows]


def _recover(source, signal):
    width = (len(source) - 1).bit_length()
    z = native_columns(width, [('S', signal)] * 2)
    return scale(_product(source, z, adjoint(source), z, source), -1)


def _bias_seed(signal):
    return [('H', signal), ('T', signal), ('H', signal),
            ('T', signal), ('H', signal)]


def _accepted(source, signal):
    code = [basis for basis in range(len(source)) if not basis >> signal & 1]
    return _block(source, code, code), code


class PrecisionSignalAndBiasTests(unittest.TestCase):
    def test_actual_signal_compression_with_mixed_clifford_walls(self):
        core, width, _ = _phase_core(3, (), 1)
        before = [('H', 2), ('CX', 0, 2), ('S', 1), ('H', 0)]
        after = [('H', 0), ('CX', 2, 1), ('H', 2), ('S', 0)]
        word = before + _word(core) + after
        source = native_columns(width, word)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(source))
        for signal in range(width):
            k, code = _accepted(source, signal)
            identity = eye(len(code))
            g = multiply(adjoint(k), k)
            m = _sub(scale(g, 2), identity)
            z = native_columns(width, [('S', signal)] * 2)
            actual_image = _product(adjoint(source), z, source)
            self.assertEqual(m, _block(actual_image, code, code))
            self.assertEqual(m, adjoint(m))
            m_squared = multiply(m, m)
            # Roots are precisely 0, +/-1 and +/-1/sqrt(2).
            self.assertEqual(_product(m, _sub(m_squared, identity),
                                      _sub(scale(m_squared, 2), identity)),
                             zeros(len(code)))
            recovered = _recover(source, signal)
            b = multiply(k, _sub(scale(identity, 3), scale(g, 4)))
            self.assertEqual(_block(recovered, code, code), b)
            rejected = [basis for basis in range(1 << width) if basis not in code]
            leakage = _block(recovered, rejected, code)
            self.assertEqual(multiply(adjoint(leakage), leakage),
                             _sub(identity, multiply(adjoint(b), b)))

        # A source which already preserves the code has exactly FJ=-QJ.
        source = native_columns(2, [('H', 0), ('T', 0), ('H', 0)])
        recovered = _recover(source, 1)
        self.assertEqual([row[:2] for row in recovered],
                         scale([row[:2] for row in source], -1))

    def test_paid_unequal_seed_has_a_nonclifford_actual_signal_image(self):
        seed = native_columns(1, _bias_seed(0))
        self.assertEqual(seed[0][0].conj() * seed[0][0], CO * Fraction(3, 4))
        self.assertEqual(seed[1][0].conj() * seed[1][0], CO * Fraction(1, 4))
        order_three = [('CX', 0, 1), ('CX', 1, 0)]
        controlled = []
        for _, control, target in order_three:
            controlled += ([('H', target)] + _ccz_7t(2, control, target)
                           + [('H', target)])
        word = _bias_seed(2) + controlled + _inverse(_bias_seed(2))
        recovered_word = (word + [('S', 2)] * 2 + _inverse(word)
                          + [('S', 2)] * 2 + word + _minus_identity(0))
        source = native_columns(3, word)
        r, identity = native_columns(2, order_three), eye(4)
        self.assertEqual(_product(r, r, r), identity)
        k, code = _accepted(source, 2)
        self.assertEqual(k, scale(add(scale(identity, 3), r), Fraction(1, 4)))
        p = _sub(identity, scale(add(add(identity, r), multiply(r, r)),
                                 Fraction(1, 3)))
        self.assertEqual(multiply(p, p), p)
        self.assertNotEqual(p, zeros(4))
        g = multiply(adjoint(k), k)
        self.assertEqual(g, _sub(identity, scale(p, Fraction(9, 16))))
        m = _sub(scale(g, 2), identity)
        # A nonzero -1/8 eigenspace is impossible for a compressed Clifford.
        self.assertEqual(multiply(m, p), scale(p, Fraction(-1, 8)))
        recovered = native_columns(3, recovered_word)
        self.assertEqual(recovered, _recover(source, 2))
        b = _block(recovered, code, code)
        self.assertEqual(multiply(multiply(adjoint(b), b), p),
                         scale(p, Fraction(175, 256)))
        for native in (word, recovered_word):
            self.assertEqual(native_columns(3, native + _inverse(native)), eye(8))
        count = lambda native: sum(gate[0] in ('T', 'TDG') for gate in native)
        self.assertEqual((count(word), count(recovered_word)), (18, 54))
        for native, cliffords, bands in ((word, 22, 14), (recovered_word, 82, 42)):
            for oriented in (native, _inverse(native)):
                self.assertLessEqual(len(oriented) - count(oriented), cliffords)
                self.assertLessEqual(len(_bands(oriented)), bands)
                self.assertTrue(all(len({gate[1] for gate in band}) == len(band)
                                    for band in _bands(oriented)))

    def test_biased_recovery_and_rejected_gram_identities_without_commutation(self):
        # The generic branches do not commute, and the formula retains them.
        u = native_columns(2, [('H', 0), ('T', 0), ('CX', 0, 1), ('S', 1)])
        v = native_columns(2, [('H', 1), ('TDG', 0), ('CX', 1, 0), ('T', 1)])
        self.assertNotEqual(multiply(u, v), multiply(v, u))
        identity = eye(4)
        r = multiply(adjoint(u), v)
        t = scale(add(r, adjoint(r)), Fraction(1, 2))
        k = scale(add(scale(u, 3), v), Fraction(1, 4))
        g = multiply(adjoint(k), k)
        self.assertEqual(g, scale(add(scale(identity, 5), scale(t, 3)), Fraction(1, 8)))
        b = multiply(k, _sub(scale(identity, 3), scale(g, 4)))
        polynomial = add(_sub(scale(identity, 3), scale(r, 7)),
                         add(scale(adjoint(r), -9), scale(multiply(r, r), -3)))
        self.assertEqual(b, scale(multiply(u, polynomial), Fraction(1, 16)))
        leakage = _sub(identity, multiply(adjoint(b), b))
        self.assertEqual(leakage,
                         scale(_product(_sub(identity, t), add(identity, t),
                                        add(identity, t)), Fraction(27, 32)))
        d = add(b, v)
        defect = multiply(adjoint(d), d)
        self.assertEqual(defect,
                         scale(_product(_sub(scale(identity, 5), scale(t, 3)),
                                        _sub(identity, t), add(identity, t)),
                               Fraction(9, 32)))
        # Exact factorization behind defect^2 <= 12*leakage on -I <= t <= I.
        positive_quadratic = add(_sub(scale(multiply(t, t), 9), scale(t, 48)),
                                 scale(identity, 103))
        self.assertEqual(_sub(scale(leakage, 12), multiply(defect, defect)),
                         scale(_product(leakage, add(identity, t), positive_quadratic),
                               Fraction(3, 32)))

        # Noncommuting branch/projector exact recovery is just the -V branch.
        r = native_columns(2, [('S', 1), ('S', 1)])
        v = multiply(u, r)
        k = scale(add(scale(u, 3), v), Fraction(1, 4))
        g = multiply(adjoint(k), k)
        b = multiply(k, _sub(scale(identity, 3), scale(g, 4)))
        p = scale(_sub(identity, r), Fraction(1, 2))
        self.assertNotEqual(multiply(v, p), multiply(p, v))
        self.assertEqual(b, scale(v, -1))
        self.assertEqual(k, multiply(scale(v, -1),
                                    _sub(scale(p, Fraction(3, 2)), identity)))


if __name__ == '__main__':
    unittest.main()
