"""Exact rational checks of the modular occupancy-clock counterexample.

These test the supplied-unitary gadget in the source comparison, not a
Clifford+T emitter. The analytical counterexample is in SOURCE_CATALOGUE.md.
"""
from fractions import Fraction as F
import unittest


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, column))
             for column in zip(*b)] for row in a]


def rotation(k, inverse=False):
    # theta = 2 atan(1/(2k)); full-unitary distance squared is 4/(4k^2+1).
    c, s = F(4*k*k-1, 4*k*k+1), F(4*k, 4*k*k+1)
    if inverse:
        s = -s
    return [[c, -s], [s, c]]


def flag_gate(r, period):
    return [[r[i//period][j//period] * (i % period == j % period)
             for j in range(2*period)] for i in range(2*period)]


def clock_gate(period):
    result = [[F(0) for _ in range(2*period)] for _ in range(2*period)]
    for flag in range(2):
        for clock in range(period):
            result[flag*period+(clock+flag) % period][flag*period+clock] = F(1)
    return result


def gadget(k, period, pulses):
    """U_k D U_(k-1) ... D U_1, retaining every rejected clock mode."""
    result, tick = identity(2*period), clock_gate(period)
    for stage in range(1, k+1):
        if stage in pulses:
            result = multiply(flag_gate(pulses[stage], period), result)
        if stage != k:
            result = multiply(tick, result)
    return result


class ModularHistoryAuditTests(unittest.TestCase):
    def test_two_pulses_return_after_four_ticks_but_not_without_wrap(self):
        for k in (5, 9, 17):
            r = rotation(k)
            self.assertEqual(multiply(rotation(k, True), r), identity(2))
            self.assertLessEqual(F(4, 4*k*k+1), F(1, k*k))
            pulses = {1: r, 5: rotation(k, True)}
            actual = gadget(k, 4, pulses)
            self.assertEqual([row[0] for row in actual], [F(1)]+[F(0)]*7)
            intended = r[0][0]**2
            defect = F(1)-intended
            self.assertEqual(defect, F(16*k*k, (4*k*k+1)**2))
            # e<3 gives a rational upper bound on the claimed e^5/(128 k^4).
            self.assertGreater(defect, F(3**5, 128*k**4))
            # More clock states than all possible ticks: genuinely nonwrapping.
            no_wrap = gadget(k, 1 << k.bit_length(), pulses)
            self.assertEqual(no_wrap[0][0], intended)

    def test_repeated_pulses_reduce_to_one_accumulated_rotation(self):
        period, repetitions = 4, 4
        k = repetitions*period+1
        r = rotation(k)
        pulses = {1+j*period: r for j in range(repetitions+1)}
        actual = gadget(k, period, pulses)
        accumulated = identity(2)
        for _ in pulses:
            accumulated = multiply(r, accumulated)
        self.assertEqual(actual, flag_gate(accumulated, period))
        self.assertNotEqual(actual[0][0], r[0][0]**len(pulses))


if __name__ == "__main__":
    unittest.main()
