"""Exact bounded evidence for one-relation projected native bands.

The phase source is evaluated on every basis input and differentiated
before its Pauli coefficients are recovered by an exact Walsh transform.
No predicted support or norm formula is used to generate those values.
The analytic matching theorem, rather than these finite cases, supplies
the arbitrary-width and robust conclusions.
"""
from fractions import Fraction
import unittest

from tests.exact_matrices import CO, CQ2, Q2, CZ, adjoint, native_columns


_ROOTS = (
    (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1),
    (-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1),
)


def _walsh(values):
    """Unnormalized exact transform, for integer four-tuples or integers."""
    result = list(values)
    step = 1
    while step < len(result):
        for start in range(0, len(result), 2 * step):
            for offset in range(step):
                i, j = start + offset, start + offset + step
                a, b = result[i], result[j]
                if isinstance(a, tuple):
                    result[i] = tuple(x + y for x, y in zip(a, b))
                    result[j] = tuple(x - y for x, y in zip(a, b))
                else:
                    result[i], result[j] = a + b, a - b
        step *= 2
    return result


def _norm_squared_numerator(value):
    """Exact |a+b*w+c*w^2+d*w^3|^2 as A+B*sqrt(2)."""
    a, b, c, d = value
    return (a*a + b*b + c*c + d*d,
            a*(b-d) + c*(b+d))


def _phase_exponents(width, axes, orientations=None, offsets=0):
    if orientations is None:
        orientations = [1] * len(axes)
    return [sum(sign * (((axis & x).bit_count() & 1)
                        ^ ((offsets >> i) & 1))
                for i, (axis, sign) in enumerate(zip(axes, orientations))) % 8
            for x in range(1 << width)]


def _observed_profile(width, axes, orientations=None, offsets=0):
    """Recover log2 Pauli support from the actual diagonal operator."""
    phase = _phase_exponents(width, axes, orientations, offsets)
    size = 1 << width
    profile = []
    for translation in range(size):
        derivative = [_ROOTS[(phase[x ^ translation] - phase[x]) % 8]
                      for x in range(size)]
        coefficients = _walsh(derivative)
        nonzero = [value for value in coefficients if value != (0, 0, 0, 0)]
        support = len(nonzero)
        if support == 0 or support & (support - 1):
            raise AssertionError('Pauli support is not a nonzero power of two')
        for value in nonzero:
            rational, irrational = _norm_squared_numerator(value)
            if irrational != 0 or rational * support != size * size:
                raise AssertionError('Exact Pauli coefficients are not flat')
        profile.append(support.bit_length() - 1)
    return profile


def _spectrum(profile):
    return {label: Fraction(value, len(profile))
            for label, value in enumerate(_walsh(profile)) if label}


def _cycle_axes(length):
    return [1 << j for j in range(length - 1)] + [(1 << (length - 1)) - 1]


def _cq_root(exponent):
    a, b, c, d = _ROOTS[exponent % 8]
    return CQ2(Q2(a, Fraction(b-d, 2)), Q2(c, Fraction(b+d, 2)))


class ProjectedBandMatchingTests(unittest.TestCase):
    def test_exact_native_lift_retains_all_columns_and_inverse(self):
        # One physical clean wire supplies the unique projected relation.
        # Every physical column is checked, including the unused clean-one code.
        width, clean = 4, 3
        compute = [('CX', j, clean) for j in range(clean)]
        band = [('T', j) for j in range(width)]
        word = compute + band + list(reversed(compute))
        actual = native_columns(width, word)
        inverse = compute + [('TDG', j) for j in range(width)] + list(reversed(compute))
        expected = [[CZ for _ in range(1 << width)] for _ in range(1 << width)]
        for x in range(1 << width):
            low = x & 7
            exponent = low.bit_count() + (((x >> clean) & 1) ^ (low.bit_count() & 1))
            expected[x][x] = _cq_root(exponent)
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(width, inverse), adjoint(actual))
        self.assertEqual(len(band), 4)
        self.assertEqual(len(compute) * 2, 6)
        # The clean-zero block is the four-parity source used below.
        phase = _phase_exponents(3, _cycle_axes(4))
        for x, exponent in enumerate(phase):
            self.assertEqual(actual[x][x], _cq_root(exponent))

    def test_exact_pauli_images_for_relation_lengths_two_through_eight(self):
        translations = 0
        for length in range(2, 9):
            width = length - 1
            axes = _cycle_axes(length)
            reference = None
            for orientations, offsets in (
                    ([1] * length, 0),
                    ([(-1) ** j for j in range(length)], 0),
                    ([1] * length, 5),
                    ([(-1) ** (j + 1) for j in range(length)], 3)):
                profile = _observed_profile(width, axes, orientations, offsets)
                if reference is None:
                    reference = profile
                    translations += len(profile)
                else:
                    self.assertEqual(profile, reference)
            for t, observed in enumerate(reference):
                active = sum((axis & t).bit_count() & 1 for axis in axes)
                all_active = active == length
                self.assertEqual(observed, active - (2 if all_active else 0))
                if length % 2:
                    self.assertFalse(all_active)
        self.assertEqual(translations, 254)

    def test_fourier_recovery_and_the_rank_three_exception(self):
        for length in range(3, 9):
            axes = _cycle_axes(length)
            spectrum = _spectrum(_observed_profile(length - 1, axes))
            if length == 4:
                self.assertEqual(spectrum, {v: Fraction(-1, 4) for v in range(1, 8)})
                self.assertNotEqual(set(spectrum), set(axes))
            else:
                recovered = {v for v, coefficient in spectrum.items()
                             if abs(coefficient) > Fraction(1, 4)}
                self.assertEqual(recovered, set(axes))
                if length % 2:
                    self.assertTrue(all(spectrum[v] == Fraction(-1, 2) for v in axes))
                else:
                    delta = Fraction(1, 1 << (length - 2))
                    self.assertTrue(all(spectrum[v] == Fraction(-1, 2) + delta
                                        for v in axes))
                    self.assertTrue(all(abs(value) == delta for v, value in spectrum.items()
                                        if v not in axes))

    def test_external_axes_are_separate_and_four_relation_does_not_recover_axes(self):
        # The external wire is independent modulo the unique relation span.
        for length in (4, 6):
            axes = _cycle_axes(length)
            outside = 1 << (length - 1)
            spectrum = _spectrum(_observed_profile(length, axes + [outside]))
            self.assertEqual(spectrum[outside], Fraction(-1, 2))
            self.assertTrue(all(value == 0 for v, value in spectrum.items()
                                if v & outside and v != outside))
        # A native binary coordinate change preserves the norm profile while
        # changing which four axes represent the exceptional rank-three span.
        original = [1, 2, 4, 7]
        changed = [1, 3, 4, 6]
        self.assertNotEqual(set(original), set(changed))
        self.assertEqual(_observed_profile(3, original), _observed_profile(3, changed))
        # A length-two repeated axis is an actual Clifford, not a new active
        # axis detectable by this non-Clifford coefficient norm invariant.
        self.assertEqual(_observed_profile(1, [1, 1]), [0, 0])
        s = native_columns(1, [('S', 0)])
        t_squared = native_columns(1, [('T', 0), ('T', 0)])
        self.assertEqual(t_squared, s)
        self.assertEqual(native_columns(1, [('T', 0), ('TDG', 0)]), [[CO, CZ], [CZ, CO]])


if __name__ == '__main__':
    unittest.main()
