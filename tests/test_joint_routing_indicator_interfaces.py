"""Literal native checks for the selected-swap/indicator interface.

The finite fixtures certify the exposed words and Clifford ledgers, not
a new sublogarithmic compiler or an unrestricted lower bound.
"""
from itertools import combinations
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_majorana_walsh_query import _inverse


def _x(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


def _phase_parity(controls, target, dagger=False):
    compute = [('CX', control, target) for control in controls]
    return compute + [('TDG' if dagger else 'T', target)] + list(reversed(compute))


def _ccz_batch(control, pairs):
    count = len(pairs)
    word = [('S', control)] * ((count // 2) % 4)
    if count % 2:
        word.append(('T', control))
    word += [('T', wire) for pair in pairs for wire in pair]
    for left, right in pairs:
        word += _phase_parity((control,), left, True)
        word += _phase_parity((control,), right, True)
    for left, right in pairs:
        word += _phase_parity((left,), right, True)
    for left, right in pairs:
        word += _phase_parity((control, left), right)
    return word


def _shears(t, a, y):
    middle = _ccz_batch(t, list(zip(a, y)))
    hg = [('H', t)]
    hs = [('H', wire) for wire in y]
    return hg + middle + hg, hs + middle + hs


def _counts(word):
    count = sum(gate[0] in ('T', 'TDG') for gate in word)
    return count, len(word) - count


def _assert_permutation(test, width, word, target):
    matrix = native_columns(width, word)
    for column in range(1 << width):
        expected = target(column)
        for row in range(1 << width):
            test.assertEqual(matrix[row][column], CO if row == expected else CZ)
    inverse = native_columns(width, _inverse(word))
    test.assertEqual(inverse, adjoint(matrix))


def _pruned_router(p):
    m = 1 << p
    word = []
    for level in range(p):
        stride = 1 << level
        pairs = [(p + base, p + base + stride)
                 for base in range(0, m, 2 * stride)]
        outside = [('CX', left, right) for left, right in pairs]
        h = [('H', left) for left, _ in pairs]
        word += outside + h + _ccz_batch(level, pairs) + h + outside
    return word


class JointRoutingIndicatorInterfacesTests(unittest.TestCase):
    def test_shared_control_shears_and_six_call_swap_on_every_column(self):
        # All x,t,Y,A inputs are arbitrary. The one-bit indicator source
        # is Clifford; using it does not remove the bilinear native words.
        x, t, y, a = 0, 1, (2, 3), (4, 5)
        indicator = _x(a[0]) + [('CX', x, a[0]), ('CX', x, a[1])]
        bg, bs = _shears(t, a, y)
        for word, bound in ((bg, 25), (bs, 27)):
            for orientation in (word, _inverse(word)):
                self.assertEqual(_counts(orientation)[0], 12)
                self.assertLessEqual(_counts(orientation)[1], bound)

        def echo(b):
            return b + indicator + _inverse(b) + _inverse(indicator)

        gather, scatter = echo(bg), echo(bs)
        swap = gather + scatter + gather

        def target(column):
            selected = y[column & 1]
            difference = ((column >> t) ^ (column >> selected)) & 1
            return column ^ (difference << t) ^ (difference << selected)

        _assert_permutation(self, 6, swap, target)
        common_i = max(_counts(indicator)[1], _counts(_inverse(indicator))[1])
        for orientation in (swap, _inverse(swap)):
            self.assertEqual(_counts(orientation)[0], 72)
            self.assertLessEqual(_counts(orientation)[1], 6 * common_i + 154)

        extraction = swap + _x(t) + _inverse(swap)
        _assert_permutation(self, 6, extraction,
                            lambda column: column ^ (1 << y[column & 1]))
        # The correct source adjoint is indispensable: skipping the last
        # swap leaves t changed on an allowed input.
        wrong = native_columns(6, swap + _x(t), columns=(0,))
        self.assertNotEqual(wrong[1 << y[0]][0], CO)

    def test_pruned_native_router_and_actual_inverse(self):
        for p in (1, 2):
            m, router = 1 << p, _pruned_router(p)
            indicator = router + _x(p) + _inverse(router)
            _assert_permutation(self, p + m, indicator,
                                lambda column: column ^ (1 << (p + (column & (m - 1)))))
            for orientation in (indicator, _inverse(indicator)):
                t, g = _counts(orientation)
                self.assertEqual(t, 12 * (m - 1) + 2)
                self.assertLessEqual(g, 28 * (m - 1) + 6 * p + 4)

    def test_translation_difference_products_do_not_flatten(self):
        for p in range(2, 7):
            m = 1 << p
            for length in range(1, p + 1):
                for directions in combinations(range(p), length):
                    vector = 1
                    for direction in directions:
                        shifted = sum(((vector >> label) & 1) << (label ^ (1 << direction))
                                      for label in range(m))
                        vector ^= shifted
                    self.assertEqual(vector.bit_count(), 1 << length)
                    self.assertNotEqual(vector, 0)


if __name__ == '__main__':
    unittest.main()
