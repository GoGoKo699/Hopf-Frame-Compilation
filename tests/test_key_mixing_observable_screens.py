"""Two exact single-column falsifiers for the paid key-mixing words.

Only r=3, x=e_1 and the two-bank input |0,0> are evaluated. Integer
Walsh transforms retain exact amplitudes; this is not an all-column
simulation or a native-source emitter. The full identities and ledgers
precede these fixtures in KEY_MIXING_OBSERVABLE_SCREENS.md.
"""
from fractions import Fraction
import unittest

from tests.test_balanced_algebra_permutation import _hrt_k


class KeyMixingObservableScreenTests(unittest.TestCase):
    def test_two_exact_single_column_counterexamples(self):
        r, w, d = 3, 7, 128
        top = 1 << (w - 1)
        k = [_hrt_k(a << 1, r) >> 1 for a in range(d)]

        def query(state, reverse_direction=False):
            result = [0] * (d * d)
            for a in range(d):
                for beta in range(d):
                    image = (((a ^ k[beta]) * d + beta)
                             if reverse_direction
                             else (a * d + (beta ^ k[a])))
                    result[image] = state[a * d + beta]
            return result

        def walsh_a(state):
            result = state.copy()
            half = 1
            while half < d:
                for base in range(0, d, 2 * half):
                    for offset in range(half):
                        lower = (base + offset) * d
                        upper = (base + offset + half) * d
                        for beta in range(d):
                            lo, hi = result[lower + beta], result[upper + beta]
                            result[lower + beta] = lo + hi
                            result[upper + beta] = lo - hi
                half *= 2
            return result

        def apply(word):
            # State amplitudes are numerator / denominator. M=H_A Q H_A
            # contributes two integer Walsh transforms and denominator d.
            numerator = [1] + [0] * (d * d - 1)
            denominator = 1
            for letter in word:
                if letter == 'Q':
                    numerator = query(numerator)
                elif letter == 'R':
                    numerator = query(numerator, reverse_direction=True)
                elif letter == 'M':
                    numerator = walsh_a(query(walsh_a(numerator)))
                    denominator *= d
                elif letter == 'X':
                    # u=iota(e_1), including every arbitrary carrier slot.
                    numerator = [numerator[(a ^ 1) * d + beta]
                                 for a in range(d) for beta in range(d)]
                else:
                    self.fail(f'Unknown exact letter {letter}')
            norm = sum(value * value for value in numerator)
            self.assertEqual(norm, denominator * denominator)
            expectations = []
            for bank in ('A', 'B'):
                signed = sum(value * value *
                             (-1 if ((index // d if bank == 'A' else index % d)
                                     & top) else 1)
                             for index, value in enumerate(numerator))
                expectations.append(Fraction(signed, norm))
            return tuple(expectations)

        # All words below are chronological. Each letter is self-adjoint;
        # S4=Q R M Q in operator order, J4=S4 X S4†.
        source = ('Q', 'M', 'R', 'Q')
        j4 = tuple(reversed(source)) + ('X',) + source
        self.assertEqual(apply(j4), (Fraction(0), Fraction(5, 16)))

        # L_A=M Q X Q M; J=L_A R L_A R in operator order.
        la = ('M', 'Q', 'X', 'Q', 'M')
        two_insertions = ('R',) + la + ('R',) + la
        self.assertEqual(apply(two_insertions), (Fraction(0), Fraction(0)))
        # Both target signs at a=beta=0 are +1. Thus these expectations
        # imply the stated basis-vector defects without testing other columns.


if __name__ == '__main__':
    unittest.main()
