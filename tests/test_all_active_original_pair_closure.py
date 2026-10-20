"""Two actual all-active pair mechanisms, not Boolean-target enumeration.

The four-wire source has the clean Pauli in the total physical axis span.
An arbitrary fifth logical wire puts it outside that span and supplies
an actual clean flip. Both sources keep the full unused clean sector;
neither is asserted to be a Boolean target.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns
from tests.test_clean_code_clifford_rigidity import _rank
from tests.test_original_middle_pair_insertions import _ghz, _pair_word
from tests.test_three_layer_native_mechanisms import _native_t_bands
from tests.test_two_layer_diagonal_rigidity import _inverse
from tests.test_two_layer_galois_return import _symplectic


class AllActiveOriginalPairClosureTests(unittest.TestCase):
    def test_original_span_dichotomy_and_actual_shared_flip(self):
        core = 4
        p = (1 << core) - 1
        for tagged in (False, True):
            with self.subTest(arbitrary_logical_tag=tagged):
                width = core + int(tagged)
                size = 1 << width
                active = [1 << (width + bit) for bit in range(core)]
                r = p << width
                b = p | ((1 << (width + core)) if tagged else 0)
                d = b ^ r
                axes = [active, [b, d], active]
                for axis_list in axes:
                    self.assertEqual(_rank(axis_list), len(axis_list))
                    self.assertTrue(all(not _symplectic(left, right, width)
                                        for left in axis_list for right in axis_list))
                self.assertEqual(b ^ d, sum(active))
                self.assertTrue(all(_symplectic(a, p, width) for a in active))
                self.assertTrue(all(_symplectic(a, b, width) for a in active))
                entire_span = active + [b, d]
                self.assertEqual(_rank(entire_span + [p]) - _rank(entire_span),
                                 int(tagged))

                # H0 is isotropic. In the untagged case its final
                # generator is the original clean Pauli; in the tagged
                # case the shared-flip branch makes no such substitution.
                r_basis = [active[0] ^ a for a in active[1:]]
                h = b ^ (active[0] if _symplectic(p, b, width) else 0)
                h0 = r_basis + [h]
                self.assertEqual(_rank(h0) - _rank(r_basis), 1)
                self.assertTrue(all(not _symplectic(left, right, width)
                                    for left in h0 for right in h0))
                self.assertTrue(all(not _symplectic(p, label, width) for label in h0))

                # CZ(0,tag) appends Z_tag to both actual pair axes.
                tag_word = ([('H', core), ('CX', 0, core), ('H', core)]
                            if tagged else [])
                pair = tag_word + _pair_word(core, 1, 1) + tag_word
                endpoint = _ghz(core)
                band = [('T', bit) for bit in range(core)]
                word = endpoint + band + pair + _inverse(band) + _inverse(endpoint)
                actual = native_columns(width, word)
                self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
                for row in range(size):
                    for column in range(size):
                        if (row ^ column) & 1:
                            self.assertEqual(actual[row][column], CZ)
                for candidate in (word, _inverse(word)):
                    counts = Counter(gate[0] for gate in candidate)
                    self.assertEqual(counts['T'] + counts['TDG'], 2 * core + 2)
                    self.assertEqual(counts['H'] + counts['CX'],
                                     6 * core - 2 + 6 * int(tagged))
                    self.assertEqual([len(band) for band in _native_t_bands(candidate)],
                                     [core, 2, core])

                if tagged:
                    # The actual interior flip Z0 X_tag commutes all
                    # physical phase axes and flips both endpoint Paulis.
                    flip = active[0] | (1 << core)
                    self.assertEqual(_symplectic(flip, p, width), 1)
                    self.assertTrue(all(not _symplectic(flip, axis, width)
                                        for axis_list in axes for axis in axis_list))
                    # Its endpoint preimage is X0 X_tag. The dressing
                    # fixes every clean-zero column, with no extra wire.
                    dress = [('CX', 0, core)]
                    dressed = native_columns(width, dress + word + dress)
                    for row in range(size):
                        for column in range(size):
                            self.assertEqual(dressed[row][column],
                                             dressed[row ^ 1][column ^ 1])
                            if not (column & 1):
                                self.assertEqual(dressed[row][column], actual[row][column])
                    self.assertEqual(native_columns(width, _inverse(dress + word + dress)),
                                     adjoint(dressed))
                else:
                    # H0's whole logical image is computational Z here.
                    self.assertEqual(h, p)
                    self.assertTrue(all(actual[row][column] == CZ
                                        for row in range(size) for column in range(size)
                                        if row != column))
                    omega = native_columns(1, [('T', 0)])[1][1]
                    self.assertEqual(actual[2][2] / actual[6][6], omega)
                    self.assertNotIn(omega, (CO, -CO))


if __name__ == '__main__':
    unittest.main()
