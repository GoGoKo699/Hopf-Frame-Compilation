"""Exact bounded audits of the complete-bank Fourier source invariant.

The proof is in READ_ONLY_TWO_BAND_SOURCE_GAP.md.  These fixtures use
Pauli labels and eighth-root phases, not a source search or a native
construction.  No dense exponentially wide bank matrix is formed.
"""
from fractions import Fraction
from itertools import product
import unittest


def _mul(left, right):
    """Multiply omega**phase X**x Z**z, retaining its literal phase."""
    lp, lx, lz = left
    rp, rx, rz = right
    phase = lp + rp + 4 * ((lz & rx).bit_count() % 2)
    return phase % 8, lx ^ rx, lz ^ rz


def _adjoint(pauli):
    phase, x, z = pauli
    return (-phase + 4 * ((x & z).bit_count() % 2)) % 8, x, z


def _add_trace(total, pauli, sign=1):
    """Accumulate the normalized trace exactly in Z[omega]."""
    phase, x, z = pauli
    if x or z:
        return
    total[phase % 4] += sign if phase < 4 else -sign


def _character(character, address):
    return -1 if (character & address).bit_count() % 2 else 1


def _moment_numerator(targets, character, order):
    """Normalized trace of (D*D)^(order/2), before division by M^order."""
    total = [0, 0, 0, 0]
    for indices in product(range(len(targets)), repeat=order):
        word = (0, 0, 0)
        address_xor = 0
        for index, address in enumerate(indices):
            factor = targets[address]
            if index % 2 == 0:
                factor = _adjoint(factor)
            word = _mul(word, factor)
            address_xor ^= address
        _add_trace(total, word, _character(character, address_xor))
    return total


class ReadOnlyTwoBandSourceGapTests(unittest.TestCase):
    def test_projective_cocycle_cancels_with_actual_conjugate_phase(self):
        # S_h has linear X/Z labels but a nontrivial scalar cocycle.
        # A_h = H S_h H.  Gamma_h = A_h tensor conjugate(S_h).
        source = []
        gamma = []
        for h in range(4):
            x, z = h & 1, (h >> 1) & 1
            phase = (2 * x + 4 * x * z) % 8
            source.append((phase, x, z))
            a_phase = (phase + 4 * x * z) % 8
            gamma.append(((a_phase - phase) % 8,
                          z | (x << 1), x | (z << 1)))
        self.assertNotEqual(_mul(source[1], source[2]),
                            _mul(source[2], source[1]))
        for h, g in product(range(4), repeat=2):
            self.assertEqual(_mul(gamma[h], gamma[g]), gamma[h ^ g])
        for word in gamma:
            self.assertEqual(_adjoint(word), word)
            self.assertEqual(_mul(word, word), (0, 0, 0))

    def test_aggregate_fourth_moment_keeps_affine_labels_and_all_phases(self):
        for rank in (2, 3):
            size = 1 << rank
            targets = []
            for h in range(size):
                # Distinct complete-bank X labels plus a linear Pauli map.
                # Both noncommuting Z terms and genuine eighth-root phases
                # are retained.  Private labels may overlap output labels.
                x, z = 1 << h, 0
                for bit in range(rank):
                    if h & (1 << bit):
                        x ^= (1 << bit) | (1 << (size + bit))
                        z ^= (1 << ((bit + 1) % size)) | (1 << bit)
                phase = (h + 3 * ((h & 1) * ((h >> 1) & 1))) % 8
                targets.append((phase, x, z))

            second = [0, 0, 0, 0]
            fourth = [0, 0, 0, 0]
            for character in range(size):
                for total, order in ((second, 2), (fourth, 4)):
                    term = _moment_numerator(targets, character, order)
                    for coordinate in range(4):
                        total[coordinate] += term[coordinate]
            self.assertEqual(second, [size ** 2, 0, 0, 0])

            constrained = [0, 0, 0, 0]
            nonzero_traces = 0
            for a, b, c in product(range(size), repeat=3):
                d = a ^ b ^ c
                word = _mul(_mul(_mul(_adjoint(targets[a]), targets[b]),
                                 _adjoint(targets[c])), targets[d])
                # Character orthogonality removes the linear perturbation;
                # the remaining independent output labels must pair.
                paired = all((a, b, c, d).count(value) % 2 == 0
                             for value in {a, b, c, d})
                self.assertEqual(word[1:] == (0, 0), paired)
                nonzero_traces += word[1:] == (0, 0)
                _add_trace(constrained, word)
            self.assertEqual(nonzero_traces, 3 * size ** 2 - 2 * size)
            self.assertEqual(fourth, [size * value for value in constrained])
            self.assertEqual(fourth[1:], [0, 0, 0])
            actual = Fraction(fourth[0], size ** 4)
            self.assertGreaterEqual(actual, 0)
            self.assertLessEqual(actual, Fraction(3 * size - 2, size ** 2))


if __name__ == '__main__':
    unittest.main()
