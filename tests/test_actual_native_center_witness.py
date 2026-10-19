"""Actual one-clean no-pair source with an easy CCZ target and surviving center.

The checks retain all physical columns, literal phases, independent native
lists, and the separately reversed word. They test a limitation of two
derivative-cancellation criteria, not a lower bound on the target.
"""
from collections import Counter
import unittest

from tests.exact_matrices import (
    CZ, CO, CI, adjoint, eye, multiply, native_columns, scale,
)
from tests.test_two_layer_galois_return import _conjugate, _inverse, _tau


def _rank(labels):
    basis = {}
    for label in labels:
        while label:
            pivot = label.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = label
                break
            label ^= basis[pivot]
    return len(basis)


def _support(operator, width):
    """Exact Pauli coefficients via the sparse tensor-Pauli trace."""
    size = 1 << width
    support = set()
    for label in range(size * size):
        x, z = label & (size - 1), label >> width
        phase = (CO, CI, -CO, -CI)[(x & z).bit_count() % 4]
        coefficient = sum(
            (phase.conj() * (-1 if (z & column).bit_count() % 2 else 1)
             * operator[column ^ x][column] for column in range(size)), CZ
        ) / size
        if coefficient != CZ:
            support.add(label)
    return support


def _source():
    # Physical order u,v,r,s; u is the sole initialized input wire.
    u, v, r, s = range(4)
    ca = [('H', v), ('CX', v, s), ('CX', r, s), ('CX', u, r)]
    cc = [('H', u), ('CX', u, r), ('CX', s, r), ('CX', v, s)]
    cb = [('CX', 3, 1), ('CX', 3, 0), ('CX', 2, 3), ('CX', 2, 0),
          ('CX', 1, 3), ('CX', 1, 2), ('CX', 2, 1), ('CX', 0, 1)]
    phases = ([('T', u), ('T', r), ('T', v), ('TDG', s)],
              [('T', q) for q in range(4)],
              [('T', v), ('T', s), ('T', u), ('TDG', r)])
    changes = (ca, cb, cc)
    bands = [change + phase + _inverse(change)
             for change, phase in zip(changes, phases)]
    c0 = [('S', u)] * 2 + [('S', r)] * 2 + [('S', s)]
    for control, target in ((u, r), (u, s), (v, r), (r, s)):
        c0 += [('H', target), ('CX', control, target), ('H', target)]
    ka = [('H', u), ('CX', u, v)]
    kc = ka + c0
    word = ka + bands[0] + bands[1] + _inverse(bands[2]) + _inverse(kc)
    return changes, phases, bands, ka, kc, c0, word


class ActualNativeCenterWitnessTests(unittest.TestCase):
    def test_full_native_source_inverse_no_flip_and_ledger(self):
        changes, phases, bands, ka_word, kc_word, c0_word, word = _source()
        ea, eb, ec = [native_columns(4, band) for band in bands]
        ka, kc, c0 = [native_columns(4, w) for w in (ka_word, kc_word, c0_word)]
        labels = ((16, 80, 2, 194), (208, 144, 224, 96),
                  (32, 160, 1, 193))
        # Each actual band is Clifford-conjugated simultaneous T/T-dagger.
        for change, phase, axes in zip(changes, phases, labels):
            coordinate = native_columns(4, change)
            self.assertEqual(_rank(axes), 4)
            for gate, axis in zip(phase, axes):
                self.assertEqual(_conjugate(adjoint(coordinate),
                                            _tau(1 << (4 + gate[1]), 4)),
                                 _tau(axis, 4))
            for left in axes:
                for right in axes:
                    self.assertEqual(multiply(_tau(left, 4), _tau(right, 4)),
                                     multiply(_tau(right, 4), _tau(left, 4)))
        pa = _tau(3, 4)  # X_u X_v.
        pc = scale(_tau(3 + (8 << 4), 4), -1)
        self.assertEqual(_conjugate(ka, _tau(16, 4)), pa)
        self.assertEqual(_conjugate(kc, _tau(16, 4)), pc)
        self.assertEqual(_conjugate(c0, pa), pc)
        a_active_word = changes[0] + phases[0][:2] + _inverse(changes[0])
        c_active_word = changes[2] + phases[2][:2] + _inverse(changes[2])
        c_silent_word = changes[2] + phases[2][2:] + _inverse(changes[2])
        a_active = native_columns(4, a_active_word)
        c_active = native_columns(4, c_active_word)
        c_silent = native_columns(4, c_silent_word)
        boolean_residual = eye(16)
        for column in range(16):
            u, v, r, s = [(column >> q) & 1 for q in range(4)]
            if (u ^ v) * r * s:
                boolean_residual[column][column] = -CO
        self.assertEqual(multiply(adjoint(c_active), multiply(eb, a_active)),
                         multiply(c0, boolean_residual))
        self.assertEqual(_conjugate(c0, c_silent), c_silent)
        j, clean_s = _conjugate(ea, pa), _conjugate(ec, pc)
        self.assertEqual(_conjugate(eb, j), clean_s)
        self.assertEqual(len(_support(j, 4)), 4)
        self.assertEqual(len(_support(clean_s, 4)), 4)
        for label in labels[1]:
            # All four middle orbits are negative singletons: m_pair=0.
            self.assertEqual(_conjugate(j, _tau(label, 4)), scale(_tau(label, 4), -1))

        actual = native_columns(4, word)
        expected = [[CZ for _ in range(16)] for _ in range(16)]
        for column in range(16):
            clean, t, r, s = [(column >> q) & 1 for q in range(4)]
            flip = clean * (r ^ s)
            expected[column ^ (flip << 1)][column] = (
                -CI if flip else (CO if not (t * r * s) else -CO))
        self.assertEqual(actual, expected)
        self.assertEqual(actual, multiply(adjoint(kc), multiply(adjoint(ec),
                         multiply(eb, multiply(ea, ka)))))
        inverse = _inverse(word)
        self.assertEqual(native_columns(4, inverse), adjoint(expected))
        for band, matrix in zip(bands, (ea, eb, ec)):
            self.assertEqual(native_columns(4, _inverse(band)), adjoint(matrix))
        self.assertEqual(Counter(g[0] for g in word),
                         {'H': 14, 'CX': 34, 'T': 8, 'TDG': 4, 'SDG': 5})
        self.assertEqual(Counter(g[0] for g in inverse),
                         {'H': 14, 'CX': 34, 'T': 4, 'TDG': 8, 'S': 5})
        for schedule, cliffords in ((word, 63), (inverse, 53)):
            expanded = []
            for gate in schedule:
                expanded.extend([('S', gate[1])] * 3 if gate[0] == 'SDG' else [gate])
            self.assertEqual(sum(g[0] not in ('T', 'TDG') for g in expanded), cliffords)
            self.assertEqual(native_columns(4, expanded),
                             expected if schedule is word else adjoint(expected))

        # Every clean-odd Pauli acquires a nonconstant XOR permutation.
        # This rules out a shared Pauli flip, not only a commuting-band flip.
        for label in range(256):
            if not label & 1:
                continue
            image = _conjugate(adjoint(actual), _tau(label, 4))
            masks = set()
            for column in range(16):
                rows = [row for row in range(16) if image[row][column] != CZ]
                self.assertEqual(len(rows), 1)
                masks.add(rows[0] ^ column)
            self.assertGreater(len(masks), 1)

    def test_surviving_center_and_failed_middle_clifford_passage(self):
        # Logical order t,r,s. E is I for r+s=0 and omega^-1 S_X for r+s=1.
        change = [('CX', 1, 0), ('CX', 2, 0), ('H', 0)]
        tx = [('H', 0), ('T', 0), ('H', 0)]
        # Conjugate Z_t to X_t Z_r Z_s using the opposite coordinate order.
        dressed = _inverse(change) + [('TDG', 0)] + change
        e_word = tx + dressed
        e = native_columns(3, e_word)
        f = eye(8)
        f[7][7] = -CO
        self.assertEqual(_conjugate(adjoint(e), f), f)
        self.assertEqual(native_columns(3, _inverse(e_word)), adjoint(e))
        r, s, h = _tau(16, 3), _tau(32, 3), _tau(8, 3)
        center = multiply(r, s)
        j_t = _conjugate(e, h)
        self.assertEqual(_conjugate(j_t, h), multiply(h, center))
        self.assertEqual(_conjugate(j_t, r), r)
        self.assertEqual(_conjugate(j_t, s), s)
        j_r = _conjugate(e, _tau(2, 3))
        self.assertEqual(_conjugate(j_r, h), multiply(_tau(9, 3), center))
        middle_image = _conjugate(f, j_r)
        cz_rs = eye(8)
        for column in (6, 7):
            cz_rs[column][column] = -CO
        diagnostic = multiply(_tau(9, 3), multiply(center, cz_rs))
        self.assertEqual(_conjugate(middle_image, h), diagnostic)
        self.assertGreater(len(_support(diagnostic, 3)), 1)
        # The actual final derivative remains an elementary CZ_(t,s).
        cz_ts = eye(8)
        for column in (5, 7):
            cz_ts[column][column] = -CO
        self.assertEqual(multiply(_tau(2, 3), _conjugate(f, _tau(2, 3))), cz_ts)

        # The eight-dimensional encoded abelian algebra has twelve Pauli
        # support labels. A fixed Clifford cannot make it a rank-three
        # Pauli algebra, whose support would contain only eight labels.
        encoded = [_conjugate(adjoint(e), _tau(z << 3, 3)) for z in range(8)]
        union = set().union(*(_support(operator, 3) for operator in encoded))
        self.assertEqual(len(union), 12)
        self.assertEqual(_rank(union), 4)
        for left in encoded:
            for right in encoded:
                self.assertEqual(multiply(left, right), multiply(right, left))
        # Linear independence follows exactly by undoing the actual encoder.
        for z, operator in enumerate(encoded):
            self.assertEqual(_conjugate(e, operator), _tau(z << 3, 3))


if __name__ == '__main__':
    unittest.main()
