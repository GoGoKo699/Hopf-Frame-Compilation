"""Exact evidence for rank-one derivative / complete-indicator reductions.

Indicator and D calls are explicit reduction oracles, not free native gates.
The shared-control Toffoli phase core, in contrast, is emitted literally.
All principal and auxiliary basis columns are retained in the fixtures.
"""
from collections import Counter
import unittest

from tests.exact_matrices import CO
from tests.test_common_probe_source import _native, _syndrome
from tests.test_joint_lookup_pauli_source import _inverse


def _phase_core(k):
    """Four physical T bands; c=0, z=1..k, R=k+1..2k."""
    z, r = list(range(1, k + 1)), list(range(k + 1, 2 * k + 1))
    power = k % 8
    control = [('S', 0)] * (power // 2) + ([('T', 0)] if power % 2 else [])
    if power == 7:
        control = [('TDG', 0)]
    first = control + [('T', j) for j in z + r]
    bands = [[g for g in first if g[0] in ('T', 'TDG')]]
    word = list(first)
    families = [([('CX', 0, j) for j in z + r], 'TDG', z + r),
                ([('CX', a, b) for a, b in zip(z, r)], 'TDG', r),
                ([g for a, b in zip(z, r) for g in [('CX', 0, b), ('CX', a, b)]],
                 'T', r)]
    for basis, sign, targets in families:
        band = [(sign, j) for j in targets]
        word += basis + band + list(reversed(basis))
        bands.append(band)
    return word, bands


def _phase_column(word, basis):
    """Exact Z8 evaluator for the literal monomial native phase core."""
    phase = 0
    for gate, *wires in word:
        if gate == 'CX':
            if basis >> wires[0] & 1:
                basis ^= 1 << wires[1]
        else:
            phase += {'T': 1, 'TDG': -1, 'S': 2}[gate] * (basis >> wires[0] & 1)
    return basis, phase % 8


def _g(bank, p, address):
    w = bank & ~1
    return (w.bit_count() ^ (w >> (address ^ _syndrome(w, p)))) & 1


def _d(state, p):
    mask = (1 << (1 << p)) - 1
    out = {}
    for basis, value in state.items():
        bank, address = (basis >> p) & mask, basis & ((1 << p) - 1)
        out[basis] = -value if (bank & 1) and _g(bank, p, address) else value
    return out


def _fwht(values):
    out = list(values)
    step = 1
    while step < len(out):
        for start in range(0, len(out), 2 * step):
            for j in range(start, start + step):
                out[j], out[j + step] = out[j] + out[j + step], out[j] - out[j + step]
        step *= 2
    return out


class RankOneIndicatorEquivalenceTests(unittest.TestCase):
    def test_four_band_native_phase_and_actual_adjoint_every_column(self):
        k = 7
        word, _ = _phase_core(k)
        inverse = _inverse(word)
        for basis in range(1 << (2 * k + 1)):
            c = basis & 1
            overlap = sum(((basis >> j) & 1) * ((basis >> (k + j)) & 1)
                          for j in range(1, k + 1))
            target = 4 * c * overlap % 8
            moved, phase = _phase_column(word, basis)
            restored, backward = _phase_column(inverse, moved)
            self.assertEqual((moved, phase), (basis, target))
            self.assertEqual((restored, (phase + backward) % 8), (basis, 0))

    def test_native_hadamard_conjugation_really_is_shared_control_toffoli(self):
        for k in (1, 3):
            core, _ = _phase_core(k)
            wall = [('H', j) for j in range(k + 1, 2 * k + 1)]
            word = wall + core + wall
            reverse = _inverse(word)
            for basis in range(1 << (2 * k + 1)):
                target = basis
                if basis & 1:
                    for j in range(1, k + 1):
                        if basis >> j & 1:
                            target ^= 1 << (k + j)
                actual = _native(word, {basis: CO})
                self.assertEqual(actual, {target: CO})
                self.assertEqual(_native(reverse, actual), {basis: CO})

    def test_all_auxiliary_bank_columns_cancel_including_arbitrary_row_zero(self):
        p, m = 3, 8
        for s in range(m):
            for z in range(1 << m):
                c, w = z & 1, z & ~1
                x = s ^ _syndrome(w, p)
                for auxiliary in range(1 << m):
                    # Chronological Omega, A, actual Omega-adjoint, actual A-adjoint.
                    first_phase = (auxiliary >> x) & 1
                    moved = auxiliary ^ (w if c else 0)
                    second_phase = (moved >> x) & 1
                    restored = moved ^ (w if c else 0)
                    parity_phase = c * (w.bit_count() & 1)
                    self.assertEqual(restored, auxiliary)
                    self.assertEqual(first_phase ^ second_phase ^ parity_phase,
                                     c * _g(z, p, s))
                    self.assertEqual(x ^ _syndrome(w, p), s)

    def test_forward_echo_and_full_bank_hadamards_give_indicator(self):
        p, m = 3, 8
        # Literal K, C2†, D†, X0, D, C2, K, Z_pi before outer bank Hs.
        shear = [('CX', p + j, i) for j in range(m) for i in range(p) if j >> i & 1]
        c2 = [('H', p)] + [('CX', p, p + j) for j in range(1, m)]
        for x in range(m):
            diagonal = []
            for bank in range(1 << m):
                basis = x | bank << p
                state = _native(shear + _inverse(c2), {basis: CO})
                state = _d(state, p)
                state = _native([('H', p), ('S', p), ('S', p), ('H', p)], state)
                state = _d(state, p)
                state = _native(c2 + shear + [('S', p + j) for j in range(m)] * 2, state)
                sign = -1 if bank >> x & 1 else 1
                self.assertEqual(state, {basis: sign * CO})
                diagonal.append(sign)
            for column in range(1 << m):
                transformed = _fwht([diagonal[row] * (-1 if (row & column).bit_count() & 1 else 1)
                                     for row in range(1 << m)])
                expected = [0] * (1 << m)
                expected[column ^ (1 << x)] = 1 << m
                self.assertEqual(transformed, expected)

    def test_nondiagonal_address_and_work_gauge_cancels_in_actual_echo(self):
        p, m = 3, 8
        work = p + m
        # J changes address and work and couples work to c; [J,X_c]=0.
        gauge = [('H', 0), ('H', work), ('CX', work, p)]
        for basis in range(1 << (p + m + 1)):
            state = _d({basis: CO}, p)  # V† = J† D†.
            state = _native(_inverse(gauge), state)
            state = _native([('H', p), ('S', p), ('S', p), ('H', p)], state)
            state = _native(gauge, state)  # V = D J.
            state = _d(state, p)
            bank = (basis >> p) & ((1 << m) - 1)
            sign = -1 if _g(bank, p, basis & (m - 1)) else 1
            self.assertEqual(state, {basis ^ (1 << p): sign * CO})

    def test_literal_band_and_total_gate_ledgers(self):
        for p in (3, 4, 5):
            m, k = 1 << p, (1 << p) - 1
            core, bands = _phase_core(k)
            counts = Counter(g[0] for g in core)
            self.assertEqual(counts['T'] + counts['TDG'], 6 * k + 1)
            self.assertEqual(counts['CX'], 10 * k)
            self.assertEqual(len(core) + 2 * k, 18 * k + 1)
            self.assertEqual(len(bands), 4)
            for band in bands:
                self.assertEqual(len({g[1] for g in band}), len(band))
            self.assertEqual(2 * (18 * k + 1) + p * m + 4 * m + 3 * k,
                             p * m + 43 * m - 37)


if __name__ == '__main__':
    unittest.main()
