"""Exact full-space p=2 separator for an enlarged co-moving alphabet.

The sum of two co-moving stages, fixed zeta and dual-zeta CNOT networks,
global selector Hadamards and the global flip has disjoint lookup-charge
spectra.  Arithmetic reduces Z[i,1/2] modulo 65537 with i mapped to 256.
Coprime reduced characteristic polynomials certify a nonzero resultant
over Q(i), hence disjoint spectra over C and no cross-charge intertwiner.

The additional zero-address fixtures check the all-bank commutator and
anticommuting-Pauli steps of the separate analytic all-p argument.  Neither
result is a normalized accepted source, uniform gap, or native filter.
"""

from collections import Counter
import unittest

from tests.test_co_moving_lookup_stage import _Matching


_PRIME = 65537
_IMAGINARY = 256
_H_SCALE = pow(4, -1, _PRIME)
_PAIRS = ((0, 1), (1, 0))


def _address_image(x, control, target):
    return x ^ (((x >> control) & 1) << target)


def _relabel(bits, control, target):
    return sum(((bits >> row) & 1) << _address_image(row, control, target)
               for row in range(4))


def _zeta(bits, transpose=False):
    return sum((sum((bits >> u) & 1 for u in range(4)
                    if (v & u == v if transpose else v & u == u)) % 2) << v
               for v in range(4))


def _zeta_literal(bits, transpose=False, p=2):
    for coordinate in range(p):
        for row in range(1 << p):
            if row >> coordinate & 1:
                continue
            control, target = row, row ^ (1 << coordinate)
            if transpose:
                control, target = target, control
            bits ^= ((bits >> control) & 1) << target
    return bits


def _xy_orbits():
    unseen = {(x, bank) for x in range(4) for bank in range(16)}
    answer = []
    while unseen:
        first = min(unseen)
        orbit = {first}
        while True:
            expanded = orbit | {
                (_address_image(x, c, t), _relabel(bank, c, t))
                for x, bank in orbit for c, t in _PAIRS}
            if expanded == orbit:
                break
            orbit = expanded
        unseen.difference_update(orbit)
        answer.append(((first[1] >> first[0]) & 1, sorted(orbit)))
    return answer


def _stage_column(x, bank, selectors, control, target):
    active = bank & selectors
    phase = sum(((x ^ row) >> control) & 1
                for row in range(4) if active >> row & 1)
    if active.bit_count() % 2:
        x = _address_image(x, control, target)
        bank = _relabel(bank, control, target)
        selectors = _relabel(selectors, control, target)
    return (x, bank, selectors), pow(_IMAGINARY, phase, _PRIME)


def _separator_matrix(orbit):
    """R_(0->1) + R_(1->0) + M + N + H_all + F, without weights."""
    index = {xy: j for j, xy in enumerate(orbit)}
    size = len(orbit) * 16
    matrix = [[0] * size for _ in range(size)]
    for j, (x, bank) in enumerate(orbit):
        for selectors in range(16):
            column = 16 * j + selectors
            for control, target in _PAIRS:
                (xx, yy, ss), coefficient = _stage_column(
                    x, bank, selectors, control, target)
                matrix[16 * index[xx, yy] + ss][column] += coefficient
            for ss in (_zeta(selectors), _zeta(selectors, True), selectors ^ 15):
                matrix[16 * j + ss][column] += 1
            for ss in range(16):
                sign = -1 if (ss & selectors).bit_count() % 2 else 1
                matrix[16 * j + ss][column] += sign * _H_SCALE
    return [[entry % _PRIME for entry in row] for row in matrix]


def _characteristic_polynomial(matrix):
    """Ascending coefficients, using upper-Hessenberg similarity over F_p."""
    a = [row[:] for row in matrix]
    n = len(a)
    for m in range(1, n):
        pivot = next((i for i in range(m, n) if a[i][m - 1]), None)
        if pivot is None:
            continue
        a[m], a[pivot] = a[pivot], a[m]
        for row in a:
            row[m], row[pivot] = row[pivot], row[m]
        inverse = pow(a[m][m - 1], -1, _PRIME)
        for i in range(m + 1, n):
            factor = a[i][m - 1] * inverse % _PRIME
            if factor:
                a[i] = [(v - factor * w) % _PRIME for v, w in zip(a[i], a[m])]
                for row in a:
                    row[m] = (row[m] + factor * row[i]) % _PRIME
    polynomials = [[1]]
    for i in range(1, n + 1):
        current = [0] + polynomials[-1][:]
        for j, value in enumerate(polynomials[-1]):
            current[j] = (current[j] - a[i - 1][i - 1] * value) % _PRIME
        subdiagonal_product = 1
        for j in range(i - 1, 0, -1):
            subdiagonal_product = subdiagonal_product * a[j][j - 1] % _PRIME
            coefficient = subdiagonal_product * a[j - 1][i - 1] % _PRIME
            for k, value in enumerate(polynomials[j - 1]):
                current[k] = (current[k] - coefficient * value) % _PRIME
        polynomials.append(current)
    return polynomials[-1]


def _polynomial_gcd(a, b):
    a, b = a[:], b[:]
    while b:
        while len(a) >= len(b):
            coefficient = a[-1] * pow(b[-1], -1, _PRIME) % _PRIME
            offset = len(a) - len(b)
            for j, value in enumerate(b):
                a[offset + j] = (a[offset + j] - coefficient * value) % _PRIME
            while a and a[-1] == 0:
                a.pop()
        a, b = b, a
    inverse = pow(a[-1], -1, _PRIME)
    return [value * inverse % _PRIME for value in a]


class CoMovingSectorSeparatorTests(unittest.TestCase):
    def test_zero_address_commutator_retains_every_background_and_phase(self):
        families = (
            (2, (((0, 1),), ((1, 0),))),
            (3, (((0, 1),), ((1, 2),), ((2, 0),))),
            (4, (((0, 1), (2, 3)), ((1, 0), (3, 2)), ((0, 2), (1, 3)))),
        )
        for p, matchings in families:
            m, full = 1 << p, (1 << (1 << p)) - 1
            if p == 2:
                banks, selector_inputs = range(16), range(16)
            else:
                even_rows = sum(1 << row for row in range(m) if row % 2 == 0)
                bases = {0, full, even_rows, 1 << (m - 1),
                         (1 << 1) | (1 << (m // 2)) | (1 << (m - 1))}
                banks = sorted(bases | {bank ^ 1 for bank in bases})
                selectors = {0, full, even_rows, full ^ even_rows}
                selectors.update(1 << row for row in range(m))
                selectors.update(full ^ (1 << row) for row in range(m))
                selector_inputs = sorted(selectors)
            for pairs in matchings:
                matching = _Matching(p, pairs)
                for bank in banks:
                    for selectors in selector_inputs:
                        for borrowed in (0, 1):
                            # Right-to-left chronology of R X_0 R^dagger X_0.
                            actual = matching.stage_literal(
                                (0, bank, selectors ^ 1, borrowed, 0), inverse=True)
                            x, y, s, bit, phase = actual
                            actual = matching.stage_literal((x, y, s ^ 1, bit, phase))
                            expected = (0, matching.relabel(bank),
                                        matching.relabel(selectors), borrowed, 0)
                            if not bank & 1:
                                expected = (0, bank, selectors, borrowed, 0)
                            self.assertEqual(actual, expected)

    def test_zeta_flip_probes_and_anticommuting_orbit_pairs(self):
        for p in range(2, 6):
            m, full = 1 << p, (1 << (1 << p)) - 1
            top = m - 1
            if p <= 3:
                selector_inputs = range(1 << m)
            else:
                selector_inputs = sorted({0, full}
                                         | {1 << row for row in range(m)}
                                         | {full ^ (1 << row) for row in range(m)})
            first = _Matching(p, ((0, 1),))
            second = _Matching(p, ((1, 0),))
            first_image, second_image = first.images[top], second.images[top]
            self.assertEqual(len({top, first_image, second_image}), 3)
            x_mask = (1 << top) | (1 << first_image)
            z_mask = (1 << top) | (1 << second_image)
            self.assertEqual((x_mask & z_mask).bit_count() % 2, 1)
            for selectors in selector_inputs:
                mfm = _zeta_literal(_zeta_literal(selectors, p=p) ^ full, p=p)
                nfn = _zeta_literal(_zeta_literal(selectors, True, p) ^ full, True, p)
                self.assertEqual(mfm, selectors ^ 1)
                self.assertEqual(nfn, selectors ^ (1 << top))
                # The two Pauli supports meet only at top.  This retains
                # their literal anticommutation sign on arbitrary inputs.
                self.assertEqual(((selectors & z_mask).bit_count()
                                  + ((selectors ^ x_mask) & z_mask).bit_count()) % 2, 1)
                for marked in (0, 1):
                    state, sign = selectors, 0
                    for operation in ('Z', 'X', 'Z', 'X'):
                        if marked and operation == 'Z':
                            sign ^= (state & z_mask).bit_count() % 2
                        elif marked:
                            state ^= x_mask
                    self.assertEqual((state, sign), (selectors, marked))

    def test_complete_stage_columns_and_literal_selector_networks(self):
        self.assertEqual(_IMAGINARY ** 2 % _PRIME, _PRIME - 1)
        for control, target in _PAIRS:
            matching = _Matching(2, ((control, target),))
            for x in range(4):
                for bank in range(16):
                    for selectors in range(16):
                        output, coefficient = _stage_column(
                            x, bank, selectors, control, target)
                        for borrowed in (0, 1):
                            literal = matching.stage_literal(
                                (x, bank, selectors, borrowed, 0))
                            self.assertEqual(output, literal[:3])
                            self.assertEqual(literal[3], borrowed)
                            self.assertEqual(coefficient, pow(_IMAGINARY, literal[4], _PRIME))
        for selectors in range(16):
            for transpose in (False, True):
                self.assertEqual(_zeta(selectors, transpose),
                                 _zeta_literal(selectors, transpose))
                self.assertEqual(_zeta(_zeta(selectors, transpose), transpose), selectors)
            for output in range(16):
                # Exact H_all M H_all column, including both 1/4 factors.
                amplitude = sum((-1) ** (
                    (v & selectors).bit_count() + (_zeta(v) & output).bit_count())
                                for v in range(16))
                self.assertEqual(amplitude, 16 * (output == _zeta(selectors, True)))

    def test_modular_characteristic_polynomial_companions(self):
        for n in range(1, 9):
            expected = [(j * j + 3 * j + n + 1) % _PRIME for j in range(n)] + [1]
            companion = [[0] * n for _ in range(n)]
            for j in range(n - 1):
                companion[j + 1][j] = 1
            for j in range(n):
                companion[j][-1] = -expected[j] % _PRIME
            for ordering in (list(range(n)), list(reversed(range(n)))):
                matrix = [[companion[i][j] for j in ordering] for i in ordering]
                self.assertEqual(_characteristic_polynomial(matrix), expected)
        self.assertEqual(_polynomial_gcd([0, 0, 1], [0, 1]), [0, 1])
        self.assertEqual(_polynomial_gcd([1, 0, 1], [0, 1]), [1])

    def test_full_two_bit_charge_spectra_are_disjoint(self):
        orbits = _xy_orbits()
        self.assertEqual(Counter(len(orbit) for _, orbit in orbits), {1: 4, 3: 12, 6: 4})
        self.assertEqual({xy for _, orbit in orbits for xy in orbit},
                         {(x, bank) for x in range(4) for bank in range(16)})
        records = []
        for charge, orbit in orbits:
            self.assertTrue(all((bank >> x) & 1 == charge for x, bank in orbit))
            matrix = _separator_matrix(orbit)
            polynomial = _characteristic_polynomial(matrix)
            self.assertEqual(len(polynomial), len(orbit) * 16 + 1)
            self.assertEqual(polynomial[-1], 1)
            self.assertEqual(polynomial[-2], -sum(matrix[i][i] for i in range(len(matrix))) % _PRIME)
            records.append((charge, orbit, polynomial))
        self.assertEqual([sum(len(orbit) * 16 for c, orbit, _ in records if c == charge)
                          for charge in (0, 1)], [512, 512])
        comparisons = 0
        for charge, _, polynomial in records:
            if charge:
                continue
            for other_charge, _, other in records:
                if other_charge:
                    self.assertEqual(_polynomial_gcd(polynomial, other), [1])
                    comparisons += 1
        self.assertEqual(comparisons, 100)


if __name__ == '__main__':
    unittest.main()
