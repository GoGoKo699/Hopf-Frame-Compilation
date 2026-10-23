"""Exact guards for paid signal interleaving and its failed lookup sectors."""
from fractions import Fraction
import unittest

import numpy as np

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, add, adjoint, eye, multiply, native_columns,
    scale, zeros,
)
from tests.test_bent_gram_source import ry_t
from tests.test_half_bilinear_gram_source import (
    _binary_image, _ccs_word, _private_toffoli, _walsh_integer,
)
from tests.test_majorana_walsh_query import _inverse
from tests.test_nilpotent_pencil_observable_source import _field_matrix, _power
from tests.test_phase_signal_query_limits import _counts
from tests.test_t_wall_xor_gram_source import _p1_mixer


def _ct_word():
    query = _private_toffoli(0, 1, 2, 5)
    return query + [('T', 2)] + query + [('TDG', 2)] + _ccs_word()


def _root_powers():
    omega = CQ2(Q2(0, Fraction(1, 2)), Q2(0, Fraction(1, 2)))
    powers = [CO]
    for _ in range(7):
        powers.append(powers[-1] * omega)
    return powers


def _p1_parts(address):
    powers, diagonal = _root_powers(), zeros(16)
    for value in range(16):
        z, t = value & 3, value >> 2
        shifted = t ^ (z & 2)
        exponent = (shifted ^ (address * (z >> 1))).bit_count() - shifted.bit_count()
        diagonal[value][value] = powers[exponent % 8]
    mixer = native_columns(4, _p1_mixer())
    e = multiply(multiply(adjoint(mixer), diagonal), mixer)
    f = multiply(native_columns(4, [('Z', 0), ('Z', 2)]), adjoint(diagonal))
    return e, f


def _control(matrix):
    size = len(matrix)
    result = zeros(2 * size)
    for i in range(size):
        result[i][i] = CO
        for j in range(size):
            result[size + i][size + j] = matrix[i][j]
    return result


def _one_clean_block(e, f):
    a, b = scale(e, CI), scale(f, CI)
    alpha = CO + CI * CQ2(Q2(-1, 1))
    beta = CQ2(Q2(1, 1)) + CI
    numerator = add(add(add(scale(eye(len(e)), alpha), scale(a, alpha.conj())),
                        scale(b, beta)), scale(multiply(b, a), -beta.conj()))
    return scale(numerator, (CO + CI) * CQ2(Fraction(1, 8)))


def _integer_receipt(p, family):
    """Exact accepted norm at x=1,z=e_p,t=0 in Z[omega].

    Two full Walsh walls give denominator 2^(2n). Coefficients are
    kept in int64; at p<=2 their bounds are far below 2^63. Convert
    to Python integers before the quadratic norm expressions.
    """
    m, n = 1 << p, p * (1 << p)
    size, dimension = 1 << n, 1 << (2 * n)
    basis = ((3, 2) if p == 1 else
             _field_matrix([[_power(u, j) for u in range(m)] for j in range(m)]))
    images = np.array([_binary_image(basis, z) for z in range(size)], dtype=np.int64)
    labels = np.arange(dimension, dtype=np.int64)
    z, t = labels & (size - 1), labels >> n
    permutation = (images[z] ^ t) | (z << n)
    weights = np.array([value.bit_count() for value in range(size)], dtype=np.int64)
    shifted = t ^ (z & ~1)
    h = weights[shifted ^ (z >> p)] - weights[shifted]
    basis_v = 1 << p
    flip = 0 if family == 'one' else (1 if family == 'two' else 1 | (1 << p))
    values = np.zeros((4, dimension), dtype=np.int64)
    values[0, basis_v ^ flip] = 1

    def phase(old, powers):
        new = np.zeros_like(old)
        for power in range(8):
            mask = powers == power
            for degree in range(4):
                target = (degree + power) % 8
                new[target % 4, mask] = ((-1 if target >= 4 else 1)
                                        * old[degree, mask])
        return new

    def polynomial(old, coefficients):
        new = np.zeros_like(old)
        for degree in range(4):
            for shift, coefficient in enumerate(coefficients):
                target = (degree + shift) % 8
                new[target % 4] += ((-1 if target >= 4 else 1)
                                    * coefficient * old[degree])
        return new

    def norm(old, denominator):
        rational = radical = 0
        for column in old.T:
            a, b, c, d = map(int, column)
            rational += a*a + b*b + c*c + d*d
            radical += a*b + b*c + c*d - a*d
        return Fraction(rational, denominator**2), Fraction(radical, denominator**2)

    routed = np.zeros_like(values)
    routed[:, permutation] = values
    values = routed
    for row in values:
        _walsh_integer(row)
    values = phase(values, h % 8)
    for row in values:
        _walsh_integer(row)
    values = values[:, permutation]
    denominator = 1 << (2 * n)
    assert norm(values, denominator) == (1, 0)
    fpower = (-h + 4 * ((z ^ t) & 1)) % 8
    fvalues = phase(values, fpower)
    basis_column = np.zeros_like(values)
    basis_column[0, basis_v] = denominator
    if family != 'one':
        accepted = (basis_column + phase(basis_column, fpower)
                    + values - fvalues)
        return norm(accepted, 4 * denominator)
    # alpha I + i conj(alpha) E + i beta F + conj(beta) F E,
    # followed by the literal prefactor (1+i)/8.
    accepted = (polynomial(basis_column, [1, 1, -1, 1])
                + polynomial(polynomial(values, [1, -1, 1, -1]), [0, 0, 1, 0])
                + polynomial(polynomial(phase(basis_column, fpower),
                                        [1, 1, 1, -1]), [0, 0, 1, 0])
                + polynomial(fvalues, [1, 1, -1, -1]))
    accepted = polynomial(accepted, [1, 0, 1, 0])
    return norm(accepted, 8 * denominator)


class SignalInterleavedKernelSourceTests(unittest.TestCase):
    def test_literal_controlled_t_returns_all_four_dirty_helpers(self):
        word = _ct_word()
        self.assertEqual(_counts(word), (70, 88))
        self.assertEqual(_counts(_inverse(word)), (70, 88))
        self.assertTrue(all(gate[1] in (2, 3, 4, 5) for gate in word
                            if gate[0] in ('T', 'TDG')))
        actual = native_columns(6, word)
        for column in range(64):
            scalar = _root_powers()[1] if column & 3 == 3 else CO
            for row in range(64):
                self.assertEqual(actual[row][column], scalar if row == column else CZ)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))

    def test_controlled_phase_keeps_coherent_address_and_every_dirty_column(self):
        # signal,z1,t0,address, then four private helpers. Unaffected
        # bank coordinates are literal identity factors in this fixture.
        mapping = {0: 0, 1: 2, 2: 4, 3: 5, 4: 6, 5: 7}
        ct = [tuple([gate[0]] + [mapping[wire] for wire in gate[1:]])
              for gate in _ct_word()]
        query = _private_toffoli(3, 1, 2, 7)
        word = _inverse(ct) + query + ct + _inverse(query)
        self.assertEqual(_counts(word), (156, 196))
        actual = native_columns(8, word)
        powers = _root_powers()
        for column in range(256):
            signal, z1, t0, address = [(column >> bit) & 1 for bit in range(4)]
            exponent = signal * ((t0 ^ (address * z1)) - t0)
            for row in range(256):
                self.assertEqual(actual[row][column],
                                 powers[exponent % 8] if row == column else CZ)

    def test_complete_signal_word_matches_polynomial_and_base_normalization(self):
        signal = 4
        rotation = native_columns(5, ry_t(signal))
        s_gate = native_columns(5, [('S', signal)])
        h_gate = native_columns(5, [('H', signal)])
        for address in (0, 1):
            e, f = _p1_parts(address)
            source = eye(32)
            for operation in (rotation, s_gate, h_gate, _control(scale(e, CI)),
                              h_gate, _control(scale(f, CI)), h_gate, rotation):
                source = multiply(operation, source)
            accepted = [row[:16] for row in source[:16]]
            expected = _one_clean_block(e, f)
            self.assertEqual(accepted, expected)
            gram = multiply(adjoint(accepted), accepted)
            if address == 0:
                for row in range(16):
                    for column in range(16):
                        equal = not ((column ^ (column >> 2)) & 1)
                        scalar = CQ2(Fraction(1, 2)) if row == column and equal else CZ
                        self.assertEqual(gram[row][column], scalar)
            else:
                self.assertEqual(gram[2][2], CQ2(Fraction(1, 2)))

    def test_exact_small_sector_receipts(self):
        expected = {
            (1, 'two'): (Fraction(1, 4), Fraction(1, 8)),
            (1, 'paired'): (Fraction(3, 8), Fraction(1, 16)),
            (1, 'one'): (Fraction(1, 2), Fraction(0)),
            (2, 'two'): (Fraction(7991, 32768), Fraction(1391, 16384)),
            (2, 'paired'): (Fraction(9305, 32768), Fraction(1, 16)),
            (2, 'one'): (Fraction(1197, 2048), Fraction(-155, 8192)),
        }
        for inputs, receipt in expected.items():
            self.assertEqual(_integer_receipt(*inputs), receipt)

    def test_three_branch_seeded_word_full_polynomial_and_trace(self):
        signal, size = 4, 16
        identity, half = eye(size), CQ2(Fraction(1, 2))
        inverse_root_two = CQ2(Q2(0, Fraction(1, 2)))
        p0 = native_columns(4, [('Z', 0), ('Z', 2)])
        projector = scale(add(identity, p0), half)
        complement = add(identity, scale(projector, -CO))
        rotation = native_columns(5, ry_t(signal))
        rotation_cube = native_columns(5, ry_t(signal, cube=True))
        signal_h = native_columns(5, [('H', signal)])
        controlled_p0 = native_columns(5, [
            ('H', 0), ('CX', signal, 0), ('H', 0),
            ('H', 2), ('CX', signal, 2), ('H', 2),
        ])
        base = multiply(multiply(rotation_cube, controlled_p0), adjoint(rotation))
        omega = _root_powers()[1]

        def trace(matrix):
            return sum((matrix[j][j] for j in range(size)), CZ) / CQ2(size)

        for address in (0, 1):
            e, f = _p1_parts(address)
            outer = multiply(multiply(signal_h, _control(f)), signal_h)
            f_plus = scale(add(identity, f), half)
            f_minus = scale(add(identity, scale(f, -CO)), half)
            delta = add(e, scale(identity, -CO))
            for seed_kind in ('zero', 'real', 'complex'):
                seed = rotation if seed_kind != 'zero' else eye(32)
                if seed_kind == 'complex':
                    seed = multiply(seed, native_columns(5, [('S', signal)]))
                word = eye(32)
                for operation in (adjoint(outer), seed, _control(e),
                                  adjoint(seed), outer, base):
                    word = multiply(operation, word)
                accepted = [row[:size] for row in word[:size]]

                # Use the accepted polynomial, independent of the 32x32 word.
                # The scalar phase in the native seed cancels between L
                # and L-dagger, so no nested square roots are needed.
                seed_small = (native_columns(1, ry_t(0))
                              if seed_kind != 'zero' else eye(2))
                if seed_kind == 'complex':
                    seed_small = multiply(seed_small, native_columns(1, [('S', 0)]))
                seed_inverse = adjoint(seed_small)
                l_op = add(scale(f_plus, seed_inverse[0][1]),
                           scale(f_minus, seed_inverse[1][1]))
                h_op = add(scale(f_minus, seed_inverse[0][1]),
                           scale(f_plus, seed_inverse[1][1]))
                coefficient = add(
                    scale(multiply(projector, add(l_op, scale(h_op, -CO))),
                          inverse_root_two),
                    multiply(complement, h_op))
                expected = scale(add(scale(projector, inverse_root_two),
                                     multiply(multiply(coefficient, delta),
                                              adjoint(l_op))), omega)
                self.assertEqual(accepted, expected)
                gram = multiply(adjoint(accepted), accepted)
                if address == 0:
                    self.assertEqual(gram, scale(projector, half))
                    continue

                mu = trace(e)
                m_op = scale(add(f, adjoint(f)), half)
                a2 = trace(multiply(m_op, m_op))
                correlation = trace(multiply(multiply(multiply(e, m_op),
                                                       adjoint(e)), m_op))
                j_moment = trace(multiply(multiply(multiply(m_op, adjoint(delta)),
                                                   p0), delta))
                if seed_kind != 'complex':
                    q, b = ((inverse_root_two, inverse_root_two)
                            if seed_kind == 'real' else (CZ, CO))
                    trace_formula = (CQ2(Fraction(1, 4)) + (CO - mu) / CQ2(4)
                                     + (q*b*j_moment - b*b*(correlation
                                         + (CO - CQ2(2)*mu)*a2)) / CQ2(8))
                    self.assertEqual(trace(gram), trace_formula)

                alpha, beta = seed_inverse[0][1], seed_inverse[1][1]
                cross = alpha.conj() * beta
                a_bloch = cross + cross.conj()
                b_bloch = (cross - cross.conj()) / CI
                c_bloch = beta.conj()*beta - alpha.conj()*alpha
                self.assertEqual(a_bloch*a_bloch + b_bloch*b_bloch
                                 + c_bloch*c_bloch, CO)
                n_op = scale(add(f, scale(adjoint(f), -CO)), -CI*half)
                axis = add(scale(m_op, c_bloch), scale(n_op, -b_bloch))
                axis_square = trace(multiply(axis, axis))
                axis_correlation = trace(multiply(multiply(multiply(e, axis),
                                                            adjoint(e)), axis))
                axis_j = trace(multiply(multiply(multiply(axis, adjoint(delta)),
                                                  p0), delta))
                general_trace = (CQ2(Fraction(1, 4)) + (CO - mu) / CQ2(4)
                                 + (a_bloch*axis_j - axis_correlation
                                    - (CO - CQ2(2)*mu)*axis_square) / CQ2(8))
                self.assertEqual(trace(gram), general_trace)
                expected_traces = {
                    'zero': CQ2(Q2(Fraction(9, 32), Fraction(-1, 64))),
                    'real': CQ2(Q2(Fraction(11, 32), Fraction(-5, 128))),
                    'complex': CQ2(Q2(Fraction(41, 128), Fraction(-1, 32))),
                }
                self.assertEqual(trace(gram), expected_traces[seed_kind])
                expected_wrong = (CQ2(Q2(Fraction(7, 16), Fraction(1, 16)))
                                  if seed_kind == 'real' else half)
                self.assertEqual(gram[2][2], expected_wrong)

    def test_three_branch_blind_sector_and_uniform_seed_certificate(self):
        # At x=1 the same F=+1, P0=+1, Q=-1 vector is protected
        # for a genuinely changed first branch as well as the actual E.
        _, f = _p1_parts(1)
        changed_e = native_columns(4, [('T', 0), ('H', 1), ('CX', 1, 3)])
        signal_h = native_columns(5, [('H', 4)])
        outer = multiply(multiply(signal_h, _control(f)), signal_h)
        sandwich = multiply(multiply(outer, _control(changed_e)), adjoint(outer))
        self.assertEqual([row[8] for row in sandwich],
                         [CO if j == 8 else CZ for j in range(32)])

        # Exact positivity certificate for every real-axis seed, p>=3.
        self.assertLess(Fraction(7, 8)**21, Fraction(1, 16))
        self.assertLess((1 + Fraction(1, 1 << 21)) / 2, Fraction(9, 16))
        self.assertEqual(Fraction(15, 16)*Fraction(7, 16)/4
                         - (Fraction(3, 4) + Fraction(1, 16))/16,
                         Fraction(53, 1024))

    def test_partial_h_full_signal_gram_and_empty_intersection(self):
        partial_h = native_columns(4, [('H', 1), ('H', 3)])
        distinguished_x = native_columns(4, [('X', 0)])
        other_x = native_columns(4, [('X', 1)])
        rotation = native_columns(5, ry_t(4))
        signal_s = native_columns(5, [('S', 4)])
        signal_h = native_columns(5, [('H', 4)])
        half = CQ2(Fraction(1, 2))
        for address in (0, 1):
            e, f = _p1_parts(address)
            changed_f = multiply(multiply(partial_h, f), partial_h)
            source = eye(32)
            for operation in (rotation, signal_s, signal_h,
                              _control(scale(e, CI)), signal_h,
                              _control(scale(changed_f, CI)), signal_h, rotation):
                source = multiply(operation, source)
            accepted = [row[:16] for row in source[:16]]
            self.assertEqual(accepted, _one_clean_block(e, changed_f))
            gram = multiply(adjoint(accepted), accepted)

            real_f = scale(add(changed_f, adjoint(changed_f)), half)
            imag_f = scale(add(changed_f, scale(adjoint(changed_f), -CO)), -CI * half)
            weighted_e = add(scale(e, CO + CI), scale(adjoint(e), CO - CI))
            expanded = add(scale(eye(16), half), scale(weighted_e, CQ2(Fraction(-1, 8))))
            expanded = add(expanded, scale(add(multiply(multiply(adjoint(e), imag_f), e),
                                              scale(imag_f, -CO)),
                                          CQ2(Q2(0, Fraction(1, 8)))))
            expanded = add(expanded, scale(add(scale(multiply(real_f, e), CO - CI),
                                              scale(multiply(adjoint(e), real_f), CO + CI)),
                                          CQ2(Fraction(1, 8))))
            self.assertEqual(gram, expanded)
            self.assertEqual(multiply(e, distinguished_x), multiply(distinguished_x, e))
            self.assertEqual(multiply(multiply(distinguished_x, changed_f), distinguished_x),
                             scale(changed_f, -CO))
            twirled = scale(add(gram, multiply(multiply(distinguished_x, gram), distinguished_x)),
                            half)
            expected_twirl = add(scale(eye(16), half), scale(weighted_e, CQ2(Fraction(-1, 8))))
            self.assertEqual(twirled, expected_twirl)
            self.assertEqual(multiply(gram, other_x), multiply(other_x, gram))
            self.assertEqual(gram[0][0], half)
            self.assertEqual(gram[2][2], half)
            if address:
                charge = native_columns(4, [('Z', 0), ('Z', 1), ('Z', 2), ('Z', 3)])
                rotated_charge = multiply(multiply(partial_h, charge), partial_h)
                support = [label for label in range(16) if f[label][label] == CO]
                self.assertEqual(support, [0, 5, 8, 13])
                self.assertFalse(any(rotated_charge[i][j] != CZ
                                     for i in support for j in support))

    def test_local_rotation_full_signal_gram_and_trace_obstruction(self):
        # The fixed-address pi/4 rotation is exact. The parameterized
        # nonzero-address echo in the proof uses paid approximations;
        # this fixture does not certify that asymptotic synthesis.
        signal_rotation = native_columns(5, ry_t(4))
        signal_s = native_columns(5, [('S', 4)])
        signal_h = native_columns(5, [('H', 4)])
        local_rotation = native_columns(4, ry_t(0))
        z0 = native_columns(4, [('Z', 0)])
        x0 = native_columns(4, [('X', 0)])
        half = CQ2(Fraction(1, 2))
        for address in (0, 1):
            e, f = _p1_parts(address)
            v = local_rotation if address else eye(16)
            changed_f = multiply(multiply(adjoint(v), f), v)
            source = eye(32)
            for operation in (signal_rotation, signal_s, signal_h,
                              _control(scale(e, CI)), signal_h,
                              _control(scale(changed_f, CI)), signal_h,
                              signal_rotation):
                source = multiply(operation, source)
            accepted = [row[:16] for row in source[:16]]
            self.assertEqual(accepted, _one_clean_block(e, changed_f))
            gram = multiply(adjoint(accepted), accepted)
            if not address:
                self.assertEqual(gram[0][0], half)
                self.assertEqual(gram[1][1], CZ)
                continue
            rest = multiply(z0, f)
            expected_f = multiply(scale(add(z0, scale(x0, -CO)),
                                        CQ2(Q2(0, Fraction(1, 2)))), rest)
            self.assertEqual(changed_f, expected_f)
            self.assertEqual(gram[2][2], CQ2(Fraction(7, 16)))
            self.assertEqual(sum((gram[k][k] for k in range(16)), CZ),
                             CQ2(Q2(6, -1)))
            self.assertEqual(sum((multiply(e, x0)[k][k] for k in range(16)), CZ), CZ)
            y0 = scale(multiply(x0, z0), CI)
            for pauli in (x0, y0, z0):
                phase = multiply(pauli, rest)
                real_part = scale(add(phase, adjoint(phase)), half)
                product = multiply(real_part, e)
                self.assertEqual(sum((product[k][k] for k in range(16)), CZ), CZ)

    def test_x1_orbit_parity_and_attained_spectral_gap(self):
        powers = _root_powers()
        for p in (1, 2, 3):
            m, n = 1 << p, p * (1 << p)
            v = sum(1 << (block * p) for block in range(m))
            w = v ^ 1
            # At x=1 the binary B query shifts field blocks by p bits.
            self.assertEqual((w >> p).bit_count(), m - 1)
            self.assertEqual((w >> p).bit_count() % 2, 1)
            k = 1 if p == 1 else 3
            z = ((1 << k) - 1) << p
            t = z & ~1  # Literal t=Pz, so a=t+Pz=0.
            a = t ^ (z & ~1)
            exponent = (a ^ (z >> p)).bit_count() - a.bit_count()
            self.assertLessEqual(k, n - p)
            self.assertEqual(exponent, k)
            eigenvalue = powers[exponent % 8]
            spectral_error = (CO - (CO + CI) * eigenvalue / CQ2(2)
                              - (CO - CI) * eigenvalue.conj() / CQ2(2)) / CQ2(4)
            expected = CQ2(Fraction(1, 4)) if p == 1 else CQ2(Q2(Fraction(1, 4), Fraction(1, 4)))
            self.assertEqual(spectral_error, expected)

    def test_exact_branch_moments_and_global_positive_floor(self):
        powers = _root_powers()
        a = CQ2(Q2(Fraction(1, 2), Fraction(1, 4)))
        for p in (1, 2):
            n, ell = p * (1 << p), p * ((1 << p) - 1)
            counts = [0] * 8
            # Literal x=1 XOR-carry phase on every pair of carrier labels.
            # No product-distribution formula is used in this enumeration.
            for z in range(1 << n):
                for t in range(1 << n):
                    shifted = t ^ (z & ~1)
                    h = (shifted ^ (z >> p)).bit_count() - shifted.bit_count()
                    counts[h % 8] += 1
            dimension = CQ2(1 << (2 * n))
            first = sum((CQ2(count) * powers[k]
                         for k, count in enumerate(counts)), CZ) / dimension
            second = sum((CQ2(count) * powers[(2 * k) % 8]
                          for k, count in enumerate(counts)), CZ) / dimension
            mu = CO
            for _ in range(ell):
                mu *= a
            rho = CQ2(Fraction(1, 1 << ell))
            self.assertEqual(first, mu)
            self.assertEqual(second, rho)
            m_square = h_square = CZ
            for k, count in enumerate(counts):
                phase = powers[k]
                real_part = (phase + phase.conj()) / CQ2(2)
                h_part = ((CO - CI) * phase + (CO + CI) * phase.conj()) / CQ2(2)
                m_square += CQ2(count) * real_part * real_part
                h_square += CQ2(count) * h_part * h_part
            self.assertEqual(m_square / dimension, (CO + rho) / CQ2(2))
            self.assertEqual(h_square / dimension, CO)

        # Exact rational certificates at the largest global moments for p>=3.
        # sqrt(2)<3/2 proves a<7/8; no floating-point gap is inferred.
        self.assertLess(Fraction(2), Fraction(3, 2) ** 2)
        self.assertLess(Fraction(7, 8) ** 21, Fraction(1, 16))
        self.assertLess((1 + Fraction(1, 1 << 21)) / 2, Fraction(3, 4) ** 2)
        mu_upper = Fraction(1, 8) + Fraction(7, 8) * Fraction(1, 16)
        self.assertEqual(mu_upper, Fraction(23, 128))
        second_upper = (1 + Fraction(1, 8) + Fraction(1, 1 << 21)) / 2
        self.assertLess(second_upper, Fraction(49, 64) ** 2)
        self.assertEqual((1 - mu_upper - Fraction(49, 64)) / 4, Fraction(7, 512))


if __name__ == '__main__':
    unittest.main()
