"""Exact guards for a normalized, but wrong-sector, dense half-phase trial.

The p=2 fixture uses integer Walsh numerators. It is finite evidence,
not an all-size source exclusion or a scalable native emitter.
"""
from fractions import Fraction
import unittest

import numpy as np

from tests.exact_matrices import (
    CI, CO, CQ2, CZ, Q2, adjoint, eye, multiply, native_columns,
)
from tests.test_gram_four_call_readout import _partial_gram_source
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts
from tests.test_nilpotent_pencil_observable_source import (
    _field_matrix, _power, _mul, _transpose,
)


def _private_toffoli(a, b, target, phase):
    controls = (a, b, target)
    vertices = (0, 1, 3, 2, 6, 7, 5, 4)
    word, old = [('H', target)], 0
    for vertex in vertices:
        difference = old ^ vertex
        if difference:
            word.append(('CX', controls[difference.bit_length() - 1], phase))
        word.append(('TDG' if vertex.bit_count() % 2 else 'T', phase))
        old = vertex
    word.extend([('CX', target, phase), ('H', target)])
    return word


def _ccs_word():
    a, b, c, d, e, phase = range(6)
    td = _private_toffoli(a, b, d, phase)
    te = _private_toffoli(a, b, e, phase)
    l = [('T', d), ('CX', c, d), ('TDG', d), ('CX', c, d)]
    ccz = [('H', e)] + _private_toffoli(c, d, e, phase) + [('H', e)]
    return td + l + td + _inverse(l) + te + ccz + te + ccz


def _binary_image(rows, value):
    return sum(((row & value).bit_count() % 2) << bit
               for bit, row in enumerate(rows))


def _walsh_integer(values):
    distance = 1
    while distance < len(values):
        view = values.reshape(-1, 2 * distance)
        left, right = view[:, :distance].copy(), view[:, distance:].copy()
        view[:, :distance] = left + right
        view[:, distance:] = left - right
        distance *= 2


def _probe_expectation(p, left_sign=1, right_sign=1, address_shear=False):
    """Exact p<=2 carrier expectation at field x=1, z=e_p, t=0.

    Four full Hadamard walls give denominator 2^(4n). Intermediate
    numerators have magnitude <=2^(4n)<=2^32, safely inside int64.
    Convert to Python integers BEFORE squaring for the exact receipt.
    """
    m, n = 1 << p, p * (1 << p)
    size, dimension = 1 << n, 1 << (2 * n)
    basis = ((3, 2) if p == 1 else
             _field_matrix([[_power(u, j) for u in range(m)] for j in range(m)]))
    images = np.array([_binary_image(basis, z) for z in range(size)], dtype=np.int64)
    labels = np.arange(dimension, dtype=np.int64)
    z, t = labels & (size - 1), labels >> n
    # At x=1, bar N_x^T z = z xor (z >> p).  This is the
    # actual right-composed shear in R_x = R Q_{bar N_x^T}.
    sheared_t = t ^ z ^ (z >> p) if address_shear else t
    permutation = (images[z] ^ sheared_t) | (z << n)
    exponent = np.zeros(dimension, dtype=np.int64)
    # At x=1, multiplication by x is the identity in the binary basis.
    for block in range(m - 1):
        for bit in range(p):
            exponent += (((z >> (block * p + bit)) & 1)
                         * ((t >> ((block + 1) * p + bit)) & 1))
    real, imag = np.zeros(dimension, dtype=np.int64), np.zeros(dimension, dtype=np.int64)
    real[1 << p] = 1

    def phase(sign):
        old_real, old_imag = real.copy(), imag.copy()
        powers = (sign * exponent) % 4
        for power in (1, 2, 3):
            mask = powers == power
            if power == 1:
                real[mask], imag[mask] = -old_imag[mask], old_real[mask]
            elif power == 2:
                real[mask], imag[mask] = -old_real[mask], -old_imag[mask]
            else:
                real[mask], imag[mask] = old_imag[mask], -old_real[mask]

    def mixer(inverse=False):
        nonlocal real, imag
        if not inverse:
            new_real, new_imag = np.zeros_like(real), np.zeros_like(imag)
            new_real[permutation], new_imag[permutation] = real, imag
            real, imag = new_real, new_imag
        _walsh_integer(real)
        _walsh_integer(imag)
        if inverse:
            real, imag = real[permutation], imag[permutation]

    # V=D^-left E^-right D^left E^right; E=R† D R.
    mixer()
    phase(right_sign)
    mixer(True)
    phase(left_sign)
    mixer()
    phase(-right_sign)
    mixer(True)
    phase(-left_sign)
    squares = [int(a) ** 2 + int(b) ** 2 for a, b in zip(real, imag)]
    denominator = 1 << (8 * n)
    norm = Fraction(sum(squares), denominator)
    numerator = sum((-1 if (int(zi) ^ int(ti)) & 1 else 1) * square
                    for zi, ti, square in zip(z, t, squares))
    return norm, Fraction(numerator, denominator)


def _fwht(a, offset, count):
    distance = 1 << offset
    for _ in range(count):
        view = a.reshape(-1, 2 * distance)
        left, right = view[:, :distance].copy(), view[:, distance:].copy()
        view[:, :distance], view[:, distance:] = left + right, left - right
        distance *= 2


def _p2_arrays():
    p, m, n = 2, 4, 8
    size, dimension = 1 << n, 1 << (2*n+p)
    labels = np.arange(dimension, dtype=np.int64)
    z, t = labels & (size-1), (labels >> n) & (size-1)
    basis = _field_matrix([[_power(u,j) for u in range(m)] for j in range(m)])
    images = np.array([_binary_image(basis, i) for i in range(size)], dtype=np.int64)
    permutation, exponent = np.empty(dimension,dtype=np.int64), np.zeros(dimension,dtype=np.int64)
    for address in range(m):
        pencil_t = _transpose(_field_matrix([[int(i==j)^(_mul(address,1) if i==j+1 else 0)
                                             for j in range(m)] for i in range(m)]))
        pimages = np.array([_binary_image(pencil_t,i) for i in range(size)],dtype=np.int64)
        start, end = address*size*size, (address+1)*size*size
        zz, tt = z[start:end], t[start:end]
        permutation[start:end] = (images[zz]^tt^pimages[zz]) | (zz << n) | (address << (2*n))
        for i,row in enumerate(pencil_t):
            row ^= 1 << i
            for j in range(n):
                if row & (1 << j):
                    exponent[start:end] += ((zz >> i)&1) * ((tt >> j)&1)
    return z,t,permutation,exponent


def _address_h_p2_probe(mode):
    assert mode in ('left','right')
    p,m,n = 2,4,8
    size,dimension = 1 << n, 1 << (2*n+p)
    z,t,permutation,exponent = _p2_arrays()
    real,imag = np.zeros(dimension,dtype=np.int64),np.zeros(dimension,dtype=np.int64)
    real[(1 << (2*n)) + (1 << p)] = 1

    def phase(sign):
        old_r,old_i = real.copy(),imag.copy()
        powers = (sign*exponent)%4
        for power in (1,2,3):
            mask = powers == power
            if power == 1: real[mask],imag[mask] = -old_i[mask],old_r[mask]
            elif power == 2: real[mask],imag[mask] = -old_r[mask],-old_i[mask]
            else: real[mask],imag[mask] = old_i[mask],-old_r[mask]

    def wall(offset,count):
        _fwht(real,offset,count)
        _fwht(imag,offset,count)

    def permute(inverse=False):
        if inverse:
            new_r,new_i = real[permutation],imag[permutation]
        else:
            new_r,new_i = np.empty_like(real),np.empty_like(imag)
            new_r[permutation],new_i[permutation] = real,imag
        real[:],imag[:] = new_r,new_i

    def mixer(inverse=False):
        if mode == 'left':
            if inverse: wall(0,2*n+p); permute(True)
            else: permute(); wall(0,2*n+p)
        elif inverse:
            wall(0,2*n); permute(True); wall(2*n,p)
        else:
            wall(2*n,p); permute(); wall(0,2*n)

    # Chronological E, D, E^dagger, D^dagger.
    mixer(); phase(1); mixer(True); phase(1)
    mixer(); phase(-1); mixer(True); phase(-1)
    # Four complete (2n+p)-wire H walls. Intermediate values <=2^36.
    # Squared norm denominator2^72: convert BEFORE multiplying.
    denominator = 1 << (4*(2*n+p))
    squares = [int(a)**2+int(b)**2 for a,b in zip(real,imag)]
    numerator = sum((1-2*((int(zz)^int(tt))&1))*s for zz,tt,s in zip(z,t,squares))
    return Fraction(sum(squares),denominator),Fraction(numerator,denominator)


def _address_h_p1_full(mode):
    assert mode in ('left','right')
    dimension = 32
    full_h = np.array([[(-1)**((a&b).bit_count()%2) for b in range(dimension)]
                      for a in range(dimension)],dtype=np.int64)
    address_h,bank_h = np.zeros((dimension,dimension),dtype=np.int64),np.zeros((dimension,dimension),dtype=np.int64)
    for a in range(dimension):
        for b in range(dimension):
            if a&15 == b&15: address_h[a,b] = (-1)**((a>>4)&(b>>4))
            if a>>4 == b>>4: bank_h[a,b] = (-1)**(((a&15)&(b&15)).bit_count()%2)
    q = np.array([((y>>4)&1)*(y&1)*((y>>3)&1) for y in range(dimension)],dtype=np.int64)
    dr,di = np.diag(1-q),np.diag(q)
    probe = np.diag([(-1)**((y&1)^((y>>2)&1)) for y in range(dimension)]).astype(np.int64)
    permutation = np.zeros((dimension,dimension),dtype=np.int64)
    for y in range(dimension):
        x,z,t = y>>4,y&3,(y>>2)&3
        ez = ((z^(z>>1))&1)|(z&2)
        t ^= ez if x else z
        permutation[(ez^t)|(z<<2)|(x<<4),y] = 1
    raw = full_h@permutation if mode=='left' else bank_h@permutation@address_h

    def multiply(a,b):
        ar,ai=a; br,bi=b
        return ar@br-ai@bi,ar@bi+ai@br

    def adjoint(a): return a[0].T,-a[1].T
    e=(raw.T@dr@raw,raw.T@di@raw)  # denominator32
    v=multiply(multiply(multiply((dr,-di),adjoint(e)),(dr,di)),e)  # denominator1024
    norm=multiply(adjoint(v),v)
    denominator=32**4
    assert np.array_equal(norm[0],np.eye(dimension,dtype=np.int64)*denominator)
    assert not norm[1].any()
    a=multiply(multiply(adjoint(v),(probe,np.zeros_like(probe))),v)
    assert a[0][9,29] == a[1][9,29] == -denominator//4
    assert a[0][18,18] == denominator and a[1][18,18] == 0
    return a,denominator


class HalfBilinearGramSourceTests(unittest.TestCase):
    def test_base_source_is_literal_h_or_x_on_both_signal_inputs(self):
        word = _partial_gram_source()
        self.assertEqual(_counts(word), (2, 19))
        self.assertEqual(_counts(_inverse(word)), (2, 19))
        actual = native_columns(3, word)
        inverse_root = CQ2(Q2(0, Fraction(1, 2)))
        for column in range(8):
            equal = (column & 1) == ((column >> 1) & 1)
            signal = (column >> 2) & 1
            for row in range(8):
                expected = CZ
                if row & 3 == column & 3:
                    if equal:
                        expected = inverse_root * (-1) ** (signal * ((row >> 2) & 1))
                    elif row == column ^ 4:
                        expected = CO
                self.assertEqual(actual[row][column], expected)

    def test_literal_52_t_ccs_returns_every_dirty_input(self):
        word = _ccs_word()
        self.assertEqual(_counts(word), (52, 68))
        self.assertTrue(all(gate[1] in (3, 4, 5) for gate in word
                            if gate[0] in ('T', 'TDG')))
        actual = native_columns(6, word)
        self.assertEqual(native_columns(6, _inverse(word)), adjoint(actual))
        for column in range(64):
            scalar = CI if column & 7 == 7 else CO
            for row in range(64):
                self.assertEqual(actual[row][column], scalar if row == column else CZ)

    def test_binary_affine_control_keeps_arbitrary_phase_helper_return(self):
        word = []
        for gate in _ccs_word():
            if gate[0] == 'CX' and gate[1] == 0:
                word.extend([gate, ('CX', 6, gate[2])])
            else:
                word.append(gate)
        self.assertEqual(_counts(word), (52, 84))
        actual = native_columns(7, word)
        self.assertEqual(native_columns(7, _inverse(word)), adjoint(actual))
        for column in range(128):
            affine = ((column >> 0) ^ (column >> 6)) & 1
            scalar = CI if affine and column & 6 == 6 else CO
            for row in range(128):
                self.assertEqual(actual[row][column], scalar if row == column else CZ)

    def test_dense_commutator_p1_exact_sector_failure_and_orientations(self):
        orientations = ((1, 1), (1, -1), (-1, 1), (-1, -1),
                        (2, 1), (1, 2), (2, 2))
        for left, right in orientations:
            norm, expectation = _probe_expectation(1, left, right)
            self.assertEqual(norm, 1)
            self.assertEqual(expectation, 1)
            # At x=1,z=e_1,t=0 the desired paired lookup character is -1.
            self.assertEqual(abs(expectation + 1) / 4, Fraction(1, 2))

    def test_dense_commutator_p2_exact_walsh_receipts(self):
        expected = {(1, 1): Fraction(19, 64), (1, -1): Fraction(39, 128),
                    (-1, 1): Fraction(39, 128), (-1, -1): Fraction(19, 64),
                    (2, 1): Fraction(3, 32), (1, 2): Fraction(1, 4),
                    (2, 2): Fraction(0)}
        for orientations, value in expected.items():
            norm, expectation = _probe_expectation(2, *orientations)
            self.assertEqual(norm, 1)
            self.assertEqual(expectation, value)
            self.assertGreaterEqual(abs(expectation + 1) / 4, Fraction(1, 4))
        self.assertEqual((1 + expected[(1, 1)]) / 4, Fraction(83, 256))

    def test_address_hadamard_complete_p1_return_and_sector_failure(self):
        # Full32-column unitarity and off-address entries are exact integers.
        for mode in ('left', 'right'):
            probe, denominator = _address_h_p1_full(mode)
            self.assertEqual(Fraction(int(probe[0][9, 29]), denominator),
                             Fraction(-1, 4))
            self.assertEqual(Fraction(int(probe[1][9, 29]), denominator),
                             Fraction(-1, 4))
            self.assertEqual(Fraction(int(probe[0][18, 18]), denominator), 1)

    def test_address_hadamard_p2_full_address_walsh_receipts(self):
        expected = {'left': Fraction(20011, 65536),
                    'right': Fraction(6913, 16384)}
        gaps = {'left': Fraction(85547, 262144),
                'right': Fraction(23297, 65536)}
        for mode, value in expected.items():
            norm, expectation = _address_h_p2_probe(mode)
            self.assertEqual((norm, expectation), (1, value))
            self.assertEqual((1 + expectation) / 4, gaps[mode])

    def test_address_sheared_p1_complete_commutator_is_identity(self):
        # z0,z1,t0,t1.  At x=1, E=bar N_1^T=[[1,1],[0,1]].
        # Emit S_1, then M, then H^4, without assuming the cancellation.
        shear = [('CX', 0, 2), ('CX', 1, 2), ('CX', 1, 3)]
        swap = []
        for z, t in ((0, 2), (1, 3)):
            swap.extend([('CX', z, t), ('CX', t, z), ('CX', z, t)])
        matrix_m = swap + [('CX', 2, 0), ('CX', 3, 0), ('CX', 3, 1)]
        mixer = native_columns(4, shear + matrix_m
                               + [('H', bit) for bit in range(4)])
        expected_mixer = native_columns(4, swap
                                        + [('H', bit) for bit in range(4)])
        self.assertEqual(mixer, expected_mixer)
        diagonal = [[(CI if column & 9 == 9 else CO)
                     if row == column else CZ
                     for column in range(16)] for row in range(16)]
        crossed = multiply(multiply(adjoint(mixer), diagonal), mixer)
        self.assertEqual(multiply(diagonal, crossed),
                         multiply(crossed, diagonal))
        commutator = multiply(multiply(multiply(adjoint(diagonal),
                                                adjoint(crossed)), diagonal),
                              crossed)
        self.assertEqual(commutator, eye(16))
        # At x=0 the phase D_0=I, so the full coherent-address V is I.
        norm, expectation = _probe_expectation(1, address_shear=True)
        self.assertEqual((norm, expectation), (1, 1))

    def test_address_sheared_p2_exact_walsh_receipts(self):
        expected = {(1, 1): Fraction(147, 512),
                    (1, -1): Fraction(129, 512),
                    (-1, 1): Fraction(129, 512),
                    (-1, -1): Fraction(147, 512),
                    (2, 1): Fraction(3, 128),
                    (1, 2): Fraction(0), (2, 2): Fraction(0)}
        for signs, value in expected.items():
            norm, expectation = _probe_expectation(2, *signs,
                                                   address_shear=True)
            self.assertEqual((norm, expectation), (1, value))
            self.assertGreaterEqual((1 + expectation) / 4, Fraction(1, 4))
        self.assertEqual((1 + expected[(1, 1)]) / 4,
                         Fraction(659, 2048))


if __name__ == '__main__':
    unittest.main()
