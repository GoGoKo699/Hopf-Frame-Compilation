"""Paid eight-wire field recoupler with four projective-line corrections.

Every column, literal phase, arbitrary borrowed input and actual inverse is
checked exactly. The exponentially many all-size lines are not a uniform
polynomial-budget decoder.
"""
from collections import Counter
from functools import lru_cache
import unittest

import numpy as np

from tests.test_recursive_field_shear_decoder import (_tower_product as mul,
    _tower_inverse as inv, _tower_trace as trace, _tower_shear,
    _recursive_decoder_word, _apply_decoder)
from tests.test_dirty_indicator_depth import _parity_basis
from tests.test_joint_lookup_pauli_source import _inverse, _bands
from tests.test_three_layer_native_mechanisms import _relative_toffoli
from tests.test_field_shear_intertwiner import _native_x, _translation
from tests.test_flat_field_bank_native import full_intertwiner
from tests.exact_matrices import Q2
from tests.test_one_band_copy_feasibility import _conjugate

DIRECTIONS = (1, 4, 5, 6, 7)


def linear_word(columns, wires):
    rows = [sum(((column >> row) & 1) << j for j, column in enumerate(columns))
            for row in range(len(wires))]
    return _parity_basis(list(wires), rows)


def fourier_word(level, offset):
    degree = 1 << (level - 1)
    # The trace-pairing matrix is symmetric and invertible over F2.
    columns = [sum(trace(level, mul(level, 1 << i, 1 << j)) << i
                   for i in range(degree)) for j in range(degree)]
    # H T = T^-T H. Native action desired H T; emit trace map then H.
    return linear_word(columns, range(offset, offset + degree)) + [
        ('H', offset + j) for j in range(degree)]


def _raw_linear_cell(direction):
    d0, d1 = direction & 3, direction >> 2
    norm = mul(2, d0, d0) ^ mul(2, d0, d1) ^ mul(2, 2, mul(2, d1, d1))
    cols_a = []
    for j in range(4):
        uv = mul(3, direction, 1 << j)
        cols_a.append((mul(3, 4, uv) >> 2)
                      | (mul(2, uv >> 2, inv(2, norm)) << 2))
    cols_b = [mul(3, direction, 1 << j) for j in range(4)]
    shift = mul(3, direction, mul(2, mul(2, d0, d1), inv(2, norm)) << 2)
    return (linear_word(cols_a, range(4, 8)) + fourier_word(2, 4)
            + linear_word(cols_b, range(4, 8)), shift)


def _cancel_clifford_pairs(word):
    out = []
    for gate in word:
        if out and gate == out[-1] and gate[0] in ('H', 'CX'):
            out.pop()
        else:
            out.append(gate)
    return out


def c3x(a, b, c, target, borrow):
    """16 T, eight bands, one returned arbitrary borrow.

    U=(-i)^(ab) Z_borrow^(ab), V=(-i)^(c target) X_borrow^(c target).
    Their literal commutator is (-1)^(ab c target) I_borrow; the scalar
    relative-Toffoli phases cancel. Target Hadamards turn C3Z into C3X.
    """
    assert len({a, b, c, target, borrow}) == 5
    h = [('H', borrow)]
    u = h + _relative_toffoli(a, b, borrow) + h
    v = _relative_toffoli(c, target, borrow)
    return _cancel_clifford_pairs([('H', target)] + u + v
                                 + _inverse(u) + _inverse(v) + [('H', target)])


def _baseline_word():
    # C1=H_F(first) A1, with the trace-pairing wall absorbed into A1.
    return [('CX',4,5),('CX',7,5),('CX',5,4),('CX',6,5),
            ('H',4),('H',5)]


def _controlled_shift(control, vector):
    if not vector:
        return []
    pivot = (vector & -vector).bit_length() - 1
    # F sends X_pivot to X_vector; exact common conjugation costs no T.
    f = [('CX', 4 + pivot, 4 + j)
         for j in range(4) if j != pivot and (vector >> j) & 1]
    borrow = 1 - control
    return _inverse(f) + c3x(2, 3, control, 4 + pivot, borrow) + f


def pauli_diagonalizer(label, pivot=None):
    x,z=label&15,label>>4
    support=x|z
    assert support
    if pivot is None:pivot=(support&-support).bit_length()-1
    assert (support>>pivot)&1
    word=[('CX',4+j,4+pivot) for j in range(4)
          if j!=pivot and (support>>j)&1]
    for j in range(4):
        if (x>>j)&1:
            word.append(('H',4+j))
            if (z>>j)&1:word.append(('S',4+j))
    return word, 4+pivot


# J(A,B;P)=T(A)T(B) P T(B)^dagger T(A)^dagger.
# The two axes commute, are independent, and both anticommute with P.
# Chronological J0,J1 equals sign*(C_D C1^dagger), with no residual Pauli.
_PAIRED_REFLECTIONS = {
    4: (((19,50,32),(156,228,112)),1),
    5: (((156,228,64),(19,50,32)),1),
    6: (((89,247,1),(19,50,13)),1),
    7: (((19,50,32),(89,247,87)),-1),
}


def pair_diagonalizer(left,right):
    first,pivot=pauli_diagonalizer(left)
    low_inverse=[(g[0],*(j-4 for j in g[1:])) for g in _inverse(first)]
    pulled=_conjugate({right:Q2(1)},low_inverse,4)
    assert len(pulled)==1
    label,coefficient=next(iter(pulled.items()))
    assert coefficient in (Q2(1),Q2(-1))
    other=(label&15)|(label>>4)
    other &= ~(1<<(pivot-4))
    assert other
    second_pivot=(other&-other).bit_length()-1
    second,target=pauli_diagonalizer(label,second_pivot)
    sign=_native_x(target) if coefficient==Q2(-1) else []
    return sign+second+first,(pivot,target)


def paired_t(left,right):
    k,wires=pair_diagonalizer(left,right)
    return _inverse(k)+[('T',wire) for wire in wires]+k


def selected_correction(direction):
    pairs,sign=_PAIRED_REFLECTIONS[direction]
    word=[]
    for left,right,probe in pairs:
        tau=paired_t(left,right)
        k,pivot=pauli_diagonalizer(probe)
        relative_p=(_inverse(k)+[('H',pivot)]
                    +_relative_toffoli(2,3,pivot)+[('H',pivot)]+k)
        word += _inverse(tau)+relative_p+tau
    # Two relative-control phases multiply to (-1)^(ab).
    if sign>0:
        word += [('H',3),('CX',2,3),('H',3)]
    return word


def line_factor(direction):
    # The unconditional baseline C1 already handles its own line.
    assert direction in _PAIRED_REFLECTIONS
    coords=linear_word([mul(3,inv(3,direction),1<<j)
                       for j in range(4)],range(4))
    negative=_native_x(2)+_native_x(3)
    selected=selected_correction(direction)
    shift=_raw_linear_cell(direction)[1]
    selected += _controlled_shift(0,shift)
    selected += _controlled_shift(1,mul(3,2,shift))
    return coords+negative+selected+negative+_inverse(coords)


@lru_cache(None)
def recoupler_word():
    cell = _baseline_word() + [g for direction in DIRECTIONS[1:]
                                      for g in line_factor(direction)]
    # (d,a)->(a,a+d), two CNOTs per coordinate.
    sum_wall = [g for j in range(4) for g in [('CX', j, 4+j), ('CX', 4+j, j)]]
    return (fourier_word(3, 4) + cell + sum_wall
            + fourier_word(2, 2) + fourier_word(2, 6))


def decoder_word():
    bottom = [g for gate in full_intertwiner()
              for g in [gate, (gate[0], *(j + 4 for j in gate[1:]))]]
    return recoupler_word() + bottom


def ledger(word):
    c = Counter(g[0] for g in word)
    return {'T': c['T'] + c['TDG'], 'bands': len(_bands(word)),
            'gates': len(word), 'width': max(j for g in word for j in g[1:]) + 1,
            'counts': dict(c)}


def ring_matrix(word, width=8):
    """Exact Z[omega]/2^den matrices, omega^4=-1, vectorized integer arithmetic.

    All H gates multiply by (omega-omega^3)/2. Normalize common factors
    two after each H. Overflow guards precede every arithmetic operation.
    """
    size = 1 << width
    a = np.zeros((4, size, size), dtype=np.int64)
    a[0] = np.eye(size, dtype=np.int64)
    den = 0
    rows = np.arange(size)
    for index, gate in enumerate(word):
        name, *wires = gate
        assert np.max(np.abs(a)) < (1 << 57), ('overflow guard', index, den)
        if name == 'CX':
            c, t = wires
            perm = rows ^ (((rows >> c) & 1) << t)
            a = a[:, perm, :]
        elif name == 'H':
            t, = wires
            lo = rows[(rows & (1 << t)) == 0]
            hi = lo ^ (1 << t)
            b = a.copy()
            b[:, lo, :] = a[:, lo, :] + a[:, hi, :]
            b[:, hi, :] = a[:, lo, :] - a[:, hi, :]
            # (omega - omega^3) b, with omega^4=-1.
            a = np.array([b[1] - b[3], b[0] + b[2],
                          b[1] + b[3], b[2] - b[0]])
            den += 1
            while den and not np.any(a & 1):
                a //= 2
                den -= 1
        else:
            assert name in ('S', 'T', 'TDG')
            t, = wires
            selected = rows[(rows & (1 << t)) != 0]
            exponent = {'S':2, 'T':1, 'TDG':7}[name]
            b = a[:, selected, :].copy()
            for power in range(4):
                target = (power + exponent) % 4
                a[target, selected, :] = ((-1) ** ((power + exponent)//4)) * b[power]
    return a, den


def equal_ring(left, right):
    a, da = left; b, db = right
    common = max(da, db)
    for value,shift in ((a,common-da),(b,common-db)):
        assert int(np.max(np.abs(value))).bit_length()+shift<63
    return np.array_equal(a * (1 << (common-da)), b * (1 << (common-db)))


def adjoint_ring(value):
    a, d = value
    return np.array([a[0].T, -a[3].T, -a[2].T, -a[1].T]), d


def _expected_nonzero():
    columns = [_apply_decoder(_recursive_decoder_word(3)[:6], {j:1})
               for j in range(256)]
    assert {n for _, n in columns} == {10}
    a = np.zeros((4,256,256),dtype=np.int64)
    for col,(entries,n) in enumerate(columns):
        for row,value in entries.items():a[0,row,col]=value
    return a,5



class ProjectiveFieldRecouplerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.recoupler = ring_matrix(recoupler_word())
        cls.recoupler_inverse = ring_matrix(_inverse(recoupler_word()))
        cls.decoder = ring_matrix(decoder_word())
        cls.decoder_inverse = ring_matrix(_inverse(decoder_word()))

    def test_field_fourier_native_walls(self):
        for level in (2,3):
            width=1 << (level-1)
            size=1 << width
            actual=ring_matrix(fourier_word(level,0),width)
            expected=np.zeros_like(actual[0])
            for row in range(size):
                for col in range(size):
                    expected[0,row,col]=(-1)**trace(level,mul(level,row,col))
            self.assertTrue(equal_ring(actual,(expected,width//2)))

    def test_every_nonzero_projective_fiber_and_displacement(self):
        low=lambda word:[(g[0],*(j-4 for j in g[1:])) for g in word]
        for direction in DIRECTIONS:
            direct,shift=_raw_linear_cell(direction)
            for scalar in (1,2,3):
                value=mul(3,scalar,direction)
                scaled,scaled_shift=_raw_linear_cell(value)
                self.assertTrue(equal_ring(ring_matrix(low(direct),4),
                                           ring_matrix(low(scaled),4)))
                self.assertEqual(scaled_shift,mul(3,scalar,shift))
        self.assertTrue(equal_ring(ring_matrix(low(_baseline_word()),4),
                                   ring_matrix(low(_raw_linear_cell(1)[0]),4)))

    def test_quartic_borrow_every_column_and_actual_inverse(self):
        word=c3x(0,1,2,3,4)
        expected=np.zeros((4,32,32),dtype=np.int64)
        for col in range(32):
            expected[0,col^(8 if col&7==7 else 0),col]=1
        for actual in (word,_inverse(word)):
            self.assertTrue(equal_ring(ring_matrix(actual,5),(expected,0)))
        self.assertEqual((ledger(word)['T'],ledger(word)['bands'],len(word)),
                         (16,8,46))

    def test_literal_paired_reflections_and_joint_physical_bands(self):
        low=lambda word:[(g[0],*(j-4 for j in g[1:])) for g in word]
        for direction,(pairs,sign) in _PAIRED_REFLECTIONS.items():
            word=[]
            for left,right,probe in pairs:
                joint,wires=pair_diagonalizer(left,right)
                self.assertNotEqual(*wires)
                for wire,label in zip(wires,(left,right)):
                    self.assertEqual(_conjugate({1<<(wire):Q2(1)},low(joint),4),
                                     {label:Q2(1)})
                tau=paired_t(left,right)
                self.assertEqual(ledger(tau)['T'],2)
                self.assertEqual(ledger(tau)['bands'],1)
                k,pivot=pauli_diagonalizer(probe)
                native_p=_inverse(k)+[('S',pivot),('S',pivot)]+k
                word += _inverse(tau)+native_p+tau
            actual,den=ring_matrix(low(word),4)
            wanted=ring_matrix(low(_inverse(_raw_linear_cell(1)[0])
                                   +_raw_linear_cell(direction)[0]),4)
            self.assertTrue(equal_ring((sign*actual,den),wanted))

    def test_complete_recoupler_nonzero_columns_and_actual_inverse(self):
        old=_expected_nonzero()
        columns=[j for j in range(256) if j&15]
        self.assertTrue(equal_ring((self.recoupler[0][:,:,columns],self.recoupler[1]),
                                   (old[0][:,:,columns],old[1])))
        self.assertTrue(equal_ring(self.recoupler_inverse,
                                   adjoint_ring(self.recoupler)))
        zero=[j for j in range(256) if not j&15]
        self.assertFalse(equal_ring((self.recoupler[0][:,:,zero],self.recoupler[1]),
                                    (old[0][:,:,zero],old[1])))

    def test_full_intertwining_all_addresses_columns_and_actual_inverse(self):
        self.assertTrue(equal_ring(self.decoder_inverse,adjoint_ring(self.decoder)))
        for address in range(8):
            query=[_tower_shear(3,address,j) for j in range(256)]
            translation=[_translation(3,address,j) for j in range(256)]
            self.assertTrue(np.array_equal(self.decoder[0][:,:,query],
                                           self.decoder[0][:,translation,:]))

    def test_emitted_native_ledgers_and_physical_bands(self):
        rows=((recoupler_word(),160,80,1104,1104),
              (decoder_word(),172,83,1234,1286))
        for word,t,bands,forward_gates,inverse_gates in rows:
            for actual,gates in ((word,forward_gates),(_inverse(word),inverse_gates)):
                counted=ledger(actual)
                self.assertEqual((counted['T'],counted['bands'],counted['gates'],
                                  counted['width']),(t,bands,gates,8))
                self.assertLessEqual({g[0] for g in actual},
                                     {'H','S','CX','T','TDG'})
                for band in _bands(actual):
                    self.assertEqual(len({g[1] for g in band}),len(band))


if __name__=='__main__':
    unittest.main()
