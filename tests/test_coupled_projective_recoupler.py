"""Coupled eight-address decoder with six returned arbitrary phase helpers.

The native word, literal phases, both actual inverses and every field address
are checked exactly. The projective predicate is finite; no uniform decoder
depth claim follows from this fixture.
"""
from functools import lru_cache
from itertools import product
import unittest

import numpy as np

from tests.test_projective_field_recoupler import (
    paired_t, pauli_diagonalizer, linear_word, _baseline_word, _raw_linear_cell,
    fourier_word, ring_matrix, equal_ring, adjoint_ring, ledger, mul, inv, trace,
    _inverse, _native_x, _relative_toffoli, _tower_shear, _translation)
from tests.test_joint_lookup_pauli_source import _bands
from tests.test_flat_field_bank_native import full_intertwiner


def affine_r91_batch(pairs, helpers):
    """Ftilde(a,b;h): four T targets h, only logical parity CNOT controls.

    A form is (linear_mask, constant). Ftilde differs from (-i)^(ab)Z_h^(ab)
    by scalar i iff both constants are one; these cancel in the commutator.
    """
    word=[]
    # Parity vertices 0, a, a+b, b, then return to zero.
    vertices=((0,0),(1,0),(1,1),(0,1))
    changes=(0,1,0,1)
    for step,(ia,ib) in enumerate(vertices):
        for (a,b),h in zip(pairs,helpers):
            parity=(ia+ib+ia*a[1]+ib*b[1])&1
            word.append(('TDG' if parity else 'T',h))
        for pair,h in zip(pairs,helpers):
            mask=pair[changes[step]][0]
            word += [('CX',j,h) for j in range(mask.bit_length()) if mask>>j&1]
    return word


def quartic_batch(forms, helpers):
    assert len(forms)==len(helpers) and len(set(helpers))==len(helpers)
    for row in forms:
        assert all((mask>>h)&1==0 for mask,const in row for h in helpers)
    u=affine_r91_batch([(r[0],r[1]) for r in forms],helpers)
    h=[('H',j) for j in helpers]
    v=h+affine_r91_batch([(r[2],r[3]) for r in forms],helpers)+h
    return u+v+_inverse(u)+_inverse(v)


def displacement_forms():
    out=[]
    for direction in (4,5,6,7):
        columns=[mul(3,inv(3,direction),1<<j) for j in range(4)]
        rows=[sum(((column>>i)&1)<<j for j,column in enumerate(columns)) for i in range(4)]
        shift=_raw_linear_cell(direction)[1]
        if not shift:continue
        for scalarbit in (0,1):
            vector=mul(3,1<<scalarbit,shift)
            out.append(((rows[2],1),(rows[3],1),(rows[scalarbit],0),(vector<<4,0)))
    assert len(out)==6
    return out


def displacement_batch_word():
    """XOR displacement on U0..3,V4..7; arbitrary helpers8..13 return."""
    forms=displacement_forms();helpers=list(range(8,14))
    wall=[('H',j) for j in range(4,8)]
    word=wall+quartic_batch(forms,helpers)+wall
    return word,{'forms':forms,'helpers':helpers,'ledger':ledger(word),
                 'inverse_ledger':ledger(_inverse(word))}


def evalform(form, logical):
    return ((form[0] & logical).bit_count() & 1) ^ form[1]


def helper_matrix(word, logical, helper):
    """Exact one-helper word at fixed principal input, including both columns."""
    local=[]
    for name,*wires in word:
        if name=='CX':
            c,t=wires
            assert c<8 and t>=8
            if t==helper and logical>>c&1:
                # Native X as HSSH.
                local += [('H',0),('S',0),('S',0),('H',0)]
        else:
            assert len(wires)==1 and wires[0]>=8
            if wires[0]==helper:local.append((name,0))
    return ring_matrix(local,1)


def relative_pair(left, right, probe, controls=(2, 3)):
    tau=paired_t(left,right)
    k,pivot=pauli_diagonalizer(probe)
    pword=(_inverse(k)+[('H',pivot)]
           +_relative_toffoli(*controls,pivot)+[('H',pivot)]+k)
    return _inverse(tau)+pword+tau


def selected_pair(direction, pair, negative_phase=False):
    coords=linear_word([mul(3,inv(3,direction),1<<j) for j in range(4)],range(4))
    negative=_native_x(2)+_native_x(3)
    inner=relative_pair(*pair)
    if negative_phase:inner += [('H',3),('CX',2,3),('H',3)]
    return coords+negative+inner+negative+_inverse(coords)


def common_a():
    # g=u2 OR u3. The inner word is i^g P^g, P=Z_V1.
    # (iP)^c(iP)^d(-iP)^(cd) = i^(c OR d) P^(c OR d).
    tau=paired_t(19,50)
    cliff=[('S',2),('S',3),('H',5),('CX',2,5),('CX',3,5),('H',5)]
    relative=[('H',5)]+_relative_toffoli(2,3,5)+[('H',5)]
    return _inverse(tau)+cliff+relative+tau


def coupled_linear_word():
    return (_baseline_word()
            +selected_pair(5,(156,228,64))
            +selected_pair(6,(89,247,44))
            +common_a()
            +selected_pair(4,(156,228,112))
            +selected_pair(7,(89,247,87),negative_phase=True))


def exterior_after():
    wall=[g for j in range(4) for g in [('CX',j,4+j),('CX',4+j,j)]]
    return wall+fourier_word(2,2)+fourier_word(2,6)


def recoupler_word():
    batch,_=displacement_batch_word()
    return fourier_word(3,4)+coupled_linear_word()+batch+exterior_after()


def decoder_word():
    bottom=[g for gate in full_intertwiner()
            for g in [gate,(gate[0],*(j+4 for j in gate[1:]))]]
    return recoupler_word()+bottom


def multiply_ring(left, right):
    a,da=left;b,db=right
    out=np.zeros_like(a)
    bound=max(1,int(np.max(abs(a))))*max(1,int(np.max(abs(b))))*a.shape[1]*4
    assert bound<1<<61
    for i in range(4):
        for j in range(4):
            if np.any(a[i]) and np.any(b[j]):
                out[(i+j)%4]+=(-1 if i+j>=4 else 1)*(a[i]@b[j])
    den=da+db
    while den and not np.any(out&1):out//=2;den-=1
    return out,den


def _field_cell(direction):
    """Independent trace-Fourier formula; entries are integer halves."""
    d0, d1 = direction & 3, direction >> 2
    norm = mul(2, d0, d0) ^ mul(2, d0, d1) ^ mul(2, 2, mul(2, d1, d1))
    out = np.zeros((4, 16, 16), dtype=np.int64)
    for v in range(16):
        uv = mul(3, direction, v)
        rhat = mul(3, 4, uv) >> 2
        s = mul(2, uv >> 2, inv(2, norm))
        for r in range(4):
            a = mul(3, direction, r | (s << 2))
            out[0, a, v] = (-1) ** trace(2, mul(2, rhat, r))
    return out, 1


def _field_shift(u):
    if not u:
        return 0
    d0, d1 = u & 3, u >> 2
    norm = mul(2, d0, d0) ^ mul(2, d0, d1) ^ mul(2, 2, mul(2, d1, d1))
    s = mul(2, mul(2, d0, d1), inv(2, norm))
    return mul(3, u, s << 2)


def _pauli_ring(label):
    x, z = label & 15, label >> 4
    out = np.zeros((4, 16, 16), dtype=np.int64)
    power = 2 * ((x & z).bit_count() % 4)
    for col in range(16):
        out[power % 4, col ^ x, col] = (-1) ** (
            power // 4 + (z & col).bit_count())
    return out, 0


def _times_i(value):
    a, d = value
    return np.array([-a[2], -a[3], a[0], a[1]]), d


def _reflection(left, right, probe):
    """Independent J=(P-iAP-iBP-ABP)/2 with integer coefficients."""
    a, b, p = _pauli_ring(left), _pauli_ring(right), _pauli_ring(probe)
    ap, bp = multiply_ring(a, p), multiply_ring(b, p)
    abp = multiply_ring(a, bp)
    assert ap[1] == bp[1] == abp[1] == 0
    return p[0] - _times_i(ap)[0] - _times_i(bp)[0] - abp[0], 1


def _identity(size):
    out = np.zeros((4, size, size), dtype=np.int64)
    out[0] = np.eye(size, dtype=np.int64)
    return out, 0


def _displacement_permutation():
    """Effective permutation derived from the emitted quartic forms."""
    rows = displacement_forms()
    permutation = []
    for basis in range(256):
        u, v = basis & 15, basis >> 4
        shift = 0
        for row in rows:
            assert all(mask < 16 for mask, _ in row[:3])
            assert row[3][0] & 15 == 0 and row[3][1] == 0
            if all(evalform(form, u) for form in row[:3]):
                shift ^= row[3][0] >> 4
        permutation.append(u | ((v ^ shift) << 4))
    return permutation


@lru_cache(None)
def logical_ring(decoder=False, inverse=False):
    """Compose the verified tensor-identity batch with actual native walls.

    Inverse=True evaluates each actual reversed native wall separately; it
    does not manufacture the answer by taking the forward matrix's adjoint.
    The private batch is self-inverse as a complete XOR operator, proved below
    for both actual native orientations on every private-helper sector.
    """
    prefix = fourier_word(3, 4) + coupled_linear_word()
    suffix = exterior_after()
    if decoder:
        suffix += [g for gate in full_intertwiner()
                   for g in [gate, (gate[0], *(j + 4 for j in gate[1:]))]]
    first, last = (prefix, suffix) if not inverse else (
        _inverse(suffix), _inverse(prefix))
    matrix, den = ring_matrix(first)
    moved = matrix[:, _displacement_permutation(), :]
    return multiply_ring(ring_matrix(last), (moved, den))


class CoupledProjectiveRecouplerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.linear = ring_matrix(coupled_linear_word())
        cls.linear_inverse = ring_matrix(_inverse(coupled_linear_word()))
        cls.recoupler = logical_ring()
        cls.recoupler_inverse = logical_ring(inverse=True)
        cls.decoder = logical_ring(decoder=True)
        cls.decoder_inverse = logical_ring(decoder=True, inverse=True)

    def test_affine_private_phase_keeps_its_literal_constant(self):
        forms = [(mask, constant) for mask in range(4) for constant in (0, 1)]
        for a, b in product(forms, repeat=2):
            word = affine_r91_batch([(a, b)], [2])
            wanted = np.zeros((4, 8, 8), dtype=np.int64)
            for col in range(8):
                ab = evalform(a, col) * evalform(b, col)
                # i^(ka kb) (-i)^(ab) (-1)^(ab h), including scalar i.
                power = (2*a[1]*b[1] - 2*ab + 4*ab*((col >> 2) & 1)) % 8
                wanted[power % 4, col, col] = (-1) ** (power // 4)
            expected = wanted, 0
            self.assertTrue(equal_ring(ring_matrix(word, 3), expected))
            self.assertTrue(equal_ring(ring_matrix(_inverse(word), 3),
                                       adjoint_ring(expected)))

    def test_quartic_affine_commutator_all_columns_and_actual_inverse(self):
        rows = [tuple((1 << j, c) for j, c in enumerate(constants))
                for constants in product((0, 1), repeat=4)]
        # Dependent, repeated, overlapping and constant affine forms are allowed.
        rows += [((3, 1), (3, 1), (5, 0), (0, 1)),
                 ((0, 0), (7, 1), (3, 0), (15, 1))]
        for row in rows:
            word = quartic_batch([row], [4])
            wanted = np.zeros((4, 32, 32), dtype=np.int64)
            for col in range(32):
                wanted[0, col, col] = -1 if all(evalform(f, col) for f in row) else 1
            for emitted in (word, _inverse(word)):
                self.assertTrue(equal_ring(ring_matrix(emitted, 5), (wanted, 0)))

    def test_private_batch_every_helper_sector_and_displacement(self):
        forms, helpers = displacement_forms(), list(range(8, 14))
        batch = quartic_batch(forms, helpers)
        actual_inverse = _inverse(batch)
        for logical in range(256):
            for row, helper in zip(forms, helpers):
                wanted = _identity(2)
                if all(evalform(f, logical) for f in row):
                    wanted = -wanted[0], wanted[1]
                for emitted in (batch, actual_inverse):
                    self.assertTrue(equal_ring(helper_matrix(emitted, logical, helper),
                                               wanted))
            u, v = logical & 15, logical >> 4
            actual = sum(all(evalform(f, logical) for f in row) for row in forms) % 2
            self.assertEqual(actual, (_field_shift(u) & v).bit_count() % 2)
        self.assertEqual(ledger(batch)['T'], 96)
        self.assertEqual(ledger(batch)['bands'], 16)
        self.assertEqual(len(batch), 360)
        for emitted in (batch, actual_inverse):
            for band in _bands(emitted):
                self.assertEqual(len({g[1] for g in band}), len(band))
                self.assertTrue(all(g[1] in helpers for g in band))

    def test_independent_paired_reflections_and_common_or_phase(self):
        ja = _reflection(19, 50, 32)
        jb64, jb112 = _reflection(156, 228, 64), _reflection(156, 228, 112)
        jc44, jc87 = _reflection(89, 247, 44), _reflection(89, 247, 87)
        for j in (ja, jb64, jb112, jc44, jc87):
            self.assertTrue(equal_ring(j, adjoint_ring(j)))
            self.assertTrue(equal_ring(multiply_ring(j, j), _identity(16)))
        for direction, first, last, sign in (
                (4, ja, jb112, 1), (5, jb64, ja, 1),
                (6, jc44, ja, 1), (7, ja, jc87, -1)):
            actual, den = multiply_ring(last, first)
            expected = multiply_ring(_field_cell(direction),
                                     adjoint_ring(_field_cell(1)))
            self.assertTrue(equal_ring((sign * actual, den), expected))
        common = ring_matrix(common_a())
        for u in range(16):
            selected = int(bool(u & 12))
            predicates = [int(mul(3, inv(3, d), u) >> 2 == 0) for d in (4,5,6,7)]
            self.assertEqual(sum(predicates) % 2, selected)
            if u:
                self.assertEqual((1j**selected)*((-1j)**sum(predicates)), 1)
            block = _times_i(ja) if selected else _identity(16)
            indices = np.arange(16)*16 + u
            actual = common[0][:, indices, :][:, :, indices], common[1]
            self.assertTrue(equal_ring(actual, block))

    def test_complete_linear_word_and_literal_zero_fiber(self):
        # The independent nonzero field formula also fixes every literal sign.
        for u in range(1, 16):
            indices = np.arange(16)*16 + u
            actual = self.linear[0][:, indices, :][:, :, indices], self.linear[1]
            self.assertTrue(equal_ring(actual, _field_cell(u)))
        zero = _field_cell(1)
        for triple in ((156,228,64), (89,247,44), (156,228,112), (89,247,87)):
            zero = multiply_ring(_reflection(*triple), zero)
        zero = -zero[0], zero[1]
        indices = np.arange(16)*16
        actual = self.linear[0][:, indices, :][:, :, indices], self.linear[1]
        self.assertTrue(equal_ring(actual, zero))
        self.assertTrue(equal_ring(self.linear_inverse, adjoint_ring(self.linear)))

    def test_complete_intertwining_all_columns_and_actual_inverses(self):
        for forward, inverse in ((self.recoupler, self.recoupler_inverse),
                                 (self.decoder, self.decoder_inverse)):
            self.assertTrue(equal_ring(inverse, adjoint_ring(forward)))
            self.assertTrue(equal_ring(multiply_ring(inverse, forward), _identity(256)))
        for address in range(8):
            source = [_tower_shear(3, address, col) for col in range(256)]
            target = [_translation(3, address, col) for col in range(256)]
            self.assertTrue(np.array_equal(self.decoder[0][:, :, source],
                                           self.decoder[0][:, target, :]))

    def test_native_ledgers_and_distinct_physical_t_targets(self):
        rows = ((coupled_linear_word(), 40, 20, 577, None, 8),
                (recoupler_word(), 136, 36, 980, 984, 14),
                (decoder_word(), 148, 39, 1110, 1166, 14))
        for word, t, bands, forward_gates, inverse_gates, width in rows:
            for emitted, gates in ((word, forward_gates), (_inverse(word), inverse_gates)):
                counted = ledger(emitted)
                self.assertEqual((counted['T'], counted['bands'], counted['width']),
                                 (t, bands, width))
                if gates is not None:
                    self.assertEqual(counted['gates'], gates)
                self.assertLessEqual({g[0] for g in emitted}, {'H','S','CX','T','TDG'})
                for band in _bands(emitted):
                    self.assertEqual(len({g[1] for g in band}), len(band))


if __name__ == '__main__':
    unittest.main()
