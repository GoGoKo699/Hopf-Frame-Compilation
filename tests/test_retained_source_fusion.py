"""Exact bounded Laurent-ring checks of retained-source logical fusion.

Coefficients are Gaussian dyadics in Z[i,1/2][z]/(z**q-1). Every int64
operation has a conservative pre-operation bound below 2**60; normalization
is exact. The source tensor is shared polynomial multiplication, not two
independent source registers. Logical transpose does not invert z.

Independent plane words test the four/eight-mode formulas. Native Clifford
matrices and source-character/physical-bank evaluations are separate small
numerical checks. These are not a scalable selector or native group emitter.
"""
from __future__ import annotations

import itertools
import unittest

import numpy as np

try:
    from .test_eight_mode_coupling import _magic_word
    from .test_operator_source_compiler import _word_matrix
except ImportError:
    from test_eight_mode_coupling import _magic_word
    from test_operator_source_compiler import _word_matrix


ATOL = 2e-11
X = np.array([[0, 1], [1, 0]], complex)
Z = np.diag([1, -1]).astype(complex)
MAGIC_NUMERATOR = np.array([[1, 0, 1j, 0], [0, 1j, 0, 1],
                            [0, 1j, 0, -1], [1, 0, -1j, 0]], complex)


class _Ring:
    """Matrix over Gaussian dyadic Laurent coefficients, source exponent first."""

    def __init__(self, real, imag=None, denominator=0):
        self.real = np.array(real, dtype=np.int64, copy=True)
        self.imag = (np.zeros_like(self.real) if imag is None
                     else np.array(imag, dtype=np.int64, copy=True))
        assert self.real.shape == self.imag.shape and self.real.ndim == 3
        self.denominator = denominator
        self.q, self.rows, self.cols = self.real.shape
        while self.denominator and not np.any(self.real % 2) and not np.any(self.imag % 2):
            self.real //= 2
            self.imag //= 2
            self.denominator -= 1
        assert self.bound() < 1 << 60

    def bound(self):
        return max(int(np.max(np.abs(self.real))), int(np.max(np.abs(self.imag))))

    def half(self):
        return _Ring(self.real, self.imag, self.denominator + 1)

    def __neg__(self):
        return _Ring(-self.real, -self.imag, self.denominator)

    def __add__(self, other):
        assert self.real.shape == other.real.shape
        power = max(self.denominator, other.denominator)
        first, second = power - self.denominator, power - other.denominator
        assert (self.bound() << first) + (other.bound() << second) < 1 << 60
        return _Ring((self.real << first) + (other.real << second),
                     (self.imag << first) + (other.imag << second), power)

    def __sub__(self, other):
        return self + -other

    def __matmul__(self, other):
        assert self.q == other.q and self.cols == other.rows
        # Each output coefficient sums at most q*cols complex products.
        assert 2 * self.q * self.cols * self.bound() * other.bound() < 1 << 60
        real = np.zeros((self.q, self.rows, other.cols), dtype=np.int64)
        imag = np.zeros_like(real)
        for a in range(self.q):
            if not np.any(self.real[a]) and not np.any(self.imag[a]):
                continue
            for b in range(self.q):
                k = (a + b) % self.q
                real[k] += self.real[a] @ other.real[b] - self.imag[a] @ other.imag[b]
                imag[k] += self.real[a] @ other.imag[b] + self.imag[a] @ other.real[b]
        return _Ring(real, imag, self.denominator + other.denominator)

    def tensor(self, other):
        assert self.q == other.q
        assert 2 * self.q * self.bound() * other.bound() < 1 << 60
        real = np.zeros((self.q, self.rows * other.rows, self.cols * other.cols), dtype=np.int64)
        imag = np.zeros_like(real)
        for a in range(self.q):
            if not np.any(self.real[a]) and not np.any(self.imag[a]):
                continue
            for b in range(self.q):
                k = (a + b) % self.q
                real[k] += np.kron(self.real[a], other.real[b]) - np.kron(self.imag[a], other.imag[b])
                imag[k] += np.kron(self.real[a], other.imag[b]) + np.kron(self.imag[a], other.real[b])
        return _Ring(real, imag, self.denominator + other.denominator)

    def transpose(self):
        return _Ring(self.real.transpose(0, 2, 1), self.imag.transpose(0, 2, 1), self.denominator)

    def dagger(self):
        inverse = (-np.arange(self.q)) % self.q
        return _Ring(self.real[inverse].transpose(0, 2, 1),
                     -self.imag[inverse].transpose(0, 2, 1), self.denominator)

    def entry(self, row, col):
        return _Ring(self.real[:, row:row + 1, col:col + 1],
                     self.imag[:, row:row + 1, col:col + 1], self.denominator)

    def evaluate(self, character):
        weights = character ** np.arange(self.q)
        return np.einsum('k,kij->ij', weights, self.real + 1j * self.imag) / (1 << self.denominator)

    def on_source(self, shift):
        result = np.zeros((self.rows * len(shift), self.cols * len(shift)), complex)
        power = np.eye(len(shift), dtype=complex)
        for k in range(self.q):
            result += np.kron(self.real[k] + 1j * self.imag[k], power)
            power = power @ shift
        return result / (1 << self.denominator)

    def __eq__(self, other):
        difference = self - other
        return not np.any(difference.real) and not np.any(difference.imag)


def _constant(q, matrix):
    matrix = np.asarray(matrix, dtype=complex)
    assert np.array_equal(matrix.real, np.rint(matrix.real))
    assert np.array_equal(matrix.imag, np.rint(matrix.imag))
    real = np.zeros((q, *matrix.shape), dtype=np.int64)
    imag = np.zeros_like(real)
    real[0], imag[0] = matrix.real.astype(np.int64), matrix.imag.astype(np.int64)
    return _Ring(real, imag)


def _identity(q, size):
    return _constant(q, np.eye(size))


def _power(q, exponent, size=1):
    real = np.zeros((q, size, size), dtype=np.int64)
    real[exponent % q] = np.eye(size, dtype=np.int64)
    return _Ring(real)


def _block_diagonal(*blocks):
    q, total = blocks[0].q, sum(block.rows for block in blocks)
    power = max(block.denominator for block in blocks)
    real, imag = np.zeros((q, total, total), np.int64), np.zeros((q, total, total), np.int64)
    offset = 0
    for block in blocks:
        assert block.rows == block.cols and block.q == q
        scale = power - block.denominator
        assert block.bound() << scale < 1 << 60
        real[:, offset:offset + block.rows, offset:offset + block.cols] = block.real << scale
        imag[:, offset:offset + block.rows, offset:offset + block.cols] = block.imag << scale
        offset += block.rows
    return _Ring(real, imag, power)


def _plane(q, size, low, high, angle):
    cosine = (_power(q, angle) + _power(q, -angle)).half()
    sine = (_constant(q, [[1j]]) @ (_power(q, angle) - _power(q, -angle))).half()
    result = _identity(q, size)
    for row, col, value in ((low, low, cosine - _identity(q, 1)),
                            (high, high, cosine - _identity(q, 1)),
                            (low, high, -sine), (high, low, sine)):
        unit = np.zeros((size, size), complex)
        unit[row, col] = 1
        result = result + _constant(q, unit).tensor(value)
    return result


def _four(q, labels):
    a, b, c = labels
    return _plane(q, 4, 2, 3, c) @ _plane(q, 4, 0, 1, b) @ _plane(q, 4, 0, 2, a)


def _factors(q, labels):
    a, b, c = labels
    ix = _identity(q, 2)
    xp, xm = (ix + _constant(q, X)).half(), (ix - _constant(q, X)).half()
    zp, zm = (ix + _constant(q, Z)).half(), (ix - _constant(q, Z)).half()
    first = (xp + _power(q, b + c, 2) @ xm) @ (zp + _power(q, a, 2) @ zm)
    second = (_power(q, -b, 2) @ xp + _power(q, -c, 2) @ xm) @ (_power(q, -a, 2) @ zp + zm)
    return first, second


def _magic_conjugate(word):
    numerator = _constant(word.q, MAGIC_NUMERATOR)
    if word.rows == 8:
        numerator = _identity(word.q, 2).tensor(numerator)
    return (numerator @ word @ numerator.dagger()).half()


def _root_bell(q, angle):
    # Direct plane conjugation is deliberately not used in this expression.
    bell = _constant(q, [[1], [0], [0], [1]])
    projector = (bell @ bell.dagger()).half()
    cosine = (_power(q, angle) + _power(q, -angle)).half()
    sine = (_constant(q, [[1j]]) @ (_power(q, angle) - _power(q, -angle))).half()
    return (_identity(q, 8)
            + (_identity(q, 2).tensor(projector)).tensor(cosine - _identity(q, 1))
            + _constant(q, [[0, -1], [1, 0]]).tensor(projector).tensor(sine))


class RetainedSourceFusionTests(unittest.TestCase):
    def test_all_q8_four_mode_label_triples_match_exact_shared_ring_factors(self):
        q = 8
        identity = _identity(q, 4)
        for labels in itertools.product(range(q), repeat=3):
            original = _four(q, labels)
            first, second = _factors(q, labels)
            self.assertEqual(_magic_conjugate(original), first.tensor(second), labels)
            self.assertEqual(original.dagger() @ original, identity, labels)

    def test_all_characters_and_q4_full_physical_source_bank(self):
        native_magic = _word_matrix(2, _magic_word())
        np.testing.assert_allclose(native_magic, MAGIC_NUMERATOR / np.sqrt(2), atol=ATOL, rtol=0)
        for labels in itertools.product(range(8), repeat=3):
            ring = _four(8, labels)
            for character in range(8):
                expected = np.eye(4)
                for (low, high), label in zip(((0, 2), (0, 1), (2, 3)), labels):
                    theta = 2 * np.pi * character * label / 8
                    gate = np.eye(4)
                    gate[np.ix_([low, high], [low, high])] = [[np.cos(theta), -np.sin(theta)],
                                                            [np.sin(theta), np.cos(theta)]]
                    expected = gate @ expected
                np.testing.assert_allclose(ring.evaluate(np.exp(-2j * np.pi * character / 8)),
                                           expected, atol=ATOL, rtol=0)
        shift = np.zeros((16, 16), complex)
        for basis in range(16):
            shift[((basis << 1) & 15) | (basis >> 3), basis] = 1
        magic = np.kron(native_magic, np.eye(16))
        for labels in ((1, 2, 3), (3, 1, 2), (1, 3, 0)):
            first, second = _factors(4, labels)
            actual = first.tensor(second).on_source(shift)
            expected = magic @ _four(4, labels).on_source(shift) @ magic.conj().T
            np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)

    def test_determinant_gauge_literal_source_phase_and_inactive_sandwich(self):
        q = 8
        for labels in ((1, 2, 3), (3, 1, 5), (7, 4, 2)):
            first, second = _factors(q, labels)
            det_first = first.entry(0, 0) @ first.entry(1, 1) - first.entry(0, 1) @ first.entry(1, 0)
            det_second = second.entry(0, 0) @ second.entry(1, 1) - second.entry(0, 1) @ second.entry(1, 0)
            self.assertEqual(det_first, _power(q, sum(labels)))
            self.assertEqual(det_second, _power(q, -sum(labels)))
            self.assertEqual(first.dagger() @ first, _identity(q, 2))
            self.assertEqual(second.dagger() @ second, _identity(q, 2))
            # Dropping the z^-b scalar from B preserves a spurious source phase.
            wrong = first.tensor(_power(q, labels[1], 2) @ second)
            correct = first.tensor(second)
            self.assertEqual(wrong, _power(q, labels[1], 4) @ correct)
            self.assertNotEqual(wrong, correct)
        # On h=0 all selected shifts are I; the fixed Clifford boundaries
        # still cancel on the full source space. No private-work emitter is modeled.
        disabled = _factors(q, (0, 0, 0))
        active = _factors(q, (1, 3, 2))
        blocks = _block_diagonal(disabled[0].tensor(disabled[1]), active[0].tensor(active[1]))
        numerator = _identity(q, 2).tensor(_constant(q, MAGIC_NUMERATOR))
        returned = (numerator.dagger() @ blocks @ numerator).half()
        self.assertEqual(returned, _block_diagonal(_identity(q, 4), _four(q, (1, 3, 2))))

    def test_eight_modes_distinct_children_and_integer_root_exponents(self):
        q, root = 8, 1
        left, right = (1, 2, 5), (3, 6, 1)
        original = _block_diagonal(_four(q, left), _four(q, right)) @ _plane(q, 8, 0, 4, root)
        a0, b0 = _factors(q, left)
        a1, b1 = _factors(q, right)
        children = _block_diagonal(a0.tensor(b0), a1.tensor(b1))
        bell_root = _root_bell(q, root)
        self.assertEqual(_magic_conjugate(original), children @ bell_root)
        self.assertEqual(original.dagger() @ original, _identity(q, 8))
        # D=SH and M have sqrt(2) denominators; D tensor M is dyadic.
        basis = _constant(q, np.kron([[1, 1], [1j, -1j]], MAGIC_NUMERATOR)).half()
        diagonal = _block_diagonal(*[_power(q, exponent) for exponent in (root, 0, 0, 0, -root, 0, 0, 0)])
        self.assertEqual(basis @ diagonal @ basis.dagger(), bell_root)
        # Four quarter-shifts cannot be interpreted as integer powers modulo8
        # at odd root label; the joint stabilizer exponent avoids this problem.
        self.assertFalse(any(4 * exponent % q == root for exponent in range(q)))
        self.assertNotEqual(_root_bell(q, 4 * root), bell_root)

    def test_fixed_basis_single_diagonal_round_does_not_cover_root_child_family(self):
        q = 8
        root = _root_bell(q, 1)
        child = _magic_conjugate(_plane(q, 8, 0, 2, 1))
        self.assertNotEqual(root @ child, child @ root)
        # A common fixed conjugation of one diagonal shift round would commute.
        # This only restricts that family of exact rewrites, not general depth.

    def test_bell_stabilizer_extraction_uses_logical_transpose_and_keeps_full_factor(self):
        q = 8
        branches = ((1, 2, 5), (3, 6, 1))
        full, local, stabilizer = [], [], []
        bell = _constant(q, [[1], [0], [0], [1]])
        for labels in branches:
            first, second = _factors(q, labels)
            value = first @ second.transpose()
            keep = second.dagger().transpose().tensor(second)
            complete = first.tensor(second)
            moved = value.tensor(_identity(q, 2))
            self.assertEqual(moved @ keep, complete)
            self.assertEqual(keep @ bell, bell)
            self.assertEqual(keep.dagger() @ keep, _identity(q, 4))
            self.assertEqual(value.dagger() @ value, _identity(q, 2))
            # Transposing the physical source shift as well would invert z.
            # The required transpose acts only on the two logical indices.
            logical_transpose = second.transpose()
            inverse = (-np.arange(q)) % q
            physical_transpose = _Ring(logical_transpose.real[inverse],
                                       logical_transpose.imag[inverse],
                                       logical_transpose.denominator)
            wrong = (first @ physical_transpose).tensor(_identity(q, 2)) @ keep
            self.assertNotEqual(wrong, complete)
            full.append(complete)
            local.append(moved)
            stabilizer.append(keep)
        complete, moved, keep = map(lambda blocks: _block_diagonal(*blocks), (full, local, stabilizer))
        root = _root_bell(q, 1)
        self.assertEqual(keep @ root, root @ keep)
        self.assertEqual(complete @ root, moved @ root @ keep)
        self.assertEqual(complete @ root @ complete.dagger(), moved @ root @ moved.dagger())
        self.assertNotEqual(complete @ root, moved @ root)  # K cannot be dropped from the full word.
        # The transported boundary can already require noncommuting SU(2)
        # factors. Its 4x4 word describes the Bell input sector only.
        first0, second0 = _factors(q, (0, 1, 0))
        first1, second1 = _factors(q, (1, 0, 0))
        value0, value1 = first0 @ second0.transpose(), first1 @ second1.transpose()
        identity2 = _identity(q, 2)
        xp, xm = (identity2 + _constant(q, X)).half(), (identity2 - _constant(q, X)).half()
        zp, zm = (identity2 + _constant(q, Z)).half(), (identity2 - _constant(q, Z)).half()
        self.assertEqual(value0, _power(q, -1, 2) @ xp + _power(q, 1, 2) @ xm)
        self.assertEqual(value1, _power(q, -1, 2) @ zp + _power(q, 1, 2) @ zm)
        commutator = value0 @ value1 - value1 @ value0
        self.assertNotEqual(commutator, _constant(q, np.zeros((2, 2))))
        numerical = commutator.evaluate(np.exp(-2j * np.pi / q))
        np.testing.assert_allclose(numerical, [[0, 1], [-1, 0]], atol=ATOL, rtol=0)  # iY.
        self.assertAlmostEqual(np.linalg.norm(numerical, 2), 1, places=11)
        transport = _block_diagonal(value0, value1) @ _plane(q, 2, 0, 1, 1).tensor(identity2)
        self.assertEqual(transport.dagger() @ transport, _identity(q, 4))
        children = _block_diagonal(first0.tensor(second0), first1.tensor(second1))
        embedding = identity2.tensor(bell)
        self.assertEqual(children @ root @ embedding, transport.tensor(identity2) @ embedding)
        self.assertNotEqual(children @ root, transport.tensor(identity2))

    def test_orthogonal_four_target_sectors_batch_stabilizers_in_two_shift_rounds(self):
        # Two target pairs: Q0 has upper!=00,lower=00; Q1,u has
        # upper=u,lower!=00. These five rank-three sectors plus vacuum
        # partition all sixteen logical modes in the ORIGINAL basis.
        q = 8
        identity4, identity16 = _identity(q, 4), _identity(q, 16)
        supports = ((4, 8, 12),) + tuple(tuple(4 * u + v for v in (1, 2, 3)) for u in range(4))
        projectors = [_constant(q, np.diag([int(basis in support) for basis in range(16)]))
                      for support in supports]
        vacuum = _constant(q, np.diag([1] + [0] * 15))
        total = vacuum
        for projector in projectors:
            total = total + projector
        self.assertEqual(total, identity16)
        zero = _constant(q, np.zeros((16, 16)))
        for i, first in enumerate(projectors):
            for j, second in enumerate(projectors):
                self.assertEqual(first @ second, first if i == j else zero)

        pair_magic = _constant(q, MAGIC_NUMERATOR)
        magic_all = _constant(q, np.kron(MAGIC_NUMERATOR, MAGIC_NUMERATOR)).half()
        h2 = np.array([[1, 1], [1, -1]])
        hadamard_all = _constant(q, np.kron(np.kron(h2, h2), np.kron(h2, h2))).half().half()
        x_basis = hadamard_all @ magic_all
        identity2 = _identity(q, 2)
        xp, xm = (identity2 + _constant(q, X)).half(), (identity2 - _constant(q, X)).half()
        zp, zm = (identity2 + _constant(q, Z)).half(), (identity2 - _constant(q, Z)).half()
        labels = ((1, 2, 5), (3, 6, 1), (5, 3, 7), (7, 4, 2), (2, 1, 6))
        serial, batched_z, batched_x = identity16, vacuum, vacuum
        transient_z, transient_x = vacuum, vacuum
        transports = []
        for row, ((a, b, c), projector) in enumerate(zip(labels, projectors)):
            bz = _power(q, -a, 2) @ zp + zm
            bx = _power(q, -b, 2) @ xp + _power(q, -c, 2) @ xm
            z_axis, x_axis = bz.dagger().tensor(bz), bx.dagger().tensor(bx)
            # This independent serial row uses the defining B factors,
            # not the selected exponent tables used by the batched word.
            local = (pair_magic.dagger() @ x_axis @ z_axis @ pair_magic).half()
            first_factor, second_factor = _factors(q, (a, b, c))
            value = first_factor @ second_factor.transpose()
            transport = (pair_magic.dagger() @ value.tensor(identity2) @ pair_magic).half()
            transports.append(transport)
            lifted = local.tensor(identity4) if row == 0 else identity4.tensor(local)
            row_completion = identity16 + (lifted - identity16) @ projector
            if row == 0:
                top_completion = row_completion
            serial = row_completion @ serial
            for axis in (z_axis, x_axis):
                original_axis = (pair_magic.dagger() @ axis @ pair_magic).half()
                lifted_axis = (original_axis.tensor(identity4) if row == 0
                               else identity4.tensor(original_axis))
                controlled_axis = identity16 + (lifted_axis - identity16) @ projector
                for sector in projectors + [vacuum]:
                    self.assertEqual(controlled_axis @ sector, sector @ controlled_axis)

            phases_z = _block_diagonal(*[_power(q, exponent) for exponent in (0, a, -a, 0)])
            phases_x = _block_diagonal(*[_power(q, exponent) for exponent in (0, b - c, c - b, 0)])
            diagonal_z = phases_z.tensor(identity4) if row == 0 else identity4.tensor(phases_z)
            diagonal_x = phases_x.tensor(identity4) if row == 0 else identity4.tensor(phases_x)
            # Rightmost Q acts BEFORE the common basis change: it models a
            # cached sector label, not a predicate recomputed in that basis.
            batched_z = batched_z + magic_all.dagger() @ diagonal_z @ magic_all @ projector
            batched_x = batched_x + x_basis.dagger() @ diagonal_x @ x_basis @ projector
            transient_z = transient_z + diagonal_z @ projector
            transient_x = transient_x + diagonal_x @ projector
        batched = batched_x @ batched_z
        self.assertEqual(serial, batched)
        self.assertEqual(batched.dagger() @ serial, identity16)
        self.assertEqual(serial.dagger(), batched.dagger())
        # Independent complete sixteen-mode frame: the upper four-mode
        # word is enabled ONLY when the lower pair is 00, followed by four
        # distinct lower-pair frames. The transported factors keep that guard.
        lower_zero = _constant(q, np.diag([1, 0, 0, 0]))
        top_guard = identity4.tensor(lower_zero)
        top_frame = identity16 + (_four(q, labels[0]).tensor(identity4) - identity16) @ top_guard
        bottom_frame = _block_diagonal(*[_four(q, label) for label in labels[1:]])
        frame = bottom_frame @ top_frame
        top_transport = identity16 + (transports[0].tensor(identity4) - identity16) @ top_guard
        bottom_transport = _block_diagonal(*transports[1:])
        self.assertEqual(frame, bottom_transport @ top_transport @ serial)
        self.assertEqual(frame, bottom_transport @ top_transport @ batched)
        # With identity lower frames, removing the top suffix guard can
        # preserve the first column while corrupting complementary columns.
        unguarded = transports[0].tensor(identity4) @ top_completion
        guarded = top_transport @ top_completion
        self.assertEqual(guarded, top_frame)
        first_column = _constant(q, np.eye(16)[:, :1])
        self.assertEqual(unguarded @ first_column, guarded @ first_column)
        complement = _constant(q, np.eye(16)[:, 1:2])  # |00>|01>.
        self.assertNotEqual(unguarded @ complement, guarded @ complement)
        wrong = (x_basis.dagger() @ transient_x @ x_basis
                 @ magic_all.dagger() @ transient_z @ magic_all)
        self.assertNotEqual(wrong, serial)
        # Arbitrary-cache inactive identity and a scalable native selector
        # implementation are analytic obligations, not simulated here.


if __name__ == '__main__':
    unittest.main()
