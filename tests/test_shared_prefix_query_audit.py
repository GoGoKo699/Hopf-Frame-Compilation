"""Bounded native audits of shared-prefix queries and their dirty offsets.

Every word here is emitted in Clifford+T on at most eight wires. Expected
actions are independent signed Boolean basis permutations, including every
dirty input and literal phase. These are small query-interface counterexamples
and reference circuits, not complete Hopf groups or scalable lookup emitters.
The clean-baseline contrast explicitly assumes its extra initialized bit.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis, _word_matrix,
    )
    from .test_t_depth import _inverse, _native_schedule, _word as _schedule_word
except ImportError:
    from test_operator_source_compiler import (
        _adjoint, _apply_native_word, _expand_toffolis, _word_matrix,
    )
    from test_t_depth import _inverse, _native_schedule, _word as _schedule_word


ATOL = 3e-11


def _bit(basis, wire):
    return (basis >> wire) & 1


def _indicator(address, outputs):
    """Exact two-output dirty indicator (1+a,a), using only Cliffords."""
    return [('X', outputs[0]), ('CX', address, outputs[0]),
            ('CX', address, outputs[1])]


def _short_echo(matrix):
    # x=0,t=1,Y=(2,3),X=(4,5),z=6. No indicator is assumed one-hot.
    middle = _expand_toffolis([
        ('CCX', 2 + i, 4 + j, 6)
        for i in range(2) for j in range(2) if matrix[i][j]
    ])
    local = _indicator(1, (4, 5))
    return middle + local + _adjoint(middle) + _adjoint(local)


def _expected(width, action):
    result = np.zeros((1 << width, 1 << width), dtype=complex)
    for basis in range(1 << width):
        result[action(basis), basis] = 1
    return result


def _two_group_action(basis):
    # The third data bit s can be an arbitrary retained source input.
    x, t = _bit(basis, 0), _bit(basis, 1)
    return basis ^ (x << 1) ^ ((x & (t ^ x)) << 2)


def _per_group_reference():
    # x,t,s,p1,p2. Each query completes while its own address is unchanged.
    q1, k1 = [('CX', 0, 3)], [('CX', 3, 1)]
    q2 = _expand_toffolis([('CCX', 0, 1, 4)])
    k2 = [('CX', 4, 2)]
    return q1 + k1 + _adjoint(q1) + _adjoint(k1) + q2 + k2 + _adjoint(q2) + _adjoint(k2)


def _dirty_c3x(controls, target, dirty):
    a, b, c = controls
    return _expand_toffolis([
        ('CCX', a, b, dirty), ('CCX', dirty, c, target),
        ('CCX', a, b, dirty), ('CCX', dirty, c, target),
    ])


def _query(reader, feature):
    return reader + feature + _adjoint(reader) + _adjoint(feature)


def _signed_expected(width, action):
    result = np.zeros((1 << width, 1 << width), dtype=complex)
    for basis in range(1 << width):
        image, sign = action(basis)
        result[image, basis] = sign
    return result


class SharedPrefixQueryAuditTests(unittest.TestCase):
    def assert_native_action(self, width, word, action):
        actual = _word_matrix(width, word)
        expected = _expected(width, action)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        return actual, expected

    def test_retained_prefix_short_echo_has_exact_dirty_residue_and_inverse(self):
        prefix = _indicator(0, (2, 3))
        for matrix in (((1, 0), (0, 1)), ((1, 1), (0, 1))):
            word = prefix + _short_echo(matrix) + _adjoint(prefix)

            def action(basis):
                x, t = _bit(basis, 0), _bit(basis, 1)
                residue = sum(_bit(basis, 2 + i) * matrix[i][t]
                              for i in range(2)) % 2
                return basis ^ ((matrix[x][t] ^ residue) << 6)

            actual, _ = self.assert_native_action(7, word, action)
            # All dirty banks return, but the output includes d^T D e_t.
            np.testing.assert_allclose(
                _apply_native_word(7, _adjoint(word), actual),
                np.eye(128), atol=ATOL, rtol=0,
            )

    def test_stale_baseline_after_changed_address_fails_with_cache_returned(self):
        matrix = ((1, 0), (0, 1))
        short, prefix = _short_echo(matrix), _indicator(0, (2, 3))
        word = short + prefix + [('X', 1)] + short + _adjoint(prefix)

        def action(basis):
            x, old = _bit(basis, 0), _bit(basis, 1)
            new = old ^ 1
            value = (_bit(basis, 2 + old) ^ _bit(basis, 2 + new)
                     ^ int(x == new))
            return basis ^ (1 << 1) ^ (value << 6)

        actual, _ = self.assert_native_action(7, word, action)
        # x=1,t_old=0,d=(1,0), arbitrary X=(1,1): the desired new row
        # equals one, but the stale baseline cancels it to zero.
        witness = 1 | (1 << 2) | (1 << 4) | (1 << 5)
        image = action(witness)
        self.assertEqual(_bit(image, 1), 1)
        self.assertEqual(_bit(image, 6), 0)
        self.assertEqual((image >> 2) & 15, (witness >> 2) & 15)
        self.assertAlmostEqual(actual[image, witness].real, 1.0, places=11)
        desired = witness ^ (1 << 1) ^ (1 << 6)
        self.assertNotEqual(image, desired)

    def test_completed_per_group_queries_use_changed_target_and_preserve_dirty_programs(self):
        word = _per_group_reference()
        actual, expected = self.assert_native_action(5, word, _two_group_action)
        # Equivalent local offset cancellation: the second primitive is a
        # program-controlled CNOT on the already changed target t.
        q1, k1 = [('CX', 0, 3)], [('CX', 3, 1)]
        q2 = [('CX', 0, 4)]
        k2 = _expand_toffolis([('CCX', 4, 1, 2)])
        local_echoes = (q1 + k1 + _adjoint(q1) + _adjoint(k1)
                        + q2 + k2 + _adjoint(q2) + _adjoint(k2))
        self.assert_native_action(5, local_echoes, _two_group_action)
        # Every input column was checked, hence arbitrary reference extensions
        # follow; this also explicitly retains two coherent reference columns.
        reference = np.arange(1, 65).reshape(32, 2).astype(complex)
        reference[:, 1] *= 1j
        reference /= np.linalg.norm(reference)
        np.testing.assert_allclose(actual @ reference, expected @ reference,
                                   atol=ATOL, rtol=0)

    def test_delayed_whole_word_echo_returns_program_but_changes_retained_source(self):
        load = [('CX', 0, 3), ('CX', 0, 4)]
        body = _expand_toffolis([('CX', 3, 1), ('CCX', 4, 1, 2)])
        word = load + body + _adjoint(load) + _adjoint(body)

        def action(basis):
            # Exact excess source X^{x*d1} from noncommuting body factors.
            return _two_group_action(basis) ^ ((_bit(basis, 0) & _bit(basis, 3)) << 2)

        actual, _ = self.assert_native_action(5, word, action)
        reference = _expected(5, _two_group_action)
        self.assertGreater(np.linalg.norm(actual - reference, 2), 1.9)
        witness = 1 | (1 << 3)  # x=1,t=s=0,d1=1,d2=0.
        self.assertEqual(_bit(action(witness), 1), 1)
        self.assertEqual(_bit(action(witness), 2), 0)
        self.assertEqual(action(witness) >> 3, witness >> 3)
        self.assertEqual(_bit(_two_group_action(witness), 2), 1)

    def test_separate_baseline_works_only_with_its_stated_clean_promise(self):
        # x,t,s,D,B. D=d is arbitrary; B=b is explicitly NOT silently clean.
        prepare = [('CX', 3, 4), ('CX', 0, 3)]
        body = _expand_toffolis([
            ('CX', 3, 1), ('CX', 4, 1),
            ('CCX', 3, 1, 2), ('CCX', 4, 1, 2),
        ])
        word = prepare + body + _adjoint(prepare)

        def action(basis):
            effective = _bit(basis, 0) ^ _bit(basis, 4)
            return basis ^ (effective << 1) ^ ((effective & (_bit(basis, 1) ^ effective)) << 2)

        actual, _ = self.assert_native_action(5, word, action)
        clean_columns = [basis for basis in range(32) if not _bit(basis, 4)]
        np.testing.assert_allclose(actual[:, clean_columns],
                                   _expected(5, _two_group_action)[:, clean_columns],
                                   atol=ATOL, rtol=0)
        witness = 1 << 4  # x=t=s=d=0,b=1: both data targets wrongly flip.
        self.assertEqual(action(witness), witness | (1 << 1) | (1 << 2))
        self.assertEqual(action(witness) >> 3, witness >> 3)

    def test_charged_fusion_preserves_literal_rotation_with_changed_second_address(self):
        # x0,x1,t,s,p,Y,d. A stores one nonlinear feature, not a full
        # four-bit indicator: the same shear algebra fits seven wires.
        feature = _expand_toffolis([('CCX', 0, 1, 5)])
        first = [('CX', 5, 4)]
        second = _expand_toffolis([('CCX', 5, 2, 4)])
        q1, q2 = _query(first, feature), _query(second, feature)
        # Controlled XZ = controlled R(pi/2), including its minus sign.
        body1 = [('H', 2), ('CX', 4, 2), ('H', 2), ('CX', 4, 2)]
        body2 = [('CX', 4, 3)]
        correction = [('X', 2)] + _dirty_c3x((0, 1, 2), 4, 6) + [('X', 2)]
        original = q1 + body1 + _adjoint(q1) + q2 + body2 + _adjoint(q2)
        opening = first + feature + _adjoint(first)
        closing = second + _adjoint(feature) + _adjoint(second)
        fused = opening + body1 + correction + body2 + closing

        def action(basis):
            e = _bit(basis, 0) & _bit(basis, 1)
            p, old = _bit(basis, 4), _bit(basis, 2)
            loaded = p ^ e
            new = old ^ loaded
            image = basis ^ (loaded << 2) ^ ((p ^ (e & new)) << 3)
            return image, (-1) ** (loaded & old)

        expected = _signed_expected(7, action)
        self.assertTrue(np.any(expected == -1))
        np.testing.assert_allclose(_word_matrix(7, original), expected, atol=ATOL, rtol=0)
        actual = _word_matrix(7, fused)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_apply_native_word(7, _adjoint(fused), actual),
                                   np.eye(128), atol=ATOL, rtol=0)
        for bad in (opening + body1 + body2 + closing,
                    opening + correction + body1 + body2 + closing):
            self.assertGreater(np.max(np.abs(_word_matrix(7, bad) - expected)), .9)

    def test_charged_fusion_reuses_one_activity_flag_with_actual_transition(self):
        # x0,x1,t,s,p,Y,d,h: no wire stores an old copy of h.
        feature = _expand_toffolis([('CCX', 0, 1, 5)])
        first = _expand_toffolis([('CCX', 7, 5, 4)])
        second = _dirty_c3x((7, 5, 2), 4, 6)
        q1, q2 = _query(first, feature), _query(second, feature)
        rotate = _expand_toffolis([('CCX', 7, 4, 2)])
        body1 = [('H', 2)] + rotate + [('H', 2)] + rotate
        body2 = _expand_toffolis([('CCX', 7, 4, 3)])
        transition = [('X', 7)]  # gamma=1; tests both 0->1 and 1->0.
        delta_reader = _adjoint(transition) + _adjoint(first) + transition + second
        correction = _query(delta_reader, feature)

        def correction_action(basis):
            e = _bit(basis, 0) & _bit(basis, 1)
            new_h, t = _bit(basis, 7), _bit(basis, 2)
            delta = e & ((new_h & (1 ^ t)) ^ 1)
            return basis ^ (delta << 4)

        self.assert_native_action(8, correction, correction_action)
        original = (q1 + body1 + _adjoint(q1) + transition
                    + q2 + body2 + _adjoint(q2))
        opening = first + feature + _adjoint(first)
        closing = second + _adjoint(feature) + _adjoint(second)
        fused = opening + body1 + transition + correction + body2 + closing

        def action(basis):
            e = _bit(basis, 0) & _bit(basis, 1)
            h, p, old = _bit(basis, 7), _bit(basis, 4), _bit(basis, 2)
            rotate_target = h & (p ^ (h & e))
            new = old ^ rotate_target
            new_h = h ^ 1
            flip_source = new_h & (p ^ (new_h & e & new))
            image = basis ^ (rotate_target << 2) ^ (flip_source << 3) ^ (1 << 7)
            return image, (-1) ** (rotate_target & old)

        expected = _signed_expected(8, action)
        np.testing.assert_allclose(_word_matrix(8, original), expected, atol=ATOL, rtol=0)
        actual = _word_matrix(8, fused)
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(_apply_native_word(8, _adjoint(fused), actual),
                                   np.eye(256), atol=ATOL, rtol=0)
        # Omitting gamma*f1 uses the new h as if it were the old h.
        wrong_correction = correction + _expand_toffolis([('CCX', 0, 1, 4)])
        wrong = opening + body1 + transition + wrong_correction + body2 + closing
        self.assertGreater(np.max(np.abs(_word_matrix(8, wrong) - expected)), .9)

    def test_affine_parity_refresh_has_two_literal_toffolis_and_eight_layers(self):
        for count, constant in ((1, 0), (2, 1), (4, 0), (4, 1)):
            # Prefix inputs, local t, program p, and one returned dirty d.
            t, p, dirty = count, count + 1, count + 2
            width = count + 3
            parity = ([('X', dirty)] if constant else []) + [
                ('CX', q, dirty) for q in range(count)
            ]
            middle = _native_schedule([('CCX', dirty, t, p)])
            schedule = ([('C', parity)] + middle + [('C', _adjoint(parity))]
                        + _inverse(middle))
            word = _schedule_word(schedule)
            self.assertEqual(sum(kind == 'T' for kind, _ in schedule), 8)
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in word), 14)
            for kind, gates in schedule:
                if kind == 'T':
                    self.assertEqual(len(gates), len({gate[1] for gate in gates}))
            self.assertEqual(len(word), 2 * (count + constant) + 2 * len(_schedule_word(middle)))

            def action(basis):
                parity_value = (sum(_bit(basis, q) for q in range(count)) + constant) % 2
                return basis ^ ((parity_value & _bit(basis, t)) << p)

            actual, expected = self.assert_native_action(width, word, action)
            reference = np.arange(1, 2 * (1 << width) + 1).reshape(1 << width, 2).astype(complex)
            reference[:, 1] *= 1j
            reference /= np.linalg.norm(reference)
            np.testing.assert_allclose(actual @ reference, expected @ reference,
                                       atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
