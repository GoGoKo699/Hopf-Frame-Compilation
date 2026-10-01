"""Actual shared-source bodies for changing-target real tree groups.

The q=2 circuits audit exact native words on every signal and dirty input.
Their coarse target errors are diagnostics, not precision certification.
Eight-mode circuits use 128-dimensional matrices; the sixteen-mode
recurrence uses 256-dimensional matrices. No target or resolvent oracle
is used to form the candidate words.
"""
from __future__ import annotations

import unittest

import numpy as np

from tests.test_one_clean_compiler import (
    _amplification_word, _coefficients, _encoded_signs, _literal_mask,
    _mask_bits, _mask_table_word, _paired_source_word, _paired_weights,
    _pauli_t_word, _sandwich_word, _scalar_word,
)
from tests.test_operator_source_compiler import (
    ATOL, I2, X, Y, Z, _adjoint, _expand_toffolis, _pauli, _word_matrix,
)


def _borrowed_x(controls, target, helper):
    """Small exact dirty echo; only three/four controls occur here."""
    if len(controls) <= 2:
        return [(('X', 'CX', 'CCX')[len(controls)], *controls, target)]
    first = _borrowed_x(controls[:-1], helper, target)
    second = [('CCX', helper, controls[-1], target)]
    return (first + second) * 2


def _small_scalar(q, mask, flag, predicate, helper):
    if len(predicate) <= 2:
        return _scalar_word(q, mask, flag, predicate, helper)
    # The sole four-control source center occurs at the sixteen-mode root.
    # Recursive borrowed echoes swap target/helper roles and restore the
    # temporary helper exactly before any source loader or query resumes.
    source = _paired_source_word(q)
    inverse = _adjoint(source)
    middle = inverse + _borrowed_x(list(predicate), 0, helper) + source
    center = _borrowed_x([*predicate, flag], 0, helper)
    one = inverse + center + source
    zero = inverse + [('X', flag)] + center + [('X', flag)] + source
    return [('H', flag)] + one + _adjoint(mask) + middle + mask + zero + [('H', flag)]


def _scalar_inside_loader(q, word):
    """Strip the common outer U/U† after commuting them through signal H."""
    source = _paired_source_word(q)
    length = len(source)
    assert word[1:1 + length] == _adjoint(source)
    assert word[-1 - length:-1] == source
    return word[:1] + word[1 + length:-1 - length] + word[-1:]


def _fixed_inside_loader(q, flag):
    """The tail commutes with P_f, leaving two transformed two-T seeds."""
    seed = [('T', 0)] + _pauli_t_word({0: 'Y', 1: 'X'}, -1)
    mask = _literal_mask(_paired_weights(q)[1])
    transformed_inverse = seed + _adjoint(mask) + _adjoint(seed)
    transformed = seed + mask + _adjoint(seed)
    return ([('H', flag), ('CX', flag, 0)] + transformed_inverse + [('X', 0)]
            + transformed + [('X', flag), ('CX', flag, 0), ('X', flag), ('H', flag)])


def _small_mask(core, address, helper, rows):
    if not address:
        return _literal_mask(rows[0])
    if len(address) <= 2:
        return _mask_table_word(core, address, helper, rows)
    # A finite bitwise query for the three-address fixture. Its expanded
    # Toffolis are charged below; no scalable lookup bound is inferred.
    masks = [_mask_bits(signs) for signs in rows]

    def query(component):
        word = []
        for row, masks_at_row in enumerate(masks):
            flips = [('X', bit) for j, bit in enumerate(address)
                     if not (row >> j) & 1]
            word += flips
            for bit in range(core):
                if (masks_at_row[component] >> bit) & 1:
                    word += _borrowed_x(address, bit, helper)
            word += flips[::-1]
        return word

    hadamards = [('H', bit) for bit in range(core)]
    return query(0) + hadamards + query(1) + hadamards


def _group(q, n, depth, angles):
    core, flag = q + 1, q + 1 + n
    target = flag - 1 - depth
    address = list(range(target + 1, flag))
    predicate = tuple(range(core, target))
    rows = [_encoded_signs(q, theta) for theta in angles]
    mask = _small_mask(core, address, target, rows)
    flips = [('X', bit) for bit in predicate]
    scalar_core = _small_scalar(q, mask, flag, predicate, target)
    scalar = flips + scalar_core + flips
    scalar_inside = flips + _scalar_inside_loader(q, scalar_core) + flips
    controlled_y = [('SDG', target), ('CX', flag, target), ('S', target)]
    middle = controlled_y + [('SDG', flag)] + scalar + [('S', flag)] + controlled_y
    middle_inside = (controlled_y + [('SDG', flag)] + scalar_inside
                     + [('S', flag)] + controlled_y)
    return dict(target=target, predicate=predicate, rows=rows, mask=mask,
                scalar=scalar, middle=middle, middle_inside=middle_inside, flips=flips)


def _body(a, b):
    """Chronological emission of (B A^dagger B A)^2 B."""
    return b + a + b + _adjoint(a) + b + a + b + _adjoint(a) + b


def _body_matrix(a, b):
    return b @ a.conj().T @ b @ a @ b @ a.conj().T @ b @ a @ b


def _amplified_matrix(qword, reflection):
    return (qword @ reflection @ qword.conj().T @ reflection @ qword
            @ reflection @ qword.conj().T @ reflection @ qword)


def _ideal_and_accepted(q, n, depth, angles, rows):
    core, flag = q + 1, q + n + 1
    target_bit = n - 1 - depth
    target = core + target_bit
    ideal = np.eye(1 << flag, dtype=complex)
    accepted = ideal.copy()
    weights, fixed = _paired_weights(q)
    for logical in range(1 << n):
        if logical & ((1 << target_bit) - 1) or (logical >> target_bit) & 1:
            continue
        row = logical >> (target_bit + 1)
        angle = angles[row]
        s, p = _coefficients(weights, fixed, rows[row])
        for dirty in range(1 << core):
            first = (logical << core) | dirty
            indices = [first, first | (1 << target)]
            ideal[np.ix_(indices, indices)] = np.cos(angle) * I2 - 1j * np.sin(angle) * Y
            accepted[np.ix_(indices, indices)] = s * I2 - 1j * p * Y
    return ideal, accepted


def _parity_word(core, n, flag):
    return [gate for bit in range(core, core + n)
            for gate in [('H', bit), ('CX', flag, bit), ('H', bit)]]


def _t_count(word):
    return sum(gate[0] in ('T', 'TDG') for gate in word)


class JointSourceBodyTests(unittest.TestCase):
    def test_middle_route_matches_retained_primitive_on_all_ports(self):
        q, n = 2, 3
        flag, width = q + n + 1, q + n + 2
        weights, fixed = _paired_weights(q)
        a_word = _scalar_word(q, _literal_mask(fixed), flag)
        a = _word_matrix(width, a_word)
        identity = np.eye(1 << width)
        signal_z = _pauli(width, {flag: Z})
        np.testing.assert_allclose(a @ a @ a, -identity, atol=ATOL, rtol=0)
        np.testing.assert_allclose(signal_z @ a @ signal_z, a.conj().T,
                                   atol=ATOL, rtol=0)
        for precision in (2, 3, 5):
            source = _paired_source_word(precision)
            seed = [('T', 0)] + _pauli_t_word({0: 'Y', 1: 'X'}, -1)
            self.assertEqual(source[:len(seed)], seed)
            fixed_mask = _literal_mask(_paired_weights(precision)[1])
            transformed = _word_matrix(
                precision + 1, source + fixed_mask + _adjoint(source))
            seed_transformed = _word_matrix(
                precision + 1, seed + fixed_mask + _adjoint(seed))
            np.testing.assert_allclose(transformed, seed_transformed, atol=ATOL, rtol=0)
            self.assertEqual(_t_count(_fixed_inside_loader(precision, precision + 1)), 8)
        for depth, angles in enumerate(((.23,), (-.41, .37), (.18, 0, .47, -.33))):
            group = _group(q, n, depth, angles)
            middle = _word_matrix(width, _expand_toffolis(group['middle']))
            new = a.conj().T @ middle @ a
            old_word = (group['flips'] + _sandwich_word(
                q, group['mask'], group['target'], flag,
                predicate=group['predicate'], helper=group['target']) + group['flips'])
            old = _word_matrix(width, _expand_toffolis(old_word))
            target = group['target']
            route = _word_matrix(width, [('H', target), ('CX', flag, target), ('H', target)])
            np.testing.assert_allclose(new, route @ old @ route, atol=ATOL, rtol=0)
            np.testing.assert_allclose(signal_z @ middle @ signal_z,
                                       middle.conj().T, atol=ATOL, rtol=0)
            ideal, accepted = _ideal_and_accepted(q, n, depth, angles, group['rows'])
            np.testing.assert_allclose(new[:1 << flag, :1 << flag], accepted,
                                       atol=ATOL, rtol=0)
            amplified = _amplified_matrix(new, signal_z)
            np.testing.assert_allclose(amplified,
                                       route @ _amplified_matrix(old, signal_z) @ route,
                                       atol=ATOL, rtol=0)
            np.testing.assert_allclose(amplified, a.conj().T @ _body_matrix(a, middle) @ a,
                                       atol=ATOL, rtol=0)
            inactive = [basis for basis in range(1 << width)
                        if any((basis >> bit) & 1 for bit in group['predicate'])]
            np.testing.assert_allclose(amplified[:, inactive], identity[:, inactive],
                                       atol=ATOL, rtol=0)
            # The current target is a genuinely used dirty center helper at
            # the root and a query selector at the last level. Its arbitrary
            # coherence returns before the two outer controlled-Y routes.
            scalar = _word_matrix(width, _expand_toffolis(group['scalar']))
            for axis in (X, Z):
                target_axis = _pauli(width, {target: axis})
                np.testing.assert_allclose(scalar @ target_axis, target_axis @ scalar,
                                           atol=ATOL, rtol=0)
            # The independent encoded compression fixes the intended angle
            # sign, while the full SU(2) completion has the same error norm
            # as its initialized signal columns.
            desired_full = np.block([[ideal, np.zeros_like(ideal)],
                                     [np.zeros_like(ideal), ideal.conj().T]])
            full_error = amplified - desired_full
            self.assertAlmostEqual(np.linalg.norm(full_error, 2),
                                   np.linalg.norm(full_error[:, :1 << flag], 2), delta=ATOL)

    def test_literal_eight_mode_joint_word_and_complex_native_baseline(self):
        q, n = 2, 3
        core, flag, width = q + 1, q + n + 1, q + n + 2
        fixed = _paired_weights(q)[1]
        a_word = _scalar_word(q, _literal_mask(fixed), flag)
        levels = ((.23,), (-.41, .37), (.18, 0, .47, -.33))
        original, bodies, inside_original, inside_bodies = [], [], [], []
        inside_a = _fixed_inside_loader(q, flag)
        desired = np.eye(1 << flag, dtype=complex)
        error_budget = 0.0
        signal_z = _pauli(width, {flag: Z})
        a = _word_matrix(width, a_word)
        other_t = 0
        for depth, angles in enumerate(levels):
            group = _group(q, n, depth, angles)
            b_word = group['middle']
            original += _amplification_word(a_word + b_word + _adjoint(a_word), flag)
            bodies += _body(a_word, b_word)
            inside_b = group['middle_inside']
            inside_original += _amplification_word(
                inside_a + inside_b + _adjoint(inside_a), flag)
            inside_bodies += _body(inside_a, inside_b)
            b_native = _expand_toffolis(b_word)
            other_t += _t_count(b_native) - 3 * 4 * q
            b = _word_matrix(width, b_native)
            amplified = _amplified_matrix(a.conj().T @ b @ a, signal_z)
            ideal, _ = _ideal_and_accepted(q, n, depth, angles, group['rows'])
            target_full = np.block([[ideal, np.zeros_like(ideal)],
                                    [np.zeros_like(ideal), ideal.conj().T]])
            error_budget += np.linalg.norm(amplified - target_full, 2)
            desired = ideal @ desired
        compressed = a_word + bodies + _adjoint(a_word)
        original_native = _expand_toffolis(original)
        compressed_native = _expand_toffolis(compressed)
        self.assertEqual(_t_count(original_native), 45 * n * 4 * q + 5 * other_t)
        self.assertEqual(_t_count(compressed_native), (27 * n + 6) * 4 * q + 5 * other_t)
        source = _paired_source_word(q)
        simplified = _expand_toffolis(
            _adjoint(source) + inside_a + inside_bodies + _adjoint(inside_a) + source)
        simplified_baseline = _expand_toffolis(_adjoint(source) + inside_original + source)
        # Once common U boundaries and the fixed-mask tail are simplified,
        # BOTH words have the same leading precision charge. The improvement
        # below is in fixed scalar work, not in the coefficient of q.
        self.assertEqual(_t_count(simplified), (40 * n + 4) * q + 8 * (4 * n + 2) + 5 * other_t)
        self.assertEqual(_t_count(simplified_baseline), (40 * n + 4) * q + 80 * n + 5 * other_t)
        actual = _word_matrix(width, compressed_native)
        np.testing.assert_allclose(actual, _word_matrix(width, original_native),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(actual, _word_matrix(width, simplified),
                                   atol=ATOL, rtol=0)
        np.testing.assert_allclose(actual, _word_matrix(width, simplified_baseline),
                                   atol=ATOL, rtol=0)
        parity = _word_matrix(width, _parity_word(core, n, flag))
        corrected = parity @ actual @ parity
        target = np.kron(I2, desired)
        self.assertLessEqual(np.linalg.norm(corrected - target, 2), error_budget + ATOL)
        self.assertGreater(np.linalg.norm(actual[1 << flag:, :1 << flag]), .01)
        # A literal complex native coarse inverse changes the target to C†W
        # without asserting that this residual is real or SO(4)-factorable.
        coarse_word = [('T', core + 2), ('H', core + 1),
                       ('CX', core + 1, core), ('TDG', core)]
        coarse = _word_matrix(width, coarse_word)
        coarse_main = _word_matrix(flag, coarse_word)
        residual_actual = coarse.conj().T @ corrected
        residual_target = np.kron(I2, coarse_main.conj().T @ desired)
        self.assertGreater(np.linalg.norm(residual_target.imag), .1)
        self.assertAlmostEqual(np.linalg.norm(residual_actual - residual_target, 2),
                               np.linalg.norm(corrected - target, 2), delta=ATOL)

    def test_sixteen_mode_next_level_keeps_the_priced_recurrence(self):
        q, n = 2, 4
        flag, width = q + n + 1, q + n + 2
        a_word = _scalar_word(q, _literal_mask(_paired_weights(q)[1]), flag)
        a = _word_matrix(width, a_word)
        signal_z = _pauli(width, {flag: Z})
        original = np.eye(1 << width, dtype=complex)
        bodies = original.copy()
        levels = ((.19,), (-.28, .33), (.41, 0, -.22, .36),
                  (.12, -.17, .21, -.26, .31, -.38, .43, -.47))
        literal_t = compressed_t = 0
        inside_a = _fixed_inside_loader(q, flag)
        inside_bodies = []
        a_t = _t_count(a_word)
        other_t = 0
        for depth, angles in enumerate(levels):
            group = _group(q, n, depth, angles)
            b_native = _expand_toffolis(group['middle'])
            b = _word_matrix(width, b_native)
            primitive = a.conj().T @ b @ a
            _, accepted = _ideal_and_accepted(q, n, depth, angles, group['rows'])
            np.testing.assert_allclose(primitive[:1 << flag, :1 << flag], accepted,
                                       atol=ATOL, rtol=0)
            original = _amplified_matrix(primitive, signal_z) @ original
            bodies = _body_matrix(a, b) @ bodies
            literal_t += 5 * (_t_count(b_native) + 2 * a_t)
            compressed_t += 5 * _t_count(b_native) + 4 * a_t
            other_t += _t_count(b_native) - 3 * 4 * q
            inside_bodies += _body(inside_a, group['middle_inside'])
        compressed_t += 2 * a_t
        np.testing.assert_allclose(original, a.conj().T @ bodies @ a, atol=ATOL, rtol=0)
        self.assertEqual(literal_t, 180 * 4 * q + 5 * other_t)
        self.assertEqual(compressed_t, 114 * 4 * q + 5 * other_t)
        source = _paired_source_word(q)
        simplified = _expand_toffolis(
            _adjoint(source) + inside_a + inside_bodies + _adjoint(inside_a) + source)
        self.assertEqual(_t_count(simplified), (40 * n + 4) * q + 8 * (4 * n + 2) + 5 * other_t)
        # Counts above include the finite three-address query and four-control
        # source center. The source recurrence is still linear in tree height;
        # no uniform-n precision-width or endpoint improvement is asserted.


if __name__ == '__main__':
    unittest.main()
