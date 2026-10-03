"""Bounded operator audits of shared-source conjugation and stale monitors.

The two-stage identity retains every source/preparation-work input at m=2.
PREP is native; masks and reflections are reduced exact operators, as in
other grouped fixtures. These checks are not a shallow-reflection emitter
or an unrestricted T-depth lower bound.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_grouped_program_prefetch import _ReducedGroup, _pack_program
    from .test_operator_source_compiler import _adjoint, _apply_native_word
except ImportError:
    from test_grouped_program_prefetch import _ReducedGroup, _pack_program
    from test_operator_source_compiler import _adjoint, _apply_native_word

ATOL = 8e-11


def _source_picture(group, columns, inverse=False):
    word = [('H', group.b)] + group.source.prep
    return _apply_native_word(group.width, _adjoint(word) if inverse else word,
                              columns)


def _middle(group, columns, program, ell, h=1, inverse=False, mask_only=False):
    """Literal diagonal mask followed by controlled K, without any PREP call."""
    enable = bool(h) & ((group.locals >> (ell + 1)) == 0)
    prefix = group.locals & ((1 << ell) - 1)
    masks = [group.masks(program, ell, int(p))[int(beta)]
             for p, beta in zip(prefix, group.branch)]
    parity = np.array([(int(core) & mask).bit_count() & 1
                       for core, mask in zip(group.core, masks)])
    phase = (-1.0) ** (enable & parity)
    if mask_only:
        return phase[:, None] * columns
    select = enable & group.branch
    target = (group.indices >> group.local[ell]) & 1
    kphase = (-1.0) ** (select & target)
    permutation = group.indices ^ (select.astype(int) << group.local[ell])
    if inverse:
        return phase[:, None] * columns[permutation] * kphase[:, None]
    return (phase[:, None] * kphase[:, None] * columns)[permutation]


def _reflection(group, columns, ell, h=1, transformed=True):
    enable = bool(h) & ((group.locals >> (ell + 1)) == 0)
    phase = np.ones(len(group.indices))
    phase[enable & (group.branch == 0) & (group.core == 0)] = -1
    result = _source_picture(group, columns, inverse=True) if transformed else columns
    result = phase[:, None] * result
    return _source_picture(group, result) if transformed else result


def _transformed_layer(group, columns, program, ell, h=1,
                       wrong_reflection=False):
    result = _middle(group, columns, program, ell, h)
    result = _reflection(group, result, ell, h, not wrong_reflection)
    result = _middle(group, result, program, ell, h, inverse=True)
    result = _reflection(group, result, ell, h, not wrong_reflection)
    result = _middle(group, result, program, ell, h)
    enable = bool(h) & ((group.locals >> (ell + 1)) == 0)
    return (-1.0) ** enable[:, None] * result


class GroupedSourceReuseTests(unittest.TestCase):
    def test_two_stage_common_conjugation_on_all_source_and_work_inputs(self):
        group = _ReducedGroup(2, 2)
        # Exact axis rows: the first stage changes the address of the second.
        program = _pack_program(2, [[(1, 0)], [(0, 1), (1, 0)]])
        inputs = np.eye(1 << group.width, dtype=complex)
        for h in (0, 1):
            expected = group.apply(inputs, program, h)
            result = _source_picture(group, inputs)
            for ell in range(2):
                result = _transformed_layer(group, result, program, ell, h)
            result = _source_picture(group, result, inverse=True)
            np.testing.assert_allclose(result, expected, atol=ATOL, rtol=0)
            if h == 0:
                np.testing.assert_allclose(result, inputs, atol=ATOL, rtol=0)

    def test_boundary_preparation_does_not_remove_conjugated_reflections(self):
        group = _ReducedGroup(2, 2)
        program = _pack_program(2, [[(1, 0)], [(0, 1), (1, 0)]])
        inputs = group.embedding()
        expected = group.apply(inputs, program)
        wrong = _source_picture(group, inputs)
        for ell in range(2):
            wrong = _transformed_layer(group, wrong, program, ell,
                                       wrong_reflection=True)
        wrong = _source_picture(group, wrong, inverse=True)
        self.assertGreater(np.linalg.norm(wrong - expected, 2), 1.0)

    def test_exact_zero_angle_exposes_stale_monitor_and_one_use_leakage(self):
        group = _ReducedGroup(3, 1)
        program = _pack_program(3, [[(0, 1)]])  # c=1, s=0 exactly.
        initialized = group.embedding()
        omega = _source_picture(group, initialized)
        projector = omega @ omega.conj().T
        identity = np.eye(len(group.indices), dtype=complex)
        mask = _middle(group, identity, program, 0, mask_only=True)
        np.testing.assert_allclose(omega.conj().T @ mask @ omega,
                                   np.eye(2) / 2, atol=ATOL, rtol=0)
        self.assertAlmostEqual(np.linalg.norm(mask @ projector - projector @ mask, 2),
                               np.sqrt(3) / 2, places=10)

        # Monitor a is the outer block: C_Pi maps zero to one precisely on Pi.
        monitor = np.block([[identity - projector, projector],
                            [projector, identity - projector]])
        monitored = monitor @ np.vstack([omega, np.zeros_like(omega)])
        acted = np.vstack([_middle(group, monitored[:len(identity)], program, 0),
                            _middle(group, monitored[len(identity):], program, 0)])
        uncomputed = monitor @ acted
        np.testing.assert_allclose(uncomputed[:len(identity)], omega / 2,
                                   atol=ATOL, rtol=0)
        rejected = uncomputed[len(identity):]
        np.testing.assert_allclose(rejected.conj().T @ rejected,
                                   3 * np.eye(2) / 4, atol=ATOL, rtol=0)

        # The one-use Q with native PREP has the same constant rejected norm;
        # its proper same-bank amplification returns the exact identity row.
        queried = group.apply_q(initialized, program, 0, h=1)
        zero = (group.branch == 0) & (group.core == 0)
        failure = queried[~zero]
        np.testing.assert_allclose(failure.conj().T @ failure,
                                   3 * np.eye(2) / 4, atol=ATOL, rtol=0)
        np.testing.assert_allclose(group.apply(initialized, program), initialized,
                                   atol=ATOL, rtol=0)


if __name__ == '__main__':
    unittest.main()
