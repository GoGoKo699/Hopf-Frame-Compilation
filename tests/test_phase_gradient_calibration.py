"""Bounded checks of the guarded binary phase-gradient calibration.

The arithmetic reuses the emitted helper-free TTK adder. Exhaustive
integer permutations include occupied programs and gradient registers;
one two-bit native matrix checks literal Clifford+T phases. Larger
kickback checks use small exact permutation matrices, not a complete
frame or lookup emitter. A deliberately coarse native product source
checks approximate source return with its actual inverse and leakage.
Asymptotic synthesis, lookup, count and depth bounds remain analytic.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_counter_dirty_indicator import (
        _add, _primitive, _read, _replace, _resources, _reverse, _serial,
    )
    from .test_operator_source_compiler import (
        H, _adjoint, _basis_action, _rotation, _word_matrix,
    )
    from .test_t_depth import _word
except ImportError:
    from test_counter_dirty_indicator import (
        _add, _primitive, _read, _replace, _resources, _reverse, _serial,
    )
    from test_operator_source_compiler import (
        H, _adjoint, _basis_action, _rotation, _word_matrix,
    )
    from test_t_depth import _word


ATOL = 8e-11


def _guarded_shift(bits):
    """Little-endian words B,F; target t and activity z above them."""
    program = list(range(bits))
    gradient = list(range(bits, 2 * bits))
    target, activity = 2 * bits, 2 * bits + 1
    complement_t = _primitive([('CX', target, wire) for wire in gradient])
    complement_z = _primitive([('CX', activity, wire) for wire in gradient])
    signed = _serial([complement_t, _add(program, gradient), complement_t])
    echoed = _serial([signed, complement_z, _reverse(signed), complement_z])
    return dict(program=program, gradient=gradient, target=target,
                activity=activity, width=2 * bits + 2,
                signed=signed, echoed=echoed)


def _gradient_preparation(bits):
    """Exact mathematical product preparation, with its common phase fixed."""
    size = 1 << bits
    walsh = np.array([[1.0]])
    for _ in range(bits):
        walsh = np.kron(H, walsh)
    return np.exp(2j * np.pi * np.arange(size) / size)[:, None] * walsh


def _shift_matrix(bits, program, activity=1, *, signed=True):
    """A preserved-program block, on target tensor gradient."""
    size = 1 << bits
    images = []
    for basis in range(2 * size):
        target, gradient = divmod(basis, size)
        direction = (-1) ** target if signed else 1
        images.append(target * size + (gradient + 2 * activity * direction * program) % size)
    return np.eye(2 * size, dtype=complex)[:, images]


def _calibration(bits, program, preparation, activity=1, *, inverse=None, signed=True):
    """Chronological PREP, (SH)^dagger, echo, SH, actual PREP^dagger."""
    size = 1 << bits
    basis = np.diag([1, 1j]) @ H
    pre = np.kron(np.eye(2), preparation)
    unpre = np.kron(np.eye(2), preparation.conj().T if inverse is None else inverse)
    return (unpre @ np.kron(basis, np.eye(size))
            @ _shift_matrix(bits, program, activity, signed=signed)
            @ np.kron(basis.conj().T, np.eye(size)) @ pre)


def _zero_gradient_embedding(bits):
    size = 1 << bits
    return np.eye(2 * size, dtype=complex)[:, [0, size]]


class PhaseGradientCalibrationTests(unittest.TestCase):
    def test_emitted_doubled_signed_shift_every_program_and_core_offset(self):
        for bits in range(1, 6):
            fixture = _guarded_shift(bits)
            circuit = fixture['echoed']
            size = 1 << bits
            for value in range(1 << fixture['width']):
                program = _read(value, fixture['program'])
                gradient = _read(value, fixture['gradient'])
                target = value >> fixture['target'] & 1
                activity = value >> fixture['activity'] & 1
                expected = _replace(value, fixture['gradient'],
                                    (gradient + 2 * activity * (-1) ** target * program) % size)
                actual = _basis_action(value, circuit[0])
                self.assertEqual(actual, expected)
                self.assertEqual(_basis_action(actual, _reverse(circuit)[0]), value)
            # Two literal TTK words and Clifford complements; no helper wire.
            self.assertEqual(_resources(circuit), (28 * (bits - 1), 16 * (bits - 1)))
            self.assertTrue(all(max(gate[1:]) < fixture['width'] for gate in circuit[0]))
            for gate in _word(circuit[1]):
                if gate[0] in ('T', 'TDG'):
                    self.assertNotIn(gate[1], (fixture['target'], fixture['activity']))

    def test_two_bit_emitted_native_word_has_literal_phase_and_actual_inverse(self):
        fixture = _guarded_shift(2)
        circuit = fixture['echoed']
        images = [_basis_action(value, circuit[0]) for value in range(1 << fixture['width'])]
        expected = np.eye(len(images), dtype=complex)[:, images]
        actual = _word_matrix(fixture['width'], _word(circuit[1]))
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(
            _word_matrix(fixture['width'], _word(_reverse(circuit)[1])),
            expected.conj().T, atol=ATOL, rtol=0)

        # Include literal source preparation, Y-basis changes, and PREP inverse.
        low, high = fixture['gradient']
        target = fixture['target']
        prep = [('H', low), ('H', high), ('S', low), ('Z', high)]
        word = (prep + [('SDG', target), ('H', target)] + _word(circuit[1])
                + [('H', target), ('S', target)] + _adjoint(prep))
        macro = _word_matrix(fixture['width'], word)
        columns, desired = [], []
        for activity in (0, 1):
            for program in range(4):
                for bit in (0, 1):
                    value = program | (bit << target) | (activity << fixture['activity'])
                    columns.append(value)
                    # theta=pi*program: odd rows are literally -I, not I.
                    desired.append(((-1) ** (activity * program)) * np.eye(64)[:, value])
        np.testing.assert_allclose(macro[:, columns], np.array(desired).T,
                                   atol=ATOL, rtol=0)

    def test_phase_gradient_kickback_is_literal_y_rotation_for_all_rows(self):
        for bits in (2, 3, 4):
            size = 1 << bits
            preparation = _gradient_preparation(bits)
            gamma = preparation[:, 0]
            np.testing.assert_allclose(gamma,
                                       np.exp(2j * np.pi * np.arange(size) / size) / np.sqrt(size),
                                       atol=ATOL, rtol=0)
            embedding = _zero_gradient_embedding(bits)
            for program in range(size):
                angle = 4 * np.pi * program / size
                actual = _calibration(bits, program, preparation) @ embedding
                np.testing.assert_allclose(actual, embedding @ _rotation(angle),
                                           atol=ATOL, rtol=0)
        # Dropping the target sign produces a scalar phase, not the real rotation.
        bits, program = 4, 1
        preparation = _gradient_preparation(bits)
        embedding = _zero_gradient_embedding(bits)
        wrong = _calibration(bits, program, preparation, signed=False) @ embedding
        self.assertGreater(np.linalg.norm(wrong - embedding @ _rotation(np.pi / 4), 2), 1)

    def test_inactive_identity_on_occupied_program_core_and_target(self):
        bits = 3
        fixture = _guarded_shift(bits)
        for value in range(1 << fixture['activity']):
            self.assertEqual(_basis_action(value, fixture['echoed'][0]), value)

        # A missing activity control fails on an occupied program at z=0.
        unconditional = _primitive([('X', wire) for wire in fixture['gradient']])
        wrong = _serial([fixture['signed'], unconditional,
                         _reverse(fixture['signed']), unconditional])
        self.assertNotEqual(_basis_action(1, wrong[0]), 1)

        # Full blocks include every gradient and target input, not just |0>.
        preparation = _gradient_preparation(bits)
        for program in range(1 << bits):
            np.testing.assert_allclose(
                _calibration(bits, program, preparation, activity=0),
                np.eye(2 << bits), atol=ATOL, rtol=0)

    def test_actual_native_product_source_inverse_bounds_error_and_leakage(self):
        bits, size = 4, 16
        exact = _gradient_preparation(bits)
        # The first two phases should be pi/8 and pi/4. T and I give
        # deliberately coarse native approximations; S and Z are exact.
        # Perturbing only bit zero would miss leakage: doubled shifts preserve
        # that bit, so its phase error alone commutes with the arithmetic.
        prep_word = [('H', bit) for bit in range(bits)]
        prep_word += [('T', 0), ('S', 2), ('Z', 3)]
        actual = _word_matrix(bits, prep_word)
        actual_inverse = _word_matrix(bits, _adjoint(prep_word))
        np.testing.assert_allclose(actual_inverse, actual.conj().T, atol=ATOL, rtol=0)
        preparation_error = np.linalg.norm(actual - exact, 2)
        self.assertGreater(preparation_error, .1)
        self.assertLess(preparation_error, 1)

        embedding = _zero_gradient_embedding(bits)
        largest_leakage = 0.0
        for program in (1, 3, 5, 11):
            got = _calibration(bits, program, actual, inverse=actual_inverse) @ embedding
            wanted = embedding @ _rotation(4 * np.pi * program / size)
            self.assertLessEqual(np.linalg.norm(got - wanted, 2), 2 * preparation_error + ATOL)
            off_zero = np.delete(got, [0, size], axis=0)
            largest_leakage = max(largest_leakage, np.linalg.norm(off_zero, 2))
        self.assertGreater(largest_leakage, .1)

        # The actual inverse gives inactive full-space return, even for an
        # imperfect source. Substituting the ideal inverse leaves a real error.
        got = _calibration(bits, 5, actual, activity=0, inverse=actual_inverse)
        np.testing.assert_allclose(got, np.eye(2 * size), atol=ATOL, rtol=0)
        wrong = _calibration(bits, 5, actual, activity=0, inverse=exact.conj().T)
        self.assertAlmostEqual(np.linalg.norm(wrong - np.eye(2 * size), 2),
                               preparation_error, places=10)


if __name__ == '__main__':
    unittest.main()
