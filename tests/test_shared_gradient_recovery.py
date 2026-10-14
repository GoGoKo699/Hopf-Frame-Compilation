"""Bounded audits of one phase-gradient source shared by a recovery word.

The two-bit arithmetic is a literal Clifford+T compute/swap/uncompute
fixture, with no claim to emit the asymptotic parallel adder. Larger
kickback checks use exact integer permutations and matrices of dimension
at most 64. A recovery-shaped word audits source reuse, actual inverses,
literal phases and live-gradient reflection masks; it is not a complete
frame compiler or an independent proof of the recovery theorem.
"""
from __future__ import annotations

import unittest

import numpy as np

try:
    from .test_counter_dirty_indicator import (
        _primitive, _read, _replace, _reverse, _serial,
    )
    from .test_operator_source_compiler import (
        H, _adjoint, _basis_action, _pauli, _rotation, _word_matrix,
    )
    from .test_phase_gradient_calibration import _gradient_preparation
    from .test_t_depth import _word
except ImportError:
    from test_counter_dirty_indicator import (
        _primitive, _read, _replace, _reverse, _serial,
    )
    from test_operator_source_compiler import (
        H, _adjoint, _basis_action, _pauli, _rotation, _word_matrix,
    )
    from test_phase_gradient_calibration import _gradient_preparation
    from test_t_depth import _word


ATOL = 1e-10
Q = np.diag([1, 1j]) @ H


def _conditional_clean_add(*, wrong_uncompute=False):
    """Two-bit B,F,C: XOR sum, swap F/C, XOR difference into C."""
    program, gradient, work = [0, 1], [2, 3], [4, 5]
    b0, b1 = program
    x0, x1 = gradient
    y0, y1 = work
    plus = [('CX', x0, y0), ('CX', b0, y0),
            ('CX', x1, y1), ('CX', b1, y1), ('CCX', x0, b0, y1)]
    # The high subtraction bit is x1 XOR b1 XOR b0 XOR (x0 AND b0).
    minus = plus + [('CX', b0, y1)]
    swap = [gate for left, right in zip(gradient, work) for gate in
            [('CX', left, right), ('CX', right, left), ('CX', left, right)]]
    circuit = _serial([_primitive(plus), _primitive(swap),
                       _primitive(plus if wrong_uncompute else minus)])
    return dict(program=program, gradient=gradient, work=work, circuit=circuit)


def _shared_echo(*, signed=True):
    fixture = _conditional_clean_add()
    target, activity = 6, 7
    complement_t = _primitive([('CX', target, wire) for wire in fixture['gradient']])
    complement_h = _primitive([('CX', activity, wire) for wire in fixture['gradient']])
    shift = fixture['circuit']
    if signed:
        shift = _serial([complement_t, shift, complement_t])
    fixture.update(target=target, activity=activity,
                   echoed=_serial([shift, complement_h, _reverse(shift), complement_h]))
    return fixture


def _table_shift(bits, table):
    """Unsigned doubled shifts, with a coherent logical address above F."""
    size = 1 << bits
    images = [address * size + (value + 2 * program) % size
              for address, program in enumerate(table) for value in range(size)]
    return np.eye(len(images), dtype=complex)[:, images]


def _y_call(bits, wire, program):
    """Two logical bits above F; no PREP or unPREP inside this call."""
    size = 1 << bits
    basis = _pauli(2, {wire: Q})
    table = [((-1) ** ((address >> wire) & 1)) * program for address in range(4)]
    lifted = np.kron(basis, np.eye(size))
    return (lifted @ _table_shift(bits, table) @ lifted.conj().T,
            _pauli(2, {wire: _rotation(4 * np.pi * program / size)}))


def _body(bits):
    """A small noncommuting native logical word with three gradient calls."""
    physical = np.eye(4 << bits, dtype=complex)
    logical = np.eye(4, dtype=complex)
    calls = [_y_call(bits, 0, 1)]
    for word in ([('H', 1)], [('CX', 0, 1)], [('T', 1)]):
        matrix = _word_matrix(2, word)
        calls.append((np.kron(matrix, np.eye(1 << bits)), matrix))
    calls.append(_y_call(bits, 1, 1))
    matrix = _word_matrix(2, [('CX', 1, 0)])
    calls.append((np.kron(matrix, np.eye(1 << bits)), matrix))
    calls.append(_y_call(bits, 0, 3))
    for native, ideal in calls:
        physical = native @ physical
        logical = ideal @ logical
    return physical, logical


def _recovery(body, reflection, rounds):
    """(-U S U^dagger S)^rounds U, retaining every literal minus sign."""
    result = body.copy()
    iterate = -body @ reflection @ body.conj().T @ reflection
    for _ in range(rounds):
        result = iterate @ result
    return result


def _coarse_native_preparation():
    # Deliberately perturb two phases: changing bit zero alone would commute
    # with every doubled shift and fail to exercise source leakage.
    word = [('H', bit) for bit in range(4)]
    word += [('T', 0), ('S', 2), ('Z', 3)]
    return _word_matrix(4, word), _word_matrix(4, _adjoint(word))


class SharedGradientRecoveryTests(unittest.TestCase):
    def test_native_compute_swap_difference_uncompute_and_actual_inverse(self):
        fixture = _conditional_clean_add()
        circuit = fixture['circuit']
        images = [_basis_action(value, circuit[0]) for value in range(64)]
        expected = np.eye(64, dtype=complex)[:, images]
        actual = _word_matrix(6, _word(circuit[1]))
        inverse = _word_matrix(6, _word(_reverse(circuit)[1]))
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(inverse, expected.conj().T, atol=ATOL, rtol=0)
        for program in range(4):
            for gradient in range(4):
                value = program | (gradient << 2)
                plus = program | (((gradient + program) % 4) << 2)
                minus = program | (((gradient - program) % 4) << 2)
                self.assertEqual(images[value], plus)
                np.testing.assert_allclose(inverse[:, value], np.eye(64)[:, minus],
                                           atol=ATOL, rtol=0)

        # A second sum is not the required subtraction uncompute.
        wrong = _conditional_clean_add(wrong_uncompute=True)['circuit']
        self.assertNotEqual(_read(_basis_action(1, wrong[0]), fixture['work']), 0)
        # The completed arithmetic promise requires C=0; the inactive echo
        # below supplies the separate guarantee for arbitrarily occupied C.
        occupied = 1 | (1 << 4)
        claimed = _replace(occupied, fixture['gradient'], 1)
        self.assertNotEqual(_basis_action(occupied, circuit[0]), claimed)

    def test_inactive_echo_returns_every_occupied_work_input(self):
        for signed in (False, True):
            fixture = _shared_echo(signed=signed)
            circuit = fixture['echoed']
            inverse = _reverse(circuit)[0]
            for value in range(256):
                actual = _basis_action(value, circuit[0])
                self.assertEqual(_basis_action(actual, inverse), value)
                active = value >> fixture['activity'] & 1
                if not active:
                    self.assertEqual(actual, value)
                elif _read(value, fixture['work']) == 0:
                    program = _read(value, fixture['program'])
                    gradient = _read(value, fixture['gradient'])
                    sign = (-1) ** (value >> fixture['target'] & 1) if signed else 1
                    expected = _replace(value, fixture['gradient'],
                                        (gradient + 2 * sign * program) % 4)
                    self.assertEqual(actual, expected)

        # Dropping the activity control changes an occupied-program H=0 input.
        fixture = _shared_echo(signed=False)
        unconditional = _primitive([('X', wire) for wire in fixture['gradient']])
        shift = fixture['circuit']
        wrong = _serial([shift, unconditional, _reverse(shift), unconditional])
        self.assertNotEqual(_basis_action(1, wrong[0]), 1)

    def test_coefficient_bank_preserves_literal_signed_sum_phase(self):
        bits, size = 4, 16
        coefficients = [1, 2]
        table = [sum((-1) ** ((address >> wire) & 1) * coefficient
                     for wire, coefficient in enumerate(coefficients)) % size
                 for address in range(4)]
        preparation = _gradient_preparation(bits)
        gamma = preparation[:, 0]
        source = np.kron(np.eye(4), gamma[:, None])
        phases = np.diag([np.exp(-4j * np.pi * program / size) for program in table])
        shifted = _table_shift(bits, table) @ source
        np.testing.assert_allclose(shifted, source @ phases, atol=ATOL, rtol=0)

        basis = _pauli(2, {0: Q, 1: Q})
        lift = np.kron(basis, np.eye(size))
        physical = lift @ _table_shift(bits, table) @ lift.conj().T
        wanted = _pauli(2, {wire: _rotation(4 * np.pi * coefficient / size)
                            for wire, coefficient in enumerate(coefficients)})
        np.testing.assert_allclose(physical @ source, source @ wanted, atol=ATOL, rtol=0)

        # Subtracting the baseline would discard a literal active-sector phase.
        # It becomes a relative phase against an inactive branch, not a gauge.
        wrong_table = [(program - table[0]) % size for program in table]
        wrong = lift @ _table_shift(bits, wrong_table) @ lift.conj().T @ source
        self.assertGreater(np.linalg.norm(wrong - source @ wanted, 2), 1)

    def test_fixed_damping_ignores_occupied_program_instead_of_xor_loading_it(self):
        bits, size, constant = 3, 8, 1
        gamma = _gradient_preparation(bits)[:, 0]
        source = np.kron(np.eye(2), gamma[:, None])
        basis = np.kron(Q, np.eye(size))
        desired = source @ _rotation(4 * np.pi * constant / size)

        def echoed_block(initial_program, *, hardwired):
            # Exact small query permutations, using the compute/swap/uncompute
            # boundaries audited natively above. No general shallow emitter
            # is inferred from these truth-table arithmetic operations.
            images = []
            for target in (0, 1):
                for gradient in range(size):
                    program, work = initial_program, 0
                    if not hardwired:
                        program ^= constant  # Load by XOR; not overwrite.
                    amount = constant if hardwired else program

                    # V_t = J_t A J_t, followed chronologically by
                    # J_H, the actual reversed V_t word, and J_H, at H=1.
                    if target:
                        gradient ^= size - 1
                    work ^= (gradient + amount) % size
                    gradient, work = work, gradient
                    work ^= (gradient - amount) % size
                    if target:
                        gradient ^= size - 1
                    gradient ^= size - 1
                    if target:
                        gradient ^= size - 1
                    work ^= (gradient - amount) % size
                    gradient, work = work, gradient
                    work ^= (gradient + amount) % size
                    if target:
                        gradient ^= size - 1
                    gradient ^= size - 1

                    if not hardwired:
                        program ^= constant  # Actual query inverse.
                    self.assertEqual(program, initial_program)
                    self.assertEqual(work, 0)
                    images.append(target * size + gradient)
            return np.eye(2 * size, dtype=complex)[:, images]

        for program in range(size):
            fixed = basis @ echoed_block(program, hardwired=True) @ basis.conj().T
            np.testing.assert_allclose(fixed @ source, desired, atol=ATOL, rtol=0)
            programmable = (basis @ echoed_block(program, hardwired=False)
                            @ basis.conj().T @ source)
            offset_angle = 4 * np.pi * (program ^ constant) / size
            np.testing.assert_allclose(programmable, source @ _rotation(offset_angle),
                                       atol=ATOL, rtol=0)
            if program == 0:
                np.testing.assert_allclose(programmable, desired, atol=ATOL, rtol=0)
            elif program == constant:
                # Occupied B=1 XOR the loaded constant 1 gives a zero angle.
                self.assertGreater(np.linalg.norm(programmable - desired, 2), 1)

    def test_one_actual_source_inverse_has_global_error_after_recovery(self):
        bits, size = 4, 16
        exact = _gradient_preparation(bits)
        actual, actual_inverse = _coarse_native_preparation()
        preparation_error = np.linalg.norm(actual - exact, 2)
        source_error = np.linalg.norm(actual[:, 0] - exact[:, 0])
        self.assertGreater(source_error, .1)
        self.assertLess(preparation_error, 1)
        embedding = np.eye(4 * size, dtype=complex)[:, np.arange(4) * size]
        logical_reflection = np.diag([-1, -1, 1, 1])
        reflection = np.kron(logical_reflection, np.eye(size))
        physical, logical = _body(bits)
        largest_leakage = 0.0
        for rounds in (1, 4, 11):
            with self.subTest(body_calls=2 * rounds + 1):
                recovery = _recovery(physical, reflection, rounds)
                desired = embedding @ _recovery(logical, logical_reflection, rounds)
                ideal = (np.kron(np.eye(4), exact.conj().T) @ recovery
                         @ np.kron(np.eye(4), exact) @ embedding)
                np.testing.assert_allclose(ideal, desired, atol=ATOL, rtol=0)
                # PREP and its actual native inverse occur only at the two ends.
                got = (np.kron(np.eye(4), actual_inverse) @ recovery
                       @ np.kron(np.eye(4), actual) @ embedding)
                self.assertLessEqual(np.linalg.norm(got - desired, 2),
                                     2 * source_error + ATOL)
                self.assertLessEqual(np.linalg.norm(got - desired, 2),
                                     2 * preparation_error + ATOL)
                off_zero = np.delete(got, np.arange(4) * size, axis=0)
                largest_leakage = max(largest_leakage, np.linalg.norm(off_zero, 2))
        self.assertGreater(largest_leakage, .1)

    def test_live_gradient_is_excluded_from_reflection_and_actual_inverse_is_used(self):
        bits, size = 4, 16
        exact = _gradient_preparation(bits)
        embedding = np.eye(4 * size, dtype=complex)[:, np.arange(4) * size]
        logical_reflection = np.diag([-1, -1, 1, 1])
        physical, logical = _body(bits)
        # Wrong mask: it tests F=0 while F holds the live Fourier source.
        wrong_reflection = np.eye(4 * size, dtype=complex)
        wrong_reflection[0, 0] = wrong_reflection[size, size] = -1
        wrong = (np.kron(np.eye(4), exact.conj().T)
                 @ _recovery(physical, wrong_reflection, 1)
                 @ np.kron(np.eye(4), exact) @ embedding)
        desired = embedding @ _recovery(logical, logical_reflection, 1)
        self.assertGreater(np.linalg.norm(wrong - desired, 2), .2)
        self.assertGreater(np.linalg.norm(np.delete(wrong, np.arange(4) * size, axis=0), 2), .1)

        actual, actual_inverse = _coarse_native_preparation()
        # On the inactive sector the echoed body is literally I on occupied
        # gradient inputs, so the two endpoint words must cancel as operators.
        inactive = np.kron(np.eye(4), actual_inverse @ actual)
        np.testing.assert_allclose(inactive, np.eye(4 * size), atol=ATOL, rtol=0)
        wrong_inverse = np.kron(np.eye(4), exact.conj().T @ actual)
        self.assertGreater(np.linalg.norm(wrong_inverse - np.eye(4 * size), 2), .1)


if __name__ == '__main__':
    unittest.main()
