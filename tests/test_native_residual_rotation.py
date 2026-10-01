"""Literal unaddressed words and certified coefficient-to-word integration.

The native checks include every input on eight wires at q=5. Their strong
oracle is the independent Majorana algebra, not the loose small-q error
constant. Fine-precision checks use both small Clifford-algebra irreps;
they do not simulate a large precision core or an addressed state compiler.
"""
from dataclasses import replace
from fractions import Fraction
import unittest

import numpy as np

from compiler_robust_hopf.native_residual_rotation import emit_residual_row, emit_rotation
from compiler_robust_hopf.residual_table_preprocessing import CoefficientInterval
from compiler_robust_hopf.rotation_programming import program_rotation


I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
H = (X + Z) / np.sqrt(2)
ATOL = 3e-11


def _point(value):
    value = Fraction(value)
    return CoefficientInterval(value, value)


def _midpoint(interval):
    return float((interval.lower + interval.upper) / 2)


def _native_matrix(gates, width, columns=None):
    """Independent elementary-gate propagation, all input columns by default."""
    size = 1 << width
    result = (np.eye(size, dtype=complex) if columns is None
              else np.array(columns, dtype=complex, copy=True))
    indices = np.arange(size)
    low = [indices[(indices & (1 << q)) == 0] for q in range(width)]
    high = [low[q] | (1 << q) for q in range(width)]
    phases = {"Z": -1, "S": 1j, "SDG": -1j,
              "T": (1 + 1j) / np.sqrt(2), "TDG": (1 - 1j) / np.sqrt(2)}
    for name, wires in gates:
        if name == "X":
            result = result[indices ^ (1 << wires[0])]
        elif name == "CX":
            control, target = wires
            result = result[indices ^ (((indices >> control) & 1) << target)]
        elif name == "H":
            q = wires[0]
            first, second = result[low[q]].copy(), result[high[q]].copy()
            result[low[q]] = (first + second) / np.sqrt(2)
            result[high[q]] = (first - second) / np.sqrt(2)
        else:
            result[high[wires[0]]] *= phases[name]
    return result


def _tensor(*factors):
    result = np.array([[1]], dtype=complex)
    for factor in factors:
        result = np.kron(result, factor)
    return result


def _majorana_triplet(program):
    """Source vectors from their mathematical coefficients, not emitted gates."""
    q, core = program.q, program.q + 1
    weights = [Fraction(0)] * (2 * core)
    weights[0] = Fraction(1, 4)
    for j in range(q):
        geometric = Fraction(1, 1 << min(j + 1, q - 1))
        weights[2 * j + 1], weights[2 * j + 2] = geometric / 2, geometric / 4
    matrices = [np.zeros((1 << core, 1 << core), dtype=complex) for _ in range(3)]
    for index, weight in enumerate(weights):
        if not weight:
            continue
        qubit, axis = divmod(index, 2)
        gamma = _tensor(*(I if bit > qubit else (X, Y)[axis] if bit == qubit else Z
                          for bit in reversed(range(core))))
        component = np.sqrt(float(weight)) * gamma
        fixed_sign = -1 if 2 <= index <= 2 * q and index % 2 == 0 else 1
        matrices[0] += component
        matrices[1] += fixed_sign * component
        matrices[2] += (1 - 2 * program.signs[index]) * component
    return matrices


def _algebra_rotation(program, triplet, axis="y"):
    source, fixed, programmed = triplet
    core_identity = np.eye(len(source))
    df = source @ fixed - core_identity / 2
    dg = source @ programmed - float(program.s) * core_identity
    whole_identity = np.eye(4 * len(source))
    left = whole_identity / 2 + _tensor(Z, X, df)
    middle = float(program.s) * whole_identity + _tensor(X, X, dg)
    primitive = left.conj().T @ middle @ left
    signs = np.tile(np.repeat([1, -1], len(source)), 2)
    amplified = primitive.copy()
    for call in (primitive.conj().T, primitive, primitive.conj().T, primitive):
        amplified = call @ (signs[:, None] * amplified)
    if axis == "z":
        conjugator = _tensor(H @ np.diag([1, -1j]), I, core_identity)
        amplified = conjugator @ amplified @ conjugator.conj().T
    return amplified


def _small_triplet(program, chirality):
    s, p = float(program.s), float(program.p)
    transverse = -p / np.sqrt(3 / 4)
    squared = 1 - s * s - transverse * transverse
    if squared < -1e-14:
        raise AssertionError("Programmed source Gram matrix is not feasible.")
    return X, X / 2 + np.sqrt(3 / 4) * Z, (
        s * X + transverse * Z + chirality * np.sqrt(max(0, squared)) * Y)


class NativeResidualRotationTests(unittest.TestCase):
    def test_literal_native_rotation_on_every_target_signal_and_core_input(self):
        program = program_rotation(_point(Fraction(3, 5)), _point(Fraction(4, 5)), 5)
        word = emit_rotation(program)
        actual = _native_matrix(word.gates, word.target + 1)
        expected = _algebra_rotation(program, _majorana_triplet(program))
        np.testing.assert_allclose(actual, expected, atol=ATOL, rtol=0)
        core_size = 1 << len(word.core)
        signal_zero = np.concatenate((np.arange(core_size), np.arange(core_size) + 2 * core_size))
        s, p = float(program.s), float(program.p)
        polynomial = 5 - 20 * (s * s + p * p) + 16 * (s * s + p * p) ** 2
        np.testing.assert_allclose(actual[np.ix_(signal_zero, signal_zero)],
                                   np.kron(polynomial * (s * I + p * X @ Z), np.eye(core_size)),
                                   atol=ATOL, rtol=0)
        flip = np.arange(len(actual)) ^ (1 << word.signal)
        np.testing.assert_allclose(actual[flip], actual[:, flip], atol=ATOL, rtol=0)
        # Literal inverse and every dirty input; no global-phase alignment.
        np.testing.assert_allclose(actual.conj().T @ actual, np.eye(len(actual)),
                                   atol=ATOL, rtol=0)
        ideal = np.kron(Fraction(3, 5) * I + Fraction(4, 5) * X @ Z,
                        np.eye(2 * core_size)).astype(complex)
        self.assertLess(np.linalg.norm(actual - ideal, ord=2), .2)

    def test_complex_residual_row_native_axes_and_identical_outer_word(self):
        row = emit_residual_row(Fraction(-3, 4), Fraction(1, 2), 5)
        outer, middle, repeated = row.rotations
        self.assertIs(outer, repeated)
        self.assertIs(row.programs[0], row.programs[2])
        self.assertEqual(row.gates, outer.gates + middle.gates + outer.gates)
        # Propagate each distinct word once, then compose. Every input is kept;
        # the emitted concatenation above has exactly this chronological action.
        z_actual = _native_matrix(outer.gates, outer.target + 1)
        y_actual = _native_matrix(middle.gates, middle.target + 1)
        z_expected = _algebra_rotation(outer.certificate, _majorana_triplet(outer.certificate), "z")
        y_expected = _algebra_rotation(middle.certificate, _majorana_triplet(middle.certificate))
        np.testing.assert_allclose(z_actual, z_expected, atol=ATOL, rtol=0)
        np.testing.assert_allclose(y_actual, y_expected, atol=ATOL, rtol=0)
        actual, expected = z_actual @ y_actual @ z_actual, z_expected @ y_expected @ z_expected
        np.testing.assert_allclose(actual, expected, atol=3 * ATOL, rtol=0)
        inverse_names = {"T": "TDG", "TDG": "T", "S": "SDG", "SDG": "S"}
        # The inverse of the entire word is obtained literally, not by
        # recompiling a negated angle or selecting a new half-phase branch.
        inverse = tuple((inverse_names.get(name, name), wires)
                        for name, wires in reversed(outer.gates))
        selected = [0, 63, 128, 255]
        np.testing.assert_allclose(
            _native_matrix(inverse, outer.target + 1, z_actual[:, selected]),
            np.eye(len(z_actual))[:, selected], atol=2 * ATOL, rtol=0)
        np.testing.assert_allclose(actual.conj().T @ actual, np.eye(len(actual)),
                                   atol=3 * ATOL, rtol=0)

    def test_fine_precision_full_signal_bound_and_residual_boundary_cases(self):
        for q in (12, 20):
            for cosine, sine in ((1, 0), (0, 1), (-1, 0), (Fraction(3, 5), Fraction(-4, 5))):
                program = program_rotation(_point(cosine), _point(sine), q)
                for axis, pauli in (("y", Y), ("z", Z)):
                    ideal = _tensor(float(cosine) * I - 1j * float(sine) * pauli, I, I)
                    for chirality in (-1, 1):
                        actual = _algebra_rotation(program, _small_triplet(program, chirality), axis)
                        self.assertLess(np.linalg.norm(actual - ideal, ord=2),
                                        float(program.operator_error_bound))
        q = 16
        tiny = Fraction(1, 1 << 70)
        for real, imag in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1),
                           (Fraction(-3, 4), tiny), (Fraction(-3, 4), -tiny)):
            row = emit_residual_row(real, imag, q)
            self.assertIs(row.rotations[0], row.rotations[2])
            result = I.copy()
            for axis, pair in zip((Z, Y, Z), row.rotation_data.rotations):
                result = (_midpoint(pair[0]) * I - 1j * _midpoint(pair[1]) * axis) @ result
            z = complex(real, imag)
            complement = np.sqrt(max(0, 1 - abs(z) ** 2))
            target = np.array([[z, -complement], [complement, z.conjugate()]])
            self.assertLessEqual(np.linalg.norm(result - target, ord=2),
                                 float(row.rotation_data.completion_error_bound))
            self.assertLess(row.operator_error_bound, Fraction(130, 1 << q))

    def test_linear_native_resource_ledger_and_rejected_forged_certificates(self):
        alphabet = {"X", "Z", "H", "S", "SDG", "T", "TDG", "CX"}
        for q in (5, 16, 64, 129):
            program = program_rotation(_point(Fraction(3, 5)), _point(Fraction(4, 5)), q)
            for axis in ("y", "z"):
                word = emit_rotation(program, axis)
                self.assertEqual(word.core, tuple(range(q + 1)))
                self.assertEqual((word.signal, word.target), (q + 1, q + 2))
                self.assertIs(word.certificate, program)
                self.assertEqual(sum(name in ("T", "TDG") for name, _ in word.gates), 180 * q)
                self.assertLessEqual(len(word.gates), 1680 * q)
                for name, wires in word.gates:
                    self.assertIn(name, alphabet)
                    self.assertEqual(len(wires), 2 if name == "CX" else 1)
                    self.assertEqual(len(wires), len(set(wires)))
                    self.assertTrue(all(0 <= wire <= word.target for wire in wires))
            row = emit_residual_row(Fraction(1, 2), Fraction(1, 4), q)
            self.assertEqual(sum(name in ("T", "TDG") for name, _ in row.gates), 540 * q)
            self.assertLessEqual(len(row.gates), 5040 * q)
        for forged in (replace(program, signs=(1 - program.signs[0],) + program.signs[1:]),
                       replace(program, operator_error_bound=Fraction(0)),
                       replace(program, s=program.s + 1)):
            with self.assertRaises(ValueError):
                emit_rotation(forged)
        with self.assertRaises(TypeError):
            emit_rotation(None)
        with self.assertRaises(ValueError):
            emit_rotation(program, "x")


if __name__ == "__main__":
    unittest.main()
