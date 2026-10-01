"""A bounded two-row residual table with one unchanged enable literal.

Every gate is elementary Clifford+T. The address and enable are logical
inputs, not initialized work; the precision core and signal are arbitrary
dirty inputs. On the inactive enable sector the whole word is exactly I.
On the active sector the certified full-operator error takes the maximum
over the two address rows, including all work and external references.

This single-address-bit construction uses Clifford mask selection and
needs no lookup bank, selector, or predicate helper. It does not implement
general addressed lookup or the complete state-preparation circuit.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .native_residual_rotation import Word, _gate, _mask, _rotation_gates, _validate_program
from .residual_table_preprocessing import ResidualRotationData, residual_rotation_data
from .rotation_programming import RotationProgram, program_rotation


@dataclass(frozen=True)
class NativeRotationTable:
    """Two active Ry/Rz rows and exact inactive identity on every input."""

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    address: int
    enable: int
    enable_value: int
    gates: Word
    programs: tuple[RotationProgram, RotationProgram]
    operator_error_bound: Fraction


@dataclass(frozen=True)
class ResidualNativeTable:
    """Three addressed rotations, with the identical outer table used twice.

    ``programs`` is ordered by Euler factor, then by address row; the
    first and last program tuples and rotation objects are identical.
    ``rotation_data`` contains the two coefficient-preprocessing results.
    The error includes preprocessing and all three native factors.
    """

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    address: int
    enable: int
    enable_value: int
    gates: Word
    rotation_data: tuple[ResidualRotationData, ResidualRotationData]
    programs: tuple[tuple[RotationProgram, RotationProgram], ...]
    rotations: tuple[NativeRotationTable, NativeRotationTable, NativeRotationTable]
    operator_error_bound: Fraction


def _enable_literal(value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value not in (0, 1):
        raise ValueError("enable_value must be the integer zero or one, not a boolean.")
    return int(value)


def _pair(values, name):
    try:
        pair = tuple(values)
    except TypeError as error:
        raise TypeError(f"{name} must contain exactly two entries.") from error
    if len(pair) != 2:
        raise ValueError(f"{name} must contain exactly two entries.")
    return pair


def _table_mask(programs: tuple[RotationProgram, RotationProgram], address: int) -> Word:
    masks = tuple(_mask(program.signs) for program in programs)
    gates = []
    # Keep all X gates before all Z gates. The actual inverse in the scalar
    # word retains any row-dependent phase of this ordered Pauli product.
    for axis in ("X", "Z"):
        first, second = ({wires[0] for name, wires in mask if name == axis} for mask in masks)
        gates.extend(_gate(axis, target) for target in sorted(first))
        for target in sorted(first ^ second):
            if axis == "X":
                gates.append(_gate("CX", address, target))
            else:
                gates.extend((_gate("H", target), _gate("CX", address, target),
                              _gate("H", target)))
    return tuple(gates)


def emit_rotation_table(
    programs: tuple[RotationProgram, RotationProgram], axis: str = "y", *, enable_value: int = 1,
) -> NativeRotationTable:
    """Emit a two-row rotation with one positive or negative enable literal.

    The layout is core 0..q, signal q+1, target q+2, address q+3,
    enable q+4. Each program must be a genuine certified output at the
    same q. Its common-angle enclosure promise remains the caller's.
    No additional work is used. Source centers alone are conditioned;
    the native Toffoli has the enable as its first control. Negative
    enable uses X wrappers around the complete amplified rotation.
    """
    programs = _pair(programs, "programs")
    enable_value = _enable_literal(enable_value)
    for program in programs:
        _validate_program(program)
    if programs[0].q != programs[1].q:
        raise ValueError("Both rotation programs must have the same q.")
    q = programs[0].q
    core, signal, target = tuple(range(q + 1)), q + 1, q + 2
    address, enable = q + 3, q + 4
    gates = _rotation_gates(q, _table_mask(programs, address), axis, enable=enable)
    if enable_value == 0:
        wrapper = (_gate("X", enable),)
        gates = wrapper + gates + wrapper
    return NativeRotationTable(q, core, signal, target, address, enable, enable_value,
                               gates, programs, max(p.operator_error_bound for p in programs))


def emit_residual_table(
    coefficients: tuple[tuple[Fraction | int, Fraction | int], ...], q: int, *, enable_value: int = 1,
) -> ResidualNativeTable:
    """Emit two residual U(z) rows, active only at the chosen enable value.

    Each (real,imag) is exact dyadic data with the same certified input
    promise as ``residual_rotation_data``. The active full-operator error
    is below 130*2**(-q); the inactive action is exactly I. Address,
    enable, target, signal and core inputs may all be coherent or entangled.
    Active U(0) is a pi/2 Ry rotation, not an inactive identity row.
    """
    enable_value = _enable_literal(enable_value)
    if isinstance(q, bool) or not isinstance(q, Integral) or q < 5:
        raise ValueError("q must be an integer at least five, not a boolean.")
    q = int(q)
    coefficients = tuple(_pair(row, "coefficient row") for row in _pair(coefficients, "coefficients"))
    if any(isinstance(value, bool) or not isinstance(value, (Fraction, Integral))
           for row in coefficients for value in row):
        raise TypeError("Coefficient entries must be exact dyadic Fractions or integers, not booleans.")
    data = tuple(residual_rotation_data(real, imag, q) for real, imag in coefficients)
    outer_programs = tuple(program_rotation(*row.rotations[0], q) for row in data)
    middle_programs = tuple(program_rotation(*row.rotations[1], q) for row in data)
    outer = emit_rotation_table(outer_programs, "z", enable_value=enable_value)
    middle = emit_rotation_table(middle_programs, "y", enable_value=enable_value)
    programs = outer_programs, middle_programs, outer_programs
    rotations = outer, middle, outer
    error = max(row.completion_error_bound for row in data) + sum(
        (rotation.operator_error_bound for rotation in rotations), Fraction(0))
    return ResidualNativeTable(
        q, outer.core, outer.signal, outer.target, outer.address, outer.enable, enable_value,
        outer.gates + middle.gates + outer.gates, data, programs, rotations, error)
