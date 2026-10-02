"""A bounded four-row native residual lookup with no lookup workspace.

Two logical address bits select one literal Pauli mask. The mask bits
have degree at most two in the addresses. For each Pauli axis, a Clifford
fanout or parity computation reduces the entire quadratic support to one
exact seven-T Toffoli or CCZ. The core inputs remain arbitrary; no target,
signal, or additional helper is used by mask selection.

Every complete mask returns both address bits. Native Toffoli expansions
temporarily change the high address, so this return is a word-boundary
identity, not a gate-by-gate control convention. All X factors precede
all Z factors, retaining the literal phases of the original masks.

Only source centers receive the separate unchanged enable condition.
The inactive sector is exactly identity; active error takes the maximum
over the four rows. This is two-address-bit lookup, not general lookup
or a complete multi-qubit state compiler.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .native_residual_rotation import (
    Word, _gate, _inverse, _mask, _rotation_gates, _toffoli, _validate_program,
)
from .native_residual_table import _enable_literal, _pair
from .residual_table_preprocessing import ResidualRotationData, residual_rotation_data
from .rotation_programming import RotationProgram, program_rotation


@dataclass(frozen=True)
class FourRowMask:
    """Literal selected mask and quadratic supports, ordered X then Z."""

    gates: Word
    quadratic_supports: tuple[tuple[int, ...], tuple[int, ...]]

    @property
    def quadratic_axis_count(self) -> int:
        return sum(bool(support) for support in self.quadratic_supports)

    @property
    def t_count(self) -> int:
        return 7 * self.quadratic_axis_count


@dataclass(frozen=True)
class NativeFourRowRotationTable:
    """Four active Ry/Rz rows and exact inactive identity on every input.

    ``addresses`` is (low,high); row index is low+2*high. ``mask_gates``
    acts only on core and addresses. The full-operator certificate covers
    arbitrary target, core, signal, address, enable, and reference inputs.
    """

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    addresses: tuple[int, int]
    enable: int
    enable_value: int
    gates: Word
    programs: tuple[RotationProgram, ...]
    operator_error_bound: Fraction
    mask_gates: Word
    quadratic_supports: tuple[tuple[int, ...], tuple[int, ...]]
    t_count: int

    @property
    def quadratic_axis_count(self) -> int:
        return sum(bool(support) for support in self.quadratic_supports)

    @property
    def nqubits(self) -> int:
        return self.q + 6

    @property
    def dirty_qubits(self) -> int:
        return self.q + 2


@dataclass(frozen=True)
class ResidualNativeFourRowTable:
    """Four U(z) completions, retaining the same outer factor twice.

    ``programs`` is ordered by Euler factor, then by low+2*high row.
    The first and last tuples and rotation objects are identical. Error
    includes coefficient preprocessing and all three native factors.
    """

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    addresses: tuple[int, int]
    enable: int
    enable_value: int
    gates: Word
    rotation_data: tuple[ResidualRotationData, ...]
    programs: tuple[tuple[RotationProgram, ...], ...]
    rotations: tuple[NativeFourRowRotationTable, ...]
    operator_error_bound: Fraction
    t_count: int

    @property
    def nqubits(self) -> int:
        return self.q + 6

    @property
    def dirty_qubits(self) -> int:
        return self.q + 2


def _four(values, name):
    try:
        rows = tuple(values)
    except TypeError as error:
        raise TypeError(f"{name} must contain exactly four entries.") from error
    if len(rows) != 4:
        raise ValueError(f"{name} must contain exactly four entries.")
    return rows


def _controlled_axis(axis: str, control: int, target: int) -> Word:
    if axis == "X":
        return (_gate("CX", control, target),)
    return (_gate("H", target), _gate("CX", control, target), _gate("H", target))


def _quadratic_axis(axis: str, support: tuple[int, ...], addresses: tuple[int, int]) -> Word:
    """Select an entire same-axis support with one exact seven-T center."""
    if not support:
        return ()
    pivot = support[0]
    # C X_p C-dagger=X_support for fanout, and
    # C Z_p C-dagger=Z_support for parity fanin. These CNOTs commute.
    conjugator = tuple(_gate("CX", pivot, target) if axis == "X"
                       else _gate("CX", target, pivot) for target in support[1:])
    center = _toffoli(*addresses, pivot)
    if axis == "Z":
        center = (_gate("H", pivot),) + center + (_gate("H", pivot),)
    return _inverse(conjugator) + center + conjugator


def _four_row_mask(signs, addresses: tuple[int, int]) -> FourRowMask:
    """Select four literal sign masks; useful independently at small size.

    Each row contains the same positive even number of exact flip bits.
    The addresses are distinct nonnegative wires beyond this mask's core.
    Row order is low+2*high. No signal or rotation target is required.
    """
    signs = tuple(tuple(row) for row in _four(signs, "sign rows"))
    length = len(signs[0])
    if length == 0 or length % 2 or any(len(row) != length for row in signs):
        raise ValueError("Sign rows must have equal positive even length.")
    if any(isinstance(bit, bool) or not isinstance(bit, Integral) or bit not in (0, 1)
           for row in signs for bit in row):
        raise ValueError("Sign rows must contain exact integer flip bits, not booleans.")
    addresses = _pair(addresses, "addresses")
    if (any(isinstance(wire, bool) or not isinstance(wire, Integral) or wire < length // 2
            for wire in addresses) or addresses[0] == addresses[1]):
        raise ValueError("Addresses must be distinct integer wires outside the mask core.")
    addresses = tuple(int(wire) for wire in addresses)
    masks = tuple(_mask(row) for row in signs)
    gates = []
    quadratic_supports = []
    # Keep the complete X part before the complete Z part, including the
    # quadratic terms. A mixed-axis shortcut can discard a row phase.
    for axis in ("X", "Z"):
        row0, row1, row2, row3 = (
            {wires[0] for name, wires in mask if name == axis} for mask in masks)
        constant, low, high = row0, row0 ^ row1, row0 ^ row2
        quadratic = tuple(sorted(row0 ^ row1 ^ row2 ^ row3))
        quadratic_supports.append(quadratic)
        gates.extend(_gate(axis, target) for target in sorted(constant))
        for address, support in zip(addresses, (low, high)):
            for target in sorted(support):
                gates.extend(_controlled_axis(axis, address, target))
        gates.extend(_quadratic_axis(axis, quadratic, addresses))
    return FourRowMask(tuple(gates), tuple(quadratic_supports))


def emit_four_row_rotation_table(
    programs: tuple[RotationProgram, ...], axis: str = "y", *, enable_value: int = 1,
) -> NativeFourRowRotationTable:
    """Emit a bounded two-address-bit table from four certified programs.

    Layout: core 0..q, signal q+1, target q+2, address low/high q+3/q+4,
    enable q+5. Programs must be genuine certified outputs at the same q;
    their common-angle enclosure promises remain external. Negative
    enable uses X wrappers around the complete amplified rotation.

    For k nonempty quadratic axis supports, the mask uses 7*k T/TDG
    gates. Its ten occurrences in the amplified rotation give the exact
    literal count 180*q+210+70*k. No added lookup work is used.
    """
    programs = _four(programs, "programs")
    enable_value = _enable_literal(enable_value)
    for program in programs:
        _validate_program(program)
    q = programs[0].q
    if any(program.q != q for program in programs):
        raise ValueError("All four rotation programs must have the same q.")
    core, signal, target = tuple(range(q + 1)), q + 1, q + 2
    addresses, enable = (q + 3, q + 4), q + 5
    mask = _four_row_mask(tuple(program.signs for program in programs), addresses)
    gates = _rotation_gates(q, mask.gates, axis, enable=enable)
    if enable_value == 0:
        wrapper = (_gate("X", enable),)
        gates = wrapper + gates + wrapper
    t_count = sum(name in ("T", "TDG") for name, _ in gates)
    if t_count != 180 * q + 210 + 70 * mask.quadratic_axis_count:
        raise ArithmeticError("The literal four-row rotation T count disagrees with its resource identity.")
    return NativeFourRowRotationTable(
        q, core, signal, target, addresses, enable, enable_value, gates, programs,
        max(program.operator_error_bound for program in programs),
        mask.gates, mask.quadratic_supports, t_count)


def emit_four_row_residual_table(
    coefficients: tuple[tuple[Fraction | int, Fraction | int], ...],
    q: int, *, enable_value: int = 1,
) -> ResidualNativeFourRowTable:
    """Emit four residual U(z) rows with exact inactive identity.

    Each (real,imag) is exact dyadic data with the approximation promise
    of ``residual_rotation_data``. The active full-operator error is
    below 130*2**(-q), including preprocessing. Active U(0) remains the
    pi/2 Ry completion, not an inactive identity row.

    If k_outer and k_middle count nonempty quadratic mask axes, the exact
    T count is 540*q+630+70*(2*k_outer+k_middle), at most 540*q+1050.
    All q+2 precision core/signal wires may start arbitrarily; the mask
    needs no additional helper and returns both addresses literally.
    """
    enable_value = _enable_literal(enable_value)
    if isinstance(q, bool) or not isinstance(q, Integral) or q < 5:
        raise ValueError("q must be an integer at least five, not a boolean.")
    q = int(q)
    coefficients = tuple(_pair(row, "coefficient row")
                         for row in _four(coefficients, "coefficients"))
    if any(isinstance(value, bool) or not isinstance(value, (Fraction, Integral))
           for row in coefficients for value in row):
        raise TypeError("Coefficient entries must be exact dyadic Fractions or integers, not booleans.")
    data = tuple(residual_rotation_data(real, imag, q) for real, imag in coefficients)
    outer_programs = tuple(program_rotation(*row.rotations[0], q) for row in data)
    middle_programs = tuple(program_rotation(*row.rotations[1], q) for row in data)
    outer = emit_four_row_rotation_table(outer_programs, "z", enable_value=enable_value)
    middle = emit_four_row_rotation_table(middle_programs, "y", enable_value=enable_value)
    programs = outer_programs, middle_programs, outer_programs
    rotations = outer, middle, outer
    error = max(row.completion_error_bound for row in data) + sum(
        (rotation.operator_error_bound for rotation in rotations), Fraction(0))
    return ResidualNativeFourRowTable(
        q, outer.core, outer.signal, outer.target, outer.addresses, outer.enable, enable_value,
        outer.gates + middle.gates + outer.gates, data, programs, rotations, error,
        sum(rotation.t_count for rotation in rotations))
