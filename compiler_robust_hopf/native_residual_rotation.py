"""Unaddressed native rotations from certified residual coefficients.

The chronological words use only elementary Clifford+T gates and no clean
work.  The precision core and synthesis signal may have arbitrary inputs.
Their approximate return is included in the full-operator certificate;
there is no initialization, measurement, reset, or source-state input.

This is one residual SU(2) row, not an addressed table or a state compiler.
Word construction uses exact classical data and linear storage in q, with
no matrix construction or numerical trigonometry.  See
docs/ONE_CLEAN_COMPILER.md sections 2--4 and 9.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .residual_table_preprocessing import ResidualRotationData, residual_rotation_data
from .rotation_programming import RotationProgram, program_rotation


Gate = tuple[str, tuple[int, ...]]
Word = tuple[Gate, ...]


@dataclass(frozen=True)
class NativeRotationWord:
    """A literal Ry/Rz word; qubit zero is the low integer bit.

    ``certificate.operator_error_bound`` bounds the operator norm from the
    intended target rotation tensored with identity on ``core`` and
    ``signal``.  No predicate helper is needed for this unaddressed word.
    """

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    gates: Word
    certificate: RotationProgram


@dataclass(frozen=True)
class ResidualNativeWord:
    """Chronological Rz, Ry, Rz completion of the supplied disk coefficient.

    The input coefficient promise is that of ``residual_rotation_data``.
    ``rotations[0] is rotations[2]``: the same literal half-phase word is
    used twice, including its program and circuit phase.
    """

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    gates: Word
    rotation_data: ResidualRotationData
    programs: tuple[RotationProgram, RotationProgram, RotationProgram]
    rotations: tuple[NativeRotationWord, NativeRotationWord, NativeRotationWord]
    operator_error_bound: Fraction


def _gate(name: str, *qubits: int) -> Gate:
    return name, tuple(qubits)


def _inverse(word: Word) -> Word:
    inverse = {"T": "TDG", "TDG": "T", "S": "SDG", "SDG": "S"}
    return tuple((inverse.get(name, name), wires) for name, wires in reversed(word))


def _plane(left: int, right: int, left_axis: str, right_axis: str, sign: int) -> Word:
    # The scalar of this Pauli-T plane cancels in U X U-dagger and in
    # its controlled versions. The word and its inverse remain literal.
    basis: Word = ()
    for qubit, axis in ((left, left_axis), (right, right_axis)):
        basis += ((_gate("SDG", qubit), _gate("H", qubit)) if axis == "Y"
                  else (_gate("H", qubit),))
    parity = (_gate("CX", left, right),)
    return (basis + parity + (_gate("T" if sign == 1 else "TDG", right),)
            + parity + _inverse(basis))


def _loader(q: int) -> Word:
    gates = [_gate("T", 0)]
    gates.extend(_plane(0, 1, "Y", "X", -1))
    for qubit in range(q - 1):
        gates.extend(_plane(qubit, qubit + 1, "X", "Y", 1))
    for qubit in range(1, q):
        gates.extend(_plane(qubit, qubit + 1, "Y", "X", -1))
    return tuple(gates)


def _mask(signs: tuple[int, ...]) -> Word:
    xs, zs = [], []
    prefix = 0
    for qubit in range(len(signs) // 2):
        first, second = signs[2 * qubit:2 * qubit + 2]
        xbit, zbit = first ^ second, first ^ prefix
        if xbit:
            xs.append(_gate("X", qubit))
        if zbit:
            zs.append(_gate("Z", qubit))
        prefix ^= xbit
    return tuple(xs + zs)


def _scalar(loader: Word, mask: Word, signal: int) -> Word:
    inverse = _inverse(loader)
    center = (_gate("CX", signal, 0),)
    source = inverse + (_gate("X", 0),) + loader
    control_one = inverse + center + loader
    control_zero = (inverse + (_gate("X", signal),) + center
                    + (_gate("X", signal),) + loader)
    return ((_gate("H", signal),) + control_one + _inverse(mask) + source
            + mask + control_zero + (_gate("H", signal),))


def emit_rotation(program: RotationProgram, axis: str = "y") -> NativeRotationWord:
    """Emit the certified unaddressed rotation, on all work input states.

    Axes are ``'y'`` and ``'z'``, with R_axis = cos(theta) I minus
    i sin(theta) Pauli_axis.  ``program`` must come from
    ``program_rotation``; its certified unit-circle input promise applies.
    The layout is core 0..q, signal q+1, target q+2.  Its q+2 dirty wires
    exclude the target, which may itself be occupied compiler work.
    """
    if not isinstance(program, RotationProgram):
        raise TypeError("program must be a RotationProgram.")
    if axis not in ("y", "z"):
        raise ValueError("axis must be 'y' or 'z'.")
    certified = program_rotation(program.cosine, program.sine, program.q)
    if program != certified:
        raise ValueError("program must match its certified coefficient encoding.")
    q = program.q
    core, signal, target = tuple(range(q + 1)), q + 1, q + 2
    loader = _loader(q)
    fixed = tuple(int(2 <= index <= 2 * q and index % 2 == 0)
                  for index in range(2 * (q + 1)))
    scalar_fixed = _scalar(loader, _mask(fixed), signal)
    scalar_program = _scalar(loader, _mask(program.signs), signal)
    cz = (_gate("H", target), _gate("CX", signal, target), _gate("H", target))
    cx = (_gate("CX", signal, target),)
    left, middle = cz + scalar_fixed + cz, cx + scalar_program + cx
    primitive = left + middle + _inverse(left)
    inverse, z = _inverse(primitive), (_gate("Z", signal),)
    # Four reflections R=-Z have no leftover scalar minus sign.
    gates = primitive + z + inverse + z + primitive + z + inverse + z + primitive
    if axis == "z":
        # B=H S-dagger maps Y to Z. Emit B-dagger, Ry, B chronologically.
        gates = ((_gate("H", target), _gate("S", target)) + gates
                 + (_gate("SDG", target), _gate("H", target)))
    return NativeRotationWord(q, core, signal, target, gates, program)


def emit_residual_row(real: Fraction | int, imag: Fraction | int, q: int) -> ResidualNativeWord:
    """Emit one literal U(z) completion from the certified dyadic input.

    The caller certifies |real+i*imag-z| <= 2**(-2*q-20), |z|<=1.
    The full-operator error includes preprocessing and all three native
    rotations, on arbitrary target/core/signal/reference inputs.  It is
    below 130*2**(-q); no global or branch-relative phase is discarded.
    """
    data = residual_rotation_data(real, imag, q)
    outer_program = program_rotation(*data.rotations[0], q)
    middle_program = program_rotation(*data.rotations[1], q)
    outer = emit_rotation(outer_program, "z")
    middle = emit_rotation(middle_program, "y")
    programs = outer_program, middle_program, outer_program
    rotations = outer, middle, outer
    error = data.completion_error_bound + sum(
        (program.operator_error_bound for program in programs), Fraction(0))
    return ResidualNativeWord(
        q, outer.core, outer.signal, outer.target,
        outer.gates + middle.gates + outer.gates, data, programs, rotations, error)
