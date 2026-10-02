"""A bounded two-qubit residual state word using four-row native tables.

The caller supplies dyadic approximations to a root a and three scaled
tails w_j=2*phi_j, with |a|**2+sum_j|w_j|**2/4=1 and |w_j|<=1. Two
native tables give ideal success amplitude one half. One exact amplitude
amplification step prepares a|00>+sum_{j=1}^3(w_j/2)|j> and returns the
two logical flags to zero, retaining the literal state phase.

The core and synthesis signal may have arbitrary inputs entangled with
an external reference. The initial-subspace reflection briefly borrows
core[0], returning it exactly on every input, without adding a wire.
Only the two system qubits and two logical flags start in zero. There is
no measurement, reset, clean precision work, or numerical angle solver.
This handles n=2 and does not implement general addressed preparation.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .native_residual_lookup import ResidualNativeFourRowTable, emit_four_row_residual_table
from .native_residual_rotation import Word, _gate, _inverse, _toffoli
from .native_residual_state import (
    NativeStateStage, ResidualStateCertificate, _controlled_h, _good_reflection,
)
from .native_residual_table import _pair
from .residual_table_preprocessing import _dyadic


@dataclass(frozen=True)
class TwoQubitResidualNativeState:
    """Literal n=2 preparation, with its actual inverse and certificates.

    Layout: core 0..q, signal q+1, target flag q+2, system low/high
    q+3/q+4, mode flag q+5. ``systems`` is ordered (low,high), and the
    integer basis index is low+2*high. ``reflection_helper`` aliases
    core[0]; it is neither initialized nor an additional work qubit.
    """

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    systems: tuple[int, int]
    mode: int
    reflection_helper: int
    coefficients: tuple[tuple[Fraction, Fraction], ...]
    tables: tuple[ResidualNativeFourRowTable, ResidualNativeFourRowTable]
    controlled_h_gates: Word
    good_reflection_gates: Word
    initial_reflection_gates: Word
    global_minus_gates: Word
    q_gates: Word
    q_inverse_gates: Word
    stages: tuple[NativeStateStage, ...]
    gates: Word
    certificate: ResidualStateCertificate

    @property
    def flag_t(self) -> int:
        return self.target

    @property
    def flag_s(self) -> int:
        return self.mode

    @property
    def nqubits(self) -> int:
        return self.certificate.nqubits

    @property
    def dirty_qubits(self) -> int:
        return self.certificate.dirty_qubits

    @property
    def t_count(self) -> int:
        return self.certificate.t_count

    @property
    def q_operator_error_bound(self) -> Fraction:
        return self.certificate.q_operator_error_bound

    @property
    def state_error_bound(self) -> Fraction:
        return self.certificate.state_error_bound

    @property
    def operator_error_bound(self) -> Fraction:
        """Error from the ideal amplified unitary, not a state reset map."""
        return self.certificate.state_error_bound


def _two_qubit_initial_reflection(
    systems: tuple[int, int], mode: int, target: int, helper: int,
) -> Word:
    """I-2|0000><0000| on logical inputs, independent of the helper.

    F toggles helper by x_low*x_high; G toggles target by helper*mode.
    Chronological F,G,F,G returns arbitrary helper and toggles target by
    x_low*x_high*mode. The Hadamards give the exact four-qubit phase
    reflection. Each of the four literal Toffolis costs seven T gates.
    """
    negative = tuple(_gate("X", wire) for wire in (*systems, mode, target))
    f = _toffoli(*systems, helper)
    g = _toffoli(helper, mode, target)
    h = (_gate("H", target),)
    return negative + h + f + g + f + g + h + _inverse(negative)


def emit_two_qubit_residual_state(
    a: tuple[Fraction | int, Fraction | int],
    tails: tuple[tuple[Fraction | int, Fraction | int], ...],
    q: int,
) -> TwoQubitResidualNativeState:
    """Emit the bounded n=2 residual preparation from exact dyadic pairs.

    External promise: the root and three tail pairs approximate a,w_j
    with Euclidean error at most delta=2**(-2*q-20) each, |w_j|<=1, and
    |a|**2+sum_j|w_j|**2/4=1. Tail order is j=1,2,3, where j=low+2*high.
    The prepared state is a|00>+sum_j(w_j/2)|j>. No coefficient is
    reconstructed or renormalized by this routine.

    Checks cover exact types, dyadic denominators, q>=5, unit-disk
    consistency, and the necessary normalization discrepancy bound
    (7/2)*delta+(7/4)*delta**2. Passing does not certify the external
    evaluation or normalization promises.

    The full state error is three times the maximum of the two table
    bounds, strictly below 390*2**(-q), including arbitrary dirty work
    and references. The word uses q+6 wires, with q+2 arbitrary dirty
    wires. If k_z,k_y count the tail table's nonempty quadratic mask
    axes, its exact T count is 3240*q+3820+210*(2*k_z+k_y), at most
    3240*q+5080. The repeated root table has no quadratic support.
    """
    if isinstance(q, bool) or not isinstance(q, Integral) or q < 5:
        raise ValueError("q must be an integer at least five, not a boolean.")
    q = int(q)
    try:
        tails = tuple(tails)
    except TypeError as error:
        raise TypeError("tails must contain exactly three coefficient pairs.") from error
    if len(tails) != 3:
        raise ValueError("tails must contain exactly three coefficient pairs.")
    inputs = (_pair(a, "a"),) + tuple(_pair(tail, "tail coefficient") for tail in tails)
    if any(isinstance(value, bool) or not isinstance(value, (Fraction, Integral))
           for pair in inputs for value in pair):
        raise TypeError("Coefficient entries must be exact dyadic Fractions or integers, not booleans.")
    coefficients = tuple(tuple(_dyadic(value) for value in pair) for pair in inputs)
    a_data, tail_data = coefficients[0], coefficients[1:]
    delta = Fraction(1, 1 << (2 * q + 20))
    tolerance = Fraction(7, 2) * delta + Fraction(7, 4) * delta * delta
    norm = sum((value * value for value in a_data), Fraction(0)) + sum(
        (value * value for tail in tail_data for value in tail), Fraction(0)) / 4
    discrepancy = abs(norm - 1)
    if discrepancy > tolerance:
        raise ValueError("The approximations are inconsistent with the promised state normalization.")

    m0 = emit_four_row_residual_table((a_data,) * 4, q, enable_value=0)
    m1 = emit_four_row_residual_table(
        ((Fraction(0), Fraction(0)),) + tail_data, q, enable_value=1)
    systems, mode, target = m0.addresses, m0.enable, m0.target
    helper = m0.core[0]
    controlled_h = _controlled_h(mode, systems[0]) + _controlled_h(mode, systems[1])
    mode_h = (_gate("H", mode),)
    q_gates = mode_h + controlled_h + m0.gates + m1.gates + mode_h
    q_inverse = _inverse(q_gates)
    good = _good_reflection(mode, target)
    initial = _two_qubit_initial_reflection(systems, mode, target, helper)
    # Chronological X,Z,X,Z retains the global -I on the complete space.
    minus = (_gate("X", target), _gate("Z", target),
             _gate("X", target), _gate("Z", target))
    stages = (
        NativeStateStage("q", q_gates), NativeStateStage("good_reflection", good),
        NativeStateStage("q_inverse", q_inverse), NativeStateStage("initial_reflection", initial),
        NativeStateStage("q", q_gates), NativeStateStage("global_minus", minus),
    )
    gates = tuple(gate for stage in stages for gate in stage.gates)
    error = max(m0.operator_error_bound, m1.operator_error_bound)
    t_count = sum(name in ("T", "TDG") for name, _ in gates)
    k_z, k_y = (rotation.quadratic_axis_count for rotation in m1.rotations[:2])
    if t_count != 3240 * q + 3820 + 210 * (2 * k_z + k_y):
        raise ArithmeticError("The literal two-qubit residual-state T count disagrees with its resource identity.")
    certificate = ResidualStateCertificate(
        delta, tolerance, discrepancy, error, 3 * error, t_count, q + 6, q + 2)
    return TwoQubitResidualNativeState(
        q, m0.core, m0.signal, target, systems, mode, helper, coefficients, (m0, m1),
        controlled_h, good, initial, minus, q_gates, q_inverse, stages, gates, certificate)
