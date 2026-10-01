"""A bounded one-qubit residual state word with two logical flags.

The caller supplies dyadic approximations to a residual root ``a`` and
scaled tail ``w`` satisfying |a|**2 + |w|**2/2 = 1 and |w| <= 1. The
target state is a|0> + (w/sqrt(2))|1>. Two native residual tables give
success amplitude exactly one half in the ideal circuit; one exact
amplification step therefore returns both logical flags to zero.

The precision core and signal start arbitrarily, including entanglement
with a reference. The target flag, system, and mode flag start in zero.
Every emitted gate is elementary Clifford+T. No predicate helper,
measurement, reset, numerical angle solver, or clean precision work is
used. This component handles n=1; it is not a general addressed compiler.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .native_residual_rotation import Word, _gate, _inverse, _toffoli
from .native_residual_table import ResidualNativeTable, _pair, emit_residual_table
from .residual_table_preprocessing import _dyadic


@dataclass(frozen=True)
class NativeStateStage:
    """One chronological stage, retaining its literal gate word."""

    name: str
    gates: Word


@dataclass(frozen=True)
class ResidualStateCertificate:
    """Conditional analytic bounds and exact resources of the literal word.

    The input approximation and normalization promises remain external.
    ``q_operator_error_bound`` bounds Q on its complete input space.
    ``state_error_bound`` includes the three Q/Q-dagger occurrences and
    bounds the preparation-vector error on zero logical inputs, including
    arbitrary dirty work and external references. It also bounds the
    operator error from the corresponding ideal amplified unitary.
    """

    input_error_bound: Fraction
    normalization_tolerance: Fraction
    normalization_discrepancy: Fraction
    q_operator_error_bound: Fraction
    state_error_bound: Fraction
    t_count: int
    nqubits: int
    dirty_qubits: int


@dataclass(frozen=True)
class ResidualNativeState:
    """The literal n=1 residual preparation and its table certificates.

    Layout: core 0..q, signal q+1, target flag q+2, system q+3, mode
    flag q+4. Only the last three logical wires have specified zero input.
    ``q_inverse_gates`` is the actual reversed, daggered ``q_gates``.
    The full flattened ``gates`` preserves the final global minus sign.
    """

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    system: int
    mode: int
    coefficients: tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]
    tables: tuple[ResidualNativeTable, ResidualNativeTable]
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


def _controlled_h(control: int, target: int) -> Word:
    """Exact controlled-H using two T gates and no discarded phase."""
    # The matrix of this chronological word is exp(i*pi/8)*Ry(pi/8).
    # Conjugating Z gives H; the two scalar phases cancel exactly.
    basis = (_gate("SDG", target), _gate("H", target), _gate("T", target),
             _gate("H", target), _gate("S", target))
    cz = (_gate("H", target), _gate("CX", control, target), _gate("H", target))
    return _inverse(basis) + cz + basis


def _good_reflection(mode: int, target: int) -> Word:
    """I-2|00><00| on the two logical flags."""
    negative = (_gate("X", mode), _gate("X", target))
    cz = (_gate("H", target), _gate("CX", mode, target), _gate("H", target))
    return negative + cz + _inverse(negative)


def _initial_reflection(system: int, mode: int, target: int) -> Word:
    """I-2|000><000| on logical inputs, with exact seven-T CCZ."""
    negative = (_gate("X", system), _gate("X", mode), _gate("X", target))
    ccz = ((_gate("H", target),) + _toffoli(system, mode, target)
           + (_gate("H", target),))
    return negative + ccz + _inverse(negative)


def emit_residual_state(
    a: tuple[Fraction | int, Fraction | int],
    w: tuple[Fraction | int, Fraction | int],
    q: int,
) -> ResidualNativeState:
    """Emit the bounded residual state component from exact dyadic data.

    External promise: the pairs approximate complex a,w with Euclidean
    error at most delta=2**(-2*q-20) each, |a|**2+|w|**2/2=1, and |w|<=1.
    The desired state is a|0>+(w/sqrt(2))|1>. These true coefficients are
    not reconstructed or normalized by this function.

    Checks cover exact types, dyadic denominators, q>=5, consistency with
    the unit disks, and the necessary normalization condition
    abs(|a_input|**2+|w_input|**2/2-1)<=3*delta+3*delta**2/2.
    Passing these checks does not establish the external promise.

    The state error is three times the maximum of the two table operator
    bounds, strictly below 390*2**(-q). Exact inactive identity makes this
    a maximum, including coherent mode inputs. The literal T/TDG count
    is 3240*q+3793, on q+5 wires of which q+2 are arbitrary dirty work.
    """
    if isinstance(q, bool) or not isinstance(q, Integral) or q < 5:
        raise ValueError("q must be an integer at least five, not a boolean.")
    q = int(q)
    inputs = (_pair(a, "a"), _pair(w, "w"))
    if any(isinstance(value, bool) or not isinstance(value, (Fraction, Integral))
           for pair in inputs for value in pair):
        raise TypeError("Coefficient entries must be exact dyadic Fractions or integers, not booleans.")
    coefficients = tuple(tuple(_dyadic(value) for value in pair) for pair in inputs)
    a_data, w_data = coefficients
    delta = Fraction(1, 1 << (2 * q + 20))
    tolerance = 3 * delta + Fraction(3, 2) * delta * delta
    norm = sum((value * value for value in a_data), Fraction(0)) + sum(
        (value * value for value in w_data), Fraction(0)) / 2
    discrepancy = abs(norm - 1)
    if discrepancy > tolerance:
        raise ValueError("The approximations are inconsistent with the promised state normalization.")

    m0 = emit_residual_table((a_data, a_data), q, enable_value=0)
    m1 = emit_residual_table(((Fraction(0), Fraction(0)), w_data), q, enable_value=1)
    system, mode, target = m0.address, m0.enable, m0.target
    controlled_h = _controlled_h(mode, system)
    mode_h = (_gate("H", mode),)
    q_gates = mode_h + controlled_h + m0.gates + m1.gates + mode_h
    q_inverse = _inverse(q_gates)
    good = _good_reflection(mode, target)
    initial = _initial_reflection(system, mode, target)
    # Chronological X,Z,X,Z is literal -I on every target input.
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
    if t_count != 3240 * q + 3793:
        raise ArithmeticError("The literal residual-state T count disagrees with its resource identity.")
    certificate = ResidualStateCertificate(
        delta, tolerance, discrepancy, error, 3 * error, t_count, q + 5, q + 2)
    return ResidualNativeState(
        q, m0.core, m0.signal, target, system, mode, coefficients, (m0, m1),
        controlled_h, good, initial, minus, q_gates, q_inverse, stages, gates, certificate)
