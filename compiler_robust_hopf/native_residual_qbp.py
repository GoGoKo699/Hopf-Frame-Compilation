"""A bounded native common-coarse QBP fixture with fine residual synthesis.

One system qubit uses the fixed gauge W'=Rz(pi/4+v) Ry(pi/4+u), where
u=2*atan(1/1024) and v=2*atan(1/256). The actual recorded coarse word is
C=Rz(pi/4) Ry(pi/4). All coefficient construction and histogram decoding
below use exact rational arithmetic; angles are never evaluated by the
emitter. The controlled observable is the literal calibrated Hadamard.

The magnitude stream prepares the common-C pair, then applies controlled
H, actual C-dagger, system H, and an X/Y branch readout. The independent
phase stream uses unbranched residual preparation, common C, plus branch,
controlled H, and Y readout. Its word has no inverse C. All intermediate
flag/dirty leakage is retained. This finite fixture is not a general Hopf
compiler, certified trigonometric evaluator, or native fine-frame emitter.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .native_branched_residual_state import BranchedResidualNativeState, emit_branched_residual_state
from .native_residual_rotation import Word, _gate, _inverse
from .native_residual_state import NativeStateStage, ResidualNativeState, _controlled_h, emit_residual_state
from .residual_table_preprocessing import CoefficientInterval, sqrt_enclosure


ANGLE_DENOMINATOR = 1024
PHASE_DENOMINATOR = 256
_Pair = tuple[Fraction, Fraction]


def _half_angle(denominator: int) -> tuple[Fraction, Fraction]:
    return (Fraction(denominator * denominator - 1, denominator * denominator + 1),
            Fraction(2 * denominator, denominator * denominator + 1))


_CU, _SU = _half_angle(ANGLE_DENOMINATOR)
_CV, _SV = _half_angle(PHASE_DENOMINATOR)
_ROOT = (_CV * _CU, _SV * _SU)
_TAIL = (_CV * _SU, _SV * _CU)
MAGNITUDE_X_WEIGHTS = (4 * (_ROOT[0] - _TAIL[0]), -4 * (_ROOT[0] + _TAIL[0]))
MAGNITUDE_Y_WEIGHTS = (4 * (_TAIL[1] - _ROOT[1]), 4 * (_ROOT[1] + _TAIL[1]))


@dataclass(frozen=True)
class ResidualQBPData:
    """Exact fixture geometry and certified dyadic synthesis inputs.

    Complex values are (real,imag) Fraction pairs. ``residual_tail`` is
    b, while ``programmed_scaled_tail`` approximates sqrt(2)*b. The two
    squared-error bounds certify the root and scaled-tail input promises.
    Target state and magnitude derivative retain the fixed mean-zero gauge.
    """

    q: int
    angle_denominator: int
    phase_denominator: int
    cos_u: Fraction
    sin_u: Fraction
    cos_v: Fraction
    sin_v: Fraction
    residual_root: _Pair
    residual_tail: _Pair
    target_state: tuple[_Pair, _Pair]
    magnitude_derivative: tuple[_Pair, _Pair]
    coarse_matrix: tuple[tuple[_Pair, _Pair], tuple[_Pair, _Pair]]
    coarse_distance_squared: Fraction
    coarse_distance_bound: Fraction
    programmed_root: _Pair
    programmed_scaled_tail: _Pair
    sqrt_two_interval: CoefficientInterval
    working_bits: int
    input_error_bound: Fraction
    coefficient_error_squared_bounds: tuple[Fraction, Fraction]
    magnitude_x_weights: tuple[Fraction, Fraction]
    magnitude_y_weights: tuple[Fraction, Fraction]


@dataclass(frozen=True)
class NativeQBPStream:
    """One complete chronological stream, without measurement or selection.

    Readouts are native basis changes; the caller measures system and
    branch at the end. Compiler flags and dirty wires are not postselected.
    ``preparation_error_bound`` covers arbitrary dirty/reference inputs;
    subsequent wrappers and the controlled observable are exact.
    """

    name: str
    stages: tuple[NativeStateStage, ...]
    gates: Word
    t_count: int
    preparation_error_bound: Fraction


@dataclass(frozen=True)
class ResidualNativeQBP:
    """Recorded actual coarse circuit, residual emitters, and both streams.

    The unbranched phase emitter has its mode wire remapped from q+4 to
    q+5 in ``phase_residual_gates``. Thus physical branch q+4 is available
    throughout, and both streams share the same q+6-wire layout. The
    separate streams reuse q+2 arbitrary dirty core/signal wires.
    """

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    system: int
    branch: int
    mode: int
    data: ResidualQBPData
    pair_residual: BranchedResidualNativeState
    phase_residual: ResidualNativeState
    phase_residual_gates: Word
    coarse_gates: Word
    coarse_inverse_gates: Word
    controlled_observable_gates: Word
    pair_preparation_gates: Word
    phase_preparation_gates: Word
    magnitude_x: NativeQBPStream
    magnitude_y: NativeQBPStream
    phase: NativeQBPStream

    @property
    def nqubits(self) -> int:
        return self.q + 6

    @property
    def dirty_qubits(self) -> int:
        return self.q + 2

    @property
    def flag_t(self) -> int:
        return self.target

    @property
    def flag_s(self) -> int:
        return self.mode

    @property
    def pair_preparation_t_count(self) -> int:
        return self.pair_residual.t_count + 2

    @property
    def phase_preparation_t_count(self) -> int:
        return self.phase_residual.t_count + 2


def _precision(q: int) -> int:
    if isinstance(q, bool) or not isinstance(q, Integral) or q < 5:
        raise ValueError("q must be an integer at least five, not a boolean.")
    return int(q)


def _round_dyadic(value: Fraction, bits: int) -> Fraction:
    """Round a rational to the nearest dyadic; ties go toward +infinity."""
    scaled = value * (1 << bits)
    integer = (2 * scaled.numerator + scaled.denominator) // (2 * scaled.denominator)
    return Fraction(integer, 1 << bits)


def residual_qbp_data(q: int) -> ResidualQBPData:
    """Evaluate this fixed fixture exactly and certify both dyadic inputs.

    With B=2*q+24, enclose sqrt(2) on the B-bit grid, multiply its
    midpoint by the exact rational b, and round both complex coordinates.
    Endpoint distances give a rational upper bound for each squared
    Euclidean coefficient error. No floating-point values enter a word.
    """
    q = _precision(q)
    bits = 2 * q + 24
    root_two = sqrt_enclosure(Fraction(2), bits)
    midpoint = (root_two.lower + root_two.upper) / 2
    root = tuple(_round_dyadic(value, bits) for value in _ROOT)
    tail = tuple(_round_dyadic(midpoint * value, bits) for value in _TAIL)
    root_error_squared = sum(((rounded - exact) ** 2
                              for rounded, exact in zip(root, _ROOT)), Fraction(0))
    tail_error_squared = sum((max(abs(rounded - root_two.lower * exact),
                                  abs(rounded - root_two.upper * exact)) ** 2
                              for rounded, exact in zip(tail, _TAIL)), Fraction(0))
    delta = Fraction(1, 1 << (2 * q + 20))
    if max(root_error_squared, tail_error_squared) > delta * delta:
        raise ArithmeticError("The fixture coefficient evaluator failed its dyadic error certificate.")
    root_norm = sum((value * value for value in _ROOT), Fraction(0))
    tail_norm = sum((value * value for value in _TAIL), Fraction(0))
    if root_norm + tail_norm != 1 or 2 * tail_norm > 1:
        raise ArithmeticError("The fixture residual violates normalization or the scaled-tail unit disk.")
    coarse_distance_squared = 2 * (1 - _CU * _CV)
    if coarse_distance_squared >= Fraction(1, 64 ** 2):
        raise ArithmeticError("The recorded coarse circuit exceeds the promised fixture distance.")

    phase_minus, phase_plus = (_CV - _SV, -_CV - _SV), (_CV - _SV, _CV + _SV)
    state = (tuple((_CU - _SU) * value / 2 for value in phase_minus),
             tuple((_CU + _SU) * value / 2 for value in phase_plus))
    derivative = (tuple(-(_CU + _SU) * value / 2 for value in phase_minus),
                  tuple((_CU - _SU) * value / 2 for value in phase_plus))
    half = Fraction(1, 2)
    coarse = (((half, -half), (-half, half)), ((half, half), (half, half)))
    return ResidualQBPData(
        q, ANGLE_DENOMINATOR, PHASE_DENOMINATOR, _CU, _SU, _CV, _SV,
        _ROOT, _TAIL, state, derivative, coarse, coarse_distance_squared, Fraction(1, 64),
        root, tail, root_two, bits, delta, (root_error_squared, tail_error_squared),
        MAGNITUDE_X_WEIGHTS, MAGNITUDE_Y_WEIGHTS)


def _stream(name: str, stages: tuple[NativeStateStage, ...], error: Fraction) -> NativeQBPStream:
    gates = tuple(gate for stage in stages for gate in stage.gates)
    return NativeQBPStream(name, stages, gates, sum(name in ("T", "TDG") for name, _ in gates), error)


def emit_residual_qbp(q: int) -> ResidualNativeQBP:
    """Emit complete magnitude X/Y and phase streams for the fixed fixture.

    Initialize the system, two compiler flags, and protocol branch to
    zero at each execution. Core and signal inputs remain arbitrary and
    may be entangled with a reference. Both preparations have complete
    error below 390*2**(-q), including all work leakage.

    The exact literal T counts are pair-residual+6 for either magnitude
    basis and plain-residual+4 for phase. There are three charged common-C
    appearances across a magnitude/phase execution pair: C in each
    preparation and actual C-dagger in magnitude readout. This provides
    no native fine-frame comparison or optimality claim.
    """
    data = residual_qbp_data(q)
    q = data.q
    a, w = data.programmed_root, data.programmed_scaled_tail
    pair = emit_branched_residual_state((((1, 0), (0, 0)), (a, w)), q)
    plain = emit_residual_state(a, w, q)
    system, branch, mode = pair.system, pair.branch, pair.mode
    # Only the old mode wire moves; the external protocol branch stays free.
    phase_residual = tuple((name, tuple(mode if wire == plain.mode else wire for wire in wires))
                           for name, wires in plain.gates)
    # Z,H is Ry(pi/4); X,TDG,X,T is literal Rz(pi/4), including its scalar.
    coarse = tuple(_gate(name, system) for name in ("Z", "H", "X", "TDG", "X", "T"))
    coarse_inverse = _inverse(coarse)
    observable = _controlled_h(branch, system)
    branch_h, system_h = (_gate("H", branch),), (_gate("H", system),)
    readout_y = (_gate("SDG", branch), _gate("H", branch))
    pair_stages = (NativeStateStage("branch_plus", branch_h),
                   NativeStateStage("pair_residual", pair.gates),
                   NativeStateStage("coarse", coarse))
    phase_stages = (NativeStateStage("phase_residual", phase_residual),
                    NativeStateStage("coarse", coarse), NativeStateStage("branch_plus", branch_h))
    pair_preparation = tuple(gate for stage in pair_stages for gate in stage.gates)
    phase_preparation = tuple(gate for stage in phase_stages for gate in stage.gates)
    magnitude_stages = pair_stages + (
        NativeStateStage("controlled_observable", observable),
        NativeStateStage("coarse_inverse", coarse_inverse),
        NativeStateStage("system_walsh", system_h),
    )
    magnitude_x = _stream("magnitude_x", magnitude_stages + (
        NativeStateStage("branch_x_readout", branch_h),), pair.state_error_bound)
    magnitude_y = _stream("magnitude_y", magnitude_stages + (
        NativeStateStage("branch_y_readout", readout_y),), pair.state_error_bound)
    phase = _stream("phase", phase_stages + (
        NativeStateStage("controlled_observable", observable),
        NativeStateStage("branch_y_readout", readout_y),), plain.state_error_bound)
    if (magnitude_x.t_count != pair.t_count + 6 or magnitude_y.t_count != pair.t_count + 6
            or phase.t_count != plain.t_count + 4):
        raise ArithmeticError("A literal QBP stream disagrees with its complete T ledger.")
    return ResidualNativeQBP(
        q, pair.core, pair.signal, pair.target, system, branch, mode, data, pair, plain,
        phase_residual, coarse, coarse_inverse, observable, pair_preparation, phase_preparation,
        magnitude_x, magnitude_y, phase)


def _counts(values, name: str) -> tuple[int, int]:
    try:
        counts = tuple(values)
    except TypeError as error:
        raise TypeError(f"{name} must contain exactly two signed integer counts.") from error
    if len(counts) != 2:
        raise ValueError(f"{name} must contain exactly two signed integer counts.")
    if any(isinstance(value, bool) or not isinstance(value, Integral) for value in counts):
        raise TypeError("Histogram counts must be exact integers, not booleans or floats.")
    return tuple(int(value) for value in counts)


def _shots(shots: int) -> int:
    if isinstance(shots, bool) or not isinstance(shots, Integral) or shots <= 0:
        raise ValueError("shots must be a positive integer, not a boolean.")
    return int(shots)


def decode_magnitude_histograms(hist_x, hist_y, shots: int) -> Fraction:
    """Return the exact empirical raw magnitude derivative for this fixture.

    Each execution contributes its branch sign to one leaf in either
    hist_x or hist_y. Basis choices must be independent fair X/Y draws;
    shots counts all magnitude executions, including signed cancellations.
    The observable coefficient norm is one. No sampling guarantee or
    change to the fixed target/coarse coefficients is performed here.
    """
    hx, hy, shots = _counts(hist_x, "hist_x"), _counts(hist_y, "hist_y"), _shots(shots)
    if sum(abs(value) for value in hx + hy) > shots:
        raise ValueError("The total absolute magnitude histogram count cannot exceed shots.")
    return (sum((weight * count for weight, count in zip(MAGNITUDE_X_WEIGHTS, hx)), Fraction(0))
            + sum((weight * count for weight, count in zip(MAGNITUDE_Y_WEIGHTS, hy)), Fraction(0))) / shots


def decode_phase_histogram(hist, shots: int) -> tuple[Fraction, Fraction]:
    """Return both exact empirical raw leaf-phase coordinates, 2*hist/shots.

    Hist records branch signs at the measured leaf in the independent Y
    phase stream. Both raw coordinates are retained; their empirical sum
    is not projected to zero. Compiler flags and dirty inputs are never
    postselected when updating either histogram.
    """
    counts, shots = _counts(hist, "hist"), _shots(shots)
    if sum(abs(value) for value in counts) > shots:
        raise ValueError("The total absolute phase histogram count cannot exceed shots.")
    return tuple(Fraction(2 * count, shots) for count in counts)
