"""Coherent selection of two bounded one-qubit residual preparations.

An arbitrary protocol branch c selects phi_c=a_c|0>+(w_c/sqrt(2))|1>.
Both coefficient pairs are normalized separately, and each w_c lies in
the unit disk. Four-row native tables address the system and branch;
one exact amplification step prepares the selected state while returning
the two logical flags to zero. Literal branch-relative phases are kept.

The system and two flags start in zero. The protocol branch, precision
core, and synthesis signal may be arbitrary and jointly entangled with a
reference. The initial-subspace reflection excludes the protocol branch.
There is no branch initialization, measurement, reset, or extra helper.
The complete word preserves c; native lookup decompositions can change
it temporarily inside individual Toffolis. This component handles one
system qubit and two branches, not general addressed state preparation.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .native_residual_lookup import ResidualNativeFourRowTable, emit_four_row_residual_table
from .native_residual_rotation import Word, _gate, _inverse
from .native_residual_state import (
    NativeStateStage, _controlled_h, _good_reflection, _initial_reflection,
)
from .native_residual_table import _pair
from .residual_table_preprocessing import _dyadic


_CoefficientPair = tuple[Fraction, Fraction]
_BranchCoefficients = tuple[_CoefficientPair, _CoefficientPair]


@dataclass(frozen=True)
class BranchedResidualStateCertificate:
    """Conditional coherent-preparation bounds and exact literal resources.

    Normalization discrepancies are recorded separately for branches
    zero and one. Passing the finite checks does not certify the external
    coefficient-evaluation or normalization promises. The Q bound is a
    full-operator bound. The state bound covers arbitrary branch/work/
    reference inputs with zero system and flags, and also bounds the
    operator error from the corresponding ideal amplified unitary.
    """

    input_error_bound: Fraction
    normalization_tolerance: Fraction
    normalization_discrepancies: tuple[Fraction, Fraction]
    q_operator_error_bound: Fraction
    state_error_bound: Fraction
    t_count: int
    nqubits: int
    dirty_qubits: int


@dataclass(frozen=True)
class BranchedResidualNativeState:
    """The literal coherent branch word and both table certificates.

    Layout: core 0..q, signal q+1, target flag q+2, system q+3,
    protocol branch q+4, mode flag q+5. The table address order is
    (system,branch), with row index system+2*branch. ``coefficients``
    records ((a_0,w_0),(a_1,w_1)) as exact complex-coordinate pairs.
    The protocol branch is a logical input, not initialized workspace.
    """

    q: int
    core: tuple[int, ...]
    signal: int
    target: int
    system: int
    branch: int
    mode: int
    coefficients: tuple[_BranchCoefficients, _BranchCoefficients]
    tables: tuple[ResidualNativeFourRowTable, ResidualNativeFourRowTable]
    controlled_h_gates: Word
    good_reflection_gates: Word
    initial_reflection_gates: Word
    global_minus_gates: Word
    q_gates: Word
    q_inverse_gates: Word
    stages: tuple[NativeStateStage, ...]
    gates: Word
    certificate: BranchedResidualStateCertificate

    @property
    def branches(self) -> tuple[_BranchCoefficients, _BranchCoefficients]:
        return self.coefficients

    @property
    def protocol_branch(self) -> int:
        return self.branch

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
        """Error from the ideal amplified unitary, not a branch reset map."""
        return self.certificate.state_error_bound


def emit_branched_residual_state(
    branches: tuple[tuple[tuple[Fraction | int, Fraction | int], ...], ...],
    q: int,
) -> BranchedResidualNativeState:
    """Emit coherent residual preparation for an arbitrary protocol branch.

    ``branches`` is exactly ((a_0,w_0),(a_1,w_1)), each entry an exact
    dyadic (real,imag) pair. External promise, separately for each c:
    these pairs approximate a_c,w_c with Euclidean error at most
    delta=2**(-2*q-20), |a_c|**2+|w_c|**2/2=1, and |w_c|<=1.
    The desired isometry maps |c>|0> to
    |c>(a_c|0>+(w_c/sqrt(2))|1>), with both flags returned to zero.

    Checks cover exact types, dyadic denominators, q>=5, unit-disk
    consistency, and each branch's necessary normalization discrepancy
    bound 3*delta+3*delta**2/2. Passing does not prove the external
    promises. Coefficients are not normalized or rephased by this API.

    The coherent state error is three times the maximum of the two
    table bounds, strictly below 390*2**(-q). No branch is projected or
    initialized. The q+6 physical wires comprise q+2 arbitrary dirty
    core/signal wires, three initialized logical inputs, and the arbitrary
    protocol branch. For tail-table quadratic axis counts k_z,k_y, the
    exact T count is 3240*q+3793+210*(2*k_z+k_y), at most 3240*q+5053.
    """
    if isinstance(q, bool) or not isinstance(q, Integral) or q < 5:
        raise ValueError("q must be an integer at least five, not a boolean.")
    q = int(q)
    inputs = tuple(
        tuple(_pair(coefficient, "coefficient")
              for coefficient in _pair(branch, "branch coefficients"))
        for branch in _pair(branches, "branches"))
    if any(isinstance(value, bool) or not isinstance(value, (Fraction, Integral))
           for branch in inputs for pair in branch for value in pair):
        raise TypeError("Coefficient entries must be exact dyadic Fractions or integers, not booleans.")
    coefficients = tuple(
        tuple(tuple(_dyadic(value) for value in pair) for pair in branch)
        for branch in inputs)
    delta = Fraction(1, 1 << (2 * q + 20))
    tolerance = 3 * delta + Fraction(3, 2) * delta * delta
    discrepancies = tuple(abs(
        sum((value * value for value in root), Fraction(0))
        + sum((value * value for value in tail), Fraction(0)) / 2 - 1)
        for root, tail in coefficients)
    if any(discrepancy > tolerance for discrepancy in discrepancies):
        raise ValueError("The approximations are inconsistent with a promised branch normalization.")

    (a0, w0), (a1, w1) = coefficients
    zero = (Fraction(0), Fraction(0))
    m0 = emit_four_row_residual_table((a0, a0, a1, a1), q, enable_value=0)
    m1 = emit_four_row_residual_table((zero, w0, zero, w1), q, enable_value=1)
    system, branch = m0.addresses
    mode, target = m0.enable, m0.target
    controlled_h = _controlled_h(mode, system)
    mode_h = (_gate("H", mode),)
    q_gates = mode_h + controlled_h + m0.gates + m1.gates + mode_h
    q_inverse = _inverse(q_gates)
    good = _good_reflection(mode, target)
    # Identity on the protocol branch is essential for coherent selection.
    initial = _initial_reflection(system, mode, target)
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
    if t_count != 3240 * q + 3793 + 210 * (2 * k_z + k_y):
        raise ArithmeticError("The literal branched residual-state T count disagrees with its resource identity.")
    certificate = BranchedResidualStateCertificate(
        delta, tolerance, discrepancies, error, 3 * error, t_count, q + 6, q + 2)
    return BranchedResidualNativeState(
        q, m0.core, m0.signal, target, system, branch, mode, coefficients, (m0, m1),
        controlled_h, good, initial, minus, q_gates, q_inverse, stages, gates, certificate)
