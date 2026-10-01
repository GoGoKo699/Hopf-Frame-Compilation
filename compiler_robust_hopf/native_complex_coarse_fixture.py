"""A bounded native integration of the gauge-fixed complex Hopf protocol.

The two logical qubits, branch, two unused compiler flags, and active dirty
helper have the same six-wire layout as ``native_coarse_fixture``. Every
phase row acts across all suffixes. Its actual coarse word is Rz(gamma) K3,
so it need not be diagonal. All selections use complete per-reflection
dirty echoes; no helper input is assumed to be zero.

This exact finite-size fallback is a phase/order/workspace check, not a
general residual-table emitter or an asymptotic resource experiment.
"""
from __future__ import annotations

from collections.abc import Sequence
from numbers import Integral

import numpy as np

from . import native_coarse_fixture as real_native
from .native_coarse_fixture import (
    ANGLE_CASES,
    BRANCH,
    DEFAULT_ANGLES,
    ELEMENTARY_GATES,
    FLAG_S,
    FLAG_T,
    HELPER,
    HIGH,
    LOW,
    NQUBITS,
    Gate,
    Word,
    commutator_reflections,
    commutator_word,
    controlled_reflection_word,
    gate_counts,
    inverse_word,
    simulate_columns,
)


PHASE_UNITS = (9, -5, 7)
PHASE_MEAN = .37 + 4 * np.pi
OBSERVABLES = ("tilted_x", "hadamard_low")
# The exact K2 trace gives ||K2-I||<1/20, hence ||K3-I||<1/200.
# Two phase depths and one real-tree perturbation contribute at most 3/200.
COARSE_DISTANCE_BOUND = 3 / 200


def _gate(name: str, *qubits: int) -> Gate:
    return name, tuple(qubits)


def _units(value: int) -> int:
    if not isinstance(value, Integral):
        raise ValueError("Phase angles must be integer multiples of pi/4.")
    # Rz(pi)=-I: the literal period is eight units, not four.
    return (int(value) + 4) % 8 - 4


def _branch_value(value: int | None) -> int | None:
    if value is not None and value not in (0, 1):
        raise ValueError("branch_value must be zero, one, or None.")
    return value


def _prefix_query(prefix: tuple[int, int] | None) -> Word:
    """XOR a zero- or one-bit logical prefix predicate into the helper."""
    if prefix is None:
        return (_gate("X", HELPER),)
    qubit, value = prefix
    if qubit not in (LOW, HIGH) or value not in (0, 1):
        raise ValueError("The bounded phase table uses at most one logical prefix bit.")
    negative = (_gate("X", qubit),) if value == 0 else ()
    return negative + (_gate("CX", qubit, HELPER),) + negative


def _echoed_reflection(
    index: int | None,
    target: int,
    prefix: tuple[int, int] | None,
    branch_value: int | None,
) -> Word:
    """Select one literal involution while returning an arbitrary helper."""
    branch_value = _branch_value(branch_value)
    if prefix is not None and prefix[0] == target:
        raise ValueError("The target cannot be its own prefix address.")
    query = _prefix_query(prefix)
    controls = ((HELPER, 1),)
    if branch_value is not None:
        controls += ((BRANCH, branch_value),)
    reflection = controlled_reflection_word(index, target, controls)
    return query + reflection + inverse_word(query) + reflection


def _selected_k3(
    target: int, prefix: tuple[int, int] | None, branch_value: int | None,
) -> Word:
    # The inherited reflection indices are in matrix order, so reverse
    # them for a chronological gate word. Echo each involution separately.
    return sum((_echoed_reflection(index, target, prefix, branch_value)
                for index in reversed(commutator_reflections())), ())


def _selected_rz(
    units: int,
    target: int,
    prefix: tuple[int, int] | None,
    branch_value: int | None,
) -> Word:
    units = _units(units)
    if units == 0:
        return ()
    # B=H S-dagger maps Y to Z. Its conjugations act only on the target
    # and cancel literally on every inactive address/branch sector.
    conjugator = (_gate("SDG", target), _gate("H", target))
    step = (_echoed_reflection(None, target, prefix, branch_value)
            + _echoed_reflection(0, target, prefix, branch_value))
    if units < 0:
        step = inverse_word(step)
    return inverse_word(conjugator) + step * abs(units) + conjugator


def phase_table_word(*, coarse: bool = False, branch_value: int | None = None) -> Word:
    """Emit D0 or its actual coarse phase layers, optionally branch-selected.

    The root row acts on HIGH with no logical prefix, for both LOW values.
    The child rows act on LOW with HIGH=0 or HIGH=1. Thus all suffix pairs
    are included. For a coarse row K3 acts first, followed by the exact
    Rz rotation; no coarse row is assumed to remain diagonal.
    """
    branch_value = _branch_value(branch_value)
    if not isinstance(coarse, (bool, np.bool_)):
        raise ValueError("coarse must be a boolean.")
    word: Word = ()
    for units, target, prefix in (
        (PHASE_UNITS[0], HIGH, None),
        (PHASE_UNITS[1], LOW, (HIGH, 0)),
        (PHASE_UNITS[2], LOW, (HIGH, 1)),
    ):
        if coarse:
            word += _selected_k3(target, prefix, branch_value)
        word += _selected_rz(units, target, prefix, branch_value)
    return word


def fine_frame_word(angles: Sequence[int] = DEFAULT_ANGLES) -> Word:
    """Exact gauge-fixed V=D0 W_R, including every literal phase."""
    return real_native.fine_frame_word(angles) + phase_table_word()


def coarse_word(angles: Sequence[int] = DEFAULT_ANGLES) -> Word:
    """Actual C=C_P C_R, with C_R=W_R E and exact dirty return."""
    return real_native.coarse_word(angles) + phase_table_word(coarse=True)


def preparation_word(
    angles: Sequence[int] = DEFAULT_ANGLES, *, initialize_branch: bool = False,
) -> Word:
    """Coherently select C|00> on branch zero and V|00> on branch one.

    Selected E precedes common W_R. At each phase row, branch-selected K3
    precedes common prefix-selected Rz. This is a full-unitary selection
    identity, including arbitrary helper and branch inputs; the state
    preparation contract only initializes the logical system.
    """
    initial = (_gate("H", BRANCH),) if initialize_branch else ()
    word = initial + real_native.selected_error_word(0) + real_native.fine_frame_word(angles)
    for units, target, prefix in (
        (PHASE_UNITS[0], HIGH, None),
        (PHASE_UNITS[1], LOW, (HIGH, 0)),
        (PHASE_UNITS[2], LOW, (HIGH, 1)),
    ):
        word += _selected_k3(target, prefix, 0) + _selected_rz(units, target, prefix, None)
    return word


def controlled_observable_word(
    angles: Sequence[int] = DEFAULT_ANGLES, observable: str = "tilted_x",
) -> Word:
    """Use the fixed-tuple O=W_R P W_R-dagger, with no D0 conjugation.

    ``tilted_x`` and the retained optional ``y`` use the inherited high-bit
    observables. ``hadamard_low`` uses P=I tensor H, selected by one exact
    controlled reflection. The same literal oracle is used by both streams
    and both protocol comparisons.
    """
    if observable != "hadamard_low":
        return real_native.controlled_observable_word(angles, observable)
    real_word = real_native.fine_frame_word(angles)
    return (inverse_word(real_word)
            + controlled_reflection_word(0, LOW, ((BRANCH, 1),)) + real_word)


def magnitude_protocol_word(
    angles: Sequence[int] = DEFAULT_ANGLES,
    observable: str = "tilted_x",
    readout: str = "X",
) -> Word:
    """Common-coarse preparation, charged observable, actual C inverse, X/Y."""
    return (preparation_word(angles, initialize_branch=True)
            + controlled_observable_word(angles, observable)
            + inverse_word(coarse_word(angles))
            + (_gate("H", LOW), _gate("H", HIGH))
            + real_native._readout_word(readout))


def phase_protocol_word(
    angles: Sequence[int] = DEFAULT_ANGLES, observable: str = "tilted_x",
) -> Word:
    """Exact V preparation and the direct one-hot Y phase-gradient stream."""
    return (fine_frame_word(angles) + (_gate("H", BRANCH),)
            + controlled_observable_word(angles, observable)
            + real_native._readout_word("Y"))


def original_magnitude_protocol_word(
    angles: Sequence[int] = DEFAULT_ANGLES, observable: str = "tilted_x",
) -> Word:
    """Original exact-frame magnitude stream with the same gauge and oracle."""
    fine = fine_frame_word(angles)
    return (fine + (_gate("H", BRANCH),)
            + controlled_observable_word(angles, observable) + inverse_word(fine)
            + (_gate("H", LOW), _gate("H", HIGH), _gate("H", BRANCH)))


def original_phase_protocol_word(
    angles: Sequence[int] = DEFAULT_ANGLES, observable: str = "tilted_x",
) -> Word:
    """The two phase streams coincide for this exact finite-size target."""
    return phase_protocol_word(angles, observable)


def _local_phase_word(units: int, *, coarse: bool) -> Word:
    """Record the literal native row before exact reflection selection."""
    units = _units(units)
    step = (_gate("Z", 0), _gate("H", 0))
    if units < 0:
        step = inverse_word(step)
    conjugator = (_gate("SDG", 0), _gate("H", 0))
    rotation = inverse_word(conjugator) + step * abs(units) + conjugator
    return (commutator_word(target=0) if coarse else ()) + rotation


def _phase_matrix(rows: np.ndarray) -> np.ndarray:
    children = np.zeros((4, 4), complex)
    children[:2, :2], children[2:, 2:] = rows[1], rows[2]
    return children @ np.kron(rows[0], np.eye(2))


def logical_blocks(angles: Sequence[int] = DEFAULT_ANGLES) -> dict[str, object]:
    """Record actual small blocks and the supplied gauge; no dense oracle needed.

    ``coarse_blocks``/``real_coarse_blocks`` are the three C_R tree rows.
    ``phase_blocks`` are the three actual all-suffix C_P rows. W denotes
    the real frame, V the exact gauge-fixed target, and C the actual coarse
    product. Tests should build their own independent logical reference.
    """
    real = real_native.logical_blocks(angles)
    phase_rows = np.asarray([simulate_columns(_local_phase_word(value, coarse=True), nqubits=1)
                             for value in PHASE_UNITS])
    fine_rows = np.asarray([simulate_columns(_local_phase_word(value, coarse=False), nqubits=1)
                            for value in PHASE_UNITS])
    phase_coarse, gauge = _phase_matrix(phase_rows), _phase_matrix(fine_rows)
    root, left, right = np.asarray(PHASE_UNITS) * np.pi / 4
    phases = PHASE_MEAN + np.array([-root - left, -root + left,
                                    root - right, root + right])
    real_blocks = np.asarray(real["coarse_blocks"])
    return {
        "W": real["W"], "W_R": real["W"], "E": real["E"], "K": real["K"],
        "C_R": real["C"], "C_P": phase_coarse, "D0": gauge,
        "V": gauge @ real["W"], "C": phase_coarse @ real["C"],
        "coarse_blocks": real_blocks, "real_coarse_blocks": real_blocks,
        "phase_blocks": phase_rows, "fine_phase_blocks": fine_rows,
        "leaf_phases": phases, "phase_mean": PHASE_MEAN,
        "phase_angles": np.asarray(PHASE_UNITS) * np.pi / 4,
    }
