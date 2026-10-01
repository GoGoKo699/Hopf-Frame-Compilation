"""A bounded, fully native two-qubit coarse-reference integration fixture.

This is an explicit finite-size example, not a general state compiler or a
resource-scaling experiment. Gate words are chronological immutable tuples
``(name, (qubits, ...))``. Qubit zero is the low integer bit. The six-wire
layout is system low/high, protocol branch, two unused compiler flags, and
one genuinely borrowed helper. Every emitted gate is in the declared
elementary Clifford+T alphabet; all inverses retain their literal phases.

The complex coarse frame is C = W E, where E acts on the root pair (0, 2).
Its SU(2) block is K3, with K0 = [T, HTH] and
K_(r+1) = [K_r, S K_r S-dagger]. The commutator estimate certifies
||K3-I|| <= 128 (1-1/sqrt(2))**8 < 1/64.
"""
from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from numbers import Integral

import numpy as np


Gate = tuple[str, tuple[int, ...]]
Word = tuple[Gate, ...]
LOW, HIGH, BRANCH, FLAG_S, FLAG_T, HELPER = range(6)
NQUBITS = 6
ELEMENTARY_GATES = frozenset({"H", "X", "Z", "S", "SDG", "T", "TDG", "CX"})
DEFAULT_ANGLES = (1, 1, 1)
ANGLE_CASES = {"balanced": (1, 1, 1), "singular": (0, 1, -1), "signed": (1, -1, 0)}
COARSE_DISTANCE_BOUND = 128 * (1 - 1 / np.sqrt(2)) ** 8


def _gate(name: str, *qubits: int) -> Gate:
    return name, tuple(qubits)


def inverse_word(word: Iterable[Gate]) -> Word:
    """The actual circuit inverse, without discarding scalar phases."""
    inverse = {"T": "TDG", "TDG": "T", "S": "SDG", "SDG": "S"}
    return tuple((inverse.get(name, name), qubits) for name, qubits in reversed(tuple(word)))


def toffoli_word(left: int, right: int, target: int) -> Word:
    """Literal 15-gate, seven-T Toffoli, with no relative-phase replacement."""
    if len({left, right, target}) != 3:
        raise ValueError("Toffoli wires must be distinct.")
    return (
        _gate("H", target), _gate("CX", right, target), _gate("TDG", target),
        _gate("CX", left, target), _gate("T", target), _gate("CX", right, target),
        _gate("TDG", target), _gate("CX", left, target), _gate("T", right),
        _gate("T", target), _gate("H", target), _gate("CX", left, right),
        _gate("T", left), _gate("TDG", right), _gate("CX", left, right),
    )


def _t_power(power: int, target: int) -> Word:
    names = ((), ("T",), ("S",), ("S", "T"), ("Z",),
             ("Z", "T"), ("SDG",), ("TDG",))[power % 8]
    return tuple(_gate(name, target) for name in names)


def reflection_word(index: int, target: int = HIGH) -> Word:
    """Literal R_j = T**j H T**(-j), in chronological order."""
    return _t_power(-index, target) + (_gate("H", target),) + _t_power(index, target)


def controlled_reflection_word(
    index: int | None,
    target: int,
    controls: Sequence[tuple[int, int]],
) -> Word:
    """Select R_index, or Z for index=None, using at most two literals.

    B = S H T H S-dagger obeys B Z B-dagger = H. Thus each selected
    reflection is a selected Z conjugated by T**j B. The conjugating T
    gates are unconditional; this does not assume native controlled T.
    """
    controls = tuple(controls)
    if len(controls) > 2 or any(value not in (0, 1) for _, value in controls):
        raise ValueError("This bounded fixture supports at most two binary controls.")
    wires = (target,) + tuple(qubit for qubit, _ in controls)
    if len(set(wires)) != len(wires):
        raise ValueError("Reflection target and controls must be distinct.")
    negative = tuple(_gate("X", qubit) for qubit, value in controls if value == 0)
    positive = tuple(qubit for qubit, _ in controls)
    if not positive:
        selected_z = (_gate("Z", target),)
    elif len(positive) == 1:
        selected_z = (_gate("H", target), _gate("CX", positive[0], target), _gate("H", target))
    else:
        selected_z = ((_gate("H", target),)
                      + toffoli_word(*positive, target)
                      + (_gate("H", target),))
    if index is None:
        middle = selected_z
    else:
        conjugator = ((_gate("SDG", target), _gate("H", target), _gate("T", target),
                       _gate("H", target), _gate("S", target))
                      + _t_power(index, target))
        middle = inverse_word(conjugator) + selected_z + conjugator
    return negative + middle + inverse_word(negative)


def commutator_word(level: int = 3, target: int = HIGH) -> Word:
    """Direct literal K_level word; K3 contains 256 T/TDG appearances."""
    if not isinstance(level, Integral) or not 0 <= level <= 3:
        raise ValueError("Only the four fixed commutators K0 through K3 are supported.")
    word = tuple(_gate(name, target) for name in ("H", "TDG", "H", "TDG", "H", "T", "H", "T"))
    for _ in range(level):
        conjugate = (_gate("SDG", target),) + word + (_gate("S", target),)
        # Matrix order K J K-dagger J-dagger, emitted rightmost first.
        word = inverse_word(conjugate) + inverse_word(word) + conjugate + word
    return word


def commutator_reflections(level: int = 3) -> tuple[int, ...]:
    """Reflection indices of K_level in matrix order, not chronological order."""
    if not isinstance(level, Integral) or not 0 <= level <= 3:
        raise ValueError("Only K0 through K3 are supported.")
    word = (1, 2, 1, 0)  # K0 = R1 R2 R1 R0, including literal phase.
    for _ in range(level):
        conjugate = tuple((index + 2) % 8 for index in word)
        word = word + conjugate + word[::-1] + conjugate[::-1]
    return word


def root_error_word() -> Word:
    """E: K3 on the low=0 root pair, identity on low=1, no work inputs used."""
    return sum((controlled_reflection_word(index, HIGH, ((LOW, 0),))
                for index in reversed(commutator_reflections())), ())


def selected_error_word(branch_value: int = 0) -> Word:
    """Branch-selected E with a complete dirty-helper echo per reflection.

    The low=0 predicate is XOR-loaded into an arbitrary helper. Two uses
    of the same involution cancel its original helper value. Completing
    this echo separately for each reflection is essential: K3 itself is
    not an involution and cannot replace one echoed reflection.
    """
    if branch_value not in (0, 1):
        raise ValueError("branch_value must be zero or one.")
    query = (_gate("X", LOW), _gate("CX", LOW, HELPER), _gate("X", LOW))
    word: Word = ()
    for index in reversed(commutator_reflections()):
        reflection = controlled_reflection_word(index, HIGH, ((BRANCH, branch_value), (HELPER, 1)))
        word += query + reflection + inverse_word(query) + reflection
    return word


def _angle_units(angles: Sequence[int]) -> tuple[int, int, int]:
    if len(angles) != 3 or any(not isinstance(value, Integral) for value in angles):
        raise ValueError("Give exactly three integer multiples of pi/4.")
    # Ry has literal period 2*pi, or eight units, not four: Ry(pi)=-I.
    # The canonical interval [-4, 3] keeps every fixture word bounded.
    return tuple((int(value) + 4) % 8 - 4 for value in angles)


def _selected_ry(units: int, target: int, control: tuple[int, int]) -> Word:
    # Ry(pi/4)=H Z, so the chronological order is selected Z then H.
    step = (controlled_reflection_word(None, target, (control,))
            + controlled_reflection_word(0, target, (control,)))
    if units < 0:
        step = inverse_word(step)
    return step * abs(units)


def fine_frame_word(angles: Sequence[int] = DEFAULT_ANGLES) -> Word:
    """Exact n=2 Hopf W; integer pi/4 units reduce modulo eight to [-4, 3]."""
    root, left, right = _angle_units(angles)
    return (_selected_ry(root, HIGH, (LOW, 0))
            + _selected_ry(left, LOW, (HIGH, 0))
            + _selected_ry(right, LOW, (HIGH, 1)))


def coarse_word(angles: Sequence[int] = DEFAULT_ANGLES) -> Word:
    """Actual native C=W E; its inverse is always inverse_word(C)."""
    return root_error_word() + fine_frame_word(angles)


def preparation_word(
    angles: Sequence[int] = DEFAULT_ANGLES, *, initialize_branch: bool = False,
) -> Word:
    """Prepare C|0> or W|0> coherently, with an arbitrary branch input.

    By default the branch is unchanged, allowing its entire two-dimensional
    input space to be audited. Setting initialize_branch=True first applies
    H to a zero branch. The two compiler flags are unused finite-size work.
    """
    initial = (_gate("H", BRANCH),) if initialize_branch else ()
    return initial + selected_error_word(0) + fine_frame_word(angles)


def controlled_observable_word(
    angles: Sequence[int] = DEFAULT_ANGLES, observable: str = "tilted_x",
) -> Word:
    """Controlled O=W P W-dagger, with P=TXT-dagger or Y on high.

    The W and single-qubit phase conjugators act unconditionally. Their
    actual inverses give literal identity on branch zero. Only the central
    X is controlled, by an ordinary CNOT.
    """
    if observable == "tilted_x":
        middle = (_gate("TDG", HIGH), _gate("CX", BRANCH, HIGH), _gate("T", HIGH))
    elif observable == "y":
        middle = (_gate("SDG", HIGH), _gate("CX", BRANCH, HIGH), _gate("S", HIGH))
    else:
        raise ValueError("observable must be 'tilted_x' or 'y'.")
    fine = fine_frame_word(angles)
    return inverse_word(fine) + middle + fine


def _readout_word(readout: str) -> Word:
    if readout == "X":
        return (_gate("H", BRANCH),)
    if readout == "Y":
        return (_gate("SDG", BRANCH), _gate("H", BRANCH))
    raise ValueError("readout must be 'X' or 'Y'.")


def corrected_protocol_word(
    angles: Sequence[int] = DEFAULT_ANGLES, observable: str = "tilted_x", readout: str = "X",
) -> Word:
    """Complete native preparation, controlled O, actual C inverse, and readout."""
    return (preparation_word(angles, initialize_branch=True)
            + controlled_observable_word(angles, observable)
            + inverse_word(coarse_word(angles))
            + (_gate("H", LOW), _gate("H", HIGH)) + _readout_word(readout))


def original_protocol_word(
    angles: Sequence[int] = DEFAULT_ANGLES, observable: str = "tilted_x",
) -> Word:
    """Exact-frame comparison with the same observable and elementary pricing."""
    fine = fine_frame_word(angles)
    return (fine + (_gate("H", BRANCH),)
            + controlled_observable_word(angles, observable) + inverse_word(fine)
            + (_gate("H", LOW), _gate("H", HIGH), _gate("H", BRANCH)))


def gate_counts(word: Iterable[Gate]) -> dict[str, int]:
    """Literal unsimplified full-word counts; no proxy or cancellation estimate."""
    by_name = Counter(name for name, _ in word)
    if set(by_name) - ELEMENTARY_GATES:
        raise ValueError("The word contains a non-elementary gate.")
    counts = {name: by_name[name] for name in sorted(ELEMENTARY_GATES)}
    counts.update(total=sum(by_name.values()), T_count=by_name["T"] + by_name["TDG"])
    counts["Clifford_count"] = counts["total"] - counts["T_count"]
    return counts


def simulate_columns(
    word: Iterable[Gate], inputs: np.ndarray | None = None, *, nqubits: int = NQUBITS,
) -> np.ndarray:
    """Propagate bounded native columns, including arbitrary helper inputs.

    The default propagates all columns on six wires. Supplying columns
    avoids building a square matrix. This diagnostic is deliberately
    limited to at most six wires and does not claim scalable simulation.
    """
    if not isinstance(nqubits, Integral) or not 1 <= nqubits <= NQUBITS:
        raise ValueError("The bounded fixture supports one through six wires.")
    size = 1 << nqubits
    result = np.eye(size, dtype=complex) if inputs is None else np.array(inputs, dtype=complex, copy=True)
    if result.ndim not in (1, 2) or result.shape[0] != size:
        raise ValueError("Input vectors or columns must have 2**nqubits rows.")
    indices = np.arange(size)
    low = {qubit: indices[(indices & (1 << qubit)) == 0] for qubit in range(nqubits)}
    high = {qubit: low[qubit] | (1 << qubit) for qubit in range(nqubits)}
    phases = {"Z": -1, "S": 1j, "SDG": -1j,
              "T": np.exp(1j * np.pi / 4), "TDG": np.exp(-1j * np.pi / 4)}
    for name, qubits in word:
        expected_arity = 2 if name == "CX" else 1
        if (name not in ELEMENTARY_GATES or len(qubits) != expected_arity
                or len(set(qubits)) != len(qubits)
                or any(not isinstance(q, Integral) or not 0 <= q < nqubits for q in qubits)):
            raise ValueError("Invalid elementary gate or wire for this register.")
        if name == "X":
            result = result[indices ^ (1 << qubits[0])]
        elif name == "CX":
            control, target = qubits
            result = result[indices ^ (((indices >> control) & 1) << target)]
        elif name == "H":
            q = qubits[0]
            first, second = result[low[q]].copy(), result[high[q]].copy()
            result[low[q]] = (first + second) / np.sqrt(2)
            result[high[q]] = (first - second) / np.sqrt(2)
        else:
            result[high[qubits[0]]] *= phases[name]
    return result


def logical_blocks(
    angles: Sequence[int] = DEFAULT_ANGLES,
) -> dict[str, np.ndarray | tuple[np.ndarray, ...]]:
    """Actual emitted small blocks; independent tests should build their own oracle."""
    fine = simulate_columns(fine_frame_word(angles), nqubits=2)
    error = simulate_columns(root_error_word(), nqubits=2)
    root = error[np.ix_((0, 2), (0, 2))]
    units = _angle_units(angles)
    rotations = []
    for value in units:
        local_word = (_gate("Z", 0), _gate("H", 0))
        if value < 0:
            local_word = inverse_word(local_word)
        rotations.append(simulate_columns(local_word * abs(value), nqubits=1))
    rotation = rotations[0]
    coarse = simulate_columns(coarse_word(angles), nqubits=2)
    return {"W": fine, "E": error, "C": coarse, "K": root,
            "fine_root": rotation, "coarse_root": rotation @ root,
            "coarse_blocks": (rotation @ root, rotations[1], rotations[2])}
