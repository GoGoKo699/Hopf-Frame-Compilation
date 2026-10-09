"""Finite boundary-feedback matrices and exact component resource ledgers.

The dense matrices expose the identities in BOUNDARY_PROPAGATION.md. They
are diagnostic constructions, not a native feedback compiler. In particular,
the two-signal block encoding retains rejected signal amplitudes.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np


def _validate_height(height: int) -> None:
    if height < 1:
        raise ValueError("height must be positive.")


def _table(table: object) -> tuple[int, np.ndarray]:
    raw = np.asarray(table)
    if np.iscomplexobj(raw):
        raise ValueError("table must contain finite real values.")
    values = np.asarray(raw, dtype=float).reshape(-1)
    size = values.size + 1
    if size < 2 or size & (size - 1):
        raise ValueError("table must have length 2**height - 1 for height>=1.")
    if not np.all(np.isfinite(values)):
        raise ValueError("table must contain finite real values.")
    return size.bit_length() - 1, values


def boundary_incidence(height: int) -> tuple[np.ndarray, ...]:
    """Return centered cuts, boundary F, triangular K, and decrement P0.

    Columns have chronological breadth-first edge order. F has normalized
    columns, and K=X^dagger (I-P0) X/2. No target table is used here.
    """

    _validate_height(height)
    size = 1 << height
    cuts = np.zeros((size, size - 1))
    boundary = np.zeros_like(cuts)
    column = 0
    for depth in range(height):
        length = 1 << (height - depth - 1)
        for prefix in range(1 << depth):
            start = (2 * prefix + 1) * length
            cuts[:, column] = -length / size
            cuts[start:start + length, column] += 1
            boundary[start + length - 1, column] = 1 / np.sqrt(2)
            boundary[start - 1, column] = -1 / np.sqrt(2)
            column += 1
    decrement = np.roll(np.eye(size), -1, axis=0)
    baseline = (np.eye(size) - decrement) / 2
    triangular = cuts.T @ baseline @ cuts
    return cuts, boundary, triangular, decrement


def coefficient_matrix(table: object) -> np.ndarray:
    """Return H=(I+2i diag(t) K)^-1 diag(t), including zero entries."""

    height, values = _table(table)
    triangular = boundary_incidence(height)[2]
    diagonal = np.diag(values)
    return np.linalg.solve(np.eye(values.size) + 2j * diagonal @ triangular,
                           diagonal)


def propagation_matrices(table: object) -> tuple[np.ndarray, np.ndarray]:
    """Return the unpadded node maps B and G, both zero on terminal nodes."""

    height, values = _table(table)
    size = 1 << height
    propagation = np.zeros((size - 1, size - 1), dtype=complex)
    complement = np.zeros_like(propagation)
    for parent in range(1, size // 2):
        value = values[parent - 1]
        g = 1 / (1 + 1j * value)
        h = value * g
        propagation[2 * parent - 1:2 * parent + 1, parent - 1] = (1j * h, g)
        complement[2 * parent - 1:2 * parent + 1, parent - 1] = (
            -g.conjugate(), -1j * h.conjugate())
    return propagation, complement


def boundary_core(table: object) -> np.ndarray:
    """Return the complete logical matrix (I-2i F H F^dagger) P0."""

    height, values = _table(table)
    _, boundary, _, decrement = boundary_incidence(height)
    coefficient = coefficient_matrix(values)
    return (np.eye(1 << height)
            - 2j * boundary @ coefficient @ boundary.T) @ decrement


def shift_dilation_permutation(height: int) -> tuple[int, ...]:
    """Return D_S destinations in physical order (port, node), dummy included.

    D_S rotates the complete (height+1)-bit word left, then interchanges
    output modes (0,0) and (1,1). Its port-zero block is the partial
    doubling map on nonzero nonterminal nodes.
    """

    _validate_height(height)
    size = 1 << height
    destinations = []
    for source in range(2 * size):
        port, node = divmod(source, size)
        destination = (node // (size // 2)) * size + (2 * node + port) % size
        if destination == 0:
            destination = size + 1
        elif destination == size + 1:
            destination = 0
        destinations.append(destination)
    return tuple(destinations)


def propagation_dilation(table: object) -> np.ndarray:
    """Return the full 2N-port U_B=(I_port tensor V) D_S.

    The top-left N block is B padded by a zero dummy row and column.
    The node pair (0,1) of V is exactly identity.
    """

    height, values = _table(table)
    size = 1 << height
    bank = np.eye(size, dtype=complex)
    for parent in range(1, size // 2):
        value = values[parent - 1]
        g = 1 / (1 + 1j * value)
        h = value * g
        bank[2 * parent:2 * parent + 2, 2 * parent:2 * parent + 2] = (
            (1j * h, -g.conjugate()), (g, -1j * h.conjugate()))
    permutation = np.zeros((2 * size, 2 * size))
    permutation[shift_dilation_permutation(height), np.arange(2 * size)] = 1
    return np.kron(np.eye(2), bank) @ permutation


def propagation_block_encoding(table: object) -> np.ndarray:
    """Return H_c Z_c C_c(U_B) H_c in order (outer, inner, node).

    Its top-left N block is (I-B)/2, with B padded by dummy node zero.
    This complete unitary does not promise that either signal returns zero.
    """

    dilation = propagation_dilation(table)
    identity = np.eye(len(dilation))
    return np.block([[identity - dilation, identity + dilation],
                     [identity + dilation, identity - dilation]]) / 2


@dataclass(frozen=True)
class EndpointBlockLedger:
    """Exact reservations; error ratio multiplies the endpoint 2^(-2^n)."""

    height: int
    logical_dimension: int
    source_width: int
    core_qubits: int
    selector_qubits: int
    predicate_helper_qubits: int
    borrowed_signal_qubits: int
    clean_signal_qubits: int
    dirty_qubits: int
    total_width: int
    sectors_per_bank: int
    variable_banks: int
    predicate_literals: int
    free_table_rows: int
    error_ratio: Fraction


def endpoint_block_ledger(height: int, calls: int = 1,
                          extra_controls: int = 0) -> EndpointBlockLedger:
    """Price repeated blocks using self-borrowed queries and occupied port a.

    ``error_ratio`` is the accumulated error of all ``calls``, divided by
    eta. Additional controls are reserved arbitrary wires, not clean flags.
    """

    if height < 3:
        raise ValueError("the endpoint block reservation requires height>=3.")
    if calls < 1 or extra_controls < 0:
        raise ValueError("calls must be positive and extra_controls nonnegative.")
    overhead = (calls - 1).bit_length()
    if overhead + extra_controls > height - 3:
        raise ValueError("repeated calls and controls exceed the dirty allocation.")
    size = 1 << height
    source_width = size + overhead + 9
    dirty = source_width + 1 + extra_controls
    return EndpointBlockLedger(
        height=height, logical_dimension=size, source_width=source_width,
        core_qubits=source_width + 1, selector_qubits=0,
        predicate_helper_qubits=0, borrowed_signal_qubits=0,
        clean_signal_qubits=2, dirty_qubits=dirty,
        total_width=height + 2 + dirty, sectors_per_bank=1,
        variable_banks=2, predicate_literals=1 + extra_controls,
        free_table_rows=size // 2,
        error_ratio=Fraction(calls * 2 * 43, 1 << (source_width - size)),
    )


def full_tree_coin_bank(table: object) -> np.ndarray:
    """Return the SU(2) child bank on all 2N heap modes, fixing (0,1)."""

    height, values = _table(table)
    bank = np.eye(2 << height, dtype=complex)
    for node, value in enumerate(values, start=1):
        g = 1 / (1 + 1j * value)
        h = value * g
        bank[2 * node:2 * node + 2, 2 * node:2 * node + 2] = (
            (1j * h, -g.conjugate()), (g, -1j * h.conjugate()))
    return bank


def prefix_encoder(table: object) -> np.ndarray:
    """Return T_n=V_n E_even(T_(n-1)) on every heap-mode column.

    The dense recurrence is an ideal matrix diagnostic, not a fine-precision
    O(N)-T emitter. The encoder fixes dummy zero and preserves node depth.
    """

    height, values = _table(table)
    dimension = 2 << height
    embedded = np.eye(dimension, dtype=complex)
    if height > 1:
        even = np.arange(0, dimension, 2)
        embedded[np.ix_(even, even)] = prefix_encoder(values[:(1 << (height - 1)) - 1])
    return full_tree_coin_bank(values) @ embedded


def canonical_reversal_permutation(height: int) -> tuple[int, ...]:
    """Reverse each canonical even-child chain, including the zero dummy."""

    _validate_height(height)
    destinations = [0]
    for node in range(1, 2 << height):
        age = (node & -node).bit_length() - 1
        odd = node >> age
        last_age = height - (odd.bit_length() - 1)
        destinations.append(odd << (last_age - age))
    return tuple(destinations)


def feedback_entrance(table: object) -> tuple[np.ndarray, np.ndarray]:
    """Return the complete P_in permutation and n-bit entrance diagonal D."""

    height, values = _table(table)
    size = 1 << height
    node_of_input = np.zeros(size, dtype=int)
    diagonal = np.ones(size, dtype=complex)
    for depth in range(height):
        for prefix in range(1 << depth):
            node = (1 << depth) + prefix
            marker = (2 * prefix + 1) << (height - depth - 1)
            node_of_input[marker - 1] = node
            value = values[node - 1]
            diagonal[marker - 1] = -(1 - 1j * value) / (1 + 1j * value)
    incoming = np.zeros((2 * size, 2 * size))
    for port in range(2):
        incoming[2 * node_of_input + 1 - port, port * size + np.arange(size)] = 1
    return incoming, np.diag(diagonal)


def boundary_feedback(table: object) -> np.ndarray:
    """Return X_a T_n R_0 P_in (I_a tensor D P0), with one encoder.

    Its initialized-port columns equal (C_n(t), 0)^T. The other columns
    specify a complete unitary; the second signal is ideally untouched.
    """

    height, values = _table(table)
    size = 1 << height
    incoming, diagonal = feedback_entrance(values)
    reversal = np.zeros((2 * size, 2 * size))
    reversal[canonical_reversal_permutation(height), np.arange(2 * size)] = 1
    flip = np.kron(np.array([[0, 1], [1, 0]]), np.eye(size))
    decrement = np.roll(np.eye(size), -1, axis=0)
    return (flip @ prefix_encoder(values) @ reversal @ incoming
            @ np.kron(np.eye(2), diagonal @ decrement))


@dataclass(frozen=True)
class CoarseEncoderLedger:
    """Canonical SU(2) full-operator encoder; ratio multiplies 2**(-s)."""

    height: int
    coarse_bits: int
    source_widths: tuple[int, ...]
    dirty_qubits: int
    total_table_rows: int
    error_ratio: Fraction


def coarse_encoder_ledger(height: int) -> CoarseEncoderLedger:
    """Return the geometric q_j=s+n-j+7 schedule and its exact error sum."""

    if height < 3:
        raise ValueError("the coarse encoder reservation requires height>=3.")
    size = 1 << height
    bits = (size + height - 1) // height
    widths = tuple(bits + height - stage + 7 for stage in range(height))
    return CoarseEncoderLedger(
        height=height, coarse_bits=bits, source_widths=widths,
        dirty_qubits=bits + height + 9, total_table_rows=size - 1,
        error_ratio=sum((Fraction(86, 1 << (q - bits)) for q in widths), Fraction()),
    )
