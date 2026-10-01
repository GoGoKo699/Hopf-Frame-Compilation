#!/usr/bin/env python3
"""Reproduce the bounded native coarse-frame integration, without shot sampling."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compiler_robust_hopf.frames import real_tree_data
from compiler_robust_hopf.coarse_frame_decoder import decode_coarse_frame_histograms
from compiler_robust_hopf.native_coarse_fixture import (
    ANGLE_CASES, BRANCH, COARSE_DISTANCE_BOUND, HELPER, NQUBITS,
    coarse_word, controlled_observable_word, corrected_protocol_word,
    fine_frame_word, gate_counts, logical_blocks, original_protocol_word,
    preparation_word, simulate_columns,
)


def _error_norm(columns: np.ndarray) -> float:
    """Only the tiny input-space Gram matrix is diagonalized."""
    return float(np.sqrt(max(0., np.linalg.eigvalsh(columns.conj().T @ columns)[-1])))


def build_report(case: str = "balanced") -> dict:
    units = ANGLE_CASES[case]
    angles = np.asarray(units) * np.pi / 4
    data = real_tree_data(angles)  # The independent four-mode analytic oracle.
    derivatives = np.asarray(data.derivatives)
    blocks = logical_blocks(units)
    fine, coarse = blocks["W"], blocks["C"]
    hadamard = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
    walsh = np.kron(hadamard, hadamard)
    corrected = (walsh @ coarse.conj().T @ derivatives.T).T
    ideal_walsh = (walsh @ fine.conj().T @ derivatives.T).T
    dimension = 1 << NQUBITS
    # Four columns preserve arbitrary branch and helper inputs, including phase.
    inputs = np.zeros((dimension, 4), dtype=complex)
    target = inputs.copy()
    for dirty in (0, 1):
        for branch in (0, 1):
            offset = (dirty << HELPER) | (branch << BRANCH)
            column = 2 * dirty + branch
            inputs[offset, column] = 1
            target[offset:offset + 4, column] = (coarse if branch == 0 else fine)[:, 0]
    prepared = simulate_columns(preparation_word(units), inputs)
    preparation_error = _error_norm(prepared - target)
    execution_inputs = inputs[:, [0, 2]]
    indices = np.arange(dimension)
    signs = 1 - 2 * ((indices >> BRANCH) & 1)
    leaves = indices & 3
    identity = np.eye(2)
    observable_results = {}
    for name in ("tilted_x", "y"):
        pauli_y = np.array([[0., -1j], [1j, 0.]])
        axis = (np.array([[0., 1.], [1., 0.]]) + pauli_y) / np.sqrt(2)
        if name == "y":
            axis = pauli_y
        observable = fine @ np.kron(axis, identity) @ fine.conj().T
        analytic = 2 * np.real(derivatives @ observable @ data.state)
        score_operators = np.zeros((3, 2, 2), dtype=complex)
        wrong_operators = np.zeros_like(score_operators)
        x_operators = None
        for setting in ("X", "Y"):
            output = simulate_columns(corrected_protocol_word(units, name, setting), execution_inputs)
            quadrature = corrected.real if setting == "X" else corrected.imag
            wrong = ideal_walsh.real if setting == "X" else ideal_walsh.imag
            terms = []
            for j in range(3):
                # Factor 1/2 averages the independently chosen X/Y settings.
                score = 8 * signs * quadrature[j, leaves]
                terms.append(.5 * output.conj().T @ (score[:, None] * output))
                bad_score = 8 * signs * wrong[j, leaves]
                wrong_operators[j] += .5 * output.conj().T @ (bad_score[:, None] * output)
            score_operators += terms
            if setting == "X":
                x_operators = np.asarray(terms)
        original = simulate_columns(original_protocol_word(units, name), execution_inputs)
        original_operators = []
        for j in range(3):
            score = 4 * signs * ideal_walsh[j, leaves].real
            original_operators.append(original.conj().T @ (score[:, None] * original))
        wanted = analytic[:, None, None] * identity
        residual = max(float(np.linalg.norm(actual - expected, ord=2))
                       for actual, expected in zip(score_operators, wanted, strict=True))
        original_residual = max(float(np.linalg.norm(actual - expected, ord=2))
                                for actual, expected in zip(original_operators, wanted, strict=True))
        observable_results[name] = {
            "analytic_gradient": analytic.tolist(),
            "decoded_gradient": (np.trace(score_operators, axis1=1, axis2=2).real / 2).tolist(),
            "dirty_score_operator_error": residual,
            "original_score_operator_error": original_residual,
            "omitted_y_bias": float(np.max(np.abs(np.trace(x_operators - wanted, axis1=1, axis2=2).real / 2))),
            "ideal_walsh_bias": float(np.max(np.abs(np.trace(wrong_operators - wanted, axis1=1, axis2=2).real / 2))),
            "controlled_observable_gates": gate_counts(controlled_observable_word(units, name)),
            "corrected_x_gates": gate_counts(corrected_protocol_word(units, name, "X")),
            "corrected_y_gates": gate_counts(corrected_protocol_word(units, name, "Y")),
            "original_gates": gate_counts(original_protocol_word(units, name)),
        }
    # A fixed signed-record fixture audits histogram reconstruction separately
    # from the exact quantum-distribution calculation. These are not sampled data.
    hx, hy, shots = [2, -1, 0, 1], [0, 1, -1, 0], 8
    histogram_gradient = decode_coarse_frame_histograms(
        angles, blocks["coarse_blocks"], hx, hy, shots)
    direct_histogram = 8 * (corrected.real @ hx + corrected.imag @ hy) / shots
    histogram_error = float(np.max(np.abs(histogram_gradient - direct_histogram)))
    distance = float(np.linalg.norm(coarse - fine, ord=2))
    largest_residual = max(preparation_error, histogram_error, *(row[key] for row in observable_results.values()
                                              for key in ("dirty_score_operator_error", "original_score_operator_error")))
    if distance > COARSE_DISTANCE_BOUND + 1e-11 or largest_residual > 3e-10:
        raise RuntimeError("Native integration exceeds its stated finite-check tolerance.")
    return {
        "case": case,
        "scope": "Exact two-qubit finite-size fallback; no general residual-table emitter or advantage claim.",
        "layout": {"system": [0, 1], "branch": BRANCH, "unused_compiler_flags": [3, 4], "dirty_helper": HELPER},
        "largest_column_batch": [dimension, 4],
        "coarse_distance": distance,
        "analytic_distance_bound": float(COARSE_DISTANCE_BOUND),
        "preparation_isometry_error": preparation_error,
        "fixed_record_histogram_error": histogram_error,
        "floating_check_tolerance": 3e-10,
        "fine_frame_gates": gate_counts(fine_frame_word(units)),
        "coarse_frame_gates": gate_counts(coarse_word(units)),
        "preparation_gates": gate_counts(preparation_word(units)),
        "observables": observable_results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=tuple(ANGLE_CASES), default="balanced")
    parser.add_argument("--format", choices=("text", "json"), default="text")
    args = parser.parse_args()
    report = build_report(args.case)
    if args.format == "json":
        print(json.dumps(report, indent=2))
        return
    print(f"Native coarse-frame QBP: {args.case} (four logical modes)")
    print(f"Coarse distance: {report['coarse_distance']:.9g}; analytic bound {report['analytic_distance_bound']:.9g} < 1/64")
    print(f"Complete preparation error: {report['preparation_isometry_error']:.3g}")
    for name, row in report["observables"].items():
        print(f"{name}: gradient {np.asarray(row['decoded_gradient']).round(9).tolist()}")
        print(f"  dirty-input score error {row['dirty_score_operator_error']:.3g}; original error {row['original_score_operator_error']:.3g}")
        print(f"  T gates: corrected {row['corrected_x_gates']['T_count']}, original {row['original_gates']['T_count']}")
        print(f"  bias from omitted Y {row['omitted_y_bias']:.6g}; ideal Walsh scores {row['ideal_walsh_bias']:.6g}")
    print(report["scope"])


if __name__ == "__main__":
    main()
