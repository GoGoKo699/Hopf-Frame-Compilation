#!/usr/bin/env python3
"""Reproduce both bounded native complex Hopf gradient streams, without sampling."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compiler_robust_hopf import native_complex_coarse_fixture as native
from compiler_robust_hopf.complex_analysis import complex_magnitude_gradient, complex_phase_gradient
from compiler_robust_hopf.complex_coarse_decoder import decode_complex_coarse_frame_histograms
from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import direct_real_frame, real_tree_data


def _error_norm(columns: np.ndarray) -> float:
    return float(np.sqrt(max(0., np.linalg.eigvalsh(columns.conj().T @ columns)[-1])))


def build_report(case: str = "balanced") -> dict:
    units = native.ANGLE_CASES[case]
    angles = np.asarray(units) * np.pi / 4
    data = real_tree_data(angles)
    real_frame = direct_real_frame(2, angles)
    blocks = native.logical_blocks(units)
    phases = blocks["leaf_phases"]
    gauge, fine, coarse = blocks["D0"], blocks["V"], blocks["C"]
    derivatives = (gauge @ np.asarray(data.derivatives).T).T
    h = np.array([[1., 1.], [1., -1.]]) / np.sqrt(2)
    walsh = np.kron(h, h)
    transformed = (walsh @ coarse.conj().T @ derivatives.T).T
    ideal = (walsh @ fine.conj().T @ derivatives.T).T
    ungauged = (walsh @ coarse.conj().T @ np.asarray(data.derivatives).T).T
    dimension = 1 << native.NQUBITS
    initial = np.zeros((dimension, 4), dtype=complex)
    target = np.zeros_like(initial)
    for dirty in (0, 1):
        for branch in (0, 1):
            offset = (dirty << native.HELPER) | (branch << native.BRANCH)
            column = 2 * dirty + branch
            initial[offset, column] = 1
            target[offset:offset + 4, column] = (coarse if branch == 0 else fine)[:, 0]
    preparation_error = _error_norm(native.simulate_columns(native.preparation_word(units), initial) - target)
    execution_inputs = initial[:, [0, 2]]
    labels = np.arange(dimension)
    signs, leaves = 1 - 2 * ((labels >> native.BRANCH) & 1), labels & 3

    def score_operator(output, weights):
        score = signs * weights[leaves]
        return output.conj().T @ (score[:, None] * output)

    def operator_error(actual, gradient):
        return max(float(np.linalg.norm(a - g * np.eye(2), ord=2))
                   for a, g in zip(actual, gradient, strict=True))

    results = {}
    for name in native.OBSERVABLES:
        if name == "tilted_x":
            axis = np.array([[0., np.exp(-1j * np.pi / 4)],
                             [np.exp(1j * np.pi / 4), 0.]])
            operator = real_frame @ np.kron(axis, np.eye(2)) @ real_frame.conj().T
        else:
            operator = real_frame @ np.kron(np.eye(2), h) @ real_frame.conj().T
        magnitude = complex_magnitude_gradient(angles, phases, operator)
        phase = complex_phase_gradient(angles, phases, operator)
        corrected = np.zeros((3, 2, 2), dtype=complex)
        ideal_scores, missing_gauge = np.zeros_like(corrected), np.zeros_like(corrected)
        x_piece = None
        for setting in ("X", "Y"):
            output = native.simulate_columns(native.magnitude_protocol_word(units, name, setting), execution_inputs)
            rows = transformed.real if setting == "X" else transformed.imag
            wrong = ideal.real if setting == "X" else ideal.imag
            absent = ungauged.real if setting == "X" else ungauged.imag
            piece = np.asarray([.5 * score_operator(output, 8 * row) for row in rows])
            corrected += piece
            ideal_scores += [.5 * score_operator(output, 8 * row) for row in wrong]
            missing_gauge += [.5 * score_operator(output, 8 * row) for row in absent]
            if setting == "X":
                x_piece = piece
        original_output = native.simulate_columns(native.original_magnitude_protocol_word(units, name), execution_inputs)
        original = []
        for node in range(1, 4):
            character = np.asarray([(-1) ** ((leaf & marker_label(node, 2)).bit_count()) for leaf in range(4)])
            original.append(score_operator(original_output, 2 * data.incoming_amplitude[node - 1] * character))
        phase_word = native.phase_protocol_word(units, name)
        original_phase_word = native.original_phase_protocol_word(units, name)
        # The exact finite-size phase preparation is the same supplied word
        # in both protocols, so propagate it once and price it in each ledger.
        if phase_word != original_phase_word:
            raise RuntimeError("This fixture's phase-stream comparison changed.")
        phase_output = native.simulate_columns(phase_word, execution_inputs)
        phase_scores = np.asarray([score_operator(phase_output, 2 * np.eye(4)[leaf]) for leaf in range(4)])
        corrected_x = native.gate_counts(native.magnitude_protocol_word(units, name, "X"))
        corrected_y = native.gate_counts(native.magnitude_protocol_word(units, name, "Y"))
        original_count = native.gate_counts(native.original_magnitude_protocol_word(units, name))
        phase_count = native.gate_counts(phase_word)
        results[name] = {
            "analytic_magnitude_gradient": magnitude.tolist(),
            "analytic_phase_gradient": phase.tolist(),
            "decoded_magnitude_gradient": (np.trace(corrected, axis1=1, axis2=2).real / 2).tolist(),
            "decoded_phase_gradient": (np.trace(phase_scores, axis1=1, axis2=2).real / 2).tolist(),
            "magnitude_dirty_score_error": operator_error(corrected, magnitude),
            "phase_dirty_score_error": operator_error(phase_scores, phase),
            "original_magnitude_dirty_score_error": operator_error(original, magnitude),
            "original_phase_dirty_score_error": operator_error(phase_scores, phase),
            "omitted_y_bias": operator_error(x_piece, magnitude),
            "ideal_walsh_bias": operator_error(ideal_scores, magnitude),
            "omitted_gauge_bias": operator_error(missing_gauge, magnitude),
            "controlled_observable_gates": native.gate_counts(native.controlled_observable_word(units, name)),
            "corrected_magnitude_x_gates": corrected_x,
            "corrected_magnitude_y_gates": corrected_y,
            "phase_gates": phase_count,
            "original_magnitude_gates": original_count,
            "original_phase_gates": native.gate_counts(original_phase_word),
            "corrected_pair_T_count": corrected_x["T_count"] + phase_count["T_count"],
            "original_pair_T_count": original_count["T_count"] + phase_count["T_count"],
        }
    hx, hy, shots = [2, -1, 0, 1], [0, 1, -1, 0], 8
    decoded = decode_complex_coarse_frame_histograms(
        angles, blocks["coarse_blocks"], blocks["phase_blocks"], phases, hx, hy, shots)
    direct = 8 * (transformed.real @ hx + transformed.imag @ hy) / shots
    histogram_error = float(np.max(abs(decoded - direct)))
    distance = float(np.linalg.norm(coarse - fine, ord=2))
    residual = max(preparation_error, histogram_error,
                   *(row[key] for row in results.values() for key in
                     ("magnitude_dirty_score_error", "phase_dirty_score_error",
                      "original_magnitude_dirty_score_error", "original_phase_dirty_score_error")))
    tolerance = 1e-9
    if distance > native.COARSE_DISTANCE_BOUND + tolerance or residual > tolerance:
        raise RuntimeError("Native complex integration exceeds its finite-check tolerance.")
    return {
        "case": case,
        "scope": "Complete exact two-qubit complex fallback; no general residual-table emitter or advantage claim.",
        "layout": {"system": [0, 1], "branch": native.BRANCH,
                   "unused_compiler_flags": [3, 4], "dirty_helper": native.HELPER},
        "largest_column_batch": [dimension, 4],
        "supplied_leaf_phases": phases.tolist(),
        "coarse_distance": distance,
        "analytic_distance_bound": native.COARSE_DISTANCE_BOUND,
        "preparation_isometry_error": preparation_error,
        "fixed_record_histogram_error": histogram_error,
        "floating_check_tolerance": tolerance,
        "fine_frame_gates": native.gate_counts(native.fine_frame_word(units)),
        "coarse_frame_gates": native.gate_counts(native.coarse_word(units)),
        "preparation_gates": native.gate_counts(native.preparation_word(units)),
        "observables": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=tuple(native.ANGLE_CASES), default="balanced")
    parser.add_argument("--format", choices=("text", "json"), default="text")
    args = parser.parse_args()
    report = build_report(args.case)
    if args.format == "json":
        print(json.dumps(report, indent=2))
        return
    print(f"Native complex Hopf QBP: {args.case} (four logical modes)")
    print(f"Coarse distance: {report['coarse_distance']:.9g}; analytic bound 3/200 < 1/64")
    print(f"Complete preparation error: {report['preparation_isometry_error']:.3g}")
    for name, row in report["observables"].items():
        print(f"{name}: magnitude {np.asarray(row['decoded_magnitude_gradient']).round(9).tolist()}")
        print(f"  phase {np.asarray(row['decoded_phase_gradient']).round(9).tolist()}")
        print(f"  dirty-input errors: magnitude {row['magnitude_dirty_score_error']:.3g}, phase {row['phase_dirty_score_error']:.3g}")
        print(f"  T gates per stream pair: corrected {row['corrected_pair_T_count']}, original {row['original_pair_T_count']}")
    print(report["scope"])


if __name__ == "__main__":
    main()
