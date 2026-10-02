#!/usr/bin/env python3
"""Report exact fixture, gate, and bias certificates for native residual QBP."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from math import sqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compiler_robust_hopf.native_residual_qbp import (
    decode_magnitude_histograms, decode_phase_histogram, emit_residual_qbp,
)


def build_report(q: int = 16) -> dict:
    """Emit finite words and exact certificates, without a statevector."""
    word = emit_residual_qbp(q)
    data = word.data
    cu, su, cv, sv = data.cos_u, data.sin_u, data.cos_v, data.sin_v
    magnitude_factor = -(cu * cu - su * su) + 4 * cu * su * cv * sv
    phase_factor = (cu * cu - su * su) * (cv * cv - sv * sv)
    streams = {}
    for stream in (word.magnitude_x, word.magnitude_y, word.phase):
        counts = dict(sorted(Counter(name for name, _ in stream.gates).items()))
        streams[stream.name] = {
            "gate_counts": counts,
            "total_gates": len(stream.gates),
            "T_count": stream.t_count,
            "preparation_error_bound": str(stream.preparation_error_bound),
            "stages": [stage.name for stage in stream.stages],
        }
    return {
        "scope": "Fixed one-qubit fine-residual QBP integration; no sampling or native fine-frame comparison.",
        "q": q,
        "physical_qubits": word.nqubits,
        "dirty_qubits": word.dirty_qubits,
        "clean_compiler_flags": [word.target, word.mode],
        "initialized_system": word.system,
        "separate_initialized_protocol_branch": word.branch,
        "half_angle_denominators": [data.angle_denominator, data.phase_denominator],
        "exact_coarse_distance_squared": str(data.coarse_distance_squared),
        "coarse_distance_bound": str(data.coarse_distance_bound),
        "coefficient_error_squared_bounds": [str(value) for value in data.coefficient_error_squared_bounds],
        "coefficient_error_promise_squared": str(data.input_error_bound ** 2),
        "ideal_magnitude_gradient": {
            "sqrt_two_factor": str(magnitude_factor),
            "decimal": sqrt(2) * float(magnitude_factor),
        },
        "ideal_phase_gradient": {
            "inverse_sqrt_two_factors": [str(phase_factor), str(-phase_factor)],
            "decimals": [float(phase_factor) / sqrt(2), -float(phase_factor) / sqrt(2)],
        },
        "streams": streams,
        "gradient_bias_bounds": {
            "magnitude_fair_XY": str(10 * word.magnitude_x.preparation_error_bound),
            "phase_vector": str(4 * word.phase.preparation_error_bound),
            "magnitude_fair_XY_decimal": float(10 * word.magnitude_x.preparation_error_bound),
            "phase_vector_decimal": float(4 * word.phase.preparation_error_bound),
        },
        "T_count_per_magnitude_phase_pair": word.magnitude_x.t_count + word.phase.t_count,
        "example_empirical_decoding": {
            "magnitude_X_counts": [2, -1], "magnitude_Y_counts": [1, 0],
            "magnitude_shots": 6,
            "magnitude": str(decode_magnitude_histograms((2, -1), (1, 0), 6)),
            "phase_counts": [2, -1], "phase_shots": 5,
            "phase": [str(value) for value in decode_phase_histogram((2, -1), 5)],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", type=int, default=16, help="Source precision, at least five (default: 16).")
    parser.add_argument("--format", choices=("text", "json"), default="text")
    args = parser.parse_args()
    if args.q < 5:
        parser.error("--q must be at least five")
    report = build_report(args.q)
    if args.format == "json":
        print(json.dumps(report, indent=2))
        return
    print(f"Native residual QBP: q={args.q}, {report['physical_qubits']} physical wires, "
          f"{report['dirty_qubits']} arbitrary dirty wires")
    print(f"Exact squared coarse distance: {report['exact_coarse_distance_squared']} < 1/4096")
    print(f"Ideal magnitude gradient: {report['ideal_magnitude_gradient']['decimal']:.12g}")
    print(f"Ideal raw phase gradients: {report['ideal_phase_gradient']['decimals']}")
    for name, stream in report["streams"].items():
        print(f"{name}: T={stream['T_count']}, gates={stream['total_gates']}")
    bounds = report["gradient_bias_bounds"]
    print(f"Certified gradient-bias bounds: fair X/Y magnitude={bounds['magnitude_fair_XY_decimal']:.9g}, "
          f"phase vector={bounds['phase_vector_decimal']:.9g}")
    print(f"T gates per magnitude/phase execution pair: {report['T_count_per_magnitude_phase_pair']}")
    print(report["scope"])


if __name__ == "__main__":
    main()
