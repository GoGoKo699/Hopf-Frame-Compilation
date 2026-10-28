"""Independent exact augmented-decoder audit in the helper computational basis.

The word here contains the specified algebraic XOR and Fourier primitives.
Their native synthesis is charged by the retained multiplication-XOR proof;
this file does not call the helper-sector evaluator or claim a Clifford-only
physical implementation.
"""
from pathlib import Path
import hashlib
import json
import os
import sys
import time


def repository_root():
    for parent in Path(__file__).resolve().parents:
        if (parent / 'WORKSPACE.md').is_file() and (parent / 'tests').is_dir():
            return parent
    configured = os.environ.get('HOPF_FRAME_REPO')
    if configured:
        return Path(configured).resolve()
    return Path('/workspace/scratch/e3eb9dbc5d0b/Hopf-Frame-Compilation-continuation')


sys.path.insert(0, str(repository_root()))
from tests.test_recursive_field_shear_decoder import (
    _tower_product as mul, _tower_trace as trace, _address_field, _tower_shear)
from tests.test_field_shear_intertwiner import _translation


def word(level, top, offset=0):
    if level == 1:
        return [('CX', offset + 1, offset)]
    ell, degree = 1 << (level - 2), 1 << (level - 1)
    helper = (1 << top) + ell - 1
    recoupler = [('F', level, offset + degree, 1),
                 ('ADD', level, offset, helper, 'before'),
                 ('F', level, offset + degree, 1 << ell),
                 ('ADD', level, offset, helper, 'after')]
    recoupler += [gate for j in range(degree)
                  for gate in [('CX', offset + j, offset + degree + j),
                               ('CX', offset + degree + j, offset + j)]]
    recoupler += [('F', level - 1, offset + ell, 1),
                  ('F', level - 1, offset + degree + ell, 1)]
    return (recoupler + word(level - 1, top, offset)
            + word(level - 1, top, offset + degree))


def apply(operations, column):
    state = {column: 1}
    normalization = 0
    for kind, *args in operations:
        if kind == 'CX':
            control, target = args
            state = {row ^ (((row >> control) & 1) << target): value
                     for row, value in state.items()}
        elif kind == 'ADD':
            level, offset, helper, which = args
            ell, degree = 1 << (level - 2), 1 << (level - 1)
            mask = (1 << degree) - 1
            low = (1 << ell) - 1
            changed = {}
            for row, value in state.items():
                u, v = (row >> offset) & mask, (row >> (offset + degree)) & mask
                u0, u1, v0, v1 = u & low, u >> ell, v & low, v >> ell
                if which == 'before':
                    # Direct relative trace, independent of the two-product
                    # expression used for the proposed native resource count.
                    delta = (mul(level, u, v) >> ell) ^ mul(level - 1, u0, u1)
                    compressed = (mul(level - 1, u0, v1 ^ u1)
                                  ^ mul(level - 1, u1, v0 ^ v1))
                    assert delta == compressed
                else:
                    assert which == 'after'
                    delta = mul(level - 1, u0, v1) ^ mul(level - 1, u1, v0)
                target = row ^ (delta << helper)
                assert target not in changed  # The whole add is a permutation.
                changed[target] = value
            state = changed
        else:
            assert kind == 'F'
            level, offset, multiplier = args
            degree = 1 << (level - 1)
            mask = ((1 << degree) - 1) << offset
            normalization += degree
            changed = {}
            for row, value in state.items():
                source = (row & mask) >> offset
                for target in range(1 << degree):
                    exponent = trace(level, mul(level, multiplier,
                                                mul(level, source, target)))
                    out = (row & ~mask) | (target << offset)
                    changed[out] = changed.get(out, 0) + (-value if exponent else value)
            state = {row: value for row, value in changed.items() if value}
    return state, normalization


def helper_phase(top, address, basis):
    exponent = 0
    for level in range(2, top + 1):
        ell = 1 << (level - 2)
        helper = (1 << top) + ell - 1
        value = (basis >> helper) & ((1 << ell) - 1)
        coefficient = _address_field(level - 1, address & ((1 << (level - 1)) - 1))
        exponent ^= trace(level - 1, mul(level - 1, coefficient, value))
    return -1 if exponent else 1


def run():
    started = time.monotonic()
    cases = []
    for top in (2, 3):
        core = 1 << top
        helpers = (1 << (top - 1)) - 1
        width, size = core + helpers, 1 << (core + helpers)
        operations = word(top, top)
        columns = []
        normalizations = set()
        for column in range(size):
            state, normalization = apply(operations, column)
            assert sum(value * value for value in state.values()) == 1 << normalization
            columns.append(state)
            normalizations.add(normalization)
        assert len(normalizations) == 1
        normalization, = normalizations
        # The primitive inverses are the same Fourier/add/CNOT operators in
        # reversed order. Check every inverse column against the transpose.
        expected_inverse = [{} for _ in range(size)]
        for column, entries in enumerate(columns):
            for row, value in entries.items():
                expected_inverse[row][column] = value
        for column in range(size):
            inverse, inverse_normalization = apply(list(reversed(operations)), column)
            assert inverse_normalization == normalization
            assert inverse == expected_inverse[column], (top, column)
        comparisons = 0
        core_mask = (1 << core) - 1
        for address in range(1 << top):
            for column in range(size):
                source_column = ((column & ~core_mask)
                                 | _tower_shear(top, address, column & core_mask))
                source_sign = helper_phase(top, address, column)
                left = {row: source_sign * value
                        for row, value in columns[source_column].items()}
                right = {}
                for row, value in columns[column].items():
                    target_row = ((row & ~core_mask)
                                  | _translation(top, address, row & core_mask))
                    right[target_row] = helper_phase(top, address, row) * value
                assert left == right, (top, address, column)
                comparisons += 1
        cases.append({'p': top, 'core_qubits': core, 'moving_helper_qubits': helpers,
                      'complete_computational_basis_columns': size,
                      'address_column_comparisons': comparisons,
                      'actual_reversed_algebraic_word_columns': size,
                      'fixed_common_normalization_exponent': normalization,
                      'all_zero_U_columns_included': True})
        print('PASS p=', top, 'columns=', size, flush=True)
    report = {'status': 'PASS', 'method': 'Exact integer amplitudes; computational helper basis',
              'scope': 'Full algebraic primitive word, not a separately emitted native multiplication circuit',
              'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'cases': cases, 'seconds': time.monotonic() - started}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
