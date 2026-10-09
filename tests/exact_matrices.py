"""Small exact matrix helpers using the retained Q(sqrt(2), i) arithmetic.

Native gate action is delegated to the existing sparse gate-word evaluator.
These helpers provide no synthesis algorithm or asymptotic resource claim.
"""
from pathlib import Path
import sys

# These retained fixtures intentionally support standalone execution. Keep
# their source/receipt hashes unchanged and restore the import path at once.
_fixture_path = str(Path(__file__).resolve().parents[1] / 'verification/fault_tolerant')
sys.path.insert(0, _fixture_path)
try:
    from exact_arithmetic import CQ2, Q2, CZ, CO, CI
    from gray_geometric import Gate, apply_word
finally:
    sys.path.pop(0)


def zeros(rows, columns=None):
    return [[CZ for _ in range(rows if columns is None else columns)]
            for _ in range(rows)]


def eye(size):
    return [[CO if i == j else CZ for j in range(size)] for i in range(size)]


def adjoint(matrix):
    return [[entry.conj() for entry in column] for column in zip(*matrix)]


def scale(matrix, scalar):
    return [[entry * scalar for entry in row] for row in matrix]


def add(left, right):
    return [[a + b for a, b in zip(row, other)]
            for row, other in zip(left, right)]


def multiply(left, right):
    result = zeros(len(left), len(right[0]))
    sparse_right = [[(j, value) for j, value in enumerate(row) if value != CZ]
                    for row in right]
    for i, row in enumerate(left):
        for k, value in enumerate(row):
            if value != CZ:
                for j, other in sparse_right[k]:
                    result[i][j] += value * other
    return result


def native_columns(width, word, columns=None):
    """Evaluate literal native words, retaining every flag/output amplitude."""
    native = []
    for name, *qubits in word:
        # Z is the exact native word SS; no quotient by scalar phase.
        if name == 'Z':
            native.extend([Gate('S', tuple(qubits))] * 2)
        else:
            native.append(Gate(name, tuple(qubits)))
    columns = range(1 << width) if columns is None else columns
    states = [apply_word(native, {column: CO}) for column in columns]
    return [[state.get(row, CZ) for state in states] for row in range(1 << width)]
