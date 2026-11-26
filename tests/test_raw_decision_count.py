"""Finite Pauli bookkeeping checks for the analytic raw-decision count bound."""

import unittest

import numpy as np

from tests.test_operator_source_compiler import _word_matrix


def _pauli_matrix(width, record):
    """Literal i**phase X(x) Z(z), with qubit zero the low integer bit."""
    x, z, phase = record
    matrix = np.zeros((1 << width, 1 << width), dtype=complex)
    for basis in range(1 << width):
        matrix[basis ^ x, basis] = (
            1j ** phase * (-1) ** ((z & basis).bit_count() % 2)
        )
    return matrix


def _pauli_record(matrix):
    """Recover a signed Pauli from an independently simulated Clifford image."""
    size = len(matrix)
    width = (size - 1).bit_length()
    x = int(np.argmax(abs(matrix[:, 0])))
    phase = min(range(4), key=lambda k: abs(matrix[x, 0] - 1j ** k))
    z = 0
    for bit in range(width):
        amplitude = matrix[x ^ (1 << bit), 1 << bit] / (1j ** phase)
        if amplitude.real < 0:
            z |= 1 << bit
    record = (x, z, phase)
    np.testing.assert_allclose(matrix, _pauli_matrix(width, record), atol=2e-12)
    return record


def _kernel_basis(labels):
    """Relations among arbitrary-width X labels, as vectors on generator indices."""
    pivots = {}
    relations = []
    for index, label in enumerate(labels):
        value, relation = label, 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = (value, relation)
                break
            other, coefficients = pivots[pivot]
            value ^= other
            relation ^= coefficients
        else:
            relations.append(relation)
    return tuple(relations)


def _descriptor(records, address_bits):
    """B, phases, ker(X), and logical Z labels; physical width is discarded."""
    gram = tuple(tuple((left[1] & right[0]).bit_count() % 2
                       for right in records) for left in records)
    phases = tuple(record[2] for record in records)
    kernel = _kernel_basis([record[0] for record in records])
    logical_z = tuple(record[1] & ((1 << address_bits) - 1)
                      for record in records)
    return gram, phases, kernel, logical_z


def _product_expectation(description, mask, address):
    gram, phases, kernel, logical_z = description
    remainder = mask
    for relation in reversed(kernel):
        if remainder & (1 << (relation.bit_length() - 1)):
            remainder ^= relation
    if remainder:
        return 0j
    active = [i for i in range(len(phases)) if mask >> i & 1]
    exponent = (sum(phases[i] for i in active)
                + 2 * sum(gram[i][j] for i in active for j in active if i < j))
    label = 0
    for i in active:
        label ^= logical_z[i]
    return 1j ** (exponent % 4) * (-1) ** ((label & address).bit_count() % 2)


def _bias(description, address_bits, pair_signal=False):
    """Conjugate P_0 by the listed rotations, optionally multiplying a final S."""
    gram, phases, _, _ = description
    stop = len(phases) - int(pair_signal)
    coefficients = {1: 1 + 0j}
    for axis in range(1, stop):
        updated = {}
        for mask, coefficient in coefficients.items():
            anticommutes = sum(gram[axis][i] + gram[i][axis]
                               for i in range(axis) if mask >> i & 1) % 2
            if anticommutes:
                updated[mask] = updated.get(mask, 0) + coefficient / np.sqrt(2)
                enlarged = mask | (1 << axis)
                # The new P_axis moves past all earlier selected generators;
                # their total commutation sign is minus in this branch.
                updated[enlarged] = (updated.get(enlarged, 0)
                                     - 1j * coefficient / np.sqrt(2))
            else:
                updated[mask] = updated.get(mask, 0) + coefficient
        coefficients = updated
    final_mask = (1 << stop) if pair_signal else 0
    return np.array([
        sum(coefficient * _product_expectation(description, mask | final_mask, x)
            for mask, coefficient in coefficients.items()).real
        for x in range(1 << address_bits)
    ])


def _native_records(width, word, output):
    """Push actual Clifford gates to the output; retain each signed T axis."""
    frame = np.eye(1 << width, dtype=complex)
    axes = []
    for gate in word:
        if gate[0] in ('T', 'TDG'):
            z = _pauli_matrix(width, (0, 1 << gate[1], 0))
            image = frame.conj().T @ z @ frame
            axes.append(_pauli_record(image if gate[0] == 'T' else -image))
        else:
            frame = _word_matrix(width, [gate]) @ frame
    initial = _pauli_record(frame.conj().T @ output @ frame)
    # The normal form is C R_t ... R_1, so the observable first meets R_t.
    return [initial, *reversed(axes)]


def _direct_observable(width, records):
    result = _pauli_matrix(width, records[0])
    for record in records[1:]:
        rotation = (np.cos(np.pi / 8) * np.eye(1 << width)
                    - 1j * np.sin(np.pi / 8) * _pauli_matrix(width, record))
        result = rotation.conj().T @ result @ rotation
    return result


class RawDecisionCountTests(unittest.TestCase):
    def test_ordered_cocycle_retains_imaginary_phases_and_dependencies(self):
        # X_0, Z_0, Y_0 are dependent as Pauli labels. Their ordered product
        # is -i I; treating all products as Hermitian would lose this sign.
        records = [(1, 0, 0), (0, 1, 0), (1, 1, 1),
                   (2, 2, 3), (1, 2, 2)]
        description = _descriptor(records, 2)
        self.assertEqual(_product_expectation(description, 7, 0), -1j)
        nonzero_values = set()
        for mask in range(1 << len(records)):
            direct = np.eye(4, dtype=complex)
            for index, record in enumerate(records):
                if mask >> index & 1:
                    direct = direct @ _pauli_matrix(2, record)
            expected = np.array([_product_expectation(description, mask, x)
                                 for x in range(4)])
            np.testing.assert_array_equal(expected, np.diag(direct))
            nonzero_values.update(value for value in expected if value)
        self.assertEqual(nonzero_values, {1, -1, 1j, -1j})

    def test_native_conjugation_and_pair_signal_match_compressed_description(self):
        words = [
            [('H', 2), ('CX', 0, 2), ('T', 2), ('H', 0),
             ('CX', 2, 1), ('TDG', 1), ('S', 2), ('H', 1),
             ('CX', 1, 2), ('T', 0), ('CX', 0, 2), ('H', 2)],
            [('S', 0), ('H', 0), ('CX', 1, 0), ('TDG', 0),
             ('H', 2), ('CX', 0, 2), ('T', 1), ('S', 2),
             ('CX', 2, 1), ('T', 2), ('H', 1), ('CX', 1, 0)],
        ]
        output = _pauli_matrix(3, (0, 4, 0))
        for word in words:
            with self.subTest(word=word):
                source = _word_matrix(3, word)
                records = _native_records(3, word, output)
                self.assertEqual(len(records), 1 + sum(g[0] in ('T', 'TDG')
                                                       for g in word))
                direct = source.conj().T @ output @ source
                np.testing.assert_allclose(_direct_observable(3, records), direct,
                                           atol=3e-12)
                np.testing.assert_allclose(_bias(_descriptor(records, 2), 2),
                                           np.diag(direct)[:4].real, atol=3e-12)
                # H = Re(V^dagger S V S) needs exactly one additional
                # generator in the description; its zero-work bias is tested.
                pair_records = records + [(0, 4, 0)]
                hermitian_pair = (direct @ output + output @ direct) / 2
                np.testing.assert_allclose(
                    _bias(_descriptor(pair_records, 2), 2, pair_signal=True),
                    np.diag(hermitian_pair)[:4].real, atol=3e-12)

    def test_invisible_work_labels_do_not_add_witness_bias_functions(self):
        records = [(0, 1, 0), (1, 1, 1), (3, 2, 1),
                   (2, 1, 0), (0, 3, 0)]
        # Extra Z labels change the physical observables, but contribute +1
        # on the zero-work witness and change neither B, ker(X), nor L.
        tags = [4, 8, 12, 0, 8]
        enlarged = [(x, z | tag, phase)
                    for (x, z, phase), tag in zip(records, tags)]
        small_description = _descriptor(records, 1)
        self.assertEqual(small_description, _descriptor(enlarged, 1))
        small = _direct_observable(2, records)
        large = _direct_observable(4, enlarged)
        np.testing.assert_allclose(np.diag(small)[:2], np.diag(large)[:2],
                                   atol=3e-12)
        np.testing.assert_allclose(_bias(small_description, 1),
                                   np.diag(large)[:2].real, atol=3e-12)
        self.assertGreater(np.linalg.norm(large - np.kron(np.eye(4), small), 2), 1)


if __name__ == '__main__':
    unittest.main()
