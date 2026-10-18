"""Finite checks for uniform transport from bounded Pauli deficiency.

The proof supplies the arbitrary-width theorem. These fixtures check signed
localization of a degenerate domain, heavy-support phase reconstruction,
and a genuinely approximate endpoint. They are endpoint-mechanism checks,
not returned-query constructions or free-source compiler implementations.
"""
from __future__ import annotations

import functools
import itertools
import math
import unittest

import numpy as np

from tests.test_outer_axis_normalizer import _tau
from tests.test_read_only_three_layer_rigidity import H, S, T, X, Z, _cx, _wire


@functools.lru_cache(maxsize=None)
def _paulis(width):
    return np.array([_tau(p, width) for p in range(1 << (2 * width))])


def _hs(matrix):
    return np.linalg.norm(matrix, "fro") / math.sqrt(len(matrix))


def _coefficients(matrix, width):
    return np.einsum("pij,ji->p", _paulis(width), matrix) / (1 << width)


def _image(unitary, label, width):
    image = unitary @ _tau(label, width) @ unitary.conj().T
    coefficients = _coefficients(image, width)
    result = int(np.argmax(abs(coefficients)))
    sign = int(coefficients[result].real < 0)
    if np.linalg.norm(image - (-1) ** sign * _tau(result, width), 2) > 2e-11:
        raise ValueError("The literal image is not a signed Pauli")
    return result, sign


def _active(label, width):
    """Local one-qubit label, with wire zero the active subsystem."""
    return (label & 1) | (((label >> width) & 1) << 1)


def _spectator(label, width):
    return label & ~(1 | (1 << width))


def _embed(label, width):
    return (label & 1) | (((label >> 1) & 1) << width)


def _basis_with_coordinates(vectors):
    basis, coordinates = [], {0: 0}
    for vector in vectors:
        if vector not in coordinates:
            bit = 1 << len(basis)
            basis.append(vector)
            coordinates.update({old ^ vector: coefficient | bit
                                for old, coefficient in list(coordinates.items())})
    return basis, coordinates


def _pair(left, right, width):
    mask = (1 << width) - 1
    return (((left & mask) & (right >> width)).bit_count()
            + ((right & mask) & (left >> width)).bit_count()) & 1


class BoundedDeficiencyTransport(unittest.TestCase):
    def _check_recovery(self, local_source, inner, target, width=3):
        """Check all probes, including the discarded coefficient tail."""
        source = _wire(local_source, 0, width)
        endpoint = source @ inner @ source.conj().T
        endpoint_error = np.linalg.norm(target - endpoint, 2)
        delta = 2 * endpoint_error
        rho = 1 / (math.sqrt(2) * 3)  # a=1: 1/[sqrt(2)(2^a+1)].
        self.assertLess(delta, rho)
        images = [_image(inner, p, width) for p in range(1 << (2 * width))]
        max_heavy = max_basis = 0
        nonreal_phase = False
        max_local_error = 0

        for probe in range(1 << (2 * width)):
            qa, qb = _active(probe, width), _spectator(probe, width)
            coefficients = _coefficients(local_source.conj().T @ _tau(qa, 1)
                                         @ local_source, 1)
            np.testing.assert_allclose(coefficients.imag, 0, atol=2e-14)
            coefficients = coefficients.real
            self.assertAlmostEqual(float(np.dot(coefficients, coefficients)), 1)
            heavy = [s for s, coefficient in enumerate(coefficients)
                     if abs(coefficient) > rho]
            self.assertTrue(heavy)
            max_heavy = max(max_heavy, len(heavy))
            result, result_sign = _image(target, probe, width)
            result_a, result_b = _active(result, width), _spectator(result, width)
            base = heavy[0]
            basis, coordinates = _basis_with_coordinates([s ^ base for s in heavy])
            max_basis = max(max_basis, len(basis))
            self.assertLessEqual(len(basis), 2)

            base_label, base_sign = images[_embed(base, width) | qb]
            self.assertEqual(_spectator(base_label, width), result_b)
            basis_images = []
            for difference in basis:
                label, sign = images[_embed(difference, width)]
                self.assertEqual(_spectator(label, width), 0)
                basis_images.append((-1) ** sign * _tau(_active(label, width), 1))

            local = np.zeros((2, 2), dtype=complex)
            direct = np.zeros((2, 2), dtype=complex)
            full_native = np.zeros_like(source)
            for s, coefficient in enumerate(coefficients):
                label, sign = images[_embed(s, width) | qb]
                full_native += coefficient * (-1) ** sign * _tau(label, width)
                if s not in heavy:
                    continue
                self.assertEqual(_spectator(label, width), result_b)
                direct_term = (-1) ** sign * _tau(_active(label, width), 1)
                direct += coefficient * direct_term
                input_product = np.eye(2, dtype=complex)
                output_product = np.eye(2, dtype=complex)
                word = coordinates[s ^ base]
                for j, difference in enumerate(basis):
                    if (word >> j) & 1:
                        input_product = input_product @ _tau(difference, 1)
                        output_product = output_product @ basis_images[j]
                input_product = input_product @ _tau(base, 1)
                alpha = np.trace(input_product.conj().T @ _tau(s, 1)) / 2
                self.assertAlmostEqual(abs(alpha), 1)
                nonreal_phase |= abs(alpha.imag) > 0.5
                output_product = (alpha * output_product
                                  @ ((-1) ** base_sign
                                     * _tau(_active(base_label, width), 1)))
                np.testing.assert_allclose(output_product, direct_term, atol=3e-13)
                local += coefficient * output_product
            np.testing.assert_allclose(local, direct, atol=3e-13)

            full_target = source.conj().T @ ((-1) ** result_sign
                                            * _tau(result, width)) @ source
            self.assertLessEqual(np.linalg.norm(full_native - full_target, 2),
                                 delta + 2e-13)
            true_local = ((-1) ** result_sign * local_source.conj().T
                          @ _tau(result_a, 1) @ local_source)
            local_error = _hs(local - true_local)
            max_local_error = max(max_local_error, local_error)
            self.assertLessEqual(local_error, delta + 2 * rho + 2e-13)
            self.assertLess(local_error, 1 / math.sqrt(2))
            candidates = [(_hs(local - (-1) ** sign * local_source.conj().T
                               @ _tau(label, 1) @ local_source), label, sign)
                          for label in range(4) for sign in range(2)]
            recovered = [item for item in candidates if item[0] < 1 / math.sqrt(2)]
            self.assertEqual(len(recovered), 1)
            self.assertEqual(recovered[0][1:], (result_a, result_sign))

        return {"endpoint_error": endpoint_error, "heavy": max_heavy,
                "basis": max_basis, "nonreal_phase": nonreal_phase,
                "local_error": max_local_error}

    def test_degenerate_domain_localizes_with_literal_signs(self):
        width = 3
        input_coordinates = _cx(1, 0, width) @ _wire(H, 2, width)
        unsigned_output = _cx(0, 2, width) @ _wire(H, 0, width)
        output_coordinates = _wire(Z, 1, width) @ unsigned_output
        active_source = _wire(T, 0, width)
        source = output_coordinates @ active_source @ input_coordinates
        domain = []
        for label in range(1 << (2 * width)):
            try:
                _image(source, label, width)
                domain.append(label)
            except ValueError:
                pass
        radical = [p for p in domain if all(_pair(p, q, width) == 0 for q in domain)]
        self.assertEqual(len(domain), 1 << (2 * width - 1))  # k=1.
        self.assertEqual(len(radical), 2)                    # r=1, a=1.
        self.assertTrue(all(p ^ q in domain for p in domain for q in domain))
        localized = output_coordinates.conj().T @ source @ input_coordinates.conj().T
        np.testing.assert_allclose(localized, active_source, atol=2e-14)
        for label in range(1 << (2 * width)):
            if _active(label, width) == 0:
                np.testing.assert_allclose(localized @ _tau(label, width),
                                           _tau(label, width) @ localized, atol=2e-14)
        wrong = unsigned_output.conj().T @ source @ input_coordinates.conj().T
        self.assertEqual(_image(wrong, 1 << 1, width), (1 << 1, 1))
        self.assertGreater(np.linalg.norm(wrong @ _tau(1 << 1, width)
                                         - _tau(1 << 1, width) @ wrong, 2), 1)

    def test_heavy_difference_basis_retains_nonreal_pauli_phases(self):
        width = 3
        y = _tau(3, 1)
        axis = (X + y + Z) / math.sqrt(3)
        angle = math.pi / 7
        local_source = math.cos(angle) * np.eye(2) + 1j * math.sin(angle) * axis
        first_band = _wire(T.conj().T, 0, width)
        second_band = _wire(H @ T, 0, width)
        records = []
        for u, v in itertools.product((0, 1), repeat=2):
            # Two actual one-band Pauli images and a middle linear Pauli.
            first = first_band @ _tau(u, width) @ first_band.conj().T
            second = second_band @ _tau(u, width) @ second_band.conj().T
            inner = first @ _tau(v << 1, width) @ second
            records.append(self._check_recovery(local_source, inner, inner, width))
        self.assertEqual(max(record["heavy"] for record in records), 3)
        self.assertEqual(max(record["basis"] for record in records), 2)
        self.assertTrue(any(record["nonreal_phase"] for record in records))
        self.assertLess(max(record["endpoint_error"] for record in records), 2e-14)
        self.assertLess(max(record["local_error"] for record in records), 2e-14)

    def test_degenerate_source_can_have_entangling_endpoint(self):
        width = 3
        cz = _wire(H, 1, width) @ _cx(0, 1, width) @ _wire(H, 1, width)
        inner = _wire(H, 2, width) @ cz
        self.assertNotEqual(_spectator(_image(inner, 1, width)[0], width), 0)
        record = self._check_recovery(T, inner, inner, width)
        self.assertEqual(record["heavy"], 2)
        self.assertEqual(record["basis"], 1)
        self.assertTrue(record["nonreal_phase"])

    def test_near_clifford_source_does_not_round_the_endpoint_exactly(self):
        width = 3
        axis = (X + math.sqrt(2) * _tau(3, 1) + math.pi * Z) / math.sqrt(3 + math.pi**2)
        first_band = _wire(T.conj().T, 0, width)
        second_band = _cx(0, 1, width) @ _wire(H @ T, 0, width)
        errors = []
        for angle in (0.01, 0.001):
            local_source = math.cos(angle) * np.eye(2) + 1j * math.sin(angle) * axis
            source = _wire(local_source, 0, width)
            records, nonclifford_endpoint = [], False
            for u, v in itertools.product((0, 1), repeat=2):
                first = first_band @ _tau(u, width) @ first_band.conj().T
                second = second_band @ _tau(v, width) @ second_band.conj().T
                inner = first @ _tau(u << (width + 1), width) @ second
                records.append(self._check_recovery(local_source, inner, inner, width))
                endpoint = source @ inner @ source.conj().T
                for probe in range(1 << (2 * width)):
                    try:
                        _image(endpoint, probe, width)
                    except ValueError:
                        nonclifford_endpoint = True
                        break
            self.assertTrue(nonclifford_endpoint)
            self.assertEqual(max(record["heavy"] for record in records), 1)
            self.assertGreater(max(record["local_error"] for record in records), 1e-5)
            error = max(record["endpoint_error"] for record in records)
            self.assertGreater(error, 1e-5)
            # The proof starts from epsilon=2 eta, hence delta=2 epsilon.
            self.assertLess(error / 2, 1 / (4 * math.sqrt(2) * 3))
            errors.append(error)
        self.assertLess(errors[1], errors[0] / 5)


if __name__ == "__main__":
    unittest.main()
