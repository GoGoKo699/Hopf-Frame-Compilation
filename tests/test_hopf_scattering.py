"""Small algebraic checks of the complete Hopf tree scattering identity.

There are at most sixteen ports plus one borrowed permutation helper.
Only the port permutations are checked as literal Clifford+T words; the
coin rotations remain ideal.  Eliminating internal ports is not a circuit.
"""
from __future__ import annotations

from math import comb
import unittest

import numpy as np

from compiler_robust_hopf.conventions import marker_label
from compiler_robust_hopf.frames import direct_real_frame, hopf_ry
from tests.test_operator_source_compiler import _expand_toffolis, _word_matrix


ATOL = 3e-12


def _scattering_ports(height, angles):
    """Return S in physical order (external x, internal heap node v).

    External x=0 enters the root's first port; external marker lambda(v)
    enters the second port of node v.  Internal v=0,1 are disconnected
    dummy ports.  Coin order is (v,p), with integer label 2*v+p.
    """
    size = 1 << height
    pin = np.zeros((2 * size, 2 * size))
    pin[2, 0] = 1
    for node in range(1, size):
        pin[2 * node + 1, marker_label(node, height)] = 1
        if node >= 2:
            pin[2 * node, size + node] = 1
    pin[0, size] = pin[1, size + 1] = 1

    coin = np.zeros_like(pin)
    coin[:2, :2] = -np.eye(2)
    for node, theta in enumerate(angles, start=1):
        coin[2 * node:2 * node + 2, 2 * node:2 * node + 2] = hopf_ry(theta).real

    # The child label is 2*v+p.  Its high bit distinguishes an internal
    # child from a leaf, with the opposite convention to the physical flag.
    pout = np.zeros_like(pin)
    for label in range(2 * size):
        pout[label ^ size, label] = 1
    scattering = pout @ coin @ pin
    d, c = scattering[:size, :size], scattering[:size, size:]
    b, a = scattering[size:, :size], scattering[size:, size:]
    return dict(pin=pin, coin=coin, pout=pout, s=scattering,
                a=a, b=b, c=c, d=d)


def _external_transfer(ports):
    a, b, c, d = (ports[key] for key in ('a', 'b', 'c', 'd'))
    return d + c @ np.linalg.solve(np.eye(len(a)) - a, b)


def _right_path_angles(height, theta):
    angles = np.zeros((1 << height) - 1)
    for depth in range(height):
        node = (1 << (depth + 1)) - 1
        angles[node - 1] = theta
    return angles


def _input_port_word(height):
    """The bit-reversal algorithm, using one dirty helper for height<=3."""
    word, helper = [], height + 1

    def swap(u, v, controls=()):
        flips = [('X', bit) for bit, value in controls if value == 0]
        word.extend(flips + [('CX', u, v)])
        control = [bit for bit, _ in controls] + [v]
        if len(control) <= 2:
            word.append(('CX' if len(control) == 1 else 'CCX', *control, u))
        else:
            assert len(control) == 3
            first, second, third = control
            word.extend([('CCX', first, second, helper),
                         ('CCX', helper, third, u)] * 2)
        word.extend([('CX', u, v)] + flips)

    # [a,x] -> [x,1-a].  The low bit is now the coin port p.
    for bit in range(height, 0, -1):
        swap(bit, bit - 1)
    word.append(('X', 0))
    for bit in range(height // 2):
        swap(1 + bit, height - bit, ((0, 1),))
    for highest in range(height):
        controls = [(0, 1), (highest + 1, 1)]
        controls += [(bit + 1, 0) for bit in range(highest + 1, height)]
        for bit in range(highest // 2):
            swap(1 + bit, highest - bit, controls)
    # Transpose only the complete coin labels 1 and 2.
    swap(0, 1, [(bit, 0) for bit in range(2, height + 1)])
    return word


class HopfScatteringTests(unittest.TestCase):
    def test_literal_port_permutations_return_arbitrary_helper(self):
        for height in (2, 3):
            ports = _scattering_ports(height, np.zeros((1 << height) - 1))
            native = _expand_toffolis(_input_port_word(height))
            np.testing.assert_allclose(_word_matrix(height + 2, native),
                                       np.kron(np.eye(2), ports['pin']), atol=ATOL, rtol=0)
            np.testing.assert_array_equal(_word_matrix(height + 2, [('X', height)]),
                                          np.kron(np.eye(2), ports['pout']))

    def test_complete_external_frame_and_unitary_physical_coin(self):
        for height in (2, 3):
            size = 1 << height
            fixtures = (
                np.linspace(.17, .83, size - 1),
                np.zeros(size - 1),
                np.full(size - 1, np.pi / 2),
                np.array([0 if node % 2 else np.pi / 2
                          for node in range(1, size)]),
                np.linspace(-1.2, .9, size - 1),
            )
            for angles in fixtures:
                with self.subTest(height=height, angles=angles):
                    ports = _scattering_ports(height, angles)
                    for name in ('pin', 'coin', 'pout', 's'):
                        operator = ports[name]
                        np.testing.assert_allclose(operator.T @ operator,
                                                   np.eye(2 * size),
                                                   atol=ATOL, rtol=0)
                    # This reference composes addressed depth rotations;
                    # it does not build or eliminate the scattering ports.
                    actual = _external_transfer(ports)
                    np.testing.assert_allclose(actual,
                                               direct_real_frame(height, angles),
                                               atol=ATOL, rtol=0)
                    np.testing.assert_allclose(actual.T @ actual, np.eye(size),
                                               atol=ATOL, rtol=0)

    def test_physical_nilpotence_and_disconnected_dummy_inverse(self):
        for height in (2, 3):
            size = 1 << height
            ports = _scattering_ports(height, np.linspace(.23, .79, size - 1))
            a, b, c, d = (ports[key] for key in ('a', 'b', 'c', 'd'))
            physical = a[2:, 2:]
            np.testing.assert_array_equal(a[:2, :2], -np.eye(2))
            np.testing.assert_array_equal(a[:2, 2:], 0)
            np.testing.assert_array_equal(a[2:, :2], 0)
            np.testing.assert_array_equal(b[:2], 0)
            np.testing.assert_array_equal(c[:, :2], 0)
            np.testing.assert_allclose(np.linalg.matrix_power(physical, height - 1),
                                       0, atol=ATOL, rtol=0)

            finite_inverse = sum(np.linalg.matrix_power(physical, power)
                                 for power in range(height - 1))
            np.testing.assert_allclose((np.eye(size - 2) - physical) @ finite_inverse,
                                       np.eye(size - 2), atol=ATOL, rtol=0)
            padded_inverse = np.zeros_like(a)
            padded_inverse[:2, :2] = .5 * np.eye(2)
            padded_inverse[2:, 2:] = finite_inverse
            np.testing.assert_allclose((np.eye(size) - a) @ padded_inverse,
                                       np.eye(size), atol=ATOL, rtol=0)
            np.testing.assert_allclose(d + c[:, 2:] @ finite_inverse @ b[2:],
                                       _external_transfer(ports), atol=ATOL, rtol=0)

            # Only the physical internal block is nilpotent.  Applying its
            # finite-series inverse formula to padded A leaves a dummy defect.
            self.assertGreaterEqual(np.linalg.norm(np.linalg.matrix_power(a, height - 1)), 1)
            invalid_inverse = sum(np.linalg.matrix_power(a, power)
                                  for power in range(height - 1))
            self.assertGreater(np.linalg.norm((np.eye(size) - a) @ invalid_inverse
                                             - np.eye(size)), 1)

    def test_right_path_has_the_analytic_degree_height_fourier_coefficient(self):
        for height in (2, 3):
            # A deterministic Fourier grid recovers every coefficient of the
            # known degree-height polynomial; this is not a random sweep.
            grid = 2 * np.pi * np.arange(2 * height + 1) / (2 * height + 1)
            amplitudes = np.array([
                _external_transfer(_scattering_ports(height,
                                                     _right_path_angles(height, theta)))[-1, 0]
                for theta in grid
            ])
            np.testing.assert_allclose(amplitudes, np.sin(grid) ** height,
                                       atol=ATOL, rtol=0)
            expected = {height - 2 * k: (-1) ** k * comb(height, k) / (2j) ** height
                        for k in range(height + 1)}
            for frequency in range(-height, height + 1):
                coefficient = np.mean(amplitudes * np.exp(-1j * frequency * grid))
                self.assertAlmostEqual(abs(coefficient - expected.get(frequency, 0)),
                                       0, delta=ATOL)
            self.assertAlmostEqual(abs(expected[height]), 2 ** (-height), delta=ATOL)


if __name__ == '__main__':
    unittest.main()
