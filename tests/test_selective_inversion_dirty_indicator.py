"""Exact native checks for selective inversion on arbitrary dirty signals.

The finite digital word is deliberately coarse and is not an accurate
indicator. These checks retain literal phases, all physical columns and
actual inverses; they do not establish the scattering formula or an
asymptotic resource bound by finite simulation.
"""
import unittest

from tests.exact_matrices import (
    CO, CZ, add, adjoint, eye, multiply, native_columns, scale, zeros,
)
from tests.test_dirty_geometric_phase_loops import _rotation_word
from tests.test_majorana_walsh_query import _inverse


def _x(wire):
    return [('H', wire), ('S', wire), ('S', wire), ('H', wire)]


def _on(word, wire):
    return [(name, wire) for name, _ in word]


def _echo(word, signals):
    z = [('S', wire) for wire in signals] * 2
    return word + z + _inverse(word) + z


def _z_pair(k, control, signal):
    return ([('CX', control, signal)]
            + _on(_rotation_word(k), signal)
            + [('CX', control, signal)])


class SelectiveInversionDirtyIndicatorTests(unittest.TestCase):
    def test_actual_adjoint_echo_norm_identities_and_literal_scalar(self):
        samples = (
            [('T', 0)],
            _x(0) + [('T', 0)],
            [('H', 0), ('T', 0), ('H', 0)],
            [('H', 0), ('T', 0), ('S', 0), ('H', 0), ('TDG', 0)],
            _rotation_word(1, 'X') + _rotation_word(1),
        )
        z = native_columns(1, [('S', 0), ('S', 0)])
        x = native_columns(1, _x(0))
        for word in samples:
            u = native_columns(1, word)
            e = native_columns(1, _echo(word, (0,)))
            expected = multiply(z, multiply(adjoint(u), multiply(z, u)))
            self.assertEqual(e, expected)
            minus = add(e, scale(eye(2), -CO))
            plus = add(e, eye(2))
            self.assertEqual(multiply(adjoint(minus), minus),
                             scale(eye(2), 4 * u[0][1].abs2()))
            self.assertEqual(multiply(adjoint(plus), plus),
                             scale(eye(2), 4 * u[0][0].abs2()))
            self.assertEqual(native_columns(1, _inverse(_echo(word, (0,)))),
                             adjoint(e))
            # A transpose word has a determinant scalar for a U(2) source.
            determinant = u[0][0] * u[1][1] - u[0][1] * u[1][0]
            transpose = [list(row) for row in zip(*u)]
            self.assertEqual(multiply(x, multiply(transpose, multiply(x, u))),
                             scale(e, determinant))
        self.assertEqual(native_columns(1, _echo(samples[0], (0,))), eye(2))
        self.assertEqual(native_columns(1, _echo(samples[1], (0,))),
                         scale(eye(2), -CO))
        # A genuine non-Clifford source need not have determinant one.
        self.assertNotEqual(native_columns(1, samples[2])[0][0], CZ)
        u = native_columns(1, samples[2])
        self.assertNotEqual(u[0][0] * u[1][1] - u[0][1] * u[1][0], CO)

    def test_complete_finite_digital_word_and_inverse_on_every_column(self):
        # p=1, m=2. Each finite step has detuning angle pi*d/2,
        # decomposed exactly as (p+1)/2 + Z_Y/2 - (-1)^u Z_x/2.
        # HTH is a native X rotation with its actual scalar retained.
        p, m, width = 1, 2, 5
        signals = (3, 4)
        pulse = []
        for u, signal in enumerate(signals):
            diagonal = _on(_rotation_word(p + 1), signal)
            diagonal += _z_pair(1, p + u, signal)
            diagonal += _z_pair(1 if u else -1, 0, signal)
            x_rotation = [('H', signal), ('T', signal), ('H', signal)]
            pulse += (diagonal + x_rotation) * 2
        word = _echo(pulse, signals)
        actual = native_columns(width, word)
        expected = zeros(1 << width)
        local = {}
        for detuning in range(p + 2):
            source = (_rotation_word(2 * detuning)
                      + [('H', 0), ('T', 0), ('H', 0)]) * 2
            local[detuning] = native_columns(1, _echo(source, (0,)))
        for column in range(1 << width):
            address = column & 1
            bank = (column >> p) & 3
            for out_signals in range(1 << m):
                row = (column & 7) | (out_signals << 3)
                amplitude = CO
                for u, signal in enumerate(signals):
                    detuning = (address ^ u).bit_count() + 1 - ((bank >> u) & 1)
                    amplitude *= local[detuning][(out_signals >> u) & 1][
                        (column >> signal) & 1]
                expected[row][column] = amplitude
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(width, _inverse(word)), adjoint(actual))
        for orientation in (word, _inverse(word)):
            self.assertEqual(sum(gate[0] in ('T', 'TDG') for gate in orientation), 8)
        target = zeros(1 << width)
        for column in range(1 << width):
            address = column & 1
            selected = (column >> (p + address)) & 1
            target[column][column] = -CO if selected else CO
        self.assertNotEqual(actual, target)


if __name__ == '__main__':
    unittest.main()
