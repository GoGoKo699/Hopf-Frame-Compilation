"""Exact checks for a broad four-path Bell/swap/relative-T source.

The accepted Gram is normalized but has the wrong complete sectors. These
fixtures check native controls, literal inverses and the coherent pair witness;
they do not infer an asymptotic source or lower bound from small matrices.
"""
import unittest

from tests.exact_matrices import (
    CI, CO, CZ, Gate, add, adjoint, apply_word, eye, multiply, scale, zeros,
)
from tests.test_corank_fourier_bell_source import _kernel, _w_word, _x
from tests.test_half_bilinear_gram_source import _private_toffoli
from tests.test_majorana_walsh_query import _inverse
from tests.test_phase_signal_query_limits import _counts
from tests.test_signal_interleaved_kernel_source import _ct_word, _root_powers


def _native(word, state):
    return apply_word([Gate(name, tuple(wires)) for name, *wires in word], state)


def _inner_word():
    # y=0,z=1; clean signals a=2,b=3; eight arbitrary private phase bits.
    # The Fredkin reuses a returned controlled-T helper at its earlier stage.
    fredkin = ([('CX', 0, 1)] + _private_toffoli(2, 1, 0, 4)
               + [('CX', 0, 1)])
    def controlled_t(target, helpers):
        mapping = (2, target, *helpers)
        return [(name, *(mapping[wire] for wire in wires))
                for name, *wires in _ct_word()]
    cy = controlled_t(0, (4, 5, 6, 7))
    cz = controlled_t(1, (8, 9, 10, 11))
    return ([('H', 2)] + fredkin + [('H', 2)] + cy + _inverse(cz)
            + [('H', 2), ('H', 3)])


def _inner_matrices(n):
    q, powers = 1 << n, _root_powers()
    dimension = q * q
    swap, phase = zeros(dimension), zeros(dimension)
    for column in range(dimension):
        y, z = column & (q - 1), column >> n
        swap[z | (y << n)][column] = CO
        phase[column][column] = powers[(y.bit_count() - z.bit_count()) % 8]
    accepted = scale(add(add(add(eye(dimension), swap), phase),
                         scale(multiply(phase, swap), -CO)), CO / 4)
    gram = add(scale(eye(dimension), CO / 4),
               scale(multiply(add(adjoint(phase), scale(phase, -CO)), swap), CO / 8))
    return accepted, gram


class MultipathBellSourceTests(unittest.TestCase):
    def test_native_private_controls_normalization_and_actual_inverse(self):
        word, expected = _inner_word(), _inner_matrices(1)[0]
        self.assertEqual(_counts(word), (148, 192))
        reverse = _inverse(word)
        self.assertEqual(_counts(reverse), (148, 192))
        # Both extreme dirty assignments, all principal/signal columns.
        # Universal private-helper return is the retained primitive theorem.
        for helper_bits in (0, 255):
            helper = helper_bits << 4
            for logical in range(16):
                column = helper | logical
                result = _native(word, {column: CO})
                self.assertTrue(all(row >> 4 == helper_bits for row in result))
                self.assertEqual(_native(reverse, result), {column: CO})
                if logical < 4:
                    accepted = {row & 3: value for row, value in result.items()
                                if not (row & 12)}
                    for row in range(4):
                        self.assertEqual(accepted.get(row, CZ), expected[row][logical])

    def test_complete_gram_pair_witness_and_bell_charge_orientation(self):
        accepted, gram = _inner_matrices(2)
        self.assertEqual(multiply(adjoint(accepted), accepted), gram)
        # y=11,z=00 and its swap have phases +i,-i. Their charge is
        # identical, but the Gram has eigenvalues zero and one half.
        u, v = 3, 12
        self.assertEqual(gram[u][u], CO / 4)
        self.assertEqual(gram[v][v], CO / 4)
        self.assertEqual(gram[u][v], -CI / 4)
        self.assertEqual(gram[v][u], CI / 4)
        for charge in (1, 2, 3):
            desired = CZ if (charge & 3).bit_count() % 2 else CO / 2
            defect = [[gram[i][j] - (desired if i == j else CZ)
                       for j in (u, v)] for i in (u, v)]
            # On this pair the spectrum is {0,-1/2} or {0,+1/2}.
            square = multiply(defect, defect)
            sign = 1 if desired == CZ else -1
            self.assertEqual(square, scale(defect, sign * CO / 2))
            self.assertEqual(defect[0][0] + defect[1][1], sign * CO / 2)
        for address in (0, 1):
            rows, vector, n = _kernel(1, address)
            self.assertEqual(n, 2)
            w2 = (_w_word(rows, [0, 1], [4, 5])
                  + _w_word(rows, [2, 3], [6, 7]))
            charge_word = [gate for offset in (0, 2) for bit in range(2)
                           if (vector >> bit) & 1 for gate in _x(offset + bit)]
            conjugated = w2 + charge_word + _inverse(w2)
            for column in (0, 1, 16, 80):
                y, z = (column >> 4) & 3, (column >> 6) & 3
                sign = -1 if ((y ^ z) & vector).bit_count() % 2 else 1
                self.assertEqual(_native(conjugated, {column: CO}),
                                 {column: sign * CO})


if __name__ == '__main__':
    unittest.main()
