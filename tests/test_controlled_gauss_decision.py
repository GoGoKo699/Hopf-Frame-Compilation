"""Two bounded exact audits of the controlled-Gauss raw decision.

The native quartic echo is checked on every logical/helper column in
Q(sqrt(2), i). A separate 32-dimensional principal circuit checks the
whole compression with integer Walsh matrices, not a large simulation
of the all-size helper allocation. Neither fixture asserts a coherent
dirty-input decision or a growing-depth bound.
"""
from collections import Counter
import unittest

import numpy as np

from tests.exact_matrices import CO, CZ, native_columns
from tests.test_operator_source_compiler import _adjoint


def _relative_phase(u, v, dirty):
    """Cyclic Gray exposure of p(d;u,v), including its literal phase."""
    return [
        ('T', dirty), ('CX', u, dirty),
        ('TDG', dirty), ('CX', v, dirty),
        ('T', dirty), ('CX', u, dirty),
        ('TDG', dirty), ('CX', v, dirty),
    ]


def _private_toffoli(u, v, target, phase_helper):
    """Eight-phase Gray cycle; only the private helper receives T."""
    word = [('H', target)]
    # 000,100,110,010,011,111,101,001,000; phases alternate in sign.
    for step, control in enumerate((u, v, u, target, u, v, u, target)):
        word += [('T' if step % 2 == 0 else 'TDG', phase_helper),
                 ('CX', control, phase_helper)]
    return word + [('H', target)]


def _walsh(width):
    size = 1 << width
    return np.array([
        [1 if (row & column).bit_count() % 2 == 0 else -1
         for column in range(size)]
        for row in range(size)
    ], dtype=np.int64)


class ControlledGaussDecisionTests(unittest.TestCase):
    def test_native_quartic_echo_returns_both_arbitrary_helpers(self):
        c, x, y, z, dirty, phase_helper = range(6)
        f = _relative_phase(y, z, dirty)
        q = _private_toffoli(c, x, dirty, phase_helper)
        word = _adjoint(f) + q + f + _adjoint(q)

        # The complete matrix includes every input of both dirty helpers.
        # Equality to phase(control bits) times identity also covers
        # arbitrary superpositions and external reference correlations.
        expected = [[CZ for _ in range(64)] for _ in range(64)]
        for basis in range(64):
            expected[basis][basis] = -CO if basis & 15 == 15 else CO
        self.assertEqual(native_columns(6, word), expected)

        counts = Counter(gate[0] for gate in word)
        self.assertEqual(counts['T'] + counts['TDG'], 24)
        self.assertEqual((counts['CX'], counts['H']), (24, 4))
        self.assertTrue(all(gate[1] in (dirty, phase_helper)
                            for gate in word if gate[0] in ('T', 'TDG')))
        self.assertTrue(all(gate[1] == dirty
                            for gate in word if gate[0] == 'H'))
        self.assertTrue(all(gate[2] in (dirty, phase_helper)
                            for gate in word if gate[0] == 'CX'))

    def test_complete_compression_diagonals_and_fixed_product_witnesses(self):
        # A(x)=[[1,0],[x,1]], b=e_0, d=e_1, so d^T A(x)^-1 b=x.
        # This small affine pencil tests the full principal compression;
        # the retained Gauss fixtures cover the prefix-tree construction.
        m, work_size, full_size = 2, 16, 32
        work_h = _walsh(2 * m)
        outer_h = _walsh(2 * m + 1)
        signal_z = np.r_[np.ones(work_size, dtype=np.int64),
                         -np.ones(work_size, dtype=np.int64)]
        zero = np.zeros(work_size, dtype=np.int64)
        plus = np.ones(work_size, dtype=np.int64)
        minus = work_h[:, 1]  # Only y_s=y_0 is in the minus state.

        for x in (0, 1):
            with self.subTest(address=x):
                signs = []
                for work in range(work_size):
                    y0, y1, z0, z1 = [(work >> bit) & 1 for bit in range(4)]
                    phase = y0 * z0 + x * y1 * z0 + y1 * z1 + y0 + z1
                    signs.append(1 - 2 * (phase % 2))
                signs = np.array(signs, dtype=np.int64)

                # Exact rational matrices: U=u_num/16 and V=v_num/32.
                # V is built from H_all, controlled D, H_all, independently
                # of the block formula for the expected pulled-back signal.
                u_num = (work_h * signs[None, :]) @ work_h
                controlled = np.r_[np.ones(work_size, dtype=np.int64), signs]
                v_num = (outer_h * controlled[None, :]) @ outer_h
                observable_num = v_num.T @ (signal_z[:, None] * v_num)
                expected = np.zeros((full_size, full_size), dtype=np.int64)
                factor = full_size ** 2 // work_size
                expected[:work_size, :work_size] = factor * u_num
                expected[work_size:, work_size:] = -factor * u_num
                np.testing.assert_array_equal(observable_num, expected)

                sigma = 1 - 2 * x
                np.testing.assert_array_equal(
                    np.diag(u_num), np.full(work_size, (1 << m) * sigma))
                self.assertEqual(np.count_nonzero(signs == 1),
                                 (work_size + (1 << m) * sigma) // 2)
                self.assertEqual(np.count_nonzero(signs == -1),
                                 (work_size - (1 << m) * sigma) // 2)

                # Normalizing either product vector multiplies both sides
                # by 1/4. The actual full source gives opposite deterministic
                # signal answers on two fixed product work inputs at both x.
                np.testing.assert_array_equal(
                    v_num @ np.r_[plus, zero],
                    full_size * np.r_[plus, zero])
                np.testing.assert_array_equal(
                    v_num @ np.r_[minus, zero],
                    full_size * np.r_[zero, minus])


if __name__ == '__main__':
    unittest.main()
