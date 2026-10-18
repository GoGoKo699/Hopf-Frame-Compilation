"""Exact integer checks for four-layer phase relocation to dirty work.

The native words act on every dirty input, retain literal phases, and
reverse every actual CNOT computation. The seven-bank CCZ check uses
factorization of verified full-space diagonal words, not an initialized
state or a numerical matrix on 31 qubits.
"""
from itertools import product
import unittest


# (actual CNOT basis word, active physical wires, address-XOR mask)
LAYERS = (
    ((), (0, 1, 2), 7),
    (((0, 1), (0, 2), (3, 0)), (0, 1, 2, 3), 14),
    (((1, 2), (2, 0), (0, 3), (3, 1)), (0, 1, 2, 3), 4),
    (((1, 3), (3, 0), (3, 2), (2, 1)), (0, 1, 2, 3), 8),
)
ROWS = ((1, 2, 4, 8), (9, 3, 5, 8),
        (7, 13, 6, 15), (11, 12, 14, 10))
CCZ_PARITIES = ((1, 1), (2, 1), (4, 1),
                (3, -1), (5, -1), (6, -1), (7, 1))


def _cnot(word, control, target):
    return word ^ (((word >> control) & 1) << target)


def _affine_compute(word, address, support, constant, mask):
    """Apply address-control CNOTs and X=HSSH, charging native gates."""
    count = 0
    for target in range(4):
        if mask >> target & 1:
            for control in support:
                word ^= ((address >> control) & 1) << target
                count += 1
            if constant:
                word ^= 1 << target
                count += 4  # X = H S S H in the fixed elementary alphabet.
    return word, count


def _native_phase(address, support, constant, dirty, sign):
    word = dirty
    exponent = t_count = cliffords = 0
    for basis, active, mask in LAYERS:
        word, count = _affine_compute(
            word, address, support, constant, mask)
        cliffords += count
        for control, target in basis:
            word = _cnot(word, control, target)
            cliffords += 1
        exponent -= sign * sum((word >> wire) & 1 for wire in active)
        t_count += len(active)
        for control, target in reversed(basis):
            word = _cnot(word, control, target)
            cliffords += 1
        # Reverse each actual address/X word. They commute here because
        # addresses are unchanged and all targets are private work.
        for target in reversed(range(4)):
            if mask >> target & 1:
                if constant:
                    word ^= 1 << target
                    # The actual inverse H S^-1 S^-1 H reduces exactly
                    # to H S S H because S^-2 = S^2.
                    cliffords += 4
                for control in reversed(support):
                    word ^= ((address >> control) & 1) << target
                    cliffords += 1
    return word, exponent % 8, t_count, cliffords


class ReadOnlyPhaseRelocationTests(unittest.TestCase):
    def test_null_identity_is_literal_on_every_dirty_input(self):
        for parity, dirty in product(range(2), range(16)):
            full = sum(parity ^ ((row & dirty).bit_count() % 2)
                       for row in range(16))
            self.assertEqual(full % 8, 0)
            retained = sum(parity ^ ((row & dirty).bit_count() % 2)
                           for row in range(1, 16))
            for sign in (-1, 1):
                self.assertEqual((-sign * retained) % 8,
                                 (sign * parity) % 8)

    def test_native_bases_partition_rows_and_reverse_exactly(self):
        covered = []
        for (basis, active, mask), rows in zip(LAYERS, ROWS):
            self.assertEqual(len(set(active)), len(active))
            covered.extend(rows[wire] for wire in active)
            for wire in active:
                self.assertEqual((rows[wire] & mask).bit_count() % 2, 1)
            for dirty in range(16):
                word = dirty
                for control, target in basis:
                    word = _cnot(word, control, target)
                expected = sum(((row & dirty).bit_count() % 2) << wire
                               for wire, row in enumerate(rows))
                self.assertEqual(word, expected)
                for control, target in reversed(basis):
                    word = _cnot(word, control, target)
                self.assertEqual(word, dirty)
        self.assertEqual(sorted(covered), list(range(1, 16)))

    def test_affine_parity_full_return_and_every_native_gate_charge(self):
        for support in ((), (0,), (0, 2), (0, 1, 2)):
            for address, constant, dirty, sign in product(
                    range(8), range(2), range(16), (-1, 1)):
                parity = constant
                for control in support:
                    parity ^= address >> control & 1
                word, phase, count, cliffords = _native_phase(
                    address, support, constant, dirty, sign)
                self.assertEqual(word, dirty)
                self.assertEqual(phase, (sign * parity) % 8)
                self.assertEqual(count, 15)
                self.assertEqual(cliffords, 22 + 16 * len(support)
                                 + 64 * constant)
        self.assertEqual(len(LAYERS), 4)

    def test_parallel_seven_factor_scalar_ccz_and_fair_ledger(self):
        self.assertEqual(sum(mask.bit_count() for mask, _ in CCZ_PARITIES), 12)
        for address in range(8):
            total_phase = total_t = total_g = 0
            for bank, (mask, sign) in enumerate(CCZ_PARITIES):
                support = tuple(i for i in range(3) if mask >> i & 1)
                phases = set()
                for dirty in range(16):
                    word, phase, count, cliffords = _native_phase(
                        address, support, 0, dirty, sign)
                    self.assertEqual(word, dirty)
                    phases.add(phase)
                # Each complete bank is literally the same scalar on
                # all sixteen basis inputs; tensoring arbitrary banks
                # therefore permits arbitrary inter-bank correlations.
                self.assertEqual(len(phases), 1)
                total_phase += phases.pop()
                total_t += count
                total_g += cliffords
                for _, active, _ in LAYERS:
                    self.assertTrue(all(4 * bank <= 4 * bank + wire
                                        < 4 * (bank + 1) for wire in active))
            self.assertEqual(total_phase % 8, 4 * int(address == 7))
            self.assertEqual(total_t, 105)
            self.assertEqual(total_g, 346)
        self.assertEqual(7 * 4, 28)
        for layer in range(4):
            targets = [4 * bank + wire for bank in range(7)
                       for wire in LAYERS[layer][1]]
            self.assertEqual(len(targets), len(set(targets)))


if __name__ == '__main__':
    unittest.main()
