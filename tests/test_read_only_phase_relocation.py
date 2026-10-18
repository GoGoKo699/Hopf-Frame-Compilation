"""Exact integer checks for four-layer phase relocation to dirty work.

The native words act on every dirty input, retain literal phases, and
reverse every actual CNOT computation. The scalar CCZ checks factor
verified diagonal words on arbitrary private banks; no initialized
work state or large dense numerical matrix is used.
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


def _hybrid_ccz(address, shared, dirty_banks, c_support=(2,)):
    """Actual four-layer 64T word, including shared-w control ordering."""
    moved = (((), -1), ((0,), 1), ((1,), 1), (c_support, 1))
    retained = (((0, 1), -1), ((0,) + c_support, -1),
                ((1,) + c_support, -1), ((0, 1) + c_support, 1))
    banks = list(dirty_banks)
    exponent = t_count = cliffords = 0
    layers = []
    for layer, (basis, active, mask) in enumerate(LAYERS):
        original = (shared, tuple(banks))
        # All bank computations read the original shared w.
        for bank, (support, _) in enumerate(moved):
            for wire in range(4):
                if mask >> wire & 1:
                    banks[bank] ^= shared << wire
                    cliffords += 1
                    for control in support:
                        banks[bank] ^= ((address >> control) & 1) << wire
                        cliffords += 1
            for control, target in basis:
                banks[bank] = _cnot(banks[bank], control, target)
                cliffords += 1
        retained_support, retained_sign = retained[layer]
        for control in retained_support:
            shared ^= address >> control & 1
            cliffords += 1
        # Physical wire 0 is w; four private banks begin at wire 1.
        targets = [0]
        exponent += retained_sign * shared
        t_count += 1
        for bank, (_, sign) in enumerate(moved):
            exponent -= sign * sum(banks[bank] >> wire & 1 for wire in active)
            t_count += len(active)
            targets.extend(1 + 4 * bank + wire for wire in active)
        # Reverse the entire actual preparation, starting with w return.
        for control in reversed(retained_support):
            shared ^= address >> control & 1
            cliffords += 1
        for bank in reversed(range(4)):
            support, _ = moved[bank]
            for control, target in reversed(basis):
                banks[bank] = _cnot(banks[bank], control, target)
                cliffords += 1
            for wire in reversed(range(4)):
                if mask >> wire & 1:
                    for control in reversed(support):
                        banks[bank] ^= ((address >> control) & 1) << wire
                        cliffords += 1
                    banks[bank] ^= shared << wire
                    cliffords += 1
        layers.append((original, (shared, tuple(banks)), tuple(targets)))
    return shared, tuple(banks), exponent % 8, t_count, cliffords, layers


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

    def test_hybrid_scalar_ccz_native_ordering_and_simultaneous_ledger(self):
        for address, shared in product(range(8), range(2)):
            # Each bank phase is additive and each bank returns at every
            # layer. Varying each bank through all16 values therefore
            # checks the factorized full-space identity, not just a
            # prepared-work subspace. Shared w is also arbitrary.
            for bank, dirty in product(range(4), range(16)):
                banks = [0] * 4
                banks[bank] = dirty
                out_w, out_banks, phase, count, gates, layers = _hybrid_ccz(
                    address, shared, banks)
                self.assertEqual((out_w, out_banks), (shared, tuple(banks)))
                self.assertEqual(phase, 4 * int(address == 7))
                self.assertEqual((count, gates), (64, 218))
                self.assertEqual(len(layers), 4)
                for before, after, targets in layers:
                    self.assertEqual(before, after)
                    self.assertEqual(len(targets), len(set(targets)))
                    self.assertTrue(all(0 <= target < 17 for target in targets))
                self.assertEqual([len(item[2]) for item in layers],
                                 [13, 17, 17, 17])

    def test_hybrid_linear_row_and_r91_clifford_charge(self):
        for row_size in (0, 1, 3):
            c_support = tuple(range(2, 2 + row_size))
            for address, shared in product(range(1 << (2 + row_size)), range(2)):
                banks = (0, 5, 10, 15)
                out_w, out_banks, phase, count, gates, _ = _hybrid_ccz(
                    address, shared, banks, c_support)
                c = sum(address >> wire & 1 for wire in c_support) % 2
                expected = 4 * (address & 1) * (address >> 1 & 1) * c
                self.assertEqual(phase, expected)
                self.assertEqual((out_w, out_banks), (shared, banks))
                self.assertEqual((count, gates), (64, 196 + 22 * row_size))
        for terms in (1, 2, 5):
            for layer in range(4):
                local = [0] + [1 + 4 * bank + wire for bank in range(4)
                               for wire in LAYERS[layer][1]]
                targets = [17 * term + wire for term in range(terms)
                           for wire in local]
                self.assertEqual(len(targets), len(set(targets)))
                self.assertTrue(all(target < 17 * terms for target in targets))


if __name__ == '__main__':
    unittest.main()
