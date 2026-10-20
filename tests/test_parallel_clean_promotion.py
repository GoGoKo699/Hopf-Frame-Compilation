"""Two bounded exact receipts for the parallel clean-source promotion.

These protect a literal native phase word and a finite packing/echo order.
The packing is specified as a basis permutation, not emitted as a native
circuit. Neither fixture supplies its synthesis cost or a depth theorem.
"""
import unittest

from tests.exact_matrices import CO, CZ, adjoint, native_columns


def inverse_native(word):
    inverse = {"H": "H", "CX": "CX", "T": "TDG", "TDG": "T"}
    return [(inverse[name], *wires) for name, *wires in reversed(word)]


def private_ccz(u, v, s, helper):
    """Literal CCZ with eight T gates, ten CNOTs, and dirty helper return."""
    gray = [("T", helper), ("CX", u, helper),
            ("TDG", helper), ("CX", v, helper),
            ("T", helper), ("CX", u, helper),
            ("TDG", helper), ("CX", v, helper)]
    return inverse_native(gray) + [("CX", s, helper)] + gray + [("CX", s, helper)]


def active_indices(label):
    # R=(A_0,t_0,A_1,t_1); bit four is the additional arbitrary z.
    return [j for j in range(2)
            if not (label >> (2*j) & 1) and (label >> (2*j+1) & 1)]


def refined_small_packing():
    """The specified k=1,h=2 permutation, without a native cost claim."""
    domain = list(range(32))
    source_one = [v for v in domain if len(active_indices(v)) == 1][:8]
    source_zero = [v for v in domain if not active_indices(v)][:8]
    target_one = [v for v in domain if (v >> 4 & 1) and (v & 1)]
    target_zero = [v for v in domain if (v >> 4 & 1) and not (v & 1)]
    image = dict(zip(source_one, target_one))
    image.update(zip(source_zero, target_zero))
    used = set(image.values())
    for parity in (0, 1):
        remaining_source = [v for v in domain if v not in image
                            and len(active_indices(v)) % 2 == parity]
        remaining_target = [v for v in domain if v not in used
                            and int(not (v >> 4 & 1) and (v & 15) < 4) == parity]
        assert len(remaining_source) == len(remaining_target)
        image.update(zip(remaining_source, remaining_target))
    return [image[v] for v in domain]


def inverse_symbolic(word):
    names = {"W": "WI", "WI": "W", "G": "G", "C": "C", "X": "X"}
    return [names[name] for name in reversed(word)]


def apply_symbolic(word, label, permutation, inverse):
    phase = 1
    for name in word:
        low, control = label & 31, label & 32
        if name == "W":
            label = permutation[low] | control
        elif name == "WI":
            label = inverse[low] | control
        elif name == "G":
            phase *= -1 if control and len(active_indices(low)) % 2 else 1
        elif name == "C":
            label ^= 1 if label >> 4 & 1 else 0
        elif name == "X":
            label ^= 16
        else:
            raise AssertionError(name)
    return label, phase


class ParallelCleanPromotionTests(unittest.TestCase):
    def test_private_quartic_phase_all_columns_and_actual_inverse(self):
        # Inputs a,b,c,d and two entirely arbitrary borrowed bits rho,w.
        rho, helper = 4, 5
        query = [("H", rho)] + private_ccz(0, 1, rho, helper) + [("H", rho)]
        phase = private_ccz(rho, 2, 3, helper)
        word = query + phase + query + phase
        expected = [[CZ for _ in range(64)] for _ in range(64)]
        for column in range(64):
            expected[column][column] = -CO if column & 15 == 15 else CO
        actual = native_columns(6, word)
        self.assertEqual(actual, expected)
        self.assertEqual(native_columns(6, inverse_native(word)), adjoint(actual))
        for direction in (word, inverse_native(word)):
            self.assertEqual(sum(g[0] in ("T", "TDG") for g in direction), 32)
            self.assertEqual(sum(g[0] not in ("T", "TDG") for g in direction), 44)
            self.assertTrue(all(g[1] == helper for g in direction
                                if g[0] in ("T", "TDG")))
            self.assertTrue(all(g[-1] in (rho, helper) for g in direction))

    def test_refined_packing_and_echo_all_basis_labels(self):
        permutation = refined_small_packing()
        self.assertEqual(sorted(permutation), list(range(32)))
        inverse = [permutation.index(v) for v in range(32)]
        for target in range(32):
            before = inverse[target]
            expected_membership = bool((target >> 4 & 1) and (target & 1)
                                       or not (target >> 4 & 1) and (target & 15) < 4)
            self.assertEqual(bool(len(active_indices(before)) % 2), expected_membership)
            if target >> 4 & 1:
                self.assertEqual(len(active_indices(before)), target & 1)

        # Chronological words: P'=W G W†; E=P'† C P' C; D=E† X E X.
        packed = ["WI", "G", "W"]
        echo = ["C"] + packed + ["C"] + inverse_symbolic(packed)
        final = ["X"] + echo + ["X"] + inverse_symbolic(echo)
        for label in range(64):
            f, z = (label >> 5) & 1, (label >> 4) & 1
            self.assertEqual(apply_symbolic(echo, label, permutation, inverse),
                             (label, (-1) ** (f*z)))
            self.assertEqual(apply_symbolic(final, label, permutation, inverse),
                             (label, (-1) ** f))
            self.assertEqual(apply_symbolic(inverse_symbolic(final), label, permutation, inverse),
                             (label, (-1) ** f))
        self.assertEqual(final.count("G"), 4)
        self.assertEqual(final.count("W") + final.count("WI"), 8)
        self.assertEqual(final.count("C"), 4)
        self.assertEqual(final.count("X"), 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
