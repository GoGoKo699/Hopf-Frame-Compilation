"""Exact scratch fixtures for the weak complete-probe Pauli-domain argument."""
from itertools import combinations, product
from pathlib import Path
import os
import sys
import unittest


def _repository_root():
    # The durable location is research/checkpoints/pr269/weak-probe. Resolve
    # that checkout first; the environment/fallback supports this scratch audit.
    for parent in Path(__file__).resolve().parents:
        if (parent / 'WORKSPACE.md').is_file() and (parent / 'tests').is_dir():
            return parent
    configured = os.environ.get('HOPF_FRAME_REPO')
    if configured:
        candidate = Path(configured).resolve()
        if (candidate / 'WORKSPACE.md').is_file():
            return candidate
    return Path('/workspace/scratch/e3eb9dbc5d0b/Hopf-Frame-Compilation-next')


sys.path.insert(0, str(_repository_root()))
from tests.test_field_decoder_pauli_domain import (
    _field_image, _hermitian_phase, _pairing, _product_phase, _rank, _span,
    _translate)
from tests.exact_matrices import Q2
from tests.test_flat_field_bank_native import full_intertwiner
from tests.test_one_band_copy_feasibility import _conjugate


def _rm_basis(p, degree):
    return [sum(int(all((x >> j) & 1 for j in subset)) << x
                for x in range(1 << p))
            for size in range(min(p, degree) + 1)
            for subset in combinations(range(p), size)]


def _gram_rank(basis, width):
    return _rank([sum(_pairing(v, w, width) << j for j, w in enumerate(basis))
                  for v in basis])


def _clifford_label(label, word, width):
    mask = (1 << width) - 1
    a, b = label & mask, label >> width
    for name, *wires in word:
        if name == 'H':
            j, = wires
            difference = ((a ^ b) >> j) & 1
            a ^= difference << j
            b ^= difference << j
        elif name == 'S':
            j, = wires
            b ^= ((a >> j) & 1) << j
        else:
            assert name == 'CX'
            control, target = wires
            a ^= ((a >> control) & 1) << target
            b ^= ((b >> target) & 1) << control
    return a | (b << width)


class WeakProbeDomainAuditTests(unittest.TestCase):
    def test_field_four_cubes_are_literal_scalars_with_helpers(self):
        degree, polynomial, helpers = 4, 0b10011, 1
        width = 2 * degree + helpers
        basis = [1 << j for j in range(2 * width)]
        probes = set(basis) | {a ^ b for a, b in combinations(basis, 2)}
        probes.add((1 << (2 * width)) - 1)
        cubes = {}
        for origin in range(8):
            for a, b in combinations(range(1, 8), 2):
                order = (origin, origin ^ a, origin ^ b, origin ^ a ^ b)
                cubes.setdefault(frozenset(order), order)
        self.assertEqual(len(cubes), 14)
        for label in probes:
            orbit = [_field_image(label, x, degree, polynomial, helpers)
                     for x in range(8)]
            self.assertTrue(all(_pairing(a, b, width) == 0
                                for a, b in combinations(orbit, 2)))
        for label, cube in product(probes, cubes.values()):
            word_label, word_phase = 0, 0
            for address in cube:
                image = _field_image(label, address, degree, polynomial, helpers)
                sign = (_hermitian_phase(label, width)
                        - _hermitian_phase(image, width)) % 4
                word_phase = (word_phase + sign
                              + _product_phase(word_label, image, width)) % 4
                word_label ^= image
            self.assertEqual(word_label, 0)
            self.assertIn(word_phase, (0, 2))

    def test_noncommuting_cube_can_have_negative_scalar(self):
        # Q_x=H^(x0), P=X, on p=2 addresses. Its symplectic entries have
        # degree one and P_x=(X,Z,X,Z) in this fixed cube order.
        # The scalar product is -I, not +I; the proof must retain phases.
        label, phase = 0, 0
        for value in (1, 2, 1, 2):
            phase = (phase + _product_phase(label, value, 1)) % 4
            label ^= value
        self.assertEqual((label, phase), (0, 2))
        self.assertEqual(_pairing(1, 2, 1), 1)
        nonlinear_character = (0, 0, 0, 1)
        requested_phase = (phase + 2 * sum(nonlinear_character)) % 4
        self.assertEqual(requested_phase, 0)
        # The exact products differ by sign despite the negative source word.
        self.assertEqual((requested_phase - phase) % 4, 2)

    def test_derivative_characterization_and_symplectic_radical(self):
        # Exhaust every p=3 Boolean mask, rather than sampling polynomials.
        for degree in (0, 1, 2, 3):
            code = _span(_rm_basis(3, degree))
            permitted = set()
            directions = list(range(1, 8))
            for mask in range(256):
                good = True
                for ds in combinations(directions, degree + 1):
                    if _rank(ds) != degree + 1:
                        continue
                    cube = _span(ds)
                    for origin in range(8):
                        if sum((mask >> (origin ^ x)) & 1 for x in cube) % 2:
                            good = False
                            break
                    if not good:
                        break
                if good:
                    permitted.add(mask)
            self.assertEqual(permitted, code)
        for p, degree, helpers in ((3, 1, 0), (4, 1, 3), (5, 1, 3), (5, 2, 1)):
            size, width = 1 << p, (1 << p) + helpers
            code_basis = _rm_basis(p, degree)
            dimension = len(code_basis)
            self.assertEqual(_rank(code_basis), dimension)
            allowed = code_basis + [1 << (width + j) for j in range(size)]
            allowed += [1 << j for j in range(size, width)]
            allowed += [1 << (width + j) for j in range(size, width)]
            self.assertEqual(_rank(allowed), size + dimension + 2 * helpers)
            self.assertEqual(_gram_rank(allowed, width), 2 * dimension + 2 * helpers)
            self.assertEqual(len(allowed) - _gram_rank(allowed, width), size - dimension)

    def test_independent_cell_rank_survives_joint_clifford_walls(self):
        for p, helpers in product(range(2, 7), (0, 3)):
            size, cells = 1 << p, 1 << (p - 2)
            width = size + helpers
            generators = []
            for cell in range(cells):
                for label in (1, 4, 8, 18):
                    a, b = label & 15, label >> 4
                    generators.append((a << (4 * cell))
                                      | (b << (width + 4 * cell)))
            generators += [1 << j for j in range(size, width)]
            generators += [1 << (width + j) for j in range(size, width)]
            original_rank = _gram_rank(generators, width)
            self.assertEqual(original_rank, 2 * cells + 2 * helpers)
            wall = [('H', 0), ('S', 1), ('CX', 0, size - 1), ('CX', size - 1, 2)]
            if helpers:
                wall += [('H', width - 1), ('CX', 1, width - 1),
                         ('CX', width - 1, 3), ('S', width - 2)]
            changed = [_clifford_label(v, wall, width) for v in generators]
            self.assertEqual(_rank(changed), 4 * cells + 2 * helpers)
            self.assertEqual(_gram_rank(changed, width), original_rank)
            allowed_rank = 2 * (p + 1) + 2 * helpers
            self.assertEqual(original_rank > allowed_rank, p >= 5)
            if p >= 5:
                # The weak-probe T-count bound alone would not reject these
                # cells; their excessive domain symplectic rank does.
                self.assertGreaterEqual(6 * cells, size - (p + 1))

    def test_exact_field_z_kernel_can_only_contain_constant_affine_masks(self):
        for p in range(2, 7):
            size = 1 << p
            constant = {0, (1 << size) - 1}
            for mask in _span(_rm_basis(p, 1)):
                orbit = {_translate(mask, x, size) for x in range(size)}
                self.assertEqual(len(orbit), 1 if mask in constant else 2)
                self.assertTrue(any(_translate(mask, x, size) == mask
                                    for x in range(1, size)))
            # Tensor-cell domain dimension exceeds the exact bound from p=3,
            # while the four-address cell saturates its permitted dimension.
            self.assertEqual(size > p + 2, p >= 3)
        self.assertFalse(any(_translate(0b10, x, 2) == 0b10
                             for x in range(1, 2)))
        for p, degree, polynomial in ((2, 2, 0b111), (3, 4, 0b10011)):
            width, addresses = 2 * degree, 1 << p
            for label in [1 << j for j in range(2 * width)]:
                orbit = {_field_image(label, x, degree, polynomial)
                         for x in range(addresses)}
                self.assertIn(len(orbit), (1, addresses))
                # An affine mask's nonconstant two-point orbit therefore
                # cannot inject into a field Pauli orbit at p>=2.
                self.assertNotEqual(len(orbit), 2)

    def test_actual_p2_decoder_saturates_exact_dimension_and_z_kernel(self):
        word = full_intertwiner()
        image_domain = {}
        for label in range(256):
            row = _conjugate({label: Q2(1)}, word, 4)
            if len(row) == 1:
                target, coefficient = next(iter(row.items()))
                self.assertIn(coefficient, (Q2(1), Q2(-1)))
                image_domain[label] = target
        self.assertEqual(len(image_domain), 16)
        self.assertEqual(_rank(image_domain), 4)
        images = set(image_domain.values())
        self.assertEqual({v & 15 for v in images}, _span(_rm_basis(2, 1)))
        self.assertEqual({v for v in images if (v & 15) == 0}, {0, 240})
        # Pure target Z pullbacks have exactly the one-dimensional all-Z
        # kernel. Three independent affine X masks plus that kernel give
        # the saturated dimension p+2=4.
        self.assertEqual(_rank(images), 2 + 2)
        for helpers in (1, 3):
            width = 4 + helpers
            lifted = [(v & 15) | ((v >> 4) << width) for v in images]
            lifted += [1 << j for j in range(4, width)]
            lifted += [1 << (width + j) for j in range(4, width)]
            self.assertEqual(_rank(lifted), 2 + 2 + 2 * helpers)

    def test_radius_two_code_bound_and_weight_four_boundary(self):
        for p in range(2, 21):
            size = 1 << p
            sphere_size = 1 + size + size * (size - 1) // 2
            self.assertEqual(sphere_size - 2 * size,
                             (size - 1) * (size - 2) // 2)
            self.assertGreater(sphere_size, 2 * size)
            if p >= 3:
                self.assertLess(4, size)
        # The dual of RM(1,p) supplies the limiting code dimensions used
        # in the argument. Exhaust the small p=3,4 masks and find an actual
        # minimum of four; the packing step cannot generically promise less.
        for p in (3, 4):
            size = 1 << p
            checks = _rm_basis(p, 1)
            code = [mask for mask in range(1 << size)
                    if all((mask & row).bit_count() % 2 == 0 for row in checks)]
            self.assertEqual(len(code), 1 << (size - p - 1))
            self.assertEqual(min(mask.bit_count() for mask in code if mask), 4)


if __name__ == '__main__':
    unittest.main()
