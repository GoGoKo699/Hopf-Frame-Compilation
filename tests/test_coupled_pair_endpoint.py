"""Exact six-wire endpoint obstruction, not an unrestricted depth bound.

Sparse Pauli arithmetic retains literal signed axes and T phases in
Q(sqrt(2), i).  Operator identities cover the full physical space; the
analytic coefficient/projector argument supplies the return-error floor.
"""
from fractions import Fraction
import itertools
import unittest

from tests.exact_matrices import CI, CO, CQ2, CZ, Q2
from tests.test_clean_code_clifford_rigidity import _rank


IDENTITY = 'IIIIII'
HALF = Fraction(1, 2)
ROOT_HALF = CQ2(Q2(0, HALF))
OMEGA = ROOT_HALF * (CO + CI)


def _add(*operators):
    result = {}
    for operator in operators:
        for axis, value in operator.items():
            result[axis] = result.get(axis, CZ) + value
    return {axis: value for axis, value in result.items() if value != CZ}


def _scale(operator, value):
    return {axis: coefficient * value for axis, coefficient in operator.items()
            if coefficient * value != CZ}


def _pauli_product(left, right):
    phase, result = CO, []
    cyclic = {'XY': 'Z', 'YZ': 'X', 'ZX': 'Y'}
    for a, b in zip(left, right):
        if a == 'I':
            result.append(b)
        elif b == 'I':
            result.append(a)
        elif a == b:
            result.append('I')
        elif a + b in cyclic:
            phase *= CI
            result.append(cyclic[a + b])
        else:
            phase *= -CI
            result.append(cyclic[b + a])
    return ''.join(result), phase


def _product(*operators):
    result = {IDENTITY: CO}
    for operator in operators:
        terms = []
        for a, ca in result.items():
            for b, cb in operator.items():
                axis, phase = _pauli_product(a, b)
                terms.append({axis: ca * cb * phase})
        result = _add(*terms)
    return result


def _adjoint(operator):
    return {axis: value.conj() for axis, value in operator.items()}


def _conjugate(unitary, operator):
    return _product(unitary, operator, _adjoint(unitary))


def _axis(axis, sign=1):
    return {axis: CO * sign}


def _t(axis, sign=1):
    return {IDENTITY: (CO + OMEGA) / 2, axis: sign * (CO - OMEGA) / 2}


def _wire(bit, axis):
    return 'I' * bit + axis + 'I' * (5 - bit)


def _gate(gate):
    kind, *bits = gate
    if kind == 'H':
        return {_wire(bits[0], 'X'): ROOT_HALF, _wire(bits[0], 'Z'): ROOT_HALF}
    if kind in ('S', 'SDG'):
        phase = CI if kind == 'S' else -CI
        return {IDENTITY: (CO + phase) / 2,
                _wire(bits[0], 'Z'): (CO - phase) / 2}
    if kind == 'CX':
        z, x = _wire(bits[0], 'Z'), _wire(bits[1], 'X')
        zx, phase = _pauli_product(z, x)
        assert phase == CO
        return {IDENTITY: CO / 2, z: CO / 2, x: CO / 2, zx: -CO / 2}
    raise ValueError(kind)


def _word_image(word, operator):
    for gate in word:
        operator = _conjugate(_gate(gate), operator)
    return operator


def _inverse(word):
    names = {'H': 'H', 'CX': 'CX', 'S': 'SDG', 'SDG': 'S'}
    return [(names[kind], *bits) for kind, *bits in reversed(word)]


def _clifford_count(word):
    return sum(3 if gate[0] == 'SDG' else 1 for gate in word)


def _label(axis):
    x = sum((letter in 'XY') << bit for bit, letter in enumerate(axis))
    z = sum((letter in 'ZY') << bit for bit, letter in enumerate(axis))
    return x | (z << 6)


B_AXES = [('XXIIII', 1), ('YYIIII', -1), ('IIZIIX', 1),
          ('IIIZIX', 1), ('IIIIZX', 1), ('IIXXXY', 1)]
D_AXES = [('ZIIIII', 1), ('IZIIII', 1), ('ZIIIIZ', 1), ('ZIYYYI', 1)]
J_CORE = _scale(_add(*[_axis(axis) for axis in
                      ('ZIIIII', 'IZIIII', 'XXIIII', 'YYIIII')]), HALF)
J_W = _scale(_add(_axis('IIIIIX'), _axis('IIIIIY')), ROOT_HALF)
J = _product(J_CORE, J_W)


def _middle_image(operator, inverse=False, b_exponents=None, d_exponents=None):
    b_exponents = (1,) * len(B_AXES) if b_exponents is None else b_exponents
    d_exponents = (1,) * len(D_AXES) if d_exponents is None else d_exponents
    assert len(b_exponents) == len(B_AXES) and len(d_exponents) == len(D_AXES)
    exponents = tuple(b_exponents) + tuple(d_exponents)
    assert set(exponents) <= {-1, 1}
    # T-dagger is the actual adjoint; replacing the signed axis would
    # introduce a literal scalar, even though conjugation cancels it.
    gates = [(_t(axis, sign) if exponent == 1 else _adjoint(_t(axis, sign)))
             for (axis, sign), exponent in zip(B_AXES + D_AXES, exponents)]
    if inverse:
        gates = [_adjoint(gate) for gate in reversed(gates)]
    for gate in gates:
        operator = _conjugate(gate, operator)
    return operator


class CoupledPairEndpointTests(unittest.TestCase):
    def test_actual_signed_native_lists_and_paid_encoders(self):
        b_encoder = [('CX', 0, 1), ('H', 0), ('CX', 0, 1)]
        b_encoder += [('H', bit) for bit in (2, 3, 4, 5)]
        for bit in (2, 3, 4):
            b_encoder += [('H', 5), ('CX', bit, 5), ('H', 5)]
        b_encoder += [('H', bit) for bit in (2, 3, 4, 5)]
        b_encoder += [('H', 5), ('SDG', 5), ('H', 5)]
        d_encoder = [('CX', 0, 5), ('CX', 0, 2), ('CX', 3, 2), ('CX', 4, 2)]
        for bit in (2, 3, 4):
            d_encoder += [('H', bit), ('S', bit)]
        for word, source_bits, axes in ((b_encoder, range(6), B_AXES),
                                        (d_encoder, (0, 1, 5, 2), D_AXES)):
            for bit, (axis, sign) in zip(source_bits, axes):
                actual = _word_image(word, _axis(_wire(bit, 'Z')))
                self.assertEqual(actual, _axis(axis, sign))
                self.assertEqual(_word_image(_inverse(word), actual),
                                 _axis(_wire(bit, 'Z')))
            self.assertEqual(_rank([_label(axis) for axis, _ in axes]), len(axes))
            for (a, _), (b, _) in itertools.combinations(axes, 2):
                self.assertEqual(_product(_axis(a), _axis(b)),
                                 _product(_axis(b), _axis(a)))
        self.assertEqual([_clifford_count(word) for word in
                          (b_encoder, _inverse(b_encoder),
                           d_encoder, _inverse(d_encoder))], [25, 23, 10, 16])
        b_labels = [_label(axis) for axis, _ in B_AXES]
        d_labels = [_label(axis) for axis, _ in D_AXES]
        self.assertEqual(_rank(b_labels + d_labels), 8)
        for center in ('ZZIIII', 'IIYYYZ'):
            self.assertEqual(_rank(b_labels + [_label(center)]), 6)
            self.assertEqual(_rank(d_labels + [_label(center)]), 4)
        self.assertEqual(_conjugate(J, _axis('XXIIII')), _axis('YYIIII'))
        self.assertEqual(_product(*[_axis(a, s) for a, s in B_AXES[2:]]),
                         _axis('IIYYYZ'))

    def test_full_probe_and_inverse_have_the_literal_four_term_witness(self):
        probe = _middle_image(J)
        self.assertEqual(_adjoint(probe), probe)
        self.assertEqual(_product(probe, probe), _axis(IDENTITY))
        self.assertEqual(_middle_image(probe, inverse=True), J)
        expected = {'YIIIIZ': CO / 2, 'YIYYYI': CO / 2,
                    'YIXXXX': -CO / 2, 'YIZZZY': -CO / 2}
        self.assertEqual(_conjugate(probe, _axis('XIIIII')), expected)
        # Orthogonal Pauli coefficients give squared normalized-HS distance
        # at least 2 - 2*(1/2) = 1 from every signed Pauli.
        self.assertEqual(sum((value.real * value.real + value.imag * value.imag
                              for value in expected.values()), Q2()), Q2(1))
        self.assertEqual(max(abs(value.real) for value in expected.values()),
                         Q2(HALF))

    def test_sector_calculation_keeps_the_flipped_character_and_pair_phases(self):
        p, c, b, h = map(_axis, ('XXIIII', 'ZIIIII', 'IIIIIZ', 'IIYYYZ'))
        f = _product(_t('ZIIIIZ'), _t('ZIYYYI'))
        e = _product(*[_t(a, s) for a, s in B_AXES[2:]])
        w = _conjugate(e, J_W)
        plus = _scale(_add(_axis(IDENTITY), h), HALF)
        minus = _scale(_add(_axis(IDENTITY), _scale(h, -1)), HALF)
        self.assertEqual(_middle_image(J), _conjugate(f, _product(p, w)))
        self.assertEqual(_conjugate(_axis('XIIIII'), f), _scale(_adjoint(f), CI))
        self.assertEqual(_product(f, f),
                         _add(_product(plus, c, b), _scale(minus, CI)))
        self.assertEqual(_conjugate(w, h), _scale(h, -1))
        self.assertEqual(_product(minus, _conjugate(w, b)),
                         _scale(_product(minus, _axis('IIZZZY')), -1))

    def test_all_core_orientations_and_mixed_tail_clifford_reduction(self):
        positive = {'YIIIIZ': CO / 2, 'YIYYYI': CO / 2,
                    'YIXXXX': -CO / 2, 'YIZZZY': -CO / 2}
        opposite = {'ZXIIII': CO / 2, 'ZXYYYZ': CO / 2,
                    'ZYXXXX': CO / 2, 'ZYZZZY': CO / 2}
        for core in itertools.product((1, -1), repeat=2):
            exponents = core + (1,) * 4
            probe = _middle_image(J, b_exponents=exponents)
            expected = (positive if core == (1, 1) else
                        _scale(positive, -1) if core == (-1, -1) else opposite)
            self.assertEqual(_conjugate(probe, _axis('XIIIII')), expected)
            self.assertEqual(_middle_image(probe, inverse=True,
                                           b_exponents=exponents), J)

        # Three flipped tails are induced by Y_z1 X_z3. Their odd
        # parity also flips dh in the conjugated D list. The reduced
        # D orientations are therefore (-,+,-,-), all paid Cliffords
        # relative to the positive D band.
        b_exponents, d_exponents = (-1, 1, -1, 1, -1, -1), (-1, 1, -1, 1)
        actual = _middle_image(J, b_exponents=b_exponents, d_exponents=d_exponents)
        reduced = _middle_image(J, b_exponents=(-1, 1, 1, 1, 1, 1))
        witness = _axis('XIIIII')
        for index in (0, 2, 3):
            axis, sign = D_AXES[index]
            s_dagger = _adjoint(_product(_t(axis, sign), _t(axis, sign)))
            reduced = _conjugate(s_dagger, reduced)
            witness = _conjugate(s_dagger, witness)
        tail_pauli = _axis('IIYIXI')
        reduced = _conjugate(tail_pauli, reduced)
        witness = _conjugate(tail_pauli, witness)
        self.assertEqual(actual, reduced)
        self.assertEqual(len(witness), 1)
        image = _conjugate(actual, witness)
        self.assertEqual(len(image), 4)
        self.assertTrue(all(value.imag == Q2() and abs(value.real) == Q2(HALF)
                            for value in image.values()))
        self.assertEqual(_middle_image(actual, inverse=True,
                                       b_exponents=b_exponents,
                                       d_exponents=d_exponents), J)


if __name__ == '__main__':
    unittest.main()
