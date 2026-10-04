"""Certified rotation coefficients for one residual SU(2) completion.

This small classical helper emits rational enclosures, not quantum gates.
The caller supplies a dyadic approximation a+ib to an unknown unit-disk z,
with certified Euclidean error at most 2**(-2*q-20). The three returned
(cosine, sine) pairs describe Rz, Ry, Rz in the Hopf convention. The two
outer pairs use the same square root of the phase, including its sign.

All decisions concern finite dyadic data. Square roots use integer isqrt;
there is no angle solver, search, floating arithmetic, or exact-zero oracle
for the original coefficient. See supplements/state_based_qbp/RESIDUAL_TABLE_PREPROCESSING.md.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from numbers import Integral


def _power_of_two(exponent: int) -> Fraction:
    return Fraction(1 << exponent) if exponent >= 0 else Fraction(1, 1 << -exponent)


@dataclass(frozen=True)
class CoefficientInterval:
    lower: Fraction
    upper: Fraction

    def __post_init__(self) -> None:
        if self.lower > self.upper:
            raise ValueError("An interval's lower endpoint must not exceed its upper endpoint.")

    @property
    def width(self) -> Fraction:
        return self.upper - self.lower


@dataclass(frozen=True)
class ResidualRotationData:
    interior: tuple[Fraction, Fraction]
    coefficient_error: Fraction
    small_radius: bool
    # Chronological Rz, Ry, Rz pairs; each pair is (cosine, sine).
    rotations: tuple[tuple[CoefficientInterval, CoefficientInterval], ...]
    working_bits: int
    target_bits: int
    completion_error_bound: Fraction


def _point(value: Fraction | int) -> CoefficientInterval:
    value = Fraction(value)
    return CoefficientInterval(value, value)


def _add(left: CoefficientInterval, right: CoefficientInterval) -> CoefficientInterval:
    return CoefficientInterval(left.lower + right.lower, left.upper + right.upper)


def _negate(value: CoefficientInterval) -> CoefficientInterval:
    return CoefficientInterval(-value.upper, -value.lower)


def _multiply(left: CoefficientInterval, right: CoefficientInterval) -> CoefficientInterval:
    products = tuple(a * b for a in (left.lower, left.upper)
                     for b in (right.lower, right.upper))
    return CoefficientInterval(min(products), max(products))


def _divide(left: CoefficientInterval, right: CoefficientInterval) -> CoefficientInterval:
    if right.lower <= 0 <= right.upper:
        raise ArithmeticError("The certified denominator interval contains zero.")
    inverse = CoefficientInterval(1 / right.upper, 1 / right.lower)
    return _multiply(left, inverse)


def _clip(value: CoefficientInterval, lower: Fraction | int, upper: Fraction | int) -> CoefficientInterval:
    return CoefficientInterval(max(value.lower, Fraction(lower)), min(value.upper, Fraction(upper)))


def sqrt_enclosure(value: Fraction, bits: int) -> CoefficientInterval:
    """Enclose a nonnegative rational square root on the 2**(-bits) grid."""
    if not isinstance(value, Fraction) or value < 0:
        raise ValueError("The radicand must be a nonnegative Fraction.")
    if not isinstance(bits, Integral) or bits < 0:
        raise ValueError("bits must be a nonnegative integer.")
    scale = 1 << int(bits)
    root = isqrt((value.numerator * scale * scale) // value.denominator)
    lower = Fraction(root, scale)
    upper = lower if lower * lower == value else Fraction(root + 1, scale)
    return CoefficientInterval(lower, upper)


def _sqrt_interval(value: CoefficientInterval, bits: int) -> CoefficientInterval:
    if value.lower < 0:
        raise ArithmeticError("The certified square-root interval is negative.")
    return CoefficientInterval(sqrt_enclosure(value.lower, bits).lower,
                               sqrt_enclosure(value.upper, bits).upper)


def _dyadic(value: Fraction | int) -> Fraction:
    if not isinstance(value, (Fraction, Integral)):
        raise TypeError("Coefficient approximations must be exact dyadic Fractions or integers.")
    value = Fraction(value)
    denominator = value.denominator
    if denominator & (denominator - 1):
        raise ValueError("Coefficient approximations must be dyadic.")
    return value


def residual_rotation_data(
    real: Fraction | int,
    imag: Fraction | int,
    q: int,
    *,
    target_bits: int | None = None,
) -> ResidualRotationData:
    """Return certified direct-programming coefficients for U(z).

    Input promise: |(real+i*imag)-z| <= delta = 2**(-2*q-20), |z|<=1.
    This promise is supplied by the caller's coefficient evaluator, not
    proved by this routine. A necessary unit-disk consistency check is
    performed. The ideal rotations enclosed here approximate U(z) within
    2**(-q-4); native rotation-synthesis error is a separate charge.
    """
    if not isinstance(q, Integral) or q < 5:
        raise ValueError("q must be an integer at least five.")
    q = int(q)
    if target_bits is None:
        target_bits = q + 20
    if not isinstance(target_bits, Integral) or target_bits < 1:
        raise ValueError("target_bits must be a positive integer.")
    target_bits = int(target_bits)
    a0, b0 = _dyadic(real), _dyadic(imag)
    delta = _power_of_two(-2 * q - 20)
    if a0 * a0 + b0 * b0 > (1 + delta) ** 2:
        raise ValueError("The approximation is inconsistent with the promised unit disk.")
    shrink = 1 - 2 * delta
    a, b = shrink * a0, shrink * b0
    radius_squared = a * a + b * b
    if not 0 <= radius_squared < 1:
        raise ArithmeticError("The certified interior shrink failed.")
    threshold = _power_of_two(-q - 6)
    working_bits = target_bits + 4 * q + 40
    if radius_squared <= threshold * threshold:
        outer = (_point(1), _point(0))
        middle = (_point(0), _point(1))
        return ResidualRotationData((a, b), delta, True, (outer, middle, outer),
                                    working_bits, target_bits, _power_of_two(-q - 4))

    radius = _clip(sqrt_enclosure(radius_squared, working_bits), 0, 1)
    sine = _clip(sqrt_enclosure(1 - radius_squared, working_bits), 0, 1)
    # The larger half-phase component is at least 1/sqrt(2). Its formula
    # avoids division by a small half-phase component at either branch cut.
    radicand = _add(_point(Fraction(1, 2)),
                     _divide(_point(abs(a)), _multiply(_point(2), radius)))
    radicand = _clip(radicand, Fraction(1, 2), 1)
    primary = _clip(_sqrt_interval(radicand, working_bits), Fraction(1, 2), 1)
    if a >= 0:
        u = primary
        v = _divide(_point(b), _multiply(_point(2), _multiply(radius, u)))
    else:
        # At finite dyadic b=0 the positive branch is selected. No question
        # about the sign or exact zero of the original unknown z is asked.
        v = primary if b >= 0 else _negate(primary)
        u = _divide(_point(b), _multiply(_point(2), _multiply(radius, v)))
    u, v = _clip(u, -1, 1), _clip(v, -1, 1)
    outer = (u, _negate(v))  # D(w)=diag(w,conj(w))=Rz with sin=-Im(w).
    middle = (radius, sine)
    rotations = outer, middle, outer
    if any(interval.width > _power_of_two(-target_bits)
           for pair in rotations for interval in pair):
        raise ArithmeticError("The fixed precision bound did not certify the requested width.")
    return ResidualRotationData((a, b), delta, False, rotations,
                                working_bits, target_bits, _power_of_two(-q - 4))
