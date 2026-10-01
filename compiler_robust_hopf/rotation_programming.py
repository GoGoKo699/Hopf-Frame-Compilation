"""Certified paired-source programming from real rotation coefficients.

The caller certifies that the supplied intervals contain the cosine and
sine of one common angle. This module checks rational types, widths and
unit-circle consistency, then emits Majorana sign bits and an exact
accepted-block error certificate. It neither finds an angle nor emits
quantum gates. See ONE_CLEAN_COMPILER.md section 3 and
RESIDUAL_TABLE_PREPROCESSING.md sections 3--4.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .residual_table_preprocessing import CoefficientInterval, sqrt_enclosure


@dataclass(frozen=True)
class RotationProgram:
    """One paired-source program; ``signs`` contains flip bits, not +/-1.

    Bits have Majorana order: head 0, odd tail 1,3,...,2q-1, even tail
    2,4,...,2q, and unused final bit 2q+1 fixed to zero. The actual
    compressed operator is ``s I + p XZ``. Cosine and sine record the
    validated input intervals intersected with [-1,1]. Error bounds are
    operator norms against their promised common unit-circle target.
    """

    q: int
    signs: tuple[int, ...]
    cosine: CoefficientInterval
    sine: CoefficientInterval
    beta: CoefficientInterval
    s: Fraction
    p: Fraction
    block_error_bound: Fraction
    operator_error_bound: Fraction


def _point(value: Fraction | int) -> CoefficientInterval:
    value = Fraction(value)
    return CoefficientInterval(value, value)


def _scale(interval: CoefficientInterval, value: Fraction | int) -> CoefficientInterval:
    values = interval.lower * value, interval.upper * value
    return CoefficientInterval(min(values), max(values))


def _add(left: CoefficientInterval, right: CoefficientInterval) -> CoefficientInterval:
    return CoefficientInterval(left.lower + right.lower, left.upper + right.upper)


def _multiply(left: CoefficientInterval, right: CoefficientInterval) -> CoefficientInterval:
    values = tuple(a * b for a in (left.lower, left.upper)
                   for b in (right.lower, right.upper))
    return CoefficientInterval(min(values), max(values))


def _midpoint(interval: CoefficientInterval) -> Fraction:
    return (interval.lower + interval.upper) / 2


def _validated_coefficient(value: CoefficientInterval, q: int) -> CoefficientInterval:
    if not isinstance(value, CoefficientInterval):
        raise TypeError("Coefficients must be CoefficientInterval instances.")
    endpoints = value.lower, value.upper
    if any(isinstance(x, bool) or not isinstance(x, (Fraction, Integral))
           for x in endpoints):
        raise TypeError("Interval endpoints must be exact rational numbers, not floats.")
    lower, upper = map(Fraction, endpoints)
    if lower > upper:
        raise ValueError("The coefficient interval is reversed.")
    if upper - lower > Fraction(1, 1 << (q + 20)):
        raise ValueError("Coefficient interval width exceeds 2**(-q-20).")
    lower, upper = max(lower, Fraction(-1)), min(upper, Fraction(1))
    if lower > upper:
        raise ValueError("The coefficient interval does not intersect [-1,1].")
    return CoefficientInterval(lower, upper)


def _squared_range(interval: CoefficientInterval) -> CoefficientInterval:
    endpoints = interval.lower ** 2, interval.upper ** 2
    lower = Fraction(0) if interval.lower <= 0 <= interval.upper else min(endpoints)
    return CoefficientInterval(lower, max(endpoints))


def _encode_tail_mean(estimate: Fraction, q: int) -> tuple[tuple[int, ...], Fraction]:
    """Clip, then round on the geometric-tail grid, ties toward larger k.

    This finite rational tie convention rounds the mean downward at an
    exact midpoint. The last repeated geometric weight is left positive
    except at the literal -1 endpoint, where all bits are one.
    """
    estimate = max(Fraction(-1), min(Fraction(1), estimate))
    index = (1 - estimate) * (1 << (q - 2))
    integer = (2 * index.numerator + index.denominator) // (2 * index.denominator)
    endpoint = 1 << (q - 1)
    if integer == endpoint:
        bits = (1,) * q
    else:
        bits = tuple((integer >> shift) & 1 for shift in range(q - 2, -1, -1)) + (0,)
    mean = 1 - Fraction(integer, 1 << (q - 2))
    return bits, mean


def program_rotation(
    cosine: CoefficientInterval,
    sine: CoefficientInterval,
    q: int,
) -> RotationProgram:
    """Program a real rotation with entirely rational error certificates.

    Require q>=5 and input widths at most 2**(-q-20). Endpoints may be
    Fractions or integers. The caller promises enclosure of the intended
    common unit-circle pair; consistency checks do not certify that
    external evaluation promise. The returned native bound refers to the
    retained five-call rotation word with an arbitrary borrowed signal.
    """
    if isinstance(q, bool) or not isinstance(q, Integral) or q < 5:
        raise ValueError("q must be an integer at least five.")
    q = int(q)
    cosine = _validated_coefficient(cosine, q)
    sine = _validated_coefficient(sine, q)
    radius_squared = _add(_squared_range(cosine), _squared_range(sine))
    if not radius_squared.lower <= 1 <= radius_squared.upper:
        raise ValueError("The coefficient intervals are inconsistent with the unit circle.")

    root_five = sqrt_enclosure(Fraction(5), q + 24)
    beta = _scale(_add(root_five, _point(-1)), Fraction(1, 4))
    desired_s = _multiply(beta, cosine)
    desired_p = _multiply(beta, sine)
    plus = _add(_scale(desired_s, Fraction(3, 4)), _scale(desired_p, Fraction(-1, 2)))
    minus = _add(_scale(desired_s, Fraction(1, 4)), _scale(desired_p, Fraction(1, 2)))
    sigma = 1 if _midpoint(plus) >= 0 else -1
    odd_mean = _add(_scale(plus, 2), _point(Fraction(-sigma, 2)))
    even_mean = _scale(minus, 4)

    # An interval of width w=2**(-q-20), beta in [0,1/3], gives
    # p+ width <=5w/3 and p- width <=w. Thus these midpoint errors
    # are far below e/4, e=2**(1-q). A mistaken head sign only occurs
    # within 1/16 of zero, retaining the documented |u|<=5/8 bound.
    unit = Fraction(1, 1 << q)
    if plus.width / 2 > Fraction(1, 16):
        raise ArithmeticError("The head-sign estimate lacks its error certificate.")
    if max(odd_mean.width, even_mean.width) / 2 > unit / 2:
        raise ArithmeticError("The tail-mean estimate lacks its error certificate.")
    odd_bits, odd_value = _encode_tail_mean(_midpoint(odd_mean), q)
    even_bits, even_value = _encode_tail_mean(_midpoint(even_mean), q)
    signs = [0] * (2 * (q + 1))
    signs[0] = int(sigma == -1)
    signs[1:2 * q:2] = odd_bits
    signs[2:2 * q + 1:2] = even_bits

    actual_plus = Fraction(sigma, 4) + odd_value / 2
    actual_minus = even_value / 4
    s = actual_plus + actual_minus
    p = -actual_plus / 2 + 3 * actual_minus / 2
    error_s = max(abs(s - desired_s.lower), abs(s - desired_s.upper))
    error_p = max(abs(p - desired_p.lower), abs(p - desired_p.upper))
    block_error = sqrt_enclosure(error_s ** 2 + error_p ** 2, q + 24).upper
    if block_error >= Fraction(5, 2) * unit:
        raise ArithmeticError("The accepted-block error exceeds the proved digit budget.")
    if block_error > beta.lower / 2 or beta.upper + block_error > Fraction(1, 2):
        raise ArithmeticError("The accepted-block error exceeds the amplification regime.")

    # ONE_CLEAN sections 4/9: initialized error <12*block_error,
    # full borrowed-signal error gains sqrt(2). Consequently it is
    # <30*sqrt(2)*2**(-q)<43*2**(-q), since 1800<43**2.
    return RotationProgram(q, tuple(signs), cosine, sine, beta, s, p,
                           block_error, 43 * unit)
