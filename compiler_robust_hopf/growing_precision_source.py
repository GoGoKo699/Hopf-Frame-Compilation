"""Exact refinement identities and symbolic growing-source accounting.

These helpers describe finite algebra and a scale schedule. They do not emit
native circuits, certify clean arena capacity, or supply precision states.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral

from .heterogeneous_precision_source import dyadic_precision


def _integer(value: int, name: str, minimum: int = 0) -> int:
    if not isinstance(value, Integral) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return int(value)


def _fraction(value: Fraction | int, name: str) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise ValueError(f"{name} must be an exact rational")
    if isinstance(value, Fraction):
        return Fraction(int(value.numerator), int(value.denominator))
    return Fraction(int(value))


@dataclass(frozen=True)
class RefinementParameters:
    """Rare-LCU scalar and normalized geometric-source survival probability."""

    beta: Fraction
    rare_probability: Fraction
    scalar_coefficient: Fraction
    survival_squared: Fraction


def refinement_parameters(beta: Fraction, rare_probability: Fraction) -> RefinementParameters:
    """Specify T=(1-p)I+p*f*X=(1-p)(I+beta*X), with |f|<=1.

    On the embedded normalized source g_(beta**2,M) tensor |0>, T has
    output sqrt(survival_squared)*g_(beta,2M). The rare register and raw
    scalar rejection flag remain distinct parts of the full unitary word.
    """
    beta = _fraction(beta, "beta")
    rare_probability = _fraction(rare_probability, "rare_probability")
    dyadic_precision(beta)
    dyadic_precision(rare_probability)
    if not 0 < beta < 1 or not 0 < rare_probability < 1 or rare_probability.numerator != 1:
        raise ValueError("require 0<beta<1 and reciprocal-power-of-two probability")
    coefficient = (1 - rare_probability) * beta / rare_probability
    if coefficient > 1:
        raise ValueError("the raw scalar coefficient must be a contraction")
    survival_squared = (1 - rare_probability) ** 2 * (1 + beta ** 2)
    return RefinementParameters(beta, rare_probability, coefficient, survival_squared)


@dataclass(frozen=True)
class GrowingSourceScale:
    """Absolute-height schedule, storing log2(q) instead of q or 2**N."""

    endpoint_height: int
    absolute_height: int
    order: int
    source_length: int
    exponent_log2: int
    rare_probability: Fraction
    rescaling_bits: int

    @property
    def integral_exponent(self) -> bool:
        """Small finite cases may require the theorem's separate fallback."""
        return self.exponent_log2 >= 0

    def refinement_steps_from(self, previous: GrowingSourceScale) -> int:
        """Number of q halvings, each doubling the source length."""
        if (self.endpoint_height != previous.endpoint_height or
                self.absolute_height < previous.absolute_height):
            raise ValueError("scales must advance within one endpoint")
        return previous.exponent_log2 - self.exponent_log2


def growing_source_scale(endpoint_height: int, absolute_height: int) -> GrowingSourceScale:
    """Set P_d=2**(2 ceil(log2(d+1))), M_d=4P_d**2, q_d=N/P_d**2.

    The common source has beta_d=2**(-q_d); coefficient spacing is q_d+B_d.
    No logical padding or sufficient arena-allocation constant is inferred.
    """
    endpoint_height = _integer(endpoint_height, "endpoint_height", 1)
    absolute_height = _integer(absolute_height, "absolute_height")
    if absolute_height > endpoint_height:
        raise ValueError("absolute_height cannot exceed endpoint_height")
    logarithm = absolute_height.bit_length()
    order = 1 << (2 * logarithm)
    return GrowingSourceScale(endpoint_height, absolute_height, order, 4 * order ** 2,
                              endpoint_height - 4 * logarithm, Fraction(1, 8 * order),
                              6 + 4 * logarithm)
