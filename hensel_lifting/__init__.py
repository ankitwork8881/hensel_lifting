"""Hensel lifting of polynomial roots and modular square roots."""

from .arith import crt, egcd, factorize, inverse, is_prime
from .poly import Polynomial, lift_roots, lift_simple_root, roots_mod_prime
from .sqrt import legendre, sqrt_mod, sqrt_mod_prime_power, tonelli_shanks

__all__ = [
    "Polynomial",
    "crt",
    "egcd",
    "factorize",
    "inverse",
    "is_prime",
    "legendre",
    "lift_roots",
    "lift_simple_root",
    "roots_mod_prime",
    "sqrt_mod",
    "sqrt_mod_prime_power",
    "tonelli_shanks",
]

__version__ = "0.1.0"
