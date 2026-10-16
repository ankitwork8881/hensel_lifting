"""Hensel lifting of polynomial roots modulo prime powers.

If f(r) = 0 (mod p) and f'(r) != 0 (mod p), Hensel's lemma gives a unique
root of f modulo every p^k congruent to r modulo p, computed by the Newton
iteration r <- r - f(r) / f'(r). Each iteration doubles the number of
correct p-adic digits.

When f'(r) = 0 (mod p) the root is *singular*: it may lift to many roots
or to none, and lifting has to branch digit by digit.
"""

from .arith import inverse


class Polynomial:
    """Integer polynomial with coefficients in increasing degree order."""

    def __init__(self, coefficients):
        coefficients = [int(c) for c in coefficients]
        while len(coefficients) > 1 and coefficients[-1] == 0:
            coefficients.pop()
        self.coefficients = coefficients or [0]

    @property
    def degree(self):
        return len(self.coefficients) - 1

    def __call__(self, x, modulus=None):
        value = 0
        for c in reversed(self.coefficients):
            value = value * x + c
            if modulus is not None:
                value %= modulus
        return value

    def derivative(self):
        return Polynomial([i * c for i, c in enumerate(self.coefficients)][1:] or [0])

    def __repr__(self):
        return f"Polynomial({self.coefficients})"


def roots_mod_prime(f, p):
    """All roots of f modulo the prime p, by evaluation."""
    return [r for r in range(p) if f(r, p) == 0]


def lift_simple_root(f, root, p, k):
    """Lift a simple root of f modulo p to the unique root modulo p^k.

    Uses quadratic (Newton) lifting: the precision goes p -> p^2 -> p^4 ...
    """
    df = f.derivative()
    if f(root, p) != 0:
        raise ValueError(f"{root} is not a root of {f} modulo {p}")
    if df(root, p) == 0:
        raise ValueError(f"{root} is a singular root modulo {p}; use lift_roots")
    target = p ** k
    r, precision = root % p, 1
    while precision < k:
        precision = min(2 * precision, k)
        m = p ** precision
        r = (r - f(r, m) * inverse(df(r, m), m)) % m
    return r % target


def lift_roots(f, p, k, roots=None):
    """All roots of f modulo p^k, lifting one p-adic digit at a time.

    Handles singular roots by trying every digit. ``roots`` may supply the
    roots modulo p; otherwise they are found by evaluation.
    """
    if k < 1:
        raise ValueError("k must be at least 1")
    current = sorted(roots_mod_prime(f, p) if roots is None else set(roots))
    modulus = p
    for _ in range(1, k):
        next_modulus = modulus * p
        current = [
            r + t * modulus
            for r in current
            for t in range(p)
            if f(r + t * modulus, next_modulus) == 0
        ]
        modulus = next_modulus
    return sorted(current)
