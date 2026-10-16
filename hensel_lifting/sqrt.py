"""Square roots modulo primes, prime powers, and arbitrary moduli.

Every function returns *all* solutions as a sorted list (empty when a is
not a square).
"""

from .arith import crt, factorize, inverse, is_prime, valuation


def legendre(a, p):
    """Legendre symbol (a / p) for an odd prime p, as -1, 0 or 1."""
    s = pow(a, (p - 1) // 2, p)
    return -1 if s == p - 1 else s


def tonelli_shanks(a, p):
    """One square root of a modulo an odd prime p, or None."""
    a %= p
    if a == 0:
        return 0
    if legendre(a, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)

    # Write p - 1 = q * 2^s with q odd.
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = next(z for z in range(2, p) if legendre(z, p) == -1)

    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        # Least i with t^(2^i) = 1.
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p
    return r


def _unit_roots_odd(u, p, k):
    """Roots of x^2 = u (mod p^k), p odd, gcd(u, p) = 1."""
    r = tonelli_shanks(u, p)
    if r is None:
        return []
    m = p ** k
    precision = 1
    while precision < k:
        precision = min(2 * precision, k)
        pm = p ** precision
        r = (r - (r * r - u) * inverse(2 * r, pm)) % pm
    return sorted({r % m, -r % m})


def _unit_roots_two(u, k):
    """Roots of x^2 = u (mod 2^k), u odd."""
    m = 1 << k
    if k == 1:
        return [1]
    if k == 2:
        return [1, 3] if u % 4 == 1 else []
    if u % 8 != 1:
        return []
    # Fix one bit per step: if r^2 != u mod 2^(j+1), flip bit j-1.
    r = 1
    for j in range(3, k):
        if (r * r - u) % (1 << (j + 1)):
            r += 1 << (j - 1)
    half = m >> 1
    return sorted({r % m, -r % m, (r + half) % m, (-r + half) % m})


def sqrt_mod_prime_power(a, p, k):
    """All x in [0, p^k) with x^2 = a (mod p^k)."""
    if k < 1 or not is_prime(p):
        raise ValueError("need a prime p and k >= 1")
    m = p ** k
    a %= m

    if a == 0:
        # x^2 = 0 mod p^k  <=>  p^ceil(k/2) | x.
        step = p ** ((k + 1) // 2)
        return list(range(0, m, step))

    e = valuation(a, p)
    if e % 2:
        return []
    # a = p^(2h) u with u a unit; x = p^h y with y^2 = u (mod p^(k-2h)),
    # and y is only determined modulo p^(k-2h).
    h = e // 2
    u = a // p ** e
    reduced = k - 2 * h
    base = _unit_roots_two(u, reduced) if p == 2 else _unit_roots_odd(u, p, reduced)
    lifted_modulus = p ** reduced
    scale = p ** h
    return sorted(
        scale * (y + t * lifted_modulus) % m
        for y in base
        for t in range(p ** h)
    )


def sqrt_mod(a, n):
    """All square roots of a modulo any n >= 1, combined with the CRT."""
    if n < 1:
        raise ValueError("modulus must be positive")
    if n == 1:
        return [0]
    parts = [
        (p ** k, sqrt_mod_prime_power(a, p, k)) for p, k in factorize(n).items()
    ]
    solutions = [0]
    modulus = 1
    for m, roots in parts:
        if not roots:
            return []
        solutions = [crt([s, r], [modulus, m]) for s in solutions for r in roots]
        modulus *= m
    return sorted(solutions)
