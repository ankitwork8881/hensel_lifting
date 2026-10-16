"""Elementary modular arithmetic used by the lifting routines."""


def egcd(a, b):
    """Return (g, s, t) with g = gcd(a, b) = s*a + t*b."""
    s0, s1, t0, t1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        s0, s1 = s1, s0 - q * s1
        t0, t1 = t1, t0 - q * t1
    return a, s0, t0


def inverse(a, m):
    """Inverse of a modulo m; raises ValueError if it does not exist."""
    g, s, _ = egcd(a % m, m)
    if g != 1:
        raise ValueError(f"{a} is not invertible modulo {m}")
    return s % m


def is_prime(n):
    """Deterministic primality test for n < 3.3e24 (Miller-Rabin)."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for p in small:
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factorize(n):
    """Prime factorisation of n >= 1 by trial division, as {p: k}."""
    factors = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def valuation(a, p):
    """Largest e with p^e | a (a != 0)."""
    e = 0
    while a % p == 0:
        a //= p
        e += 1
    return e


def crt(residues, moduli):
    """Combine x = r_i (mod m_i) for pairwise coprime m_i into x mod prod(m_i)."""
    x, m = 0, 1
    for r, mi in zip(residues, moduli):
        # Solve x + m*t = r (mod mi) for t.
        t = (r - x) * inverse(m, mi) % mi
        x += m * t
        m *= mi
    return x % m
