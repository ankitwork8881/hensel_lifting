import math
import random

import pytest

from hensel_lifting import (
    Polynomial,
    crt,
    egcd,
    factorize,
    inverse,
    is_prime,
    legendre,
    lift_roots,
    lift_simple_root,
    sqrt_mod,
    sqrt_mod_prime_power,
    tonelli_shanks,
)

SMALL_PRIMES = [2, 3, 5, 7, 11, 13]


def brute_sqrt(a, m):
    return [x for x in range(m) if (x * x - a) % m == 0]


# --- arithmetic -------------------------------------------------------------

def test_egcd_and_inverse():
    rng = random.Random(0)
    for _ in range(200):
        a, b = rng.randrange(1, 10 ** 6), rng.randrange(1, 10 ** 6)
        g, s, t = egcd(a, b)
        assert g == math.gcd(a, b) == s * a + t * b
        if g == 1:
            assert a * inverse(a, b) % b == 1
    with pytest.raises(ValueError):
        inverse(6, 9)


def test_is_prime_and_factorize():
    sieve = [n for n in range(2, 2000) if all(n % d for d in range(2, int(n ** 0.5) + 1))]
    assert [n for n in range(2000) if is_prime(n)] == sieve
    assert is_prime(2 ** 61 - 1) and not is_prime(2 ** 61 + 1)
    for n in range(1, 3000):
        assert math.prod(p ** k for p, k in factorize(n).items()) == n


def test_crt():
    moduli = [4, 9, 25, 7]
    for x in range(0, 6300, 37):
        assert crt([x % m for m in moduli], moduli) == x


# --- polynomial lifting ------------------------------------------------------

def brute_poly_roots(f, m):
    return [x for x in range(m) if f(x, m) == 0]


@pytest.mark.parametrize("p", [3, 5, 7])
@pytest.mark.parametrize("k", [1, 2, 3, 4])
def test_lift_roots_matches_brute_force(p, k):
    rng = random.Random(p * 10 + k)
    for _ in range(10):
        f = Polynomial([rng.randrange(-20, 20) for _ in range(rng.randrange(2, 5))])
        assert lift_roots(f, p, k) == brute_poly_roots(f, p ** k)


def test_lift_simple_root_is_unique_lift():
    f = Polynomial([-2, 0, 1])  # x^2 - 2, simple roots 3, 4 mod 7
    for r in (3, 4):
        lifted = lift_simple_root(f, r, 7, 20)
        assert f(lifted, 7 ** 20) == 0 and lifted % 7 == r


def test_lift_simple_root_rejects_bad_input():
    f = Polynomial([0, 0, 1])  # x^2: root 0 mod p is singular
    with pytest.raises(ValueError):
        lift_simple_root(f, 0, 5, 3)
    with pytest.raises(ValueError):
        lift_simple_root(Polynomial([-2, 0, 1]), 1, 7, 3)


def test_singular_root_can_branch_or_vanish():
    # x^2 + 1 mod 2: root 1 is singular; mod 4 it has no lifts.
    assert lift_roots(Polynomial([1, 0, 1]), 2, 2) == []
    # x^2 - 9 mod 3^k: root 0 is singular and branches.
    assert lift_roots(Polynomial([-9, 0, 1]), 3, 3) == brute_poly_roots(
        Polynomial([-9, 0, 1]), 27
    )


def test_polynomial_basics():
    f = Polynomial([1, 2, 3, 0, 0])
    assert f.degree == 2
    assert f(2) == 17 and f(2, 5) == 2
    assert f.derivative().coefficients == [2, 6]


# --- square roots ------------------------------------------------------------

def test_tonelli_shanks_all_residues():
    for p in [p for p in range(3, 400) if is_prime(p)]:
        for a in range(p):
            r = tonelli_shanks(a, p)
            if legendre(a, p) == -1:
                assert r is None
            else:
                assert r * r % p == a


@pytest.mark.parametrize("p", SMALL_PRIMES)
def test_sqrt_mod_prime_power_exhaustive(p):
    k = 1
    while p ** k <= 1500:
        m = p ** k
        for a in range(m):
            assert sqrt_mod_prime_power(a, p, k) == brute_sqrt(a, m), (a, p, k)
        k += 1


def test_sqrt_mod_large_prime_power():
    p, k = 1_000_003, 6
    m = p ** k
    for x in (2, 5, 123_456, 7 * p):
        target = x * x % m
        roots = sqrt_mod_prime_power(target, p, k)
        assert x % m in roots
        assert all(r * r % m == target for r in roots)


def test_root_counts_for_units():
    # Odd p: 0 or 2 roots; p = 2, k >= 3: 0 or 4 roots.
    for a in range(1, 200, 2):
        assert len(sqrt_mod_prime_power(a, 2, 10)) in (0, 4)
        if a % 3:
            assert len(sqrt_mod_prime_power(a, 3, 6)) in (0, 2)


def test_sqrt_mod_composite_exhaustive():
    for n in range(1, 400):
        for a in range(n):
            assert sqrt_mod(a, n) == brute_sqrt(a, n), (a, n)


def test_invalid_arguments():
    with pytest.raises(ValueError):
        sqrt_mod_prime_power(4, 9, 2)
    with pytest.raises(ValueError):
        sqrt_mod_prime_power(4, 3, 0)
    with pytest.raises(ValueError):
        sqrt_mod(1, 0)
