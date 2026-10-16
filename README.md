# hensel_lifting

Pure-Python (3.9+, no dependencies) implementation of **Hensel lifting**: turning roots of a polynomial modulo a prime `p` into roots modulo `p^k`. On top of that, the library computes **every square root modulo any integer n**.

## Background

If `f(r) ≡ 0 (mod p)` and `f'(r) ≢ 0 (mod p)`, Hensel's lemma says `r` lifts to exactly one root modulo every `p^k`. That root comes from the p-adic Newton iteration

```
r ← r − f(r) · f'(r)⁻¹   (mod p^(2j))
```

which doubles the number of correct p-adic digits per step.

When `f'(r) ≡ 0 (mod p)` the root is singular. It can split into many roots or disappear, so lifting has to try each next digit. For example, `x² − 9` has many roots modulo `3^k`, while `x² + 1` has a root mod 2 but none mod 4.

For square roots, `x² ≡ a (mod p^k)` needs separate handling of:

| Case | Method |
|---|---|
| `p` odd, `p ∤ a` | Tonelli–Shanks mod `p`, then quadratic Newton lifting; 0 or 2 roots |
| `p = 2`, `a` odd | needs `a ≡ 1 (mod 8)` when `k ≥ 3`; one root built bit by bit, then `±r`, `±r + 2^(k−1)`; 0 or 4 roots |
| `a = p^(2h)·u` | write `x = p^h·y` with `y² ≡ u (mod p^(k−2h))`; each `y` gives `p^h` roots mod `p^k` |
| odd valuation | no roots |
| `a ≡ 0` | `x` is any multiple of `p^⌈k/2⌉` |

For composite `n`, the roots for each prime power are combined with the Chinese remainder theorem.

## Usage

```python
from hensel_lifting import Polynomial, lift_roots, lift_simple_root, sqrt_mod

f = Polynomial([-2, 0, 1])            # x^2 - 2
lift_simple_root(f, 3, 7, 10)         # the 7-adic sqrt(2) to 10 digits
lift_roots(Polynomial([-9, 0, 1]), 3, 3)   # all roots mod 27, singular branch included

sqrt_mod(4, 2 ** 5)                   # [2, 6, 10, 14, 18, 22, 26, 30]
sqrt_mod(1, 105)                      # 8 roots, combined by CRT
```

| Function | Returns |
|---|---|
| `lift_simple_root(f, r, p, k)` | the unique lift of a simple root |
| `lift_roots(f, p, k)` | all roots mod `p^k`, singular ones included |
| `tonelli_shanks(a, p)` | one square root mod an odd prime, or `None` |
| `sqrt_mod_prime_power(a, p, k)` | all square roots mod `p^k` |
| `sqrt_mod(a, n)` | all square roots mod any `n ≥ 1` |
| `egcd`, `inverse`, `crt`, `is_prime`, `factorize` | supporting arithmetic |

## Tests

```bash
python -m pip install -e ".[test]"
python -m pytest
```

Every result is compared against exhaustive search:

- all `a` modulo every prime power up to 1500 for `p ≤ 13`;
- all `a` modulo every `n < 400`;
- random polynomials modulo `3^k`, `5^k`, `7^k`.

Further tests cover root counts (0/2 for odd primes, 0/4 for powers of 2), lifts modulo `(10^6+3)^6`, and invalid input.

## Credits

The original idea came from [lapets/hensel](https://github.com/lapets/hensel) by Andrei Lapets (MIT), which lifts a single square root modulo a prime power. This repository is an independent rewrite with a different design and a wider scope: general polynomials, singular roots, all roots, `p = 2`, and composite moduli.

## License

MIT. See `LICENSE`.
