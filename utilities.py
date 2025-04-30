#!/usr/bin/env python3

from sage.arith.misc import gcd
from sage.rings.integer_ring import ZZ
from sage.structure.factorization import Factorization
from sage.rings.finite_rings.integer_mod import Mod

from const_precomp import small_sos


def remove_common_factors(n, r):
    r"""Remove all primes dividing r from n

    If n = \prod_i p_i^e_i and r = \prod_i p_i^f_i, returns
    tuple
        (remainder, common_factors)
    where
      common_factors = \prod_i p_i^{min(e_i, f_i)}
    and
      remainder = n / common_factors
    """

    common_factors = 1
    remainder = n

    while (g := gcd(remainder, r)) != 1:
        common_factors *= g
        remainder /= g

    assert common_factors * remainder == n

    return remainder, common_factors


def remove_primes(n, primes, odd_parity_only=True):
    """Remove all powers of primes in primes from n

    Input:
      - n: int/Integer
      - primes: List of prime numbers
      - odd_parity_only: (optional) For every prime p in primes, remove only as
        many powers of p until p has even multiplicity in n

    Output:
      - Tuple (remainder, prime_powers_contained)

        Where remainder is n with all the power of primes removed, and
        prime_powers_contained the product of all primes and their
        multiplicities in n. In other words

           n == remainder * prime_powers_contained
    """
    remainder = ZZ(n)
    prime_powers_contained = 1

    for prime in primes:
        val = remainder.valuation(prime)

        if odd_parity_only:
            val = val % 2

        prime_powers_contained *= prime ** val
        remainder /= prime ** val

    assert prime_powers_contained * remainder == n

    return remainder, prime_powers_contained


def two_squares_factored(factors):
    """
    This is the function `two_squares` from sage, except we give it the
    factorisation of n already.
    """
    F = Factorization(factors)
    for p, e in F:
        if e % 2 == 1 and p % 4 == 3:
            raise ValueError("%s is not a sum of 2 squares" % n)

    n = F.expand()
    if n == 0:
        return (0, 0)
    a = ZZ.one()
    b = ZZ.zero()
    for p, e in F:
        if p == 1:
            continue
        if e >= 2:
            m = p ** (e // 2)
            a *= m
            b *= m
        if e % 2 == 1:
            if p == 2:
                # (a + bi) *= (1 + I)
                a, b = a - b, a + b
            else:  # p = 1 mod 4
                if p in small_sos:
                    r, s = small_sos[p]
                else:
                    # Find a square root of -1 mod p.
                    # If y is a non-square, then y^((p-1)/4) is a square root of -1.
                    y = Mod(2, p)
                    while True:
                        s = y ** ((p - 1) / 4)
                        if not s * s + 1:
                            s = s.lift()
                            break
                        y += 1
                    # Apply Cornacchia's algorithm to write p as r^2 + s^2.
                    r = p
                    while s * s > p:
                        r, s = s, r % s
                    r %= s

                # Multiply (a + bI) by (r + sI)
                a, b = a * r - b * s, b * r + a * s

    a = a.abs()
    b = b.abs()
    assert a * a + b * b == n
    if a <= b:
        return (a, b)
    else:
        return (b, a)


def on_surface(E):
    try:
        E.torsion_basis(2)
    except ValueError:
        return False
    try:
        E.torsion_basis(4)
    except ValueError:
        return True

    return False
