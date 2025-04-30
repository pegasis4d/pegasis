#!/usr/bin/env python3

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from itertools import product
from random import choice
from textwrap import dedent

import sage.all
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.arith.misc import is_prime
from sage.arith.misc import kronecker_symbol
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.mod_poly import classical_modular_polynomial
from sage.sets.primes import Primes
from sage.rings.fast_arith import prime_range

# Arguments
parser = ArgumentParser(
    formatter_class=RawDescriptionHelpFormatter,
    epilog=dedent(
        """
    CONTROLLING ELKIES PRIMES
    So that we have small elkies primes, we search for a prime so that there
    are ELKIES elkies primes of size at most MAX_ELKIES.

    FINDING STARTING CURVE
    To find a starting curve, we must be sufficiently far away from j=1728.
    We achieve this by doing a random walk starting at j=1728*.
    For each elkies prime up to WALK_MAX_ELKIES we do WALK_STEPS number of
    steps of that size.

    *There is a small caveat: we want a curve on the surface, so our first step
    is a 2-isogeny step to land on the surface.
"""
    ),
)
parser.add_argument("-p", "--prime", type=int, help="Prime size in bits (Default: 512 bits)", default=512)
parser.add_argument("-e", "--elkies", type=int, help="Minmal number of small elkies primes", default=4)
parser.add_argument("-m", "--max-elkies", type=int, help="Maximal elkies primes to consider small", default=23)
parser.add_argument("-s", "--walk-steps", type=int, help="Steps per elkies prime in random walk", default=20)
parser.add_argument("-w", "--walk-max-elkies", type=int, help="Largest elkies prime used in random walk", default=None)
args = parser.parse_args()


print(f":: Finding parameters for")
print(f"  => log(p) = {args.prime}")
print(f"  => With at least {args.elkies} odd elkies primes of size at most {args.max_elkies}")
print(f"  => Doing random walks of length {args.walk_steps} for every prime up to {args.prime}")
print()


def small_elkies_primes(p):
    elkies_primes = []
    for ell in Primes(proof=False):
        if ell == 2:
            continue
        if ell >= args.max_elkies:
            return elkies_primes
        if kronecker_symbol(-p, ell) == 1:
            elkies_primes += [ell]
    return []


def random_elkies(j, ell):
    return choice(list(set(classical_modular_polynomial(ell, j=j).roots(GF(p), multiplicities=False)) - set([j])))


def walk(j, p):
    # Do some walking in the graph, to get away from j=1728
    for ell in prime_range(3, args.walk_max_elkies):
        print(f"Walking {args.walk_steps} times for prime {ell}")
        for _ in range(args.walk_steps):
            if kronecker_symbol(-p, ell) == 1:
                j = random_elkies(j, ell)
    return j


p = None
for e, f in product(range(args.prime, 2 * args.prime), range(1, 10 * args.prime)):
    p = f * 2**e - 1
    if is_prime(p):
        j = GF(p)(1728)
        if len(small_elkies_primes(p)) > args.elkies:
            if not on_surface(EllipticCurve(j=j)):
                j = random_elkies(j, 2)
                assert on_surface(EllipticCurve(j=j))
                break

if p is None:
    print("Failure")
    exit()

print(f":: Found prime p = {f} * 2**{e} - 1")
print(f"  => Corresponding elkies primes = {small_elkies_primes(p)}")
print()

if args.walk_max_elkies is None:
    args.walk_max_elkies = args.prime

A = None

while True:
    print()
    j = walk(j, p)
    for E in EllipticCurve(j=j).twists():
        E = E.montgomery_model()

        assert E.is_supersingular()
        assert E.base_field() == GF(E.base_field().characteristic())
        assert on_surface(E)

        A = E.a2()

        print()
        print(f"self.f = {f}")
        print(f"self.e = {e}")
        print(f"self.p = self.f * 2**self.e - 1")
        print(f"self.A = {A}")
        print(f"self.allowed_primes = {small_elkies_primes(p)}")
