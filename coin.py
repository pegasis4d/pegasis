#!/usr/bin/env python3

from copy import copy

from sage.arith.misc import is_pseudoprime, xgcd
from sage.functions.other import floor, ceil
from sage.misc.functional import isqrt, log
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.misc.misc_c import prod

from utilities import remove_common_factors, two_squares_factored, remove_primes

# Primes 1 mod 4 / 3 mod 4
from const_precomp import good_prime_prod_10k, bad_prime_prod_10k

good_prime_prod = good_prime_prod_10k
bad_prime_prod = bad_prime_prod_10k


def coins(n1, n2, e, exp_bound=None):
    """Yield pairs (u, v) s.t. u*n1 + v*n2 = 2^t for t <= e

    Arguments:
      - n1, n2:               Integers
      - e:                    Largest 2-power exponent
      - exp_bound (optional): Stop searching when t > t_min + exp_bound.
                              Where t_min is the smallest 2-power for which
                              there is a solution.
    """

    g, a, b = xgcd(n1, n2)

    if g != 1:
        raise ValueError(f"{n1 = } and {n2 = } are not coprime")

    # Smallest t for which we have found a solution
    t_min = None

    for t in range(floor(log(min(n1, n2), 2)), e + 1):

        if t_min and exp_bound and t > t_min + exp_bound:
            break

        N = 2**t

        # Find conditions for u > 0 and v > 0 and solve for k
        # (using definitions as below)
        for k in range(ceil(-N * a / n2), floor(N * b / n1) + 1):

            if not t_min:
                # Found minimal solution
                t_min = t

            # u * n1 + v * n2 = N * (a * n1 + b * n2) = g * N = N
            u = N * a + k * n2
            v = N * b - k * n1

            if u % 2 == 0 and v % 2 == 0:
                # We have seen this solution for a previous power of 2
                continue

            assert u > 0
            assert v > 0
            assert u * n1 + v * n2 == 2**t

            yield u, v, t


def is_sum_of_2_squares(numbers, early_abort=False):
    """Determine whether every number in numbers is sum of (two) squares

    We know that n can be written as the sum of two squares if and only if every
    prime p dividing n that is 3 mod 4 divides n with even multiplicity

    If any of the numbers cannot be written as sum of squares, we try to abort
    as soon as possible.

    Input:
      - numbers: List of ints/Integers

        Note: Also accepts a single int/Integer n (you do not need to cast
        single element to singleton list [n])

    Output:
      - List [(success, pseudo_factorization), ...] one entry for every number
        in numbers

        Note: If numbers is  just an int, or Integer, then just a tuple
            (success, partial_factorization)
        is returned, not a singleton
            [(success, pseudo_factorization)].

        For the i-th entry (success, pseudo_factorization) in output list

        - success (True/False) tells us whether i-th number in numbers can be
          written as sum of two squares

        - pseudo_factorization is a list of tuples [(p_1, e_1), (p_2, e_2), ...]
          where p_j divides the i-th number in numbers with multiplicity e_j

          Every p_j is a prime number, with the exception of the last one, which
          is only checked for pseudoprimality
    """

    single = False
    if isinstance(numbers, (int, Integer)):
        single = True
        numbers = [ZZ(numbers)]

    # List of pairs (status, factorizations) for every number
    #   - status is one of [-1, 0, 1], where
    #     * -1 = do not know
    #     *  0 = failure
    #     *  1 = success
    #
    #     we need three status codes (more than success/failure), becuase we
    #     want to be able to abort early and so we need an "I don't know" option
    #
    #   - factorizations is a list of (partial) factorizations of the numbers
    #     If there is an early abort, it is possible that the factorization is
    #     only partial

    result = [[-1, []] for _ in numbers]

    # Numbers with factors removed (Iteratively updated throughout)
    # Must perform copy to prevent changing numbers outside function call as
    # python essentially passes mutable types by reference
    reduced_numbers = copy(numbers)

    for i, n in enumerate(reduced_numbers):
        # Any power of two is the sum of two squares, we know this for free
        v2 = n.valuation(2)
        remainder = n // (2**v2)
        # Update factorization
        result[i][1] += [(2, v2)]

        if remainder % 4 == 3:
            # n cannot be written as sum of two squares
            # Tell later steps not to operate further
            result[i][0] = 0

            if early_abort:
                # If there is only one number, do not return list
                return result[0] if single else result

        # Update list of numbers
        reduced_numbers[i] = remainder

    for i, n in enumerate(reduced_numbers):
        if result[i] == 0:
            # Previous step has determined that this number cannot be written as
            # sum of two squares
            continue

        # Get product of all "bad" primes that divide n
        remainder, bad_part = remove_common_factors(n, bad_prime_prod)

        # Check whether bad_part is a square
        bad_part_isqrt = isqrt(bad_part)
        if bad_part_isqrt**2 == bad_part:
            # Update the factorization
            result[i][1] += [(bad_part_isqrt, 2)]
        else:
            # n cannot be written as sum of two squares
            # Tell later steps not to operate further
            result[i][0] = 0

            if early_abort:
                # If there is only one number, do not return list
                return result[0] if single else result

        # Update list of numbers
        reduced_numbers[i] = remainder

    for i, n in enumerate(reduced_numbers):
        if result[i] == 0:
            # Previous step has told us, that this number is not sum of two
            # squares we skip further processing
            continue

        # Get product of all "good" primes that divide n
        remainder, good_part = remove_common_factors(n, good_prime_prod)

        # Update factorization
        # Cast Factorization instance to list so we can append
        result[i][1] += list(ZZ(good_part).factor()) + [(remainder, 1)]

        # If not successful
        if not is_pseudoprime(remainder):
            # n cannot be written as sum of two squares
            result[i][0] = 0

            if early_abort:
                # If there is only one number, do not return list
                return result[0] if single else result
        else:
            # n can be written as sum of two squares
            result[i][0] = 1

        # Update list of numbers
        reduced_numbers[i] = remainder

    # Verify factorizations are correct
    # Disabled when python is called with -O flag
    if __debug__:
        for i, (success, factorization) in enumerate(result):
            if success == 1:
                assert prod([p**e for p, e in factorization]) == numbers[i]

    # If there is only one number, do not return list
    return result[0] if single else result


def sum_of_2_squares(numbers, early_abort=False):
    """Write every number in numbers as sum of two squares if possible

    If any of the numbers cannot be written as sum of squares, we try to abort
    as soon as possible.

    See also: is_sum_of_2_squares

    Input:
      - numbers: List of ints/Integers

        Note: Also accepts a single int/Integer n (you do not need to cast
        single element to singleton list [n])

    Output:
      - Tuple (success, list_of_pairs)

        Where
        - success (True/False) tells us whether

      - List [(success, pseudo_factorization), ...] one entry for every number
        in numbers

        Note: If numbers is  just an int, or Integer, then just a tuple
            (success, partial_factorization)
        is returned, not a singleton
            [(success, pseudo_factorization)].

        For the i-th entry (success, pseudo_factorization) in output list

        - success (True/False) tells us whether i-th number in numbers can be
          written as sum of two squares

        - pseudo_factorization is a list of tuples [(p_1, e_1), (p_2, e_2), ...]
          where p_j divides the i-th number in numbers with multiplicity e_j

          Every p_j is a prime number, with the exception of the last one, which
          is only checked for pseudoprimality
    """

    single = False
    if isinstance(numbers, (int, Integer)):
        single = True
        numbers = [ZZ(numbers)]

    are_sos = is_sum_of_2_squares(numbers, early_abort=early_abort)

    result = []
    for success, factorization in are_sos:
        if success == 1:
            x, y = two_squares_factored(factorization)
            result += [(True, x, y)]
        else:
            result += [(False, 0, 0)]

    # If there is only one number, do not return list
    return result[0] if single else result


def sos_coins(n1, n2, e, allowed_primes, n_squares=1, exp_bound=None):
    """Yield tuples to solve norm equation u * n1 + v * n2 = 2 ** t

    Input:
      - n1, n2:         int/Integer
      - e:              Maximal two-exponent t to return for
      - n_squares:      How many of u and v must be sum of squares (1 or 2,
                        default 1)
      - allowed primes: List of primes that can be removed from u and v
                        when checking for sum of squares
                        e.g. u = 3 * (x_u ** 2 + y_y ** 2)
      - exp_bound: (optional) Argument to pass to coin_pair_generator

    Output:
      Tuple ((success_u, x_u, y_u, g_u), (success_v, x_v, y_v, g_v), t)

      Such that with

          u = | g_u * (x_u ** 2 + y_u ** 2)  if u/g_u is sum of two squares
              | u                            else

          v = | g_v * (x_v ** 2 + y_v ** 2)  if v/g_v is sum of two squares
              | v                            else

      we have

          u * n1 + v * n2 = 2^t

      with t <= e

      In both cases, g_{u, v} is a product of primes in allowed_primes with
      multiplicity at most 1

      If n_squares = 1, it is possible that only one of u/g_u, v/g_v is sum of
      two squares.
      If n_squares = 2, both must be.
    """

    early_abort = True if n_squares == 2 else False

    for u, v, t in coins(n1, n2, e, exp_bound):

        u, g_u = remove_primes(u, allowed_primes, odd_parity_only=True)
        v, g_v = remove_primes(v, allowed_primes, odd_parity_only=True)

        (success_u, x_u, y_u), (success_v, x_v, y_v) = sum_of_2_squares([u, v], early_abort=early_abort)

        if n_squares == 1 and (success_u or success_v):
            raise NotImplementedError("n_squares = 1 is not implemented yet")

        elif n_squares == 2 and (success_u and success_v):
            # In the n_squares = 2 case, we know that both u, v are sums of two
            # squares, so success_{u, v} are redudant
            # ...however, when we implement n_squares = 1, the caller will want
            # to know which one of u, v are sum of two squares, and so we will
            # need some sort of success variable
            # To keep a stable api, we leave the succcess_{u, v} variables here
            # for now
            yield (success_u, x_u, y_u, g_u), (success_v, x_v, y_v, g_v), t
