#!/usr/bin/env python3

import itertools as itt
from random import randint

from sage.arith.misc import next_prime
from sage.matrix.constructor import Matrix
from sage.misc.functional import isqrt
from sage.arith.misc import kronecker
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer_ring import ZZ


def element_norm(el, p):
    return el[0] ** 2 + p * el[1] ** 2


def element_multiply(el1, el2, p):
    a1, b1 = el1
    a2, b2 = el2
    return [a1 * a2 - p * b1 * b2, a1 * b2 + b1 * a2]


def element_scale(el, scal):
    a, b = el
    return [a / scal, b / scal]


def _norm(v):
    return sum([x**2 for x in v])


def _add(a, b):
    return [ia + ib for ia, ib in zip(a, b)]


def _mul(v, n):
    return [x * n for x in v]


def _is_zero(v):
    return all([x == 0 for x in v])


def short_vectors(L, B, cf_b=15):
    """
    Enumerate vectors in L with norm bounded by B.
    Tries linear combinations with coefficients up to cf_b
    """
    sv = {}
    assert len(L) == 2  # Avoiding sign repeats
    for cff in itt.product(range(-cf_b, cf_b), range(0, cf_b)):
        # for cff in itt.product(range(-cf_b, cf_b), repeat=len(L)):
        v = [0 for _ in range(len(L[0]))]
        for i in range(len(L)):
            v = _add(v, _mul(L[i], cff[i]))
        nv = _norm(v)
        if nv < B and nv != 0 and not nv in sv:
            sv[nv] = v
    sv = [(k, sv[k]) for k in sv]
    sv.sort()
    return sv


def short_ideals(ideal, ideal_norm, norm_bound, p, spr=None, cf_b=15):
    """
    For a given ideal returns the shortest vectors (or equivalently the smallest
    equivalent ideals).

    Input:
    - I: an ideal as list of generators, where the element [a1, a2] is a1 + i*a2
    - N: the norm of I
    - B: bound on the (reduced) norm of the elements returned
    - p: the field prime
    - spr: isqrt(p) if precomputed # TODO: probably remove
    - cf_b: number of combinations to try in ShortVectors
    Output:
    - res: a list of pairs (n, [a, b, c]) where [a, b] are generators (as list) of
        and ideal I1 equivalent to I, n is the norm of I1 and c is the short element
        in I sending I to I1
    """
    gens = [list(gen) for gen in ideal.gens()]

    if not spr:
        spr = isqrt(p)

    L = Matrix(ZZ, 2, [gens[0][0], spr * gens[0][1], gens[1][0], spr * gens[1][1]])

    L = L.LLL()
    L = [L[0], L[1]]
    res = []
    sv2 = short_vectors(L, ideal_norm * norm_bound, cf_b=cf_b)
    for _, sh in sv2:
        cshel = [sh[0], -sh[1] / spr]
        if sh[1] % spr != 0:
            print("Non divisible")
        idl = [element_scale(element_multiply(gens[0], cshel, p), ideal_norm), element_scale(element_multiply(gens[1], cshel, p), ideal_norm), cshel]
        res.append((element_norm(cshel, p) / ideal_norm, idl))
    return res


def reduced_basis(I_basis):
    # LLL is way overkill here, change later
    def _matrix_to_gens(M, B):
        return [sum(c * g for c, g in zip(row, B)) for row in M]

    def gram_matrix(basis):
        M = []
        for a in basis:
            M.append([(a * b.conjugate()).trace() for b in basis])
        return Matrix(QQ, M)

    G = gram_matrix(I_basis)
    U = G.LLL_gram().transpose()
    reduced_basis_elements = _matrix_to_gens(U, I_basis)
    return reduced_basis_elements


def principal_generator(frak_a):
    N = frak_a.norm()

    assert len(frak_a.gens()) == 2, "The ideal has only one generator (which is not really an issue)"
    for gen in reduced_basis(frak_a.gens()):
        if gen.norm() == N:
            return gen
    assert False, "WARNING: Not a principal ideal"


def conjugate(order, frak_a):
    a, b = frak_a.gens_two()
    return order * a.conjugate() + order * b.conjugate()


def random_degree_one_ideal(order, gen, p):
    ell = next_prime(randint(-order.discriminant(), -10 * order.discriminant()))
    while kronecker(-p, ell) != 1:
        ell = next_prime(randint(-order.discriminant(), -10 * order.discriminant()))

    lam = Integer(GF(ell)(-p).sqrt())

    frak_ell = order * ell + order * (gen - lam)

    assert frak_ell.norm() == ell
    return ell, frak_ell


def random_fixed_prime_ideal(order, gen, p, ell):
    lam = Integer(GF(ell)(-p).sqrt())

    frak_ell = order * ell + order * (gen - lam)

    assert frak_ell.norm() == ell
    return ell, frak_ell


def ideal_to_sage(I, O):
    """
    Converts an ideal as a list of lists into
    an ideal in the order O
    """
    I[0] = [QQ(i) for i in I[0]]
    I[1] = [QQ(i) for i in I[1]]

    pi = O.gens()[1]

    g1 = I[0][0] + pi * I[0][1]
    g2 = I[1][0] + pi * I[1][1]
    return O * g1 + O * g2
