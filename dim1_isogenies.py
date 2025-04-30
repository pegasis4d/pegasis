#!/usr/bin/env python3

import logging

from xonly import xPoint, isWeierstrass, translate_by_T
from basis_sampling import find_Ts
from elkies import Elkies

from sage.arith.misc import factor
from sage.rings.integer import Integer
from sage.arith.misc import kronecker_symbol
from sage.rings.finite_rings.finite_field_constructor import GF

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
logger_sh.setLevel(logging.WARNING)
formatter = logging.Formatter("%(name)s [%(levelname)s] %(message)s")
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)


def ideal_above_2_action(E, P, Q, basis_order, ideal, w):
    r"""Computes the action of a power of an ideal above 2 using Velu

    Input:
      - E: Elliptic curve over Fp
      - P, Q: Basis of E[2^basis_order],
              Here, P is in E(Fp) and Q is in Et(Fp) for a quadratic twist Et
      - basis_order: Such that 2^basis_order is the order of P, Q
      - ideal: Power of an ideal above 2
    Output:
      Tuple (phi, E_ideal)

      - phi: The isogeny corresponding to ideal
      - E_ideal: The curve obtained by the action of the ideal on E
    """
    _, v2 = factor(ideal.norm())[0]

    left = 0
    right = 0
    if (w - 1) / 2 in ideal:
        K = P
        left = v2
    else:
        assert (w + 1) / 2 in ideal
        K = Q
        right = v2

    K = K.xMUL(2 ** (basis_order - v2))

    assert K.xMUL(2 ** (v2 - 1))
    assert not K.xMUL(2**v2)

    phi_2 = K.xMUL(2 ** (v2 - 1)).xISOG(E, 2)
    phi = phi_2

    for step in range(1, v2):
        E = phi.codomain()
        K = K.push(phi_2)

        assert K.xMUL(2 ** (v2 - step - 1))
        assert not K.xMUL(2 ** (v2 - step))

        phi_2 = K.xMUL(2 ** (v2 - step - 1)).xISOG(E, 2)
        phi = phi_2 * phi

    E_out = phi.codomain().montgomery_model()

    return phi.codomain().isomorphism_to(E_out) * phi, E_out, left, right


def small_prime_ideal_action(E, ell, lam=None, prev_j=None):
    r"""Computes action of a ideal with prime norm using Elkies

    Input:
      - E:      Elliptic curve over Fp IN WEIERSTRASS FORM
      - ell:    A small prime, split in the quadratic order
      - lam:    (optional) Eigenvalue determining the ideal
      - prev_j: (optional) The j-invariant of the WRONG ell-isogenous curve

    Both lam and prev_j are passed to Elkies

    Output:
      - phi: Isogeny corresponding to an ideal above ell
    """
    # Fix later, but Elkies uses short Weierstrass
    # Change ring stuff is messy, but needed to actually compute the rational model
    assert isWeierstrass(E)

    h = Elkies(E, ell, E.base_field(), lam=lam, prev_j=prev_j)
    phi = E.isogeny(h)

    return phi


def smooth_ideal_action(E, norm, ideal, w, max_order):
    r"""Computes action of a an ideal with odd smooth norm using successive Elkies

    Input:
      - E:     Elliptic curve over Fp
      - norm:  Norm of ideal
      - ideal: An ideal

    Output: Tuple (isogenies, E_ideal)
      - isogenies: A list of isogenies, forming a chain from E to E_out
      - E_ideal:   The curve obtained by the action of the ideal on E
    """

    p = E.base_field().characteristic()

    isogenies = []
    E1 = E.short_weierstrass_model()
    isogenies += [E.isomorphism_to(E1)]

    if norm > 1:
        for ell, e in factor(norm):

            lam = Integer(GF(ell)(-p).sqrt())
            if not (w - lam) in (ideal + ell * max_order):
                # In this case, we want the other one
                lam = ell - lam

            assert (w - lam) in (ideal + ell * max_order)

            prev_j = None
            for _ in range(e):
                logger.debug(f"Starting Elkies step of degree {ell}")
                phi = small_prime_ideal_action(E1, ell, lam, prev_j)
                prev_j = E1.j_invariant()
                E1 = phi.codomain()
                isogenies.append(phi)

    E_ideal = isogenies[-1].codomain().montgomery_model()
    isogenies += [isogenies[-1].codomain().isomorphism_to(E_ideal)]

    return isogenies, E_ideal


def random_smooth_isogeny(E, g):
    r"""
    Returns a random isogeny of smooth degree g from E
    Input:
        - E: Elliptic curve over Fp
        - g: Smooth number
    Output:
        - phi : An isogeny from E of degree g
        - E_out : phi(E)
    """
    p = E.base_field().characteristic()
    # Random smooth odd norm isogeny
    isogs = []
    if g > 1:
        E1 = E.short_weierstrass_model()
        pre_isom = E.isomorphism_to(E1)
        isogs.append(pre_isom)
        logger.debug(f"Extra isogeny for g = {factor(g)}")
        # Do elkies first, to stay over Fp
        g_elkies = []
        for ell, e in factor(g):
            if kronecker_symbol(-p, ell) == 1:
                g_elkies.append((ell, e))
        for ell, e in g_elkies:
            prev_j = None
            for _ in range(e):
                phi_ell = small_prime_ideal_action(E1, ell, prev_j=prev_j)
                prev_j = E1.j_invariant()
                E1 = phi_ell.codomain()
                isogs.append(phi_ell)

        E_out = isogs[-1].codomain().montgomery_model()
        post_isom = isogs[-1].codomain().isomorphism_to(E_out)
        isogs.append(post_isom)

        return isogs, E_out
    else:
        return [], E


def eval_endomorphism(rho, P, twist, max):
    r"""
    Evaluates an element of End(E) on a point P.
    P is assumed to be Fp rational on either E or on a twist
    """
    a, b = list(rho)  # write as a + b*pi
    n = rho.denominator()
    if twist:
        # a + b*pi = a - b, since pi(P) = -P
        m = Integer(a - b)
    else:
        # a + b*pi = a + b, since pi(P) = P
        m = Integer(a + b)
    if max and (n == 2):  # See Appendix D.2
        mP = P.xMUL(m)
        E = P.curve

        T0 = find_Ts(E, only_T0=True)

        # Already know which point to choose, Remark D.2
        assert (P.X - T0.x()).is_square() != twist

        T = xPoint(T0.x(), E)

        return translate_by_T(mP, T)

    else:
        return P.xMUL(m)
