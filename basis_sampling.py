#!/usr/bin/env python3

import logging
from time import time

from xonly import xPoint, MontgomeryA

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
logger_sh.setLevel(logging.WARNING)
formatter = logging.Formatter("%(name)s [%(levelname)s] %(message)s")
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)


def find_Ts(E, only_T0 = False):
    r"""
    Given a curve E, finds and marks the non-trivial
    2-torsion points according to Lemma D.1
    """
    A = MontgomeryA(E)
    F =  E.base_field()
    R = F["X"]
    X = R.gens()[0]
    f = X**2 + A*X + 1
    logger.debug("     > root finding")
    lam1, lam2 = f.roots(multiplicities=False)

    R1 = E(lam1, 0)
    R2 = E(lam2, 0)
    R3 = E(0, 0)

    #Find T0
    Rs = [R1, R2]
    for T in Rs:
        if T.tate_pairing(T, 2, 1) != 1:
            T0 = T
            Rs.remove(T)
            break

    assert T0
    if only_T0:
        return T0

    assert T0.tate_pairing(T0, 2, 1) == -1
    Rs.append(R3)
    for T in Rs:
        if T.tate_pairing(T0, 2, 1) == 1:
            Tm1 = T
            Rs.remove(T)
            break

    assert Tm1
    T1 = Rs[0]
    assert T1.tate_pairing(T0, 2, 1) != 1

    return T0, Tm1, T1


def TwoTorsBasis(E, e):
    r"""
    Fast sampling of a basis P, Q of E[2**e], such that x(P) and x(Q) are both defined over Fp
    Input:
        - E: Elliptic curve over Fp
        - e: Exponent
    Output:
        - P, Q: Basis of E[2**e] so that P is in E(Fp), and Q is in E^t(Fp) for an Fp-twist of E.
    """
    logger.debug("     > Finding TwoTorsionBasis")
    tstart = time()

    T0, Tm1, T1 = find_Ts(E)

    A = MontgomeryA(E)
    F =  E.base_field()
    R = F["X"]
    X = R.gens()[0]
    f = X**2 + A*X + 1

    logger.debug(f"     > Done, have used {time()-tstart} sec")

    logger.debug("     > sample xP")
    xT0 = T0.x()
    xP = xT0 + F.random_element()**2
    while not (f(xP)*xP).is_square() or (xP-Tm1.x()).is_square():
        xP = xT0 + F.random_element()**2

    logger.debug("     > sample xQ")
    xQ = xT0 - F.random_element()**2
    while (f(xQ)*xQ).is_square() or not ((xQ-T1.x()).is_square()):
        xQ = xT0 - F.random_element()**2

    logger.debug(f"    > Done, have used {time()-tstart} sec")
    #print(T0.tate_pairing(E.lift_x(xP), 2, 1))

    P = xPoint(xP, E)
    Q = xPoint(xQ, E)

    assert (E.base_field().characteristic()+1) % 2**(e+1) == 0
    cofac = (E.base_field().characteristic()+1)/2**(e+1)
    P = P.xMUL(cofac)
    Q = Q.xMUL(cofac)

    assert P.xMUL(2**(e-1))
    assert not P.xMUL(2**e)
    assert Q.xMUL(2**(e-1))
    assert not Q.xMUL(2**e)
    assert Q.xMUL(2**(e-1)) != P.xMUL(2**(e-1))

    logger.debug(f"    > Total time for 2-tors basis finding: {time()-tstart}")

    return P, Q
