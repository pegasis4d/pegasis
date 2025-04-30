#!/usr/bin/env python3

from time import time
import logging

from sage.arith.misc import valuation, inverse_mod
from sage.calculus.var import var
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ

from ideals import ideal_to_sage, conjugate, principal_generator
from basis_sampling import TwoTorsBasis
from dim1_isogenies import ideal_above_2_action, smooth_ideal_action, eval_endomorphism, random_smooth_isogeny

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter("%(name)s [%(levelname)s] %(message)s")
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)


def ideal_to_kernel(E, uv, frak_bc, max_order, w, e):
    r"""
    Following sec. 5.4, prepare the data necessary for the 4D isogeny.
    Input:
        - E: starting curve
        - uv: the decomposition of u and v as cof + sum of squares
        - frak_bc: the ideals b and c as output of UVSolver
        - e: the torsion we want to use

    Output:
        - N1, N2, x_u, y_u, g_u, x_v, y_u, g_v e s.t.
            N1*g_u*(x_u^2 + y_u^2) + N2*g_v*(x_v^2 + y_v^2) == 2^e
        - P_u, Q_u, P_v, Q_v: image points on E_u, Ev as described
            in sec 5.4
    """
    p = E.base_field().characteristic()
    Fp2 = GF((p, 2), name="i", modulus=var("x") ** 2 + 1)

    # First recover the numbers
    tstart = time()
    if type(uv[0]) == tuple:
        (x_u, y_u, g_u) = uv[0]
    else:
        x_u = None
        y_u = None
        g_u = 1

    if type(uv[1]) == tuple:
        (x_v, y_v, g_v) = uv[1]
    else:
        x_v = None
        y_v = None
        g_v = 1

    # (x_u, y_u, g_u), (x_v, y_v, g_v) = uv
    b_list, c_list = frak_bc

    N1, b_cof, b_twopower, frak_b_list = b_list
    N2, c_cof, c_twopower, frak_c_list = c_list

    # assert N1 * g_u * (x_u**2 + y_u**2) + N2 * g_v * (x_v**2 + y_v**2) == 2**e

    def div_by_n(ideal, n):
        generators = [gi / n for gi in ideal.gens()]
        return sum([max_order * gi for gi in generators])

    frak_b = ideal_to_sage(frak_b_list, max_order)
    frak_c = ideal_to_sage(frak_c_list, max_order)

    if b_twopower > 1:  # factors through two
        frak_b = div_by_n(frak_b, 2)
    if c_twopower > 1:
        frak_c = div_by_n(frak_c, 2)

    # The part thats the same for b and c
    frak_bc_same = frak_b + frak_c

    # We only want the "unique" part of b and c
    frak_b = div_by_n(frak_b * conjugate(max_order, frak_bc_same), frak_bc_same.norm())
    frak_c = div_by_n(frak_c * conjugate(max_order, frak_bc_same), frak_bc_same.norm())

    frak_b_2 = frak_b + max_order * ((2 ** frak_b.norm().valuation(2)))
    frak_c_2 = frak_c + max_order * ((2 ** frak_c.norm().valuation(2)))

    odd_part = Integer(frak_bc_same.norm()).prime_to_m_part(2)
    pow2_part = frak_bc_same.norm() / odd_part
    frak_bc_same_2 = frak_bc_same + max_order * pow2_part
    frak_bc_same_odd = frak_bc_same + max_order * odd_part

    b_cof = b_cof / odd_part
    c_cof = c_cof / odd_part

    # First step to the starting curve given by same direction steps
    ee = valuation(p + 1, 2) - 1
    logger.debug("Walking to starting point (i.e. steps that are same for b and c)...")
    if frak_bc_same.norm() > 1:
        if pow2_part > 1:
            P_temp, Q_temp = TwoTorsBasis(E, ee)
            # P_temp, Q_temp = TwoTorsBasis(E, two_pow)
            _, E, _, _ = ideal_above_2_action(E, P_temp, Q_temp, ee, frak_bc_same_2, w)
        if odd_part > 1:
            _, E = smooth_ideal_action(E, odd_part, frak_bc_same_odd, w, max_order)
    logger.debug(f"Done! Have used {time() - tstart}")
    P, Q = TwoTorsBasis(E, ee)
    # P, Q = TwoTorsBasis(E, ee)

    # Doing the Elkies steps for b
    logger.debug("odd Elkies steps for b")
    isogs_cof, E1 = smooth_ideal_action(E, b_cof, frak_b, w, max_order)
    P1, Q1 = P, Q
    for phi in isogs_cof:
        P1 = P1.push(phi)
        Q1 = Q1.push(phi)
    logger.debug(f"Done! Have used {time() - tstart}")

    # Even Elkies steps for b
    if frak_b_2.norm() > 1:
        phi_2, E1, left_b, right_b = ideal_above_2_action(E1, P1, Q1, ee, frak_b_2, w)
        P1 = P1.push(phi_2)
        Q1 = Q1.push(phi_2)
    else:
        left_b = 0
        right_b = 0

    # Adjust the order
    P1, Q1 = P1.xMUL(2 ** (ee - left_b - (e + 2))), Q1.xMUL(2 ** (ee - right_b - (e + 2)))

    assert P1.xMUL(2 ** (e + 1))
    assert not P1.xMUL(2 ** (e + 2))
    assert Q1.xMUL(2 ** (e + 1))
    assert not Q1.xMUL(2 ** (e + 2))
    # At this point, P1, Q1 are the "starting points"

    logger.debug("Doing g_u isogeny")
    g_u_isogs, Eu = random_smooth_isogeny(E1, g_u)
    Pu, Qu = P1, Q1
    for phi_gu in g_u_isogs:
        Pu = Pu.push(phi_gu)
        Qu = Qu.push(phi_gu)
    logger.debug(f"Done! Have used {time() - tstart}")

    # Now the second part
    # Evaluating the endomorphism
    logger.debug("Evaluating the endomorphism frak_b*frac_c_bar")
    rho_bc = principal_generator(frak_b * conjugate(max_order, frak_c))
    rho_bc_val_2 = rho_bc.norm().valuation(2)

    # assert (frak_b*frak_c).norm() == rho_bc.norm()
    # assert rho_bc in frak_b
    # assert rho_bc.conjugate() in frak_c

    # adjust_2 = 0
    # if rho_bc_val_2 > 0: #Here we really need to clear denominators to evaluate the endomorphism
    #    rho_bc *= 2
    #    adjust_2 = 1

    # print("!!!!!!!!!")
    # print(ee-(rho_bc_val_2))
    # print(e+2)
    # assert ee-(rho_bc_val_2 + adjust_2) >= e+2, "NotImplemented: We can drop adjust_2 if we allow extension-fields"
    # print(f"Evaluating {rho_bc}")
    max = ee - rho_bc_val_2 == e + 2

    Pc = P
    Qc = Q
    # Then apply endomorphism
    Pc = eval_endomorphism(rho_bc, Pc, False, max)
    Qc = eval_endomorphism(rho_bc, Qc, True, max)

    # Act with 2-ideals above c
    logger.debug("Acting with 2-part of frak_c")
    E2 = E
    frak_c_2_norm = frak_c_2.norm()
    if frak_c_2_norm > 1:
        phi_2, E2, _, _ = ideal_above_2_action(E2, P, Q, ee, frak_c_2, w)
        Pc = Pc.push(phi_2)
        Qc = Qc.push(phi_2)
        c_2_part = ZZ(frak_c_2_norm).valuation(2)
    else:
        c_2_part = 0

    # Finally, odd part of c
    isogs_cof, E2 = smooth_ideal_action(E2, c_cof, frak_c, w, max_order)
    for phi in isogs_cof:
        Pc = Pc.push(phi)
        Qc = Qc.push(phi)

    inverse_norms = inverse_mod(c_cof, 2**ee)
    assert (inverse_norms * c_cof) % 2**ee == 1

    logger.debug(f"Done! Have used {time() - tstart}")
    P2 = Pc.xMUL(2 ** (ee - c_2_part - left_b - (e + 2)) * inverse_norms)
    Q2 = Qc.xMUL(2 ** (ee - c_2_part - right_b - (e + 2)) * inverse_norms)

    assert P2.xMUL(2 ** (e + 1))
    assert not P2.xMUL(2 ** (e + 2))
    assert Q2.xMUL(2 ** (e + 1))
    assert not Q2.xMUL(2 ** (e + 2))

    logger.debug("Doing g_v isogeny")
    Pv, Qv = P2, Q2

    g_v_isogs, Ev = random_smooth_isogeny(E2, g_v)

    for phi_gv in g_v_isogs:
        Pv = Pv.push(phi_gv)
        Qv = Qv.push(phi_gv)
    logger.debug(f"Done! Have used {time() - tstart}")

    Ev = Pv.curve.change_ring(Fp2)
    Eu = Pu.curve.change_ring(Fp2)

    Pv = Ev.lift_x(Pv.X)
    Qv = Ev.lift_x(Qv.X)

    Pu = Eu.lift_x(Pu.X)
    Qu = Eu.lift_x(Qu.X)

    # Make sure to lift correctly
    if Pu.weil_pairing(Qu, 2 ** (e + 2)) ** (g_v * N1 * N2) != Pv.weil_pairing(Qv, 2 ** (e + 2)) ** (g_u):
        Qv = -Qv

    # Double check full order
    assert Pu.weil_pairing(Qu, 2 ** (e + 2)) ** (2 ** (e + 2)) == 1
    assert Pu.weil_pairing(Qu, 2 ** (e + 2)) ** (2 ** (e + 1)) != 1
    assert Pv.weil_pairing(Qv, 2 ** (e + 2)) ** (2 ** (e + 2)) == 1
    assert Pv.weil_pairing(Qv, 2 ** (e + 2)) ** (2 ** (e + 1)) != 1
    # Weil pairing check
    assert Pu.weil_pairing(Qu, 2 ** (e + 2)) ** (g_v * N1 * N2) == Pv.weil_pairing(Qv, 2 ** (e + 2)) ** (g_u)

    return N1, N2, x_u, y_u, g_u, x_v, y_v, g_v, Pu, Qu, Pv, Qv
