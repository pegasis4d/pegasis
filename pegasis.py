from time import time
import logging

from sage.all import *

from xonly import xPoint, random_xPoint, MontgomeryA, isWeierstrass, translate_by_T
from lcsidh import UVSolver
from uv_params import UV_params
from ideals import ideal_to_sage, PrincipalGenerator, Conjugate, RandomDegreeOnePrimeIdeal
from elkies import Elkies

from theta_lib.isogenies.Kani_clapoti import KaniClapotiIsog

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

class PEGASIS:
    def __init__(self, level, params=None):
        r"""
        Main Clapoti class. Takes as input the security level
        (500|1000|1500|2000|4000) and generates the necessary parameters and
        structures to solve the norm equation and run the Clapoti algorithm.
        Additional parameters can be supplied as a dictionary.
        """
        logger.debug("Initialising")
        uv_params = UV_params(level, params)
        self.uv_params = uv_params

        self.p = uv_params.p
        self.Fp = uv_params.Fp
        self.Fp2 = uv_params.Fp2
        self.E_start = uv_params.E_start

        self.K = uv_params.K
        self.w = uv_params.w
        self.max_order = uv_params.max_order
        self.order = uv_params.order

        def eval_frob(P):
            E = P.curve()
            frob = E.base_field().frobenius_endomorphism()
            return E([frob(c) for c in P.xy()])

        def eval_x_frob(xP):
            frob = xP.curve.base_field().frobenius_endomorphism()
            return xPoint(frob(xP.X), xP.curve)

        def eval_x_frob2(xP):
            frob = xP.curve.frobenius_endomorphism()
            return xPoint(frob(frob(xP.X)), xP.curve)

        self.pi = eval_frob
        self.pi_x = eval_x_frob
        self.pi2_x = eval_x_frob2

        for ell in self.uv_params.allowed_primes:
            if kronecker_symbol(-self.p, ell) == 1:
                self.lam = Integer(GF(ell)(-self.p).sqrt())
                self.small_ell = ell
                self.frak_ell = self.order*ell + self.order*(self.w-self.lam)
                assert self.frak_ell.norm() == ell
                break

    def action(self, E, frak_a, verb_return = False):
        r"""
        Function for computes the class action of any ideal.
        Input:
            - E: Elliptic curve E over Fp
            - frak_a: An ideal of End(E)
        Output:
            - frak_a * E
            - Optional: Extra data for measuring performance
        """
        compensate = 0
        logger.info("")
        logger.info("="*50)
        logger.info("Starting UV Solver....")
        logger.info("="*50)

        # Step 1: find u and v
        _t0 = time()
        _rerand_cnt = 0

        uv_output = None
        while not uv_output:
            uv_output = UVSolver(
                    [list(beta) for beta in frak_a.gens()], frak_a.norm(),
                    self.uv_params)
            if not uv_output:
                _rerand_cnt += 1
                frak_a = frak_a*self.frak_ell
                compensate += 1
                logger.debug("Rerandomising!")
                continue

        _t1 = time()
        logger.info(f"Used {_t1 - _t0} seconds")

        # Step 2: compute the kernel
        u, v, b_list, c_list, sig_1, sig_2, t = uv_output
        _n1, _cof1, _v21, _gens1 = b_list
        _n2, _cof2, _v22, _gens2 = c_list

        if compensate > 0:
            E = E.short_weierstrass_model()
            prev_j = None
            for _ in range(compensate):
                logger.debug("Compensation step")
                E_prev = E
                E = self.small_prime_ideal_action(E, self.small_ell, self.small_ell - self.lam, prev_j).codomain()
                prev_j = E_prev.j_invariant()
            E = E.montgomery_model()

        logger.info("")
        logger.info("="*50)
        logger.info("Starting Ideal Action....")
        logger.info("="*50)
        N1, N2, x_u, y_u, g_u, x_v, y_v, g_v, Pu, Qu, Pv, Qv = self.ideal_to_kernel(E, [u,v], [b_list, c_list], t)

        _t2 = time()
        logger.info(f"Used {_t2 - _t1} seconds")

        # Sanity checks
        assert all([not R*2**(t+2) for R in [Pu, Qu, Pv, Qv]])
        assert all([bool(R*2**(t+1)) for R in [Pu, Qu, Pv, Qv]])
        assert Pu.weil_pairing(Qu, 2**(t+2))**(g_v*N1*N2) == Pv.weil_pairing(Qv, 2**(t+2))**(g_u)
        assert N1*g_u*(x_u**2 + y_u**2) + N2*g_v*(x_v**2 + y_v**2) == 2**t

        # Step 3 : 4d iso
        logger.info("")
        logger.info("="*50)
        logger.info(f"Starting 4D isog")
        logger.info("="*50)
        logger.debug("Starting product: ")
        logger.debug(Pu.curve())
        logger.debug(Qv.curve())

        F = KaniClapotiIsog([Pu,Qu,Pv,Qv],[g_u,x_u,y_u,g_v,x_v,y_v,N1,N2,t])
        logger.debug(f"Change of coordinates: {F.iso_type} Twist={F.twist}")
        _t3 = time()

        logger.info(f"Used {_t3 - _t2} seconds")
        logger.info(f"Total time: {_t3 - _t0} seconds")

        stats = {
            "time_step1": _t1 - _t0,
            "time_step2": _t2 - _t1,
            "time_step3": _t3 - _t2,
            "time_total": _t3 - _t0,
            "rerandomizations": _rerand_cnt,
            "u": u,
            "v": v,
            "g_u": g_u,
            "x_u": x_u,
            "y_u": y_u,
            "g_v": g_v,
            "x_v": x_v,
            "y_v": y_v,
            "N1": N1,
            "N2": N2,
            "_cof1": _cof1,
            "_cof2": _cof2,
            "_v21": _v21,
            "_v22": _v22,
            "t": t
        }

        if verb_return:
            return F.Ea.change_ring(self.Fp), stats

        return F.Ea.change_ring(self.Fp)


    def ideal_to_kernel(self, E, uv, frak_bc, e):
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
            return sum([self.max_order * gi for gi in generators])

        frak_b = ideal_to_sage(frak_b_list, self.max_order)
        frak_c = ideal_to_sage(frak_c_list, self.max_order)

        if b_twopower > 1: # factors through two
            frak_b = div_by_n(frak_b, 2)
        if c_twopower > 1:
            frak_c = div_by_n(frak_c, 2)

        #The part thats the same for b and c
        frak_bc_same = frak_b + frak_c

        #We only want the "unique" part of b and c
        frak_b = div_by_n(frak_b * Conjugate(self.max_order, frak_bc_same), frak_bc_same.norm())
        frak_c = div_by_n(frak_c * Conjugate(self.max_order, frak_bc_same), frak_bc_same.norm())

        frak_b_2 = frak_b + self.max_order*((2**frak_b.norm().valuation(2)))
        frak_c_2 = frak_c + self.max_order*((2**frak_c.norm().valuation(2)))

        odd_part = Integer(frak_bc_same.norm()).prime_to_m_part(2)
        pow2_part = frak_bc_same.norm()/odd_part
        frak_bc_same_2 = frak_bc_same + self.max_order*pow2_part
        frak_bc_same_odd = frak_bc_same + self.max_order*odd_part

        b_cof = b_cof/odd_part
        c_cof = c_cof/odd_part

        #First step to the starting curve given by same direction steps
        ee = valuation(self.p+1, 2) - 1
        logger.debug("Walking to starting point (i.e. steps that are same for b and c)...")
        if frak_bc_same.norm() > 1:
            if pow2_part > 1:
                P_temp, Q_temp = self.TwoTorsBasis(E, ee)
                #P_temp, Q_temp = TwoTorsBasis(E, two_pow)
                _, E = self.ideal_above_2_action(E, P_temp, Q_temp, ee, frak_bc_same_2)
            if odd_part > 1:
                _, E = self.smooth_ideal_action(E, odd_part, frak_bc_same_odd)
        logger.debug(f"Done! Have used {time() - tstart}")
        P, Q = self.TwoTorsBasis(E, ee)
        #P, Q = TwoTorsBasis(E, ee)

        # Doing the Elkies steps for b
        logger.debug("odd Elkies steps for b and conjugate of c")
        # These are really different directions, so we can
        isogs_cof, E1 = self.smooth_ideal_action(E, b_cof*c_cof, frak_b*Conjugate(self.max_order, frak_c))
        P1, Q1 = P, Q
        for phi in isogs_cof:
            P1 = P1.push(phi)
            Q1 = Q1.push(phi)
        logger.debug(f"Done! Have used {time() - tstart}")

        #Adjust the order
        P1, Q1 = P1.xMUL(2**(ee - (e+2))), Q1.xMUL(2**(ee - (e+2)))

        assert P1.xMUL(2**(e+1))
        assert not P1.xMUL(2**(e+2))
        assert Q1.xMUL(2**(e+1))
        assert not Q1.xMUL(2**(e+2))
        #At this point, P1, Q1 are the "starting points"

        logger.debug("Doing g_u isogeny")
        g_u_isogs, Eu = self.random_smooth_isogeny(E1, g_u)
        Pu, Qu = P1, Q1
        for phi_gu in g_u_isogs:
            Pu = Pu.push(phi_gu)
            Qu = Qu.push(phi_gu)
        logger.debug(f"Done! Have used {time() - tstart}")

        # Now the second part
        # Evaluating the endomorphism
        logger.debug("Evaluating the endomorphism frak_b*frac_c_bar")
        rho_bc = PrincipalGenerator(frak_b*Conjugate(self.max_order, frak_c))
        rho_bc_val_2 = rho_bc.norm().valuation(2)

        #assert (frak_b*frak_c).norm() == rho_bc.norm()
        #assert rho_bc in frak_b
        #assert rho_bc.conjugate() in frak_c

        #adjust_2 = 0
        #if rho_bc_val_2 > 0: #Here we really need to clear denominators to evaluate the endomorphism
        #    rho_bc *= 2
        #    adjust_2 = 1

        #print("!!!!!!!!!")
        #print(ee-(rho_bc_val_2))
        #print(e+2)
        #assert ee-(rho_bc_val_2 + adjust_2) >= e+2, "NotImplemented: We can drop adjust_2 if we allow extension-fields"
        #print(f"Evaluating {rho_bc}")
        max = (ee - rho_bc_val_2 == e+2)
        Pc = self.eval_endomorphism(rho_bc, P, False, max)
        Qc = self.eval_endomorphism(rho_bc, Q, True, max)

        logger.debug(f"Done! Have used {time() - tstart}")

        # Act with 2-ideals right away.
        logger.debug("Acting with 2-part of frak_b and frak_c")
        #print(f"!!!! {frak_c_2.norm(), frak_b_2.norm()}")
        #print(f"!!!? {rho_bc_val_2}")
        if frak_c_2.norm() > 1 or frak_b_2.norm() > 1:
            frak_bc_2 = frak_c_2*Conjugate(self.max_order, frak_b_2)
            phi_2, E = self.ideal_above_2_action(E, P, Q, ee, frak_bc_2)
            Pc = Pc.push(phi_2)
            Qc = Qc.push(phi_2)

        logger.debug(f"Done! Have used {time() - tstart}")
        P2 = Pc.xMUL(2**(ee - (e+2) - (rho_bc_val_2)))
        Q2 = Qc.xMUL(2**(ee - (e+2) - (rho_bc_val_2)))

        assert P2.xMUL(2**(e+1))
        assert not P2.xMUL(2**(e+2))
        assert Q2.xMUL(2**(e+1))
        assert not Q2.xMUL(2**(e+2))
        E2 = E

        logger.debug("Doing g_v isogeny")
        Pv, Qv = P2, Q2

        g_v_isogs, Ev = self.random_smooth_isogeny(E2, g_v)

        for phi_gv in g_v_isogs:
            Pv = Pv.push(phi_gv)
            Qv = Qv.push(phi_gv)
        logger.debug(f"Done! Have used {time() - tstart}")

        Ev = Pv.curve.change_ring(self.Fp2)
        Eu = Pu.curve.change_ring(self.Fp2)

        Pv = Ev.lift_x(Pv.X)
        Qv = Ev.lift_x(Qv.X)

        Pu = Eu.lift_x(Pu.X)
        Qu = Eu.lift_x(Qu.X)

        #Make sure to lift correctly
        if Pu.weil_pairing(Qu, 2**(e+2))**(g_v*N1*N2) != Pv.weil_pairing(Qv, 2**(e+2))**(g_u):
            Qv = -Qv

        # Double check full order
        assert Pu.weil_pairing(Qu, 2**(e+2))**(2**(e+2)) == 1
        assert Pu.weil_pairing(Qu, 2**(e+2))**(2**(e+1)) != 1
        assert Pv.weil_pairing(Qv, 2**(e+2))**(2**(e+2)) == 1
        assert Pv.weil_pairing(Qv, 2**(e+2))**(2**(e+1)) != 1
        # Weil pairing check
        assert Pu.weil_pairing(Qu, 2**(e+2))**(g_v*N1*N2) == Pv.weil_pairing(Qv, 2**(e+2))**(g_u)

        return N1, N2, x_u, y_u, g_u, x_v, y_v, g_v, Pu, Qu, Pv, Qv


    def sample_ideal(self):
        r"""
        Samples a random ideal of End(E)
        """
        N_a, frak_a = RandomDegreeOnePrimeIdeal(self.order, self.w, self.p)
        return frak_a

    def ideal_above_2_action(self, E, P, Q, basis_ord, frak_a):
        r"""
        Computes the action of a power of an ideal above 2 using Velu.
        Input:
            - E: Elliptic curve over Fp
            - P, Q: Basis of E[2^basis_ord], where P is in E(Fp) and Q is in Et(Fp) for a quadratic twist Et
            - basis_ord: The order of P, Q (log_2)
            - frak_a: A power of an ideal above 2
        Output:
            - phi : The isogeny corresponding to frak_a
            - frak_a * E
        """
        _, val_2_frak_a = factor(frak_a.norm())[0]

        if (self.w - 1)/2 in frak_a:
            K = P
        else:
            assert (self.w + 1)/2 in frak_a
            K = Q

        K = K.xMUL(2**(basis_ord - val_2_frak_a))
        assert K.xMUL(2**(val_2_frak_a - 1))
        assert not K.xMUL(2**val_2_frak_a)

        phi_2 = K.xMUL(2**(val_2_frak_a-1)).xISOG(E, 2)
        phi = phi_2
        for steps in range(1, val_2_frak_a):
            E = phi.codomain()
            K = K.push(phi_2)
            assert K.xMUL(2**(val_2_frak_a-1-steps))
            assert not K.xMUL(2**(val_2_frak_a-steps))
            phi_2 = K.xMUL(2**(val_2_frak_a-1-steps)).xISOG(E, 2)
            phi = phi_2 * phi

        E_out = phi.codomain().montgomery_model()

        return phi.codomain().isomorphism_to(E_out) * phi, E_out


    def small_prime_ideal_action(self, E, ell, lam=None, prev_j=None):
        r"""
        Computes the action of a an ideal above ell using Elkies algorithm.
        Input:
            - E: Elliptic curve over Fp IN WEIERSTRASS FORM
            - ell: A small prime, split in End(E)
            - lam (optional): The eigenvalue determining the ideal frak_l = (ell, pi - lam)
            - prev_j (optional): The j-invariant of the WRONG ell-isogenous curve
        Output:
            - phi : The isogeny corresponding to an ideal above ell
        """
        # Fix later, but Elkies uses short Weierstrass
        # Change ring stuff is messy, but needed to actually compute the rational model
        assert isWeierstrass(E)
        h = Elkies(E, ell, self.Fp, lam=lam, prev_j=prev_j)
        phi = E.isogeny(h)
        return phi


    def small_prime_degree_isogeny(self, E, ell):
        r"""
        Computes an (unspecified) isogeny of degree ell
        Input:
            - E: Elliptic curve over Fp
            - ell: A small prime, not necessarily split in End(E)
        Output:
            - phi : An isogeny from E of degree ell
        """
        # Just a random isogeny of degree ell
        # Not (necessarily) oriented (if prime is Elkies, we should use the other)
        Fbig, k, ontwist = self.ell_to_Fpk[ell]
        Ebig = E.change_ring(Fbig)
        if ontwist:
            Ebig_ord = self.p**k + (-1)**k
        else:
            Ebig_ord = self.p**k - (-1)**k
        assert not Ebig_ord % ell, f"ehh {(self.p**k + (-1)**k) % ell}, {(self.p**k - (-1)**k) % ell}"

        cof = Ebig_ord//ell

        # Picking point smarter makes this easier
        K = None
        while not K:
            K = random_xPoint(Ebig, Fbig, ontwist)
            # TODO: Clear some cofactors by working in trace 0 part
            K = K.xMUL(cof)

        phi = K.xISOG(E, ell)
        return phi


    def _adjust_2_tors_order(self, P, Q, e, return_scalars = False):
        r"""
        Multiplies P, Q by suitable cofactors so they form a basis of E[2**e]
        Input:
            - P, Q, points of order at least 2**e, so that E[2**e] \subseteq <P, Q>
            - e: Exponent of 2.
        Output:
            - Pm, Qm: Basis of E[2**e]
            - Optional: The scalars a, b so that Pm = a*P, Qm = b*Q
        """
        assert P.xMUL((2**(e-1)))
        assert Q.xMUL((2**(e-1)))
        P_small = P.xMUL(2**e)
        scal_1 = 1
        scal_2 = 1
        while P_small:
            P = P.xMUL(2)
            P_small = P_small.xMUL(2)
            scal_1 *= 2
        Q_small = Q.xMUL(2**e)
        while Q_small:
            Q = Q.xMUL(2)
            Q_small = Q_small.xMUL(2)
            scal_2 *= 2
        if return_scalars:
            return P, Q, scal_1, scal_2
        return P, Q


    def smooth_ideal_action(self, E, cof, frak):
        r"""
        Computes the action of a an ideal of odd, smooth norm using Elkies algorithm repeatedly.
        Input:
            - E: Elliptic curve over Fp
            - cof: norm of frak
            - frak: An ideal of End(E)
        Output:
            - phi : The isogeny corresponding to frak
            - frak * E
        """
        # Evaluate a smooth, odd normed ideal
        isogs = []
        E1 = E.short_weierstrass_model()
        pre_isom = E.isomorphism_to(E1)
        isogs.append(pre_isom)
        if cof > 1:
            for ell_cof, e_cof in factor(cof):
                lam = Integer(GF(ell_cof)(-self.p).sqrt())
                if not self.w - lam in frak + ell_cof*self.max_order:
                    # In this case, we want the other one
                    lam = ell_cof - lam
                assert self.w - lam in frak + ell_cof*self.max_order
                prev_j = None
                for _ in range(e_cof):
                    logger.debug(f'Starting Elkies step of degree {ell_cof}')
                    phi_part = self.small_prime_ideal_action(E1, ell_cof, lam, prev_j)
                    prev_j = E1.j_invariant()
                    E1 = phi_part.codomain()
                    isogs.append(phi_part)

        E_out = isogs[-1].codomain().montgomery_model()
        post_isom = isogs[-1].codomain().isomorphism_to(E_out)
        isogs.append(post_isom)
        return isogs, E_out


    def random_smooth_isogeny(self, E, g):
        r"""
        Returns a random isogeny of smooth degree g from E
        Input:
            - E: Elliptic curve over Fp
            - g: Smooth number
        Output:
            - phi : An isogeny from E of degree g
            - E_out : phi(E)
        """
        # Random smooth odd norm isogeny
        isogs = []
        if g > 1:
            E1 = E.short_weierstrass_model()
            pre_isom = E.isomorphism_to(E1)
            isogs.append(pre_isom)
            logger.debug(f"Extra isogeny for g = {factor(g)}")
            # Do elkies first, to stay over Fp
            g_elkies = []
            g_non_elkies = []
            for ell, e in factor(g):
                if kronecker_symbol(-self.p, ell) == 1:
                    g_elkies.append((ell, e))
                else:
                    g_elkies.append((ell, e))
            for ell, e in g_elkies:
                prev_j = None
                for _ in range(e):
                    phi_ell = self.small_prime_ideal_action(E1, ell, prev_j = prev_j)
                    prev_j = E1.j_invariant()
                    E1 = phi_ell.codomain()
                    isogs.append(phi_ell)

            E_out = isogs[-1].codomain().montgomery_model()
            post_isom = isogs[-1].codomain().isomorphism_to(E_out)
            isogs.append(post_isom)

            for ell, e in g_non_elkies:
                for _ in range(e):
                    phi_ell = self.small_prime_degree_isogeny(E_out, ell)
                    E_out = phi_ell.codomain()
                    isogs.append(phi_ell)

            return isogs, E_out
        else:
            return [], E


    def eval_endomorphism(self, rho, P, twist, max):
        r"""
        Evaluates an element of End(E) on a point P.
        P is assumed to be Fp rational on either E or on a twist
        """
        a, b = list(rho) #write as a + b*pi
        n = rho.denominator()
        if twist:
            #a + b*pi = a - b, since pi(P) = -P
            m = Integer(a - b)
        else:
            #a + b*pi = a + b, since pi(P) = P
            m = Integer(a + b)
        if max and (n == 2): #See Appendix D.2
            mP = P.xMUL(m)
            E = P.curve

            T0 = self.find_Ts(E, only_T0 = True)

            #Already know which point to choose, Remark D.2
            assert (P.X - T0.x()).is_square() != twist

            T = xPoint(T0.x(), E)

            return translate_by_T(mP, T)

        else:
            return P.xMUL(m)


    def TwoTorsBasis(self, E, e):
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

        T0, Tm1, T1 = self.find_Ts(E)

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

        assert (self.p+1) % 2**(e+1) == 0
        cofac = (self.p+1)/2**(e+1)
        P = P.xMUL(cofac)
        Q = Q.xMUL(cofac)

        assert P.xMUL(2**(e-1))
        assert not P.xMUL(2**e)
        assert Q.xMUL(2**(e-1))
        assert not Q.xMUL(2**e)
        assert Q.xMUL(2**(e-1)) != P.xMUL(2**(e-1))

        logger.debug(f"    > Total time for 2-tors basis finding: {time()-tstart}")

        return P, Q

    def find_Ts(self, E, only_T0 = False):
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
