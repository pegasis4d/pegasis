from time import time
import logging

from sage.rings.integer import Integer
from sage.arith.misc import kronecker_symbol
from sage.rings.finite_rings.finite_field_constructor import GF

from xonly import xPoint
from lcsidh import UVSolver
from uv_params import UV_params
from ideals import random_degree_one_ideal
from dim1_isogenies import small_prime_ideal_action
from ideals import ideal_to_sage
from ideal_to_kernel import ideal_to_kernel

#from theta_lib.isogenies.Kani_clapoti import KaniClapotiIsog
from Theta_dim4.Theta_dim4_sage.pkg.isogenies.Kani_clapoti import KaniClapotiIsog

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter("%(name)s [%(levelname)s] %(message)s")
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)


class PEGASIS:
    def __init__(self, level, params=None):
        r"""
        Main Clapoti class. Takes as input the security level
        (100|500|1000|1500|2000|4000) and generates the necessary parameters and
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
                self.frak_ell = self.order * ell + self.order * (self.w - self.lam)
                assert self.frak_ell.norm() == ell
                break

    def action(self, E, frak_a, verb_return=False):
        r"""
        Function for computes the class action of any ideal.
        Input:
            - E: Elliptic curve E over Fp
            - frak_a: An ideal of End(E)
        Output:
            - frak_a * E
            - Optional: Extra data for measuring performance
        """
        rerandomizations = 0
        logger.info("")
        logger.info("=" * 50)
        logger.info("Starting UV Solver....")
        logger.info("=" * 50)

        # Step 1: find u and v
        _t0 = time()

        uv_output = None
        while not uv_output:
            uv_output = UVSolver(frak_a, frak_a.norm(), self.uv_params)
            if not uv_output:
                rerandomizations += 1
                frak_a = frak_a * self.frak_ell
                logger.debug("Rerandomising!")

        _t1 = time()
        logger.info(f"Used {_t1 - _t0} seconds")

        # Step 2: compute the kernel
        u, v, b_list, c_list, sig_1, sig_2, t = uv_output
        _n1, _cof1, _v21, _gens1 = b_list
        _n2, _cof2, _v22, _gens2 = c_list

        N1, b_cof, b_twopower, frak_b_list = b_list
        N2, c_cof, c_twopower, frak_c_list = c_list
        frak_b = ideal_to_sage(frak_b_list, self.order)
        frak_c = ideal_to_sage(frak_c_list, self.order)
        assert frak_b.is_equivalent(frak_a)
        assert frak_c.is_equivalent(frak_a)

        E = E.short_weierstrass_model()
        prev_j = None
        for _ in range(rerandomizations):
            logger.debug("Compensation step")
            E_prev = E
            E = small_prime_ideal_action(E, self.small_ell, self.small_ell - self.lam, prev_j).codomain()
            prev_j = E_prev.j_invariant()
        E = E.montgomery_model()

        logger.info("")
        logger.info("=" * 50)
        logger.info("Starting Ideal Action....")
        logger.info("=" * 50)
        N1, N2, x_u, y_u, g_u, x_v, y_v, g_v, Pu, Qu, Pv, Qv = ideal_to_kernel(
            E, [u, v], [b_list, c_list], self.max_order, self.w, t
        )

        _t2 = time()
        logger.info(f"Used {_t2 - _t1} seconds")

        # Sanity checks
        assert all([not R * 2 ** (t + 2) for R in [Pu, Qu, Pv, Qv]])
        assert all([bool(R * 2 ** (t + 1)) for R in [Pu, Qu, Pv, Qv]])
        assert Pu.weil_pairing(Qu, 2 ** (t + 2)) ** (g_v * N1 * N2) == Pv.weil_pairing(Qv, 2 ** (t + 2)) ** (g_u)
        assert N1 * g_u * (x_u**2 + y_u**2) + N2 * g_v * (x_v**2 + y_v**2) == 2**t

        # Step 3 : 4d iso
        logger.info("")
        logger.info("=" * 50)
        logger.info(f"Starting 4D isog")
        logger.info("=" * 50)
        logger.debug("Starting product: ")
        logger.debug(Pu.curve())
        logger.debug(Qv.curve())

        F = KaniClapotiIsog([Pu, Qu, Pv, Qv], [g_u, x_u, y_u, g_v, x_v, y_v, N1, N2, t])
        logger.debug(f"Change of coordinates: {F.iso_type} Twist={F.twist}")
        _t3 = time()

        logger.info(f"Used {_t3 - _t2} seconds")
        logger.info(f"Total time: {_t3 - _t0} seconds")

        stats = {
            "time_step1": _t1 - _t0,
            "time_step2": _t2 - _t1,
            "time_step3": _t3 - _t2,
            "time_total": _t3 - _t0,
            "rerandomizations": rerandomizations,
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
            "t": t,
        }

        if verb_return:
            return F.Ea.change_ring(self.Fp), stats

        return F.Ea.change_ring(self.Fp)

    def sample_ideal(self):
        r"""Sample a random ideal of the class group
        """
        _, ideal = random_degree_one_ideal(self.order, self.w, self.p)
        return ideal
