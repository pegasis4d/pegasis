
### Directly copied from ApresSQI

#######################################################
#                                                     #
#         extensive usage from DeuringFTP             #
#  https://github.com/friends-of-quaternions/deuring  #
#                                                     #
#######################################################

"""
MIT Licence

Copyright (c) 2023 J. K. Eriksen, L. Panny, J. Sotáková, M. Veroni

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.
"""

from sage.all import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism

class xPoint:
    r"""
    Class for x-only arithmetic on Montgomery curves.
    """
    def __init__(self, X, E):
        #assert isMontgomery(E) This assertion is not true anymore, because we use intermidiate weierstrass curves 
        k = E.base_field()
        self.X = k(X) if X is not None else None
        self.curve = E

    def __repr__(self):
        return f"Point with X-coordinate {self.X} on {self.curve}"

    def __bool__(self):
        return self.X is not None
    
    def __eq__(self, other):
        assert isinstance(other, xPoint)
        return self.X == other.X and self.curve == other.curve
    
    def on_twist(self):
        y_sqr = self.X**3 + MontgomeryA(self.curve)*self.X**2 + self.X
        return not is_square(y_sqr)
    
    def xDBL(self):
        if not self:
            return self
        a = MontgomeryA(self.curve)
        X = self.X
        XX = X**2
        Z3 = 4*X*(XX + a*X + 1)
        if Z3 == 0:
            return xPoint(None, self.curve)
        X3 = (XX - 1)**2
        return xPoint(X3/Z3, self.curve)
    
    def eval_isomorphism_x(self, isomorphism, x=None):
        #return isomorphism.x_rational_map()(self.X)
        if x is None:
            x=self.X
        return (x - isomorphism.r) / isomorphism.u**2
    
    def push(self, isogeny):
        r"""
        Given an isogeny phi (where phi is computed as a composition of isogenies
        computed with Kohel's algorithm, or a composition of 2 and 4 isogenies), 
        returns phi(self)
        """

        newX = self.X

        if type(isogeny) is WeierstrassIsomorphism:
            newX = self.eval_isomorphism_x(isogeny, newX)
        else:

            if type(isogeny) is not EllipticCurveHom_composite:
                isogeny = EllipticCurveHom_composite.from_factors([isogeny])

            for isopart in isogeny.factors():
                #assert isinstance(isopart, EllipticCurveIsogeny)
                if isopart._EllipticCurveIsogeny__algorithm == 'velu_sqrt':
                    if (isom := isopart._pre_iso):
                        newX = self.eval_isomorphism_x(isom, newX)
                    
                    newX = isopart._raw_eval(newX)
                    if not newX:
                        verbose("Point seems to be in the kernel")
                        return xPoint(None, isogeny.codomain())
                    
                    if (isom := isopart._post_iso):
                        newX = self.eval_isomorphism_x(isom, newX)

                elif isopart._EllipticCurveIsogeny__algorithm == 'kohel':

                    if (isom := isopart._EllipticCurveIsogeny__pre_isomorphism):
                        newX = self.eval_isomorphism_x(isom, newX)

                    phi = isopart._EllipticCurveIsogeny__phi
                    psi = isopart._EllipticCurveIsogeny__psi
                    try:
                        newX = phi(newX) / psi(newX)**2
                    except ZeroDivisionError:
                        verbose("Point seems to be in the kernel")
                        newX = None
                        return xPoint(None, isogeny.codomain())

                    if (isom := isopart._EllipticCurveIsogeny__post_isomorphism):
                        newX = self.eval_isomorphism_x(isom, newX)

                elif 4 % isopart.degree() == 0:
                    try:
                        newX = isopart.x_rational_map()(newX)
                    except ZeroDivisionError:
                        verbose("Point seems to be in the kernel")
                        return xPoint(None, isogeny.codomain())
                else:
                    print(isopart)
                    print(isopart.degree())
                    assert False, "Cannot evaluate in this type of isogeny"

        new_curve = isogeny.codomain().base_extend(self.curve.base_field())
        return xPoint(newX, new_curve)

    def xMUL(self, n):
        """
        Given an integer n, computes [n]self
        """
        n = ZZ(n).abs()
        if n == 0:
            return xPoint(None, self.curve)
        if n == 1:
            return self
        R0, R1, diff = self, self.xDBL(), self
        for i in [int(b) for b in bin(n)[3:]]:
            R0pR1 = xADD(R0, R1, diff)
            diff = xADD(R0, R1, R0pR1)
            if i == 0:
                R0, R1 = R0.xDBL(), R0pR1
            if i == 1:
                R0, R1 = R0pR1, R1.xDBL()
        return R0
    
    def x_multiples(self):
        r"""
        Given an xPoint P, 
        returns all the roots of the kernel polynomial of <P>.
        """
        if not self:
            return []
        Ps = [self]
        P = self.xDBL()
        if not P:
            xs = [P.X for P in Ps]
            return xs
        while P not in Ps[-2:]:
            Ps.append(P)
            P = xADD(Ps[-1], Ps[0], Ps[-2])
        xs = [P.X for P in Ps]
        return xs
    
    def kernel_polynomial(self, E, l):
        r"""
        Given that self has order l,
        returns the kernel polynomial of <P>
        """
        x = self.X
        F2 = E.base_field()
        Fbig = x.parent()

        ext = Fbig.over(F2)

        R,X = F2['X'].objgen()
        assert E.base_extend(Fbig) == self.curve

        assert l.is_prime()
        if l <= 3:
            return R([-x, 1])

        try:
            X_in_F2 = F2(x)
        except ValueError:
            pass
        else:
            return prod(X - xx for xx in xPoint(X_in_F2, E).x_multiples())

        #if E.frobenius() not in ZZ:
        #    raise NotImplementedError

        def minpoly(elt):
            return ntl_funs(F2).minpoly_mod(ext(elt).vector(), ext.modulus())

        fs = [minpoly(x)]

        k = fs[0].degree()
        m = (l-1) // (2*k)

        assert k > 1    # handled above

        from sage.schemes.elliptic_curves.isogeny_small_degree import _least_semi_primitive
        a = _least_semi_primitive(l)

        xi = xPoint(x, E.change_ring(Fbig))
        for _ in range(1, m):
            xi = xi.xMUL(a)
            fs.append(minpoly(xi.X))

        return prod(fs)
    
    def xISOG(self, E, l):
        r"""
        Given a x = x(P) of a point P on E of order l,
        computes the separable isogeny with <P> as kernel
        """
        l = Integer(l)
        h = self.kernel_polynomial(E, l)
        return E.isogeny(h, model='montgomery')

        

def xADD(P, Q, PmQ):
    r"""
    Given a x(P), x(Q) and x(P-Q),
    computes x(P+Q)
    """
    if not P:
        return Q
    if not Q:
        return P
    if not PmQ:
        return P.xDBL()
    A = P.X + 1
    B = P.X - 1
    C = Q.X + 1
    D = Q.X - 1
    DA = D*A
    CB = C*B
    X5 = (DA + CB)**2
    Z5 = PmQ.X*(DA - CB)**2
    if Z5 == 0:
        return xPoint(None, P.curve)
    return xPoint(X5/Z5, P.curve)

def xDBLMUL(m, n, P, Q, PmQ):
    r"""
    Algorithm 9 from eprint2017/212
    Given a x(P), x(Q) and x(P-Q), and scalars m,n
    computes x([m]P+[n]Q)
    """
    if m == 0:
        return Q.xMUL(n)
    if n == 0:
        return P.xMUL(m)
    s0, s1 = Integer(m), Integer(n)
    x0, x1, xmin = P, Q, PmQ
     
    while s0 != 0:

        if s1 < s0:
            s0, s1 = s1, s0
            x0, x1 = x1, x0

        if s1 <= 4*s0:
            s1 = s1 - s0
            tmp = xADD(x1, x0, xmin)
            xmin = x0
            x0 = tmp
        elif s0 % 2 == s1 % 2:
            s1 = (s1-s0)//2
            x0 = xADD(x1, x0, xmin)
            x1 = x1.xDBL()
        elif s1 % 2 == 0:
            s1 = s1//2
            xmin = xADD(x1, xmin, x0)
            x1 = x1.xDBL()
        else:
            s0 = s0//2
            xmin = xADD(x0, xmin, x1)
            x0 = x0.xDBL()

    while s1 % 2 == 0:
        s1 = s1//2
        x1 = x1.xDBL()

	#extra step with tripling as described in eprint2017/212
    while s1 % 3 == 0:
        s1 = s1//3
        x1 = x1.xMUL(3)

    if s1 > 1:
        x1 = x1.xMUL(s1)

    return x1

################################################################

import ast
class _ntl_funs:
    r"""
    An object encapsulating the NTL context for a given finite
    field F, such that polynomials in F[X] can be converted to
    NTL using .ntlify() and minimal polynomials in F[X]/f can
    be computed using .minpoly_mod().
    """
    def __init__(self, F):
        self.ctx = ntl.ZZ_pEContext(ntl.ZZ_pX(F.modulus().list(), F.characteristic()))
        self.F = F
        self.R = F['X']

    def ntlify(self, poly):
        try:
            poly = poly.vector()
        except AttributeError:
            pass
        assert poly.base_ring() == self.F
        return ntl.ZZ_pEX([ntl.ZZ_pE(c.list(), self.ctx) for c in poly])

    def minpoly_mod(self, elt, mod):
        ntl_mod = self.ntlify(mod)
        ntl_elt = self.ntlify(elt) % ntl_mod
        ntl_res = ntl_elt.minpoly_mod(ntl_mod)
        return self.R(ast.literal_eval(str(ntl_res).replace(' ',',')))

@cached_function
def ntl_funs(F2):
    r"""
    Caching helper function to construct the _ntl_funs object
    for a given base field F2 of degree 2 over its prime field.
    """
    assert F2.degree() == 2
    return _ntl_funs(F2)

def isMontgomery(E):
    r"""
    Given a curve E,
    returns whether E is in Montgomery form
    """
    a1, a2, a3, a4, a6 = E.a_invariants()
    return a1 == 0 and a3 == 0 and a6 == 0 and a4 == 1

def isWeierstrass(E):
    r"""
    Given a curve E,
    returns whether E is in Montgomery form
    """
    a1, a2, a3, a4, a6 = E.a_invariants()
    return a1 == 0 and a2 == 0 and a3 == 0

def MontgomeryA(E):
    r"""
    Given a Montgomery curve E_A,
    returns the coefficient A
    """
    assert isMontgomery(E)
    return E.a2()

def random_xPoint(E, F, ontwist):
    P = xPoint(F.random_element(), E)
    while P.on_twist() != ontwist:
        P = xPoint(F.random_element(), E)
    return P

def PointDiff(P, Q):
    E = P.curve
    assert isMontgomery(E)
    A = E.a2()
    Px, Qx = P.X, Q.X

    if Px == Qx:
        return xPoint(None, E)
        
    PmQZ = Px - Qx
    t2 = Px*Qx
    t3 = t2 - 1
    t0 = PmQZ*t3
    PmQZ = PmQZ**2
    t0 = t0**2
    t1 = t2 + 1
    t3 = Px + Qx
    t1 = t1*t3
    t2 = t2*A
    t2 = 2*t2
    t1 = t1 + t2
    t2 = t1**2
    t0 = t2-t0
    assert t0.is_square()
    t0 = t0.sqrt()
    PmQX = t0 + t1

    return xPoint(PmQX/PmQZ, E)


### Not needed?
def translate_by_T(P, T):
    assert P.curve == T.curve
    A = T.X
    X = P.X

    #projectively (X, Z) + (A, B) -> (AX - BZ, BX - AZ)
    if A == 0:
        return xPoint(1/X, P.curve)
    else:
        return xPoint((A*X - 1)/(X - A), P.curve)
