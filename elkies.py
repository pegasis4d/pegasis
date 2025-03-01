from sage.all import *
import random
from xonly import xPoint

def Elkies(E, ell, Fp, lam = None, prev_j = None):
    """
    Elkies algorithm for computing isogenies as detailed in [BSS - Elliptic
    Curves in Cryptography, Chapter VII].
    Input:
        - E: starting curve
        - ell: a split prime
        - Fp: base field
        - lam: if looking for a specific Elkies isogeny, optionally provide the
          corresponding eigenvalue
        - prev_j: if the curve is the codomain of another ell-isogeny,
          optionally provide the j invariant of the previous curve to avoid
          backtracking
    Output:
        - the kernel polynomial hc
    """
    #TODO maybe we want montomery curves later
    p = Fp.characteristic()
    assert ell > 2
    d = (ell-1)/2
    j = Fp(E.j_invariant())
    _, _, _, a, b = [Fp(c) for c in E.a_invariants()]

    R = Fp['X,Y']
    X, Y = R._first_ngens(2)

    E4 = -48*a
    E6 = 864*b

    jp = -E6*j/E4

    phi_ell = R(classical_modular_polynomial(ell))
    X, Y = phi_ell.variables()
    phi_ellx = derivative(phi_ell, X)
    phi_ellxx = derivative(phi_ellx, X)
    phi_ellxy = derivative(phi_ellx, Y)

    phi_elly = derivative(phi_ell, Y)
    phi_ellyy = derivative(phi_elly, Y)

    phi_ell_eval = phi_ell(X, j).univariate_polynomial()
    if prev_j:
        x = phi_ell_eval.variables()[0]
        phi_ell_eval = phi_ell_eval // (x - prev_j)
    # print(f"Finding j: {time()-tstart}")
    if prev_j:
        j_ells = [phi_ell_eval.any_root()] #There is only one root in this case!
        #assert len(j_ells) == 1, f"j_ells was {j_ells}, and prev_j was {prev_j}"
    else:
        j_ells = [r for r,e in phi_ell_eval.roots()]
        assert len(j_ells) == 2, f"j_ells was {j_ells}"

    #print(j_ells)
    jt = j_ells[0]

    assert jt != 0 and jt != 1728, "encountered E0"

    hc = _derive_hc(Fp, a, b, ell, d, j, jt, jp, E4, E6, phi_ellx, phi_ellxx,
                    phi_ellxy, phi_elly, phi_ellyy)

    if lam and not prev_j:
        if not CheckElkies(E, ell, hc, lam):
            #Do the other one
            assert j_ells[1] != 0 and j_ells[1] != 1728, "encountered E0"
            hc = _derive_hc(Fp, a, b, ell, d, j, j_ells[1], jp, E4, E6, phi_ellx, phi_ellxx, phi_ellxy, phi_elly, phi_ellyy)
    return hc

def WeierstrassP(a,b,d):
    # Weierstrass p up to degree d in w=z^2
    coefs = [-a/5, -b/7]
    for k in range(2, d):
        ck = 0
        for j in range (k-1):
           ck = ck + coefs[j]*coefs[k-2-j]
        ck = ck*3/((k-1)*(2*(k+1)+3))
        coefs.append(ck)
    return coefs


def _derive_hc(Fp, a, b, ell, d, j, jt, jp, E4, E6, phi_ellx, phi_ellxx, phi_ellxy, phi_elly, phi_ellyy):
    jtp = -jp*phi_ellx(j, jt)/(ell*phi_elly(j, jt))
    at = -jtp**2/(48*jt*(jt-1728))
    bt = -jtp**3/(864*jt**2*(jt-1728))

    E4t = -48*at
    E6t = 864*bt

    fracjs = -(jp**2*phi_ellxx(j, jt) + 2*ell*jp*jtp*phi_ellxy(j, jt) + (ell*jtp)**2*phi_ellyy(j, jt))/ (jp*phi_ellx(j, jt))
    p1 = ell*fracjs/2 + (Fp(ell)/4)*(E4**2/E6 - ell*E4t**2/E6t) + (Fp(ell)/3)*(E6/E4 - ell*E6t/E4t)

    c = WeierstrassP(a,b,d)
    ckst = WeierstrassP(ell**4*at, ell**6*bt, d)
    ct = [ckst[i] - ell*c[i] for i in range(len(c))] # difference in formula VII.23

    # Computing the coefficients of VII.23 and store as A[i]
    Fpw, w = PowerSeriesRing(Fp, 'w', default_prec=d+1).objgen()
    evp = -(p1/2)*w - sum(w**(i+1)*ct[i-1] / ((2*i+1)*(2*i+2)) for i in range(1,d))
    exp_evp = evp.exp(d+1)
    A = exp_evp.coefficients()
    C = sum(c[i-1]*w**i for i in range(1, d+1))

    # Computing all powers of C starting with zeroth power
    Cpow = [Fpw(1), C];
    for i in range(2, d+1):
        Cpow.append(C*Cpow[-1])

    # Now doing recurrence relation VII.24
    hc = [1, -p1/2]
    for i in range(2, d+1):
        newcoeff = A[i]
        for k in range(1, i+1):
            insum = 0
            for j in range(k+1):
                insum += Fp(binomial(d-i+k, k-j))*Cpow[k-j][j]
            newcoeff -= insum*hc[i-k]
        hc.append(newcoeff)

    Rx, x = PolynomialRing(Fp, 'x').objgen()
    hc = Rx(hc[::-1])
    return hc

def CheckElkies(E, ell, h, lam):
    p = E.base_field().characteristic()
    if kronecker_symbol(-p, ell) != 1:
        assert False, "not Elkies"

    _, _, _, a, b = E.a_invariants()
    f = h.parent()([b, a, 0, 1])

    B = pow(f, (p-1)/2, h)
    check_wrong = False
    if lam > ell//2:
        lam = ell-lam
        check_wrong = True

    if lam == 1:
        if check_wrong:
            return B != 1
        else:
            return B == 1

    # Stupid way for now, no sage function for directly computing mod h
    RR = h.parent().quotient_ring(h)
    y_coord = E.multiplication_by_m(lam)[1]
    y_coord_num = y_coord.numerator()
    x, y = y_coord_num.variables()
    y_coord_den = y_coord.denominator()

    if check_wrong:
        return RR(B)*RR(y_coord_den) != RR(y_coord_num(x, 1).univariate_polynomial())
    else:
        return RR(B)*RR(y_coord_den) == RR(y_coord_num(x, 1).univariate_polynomial())


def NonElkiesIsogeny(E, ell, Fp2):
    # For u, v, when not elkies prime.
    # Make that either ell divides p+1, the other way works...

    E = E.change_ring(Fp2)
    #print(factor(E.division_polynomial(ell)))
    # Make sure these are the only cases...
    h_ell = factor(E.division_polynomial(ell))[0][0]

    m = (ell-1)//(2*h_ell.degree())
    if m == 1:
        return h_ell

    #Deuring for the people <3
    from sage.schemes.elliptic_curves.isogeny_small_degree import _least_semi_primitive
    a = _least_semi_primitive(ell)
    fs = [h_ell]
    Fbig = Fp2.extension(h_ell)
    x = Fbig.gens()[0]
    xi = xPoint(x, E.change_ring(Fbig))
    for _ in range(1, m):
        xi = xi.mul(a)
        fs.append(xi.X.minpoly())

    h_ell = prod(fs)
    assert h_ell.degree() == (ell-1)/2, f"Degree : {h_ell.degree()}, shouldve been {(ell-1)/2}"
    return h_ell



### This is way slower for now :((
def ElkiesDirect(E, ell, lam = None):
    psi_x = E.division_polynomial(ell)
    Fp = E.base_field()
    p = Fp.characteristic()
    if not lam:
        F_ell = GF(ell)
        lam = F_ell(-p).sqrt()

    R = psi_x.parent()
    x = psi_x.variables()[0]
    RR = R.quotient_ring(psi_x)

    X = RR._first_ngens(1)[0]

    check_wrong = False
    if lam > ell//2:
        lam = ell-lam
        check_wrong = True

    if lam == 1:
        h_sq_bar = X**p - X
    else:
        mult_by_lam = E.multiplication_by_m(lam)[0]
        mult_by_lam_num = mult_by_lam.numerator().univariate_polynomial()
        mult_by_lam_denom = mult_by_lam.denominator().univariate_polynomial()
        h_sq_bar = X**p*mult_by_lam_denom(X) - mult_by_lam_num(X)

    h_sq = 0
    x_pow = x**0
    for ai in list(h_sq_bar):
        h_sq += ai*x_pow
        x_pow *= x
    h_sq = gcd(psi_x, h_sq)
    #print(h_sq)
    #print(factor(h_sq))

    #print(Elkies(E, ell, Fp))


