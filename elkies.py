#!/usr/bin/env python3

import sage.all
from sage.schemes.elliptic_curves.mod_poly import classical_modular_polynomial
from sage.calculus.functional import derivative
from sage.rings.power_series_ring import PowerSeriesRing
from sage.arith.misc import kronecker_symbol
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.functions.other import binomial


def Elkies(E, ell, Fp, lam=None, prev_j=None):
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
    # TODO maybe we want montomery curves later
    p = Fp.characteristic()
    assert ell > 2
    d = (ell - 1) / 2
    j = Fp(E.j_invariant())
    _, _, _, a, b = [Fp(c) for c in E.a_invariants()]

    R = Fp["X,Y"]
    X, Y = R._first_ngens(2)

    E4 = -48 * a
    E6 = 864 * b

    jp = -E6 * j / E4

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
        j_ells = [phi_ell_eval.any_root()]  # There is only one root in this case!
        # assert len(j_ells) == 1, f"j_ells was {j_ells}, and prev_j was {prev_j}"
    else:
        j_ells = [r for r, e in phi_ell_eval.roots()]
        assert len(j_ells) == 2, f"j_ells was {j_ells}"

    # print(j_ells)
    jt = j_ells[0]

    assert jt != 0 and jt != 1728, "encountered E0"

    hc = _derive_hc(Fp, a, b, ell, d, j, jt, jp, E4, E6, phi_ellx, phi_ellxx, phi_ellxy, phi_elly, phi_ellyy)

    if lam and not prev_j:
        if not CheckElkies(E, ell, hc, lam):
            # Do the other one
            assert j_ells[1] != 0 and j_ells[1] != 1728, "encountered E0"
            hc = _derive_hc(
                Fp, a, b, ell, d, j, j_ells[1], jp, E4, E6, phi_ellx, phi_ellxx, phi_ellxy, phi_elly, phi_ellyy
            )
    if lam:
        assert CheckElkies(E, ell, hc, lam)
    return hc


def WeierstrassP(a, b, d):
    # Weierstrass p up to degree d in w=z^2
    coefs = [-a / 5, -b / 7]
    for k in range(2, d):
        ck = 0
        for j in range(k - 1):
            ck = ck + coefs[j] * coefs[k - 2 - j]
        ck = ck * 3 / ((k - 1) * (2 * (k + 1) + 3))
        coefs.append(ck)
    return coefs


def _derive_hc(Fp, a, b, ell, d, j, jt, jp, E4, E6, phi_ellx, phi_ellxx, phi_ellxy, phi_elly, phi_ellyy):
    jtp = -jp * phi_ellx(j, jt) / (ell * phi_elly(j, jt))
    at = -(jtp**2) / (48 * jt * (jt - 1728))
    bt = -(jtp**3) / (864 * jt**2 * (jt - 1728))

    E4t = -48 * at
    E6t = 864 * bt

    fracjs = -(
        jp**2 * phi_ellxx(j, jt) + 2 * ell * jp * jtp * phi_ellxy(j, jt) + (ell * jtp) ** 2 * phi_ellyy(j, jt)
    ) / (jp * phi_ellx(j, jt))
    p1 = (
        ell * fracjs / 2
        + (Fp(ell) / 4) * (E4**2 / E6 - ell * E4t**2 / E6t)
        + (Fp(ell) / 3) * (E6 / E4 - ell * E6t / E4t)
    )

    c = WeierstrassP(a, b, d)
    ckst = WeierstrassP(ell**4 * at, ell**6 * bt, d)
    ct = [ckst[i] - ell * c[i] for i in range(len(c))]  # difference in formula VII.23

    # Computing the coefficients of VII.23 and store as A[i]
    Fpw, w = PowerSeriesRing(Fp, "w", default_prec=d + 1).objgen()
    evp = -(p1 / 2) * w - sum(w ** (i + 1) * ct[i - 1] / ((2 * i + 1) * (2 * i + 2)) for i in range(1, d))
    exp_evp = evp.exp(d + 1)
    A = exp_evp.coefficients()
    C = sum(c[i - 1] * w**i for i in range(1, d + 1))

    # Computing all powers of C starting with zeroth power
    Cpow = [Fpw(1), C]
    for i in range(2, d + 1):
        Cpow.append(C * Cpow[-1])

    # Now doing recurrence relation VII.24
    hc = [1, -p1 / 2]
    for i in range(2, d + 1):
        newcoeff = A[i]
        for k in range(1, i + 1):
            insum = 0
            for j in range(k + 1):
                insum += Fp(binomial(d - i + k, k - j)) * Cpow[k - j][j]
            newcoeff -= insum * hc[i - k]
        hc.append(newcoeff)

    Rx, x = PolynomialRing(Fp, "x").objgen()
    hc = Rx(hc[::-1])
    return hc


def CheckElkies(E, ell, kernel_polynomial, lam):
    r"""Given a kernel polynomial, verify it corresponds to the correct eigenvalue

    If the multiplication-by-lambda map has the following standard form in
    rational maps (c.f. Sutherland's lectures)

        [\lambda] = (u(x)/v(x), r(x, y)/s(x))

    then the eigenvalue is correct, if

        \pi(P) = \lambda P

    on all points in the kernel of the isogeny defined by kernel_polynomial.

    Note that r(x, y) = r(x, 1) * y, by the standard form of isogenies. So, if
    P = (x, y), this is equivalent to

        (x^p, y^p) = (u(x)/v(x), r(x, 1)/s(x) * y)

    for points in ker(\varphi)

    Verifying the first component is easy. To verify the second, we note that

            y^p = r(x, 1)/s(x) * y
        <=> y^{p-1} = r(x, 1)/s(x)
        <=> f(x)^{(p-1)/2} = r(x, 1)/s(x)
        <=> f(x)^{(p-1)/2} * s(x) = r(x, 1)

    where f(x) is the defining equation of the curve E: y^2 = f(x).
    """

    p = E.base_field().characteristic()
    E = E.short_weierstrass_model()

    assert kronecker_symbol(-p, ell) == 1, f"{ell} not Elkies"

    # Defining equation of E: y^2 + g(x)y = f(x)
    f, g = E.hyperelliptic_polynomials()

    # Must be true, because E is in Weierstrass form
    assert g == 0

    # Disabled when python is called with -O flag
    if __debug__:
        # @todo: Can be deleted later. Sanity check.
        # Verify that this is really the same as previous method of getting the
        # defining equation
        # This must be equal, because E is in short weierstrass form
        assert f == kernel_polynomial.parent()([E.a6(), E.a4(), 0, 1])

    # For efficiency: replace lambda with -lambda if -lambda has smaller
    # absolute value
    # If we switch, then we need to multiply the isogeny with -1
    # (which is multiplication by -1 on the y-coordinate)

    if lam > ell - lam:
        lam = ell - lam
        sign = -1
    else:
        sign = 1

    if lam == 1:
        Y = pow(f, (p - 1) / 2, kernel_polynomial)
        return sign * Y == 1

    # Build extension over which the x-coordinates of the kernel are defined
    extension = kernel_polynomial.parent().quotient_ring(kernel_polynomial)

    # Get rational functions of multiplication-by-lambda
    # x = u/v, y = r/x as in the description in the docstring
    x, y = E.multiplication_by_m(lam)[:2]
    u = extension(x.numerator())
    v = extension(x.denominator())
    # Overwrite r(x, y) with r(x, 1)
    r = y.numerator()
    r = extension(r(r.variables()[0], 1).univariate_polynomial())
    s = extension(y.denominator())

    # x^p in the extension
    Xp = extension(pow(kernel_polynomial.parent().gens()[0], p, kernel_polynomial))

    # Verify u(x)/v(x) = x^p
    if u != Xp * v:
        return False

    Y = extension(pow(f, (p - 1) / 2, kernel_polynomial))

    # Verify y^{p-1} = f(x)^{(p-1)/2} = sign * r(x, 1)/s(x)
    return Y * s == sign * r
