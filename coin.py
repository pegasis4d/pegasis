from sage.all import (
    xgcd,
    ceil, floor,
    log,
)
from sage.all import gcd, is_pseudoprime, two_squares, Factorization, ZZ
from sage.all import isqrt
from sage.rings.finite_rings.integer_mod import Mod

small_squares = {
    5: (1, 2), 13: (2, 3), 17: (1, 4), 29: (2, 5), 37: (1, 6), 41: (4, 5),
    53: (2, 7), 61: (5, 6), 73: (3, 8), 89: (5, 8), 97: (4, 9),
    101: (1, 10), 109: (3, 10), 113: (7, 8), 137: (4, 11), 149: (7, 10),
    157: (6, 11), 173: (2, 13), 181: (9, 10), 193: (7, 12), 197: (1, 14),
    229: (2, 15), 233: (8, 13), 241: (4, 15), 257: (1, 16), 269: (10, 13),
    277: (9, 14), 281: (5, 16), 293: (2, 17), 313: (12, 13), 317: (11, 14),
    337: (9, 16), 349: (5, 18), 353: (8, 17), 373: (7, 18), 389: (10, 17),
    397: (6, 19), 401: (1, 20), 409: (3, 20), 421: (14, 15), 433: (12, 17),
    449: (7, 20), 457: (4, 21), 461: (10, 19), 509: (5, 22), 521: (11, 20),
    541: (10, 21), 557: (14, 19), 569: (13, 20), 577: (1, 24), 593: (8, 23),
    601: (5, 24), 613: (17, 18), 617: (16, 19), 641: (4, 25), 653: (13, 22),
    661: (6, 25), 673: (12, 23), 677: (1, 26), 701: (5, 26), 709: (15, 22),
    733: (2, 27), 757: (9, 26), 761: (19, 20), 769: (12, 25), 773: (17, 22),
    797: (11, 26), 809: (5, 28), 821: (14, 25), 829: (10, 27), 853: (18, 23),
    857: (4, 29), 877: (6, 29), 881: (16, 25), 929: (20, 23), 937: (19, 24),
    941: (10, 29), 953: (13, 28), 977: (4, 31), 997: (6, 31), 1009: (15, 28),
    1013: (22, 23), 1021: (11, 30), 1033: (3, 32), 1049: (5, 32), 1061: (10, 31),
    1069: (13, 30), 1093: (2, 33), 1097: (16, 29), 1109: (22, 25), 1117: (21, 26),
    1129: (20, 27), 1153: (8, 33), 1181: (5, 34), 1193: (13, 32), 1201: (24, 25),
    1213: (22, 27), 1217: (16, 31), 1229: (2, 35), 1237: (9, 34), 1249: (15, 32),
    1277: (11, 34), 1289: (8, 35), 1297: (1, 36), 1301: (25, 26), 1321: (5, 36),
    1361: (20, 31), 1373: (2, 37), 1381: (15, 34), 1409: (25, 28), 1429: (23, 30),
    1433: (8, 37), 1453: (3, 38), 1481: (16, 35), 1489: (20, 33), 1493: (7, 38),
    1549: (18, 35), 1553: (23, 32), 1597: (21, 34), 1601: (1, 40), 1609: (3, 40),
    1613: (13, 38), 1621: (10, 39), 1637: (26, 31), 1657: (19, 36), 1669: (15, 38),
    1693: (18, 37), 1697: (4, 41), 1709: (22, 35), 1721: (11, 40), 1733: (17, 38),
    1741: (29, 30), 1753: (27, 32), 1777: (16, 39), 1789: (5, 42), 1801: (24, 35),
    1861: (30, 31), 1873: (28, 33), 1877: (14, 41), 1889: (17, 40), 1901: (26, 35),
    1913: (8, 43), 1933: (13, 42), 1949: (10, 43), 1973: (23, 38), 1993: (12, 43),
    1997: (29, 34), 2017: (9, 44), 2029: (2, 45), 2053: (17, 42), 2069: (25, 38),
    2081: (20, 41), 2089: (8, 45), 2113: (32, 33), 2129: (23, 40), 2137: (29, 36),
    2141: (5, 46), 2153: (28, 37), 2161: (15, 44), 2213: (2, 47), 2221: (14, 45),
    2237: (11, 46), 2269: (30, 37), 2273: (8, 47), 2281: (16, 45), 2293: (23, 42),
    2297: (19, 44), 2309: (10, 47), 2333: (22, 43), 2341: (15, 46), 2357: (26, 41),
    2377: (21, 44), 2381: (34, 35), 2389: (25, 42), 2393: (32, 37), 2417: (4, 49),
    2437: (6, 49), 2441: (29, 40), 2473: (13, 48), 2477: (19, 46), 2521: (35, 36),
    2549: (7, 50), 2557: (21, 46), 2593: (17, 48), 2609: (20, 47), 2617: (4, 51),
    2621: (11, 50), 2633: (28, 43), 2657: (16, 49), 2677: (34, 39), 2689: (33, 40),
    2693: (22, 47), 2713: (3, 52), 2729: (5, 52), 2741: (25, 46), 2749: (30, 43),
    2753: (7, 52), 2777: (29, 44), 2789: (17, 50), 2797: (14, 51), 2801: (20, 49),
    2833: (23, 48), 2837: (34, 41), 2857: (16, 51), 2861: (19, 50), 2897: (31, 44),
    2909: (10, 53), 2917: (1, 54), 2953: (12, 53), 2957: (29, 46), 2969: (37, 40),
    3001: (20, 51), 3037: (11, 54), 3041: (4, 55), 3049: (32, 45), 3061: (6, 55),
    3089: (8, 55), 3109: (30, 47), 3121: (39, 40), 3137: (1, 56), 3169: (12, 55),
    3181: (34, 45), 3209: (20, 53), 3217: (9, 56), 3221: (14, 55), 3229: (27, 50),
    3253: (2, 57), 3257: (11, 56), 3301: (30, 49), 3313: (8, 57), 3329: (25, 52),
    3361: (15, 56), 3373: (3, 58), 3389: (5, 58), 3413: (7, 58), 3433: (27, 52),
    3449: (40, 43), 3457: (39, 44), 3461: (31, 50), 3469: (38, 45), 3517: (6, 59),
    3529: (35, 48), 3533: (13, 58), 3541: (25, 54), 3557: (34, 49), 3581: (10, 59),
    3593: (28, 53), 3613: (42, 43), 3617: (41, 44), 3637: (39, 46), 3673: (37, 48),
    3677: (14, 59), 3697: (36, 49), 3701: (26, 55), 3709: (30, 53), 3733: (22, 57),
    3761: (25, 56), 3769: (13, 60), 3793: (33, 52), 3797: (41, 46), 3821: (10, 61),
    3833: (32, 53), 3853: (3, 62), 3877: (31, 54), 3881: (20, 59), 3889: (17, 60),
    3917: (14, 61), 3929: (35, 52), 3989: (25, 58), 4001: (40, 49), 4013: (13, 62),
    4021: (39, 50), 4049: (32, 55), 4057: (24, 59), 4073: (37, 52), 4093: (27, 58),
    4129: (23, 60), 4133: (17, 62), 4153: (43, 48), 4157: (26, 59), 4177: (9, 64),
    4201: (40, 51), 4217: (11, 64), 4229: (2, 65), 4241: (4, 65), 4253: (38, 53),
    4261: (6, 65), 4273: (32, 57), 4289: (8, 65), 4297: (24, 61), 4337: (44, 49),
    4349: (43, 50), 4357: (1, 66), 4373: (23, 62), 4397: (26, 61), 4409: (40, 53),
    4421: (14, 65), 4441: (29, 60), 4457: (19, 64), 4481: (16, 65), 4493: (2, 67),
    4513: (47, 48), 4517: (46, 49), 4549: (18, 65), 4561: (31, 60), 4597: (41, 54),
    4621: (30, 61), 4637: (34, 59), 4649: (5, 68), 4657: (39, 56), 4673: (7, 68),
    4721: (25, 64), 4729: (45, 52), 4733: (37, 58), 4789: (42, 55), 4793: (13, 68),
    4801: (24, 65), 4813: (18, 67), 4817: (41, 56), 4861: (10, 69), 4877: (34, 61),
    4889: (20, 67), 4909: (3, 70), 4933: (33, 62), 4937: (29, 64), 4957: (14, 69),
    4969: (37, 60), 4973: (22, 67), 4993: (32, 63), 5009: (28, 65), 5021: (11, 70),
    5077: (6, 71), 5081: (40, 59), 5101: (50, 51), 5113: (48, 53), 5153: (23, 68),
    5189: (17, 70), 5197: (29, 66), 5209: (5, 72), 5233: (7, 72), 5237: (14, 71),
    5261: (19, 70), 5273: (28, 67), 5281: (41, 60), 5297: (16, 71), 5309: (50, 53),
    5333: (2, 73), 5381: (34, 65), 5393: (8, 73), 5413: (38, 63), 5417: (44, 59),
    5437: (26, 69), 5441: (20, 71), 5449: (43, 60), 5477: (1, 74), 5501: (5, 74),
    5521: (36, 65), 5557: (9, 74), 5569: (40, 63), 5573: (47, 58), 5581: (35, 66),
    5641: (4, 75), 5653: (18, 73), 5657: (44, 61), 5669: (38, 65), 5689: (8, 75),
    5693: (43, 62), 5701: (15, 74), 5717: (26, 71), 5737: (51, 56), 5741: (29, 70),
    5749: (50, 57), 5801: (5, 76), 5813: (22, 73), 5821: (14, 75), 5849: (35, 68),
    5857: (9, 76), 5861: (31, 70), 5869: (45, 62), 5881: (16, 75), 5897: (11, 76),
    5953: (52, 57), 5981: (50, 59), 6029: (10, 77), 6037: (41, 66), 6053: (47, 62),
    6073: (12, 77), 6089: (40, 67), 6101: (25, 74), 6113: (28, 73), 6121: (45, 64),
    6133: (7, 78), 6173: (53, 58), 6197: (34, 71), 6217: (21, 76), 6221: (50, 61),
    6229: (30, 73), 6257: (4, 79), 6269: (37, 70), 6277: (6, 79), 6301: (26, 75),
    6317: (29, 74), 6329: (20, 77), 6337: (36, 71), 6353: (32, 73), 6361: (40, 69),
    6373: (17, 78), 6389: (55, 58), 6397: (54, 59), 6421: (39, 70), 6449: (7, 80),
    6469: (50, 63), 6473: (43, 68), 6481: (9, 80), 6521: (11, 80), 6529: (48, 65),
    6553: (37, 72), 6569: (13, 80), 6577: (4, 81), 6581: (41, 70), 6637: (54, 61),
    6653: (53, 62), 6661: (10, 81), 6673: (52, 63), 6689: (17, 80), 6701: (35, 74),
    6709: (25, 78), 6733: (3, 82), 6737: (31, 76), 6761: (19, 80), 6781: (34, 75),
    6793: (48, 67), 6829: (30, 77), 6833: (47, 68), 6841: (21, 80), 6857: (56, 61),
    6869: (55, 62), 6917: (26, 79), 6949: (15, 82), 6961: (20, 81), 6977: (44, 71),
    6997: (39, 74), 7001: (35, 76), 7013: (17, 82), 7057: (1, 84), 7069: (38, 75),
    7109: (47, 70), 7121: (55, 64), 7129: (27, 80), 7177: (11, 84), 7193: (52, 67),
    7213: (18, 83), 7229: (2, 85), 7237: (26, 81), 7253: (23, 82), 7297: (39, 76),
    7309: (35, 78), 7321: (60, 61), 7333: (58, 63), 7349: (25, 82), 7369: (12, 85),
    7393: (47, 72), 7417: (19, 84), 7433: (53, 68), 7457: (41, 76), 7477: (9, 86),
    7481: (16, 85), 7489: (33, 80), 7517: (11, 86), 7529: (40, 77), 7537: (36, 79),
    7541: (50, 71), 7549: (18, 85), 7561: (44, 75), 7573: (2, 87), 7577: (59, 64),
    7589: (58, 65), 7621: (15, 86), 7649: (55, 68), 7669: (10, 87), 7673: (28, 83),
    7681: (25, 84), 7717: (34, 81), 7741: (46, 75), 7753: (3, 88), 7757: (19, 86),
    7789: (30, 83), 7793: (7, 88), 7817: (61, 64), 7829: (50, 73), 7841: (40, 79),
    7853: (58, 67), 7873: (57, 68), 7877: (49, 74), 7901: (26, 85), 7933: (43, 78),
    7937: (4, 89), 7949: (35, 82), 7993: (53, 72), 8009: (28, 85), 8017: (31, 84),
    8053: (22, 87), 8069: (62, 65), 8081: (41, 80), 8089: (60, 67), 8093: (37, 82),
    8101: (1, 90), 8117: (14, 89), 8161: (40, 81), 8209: (55, 72), 8221: (11, 90),
    8233: (48, 77), 8237: (29, 86), 8269: (13, 90), 8273: (23, 88), 8293: (47, 78),
    8297: (4, 91), 8317: (6, 91), 8329: (52, 75), 8353: (28, 87), 8369: (25, 88),
    8377: (51, 76), 8389: (17, 90), 8429: (50, 77), 8461: (19, 90), 8501: (55, 74),
    8513: (7, 92), 8521: (36, 85), 8537: (16, 91), 8573: (43, 82), 8581: (65, 66),
    8597: (26, 89), 8609: (47, 80), 8629: (23, 90), 8641: (60, 71), 8669: (38, 85),
    8677: (46, 81), 8681: (20, 91), 8689: (15, 92), 8693: (58, 73), 8713: (8, 93),
    8737: (41, 84), 8741: (50, 79), 8753: (17, 92), 8761: (56, 75), 8821: (30, 89),
    8837: (1, 94), 8849: (65, 68), 8861: (5, 94), 8893: (53, 78), 8929: (60, 73),
    8933: (47, 82), 8941: (29, 90), 8969: (35, 88), 9001: (51, 80), 9013: (38, 87),
    9029: (2, 95), 9041: (4, 95), 9049: (20, 93), 9109: (55, 78), 9133: (22, 93),
    9137: (64, 71), 9157: (54, 79), 9161: (44, 85), 9173: (62, 73), 9181: (30, 91),
    9209: (53, 80), 9221: (14, 95), 9241: (5, 96), 9257: (59, 76), 9277: (21, 94),
    9281: (16, 95), 9293: (58, 77), 9337: (11, 96), 9341: (46, 85), 9349: (18, 95),
    9377: (56, 79), 9397: (66, 71), 9413: (2, 97), 9421: (45, 86), 9433: (28, 93),
    9437: (34, 91), 9461: (25, 94), 9473: (8, 97), 9497: (61, 76), 9521: (40, 89),
    9533: (53, 82), 9601: (24, 95), 9613: (3, 98), 9629: (5, 98), 9649: (57, 80),
    9661: (69, 70), 9677: (29, 94), 9689: (35, 92), 9697: (56, 81), 9721: (64, 75),
    9733: (18, 97), 9749: (55, 82), 9769: (45, 88), 9781: (41, 90), 9817: (4, 99),
    9829: (15, 98), 9833: (37, 92), 9857: (44, 89), 9901: (10, 99), 9929: (52, 85),
    9941: (70, 71), 9949: (43, 90), 9973: (57, 82)}

# Primes 1 mod 4 / 3 mod 4
from const_precomp import (
        good_prime_prod_10k, bad_prime_prod_10k,
        good_prime_prod_100k, bad_prime_prod_100k,
        good_prime_prod_500k, bad_prime_prod_500k
)
good_prime_prod = good_prime_prod_10k
bad_prime_prod = bad_prime_prod_10k

def coin_pair_generator(n1, n2, e, exp_bound=None):
    """
    Find all pairs (u, v) s.t. u*n1 + v*n2 = 2^t for t <= e.
    If exp_bound is set, stop returning solution for t > t0 + exp_bound, where
    t0 is the minimum t for which a solution is found.
    """
    one, a, b = xgcd(n1, n2)
    # Assuming (a, b) = 1 - if gcd is 2^i we can divide before
    if one != 1:
        raise ValueError(f"{n1 = } and {n2 = } are not coprime")
    a, b, n1, n2 = map(int, [a, b, n1, n2])

    t0 = None
    for t in range(floor(log(min(n1, n2), 2)), e+1):
        if exp_bound and t0 and t > t0 + exp_bound:
            break
        l = 2**t
        if l < max(n1, n2):
            continue
        if a < 0:
            k = ((-l*a - 1) // n2) + 1 # = ceil(-l*a / n2)
            k_max = l*b // n1 # = floor(l*b / n1)
            while k <= k_max:
                u = l*a + k*n2
                v = l*b - k*n1

                k += 1
                if u % 2 == v % 2 == 0:
                    continue # Skip if both even (already included)
                if not t0:
                    t0 = t
                yield (ZZ(u), ZZ(v), ZZ(t))

        elif b < 0:
            k = ((-l*b - 1) // n1) + 1  # = ceil(-l*b / n1)
            k_max = l*a // n2 # = floor(l*a / n2)
            while k <= k_max:
                u = l*a - k*n2
                v = l*b + k*n1

                k += 1
                if u % 2 == v % 2 == 0:
                    continue
                if not t0:
                    t0 = t
                yield (ZZ(u), ZZ(v), ZZ(t))


def two_squares_factored(factors):
    """
    This is the function `two_squares` from sage, except we give it the
    factorisation of n already.
    """
    F = factors
    for (p,e) in F:
        if e % 2 == 1 and p % 4 == 3:
            raise ValueError("%s is not a sum of 2 squares"%n)

    n = factors.expand()
    if n == 0:
        return (0, 0)
    a = ZZ.one()
    b = ZZ.zero()
    for (p,e) in F:
        if p == 1:
            continue
        if e >= 2:
            m = p ** (e//2)
            a *= m
            b *= m
        if e % 2 == 1:
            if p == 2:
                # (a + bi) *= (1 + I)
                a,b = a - b, a + b
            else:  # p = 1 mod 4
                if p in small_squares:
                    r,s = small_squares[p]
                else:
                    # Find a square root of -1 mod p.
                    # If y is a non-square, then y^((p-1)/4) is a square root of -1.
                    y = Mod(2,p)
                    while True:
                        s = y**((p-1)/4)
                        if not s*s + 1:
                            s = s.lift()
                            break
                        y += 1
                    # Apply Cornacchia's algorithm to write p as r^2 + s^2.
                    r = p
                    while s*s > p:
                        r,s = s, r % s
                    r %= s

                # Multiply (a + bI) by (r + sI)
                a,b = a*r - b*s, b*r + a*s

    a = a.abs()
    b = b.abs()
    assert a*a + b*b == n
    if a <= b:
        return (a,b)
    else:
        return (b,a)

def rep_gcd(a, b):
    """
    Given a and b returns (a1, g) where a1*g = a and g contains
    the factors in common between a and b
    """
    out = 1
    g = gcd(a, b)
    while g != 1:
        out *= g
        a /= g
        g = gcd(a, g)

    return a, out

def sum_of_squares_friendly(n):
    """
    We can write any n = x^2 + y^2 providing that there
    are no prime power factors p^k | n such that
    p = 3 mod 4 and k odd.
    """
    # We consider the odd part of n and try and determine if there are bad factors
    n_val = n.valuation(2)
    n_odd = n // (2**n_val)
    fact = [(2, n_val)]

    if n_odd % 4 == 3:
        return False, fact

    n_odd, bad_cof = rep_gcd(n_odd, bad_prime_prod)

    sbf = isqrt(bad_cof)
    if sbf**2 == bad_cof:
        fact.append((sbf, 2))
    else:
        return False, fact

    # Good primes 1 mod 4
    n_odd, good_cof = rep_gcd(n_odd, good_prime_prod)
    good_cof = ZZ(good_cof)

    if n_odd == 1:
        return True, Factorization([*fact, *good_cof.factor()])

    else:
        return is_pseudoprime(n_odd), Factorization([*fact, *good_cof.factor(), (n_odd, 1)])

def sum_of_squares_friendly_pair(n1, n2):
    """
    sum_of_squares_friendly applied directly to a pair
    to minimize primality testing with early rejection
    """
    # We consider the odd part of n and try and determine if there are bad factors
    n1_val = n1.valuation(2)
    n1_odd = n1 // (2**n1_val)
    fact1 = [(2, n1_val)]
    n2_val = n2.valuation(2)
    n2_odd = n2 // (2**n2_val)
    fact2 = [(2, n2_val)]

    if n1_odd % 4 == 3 or n2_odd % 4 == 3:
        return False, fact, False, fact

    n1_odd, bad_cof1 = rep_gcd(n1_odd, bad_prime_prod)
    sbf1 = isqrt(bad_cof1)
    if sbf1**2 == bad_cof1:
        fact1.append((sbf1, 2))
    else:
        return False, fact1, False, fact2

    n2_odd, bad_cof2 = rep_gcd(n2_odd, bad_prime_prod)
    sbf2 = isqrt(bad_cof2)
    if sbf2**2 == bad_cof2:
        fact2.append((sbf2, 2))
    else:
        return False, fact1, False, fact2

    # Good primes 1 mod 4
    n1_odd, good_cof1 = rep_gcd(n1_odd, good_prime_prod)
    good_cof1 = ZZ(good_cof1)
    if n1_odd != 1 and not is_pseudoprime(n1_odd):
        return False, fact1, False, fact2

    n2_odd, good_cof2 = rep_gcd(n2_odd, good_prime_prod)
    good_cof2 = ZZ(good_cof2)
    if n2_odd != 1 and not is_pseudoprime(n2_odd):
        return False, fact1, False, fact2

    return True, Factorization([*fact1, *good_cof1.factor(), (n1_odd, 1)]), True, Factorization([*fact2, *good_cof2.factor(), (n2_odd, 1)])


def sum_of_squares(n):
    """
    Attempts to compute x,y such that n = x^2 + y^2
    """
    n = ZZ(n)

    b, fact = sum_of_squares_friendly(n)
    if not b:
        return []
    return two_squares_factored(fact)


def sum_of_squares_pair(n1, n2):
    """
    Sum of squares for both n1 and n2
    """
    n1 = ZZ(n1)
    n2 = ZZ(n2)

    b1, fact1, b2, fact2 = sum_of_squares_friendly_pair(n1, n2)
    if not (b1 and b2):
        return [], []

    return two_squares_factored(fact1), two_squares_factored(fact2)


def sos_pair_generator(n1, n2, e, n_squares=1, exp_bound=None):
    """
    Compute pairs (u, v) such that u*n1 + v*n2 = 2^t, with t <= e
    and u or v (or both) sum of squares.
    Input:
        - n1 and n2: starting values
        - e: target exponent of 2
        - n_squares: how many of u and v must be sum of squares
            (1 or 2, default 1)
        - exp_bound: (optional) exp_bound argument to pass to
            coin_pair_generator
    Output: a triple (u, v, t) where
        - t is such that u*n1 + v*n2 = 2^t, t <= e
        - the non s.o.s. among u and v (if any) is an integer
        - the s.o.s. are returned as a pair (x, y) such that
            x^2 + y^2 = u
    """

    for u, v, t in coin_pair_generator(n1, n2, e, exp_bound):
        if n_squares == 1:
            dec = sum_of_squares(u)
            if dec != []:
                yield [dec, v, t]
            dec = sum_of_squares(v)
            if dec != []:
                yield [u, dec, t]

        elif n_squares == 2:
            dec_u = sum_of_squares(u)
            if dec_u == []:
                continue
            dec_v = sum_of_squares(v)
            if dec_v != []:
                yield [dec_u, dec_v, t]

        else:
            raise ValueError("n_squares must be 1 or 2")

    return false


def remove_primes(u, primes):
    """
    Remove odd exponents of primes from u.
    Returns (x, g) such that x*g = u, and g contains
    the odd exponents of `primes` in u.
    """
    g = 1
    for pi in primes:
        val_i = u.valuation(pi) % 2
        g *= pi**val_i
        u /= pi**val_i

    return u, g

def xsos_pair_generator(n1, n2, e, allowed_primes, n_squares=1, exp_bound=None):
    """
    Compute pairs (u, v) such that u*n1 + v*n2 = 2^t, with t <= e and u or v
    (or both) sum of squares.
    Input:
        - n1 and n2: starting values
        - e: target exponent of 2
        - n_squares: how many of u and v must be sum of squares (1 or 2,
          default 1)
        - allowed primes: a list of primes that can be removed from u and v
          when checking for sum of squares, e.g.  u = u'*3 with u' sum of
          squares and 3 allowed prime
        - exp_bound: (optional) argument to pass to coin_pair_generator
    Output: a triple (u, v, t) where
        - t is such that u*n1 + v*n2 = 2^t, t <= e
        - the non s.o.s. among u and v (if any) is an integer
        - the s.o.s. are returned as a triple (x, y, g) such that g*(x^2 + y^2)
          = u
    """
    for u, v, t in coin_pair_generator(n1, n2, e, exp_bound):
        if n_squares == 1:
            uu, gu = remove_primes(u, allowed_primes)
            dec = sum_of_squares(uu)
            if dec != []:
                x, y = dec
                yield [(x, y, gu), v, t]

            vv, gv = remove_primes(v, allowed_primes)
            dec = sum_of_squares(vv)
            if dec != []:
                x, y = dec
                yield [u, (x, y, gv), t]

        elif n_squares == 2:
            uu, gu = remove_primes(u, allowed_primes)
            vv, gv = remove_primes(v, allowed_primes)

            if (uu % 4 == 3) or (vv % 4 == 3):
                continue

            dec_u, dec_v = sum_of_squares_pair(uu, vv)
            if dec_u == [] or dec_v == []:
                continue

            xu, yu = dec_u
            xv, yv = dec_v
            yield [(xu, yu, gu), (xv, yv, gv), t]

        else:
            raise ValueError("n_squares must be 1 or 2")

    return False


if __name__ == '__main__':
    from sage.all import set_random_seed, randint
    _r = randint(1, 1000)
    set_random_seed(_r)
    print(f'random seed: {_r}')

    # benchmarking for p = 2**lambda
    lb = 1500
    k = lb//2
    d1 = randint(2**k-2**(k//2), 2**k+2**(k//2))
    d2 = d1
    while gcd(d1, d2) != 1:
        d2 = randint(2**k-2**(k-1), 2**k+2**(k-1))
    print(f'{d1*d2}\n{2**lb}')

    import time
    tic = time.time()
    res = sos_pair_generator(d1, d2, lb, n_squares=1)
    print(f'time 1: {time.time() - tic:.3f}')
    print(f'{res = }')

    if res:
        u, v, t = res
        if not u in zz:
            u = u[0]**2 + u[1]**2
        if not v in zz:
            v = v[0]**2 + v[1]**2
        assert u*d1 + v*d2 == 2**t and t <= lb

