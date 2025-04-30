from sage.all import *
from pegasis import PEGASIS, Conjugate

EGA = PEGASIS(500)
order = EGA.order
generator = EGA.w

e = 245
ell = next_prime(randint(0, 2**e))
while kronecker_symbol(-EGA.p, ell) != 1:
    ell = next_prime(randint(0, 2**e))

ideal = ell * order + (generator - Integer(GF(ell)(-EGA.p).sqrt())) * order
assert ideal.norm() == ell
#ideal = EGA.sample_ideal()



E = EGA.action(EGA.E_start, ideal)

print("DONE WITH FIRST.....")
#alpha = frak_a.random_element()
#while not is_pseudoprime(ZZ(alpha.norm()/frak_a.norm())):
#    alpha = frak_a.random_element()
#gen_1, gen_2 = frak_a.gens_two()
#frak_a = EGA.order*(gen_1*alpha.conjugate()/frak_a.norm()) + EGA.order*(gen_2*alpha.conjugate()/frak_a.norm())
#assert is_pseudoprime(frak_a.norm())

E2 = EGA.action(E, Conjugate(EGA.order, ideal))

print("-------------")
print(f"original: {EGA.E_start.j_invariant()}")
print(f"new: {E2.j_invariant()}")

"""
E2 = E2.short_weierstrass_model()
prev_j = None
F = E2
for idx in range(5):
    F = EGA.small_prime_ideal_action(F, EGA.small_ell, EGA.lam, prev_j).codomain()
    print(f"left_{idx}: {F.j_invariant()}")
    if truth == F.j_invariant():
        print("!!!!!!!!!!! WOW")

F = E2
for idx in range(5):
    F = EGA.small_prime_ideal_action(F, EGA.small_ell, EGA.small_ell-EGA.lam, prev_j).codomain()
    print(f"right_{idx}: {F.j_invariant()}")
    if truth == F.j_invariant():
        print("!!!!!!!!!!! WOW")
"""