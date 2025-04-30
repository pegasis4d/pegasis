from sage.all import *

import coin
import const_precomp
from ideals import ideal_to_sage

class UV_params:
    """
    Parameter set for running the algorithm and finding solutions to the coin
    equation. By default, takes as input the securty level
    (500|1000|1500|2000|4000) and proceed to setup the predefined parameters
    for such level. Optionally, different parameters can be provided as a
    dictionary.

    Parameters:
    - p: the underlying prime
    - O: them maximal order Z[(pi + 1)/2]
    - uv_e: the 2-torsion available for solving the norm equation; by default
      is v_2(p+1) - 2;
    - spr: precomputed rounded square root of p;
    - allowed_primes: small elkies primes that can be removed from the norm
      equation and computed separately; a different (heuristically) optimal set
      is available for each level; to update the list, and all the associated
      parameters, use the functions `set_allowed_primes` and
      `add_allowed_prime` described below;
    - allowed_primes_prod: the product of all allowed_primes;
    - sos_allowed_primes: primes that can be removed from the sum of squares
      step; consists of the primes 3 mod 4 in allowed_primes;
    - norm_bound: scaling factor used to bound the norms of the short vectors
      sampled from the ideal lattice; only vectors with norm smaller than
      norm_bound*spr are considered;
    - norm_cff_bound: short vectors are computed as linear combinations of an
      LLL basis; this value bound the coefficients of the combinations; by
      default is the (rounded) square root of norm_bound;
    - comb_bound: the number of short vectors that are combined to try to find
      a solution to the coin equation; by default is 10, meaning 100
      combinations are tried;
    - n_squares: how many between u and v must be sums of squares; default 2;
    - sol_bound: how many solutions to try before returning the best one found;
      default 1;

    Methods:
    - set_allowed_primes(ls): takes in input a list `ls` and set
      `self.allowed_primes` to the Elkies primes contained in `ls`; updates the
      related parameters accordingly;
    - add_allowed_prime(ell): if the given prime `ell` is Elkies, it is added
      to `self.allowed_primes` and the function returns True; otherwise, it
      returns False.
    """

    def __init__(self, level, params=None):
        level = int(level)
        if level == 100:
            self.init_100()
        elif level == 500:
            self.init_500()
        elif level == 1000:
            self.init_1000()
        elif level == 1500:
            self.init_1500()
        elif level == 2000:
            self.init_2000()
        elif level == 4000:
            self.init_4000()
        else:
            raise ValueError('unknown level')

        # Arithmetic
        self.Fp = GF(self.p)
        self.Fp2 = GF((self.p, 2), name='i', modulus=var('x')**2 + 1)
        self.K = NumberField(name="pi", polynomial = var('x')**2 + self.p)
        self.w = self.K.gens()[0]
        self.max_order = self.K.maximal_order()
        self.order = self.K.order_of_conductor(2)

        # Starting curve
        E_start = EllipticCurve(self.Fp, [0, self.A, 0, 1, 0])
        self.E_start = E_start

        # UV parameters
        self.uv_e = self.e - 3
        self.spr = isqrt(self.p)
        p_3mod4 = list(filter(lambda x: x%4 == 3, self.allowed_primes))
        self.sos_allowed_primes = p_3mod4
        self.allowed_primes_prod = prod(self.allowed_primes)
        self.norm_bound = ZZ(self.p).nbits()
        self.norm_cff_bound = isqrt(self.norm_bound) + 1
        self.comb_bound = 10
        self.n_squares = 2
        self.sol_bound = 1

        self.two_left = ideal_to_sage([[2, 0], [-1 / 2, 1 / 2]], self.max_order)
        self.two_right = ideal_to_sage([[2, 0], [1 / 2, 1 / 2]], self.max_order)

        # Optional parameter update
        if params:
            for p_key, p_val in params.items():
                setattr(self,p_key,p_value)

    def init_100(self):
        self.f = 77
        self.e = 100
        self.p = self.f * 2**self.e - 1
        self.A = 86576444069281248423336823187435
        allowed_primes = [5, 7, 11]
        self.allowed_primes = [ZZ(li) for li in allowed_primes]

    def init_500(self):
        self.f = 33
        self.e = 503
        self.p = self.f * 2**self.e - 1
        self.A = 846923981206860774667188923127927404866138431504857441785556155002892474407686465068998563030997400041497236660113832719116298437902864009587446028542241
        allowed_primes = [3, 7, 11, 13]
        self.allowed_primes = [ZZ(li) for li in allowed_primes]

    def init_1000(self):
        self.f = 15
        self.e = 1004
        self.p = self.f * 2**self.e - 1
        self.A = 2447932165031884901753526747235802754200581132492013990029309589916334680445261935406178846882142279733036010435356728900482083533977052489550706097041090570437350041707592382243329641233414772511433265807801169694848478362818691704488235495552842938554477902999119670514198118427921959251521459174283757
        allowed_primes = [3, 5, 7, 11]
        self.allowed_primes = [ZZ(li) for li in allowed_primes]

    def init_1500(self):
        self.f = 9
        self.e = 1551
        self.p = self.f * 2**self.e - 1
        self.A = 470690138325573302381284028925922066259197385545434613961546293117200415888454393210985745811672527707750649809464293986828934656433695254827553343115710988038803620829829607213986283428022502352610750230825225753910990410861379903133396479857538580518638524776922388484859923418301448669856130467309798601801912903218678382572894516230515348525582714674329329665444893417613957336453258663721257841338989369604936267692408989319582217755257612152597484632197929506422
        allowed_primes = [3, 5, 11]
        self.allowed_primes = [ZZ(li) for li in allowed_primes]
        coin.good_prime_prod = const_precomp.good_prime_prod_100k
        coin.bad_prime_prod = const_precomp.bad_prime_prod_100k

    def init_2000(self):
        self.f = 51
        self.e = 2026
        self.p = self.f * 2**self.e - 1
        self.A = 81522175295034447595829362752089781064926402323585792843903881571299469475765731013589255502541646344652265632612080053850583985808199619939412577427961666539089410884159333672468730652939947531656027845338622561223839271607859219023997914401035458486761172816051586319123853202351034903432110329760972693423592535540529210197485000054794934251253577983683962395651835481200021277459090885409493264948676810331609354889543129654381859017748250839312739412603528665722325984598111481835819726605526313101275796618189346971789271155935744998546622047982841689962568590619138079006213242588187631607235805859652542
        allowed_primes = [3, 7, 11, 17]
        self.allowed_primes = [ZZ(li) for li in allowed_primes]
        coin.good_prime_prod = const_precomp.good_prime_prod_100k
        coin.bad_prime_prod = const_precomp.bad_prime_prod_100k

    def init_4000(self):
        self.f = 63
        self.e = 4084
        self.p = self.f * 2**self.e - 1
        self.A = 3861182440791393427115834319919792363082397953235366672989642232032664601358637530587592475724492228157994152478762370779212684014665350275282740708992301403517860204937961301178353993394123923639955895951016030909047562464435151058330155507538738419078593938731652956043256633693827720453080631845604091259290706126622597714437085574776283248730437225629704782127771641525437187101408298218133833634726212421345805005641803578982727095077266614175977073363364911048284637511496334424899430835066269145012363608931646817082193451123867313226197275905378376442493746564842029440001132446848573232704980649005361798250530086069725335482701101306335776380418104681645935311890062938192551774656426376713685296742945987418112888455314815190689399358568871101254106213839110567602195750401800211431397684678743127540890669195457658925252532743828619079040054116684464335116945294570022631103926522225020600966809347462518428416543333615647560326220108209234883366833121125898408957117320377370772129934485893908787438550276617471673661822552984223791369029779538709395655155777458500406072451055588404055831790775259870437154343571081246432780811947116682159194253450328736841478964648801067542070339026439293541698676665959781451048203
        allowed_primes = [3, 7, 11, 17, 19]
        self.allowed_primes = [ZZ(li) for li in allowed_primes]
        coin.good_prime_prod = const_precomp.good_prime_prod_100k
        coin.bad_prime_prod = const_precomp.bad_prime_prod_100k

    def set_allowed_primes(self, allowed_primes):
        """
        Set the list of allowed primes to a given list by filtering non-Elkies
        primes, and update parameters accordingly.
        """
        self.allowed_primes = list(filter(
            lambda ell: kronecker(-self.p, ell) == 1, allowed_primes))
        self.allowed_primes_prod = prod(self.allowed_primes)
        p_3mod4 = list(filter(lambda x: x%4 == 3, self.allowed_primes))
        self.sos_allowed_primes = p_3mod4

    def add_allowed_prime(self, ell):
        """
        Check if ell is an elkies prime; if so, add it to the list of allowed
        primes, update the related parameters accordingly, and return True;
        otherwise, return False.
        """
        if ell in self.allowed_primes or kronecker(-self.p, ell) != 1:
            return False
        self.allowed_primes.append(ell)
        self.allowed_primes_prod *= ell
        if ell%4 == 3:
            self.sos_allowed_primes.append(ell)
        return True
