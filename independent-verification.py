#!/usr/bin/env python3

from logging import getLogger

from sage.arith.misc import kronecker_symbol
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.schemes.elliptic_curves.mod_poly import classical_modular_polynomial
from sage.rings.fast_arith import prime_range

from multiprocessing import Process, Queue, cpu_count

from pegasis import PEGASIS

pegasis_logger = getLogger("pegasis")
pegasis_logger.setLevel("WARNING")

SENTINEL = None


def test(ell):
    EGA = PEGASIS(500)

    ideal = ell * EGA.order + (EGA.w - Integer(GF(ell)(-EGA.p).sqrt())) * EGA.order

    E = EGA.action(EGA.E_start, ideal)

    modular_polynomial = classical_modular_polynomial(ell, j=EGA.E_start.j_invariant())

    return (modular_polynomial(E.j_invariant()) == 0, ell)


def __test(q_in, q_out):
    while True:
        ell = q_in.get()

        if ell is SENTINEL:
            q_out.put((SENTINEL, SENTINEL))
            return

        q_out.put(test(ell))


if __name__ == "__main__":
    EGA = PEGASIS(500)
    good_ells = [ell for ell in prime_range(2, 2**16) if kronecker_symbol(-EGA.p, ell) == 1]

    q_in = Queue()
    q_out = Queue()

    for ell in good_ells:
        q_in.put(ell)

    processes = [Process(target=__test, args=(q_in, q_out)) for _ in range(cpu_count() - 1)]

    for process in processes:
        process.start()

    sentinels = 0
    finished = False
    tries = 0
    success = 0
    max_ell_computed = 0

    while sentinels < cpu_count() - 1:

        result, ell = q_out.get()

        if result is SENTINEL:
            sentinels += 1
            continue

        tries += 1
        if result:
            success += 1

        max_ell_computed = max(max_ell_computed, ell)

        print(f"Success rate {success / tries * 100:.1f}% of {tries:>4} attempts (Max ell: {max_ell_computed})")
