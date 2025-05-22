#!/usr/bin/env python3

import argparse
import itertools
from time import time

from pegasis import PEGASIS
#from theta_lib.utilities.strategy import precompute_strategy_with_first_eval_and_splitting as precompute
from Theta_dim4.Theta_dim4_sage.pkg.utilities.strategy import precompute_strategy_with_first_eval as precompute


def precompute_strategies(emax):
    erange = range(emax - 10, emax + 1)
    mrange = range(1, 6)
    precompute.precompute(itertools.product(erange, mrange), 4)


if __name__ == "__main__":
    """
    Timing for clapotis. Usage:
    sage --python -O timings.py
        -p P        level (500|1000...)
        --nruns N    number of runs
        --runname R  optional run identifier
        --precompute precompute strategies
    """

    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", type=int, help="Prime size", choices=(500, 1000, 1500, 2000, 4000), default=500)
    parser.add_argument("--nruns", type=int, help="Number of runs", default=5)
    parser.add_argument("--runname", type=str, help="Run name", default="run")
    parser.add_argument("--precompute", action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    assert args.p in [500, 1000, 1500, 2000, 4000]

    EGA = PEGASIS(args.p)
    if args.precompute:
        precompute_start = time()
        print(f":: Precomputing 4d strategies (p={args.p})")
        precompute_strategies(EGA.uv_params.uv_e)
        print()
        print(f"Done. Took {time() - precompute_start:.1f}s.")
        print()

    total_times = [0, 0, 0]
    start = time()

    E = EGA.E_start

    for run in range(1, args.nruns + 1):
        elapsed_time = time() - start

        if run > 1:
            # When run has value n, n-1 runs have been executed
            eta = (args.nruns - (run - 1)) * elapsed_time / (run - 1)
            print()
            print(f":: Run {run}/{args.nruns} (Estimated remaining time: {eta:.3f}s)")
        else:
            print(f":: Run {run}/{args.nruns}")

        frak_a = EGA.sample_ideal()
        E, stats = EGA.action(E, frak_a, verb_return=True)

        total_times[0] += stats["time_step1"]
        total_times[1] += stats["time_step2"]
        total_times[2] += stats["time_step3"]

        average_t1 = total_times[0] / run
        average_t2 = total_times[1] / run
        average_t3 = total_times[2] / run
        average = average_t1 + average_t2 + average_t3

        print()
        print(f"=> Average timings for {run} runs:")
        print(f"  - Step 1 (Finding uv):   {average_t1:>6.3f}s")
        print(f"  - Step 2 (Elkies steps): {average_t2:>6.3f}s")
        print(f"  - Step 3 (4d Isogeny):   {average_t3:>6.3f}s")
        print(f"  - Total:                 {average:>6.3f}s")
