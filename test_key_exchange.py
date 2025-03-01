#!/usr/bin/env python3
from argparse import ArgumentParser, RawTextHelpFormatter
from datetime import datetime, timedelta

from pegasis import PEGASIS

description = """
Testing the key exchange from the pegasis action.

The test does the following

  1) Randomly samples ideals a, b
  2) Computes Ea, Eb as the action of a and b on the starting curve E
  3) Computes Eab, Eba as the action of b on Ea and a on Eb
  4) Verifies that Eab == Eba
"""


def action(EGA, E, ideal):
    global stats
    start = datetime.now()
    E_ideal = EGA.action(E, ideal)
    end = datetime.now()
    stats["time"] += (end - start) / timedelta(microseconds=1) / 10**6
    stats["actions"] += 1
    return E_ideal


if __name__ == "__main__":
    print("Beginning tests...")
    parser = ArgumentParser(prog="test_key_exchange.py", description=description, formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "--loglevel",
        default="INFO",
        help="Set loglevel to one of 'DEBUG', 'INFO' (default), 'WARNING', 'ERROR', 'CRITICAL'",
        choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
    )
    parser.add_argument("-p", type=int, help="Prime size", choices=(500, 1000, 1500, 2000, 4000), default=500)
    args = parser.parse_args()

    stats = {
        # Tries
        "tries": 0,
        # Successes
        "success": 0,
        # Number of action computations done so far
        "actions": 0,
        # Total action time (seconds)
        "time": 0.0,
    }

    while True:
        stats["tries"] += 1

        # Get two instances of PEGASIS to prove that there is no stateful
        # information in the class that could interfere
        alice = PEGASIS(args.p)
        bob = PEGASIS(args.p)

        assert alice.E_start == bob.E_start

        alice_seckey = alice.sample_ideal()
        bob_seckey = bob.sample_ideal()

        alice_pubkey = action(alice, alice.E_start, alice_seckey)
        bob_pubkey = action(bob, bob.E_start, bob_seckey)

        alice_shared_secret = action(alice, bob_pubkey, alice_seckey)
        bob_shared_secret = action(bob, alice_pubkey, bob_seckey)

        if alice_shared_secret != bob_shared_secret:
            status = "Failure"
        else:
            status = "Success!"
            stats["success"] += 1

        _log = []
        _log += [f"#{stats['tries']:>3} {status}"]
        _log += [f"Success rate: {stats['success'] / stats['tries'] * 100:>4.1f}%"]
        _log += [f"Actions computed: {stats['actions']:>4}"]
        _log += [f"Time per action: {stats["time"] / stats["actions"]:4>.2f}s"]

        print(" | ".join(_log))
