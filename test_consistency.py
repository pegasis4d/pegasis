#!/usr/bin/env python3
from argparse import ArgumentParser, RawTextHelpFormatter
import re

import sage.all
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from pegasis import PEGASIS

from sage.rings.integer import Integer

description = """
Verifies the outptus of the action against previously collected values.

This is useful to verify that the values computed by pegasis are consistent during development.
"""

parser = ArgumentParser(prog="test_consistency.py", description=description, formatter_class=RawTextHelpFormatter)
parser.add_argument(
    "--loglevel",
    default="INFO",
    help="Set loglevel to one of 'DEBUG', 'INFO' (default), 'WARNING', 'ERROR', 'CRITICAL'",
    choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
)
parser.add_argument(
    "--verify",
    default=False,
    metavar="file",
    help="Select file containing test vectors",
)
parser.add_argument(
    "--generate",
    type=int,
    default=0,
    metavar="N",
    help="Generate N test vectors",
)
parser.add_argument(
    "--save",
    default=False,
    metavar="file",
    help="Select file to save test vectors to",
)
parser.add_argument("-p", type=int, help="Prime size", choices=(500, 1000, 1500, 2000, 4000), default=500)
args = parser.parse_args()


EGA = PEGASIS(args.p)


def montgomery(A):
    return EllipticCurve(EGA.Fp, [0, A, 0, 1, 0])


def generate(nvectors):
    E = EGA.E_start
    vectors = ["# Created with `consistency.py`. Verify with `consistency.py --verify <path to file>`"]

    for _ in range(nvectors):
        ideal = EGA.sample_ideal()
        F = EGA.action(E, ideal)

        # Unfortunately sage does not provide an easily parseable __repr__ class for ideals
        # The following diy serialisation is a little bit brittle
        #
        # - .gens_two() will return two generators as elements of the underlying order
        # - This order is defined inside the number field Q(sqrt(-p))
        # - More precisely, this is defined in pegasis as
        #    self.K = NumberField(name="pi", polynomial = var('x')**2 + self.p)
        # - Elements of K are a list [a, b] representing the element a + pi*b
        # - We print a and b into the test_vectors file

        vector = []
        vector += [f"({x}, {y})" for (x, y) in [list(g) for g in ideal.gens_two()]]
        vector += [f"{E.a2()}"]
        vector += [f"{F.a2()}"]

        vectors += [" ".join(vector)]

        E = F
        print(f"[{_+1}/{nvectors}] Generated vector")

    return vectors


def verify(vectors):

    nvectors = len(vectors)
    success = 0

    for nvector, vector in enumerate(vectors):
        # re.match('([0-9]*)')
        # python's re does not support repeating capture groups, must repeat manually
        pattern = "^"
        # First ideal generator (two coordinates)
        pattern = r"^\(([-0-9]*),\s*([-0-9]*)\)"
        # Second ideal generator (two coordinates)
        pattern += r"\s*\(([-0-9]*),\s*([-0-9]*)\)"
        # Starting Montgomery coefficient
        pattern += r"\s*([-0-9]*)"
        # Ending Montgomery coefficient
        pattern += r"\s*([-0-9]*)"
        pattern += "$"

        matches = re.match(pattern, vector)
        if matches is None:
            print("Parsing error of test vector")
            print(vector)
            continue

        # - Coefficients g11, g12 of the first generator of the ideal
        # - Coefficients g21, g22 of the second generator of the ideal
        # - A Montgomery coefficient of the starting curve
        # - B Montgomery coefficient of the ending curve
        g11, g12, g21, g22, A, B = matches.groups()

        pi = EGA.order.gens()[1]
        g1 = Integer(g11) + Integer(g12) * pi
        g2 = Integer(g21) + Integer(g22) * pi
        ideal = EGA.order * g1 + EGA.order * g2
        A = EGA.Fp(A)
        B = EGA.Fp(B)

        E = montgomery(A)
        F = EGA.action(E, ideal)

        if F.a2() == B:
            success += 1
            print(f"{nvector+1}/{nvectors} (Succes rate: {success/(nvector+1)*100:.2f}) Success!")
        else:
            print(f"{nvector+1}/{nvectors} Failure!")


if args.generate and args.verify:
    print("--generate and --verify are mutually exclusive")
    exit(1)

if args.generate:
    print("Generating test vectors")
    if args.save:
        print(f"Will save vectors to file '{args.save}'")
        # Test that we can actually write before generating vectors
        with open(args.save, "w") as file:
            pass
    else:
        print("Will print test vectors to stdout redirect to file or copy-paste manually")

    vectors = generate(args.generate)

    if args.save:
        with open(args.save, "w") as file:
            vectors = [vector + "\n" for vector in vectors]
            file.writelines(vectors)
    else:
        for vector in vectors:
            print(vector)

if args.verify:
    with open(args.verify, "r") as vectors_file:
        vectors = vectors_file.readlines()
        # Ignore comment lines
        vectors = [vector for vector in vectors if not vector[0] == "#"]
    verify(vectors)
