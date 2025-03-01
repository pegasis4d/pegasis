# PEGASIS

## Examples

This code is compatible with `sage 10.5`.

The general API is as follows

```python
from pegasis import PEGASIS
EGA = PEGASIS(500) # 1000 | 1500 | 2000 | 4000
a = EGA.sample_ideal()
Ea = EGA.action(EGA.E_start, a)
```

### Key Exchange

Run (and check) a CSIDH key exchange

```python
from pegasis import PEGASIS

lvl = 500 # 1000 | 1500 ...

# Two distinct instances: there is no shared stateful information
Alice = PEGASIS(lvl)
Bob = PEGASIS(lvl)

assert Alice.E_start == Bob.E_start

alice_seckey = Alice.sample_ideal()
alice_pubkey = Alice.action(Alice.E_start, alice_seckey)

bob_seckey = Bob.sample_ideal()
bob_pubkey = Bob.action(Bob.E_start, bob_seckey)

alice_shared_secret = Alice.action(bob_pubkey, alice_seckey)
bob_shared_secret = Bob.action(alice_pubkey, bob_seckey)

assert alice_shared_secret == bob_shared_secret
```

### Timings

Timings contained in the paper may be verified by running
```bash
$ python -O time_pegasis.py [-p {500,1000,1500,2000,4000}] [--nruns <num>] [--precompute | --no-precompute]
```

- `-p` set the underlying prime size
- `--nruns` set how many runs to execute
- `--precompute`, pass this option to precompute strategies for 4D isogeny
  computations before beginning

**The `-O` switch must be passed to python, else expensive assertions will
affect timings**

Defaults to `-p 500 --nrun 5 --no-precompute`

### Testing

To test whether the library computes key exchanges correctly, run
`test_key_exchange.py`. It performs the key exchange example above in an
infinite loop and prints some crude timing information whilst doing so.

(For development) To test whether the action is computing the same outputs every
time, use `test_consistency.py`. It allows you to generate test vectors as
follows

```bash
$ python -O test_consistency.py --generate 10 --save /tmp/test_vectors.txt
```

then after some changes have been made, you can execute

```bash
$ python -O test_consistency.py --verify /tmp/test_vectors.txt
```

We have supplied test vectors in the repository `test_vectors.txt`.

## Additional information

- The code is a proof of concept, and as such contains various unnecessary
  assertions that may slow down the computation; to ignore them, use the
  `-O` flag (e.g. `sage --python -O script.py`)
- additional information may be printed by adding the following snippet at the
  beginning of the code:

```python
import logging
logging.getLogger("pegasis").setLevel(logging.INFO) #or logging.DEBUG
logging.getLogger("lcsidh").setLevel(logging.INFO)
```

## Project structure

- `theta_lib`: code for computing with 4D isogenies;
- `pegasis.py`: contains `PEGASIS`, the main class running the algorithm
- `coin.py`: utilities to solve the norm equation with u and v sums of squares
- `const_precomp.py`: precomputed constants for trial division
- `elkies.py`: Elkies algorithm for computing rational isogenies
- `ideal.py`: code for handling fractional ideals
- `lcsidh.py`: solve the norm equation for a given ideal
- `xonly.py`: code for x-only arithmetic on Montgomery curves
- `uv_params.py`: contains the parameters for solving the norm equation
  ("finding uv")
- `test_vectors.txt`: contains test vectors for verification
