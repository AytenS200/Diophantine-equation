# diophantine

A Diophantine equation solver written from scratch in pure Python. No solver or math libraries are used — not even `math`, as `isqrt` is implemented here with Newton's method; only the standard library (`fractions`, `itertools`, `heapq`, `random`, `re`, `argparse`, `json`, `sys`). Every answer comes with a step-by-step explanation and, when there are no solutions, a reason.

Requires Python 3.8 or newer.

## Installation

From the project folder:

```bash
pip install .            # installs the `diophantine` command and the importable package
pip install -e .         # editable install, for development
```

From a built distribution:

```bash
pip install diophantine-1.0.0-py3-none-any.whl
```

Without installing, you can also run it from the project folder with `python -m diophantine ...`.

## Terminal usage

```bash
diophantine "3x + 5y = 7"
diophantine "x^2 - 61y^2 = 1" --limit 5
diophantine "1/x + 1/y + 1/z = 1/2" --domain P
diophantine "x + y + z = 6; x - y = 1"
diophantine "x ≡ 2 (mod 3); x ≡ 3 (mod 5)"
diophantine "frobenius 6 9 20"
diophantine "x^2 + y^2 = 325" --json        # machine-readable output
diophantine                                 # interactive mode
diophantine < equations.txt                 # one equation per line
```

| Option | Meaning |
|---|---|
| `-d, --domain Z\|N\|P` | all integers (default), non-negative, positive |
| `-l, --limit N` | examples shown per infinite family (default 10) |
| `--json` | print the result as JSON |
| `-V, --version` | show the version |

In interactive mode, `:domain P`, `:limit 5` and `:json on` change the settings.

## Library usage

```python
from diophantine import solve, ParseError

r = solve("x^2 - 2y^2 = 1", domain="P")

print(r)                         # full step-by-step report (same as the terminal)
r.status                         # 'none' | 'finite' | 'infinite' | 'identity' | 'info' | 'error'
r.vars                           # ['x', 'y'] — tuple order of every solution
r.solutions                      # finite solutions as tuples
r.families                       # infinite families (parametric / recursive / generator)
r.is_solvable                    # True / False (None when inconclusive: see below)
r.partial                        # True if completeness is not guaranteed
r.proof                          # reason when there are no solutions
r.trace                          # solution steps as a list of strings
list(r.iter_solutions(4))        # [(3, 2), (17, 12), (99, 70), (577, 408)]
r.to_dict()                      # JSON-friendly dict (the trace is under the "steps" key)

try:
    solve("2^x = 8")
except ParseError as e:
    print(e)                     # exponential equations are rejected with a clear message
```

`ParseError` is raised only for input the parser rejects outright. When an equation parses but falls outside what the solvers cover — a non-linear system, or a congruence system whose search space is too large — `solve` still returns a `Result`, with `status = 'error'`, `is_solvable = None` and a note explaining why. `is_solvable` is also `None` for a partial result and for the `'info'` status; it is `False` only when insolvability has been proved.

Each family also exposes `samples(n, domain)` and `describe(variable_names)`. Parametric families have `exprs` (polynomials in the parameters) and `value({"t": 5})`; recursive families have `base` and the forward/backward step matrices `Mf, cf, Mb, cb`.

The number-theory toolkit can be used on its own: `from diophantine import numtheory`, then `numtheory.factorint`, `pell_fundamental`, `sqrt_mod_all`, `crt_pair`, `lll` and so on.

### Input syntax

Implicit multiplication (`3x`, `xy`, `2(x+1)`), powers (`x^2`, `x**2`, `x²`), fractions (`1/x`), indexed variables (`x1`, `x_2`) and the symbols `−`, `·`, `×` are all accepted. Separate system equations with `;` or newlines. Congruences may be written `3x ≡ 4 (mod 7)` or `3x = 4 mod 7`.

## How it works

The pipeline (`core.py`) runs these stages in order:

1. Parsing into an AST (`parser.py`).
2. Normalisation into a canonical integer polynomial, clearing denominators (`poly.py`).
3. A gcd divisibility check.
4. A modular obstruction search over the primes and prime powers up to 64, as far as the budget allows — it stops once modulus^variables would exceed 40 000. If an obstruction is found, it proves there are no solutions.
5. Classification and dispatch to the most specific matching solver.
6. Verification of every solution against the original equation.
7. The domain filter (Z / N / P), including integer intervals for parameters (`intervals.py`).

## Solvers (`diophantine/solvers/`)

| Code | Equation type | Method |
|---|---|---|
| A1/A2 | Linear equation | Extended Euclid; Hermite normal form + LLL; complete enumeration in N/P when coefficients share a sign |
| A3 | Linear system | Hermite normal form via column operations |
| A4 | Congruences | Modular inverse, generalised Chinese remainder theorem |
| A5 | Frobenius coin problem | Dijkstra on the residue graph |
| A6 | x² + y² = z² | Euclid's parametrisation |
| A7 | axy + bx + cy + d = 0 | (ax + c)(ay + b) = bc − ad, divisor enumeration |
| A8 | Egyptian fractions (≥ 3 terms) | Recursion bounded by the largest term |
| A9 | x² + y² = n | Fermat–Euler criterion, Cornacchia, Gaussian integers |
| B1 | x² − Dy² = 1 | Continued fraction of √D |
| B2 | x² − Dy² = −1 | Parity of the period |
| B3 | x² − Dy² = N | LMM (Lagrange–Matthews–Mollin) algorithm |
| B4 | General two-variable quadratic | Δ<0: finite interval; Δ=0: unimodular substitution; Δ>0 square: divisors; Δ>0 non-square: Pell reduction with periodic integrality |
| B5 | Sums of 3, 4, … k squares | Legendre, Lagrange, canonical enumeration, Jacobi check |
| B6 | ax² + by² + cz² = 0 | Normalisation, Legendre criterion, Holzer search, parametrisation |
| B7 | x² + y² + z² = kxyz | Hurwitz's theorem, Vieta (Markov) tree |
| U | One-variable polynomial | Rational root theorem |
| F | Unrecognised type | Bounded search, flagged as a partial result |

The number-theory toolkit (`numtheory.py`) provides the Euclidean algorithm, a Newton `isqrt`, Miller–Rabin, Pollard–Brent, Tonelli–Shanks, Hensel lifting, CRT, continued fractions and LLL. Everything is deterministic: Miller–Rabin uses the first twelve primes as fixed bases, which is a proof below 3.3·10²⁴ and a strong probable-prime test above it, and Pollard–Brent seeds its own generator with the number being factored, so the same input always takes the same path.

## Tests

```bash
python tests/test_solver.py
```

The suite has 129 tests and takes about 3–4 minutes. Each equation's answer is compared with brute-force search inside a box, checking both soundness (every returned solution satisfies the equation) and completeness (no solution in the box is missed). The suite includes 60 randomly generated two-variable quadratic equations, drawn from a fixed seed so that every run is reproducible.

## Limitations

- **Level-3 equations.** Exponential equations (`2^x + 3^y = z^2`) are rejected by the parser. Elliptic curves (`y² = x³ + k`) and Thue equations fall back to bounded search and are flagged as partial.
- **Non-linear systems.** These are not supported: `solve` returns a result with `status = 'error'` and a note, rather than raising.
- **Very large parameters.** For example, a B4 elliptic interval above 5·10⁶ or a parabolic modulus above 10⁷ gives a partial result, and the output says so.

## Licence

MIT — see [LICENSE](LICENSE).

Written by Aytən Süleymanlı.
