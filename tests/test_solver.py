import random
import sys
import time
from itertools import product

sys.path.insert(0, ".")
from diophantine import solve
from diophantine.parser import parse_equation
from diophantine.solvers.squares import expand_variants


def brute(eq, names, B, domain="Z"):
    lo = {"Z": -B, "N": 0, "P": 1}[domain]
    poly = eq.integer_poly()
    if len(names) == 2 and poly.degree_in(names[1]) <= 2:
        from diophantine.numtheory import isqrt
        x, y = names
        groups = {}
        for m, c in poly.terms.items():
            d = dict(m)
            groups.setdefault(d.get(y, 0), []).append((int(c), d.get(x, 0)))
        out = set()
        for xv in range(lo, B + 1):
            cs = [sum(c * xv ** ex for c, ex in groups.get(k, [])) for k in range(3)]
            c0, c1, c2 = cs
            cands = []
            if c2 == 0:
                if c1 == 0:
                    cands = range(lo, B + 1) if c0 == 0 else []
                elif c0 % c1 == 0:
                    cands = [-c0 // c1]
            else:
                disc = c1 * c1 - 4 * c2 * c0
                if disc >= 0 and isqrt(disc) ** 2 == disc:
                    s = isqrt(disc)
                    cands = [n // (2 * c2) for n in (-c1 + s, -c1 - s) if n % (2 * c2) == 0]
            for yv in cands:
                if lo <= yv <= B and eq.check({x: xv, y: yv}):
                    out.add((xv, yv))
        return out
    return {t for t in product(range(lo, B + 1), repeat=len(names))
            if eq.check(dict(zip(names, t)))}


def produced(res, B, domain):
    out = set()
    sols = res.solutions
    if res.canonical:
        sols = [v for t in sols for v in expand_variants(t, domain)]
    out |= set(sols)
    for fam in res.families:
        if fam.kind == "param":
            if len(fam.params) == 1:
                p = fam.params[0]
                for t in range(-4000, 4001):
                    if t in fam.excluded:
                        continue
                    out.add(fam.value({p: t}))
            else:
                R = 12
                for vals in product(range(-R, R + 1), repeat=len(fam.params)):
                    out.add(fam.value(dict(zip(fam.params, vals))))
        elif fam.kind == "recursive":
            out |= set(fam.terms(40))
        else:
            out |= set(fam.samples(20000, domain))
    return {t for t in out if all(abs(v) <= B for v in t)}


FAILS = []


def check(text, B, domain="Z", complete=True):
    eq = parse_equation(text)
    res = solve(text, domain=domain)
    names = res.vars
    br = brute(eq, names, B, domain) if names else set()
    pr = produced(res, B, domain)
    bad = [t for t in pr if not eq.check(dict(zip(names, t)))]
    missing = br - pr if complete else set()
    if res.status == "none" and br:
        missing = br
    ok = not bad and not missing and not any("INTERNAL ERROR" in n for n in res.notes)
    tag = "OK  " if ok else "FAIL"
    print(f"{tag} [{res.code:6}] {text:42} domain={domain} status={res.status:8} in_box={len(br)}")
    if not ok:
        FAILS.append(text)
        if bad:
            print("     WRONG solutions:", sorted(bad)[:8])
        if missing:
            print("     MISSING:", sorted(missing)[:8])
        for n in res.notes:
            print("     note:", n)


CASES = [
    ("3x + 5y = 7", 60, "Z"), ("6x + 9y = 4", 30, "Z"), ("-4x + 6y = 10", 60, "Z"),
    ("3x + 5y = 47", 47, "N"), ("3x + 5y = 47", 47, "P"), ("6x + 9y + 20z = 43", 43, "N"),
    ("6x + 9y + 20z = 44", 44, "N"), ("2x + 3y - 5z = 1", 12, "Z"), ("4x+6y+10z=8", 10, "Z"),
    ("x^3 - 6x^2 + 11x - 6 = 0", 50, "Z"), ("2x^3 + x^2 = 0", 50, "Z"),
    ("xy = 12", 40, "Z"), ("xy + 2x + 3y = 20", 80, "Z"), ("2xy - x + 3y = 7", 80, "Z"),
    ("1/x + 1/y = 1/6", 100, "Z"), ("1/x + 1/y = 1/6", 100, "P"), ("x/y = 2", 20, "Z"),
    ("1/x + 1/y + 1/z = 1", 45, "P"), ("1/x + 1/y + 1/z = 1/2", 45, "P"),
    ("2/x + 3/y + 1/z = 2", 40, "P"),
    ("x^2 + y^2 = 325", 30, "Z"), ("x^2 + y^2 = 21", 10, "Z"), ("x^2 + y^2 = 0", 5, "Z"),
    ("x^2 + y^2 = 2", 5, "Z"), ("x^2 + y^2 = 1105", 40, "N"),
    ("x^2 + y^2 + z^2 = 27", 8, "Z"), ("x^2 + y^2 + z^2 = 28", 8, "Z"),
    ("a^2 + b^2 + c^2 + d^2 = 30", 6, "Z"), ("a^2 + b^2 + c^2 + d^2 = 30", 6, "P"),
    ("x^2 + y^2 = z^2", 25, "Z"), ("x^2 + y^2 = z^2", 30, "P"),
    ("x^2 - 2y^2 = 1", 1000, "Z"), ("x^2 - 61y^2 = 1", 1000, "Z"), ("x^2 - 2y^2 = -1", 1000, "Z"),
    ("x^2 - 3y^2 = -1", 300, "Z"), ("x^2 - 34y^2 = -1", 2000, "Z"),
    ("x^2 - 5y^2 = -4", 1000, "Z"), ("x^2 - 7y^2 = 9", 2000, "Z"), ("x^2 - 6y^2 = 10", 2000, "Z"),
    ("y^2 - 13x^2 = 12", 3000, "Z"), ("x^2 - 10y^2 = -9", 3000, "Z"), ("x^2 - 2y^2 = 0", 100, "Z"),
    ("x^2 - 11y^2 = 25", 3000, "Z"), ("x^2 - 7y^2 = 9", 2000, "P"),
    ("x^2 + xy + y^2 = 49", 20, "Z"), ("2x^2 + 3y^2 - x + y = 100", 20, "Z"),
    ("x^2 - 2xy + y^2 + x - y = 6", 300, "Z"), ("x^2 + 2xy + y^2 - 4x + 2y = 3", 300, "Z"),
    ("4x^2 - 12xy + 9y^2 + 2x + 5y = 7", 300, "Z"), ("x^2 + y = 5", 300, "Z"),
    ("y^2 = 3x + 1", 300, "Z"), ("x^2 - y^2 = 45", 60, "Z"), ("x^2 - 4y^2 + x = 0", 200, "Z"),
    ("x^2 - y^2 + 2x - 4y = 3", 200, "Z"),
    ("x^2 + 3xy + y^2 = 11", 2000, "Z"), ("2x^2 - 7y^2 + x + 3y = 5", 3000, "Z"),
    ("a^2 + b^2 - 3ab - 3 = 0", 2000, "Z"),
    ("x^2 + y^2 + 1 = 3xy", 2000, "Z"),
    ("3x^2 - 5xy + y^2 + 2x - 3y = 10", 3000, "Z"),
    ("x^2 + y^2 + z^2 = 3xyz", 40, "P"), ("x^2 + y^2 + z^2 = xyz", 40, "P"),
    ("x^2 + y^2 + z^2 = 5xyz", 20, "Z"),
    ("y^2 = x^3 + 7", 30, "Z"),
]


def legendre_cases():
    for text, B in [("x^2 + y^2 - 3z^2 = 0", 30), ("2x^2 + 3y^2 - 5z^2 = 0", 25),
                    ("x^2 + y^2 - 2z^2 = 0", 30), ("3x^2 - 7y^2 + 11z^2 = 0", 20),
                    ("5x^2 - 6y^2 + 15z^2 = 0", 20), ("x^2 + 2y^2 - 6z^2 = 0", 30)]:
        eq = parse_equation(text)
        res = solve(text)
        names = res.vars
        br = {t for t in brute(eq, names, B) if any(t)}
        smp = set()
        for fam in res.families:
            smp |= set(fam.samples(30, "Z"))
        bad = [t for t in smp if not eq.check(dict(zip(names, t)))]
        claims_nontrivial = bool(res.families)
        ok = not bad and claims_nontrivial == bool(br)
        print(f"{'OK  ' if ok else 'FAIL'} [B6    ] {text:42} nontrivial={claims_nontrivial} "
              f"in_box={len(br)}")
        if not ok:
            FAILS.append(text)


def random_bq(n=60, seed=1):
    rnd = random.Random(seed)
    for _ in range(n):
        a, b, c, d, e, f = (rnd.randint(-6, 6) for _ in range(6))
        if a == b == c == 0:
            continue
        text = f"{a}x^2 + {b}xy + {c}y^2 + {d}x + {e}y + {f} = 0".replace("+ -", "- ")
        check(text, 150)


if __name__ == "__main__":
    t0 = time.time()
    for text, B, dom in CASES:
        check(text, B, dom, complete=not text.startswith("y^2 = x^3"))
    legendre_cases()
    print("--- random two-variable quadratic equations ---")
    random_bq()
    print(f"\n{len(FAILS)} failed test(s); time {time.time() - t0:.1f} s")
    sys.exit(1 if FAILS else 0)
