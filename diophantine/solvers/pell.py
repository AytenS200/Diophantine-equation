from ..numtheory import (sqrt_cf, pell_fundamental, pell_neg_fundamental, isqrt,
                         factorint, sqrt_mod_all, is_square)
from ..poly import M
from ..result import RecursiveFamily

CODE = "B1–B3"


def _lmm(D, s, m, Q0, P0, neg):
    P, Q = P0, Q0
    A1, A2, B1, B2 = 1, 0, 0, 1
    seen = set()
    for _ in range(2_000_000):
        if (P, Q) in seen:
            return None
        seen.add((P, Q))
        a = (P + s) // Q if Q > 0 else (P + s + 1) // Q
        A1, A2 = a * A1 + A2, A1
        B1, B2 = a * B1 + B2, B1
        G = Q0 * A1 - P0 * B1
        P = a * Q - P
        Q = (D - P * P) // Q
        if Q in (1, -1):
            val = G * G - D * B1 * B1
            if val == m:
                return G, B1
            if val == -m:
                if neg:
                    t, u = neg
                    return G * t + D * B1 * u, G * u + B1 * t
                return None
    return None


def same_class(a, b, D, N):
    n = abs(N)
    return (a[0] * b[0] - D * a[1] * b[1]) % n == 0 and (a[0] * b[1] - b[0] * a[1]) % n == 0


def gen_pell_reps(D, N):
    neg = pell_neg_fundamental(D)
    s = isqrt(D)
    fs = [1]
    for p, e in factorint(N).items():
        fs = [f * p ** k for f in fs for k in range(e // 2 + 1)]
    reps = []
    for f in sorted(fs):
        m = N // (f * f)
        am = abs(m)
        for z in ([0] if am == 1 else sqrt_mod_all(D % am, am)):
            if 2 * z > am:
                z -= am
            sol = _lmm(D, s, m, am, z, neg)
            if sol:
                cand = (f * sol[0], f * sol[1])
                assert cand[0] ** 2 - D * cand[1] ** 2 == N
                if not any(same_class(cand, r, D, N) for r in reps):
                    reps.append(cand)
    return reps


def _orient(ctx):
    x, y = ctx.vars
    if set(ctx.ip) - {M((x, 2)), M((y, 2)), ()}:
        return None
    cx, cy = ctx.ip.get(M((x, 2)), 0), ctx.ip.get(M((y, 2)), 0)
    if cx * cy >= 0:
        return None
    const = ctx.ip.get((), 0)
    for u, cu, v, cv in ((x, cx, y, cy), (y, cy, x, cx)):
        if abs(cu) == 1:
            D, N = -cv * cu, -const * cu
            if D > 0 and not is_square(D):
                return u, v, D, N
    return None


def match(ctx):
    if len(ctx.vars) != 2 or ctx.poly.degree() != 2:
        return False
    o = _orient(ctx)
    if o:
        ctx.data["pell"] = o
    return bool(o)


def _perm(Mx):
    return [[Mx[1][1], Mx[1][0]], [Mx[0][1], Mx[0][0]]]


def solve(ctx):
    u, v, D, N = ctx.data["pell"]
    kind = {1: "B1", -1: "B2"}.get(N, "B3")
    title = {"B1": "Pell equation", "B2": "Negative Pell equation",
             "B3": "Generalised Pell equation"}[kind]
    r = ctx.result(kind, f"{title}: {u}² − {D}{v}² = {N}")
    swap = ctx.vars.index(u) == 1

    def tup(a, b):
        return (b, a) if swap else (a, b)

    if N == 0:
        ctx.step(f"√{D} is irrational ⇒ {u}² = {D}{v}² only when {u} = {v} = 0")
        r.solutions = [(0, 0)]
        return r
    a0, per = sqrt_cf(D)
    shown = ", ".join(map(str, per[:30])) + (", …" if len(per) > 30 else "")
    ctx.step(f"√{D} = [{a0}; ({shown})] - continued fraction with period length L = {len(per)}")
    x1, y1 = pell_fundamental(D)
    ctx.step(f"Fundamental solution of x² − {D}y² = 1 (L {'even' if len(per) % 2 == 0 else 'odd'} ⇒ "
             f"convergent {'p_{L−1}' if len(per) % 2 == 0 else 'p_{2L−1}'}): ε = {x1} + {y1}√{D}")
    Mf, Mb = [[x1, D * y1], [y1, x1]], [[x1, -D * y1], [-y1, x1]]
    if swap:
        Mf, Mb = _perm(Mf), _perm(Mb)
    if N == 1:
        bases = [(1, 0), (-1, 0)]
        ctx.step("All solutions are ±ε^k (k ∈ ℤ): multiplying (x + y√D) by ε gives a new solution")
    elif N == -1:
        neg = pell_neg_fundamental(D)
        if not neg:
            return r.none(f"The period of √{D} has even length L = {len(per)}; in that case "
                          f"x² − {D}y² = −1 has no solutions.")
        ctx.step(f"L is odd ⇒ fundamental solution p_(L−1)/q_(L−1) = ({neg[0]}, {neg[1]}); "
                 f"all solutions are ±η^(2k+1) with η² = ε")
        bases = [neg, (-neg[0], -neg[1])]
    else:
        ctx.step(f"LMM algorithm: for every f with f² | {N} and every z with z² ≡ {D} (mod |{N}/f²|), "
                 f"expand (z + √{D})/|m| as a continued fraction until |Q| = 1")
        reps = gen_pell_reps(D, N)
        if not reps:
            return r.none(f"The LMM algorithm found no solution class: x² − {D}y² = {N} "
                          f"has no integer solutions.")
        ctx.step(f"Found {len(reps)} solution class(es), representatives: "
                 + ", ".join(f"({a}, {b})" for a, b in reps))
        bases = reps + [(-a, -b) for a, b in reps]
        ctx.step("Each class gives an infinite family ±(x0 + y0√D)·ε^k")
    for b in bases:
        r.families.append(RecursiveFamily(tup(*b), Mf, [0, 0], Mb, [0, 0]))
    return r
