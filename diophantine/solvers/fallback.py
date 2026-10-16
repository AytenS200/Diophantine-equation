from itertools import product

from ..numtheory import isqrt

CODE, TITLE = "F", "Unrecognised type - bounded search"

MODULI = (2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23, 25, 27, 29, 31, 32,
          37, 41, 43, 47, 49, 53, 59, 61, 64)


def _compile(ip, names):
    idx = {v: i for i, v in enumerate(names)}
    return [(c, [(idx[v], e) for v, e in m]) for m, c in ip.items()]


def modular_obstruction(ip, names, budget=40000):
    n = len(names)
    comp = _compile(ip, names)
    for m in MODULI:
        if m ** n > budget:
            break
        pw = [[pow(r, e, m) for r in range(m)] for e in range(0, 1 + max(
            (e for _, mon in comp for _, e in mon), default=1))]
        seen, hit = set(), False
        for vals in product(range(m), repeat=n):
            s = 0
            for c, mon in comp:
                t = c
                for i, e in mon:
                    t *= pw[e][vals[i]]
                s += t
            s %= m
            if s == 0:
                hit = True
                break
            seen.add(s)
        if not hit:
            return m, sorted(seen)
    return None


def solve(ctx):
    r = ctx.result(CODE, TITLE)
    r.partial = True
    names = ctx.vars
    degs = {v: ctx.poly.degree_in(v) for v in names}
    sv = min(names, key=lambda v: degs[v]) if len(names) > 1 else names[0]
    sv = sv if degs[sv] <= 2 else None
    loop = [v for v in names if v != sv]
    ctx.step("This equation does not belong to any level-1 or level-2 family")
    ctx.step("A modular obstruction was searched for and not found")
    groups = {}
    for m, c in ctx.ip.items():
        e = dict(m).get(sv, 0) if sv else 0
        rest = tuple((v, ee) for v, ee in m if v != sv)
        groups.setdefault(e, []).append((c, rest))
    idx = {v: i for i, v in enumerate(loop)}
    budget = 300000
    B = max(1, int((budget ** (1 / max(1, len(loop))) - 1) / 2)) if loop else 0
    if ctx.domain == "Z":
        rng = range(-B, B + 1)
    else:
        rng = range(1 if ctx.domain == "P" else 0, 2 * B + 2)
    ctx.step(f"Bounded search: {', '.join(loop) or '—'} ∈ [{rng.start}, {rng.stop - 1}]"
             + (f", {sv} is solved exactly each time (degree ≤ 2)" if sv else ""))
    sols, inf_note = set(), False
    for vals in product(rng, repeat=len(loop)):
        cs = {}
        for e, terms in groups.items():
            s = 0
            for c, rest in terms:
                t = c
                for v, ee in rest:
                    t *= vals[idx[v]] ** ee
                s += t
            cs[e] = s
        if sv is None:
            if cs.get(0, 0) == 0:
                sols.add(tuple(vals))
            continue
        c2, c1, c0 = cs.get(2, 0), cs.get(1, 0), cs.get(0, 0)
        roots = []
        if c2 == 0:
            if c1 == 0:
                if c0 == 0:
                    inf_note = True
            elif c0 % c1 == 0:
                roots = [-c0 // c1]
        else:
            disc = c1 * c1 - 4 * c2 * c0
            if disc >= 0 and isqrt(disc) ** 2 == disc:
                s = isqrt(disc)
                roots = [num // (2 * c2) for num in {-c1 + s, -c1 - s} if num % (2 * c2) == 0]
        for root in roots:
            full = dict(zip(loop, vals))
            full[sv] = root
            sols.add(tuple(full[v] for v in names))
    r.solutions = sorted(sols)
    if inf_note:
        r.notes.append(f"For some values {sv} is arbitrary (an infinite family).")
    r.notes.append("Solutions were only searched within the stated range; this method cannot decide "
                   "whether larger solutions exist (e.g. level 3: elliptic curves, Thue equations).")
    if not sols:
        r.proof = "No solution in the searched range (this does not prove that none exist)."
    return r
