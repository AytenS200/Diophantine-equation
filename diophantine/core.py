import re

from . import intervals as iv
from .numtheory import gcd_list, signed_divisors
from .parser import integerize
from .parser import parse_equation, ParseError
from .result import Result, in_domain
from .solvers import (univariate, linear, egyptian, pythagorean, squares, pell, bilinear,
                      markov, legendre, binary_quadratic, fallback, congruence, frobenius)

SOLVERS = [
    univariate,
    linear,
    egyptian,
    pythagorean,
    squares,
    pell,
    bilinear,
    markov,
    legendre,
    binary_quadratic,
]


class Ctx:
    def __init__(self, eq, poly, domain, limit):
        self.eq, self.poly = eq, poly
        self.ip = {m: int(c) for m, c in poly.terms.items()}
        self.vars = poly.variables()
        self.domain, self.limit = domain, limit
        self.trace, self.data = [], {}

    def step(self, msg):
        self.trace.append(msg)

    def result(self, code, title):
        r = Result(self.eq.text, self.vars, self.domain, code, title, self.limit)
        r.trace = self.trace
        r.normal_form = str(self.poly)
        return r


def solve(text, domain="Z", limit=10):
    domain = domain.upper()
    if domain not in ("Z", "N", "P"):
        raise ValueError("domain must be 'Z', 'N' or 'P'")
    raw = text.strip()
    m = re.match(r"^\s*frobenius\s*\(?\s*([-\d\s,]+?)\s*\)?\s*$", raw, re.I)
    if m:
        return frobenius.solve([int(v) for v in re.split(r"[\s,]+", m.group(1)) if v], limit)
    pieces = [p.strip() for p in re.split(r"[;\n]", raw) if p.strip()]
    if not pieces:
        raise ParseError("Empty input")
    if any(re.search(r"\bmod\b", p, re.I) for p in pieces):
        return congruence.solve(pieces, domain, limit)
    eqs = [parse_equation(p) for p in pieces]
    if len(eqs) > 1:
        return linear.solve_system(eqs, domain, limit)
    eq = eqs[0]
    ctx = Ctx(eq, eq.integer_poly(), domain, limit)
    res = _dispatch(ctx)
    _finalize(ctx, res)
    lost = set(eq.variables()) - set(ctx.vars)
    if lost and res.status not in ("identity",):
        res.notes.append(f"{', '.join(sorted(lost))} cancels out of the equation and may take any value.")
    return res


def _dispatch(ctx):
    poly = ctx.poly
    if poly.is_zero():
        r = ctx.result("—", "Eynilik")
        r.status = "identity"
        return r
    if not ctx.vars:
        r = ctx.result("—", "Constant equation")
        return r.none(f"The equation reduces to {poly} = 0, which is false.")
    nonconst = [c for m, c in ctx.ip.items() if m]
    g, c0 = gcd_list(nonconst), ctx.ip.get((), 0)
    if g > 1 and c0 % g:
        ctx.step(f"The gcd of the variable coefficients is {g}, but the constant term {c0} is not divisible by {g}")
        r = ctx.result("M", "Divisibility obstruction")
        return r.none(f"For integer values the variable part is always divisible by {g}, "
                      f"but {-c0} is not.")
    if poly.degree() >= 2 and len(ctx.vars) <= 4:
        obs = fallback.modular_obstruction(ctx.ip, ctx.vars)
        if obs:
            m, vals = obs
            ctx.step(f"All {m}^{len(ctx.vars)} residue combinations modulo {m} were checked")
            r = ctx.result("M", "Modular obstruction")
            return r.none(f"Modulo {m} the left side only takes the residues "
                          f"{{{', '.join(map(str, vals))}}} and is never 0.")
    for S in SOLVERS:
        if S.match(ctx):
            r = S.solve(ctx)
            if r is not None:
                return r
    return fallback.solve(ctx)


def _int_roots(p, t):
    if p.is_zero():
        return None
    cs = [int(c) for c in integerize(p).univariate_coeffs(t)] if not p.is_const() else [1]
    low = next(k for k, c in enumerate(cs) if c)
    roots = {0} if low else set()
    if len(cs) - 1 > low:
        for d in signed_divisors(cs[low]):
            v = 0
            for c in reversed(cs[low:]):
                v = v * d + c
            if v == 0:
                roots.add(d)
    return roots


def _apply_denominators(ctx, r):
    dens = ctx.eq.denominators()
    if not dens:
        return
    kept = []
    for fam in r.families:
        if fam.kind == "param" and len(fam.params) == 1:
            p = fam.params[0]
            sub = dict(zip(r.vars, fam.exprs))
            dead = False
            for D in dens:
                roots = _int_roots(D.substitute(sub), p)
                if roots is None:
                    dead = True
                    break
                fam.excluded |= roots
            if dead:
                continue
        kept.append(fam)
    r.families = kept


def _finalize(ctx, r):
    if r.status == "identity":
        return
    eq, names = ctx.eq, r.vars
    _apply_denominators(ctx, r)

    def ok(t):
        return eq.check(dict(zip(names, t)))

    good, dropped = [], 0
    for t in dict.fromkeys(r.solutions):
        if ok(t):
            good.append(t)
        else:
            dropped += 1
    if dropped:
        r.notes.append(f"{dropped} candidate(s) failed the original equation (e.g. a denominator vanishes) and were discarded.")
    r.solutions = good
    if ctx.domain != "Z" and not r.domain_applied:
        r.solutions = [t for t in r.solutions if in_domain(t, ctx.domain)]
        kept = []
        for fam in r.families:
            if fam.kind == "param" and len(fam.params) == 1 and all(e.degree() <= 2 for e in fam.exprs):
                p = fam.params[0]
                sets = [iv.poly_int_set(e.univariate_coeffs(p) if not e.is_const()
                                        else [e.const_value()], ctx.domain == "P")
                        for e in fam.exprs]
                ivs = iv.intersect_all(sets)
                if not ivs:
                    continue
                if iv.is_bounded(ivs) and iv.size(ivs) <= 100000:
                    for t in iv.iterate(ivs):
                        if t in fam.excluded:
                            continue
                        tup = fam.value({p: t})
                        if ok(tup) and tup not in r.solutions:
                            r.solutions.append(tup)
                    continue
                fam.ranges[p] = ivs
            kept.append(fam)
        r.families = kept
        before = len(r.families)
        r.families = [f for f in r.families if f.kind != "recursive" or f.samples(1, ctx.domain)]
        if len(r.families) < before:
            r.notes.append(f"{before - len(r.families)} famil(y/ies) have no member in the chosen domain "
                           "(signs are constant along the orbit) and are not shown.")
    r.solutions.sort(key=lambda t: (sum(abs(v) for v in t), t))
    for fam in r.families:
        for t in fam.samples(4, "Z"):
            res = ok(t)
            if res is None:
                r.notes.append(f"Some family members (e.g. {t}) make a denominator zero and are excluded.")
                break
            if not res:
                r.notes.append(f"INTERNAL ERROR: {t} does not satisfy the equation - please report this.")
                break
    r.finalize_status()
    if r.status == "none" and not r.proof:
        r.proof = "No solution satisfies the domain constraint."
