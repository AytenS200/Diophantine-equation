from ..numtheory import signed_divisors
from ..poly import M

CODE, TITLE = "U", "Polynomial equation in one variable"


def match(ctx):
    return len(ctx.vars) == 1 and ctx.poly.degree() >= 2


def solve(ctx):
    x = ctx.vars[0]
    deg = ctx.poly.degree()
    cs = [ctx.ip.get(M((x, k)), 0) for k in range(deg + 1)]
    r = ctx.result(CODE, f"{TITLE} (degree {deg})")
    low = next(k for k, c in enumerate(cs) if c)
    roots = set()
    if low:
        roots.add(0)
        ctx.step(f"{x}^{low} is a common factor ⇒ {x} = 0 is a root")
    c0, cn = cs[low], cs[-1]
    ctx.step(f"Rational root theorem: every integer root divides {c0}")
    for d in signed_divisors(c0):
        v = 0
        for c in reversed(cs[low:]):
            v = v * d + c
        if v == 0:
            roots.add(d)
    ctx.step(f"{len(signed_divisors(c0))} candidates checked with Horner's scheme")
    if not roots:
        return r.none(f"No divisor of {c0} is a root.")
    r.solutions = [(v,) for v in sorted(roots)]
    return r
