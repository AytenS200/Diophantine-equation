from ..numtheory import signed_divisors
from ..poly import Poly, M
from ..result import ParamFamily

CODE, TITLE = "A7", "Bilinear equation (factoring trick)"


def match(ctx):
    if len(ctx.vars) != 2 or ctx.poly.degree() != 2:
        return False
    x, y = ctx.vars
    allowed = {M((x, 1), (y, 1)), M((x, 1)), M((y, 1)), ()}
    return M((x, 1), (y, 1)) in ctx.ip and set(ctx.ip) <= allowed


def bilinear_core(a, b, c, d, step):
    N = b * c - a * d
    X, Y = Poly.var("x"), Poly.var("y")
    step(f"Multiply by {a} and factor: ({a * X + c})({a * Y + b}) = "
         f"bc - ad = {N}")
    sols, fams = set(), []
    if N != 0:
        ds = signed_divisors(N)
        step(f"Checking the {len(ds)} signed divisors D of {N}: ax + c = D, ay + b = N/D")
        for D in ds:
            if (D - c) % a == 0 and (N // D - b) % a == 0:
                sols.add(((D - c) // a, (N // D - b) // a))
        if not sols:
            return sols, fams, (f"No divisor pair of ({a}x + {c})({a}y + {b}) = {N} "
                                f"yields integer x, y.")
        return sols, fams, None
    t = Poly.var("t")
    step("The right side is 0 ⇒ one factor must vanish: ax + c = 0 or ay + b = 0")
    if c % a == 0:
        fams.append(ParamFamily([Poly.const(-c // a), t], ["t"]))
    if b % a == 0:
        fams.append(ParamFamily([t, Poly.const(-b // a)], ["t"]))
    if not fams:
        return sols, fams, f"Neither {-c}/{a} nor {-b}/{a} is an integer."
    return sols, fams, None


def solve(ctx):
    x, y = ctx.vars
    a = ctx.ip[M((x, 1), (y, 1))]
    b, c, d = ctx.ip.get(M((x, 1)), 0), ctx.ip.get(M((y, 1)), 0), ctx.ip.get((), 0)
    r = ctx.result(CODE, TITLE)
    ctx.step(f"Form axy + bx + cy + d = 0 with a={a}, b={b}, c={c}, d={d}  ({x} → x, {y} → y)")
    sols, fams, proof = bilinear_core(a, b, c, d, ctx.step)
    if proof:
        return r.none(proof)
    r.solutions, r.families = sorted(sols), fams
    return r
