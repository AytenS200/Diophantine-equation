from ..result import GeneratorFamily, in_domain
from ..poly import M
from .squares import two_squares_all

CODE, TITLE = "A6", "Pythagorean equation x² + y² = z²"


def match(ctx):
    if len(ctx.vars) != 3 or ctx.poly.degree() != 2 or len(ctx.ip) != 3:
        return False
    sq = {v: ctx.ip.get(M((v, 2)), 0) for v in ctx.vars}
    vals = sorted(sq.values())
    if vals == [-1, 1, 1]:
        z = next(v for v in ctx.vars if sq[v] == -1)
    elif vals == [-1, -1, 1]:
        z = next(v for v in ctx.vars if sq[v] == 1)
    else:
        return False
    ctx.data["hyp"] = z
    return True


def solve(ctx):
    z = ctx.data["hyp"]
    names = ctx.vars
    zi = names.index(z)
    legs = [i for i in range(3) if i != zi]
    x, y = names[legs[0]], names[legs[1]]
    r = ctx.result(CODE, f"Pythagorean equation {x}² + {y}² = {z}²")
    ctx.step(f"In a primitive solution exactly one of {x}, {y} is even; for it, "
             f"{z}² − {x}² = ({z} − {x})({z} + {x}) forces both factors to be squares (up to 2)")
    ctx.step("This yields Euclid's parametrisation (m² − n², 2mn, m² + n²), "
             "m > n ≥ 0, gcd(m, n) = 1, m − n odd; then scale by k")

    def assemble(p, q, c):
        t = [0, 0, 0]
        t[legs[0]], t[legs[1]], t[zi] = p, q, c
        return tuple(t)

    def gen(n, domain):
        out = []
        if domain != "P":
            out.append((0, 0, 0))
        c = 0
        while len(out) < n and c < 100000:
            c += 1
            batch = []
            for p, q in two_squares_all(c * c):
                for cs in ((c, -c) if domain == "Z" else (c,)):
                    t = assemble(p, q, cs)
                    if in_domain(t, domain):
                        batch.append(t)
            out += sorted(batch, key=lambda t: (sum(v < 0 for v in t), 0 in t, t))
        return out[:n]

    lines = [
        f"{x} = k(m² − n²),  {y} = 2kmn,  {z} = k(m² + n²)   (or with {x} and {y} swapped)",
        "m > n ≥ 0, gcd(m, n) = 1, m − n odd, k ≥ 0; every coordinate may change sign",
    ]
    r.families.append(GeneratorFamily(lines, gen))
    r.domain_applied = True
    return r
