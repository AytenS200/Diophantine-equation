from fractions import Fraction

from ..numtheory import gcd, gcd_list, lcm, isqrt, factorint, squarefree_decomp, is_qr_prime
from ..poly import Poly, M
from ..result import GeneratorFamily, in_domain

CODE, TITLE = "B6", "Legendre equation ax² + by² + cz² = 0"


def match(ctx):
    if len(ctx.vars) != 3 or ctx.poly.degree() != 2:
        return False
    sq = {M((v, 2)) for v in ctx.vars}
    return set(ctx.ip) == sq


def solve(ctx):
    names = ctx.vars
    coef = [ctx.ip[M((v, 2))] for v in names]
    r = ctx.result(CODE, TITLE)
    ctx.step(f"Coefficients: ({', '.join(map(str, coef))})")
    scale = [Fraction(1)] * 3
    while True:
        changed = False
        g = gcd_list(coef)
        if g > 1:
            coef = [c // g for c in coef]
        for i in range(3):
            sf, s = squarefree_decomp(coef[i])
            if s > 1:
                coef[i] = sf
                scale[i] /= s
                changed = True
        for i in range(3):
            for j in range(i + 1, 3):
                g2 = gcd(coef[i], coef[j])
                if g2 > 1 and not changed:
                    p = min(factorint(g2))
                    k = 3 - i - j
                    coef[i] //= p
                    coef[j] //= p
                    coef[k] *= p
                    scale[k] *= p
                    changed = True
        if not changed:
            break
    ctx.step(f"Normalised form (squarefree, pairwise coprime): "
             f"{coef[0]}X² + {coef[1]}Y² + {coef[2]}Z² = 0")
    trivial = r
    if all(c > 0 for c in coef) or all(c < 0 for c in coef):
        r.solutions = [(0, 0, 0)]
        ctx.step("All coefficients have the same sign ⇒ the left side vanishes only at zero")
        return r
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        mod = abs(coef[i])
        val = -coef[j] * coef[k]
        for p in factorint(mod):
            if not is_qr_prime(val, p):
                ctx.step(f"Legendre criterion: {val} is not a quadratic residue modulo {p}")
                r.solutions = [(0, 0, 0)]
                r.notes.append(f"Only the trivial solution: {val} is not a quadratic residue mod {p} "
                               f"(a prime divisor of |{coef[i]}|), so no non-trivial solution exists.")
                return r
    ctx.step("Legendre criterion holds ⇒ a non-trivial solution exists")
    A = coef
    bounds = [isqrt(abs(A[1] * A[2])), isqrt(abs(A[0] * A[2])), isqrt(abs(A[0] * A[1]))]
    k = max(range(3), key=lambda t: bounds[t])
    i, j = [t for t in range(3) if t != k]
    ctx.step(f"Holzer bounds: |X| ≤ {bounds[0]}, |Y| ≤ {bounds[1]}, |Z| ≤ {bounds[2]} - searching")
    if (bounds[i] + 1) * (bounds[j] + 1) > 6_000_000:
        r.partial = True
        r.solutions = [(0, 0, 0)]
        r.notes.append("A solution exists (Legendre), but the coefficients are too large for the Holzer search.")
        return r
    P0 = None
    for vi in range(bounds[i] + 1):
        for vj in range(bounds[j] + 1):
            if vi == 0 and vj == 0:
                continue
            num = -(A[i] * vi * vi + A[j] * vj * vj)
            if num % A[k]:
                continue
            w = num // A[k]
            if w >= 0 and isqrt(w) ** 2 == w:
                P0 = [0, 0, 0]
                P0[i], P0[j], P0[k] = vi, vj, isqrt(w)
                break
        if P0:
            break
    assert P0, "Holzer's theorem guarantees a solution here"
    L = 1
    for s in scale:
        L = lcm(L, s.denominator)

    def to_orig(vec):
        return [scale[t] * vec[t] * L for t in range(3)]

    o0 = [int(v) for v in to_orig(P0)]
    g0 = gcd_list(o0)
    o0 = tuple(v // g0 for v in o0)
    ctx.step(f"Solution of the normalised form: {tuple(P0)} ⇒ in the original variables: {o0}")
    k0 = next(t for t in range(3) if P0[t] != 0)
    others = [t for t in range(3) if t != k0]
    mP, nP = Poly.var("m"), Poly.var("n")
    w = [Poly(), Poly(), Poly()]
    w[others[0]], w[others[1]] = mP, nP
    Qw = sum((A[t] * w[t] * w[t] for t in range(3)), Poly())
    Bw = sum((A[t] * P0[t] * w[t] for t in range(3)), Poly())
    vec = [Qw * P0[t] - 2 * Bw * w[t] for t in range(3)]
    V = [vec[t] * (scale[t] * L) for t in range(3)]
    ctx.step("Every line through P0 with rational slope meets the conic in one more point: "
             "v = Q(w)·P0 − 2B(P0, w)·w - this gives every rational solution")

    def prim(m, n):
        vals = [int(e.evaluate({"m": m, "n": n})) for e in V]
        g = gcd_list(vals)
        return tuple(abs(v) // g for v in vals) if g else None

    def gen(cnt, domain):
        found = set()
        if domain != "P":
            found.add((0, 0, 0))
        found.add(tuple(abs(v) for v in o0))
        R = 1
        while len(found) < cnt + 5 and R < 60:
            for m in range(-R, R + 1):
                for n in range(-R, R + 1):
                    if max(abs(m), abs(n)) == R and gcd(m, n) == 1:
                        t = prim(m, n)
                        if t:
                            found.add(t)
            R += 1
        res = sorted((t for t in found if in_domain(t, domain)), key=lambda t: (max(t), t))
        return res[:cnt]

    lines = [f"{v} = {e}" for v, e in zip(names, V)]
    lines.append("where gcd(m, n) = 1; divide the triple by its gcd to get a primitive solution")
    lines.append("all solutions: integer multiples k of primitive solutions; coordinate signs are arbitrary")
    r.families.append(GeneratorFamily(lines, gen))
    r.domain_applied = True
    return r
