from fractions import Fraction

from ..numtheory import egcd, gcd, isqrt, is_square, signed_divisors, pell_fundamental
from ..intervals import poly_int_set
from ..poly import Poly, M
from ..result import ParamFamily, RecursiveFamily
from .bilinear import bilinear_core
from .pell import gen_pell_reps

CODE, TITLE = "B4", "General two-variable quadratic equation"


def match(ctx):
    return len(ctx.vars) == 2 and ctx.poly.degree() == 2


class _Out:
    def __init__(self):
        self.sols, self.fams, self.proof, self.partial, self.notes = set(), [], None, False, []
        self.kind = ""


def linear2(A, B, C):
    if A == 0 and B == 0:
        return None
    g, u, v = egcd(A, B)
    if C % g:
        return None
    k, t = C // g, Poly.var("t")
    return ParamFamily([u * k + (B // g) * t, v * k - (A // g) * t], ["t"])


def _swap(out):
    o = _Out()
    o.sols = {(y, x) for x, y in out.sols}
    o.proof, o.partial, o.notes, o.kind = out.proof, out.partial, out.notes, out.kind
    for f in out.fams:
        if f.kind == "param":
            o.fams.append(ParamFamily(list(reversed(f.exprs)), f.params, f.conditions, f.ranges))
        else:
            P = lambda Mx: [[Mx[1][1], Mx[1][0]], [Mx[0][1], Mx[0][0]]]
            o.fams.append(RecursiveFamily((f.base[1], f.base[0]), P(f.Mf), f.cf[::-1],
                                          P(f.Mb), f.cb[::-1]))
    return o


def solve_bq(a, b, c, d, e, f, step):
    if a == 0 and c == 0:
        o = _Out()
        o.kind = "biliniyar"
        o.sols, o.fams, o.proof = bilinear_core(b, d, e, f, step)
        return o
    if a == 0:
        step("a = 0, so the roles of x and y are swapped")
        return _swap(solve_bq(c, b, a, e, d, f, step))
    D = b * b - 4 * a * c
    step(f"Diskriminant Δ = b² − 4ac = {b}² − 4·{a}·{c} = {D}")
    if D < 0:
        return _elliptic(a, b, c, d, e, f, D, step)
    if D == 0:
        return _parabolic(a, b, c, d, e, f, step)
    if is_square(D):
        return _hyper_square(a, b, c, d, e, f, D, step)
    return _hyper_pell(a, b, c, d, e, f, D, step)


def _elliptic(a, b, c, d, e, f, D, step):
    o = _Out()
    o.kind = "elliptik (Δ < 0)"
    q = [e * e - 4 * c * f, 2 * (b * e - 2 * c * d), D]
    step(f"The discriminant of the quadratic in y must satisfy {D}x² + {q[1]}x + {q[0]} ≥ 0; "
         f"the leading coefficient is negative ⇒ x lies in a finite interval")
    ivs = poly_int_set(q, strict=False)
    if not ivs:
        o.proof = "The discriminant is negative for every integer x - the curve has no real points."
        return o
    lo, hi = ivs[0][0], ivs[-1][1]
    step(f"x ∈ [{lo}, {hi}] - for each x check whether the discriminant is a perfect square")
    if hi - lo > 5_000_000:
        o.partial = True
        o.notes.append(f"The interval is very large ({hi - lo + 1}); only the first 5·10⁶ values were checked.")
        hi = lo + 5_000_000
    for x in range(lo, hi + 1):
        B_, C_ = b * x + e, a * x * x + d * x + f
        disc = B_ * B_ - 4 * c * C_
        if disc < 0:
            continue
        s = isqrt(disc)
        if s * s != disc:
            continue
        for num in (-B_ + s, -B_ - s):
            if num % (2 * c) == 0:
                o.sols.add((x, num // (2 * c)))
    if not o.sols:
        o.proof = f"No x in [{lo}, {hi}] gives an integer y."
    return o


def _parabolic(a, b, c, d, e, f, step):
    o = _Out()
    o.kind = "parabolik (Δ = 0)"
    g = gcd(a, c) * (1 if a > 0 else -1)
    alpha, gamma = isqrt(a // g), isqrt(c // g)
    beta = gamma if b * g >= 0 else -gamma
    assert g * alpha ** 2 == a and g * beta ** 2 == c and 2 * g * alpha * beta == b
    _, u, v = egcd(alpha, beta)
    p, h = d * u + e * v, e * alpha - d * beta
    step(f"The quadratic part is a perfect square: {a}x² + {b}xy + {c}y² = {g}·({alpha}x + {beta}y)²")
    step(f"Unimodular substitution: w = {alpha}x + {beta}y, z = {-v}x + {u}y  "
         f"⇒  x = {u}w − {beta}z, y = {v}w + {alpha}z")
    step(f"Equation: {g}w² + {p}w + {f} + {h}·z = 0")
    t = Poly.var("t")
    if h == 0:
        disc = p * p - 4 * g * f
        ws = set()
        if disc >= 0 and is_square(disc):
            s = isqrt(disc)
            ws = {num // (2 * g) for num in (-p + s, -p - s) if num % (2 * g) == 0}
        step(f"The coefficient of z is 0 ⇒ w must be an integer root of {g}w² + {p}w + {f} = 0")
        if not ws:
            o.proof = f"{g}w² + {p}w + {f} = 0 has no integer root."
        for w in sorted(ws):
            o.fams.append(ParamFamily([u * w - beta * t, v * w + alpha * t], ["t"]))
        return o
    H = abs(h)
    step(f"For z to be an integer we need {g}w² + {p}w + {f} ≡ 0 (mod {H})")
    if H > 10_000_000:
        o.partial = True
        o.notes.append("The modulus is too large; the residue search was skipped.")
        return o
    res = [r0 for r0 in range(H) if (g * r0 * r0 + p * r0 + f) % H == 0]
    if not res:
        o.proof = f"The congruence {g}w² + {p}w + {f} ≡ 0 (mod {H}) has no solution."
        return o
    step(f"Admissible residues: w ≡ {', '.join(map(str, res[:20]))}{' …' if len(res) > 20 else ''} "
         f"(mod {H}); for each, substitute w = r + {H}t")
    for r0 in res:
        w = r0 + H * t
        z = (g * w * w + p * w + f) * Fraction(-1, h)
        o.fams.append(ParamFamily([u * w - beta * z, v * w + alpha * z], ["t"]))
    return o


def _back(a, b, d, D, h0, X, Y):
    if (Y + h0) % D:
        return None
    y = (Y + h0) // D
    num = X - b * y - d
    if num % (2 * a):
        return None
    return num // (2 * a), y


def _reduction_steps(a, b, d, e, f, D, step):
    h0 = 2 * a * e - b * d
    K = h0 * h0 + D * (4 * a * f - d * d)
    step(f"X = 2ax + by + d = {2 * a}x + {b}y + {d},  Y = Δy − (2ae − bd) = {D}y − {h0}")
    step(f"Multiplying by 4aΔ and completing squares: Y² − {D}X² = {K}")
    return h0, K


def _hyper_square(a, b, c, d, e, f, D, step):
    o = _Out()
    o.kind = "hiperbolik, Δ tam kvadrat"
    s = isqrt(D)
    h0, K = _reduction_steps(a, b, d, e, f, D, step)
    step(f"Δ = {s}² ⇒ (Y − {s}X)(Y + {s}X) = {K}")
    if K != 0:
        for d1 in signed_divisors(K):
            d2 = K // d1
            if (d1 + d2) % 2 or (d2 - d1) % (2 * s):
                continue
            xy = _back(a, b, d, D, h0, (d2 - d1) // (2 * s), (d1 + d2) // 2)
            if xy:
                o.sols.add(xy)
        step(f"All divisor pairs of {K} checked; mapping back requires integrality")
        if not o.sols:
            o.proof = "No divisor pair gives an integer (x, y)."
        return o
    step("K = 0 ⇒ Y = ±sX: two linear equations")
    for sg in (1, -1):
        fam = linear2(-2 * a * sg * s, D - sg * s * b, h0 + sg * s * d)
        if fam:
            o.fams.append(fam)
    if not o.fams:
        o.proof = "Neither linear equation has an integer solution."
    return o


def _matmul(A, B, mod=None):
    R = [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
         [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]
    if mod:
        R = [[v % mod for v in row] for row in R]
    return R


def _matpow(A, k):
    R = [[1, 0], [0, 1]]
    while k:
        if k & 1:
            R = _matmul(R, A)
        A = _matmul(A, A)
        k >>= 1
    return R


def _hyper_pell(a, b, c, d, e, f, D, step):
    o = _Out()
    o.kind = "hiperbolik, Pell tipli"
    h0, K = _reduction_steps(a, b, d, e, f, D, step)
    if K == 0:
        step(f"Δ is not a square ⇒ Y² = {D}X² only for X = Y = 0")
        xy = _back(a, b, d, D, h0, 0, 0)
        if xy:
            o.sols.add(xy)
        else:
            o.proof = "The point X = Y = 0 does not correspond to an integer (x, y)."
        return o
    reps = gen_pell_reps(D, K)
    if not reps:
        o.proof = f"Y² − {D}X² = {K} has no solutions (LMM algorithm)."
        return o
    step(f"LMM gives {len(reps)} solution class(es) of Y² − {D}X² = {K}: "
         + ", ".join(f"(Y, X) = ({p}, {q})" for p, q in reps))
    p1, q1 = pell_fundamental(D)
    S = [[p1, D * q1], [q1, p1]]
    Si = [[p1, -D * q1], [-q1, p1]]
    Mod = abs(2 * a * D)
    P, cur = 1, [[v % Mod for v in row] for row in S]
    ident = [[1 % Mod, 0], [0, 1 % Mod]]
    while cur != ident:
        cur = _matmul(cur, S, Mod)
        P += 1
        if P > 2_000_000:
            o.partial = True
            o.notes.append("The order of the Pell matrix modulo M is too large.")
            return o
    step(f"Pell unit ε = {p1} + {q1}√{D}; its matrix has order {P} modulo {Mod} "
         f"⇒ the integrality pattern repeats with period {P}")
    SP, SPi = _matpow(S, P), _matpow(Si, P)

    def conj(SM):
        def T(x, y):
            Y, X = D * y - h0, 2 * a * x + b * y + d
            Y2, X2 = SM[0][0] * Y + SM[0][1] * X, SM[1][0] * Y + SM[1][1] * X
            yy = Fraction(Y2 + h0, D)
            xx = (X2 - b * yy - d) / Fraction(2 * a)
            return xx, yy
        c0 = T(0, 0)
        e1, e2 = T(1, 0), T(0, 1)
        L = [[e1[0] - c0[0], e2[0] - c0[0]], [e1[1] - c0[1], e2[1] - c0[1]]]
        for v in (L[0][0], L[0][1], L[1][0], L[1][1], c0[0], c0[1]):
            assert v.denominator == 1
        return [[int(v) for v in row] for row in L], [int(c0[0]), int(c0[1])]

    Tf, cf = conj(SP)
    Tb, cb = conj(SPi)
    bases = []
    for Y0, X0 in reps:
        for sg in (1, -1):
            Y, X = sg * Y0, sg * X0
            for _ in range(P):
                xy = _back(a, b, d, D, h0, X, Y)
                if xy and xy not in bases:
                    bases.append(xy)
                Y, X = S[0][0] * Y + S[0][1] * X, S[1][0] * Y + S[1][1] * X
    if not bases:
        o.proof = (f"Y² − {D}X² = {K} has solutions, but none maps back to "
                   f"an integer (x, y) (one full period was checked).")
        return o
    step(f"Within one period {len(bases)} integer starting point(s) were found; for each, "
         f"the step corresponding to ε^{P} is an integer affine map")
    for bs in bases:
        o.fams.append(RecursiveFamily(bs, Tf, cf, Tb, cb))
    return o


def solve(ctx):
    x, y = ctx.vars
    g = lambda *pr: ctx.ip.get(M(*pr), 0)
    a, b, c = g((x, 2)), g((x, 1), (y, 1)), g((y, 2))
    d, e, f = g((x, 1)), g((y, 1)), g()
    ctx.step(f"Coefficients: a={a}, b={b}, c={c}, d={d}, e={e}, f={f}  ({x} → x, {y} → y)")
    out = solve_bq(a, b, c, d, e, f, ctx.step)
    r = ctx.result(CODE, f"{TITLE} — {out.kind}")
    r.partial = out.partial
    r.notes += out.notes
    if out.proof and not out.sols and not out.fams:
        return r.none(out.proof)
    r.solutions = sorted(out.sols)
    r.families = out.fams
    return r
