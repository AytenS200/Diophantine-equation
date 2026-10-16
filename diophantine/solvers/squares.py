from itertools import permutations

from ..numtheory import factorint, isqrt, sqrt_mod_prime, factor_str, divisors
from ..poly import M

CODE, TITLE = "A9/B5", "Sum of squares"


def prime_two_squares(p):
    if p == 2:
        return 1, 1
    r = sqrt_mod_prime(p - 1, p)
    a, b = p, r
    lim = isqrt(p)
    while b > lim:
        a, b = b, a % b
    y2 = p - b * b
    y = isqrt(y2)
    if y * y == y2:
        return b, y
    for x in range(isqrt(p) + 1):
        y2 = p - x * x
        if y2 >= 0 and isqrt(y2) ** 2 == y2:
            return x, isqrt(y2)
    raise ArithmeticError


def _gmul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def _gpow(z, e):
    r = (1, 0)
    for _ in range(e):
        r = _gmul(r, z)
    return r


def two_squares_obstruction(n):
    for p, e in factorint(n).items():
        if p % 4 == 3 and e % 2:
            return p, e
    return None


def two_squares_all(n):
    if n < 0:
        return []
    if n == 0:
        return [(0, 0)]
    base, choices = (1, 0), []
    for p, e in factorint(n).items():
        if p == 2:
            base = _gmul(base, _gpow((1, 1), e))
        elif p % 4 == 3:
            if e % 2:
                return []
            base = _gmul(base, (p ** (e // 2), 0))
        else:
            a, b = prime_two_squares(p)
            choices.append([_gmul(_gpow((a, b), j), _gpow((a, -b), e - j)) for j in range(e + 1)])
    combos = [base]
    for opts in choices:
        combos = [_gmul(z, w) for z in combos for w in opts]
    res = set()
    for z in combos:
        for u in ((1, 0), (0, 1), (-1, 0), (0, -1)):
            res.add(_gmul(z, u))
    return sorted(res)


def _canon(k, n, maxv):
    if k == 1:
        s = isqrt(n)
        if s * s == n and s <= maxv:
            yield (s,)
        return
    if k == 2:
        for x, y in two_squares_all(n):
            if 0 <= x <= y <= maxv:
                yield (x, y)
        return
    v = min(maxv, isqrt(n))
    while v >= 0 and k * v * v >= n:
        for rest in _canon(k - 1, n - v * v, v):
            yield rest + (v,)
        v -= 1


def _fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def count_variants(t, domain):
    k = len(t)
    perms = _fact(k)
    for v in set(t):
        perms //= _fact(t.count(v))
    nz = sum(1 for v in t if v)
    if domain == "Z":
        return perms * 2 ** nz
    if domain == "N":
        return perms
    return perms if nz == k else 0


def expand_variants(t, domain):
    out = set()
    for p in set(permutations(t)):
        if domain == "Z":
            signs = [(1, -1) if v else (1,) for v in p]
            combos = [()]
            for s in signs:
                combos = [c + (x,) for c in combos for x in s]
            for c in combos:
                out.add(tuple(a * b for a, b in zip(p, c)))
        elif domain == "N" or all(p):
            out.add(p)
    return out


def match(ctx):
    if len(ctx.vars) < 2 or ctx.poly.degree() != 2:
        return False
    sq = {M((v, 2)) for v in ctx.vars}
    if not set(ctx.ip) <= sq | {()} or not sq <= set(ctx.ip):
        return False
    cs = {ctx.ip[m] for m in sq}
    if cs not in ({1}, {-1}):
        return False
    s = cs.pop()
    ctx.data["sumsq_n"] = -ctx.ip.get((), 0) * s
    return True


def solve(ctx):
    k, n = len(ctx.vars), ctx.data["sumsq_n"]
    r = ctx.result("A9" if k == 2 else "B5", f"Sum of {k} squares = {n}")
    if n < 0:
        return r.none("A sum of squares cannot be negative.")
    if n > 0:
        ctx.step(f"Faktorizasiya: {n} = {factor_str(n)}")
    if k == 2:
        obs = two_squares_obstruction(n)
        if obs:
            p, e = obs
            return r.none(f"The prime {p} ≡ 3 (mod 4) divides {n} to an odd power ({e}); "
                          f"by the Fermat-Euler theorem such a number is not a sum of two squares.")
        for p in factorint(n):
            if p % 4 == 1:
                a, b = prime_two_squares(p)
                ctx.step(f"{p} ≡ 1 (mod 4) ⇒ Cornacchia: {p} = {a}² + {b}² = ({a}+{b}i)({a}−{b}i)")
        ctx.step("Choosing each Gaussian prime or its conjugate and multiplying by the "
                 "4 units (±1, ±i) produces every representation")
        sols = two_squares_all(n)
        if n:
            d1 = sum(1 for d in divisors(n) if d % 4 == 1)
            d3 = sum(1 for d in divisors(n) if d % 4 == 3)
            ctx.step(f"Yoxlama: r₂(n) = 4·(d₁ − d₃) = 4·({d1} − {d3}) = {4 * (d1 - d3)} = {len(sols)}")
        r.solutions = sols
        return r
    if k == 3:
        m, a = n, 0
        while m and m % 4 == 0:
            m //= 4
            a += 1
        if n and m % 8 == 7:
            return r.none(f"{n} = 4^{a}·(8·{m // 8} + 7); by the Legendre-Gauss three-square theorem "
                          f"such numbers are not sums of three squares.")
        ctx.step("Legendre condition holds: n ≠ 4ᵃ(8b + 7) ⇒ solutions exist")
    if k >= 4:
        ctx.step("Lagrange's theorem: every natural number is a sum of four squares ⇒ solutions exist")
    ctx.step("Canonical enumeration: 0 ≤ x₁ ≤ … ≤ xₖ, the largest element lies in √(n/k) ≤ xₖ ≤ √n, "
             "the rest recursively; the last step uses the two-square algorithm")
    reps, cap, trunc = [], 5000, False
    for t in _canon(k, n, isqrt(n)):
        if ctx.domain == "P" and 0 in t:
            continue
        reps.append(t)
        if len(reps) >= cap:
            trunc = True
            break
    total = sum(count_variants(t, ctx.domain) for t in reps)
    if k == 4 and ctx.domain == "Z" and n:
        r4 = 8 * sum(d for d in divisors(n) if d % 4)
        ctx.step(f"Check (Jacobi's formula): r₄({n}) = 8·Σ(d | n, 4 ∤ d) d = {r4}")
        if trunc:
            total = r4
    r.domain_applied = True
    if not reps:
        return r.none("No representation satisfies the domain constraint.")
    if total <= 300 and not trunc:
        allv = set()
        for t in reps:
            allv |= expand_variants(t, ctx.domain)
        r.solutions = sorted(allv)
    else:
        r.solutions = reps
        r.canonical = True
        r.total_count = total
        if trunc:
            r.notes.append(f"Canonical list truncated at {cap}.")
    return r
