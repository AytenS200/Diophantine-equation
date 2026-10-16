from ..numtheory import egcd, gcd, lll, size_reduce
from ..poly import Poly, M
from ..result import ParamFamily, Result

CODE, TITLE = "A1/A2", "Linear Diophantine equation"


def match(ctx):
    return ctx.poly.degree() == 1


def _lin_str(coeffs, names, c):
    p = Poly()
    for a, v in zip(coeffs, names):
        p = p + a * Poly.var(v)
    return f"{p} = {c}"


def hnf_solve(A, b):
    m, n = len(A), len(A[0])
    H = [list(r) for r in A]
    V = [[int(i == j) for j in range(n)] for i in range(n)]
    rows = []
    k = 0
    for i in range(m):
        if k < n:
            for j in range(k + 1, n):
                if H[i][j]:
                    p, q = H[i][k], H[i][j]
                    g, u, v = egcd(p, q)
                    pg, qg = p // g, q // g
                    for Mx in (H, V):
                        for row in Mx:
                            ck, cj = row[k], row[j]
                            row[k] = u * ck + v * cj
                            row[j] = -qg * ck + pg * cj
            if H[i][k]:
                rows.append((i, k))
                k += 1
                continue
        rows.append((i, None))
    y = [0] * n
    for i, col in rows:
        lim = k if col is None else col
        rem = b[i] - sum(H[i][j] * y[j] for j in range(lim))
        if col is None:
            if rem:
                return None, i
        else:
            if rem % H[i][col]:
                return None, i
            y[col] = rem // H[i][col]
    x0 = [sum(V[r][j] * y[j] for j in range(k)) for r in range(n)]
    basis = [[V[r][j] for r in range(n)] for j in range(k, n)]
    return (x0, basis), None


def _family_from(x0, basis):
    if basis:
        basis = lll(basis)
        x0 = size_reduce(x0, basis)
    params = ["t"] if len(basis) == 1 else [f"t{i + 1}" for i in range(len(basis))]
    exprs = []
    for r in range(len(x0)):
        e = Poly.const(x0[r])
        for p, vec in zip(params, basis):
            e = e + vec[r] * Poly.var(p)
        exprs.append(e)
    return ParamFamily(exprs, params)


def solve(ctx):
    names = ctx.vars
    a = [ctx.ip.get(M((v, 1)), 0) for v in names]
    c = -ctx.ip.get((), 0)
    if len(a) == 1:
        r = ctx.result("A1", "Linear equation in one variable")
        ctx.step(f"{_lin_str(a, names, c)}  ⇒  {names[0]} = {c}/{a[0]}")
        if c % a[0]:
            return r.none(f"{c} is not divisible by {a[0]}, so there is no integer solution.")
        r.solutions = [(c // a[0],)]
        return r
    if ctx.domain != "Z" and (all(v > 0 for v in a) or all(v < 0 for v in a)):
        return _bounded(ctx, a, c)
    if len(a) == 2:
        return _two(ctx, a, c)
    return _multi(ctx, a, c)


def _two(ctx, a, c):
    (x, y), (A, B) = ctx.vars, a
    r = ctx.result("A1", "Linear equation in two variables: ax + by = c")
    g, u, v = egcd(A, B)
    ctx.step(f"Equation: {_lin_str(a, ctx.vars, c)}")
    ctx.step(f"Extended Euclidean algorithm: gcd({A}, {B}) = {g}, "
             f"since {A}·({u}) + {B}·({v}) = {g}")
    if c % g:
        return r.none(f"The left side is always divisible by gcd({A}, {B}) = {g}, but {c} is not.")
    k = c // g
    x0, y0 = u * k, v * k
    ctx.step(f"Multiplying by {k} gives a particular solution: ({x}, {y}) = ({x0}, {y0})")
    dx, dy = B // g, -(A // g)
    s = -(x0 // abs(dx)) * (1 if dx > 0 else -1)
    x0, y0 = x0 + dx * s, y0 + dy * s
    ctx.step(f"General solution: {x} = {x0} + ({dx})t, {y} = {y0} + ({dy})t; "
             f"step (b/g, -a/g) = ({dx}, {dy})")
    t = Poly.var("t")
    r.families.append(ParamFamily([x0 + dx * t, y0 + dy * t], ["t"]))
    return r


def _multi(ctx, a, c):
    r = ctx.result("A2", f"Linear equation in {len(a)} variables")
    ctx.step(f"Equation: {_lin_str(a, ctx.vars, c)}")
    sol, _ = hnf_solve([a], [c])
    g = 0
    for v in a:
        g = gcd(g, v)
    ctx.step(f"gcd(coefficients) = {g} divides {c} ⇒ solutions exist")
    if sol is None:
        return r.none(f"gcd(coefficients) = {g} does not divide {c}.")
    x0, basis = sol
    ctx.step("Column operations build the Hermite normal form A·V = (g, 0, …, 0); "
             f"the last {len(basis)} columns of V span the homogeneous solutions; "
             "the basis is shortened with LLL")
    r.families.append(_family_from(x0, basis))
    if ctx.domain != "Z":
        r.notes.append("Coefficients have mixed signs, so the domain-restricted solution set may "
                       "be infinite; examples are filtered by the domain.")
    return r


def _bounded(ctx, a, c):
    r = ctx.result("A2", "Linear equation - non-negative/positive solutions (finite enumeration)")
    if a[0] < 0:
        a, c = [-v for v in a], -c
    off = 1 if ctx.domain == "P" else 0
    cc = c - off * sum(a)
    ctx.step("All coefficients share a sign ⇒ every variable is bounded, so there are finitely many solutions")
    if off:
        ctx.step(f"For positivity substitute v = v' + 1 for each variable: the right side becomes {cc}")
    if cc < 0:
        return r.none(f"Even the smallest possible value exceeds the right side ({c - cc} > {c}).")
    n = len(a)
    sg = [0] * (n + 1)
    for i in reversed(range(n)):
        sg[i] = gcd(a[i], sg[i + 1])
    sols, cap, state = [], 20000, {"trunc": False}

    def rec(i, rem, cur):
        if state["trunc"] or rem % sg[i]:
            return
        if i == n - 1:
            if rem % a[i] == 0:
                sols.append(tuple(v + off for v in cur + [rem // a[i]]))
                if len(sols) >= cap:
                    state["trunc"] = True
            return
        for v in range(rem // a[i] + 1):
            rec(i + 1, rem - a[i] * v, cur + [v])
            if state["trunc"]:
                return

    rec(0, cc, [])
    r.solutions = sols
    r.domain_applied = True
    if state["trunc"]:
        if cc <= 3_000_000:
            ways = [0] * (cc + 1)
            ways[0] = 1
            for ai in a:
                for s in range(ai, cc + 1):
                    ways[s] += ways[s - ai]
            r.total_count = ways[cc]
        r.notes.append(f"The list was truncated at {cap} solutions.")
    if not sols:
        return r.none("No non-negative combination reaches the right side "
                      "(all possibilities were checked).")
    return r


def solve_system(eqs, domain, limit):
    names = sorted({v for e in eqs for v in e.variables()})
    rows, rhs = [], []
    for e in eqs:
        p = e.integer_poly()
        if p.degree() > 1:
            r = Result("; ".join(x.text for x in eqs), names, domain, "—", "System of equations", limit)
            r.status = "error"
            r.notes.append("Only linear systems are supported.")
            return r
        rows.append([int(p.terms.get(M((v, 1)), 0)) for v in names])
        rhs.append(-int(p.terms.get((), 0)))
    r = Result("; ".join(x.text for x in eqs), names, domain, "A3", "System of linear equations", limit)
    r.trace.append(f"{len(rows)} equations, {len(names)} variables; written as A·x = b")
    r.trace.append("Column operations (extended Euclid at each step) produce the echelon form "
                   "A·V = H with V unimodular")
    sol, bad = hnf_solve(rows, rhs)
    if sol is None:
        r.none(f"Equation #{bad + 1} cannot be satisfied in integers together with the previous "
               "ones (a divisibility/consistency condition fails in echelon form).")
        r.finalize_status()
        return r
    x0, basis = sol
    if not basis:
        r.solutions = [tuple(x0)]
        r.trace.append("The homogeneous system has only the zero solution ⇒ the solution is unique")
    else:
        r.trace.append(f"Number of free parameters: {len(basis)}")
        r.families.append(_family_from(x0, basis))
    for e in eqs:
        for t in (r.solutions or r.families[0].samples(3, "Z")):
            assert e.check(dict(zip(names, t)))
    if domain != "Z":
        r.solutions = [t for t in r.solutions if all(v >= (1 if domain == "P" else 0) for v in t)]
    r.finalize_status()
    return r
