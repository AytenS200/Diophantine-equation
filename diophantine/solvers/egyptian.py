from fractions import Fraction

from ..parser import ast_eval, ast_vars

CODE, TITLE = "A8", "Fraction equation (Egyptian fractions)"


class _NoMatch(Exception):
    pass


def _const_times_var(node):
    if node[0] == "var":
        return Fraction(1), node[1]
    if node[0] == "mul":
        l, r = node[1], node[2]
        if not ast_vars(l) and r[0] == "var":
            return ast_eval(l, {}), r[1]
        if not ast_vars(r) and l[0] == "var":
            return ast_eval(r, {}), l[1]
    raise _NoMatch


def _collect(node, sign, acc):
    k = node[0]
    if k in ("add", "sub"):
        _collect(node[1], sign, acc)
        _collect(node[2], sign if k == "add" else -sign, acc)
        return
    if k == "neg":
        _collect(node[1], -sign, acc)
        return
    if not ast_vars(node):
        acc["const"] += sign * ast_eval(node, {})
        return
    if k == "div" and not ast_vars(node[1]):
        coef, var = _const_times_var(node[2])
        acc["vars"][var] = acc["vars"].get(var, 0) + sign * ast_eval(node[1], {}) / coef
        return
    if k == "pow" and node[2] == -1 and node[1][0] == "var":
        v = node[1][1]
        acc["vars"][v] = acc["vars"].get(v, 0) + sign
        return
    raise _NoMatch


def match(ctx):
    acc = {"const": Fraction(0), "vars": {}}
    try:
        _collect(ctx.eq.lhs, 1, acc)
        _collect(ctx.eq.rhs, -1, acc)
    except (_NoMatch, ZeroDivisionError):
        return False
    vs = acc["vars"]
    if len(vs) < 3 or any(c == 0 for c in vs.values()) or sorted(vs) != ctx.vars:
        return False
    ctx.data["egypt"] = (vs, -acc["const"])
    return True


def _floor(q):
    return q.numerator // q.denominator


def _ceil(q):
    return -((-q.numerator) // q.denominator)


def solve_positive(items, Q, cap=200000):
    results = set()

    def rec(remaining, q, maxterm, assign):
        if len(results) >= cap or q <= 0:
            return
        if len(remaining) == 1:
            var, c = remaining[0]
            v = c / q
            if v.denominator == 1 and v >= 1:
                assign[var] = int(v)
                results.add(tuple(sorted(assign.items())))
                del assign[var]
            return
        m = len(remaining)
        for idx, (var, c) in enumerate(remaining):
            lo = _floor(c / q) + 1
            if maxterm is not None:
                lo = max(lo, _ceil(c / maxterm))
            hi = _floor(m * c / q)
            rest = remaining[:idx] + remaining[idx + 1:]
            for v in range(lo, hi + 1):
                assign[var] = v
                term = c / v
                rec(rest, q - term, term, assign)
                del assign[var]

    rec(sorted(items.items()), Q, None, {})
    return results


def solve(ctx):
    vs, Q = ctx.data["egypt"]
    r = ctx.result(CODE, f"{TITLE}, {len(vs)} variables")
    terms = " + ".join(f"{c}/{v}" for v, c in sorted(vs.items()))
    ctx.step(f"Equation: {terms} = {Q}")
    if all(c < 0 for c in vs.values()):
        vs, Q = {v: -c for v, c in vs.items()}, -Q
    if not all(c > 0 for c in vs.values()) or Q <= 0:
        return None
    k = len(vs)
    ctx.step(f"In a positive solution the largest term c/x is ≥ Q/{k} ⇒ that variable "
             f"satisfies c/Q < x ≤ {k}·c/Q; recurse on the remaining terms")
    found = solve_positive(vs, Q)
    names = ctx.vars
    r.solutions = sorted(tuple(dict(t)[v] for v in names) for t in found)
    r.domain_applied = True
    if ctx.domain == "Z":
        r.partial = True
        r.solutions = r.solutions
        r.notes.append("The complete set of positive solutions was found. If negative values are "
                       f"allowed the equation usually has infinitely many solutions (e.g. pairs c/x + (−c/x) = 0), "
                       "so over ℤ the result is partial. Use domain P for a complete answer.")
    if not r.solutions:
        return r.none("All possible positive values were checked - there is no solution.")
    return r
