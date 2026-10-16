import re
from itertools import product

from ..numtheory import gcd, modinv, crt_pair
from ..parser import parse_equation, parse_expr, normalize_text, ast_vars, ast_eval, ParseError
from ..poly import M
from ..result import Result, ParamFamily
from .linear import hnf_solve, _family_from

_PAT = re.compile(r"^(?P<l>.+?)(?:≡|==|=)(?P<r>.+?)\(?\s*mod\s*(?P<m>[^()]+?)\s*\)?\s*$", re.I)


def _parse(piece):
    m = _PAT.match(normalize_text(piece))
    if not m:
        raise ParseError(f"Congruence not recognised: '{piece}'. Format: 3x ≡ 4 (mod 7)")
    mod_ast = parse_expr(m.group("m"))
    if ast_vars(mod_ast):
        raise ParseError("The modulus must be a constant")
    mod = ast_eval(mod_ast, {})
    if mod.denominator != 1 or mod <= 0:
        raise ParseError("The modulus must be a positive integer")
    eq = parse_equation(f"{m.group('l')} = {m.group('r')}")
    p = eq.raw_poly()
    if any(c.denominator != 1 for c in p.terms.values()):
        raise ParseError("Fractional coefficients are not supported in congruences")
    return eq, {k: int(c) for k, c in p.terms.items()}, int(mod)


def _single_classes(ip, x, m, step):
    deg = max((e for mono in ip for _, e in mono), default=0)
    a = ip.get(M((x, 1)), 0)
    c0 = ip.get((), 0)
    if deg <= 1:
        b = -c0
        g = gcd(a, m)
        if a % m == 0:
            ok = b % m == 0
            step(f"{a}{x} ≡ {b} (mod {m}): the left side is always 0 ⇒ " + ("every x is a solution" if ok else "contradiction"))
            return [(0, 1)] if ok else []
        if b % g:
            step(f"{a}{x} ≡ {b} (mod {m}): gcd({a}, {m}) = {g} ∤ {b} ⇒ no solution")
            return []
        mm = m // g
        r = (b // g) * modinv(a // g, mm) % mm
        step(f"{a}{x} ≡ {b} (mod {m}): gcd = {g}; divide by {g} and multiply by the inverse ⇒ "
             f"{x} ≡ {r} (mod {mm})")
        return [(r, mm)]
    if m > 1_000_000:
        raise ParseError("For a non-linear congruence the modulus must be ≤ 10^6")
    cs = [ip.get(M((x, k)), 0) for k in range(deg + 1)]
    res = []
    for r0 in range(m):
        v = 0
        for c in reversed(cs):
            v = (v * r0 + c) % m
        if v == 0:
            res.append((r0, m))
    step(f"Degree {deg}: all {m} residues checked ⇒ {len(res)} solution(s)")
    return res


def solve(pieces, domain, limit):
    items = [_parse(p) for p in pieces]
    names = sorted({v for eq, _, _ in items for v in eq.variables()})
    text = "; ".join(pieces)
    r = Result(text, names, domain, "A4", "Congruence" if len(items) == 1
               else "System of congruences - Chinese remainder theorem", limit)
    step = r.trace.append
    if len(names) == 1:
        x = names[0]
        classes = [(0, 1)]
        for eq, ip, m in items:
            cl = _single_classes(ip, x, m, step)
            classes = [c for c in (crt_pair(r1, m1, r2, m2) for r1, m1 in classes for r2, m2 in cl) if c]
            if not classes:
                r.none(f"'{eq.text} (mod {m})' is incompatible with the previous conditions." if len(items) > 1
                       else "The congruence has no solution.")
                r.finalize_status()
                return r
        Mod = classes[0][1]
        rs = sorted({c[0] for c in classes})
        if len(items) > 1:
            step(f"Generalised CRT: write x = r₁ + m₁k, substitute into the next condition, "
                 f"combine pairwise ⇒ mod {Mod}")
        r.classes = (x, Mod, rs)
        r.finalize_status()
        return r
    if len(items) > 1:
        r.status = "error"
        r.notes.append("Systems of multivariate congruences are not supported.")
        return r
    eq, ip, m = items[0]
    deg = max((sum(e for _, e in mono) for mono in ip), default=0)
    if deg <= 1:
        a = [ip.get(M((v, 1)), 0) for v in names]
        c = -ip.get((), 0)
        step(f"The congruence becomes a linear Diophantine equation: Σaᵢxᵢ − {m}k = {c}")
        sol, _ = hnf_solve([a + [-m]], [c])
        if sol is None:
            r.none(f"gcd(coefficients, {m}) does not divide {c}.")
        else:
            x0, basis = sol
            fam = _family_from(x0[:-1], [b[:-1] for b in basis])
            r.families.append(fam)
        r.finalize_status()
        return r
    n = len(names)
    if m ** n > 2_000_000:
        r.status = "error"
        r.notes.append("The residue space is too large.")
        return r
    comp = [(c, [(names.index(v), e) for v, e in mono]) for mono, c in ip.items()]
    tups = []
    for vals in product(range(m), repeat=n):
        s = 0
        for c, mon in comp:
            t = c
            for i, e in mon:
                t *= vals[i] ** e
            s += t
        if s % m == 0:
            tups.append(vals)
    step(f"All {m ** n} residue combinations checked")
    if tups:
        r.residue_tuples = (m, tups)
    else:
        r.none(f"No residue combination modulo {m} works.")
    r.finalize_status()
    return r
