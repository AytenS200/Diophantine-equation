from fractions import Fraction
from itertools import product

from . import intervals as iv
from .poly import Poly

DOMAINS = {
    "Z": "all integers (ℤ)",
    "N": "non-negative integers (≥ 0)",
    "P": "positive integers (≥ 1)",
}


def in_domain(tup, domain):
    if domain == "N":
        return all(v >= 0 for v in tup)
    if domain == "P":
        return all(v > 0 for v in tup)
    return True


def fmt_tuple(t):
    return "(" + ", ".join(str(v) for v in t) + ")"


class ParamFamily:
    kind = "param"

    def __init__(self, exprs, params, conditions=None, ranges=None):
        self.exprs = [e if isinstance(e, Poly) else Poly.const(e) for e in exprs]
        self.params = list(params)
        self.conditions = list(conditions or [])
        self.ranges = dict(ranges or {})
        self.excluded = set()

    def value(self, assignment):
        out = []
        for e in self.exprs:
            v = e.evaluate(assignment)
            if v.denominator != 1:
                raise ArithmeticError("parametric expression produced a non-integer value")
            out.append(int(v))
        return tuple(out)

    def _single_values(self, p, cap):
        rng = self.ranges.get(p)
        if rng and iv.is_bounded(rng):
            yield from iv.iterate(rng)
            return
        start = 0
        if rng:
            for lo, hi in rng:
                if lo is not None and hi is None:
                    start = lo
                elif hi is not None and lo is None:
                    start = hi
        for k in range(cap):
            for t in ((start,) if k == 0 else (start + k, start - k)):
                if rng is None or iv.contains(rng, t):
                    yield t

    def samples(self, n, domain):
        out, seen = [], set()
        if len(self.params) == 1:
            p = self.params[0]
            for t in self._single_values(p, 20000):
                if t in self.excluded:
                    continue
                tup = self.value({p: t})
                if in_domain(tup, domain) and tup not in seen:
                    seen.add(tup)
                    out.append(tup)
                    if len(out) >= n:
                        break
            return out
        budget = 200000
        for R in range(0, 60):
            for vals in product(range(-R, R + 1), repeat=len(self.params)):
                if max((abs(v) for v in vals), default=0) != R:
                    continue
                budget -= 1
                tup = self.value(dict(zip(self.params, vals)))
                if in_domain(tup, domain) and tup not in seen:
                    seen.add(tup)
                    out.append(tup)
                    if len(out) >= n:
                        return out
            if budget <= 0:
                break
        return out

    def describe(self, names):
        lines = [f"{v} = {e}" for v, e in zip(names, self.exprs)]
        if len(self.params) == 1:
            p = self.params[0]
            lines.append(iv.describe(self.ranges[p], p) if p in self.ranges else f"{p} ∈ ℤ")
            if self.excluded:
                lines.append(f"{p} ≠ " + ", ".join(map(str, sorted(self.excluded)))
                             + "  (a denominator vanishes at these values)")
        elif self.params:
            lines.append(", ".join(self.params) + " ∈ ℤ")
        lines += self.conditions
        return lines


class RecursiveFamily:
    kind = "recursive"

    def __init__(self, base, Mf, cf, Mb, cb, note=None):
        self.base = tuple(base)
        self.Mf, self.cf, self.Mb, self.cb = Mf, cf, Mb, cb
        self.note = note

    @staticmethod
    def _apply(Mx, c, v):
        return tuple(sum(Mx[i][j] * v[j] for j in range(len(v))) + c[i] for i in range(len(v)))

    def terms(self, k):
        out = [self.base]
        v = self.base
        for _ in range(k):
            v = self._apply(self.Mf, self.cf, v)
            out.append(v)
        v = self.base
        for _ in range(k):
            v = self._apply(self.Mb, self.cb, v)
            out.append(v)
        return out

    def samples(self, n, domain):
        seen, res = set(), []
        for t in self.terms(max(25, n)):
            if in_domain(t, domain) and t not in seen:
                seen.add(t)
                res.append(t)
        res.sort(key=lambda t: (max(abs(v) for v in t), sum(v < 0 for v in t), t))
        return res[:n]

    def _lin(self, row, c, names):
        p = Poly.const(c)
        for coef, name in zip(row, names):
            p = p + coef * Poly.var(name)
        return str(p)

    def describe(self, names):
        lines = [f"Start:    ({', '.join(names)}) = {fmt_tuple(self.base)}"]
        lines.append("Forward:  " + ",  ".join(
            f"{names[i]}' = {self._lin(self.Mf[i], self.cf[i], names)}" for i in range(len(names))))
        lines.append("Backward: " + ",  ".join(
            f"{names[i]}' = {self._lin(self.Mb[i], self.cb[i], names)}" for i in range(len(names))))
        lines.append("(applying the step any number of times forward and backward yields every member)")
        if self.note:
            lines.append(self.note)
        return lines


class GeneratorFamily:
    kind = "generator"

    def __init__(self, lines, gen):
        self.lines, self.gen = list(lines), gen

    def samples(self, n, domain):
        return self.gen(n, domain)

    def describe(self, names):
        return self.lines


class Result:
    def __init__(self, equation, variables, domain, code, title, limit=10):
        self.equation = equation
        self.vars = list(variables)
        self.domain = domain
        self.code, self.title = code, title
        self.limit = limit
        self.status = None
        self.solutions, self.families = [], []
        self.proof = None
        self.notes, self.trace, self.extra_lines = [], [], []
        self.normal_form = None
        self.partial = False
        self.domain_applied = False
        self.canonical = False
        self.total_count = None
        self.classes = None
        self.residue_tuples = None

    def none(self, proof):
        self.proof = proof
        self.status = "none"
        return self

    def finalize_status(self):
        if self.status in ("identity", "info", "error"):
            return
        if self.classes or self.families:
            self.status = "infinite"
        elif self.solutions or self.residue_tuples:
            self.status = "finite"
        else:
            self.status = "none"

    def __str__(self):
        L = ["═" * 64]
        if self.equation:
            L.append(f"Equation:      {self.equation}")
        if self.normal_form:
            L.append(f"Normal form:   {self.normal_form} = 0")
        L.append(f"Type:          [{self.code}] {self.title}")
        if self.status != "info":
            L.append(f"Domain:        {DOMAINS[self.domain]}")
        if self.trace:
            L += ["", "Solution steps:"]
            L += [f"  {i}. {s}" for i, s in enumerate(self.trace, 1)]
        L.append("")
        vt = "(" + ", ".join(self.vars) + ")"
        if self.partial:
            L.append("⚠ PARTIAL RESULT - completeness is not guaranteed for this case")
        st = self.status
        if st == "identity":
            L.append("RESULT: IDENTITY - every integer assignment is a solution "
                     "(except those making a denominator zero)")
        elif st == "none":
            L.append("RESULT: NO SOLUTIONS")
            if self.proof:
                L.append(f"Reason: {self.proof}")
        elif st == "finite":
            if self.residue_tuples:
                Mod, tups = self.residue_tuples
                L.append(f"RESULT: {len(tups)} residue combinations modulo {Mod}")
                L += [f"  {vt} ≡ {fmt_tuple(t)} (mod {Mod})" for t in tups[:200]]
            else:
                cnt = self.total_count if self.total_count is not None else len(self.solutions)
                L.append(f"RESULT: finitely many solutions - {cnt} in total")
                if self.canonical:
                    L.append("  (canonical representatives shown; all other solutions are their "
                             "permutations and sign changes)")
                show = self.solutions[:max(self.limit * 10, 100)]
                L += [f"  {vt} = {fmt_tuple(t)}" for t in show]
                if len(self.solutions) > len(show):
                    L.append(f"  … and {len(self.solutions) - len(show)} more in the list")
        elif st == "infinite":
            L.append("RESULT: INFINITELY MANY SOLUTIONS")
            if self.classes:
                var, Mod, rs = self.classes
                L.append(f"  {var} ≡ " + ", ".join(map(str, rs)) + f"  (mod {Mod})")
                L.append(f"  i.e. {var} = r + {Mod}t, t ∈ ℤ")
                vals = sorted({r + Mod * k for r in rs for k in range(-3, 6)})
                vals = [v for v in vals if in_domain((v,), self.domain)]
                L.append("  Examples: " + ", ".join(map(str, vals[: self.limit])))
            if self.solutions:
                L.append("  Additional isolated solutions:")
                L += [f"    {vt} = {fmt_tuple(t)}" for t in self.solutions]
            for i, fam in enumerate(self.families, 1):
                L.append(f"  Family {i}:" if len(self.families) > 1 else "  General solution:")
                L += ["    " + line for line in fam.describe(self.vars)]
                smp = fam.samples(self.limit, self.domain)
                if smp:
                    L.append("    Examples: " + ", ".join(fmt_tuple(t) for t in smp))
                else:
                    L.append("    (no example in this domain among the first members of the family)")
        L += self.extra_lines
        L += [f"Note: {n}" for n in self.notes]
        return "\n".join(L)

    def __repr__(self):
        return (f"<Result {self.code} status={self.status!r} vars={self.vars} "
                f"solutions={len(self.solutions)} families={len(self.families)}>")

    def explain(self):
        return str(self)

    @property
    def is_solvable(self):
        if self.status in ("finite", "infinite", "identity"):
            return True
        if self.partial:
            return None
        return False if self.status == "none" else None

    @property
    def is_finite(self):
        return self.status in ("finite", "none")

    def iter_solutions(self, n=None, domain=None):
        domain = domain or self.domain
        count = 0
        for t in self.solutions:
            if n is not None and count >= n:
                return
            yield t
            count += 1
        if self.classes:
            var, Mod, rs = self.classes
            k = 0
            while n is None or count < n:
                for r in rs:
                    for v in ((r + Mod * k,) if k == 0 else (r + Mod * k, r - Mod * k)):
                        if in_domain((v,), domain) and (n is None or count < n):
                            yield (v,)
                            count += 1
                k += 1
                if n is None and k > 50:
                    return
        per = 50 if n is None else max(1, n)
        for fam in self.families:
            for t in fam.samples(per, domain):
                if n is not None and count >= n:
                    return
                yield t
                count += 1

    def to_dict(self, samples=10):
        fams = []
        for f in self.families:
            d = {"kind": f.kind, "description": f.describe(self.vars),
                 "examples": [list(t) for t in f.samples(samples, self.domain)]}
            if f.kind == "param":
                d["expressions"] = {v: str(e) for v, e in zip(self.vars, f.exprs)}
                d["parameters"] = list(f.params)
            elif f.kind == "recursive":
                d.update(base=list(f.base), forward=[f.Mf, f.cf], backward=[f.Mb, f.cb])
            fams.append(d)
        out = {
            "equation": self.equation, "variables": list(self.vars), "domain": self.domain,
            "type": self.code, "title": self.title, "status": self.status,
            "partial": self.partial, "solutions": [list(t) for t in self.solutions],
            "total_count": self.total_count, "canonical": self.canonical,
            "families": fams, "proof": self.proof, "steps": list(self.trace),
            "notes": list(self.notes),
        }
        if self.classes:
            var, Mod, rs = self.classes
            out["congruence_classes"] = {"variable": var, "modulus": Mod, "residues": rs}
        if self.residue_tuples:
            out["residue_tuples"] = {"modulus": self.residue_tuples[0],
                                     "tuples": [list(t) for t in self.residue_tuples[1]]}
        if self.extra_lines:
            out["info"] = list(self.extra_lines)
        return out
