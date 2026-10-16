from fractions import Fraction

_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def M(*pairs):
    d = {}
    for v, e in pairs:
        d[v] = d.get(v, 0) + e
    return tuple(sorted((v, e) for v, e in d.items() if e))


def _mono_mul(m1, m2):
    if not m1:
        return m2
    if not m2:
        return m1
    d = dict(m1)
    for v, e in m2:
        d[v] = d.get(v, 0) + e
    return tuple(sorted(d.items()))


def mono_deg(m):
    return sum(e for _, e in m)


def fmt_num(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


class Poly:
    __slots__ = ("terms",)

    def __init__(self, terms=None):
        self.terms = {}
        if terms:
            for m, c in terms.items():
                if c != 0:
                    self.terms[m] = Fraction(c)

    @classmethod
    def const(cls, c):
        return cls({(): c})

    @classmethod
    def var(cls, name):
        return cls({((name, 1),): 1})

    @staticmethod
    def _co(o):
        return o if isinstance(o, Poly) else Poly.const(o)

    def __add__(self, o):
        o = Poly._co(o)
        t = dict(self.terms)
        for m, c in o.terms.items():
            t[m] = t.get(m, 0) + c
        return Poly(t)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.terms.items()})

    def __sub__(self, o):
        return self + (-Poly._co(o))

    def __rsub__(self, o):
        return Poly._co(o) - self

    def __mul__(self, o):
        o = Poly._co(o)
        t = {}
        for m1, c1 in self.terms.items():
            for m2, c2 in o.terms.items():
                m = _mono_mul(m1, m2)
                t[m] = t.get(m, 0) + c1 * c2
        return Poly(t)

    __rmul__ = __mul__

    def __pow__(self, k):
        if k < 0:
            raise ValueError("a polynomial cannot be raised to a negative power")
        r, b = Poly.const(1), self
        while k:
            if k & 1:
                r = r * b
            b = b * b
            k >>= 1
        return r

    def __eq__(self, o):
        return isinstance(o, Poly) and self.terms == o.terms

    __hash__ = None

    def is_zero(self):
        return not self.terms

    def is_const(self):
        return all(m == () for m in self.terms)

    def const_value(self):
        return self.terms.get((), Fraction(0))

    def variables(self):
        return sorted({v for m in self.terms for v, _ in m})

    def degree(self):
        return max((mono_deg(m) for m in self.terms), default=0)

    def degree_in(self, v):
        return max((e for m in self.terms for w, e in m if w == v), default=0)

    def univariate_coeffs(self, v):
        d = self.degree_in(v)
        cs = [Fraction(0)] * (d + 1)
        for m, c in self.terms.items():
            e = dict(m).get(v, 0)
            if len(m) > (1 if e else 0):
                raise ValueError("polynomial is not univariate")
            cs[e] += c
        return cs

    def evaluate(self, values):
        s = Fraction(0)
        for m, c in self.terms.items():
            t = c
            for v, e in m:
                t *= Fraction(values[v]) ** e
            s += t
        return s

    def substitute(self, mapping):
        res = Poly()
        for m, c in self.terms.items():
            t = Poly.const(c)
            for v, e in m:
                t = t * (mapping[v] if v in mapping else Poly.var(v)) ** e
            res = res + t
        return res

    def __str__(self):
        if not self.terms:
            return "0"
        single = all(len(v) == 1 for m in self.terms for v, _ in m)
        sep = "" if single else "·"
        items = sorted(self.terms.items(),
                       key=lambda mc: (-mono_deg(mc[0]), tuple((v, -e) for v, e in mc[0])))
        out = []
        for i, (m, c) in enumerate(items):
            neg = c < 0
            a = -c if neg else c
            if m == ():
                body = fmt_num(a)
            else:
                ms = sep.join(v + (str(e).translate(_SUP) if e != 1 else "") for v, e in m)
                if a == 1:
                    body = ms
                else:
                    body = fmt_num(a) + ("" if a.denominator == 1 and single else "·") + ms
            if i == 0:
                out.append(("-" if neg else "") + body)
            else:
                out.append((" - " if neg else " + ") + body)
        return "".join(out)

    __repr__ = __str__


class Rat:
    __slots__ = ("num", "den")

    def __init__(self, num, den=None):
        if den is None:
            den = Poly.const(1)
        if den.is_zero():
            raise ZeroDivisionError("division by zero")
        if den.is_const():
            num = num * (Fraction(1) / den.const_value())
            den = Poly.const(1)
        self.num, self.den = num, den

    def __add__(self, o):
        if self.den == o.den:
            return Rat(self.num + o.num, self.den)
        return Rat(self.num * o.den + o.num * self.den, self.den * o.den)

    def __neg__(self):
        return Rat(-self.num, self.den)

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        return Rat(self.num * o.num, self.den * o.den)

    def __truediv__(self, o):
        if o.num.is_zero():
            raise ZeroDivisionError("division by zero")
        return Rat(self.num * o.den, self.den * o.num)

    def __pow__(self, k):
        if k >= 0:
            return Rat(self.num ** k, self.den ** k)
        if self.num.is_zero():
            raise ZeroDivisionError("negative power of zero")
        return Rat(self.den ** (-k), self.num ** (-k))
