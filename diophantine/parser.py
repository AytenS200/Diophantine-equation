import re
from fractions import Fraction

from .poly import Poly, Rat


class ParseError(ValueError):
    pass


_SUPMAP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
           "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9"}


def normalize_text(s):
    for a, b in (("−", "-"), ("–", "-"), ("—", "-"), ("·", "*"), ("×", "*"),
                 ("⋅", "*"), ("÷", "/"), ("**", "^")):
        s = s.replace(a, b)
    out, i = [], 0
    while i < len(s):
        if s[i] in _SUPMAP:
            digits = ""
            while i < len(s) and s[i] in _SUPMAP:
                digits += _SUPMAP[s[i]]
                i += 1
            out.append("^" + digits)
        else:
            out.append(s[i])
            i += 1
    return "".join(out)


def tokenize(s):
    toks, i, n = [], 0, len(s)
    while i < n:
        ch = s[i]
        if ch.isspace():
            i += 1
        elif ch.isdigit() or (ch == "." and i + 1 < n and s[i + 1].isdigit()):
            j = i
            while j < n and (s[j].isdigit() or s[j] == "."):
                j += 1
            try:
                toks.append(("NUM", Fraction(s[i:j])))
            except ValueError:
                raise ParseError(f"Invalid number: '{s[i:j]}'")
            i = j
        elif ch.isalpha():
            j = i + 1
            while j < n and (s[j].isdigit() or s[j] == "_"):
                j += 1
            toks.append(("VAR", s[i:j]))
            i = j
        elif ch in "+-*/^()":
            toks.append(("OP", ch))
            i += 1
        else:
            raise ParseError(f"Unknown character: '{ch}'")
    out = []
    for t in toks:
        if out:
            p = out[-1]
            left = p[0] in ("NUM", "VAR") or p == ("OP", ")")
            right = t[0] == "VAR" or t == ("OP", "(") or (t[0] == "NUM" and p[0] != "NUM")
            if left and right:
                out.append(("OP", "*"))
        out.append(t)
    return out


class _Parser:
    def __init__(self, toks):
        self.t, self.i = toks, 0

    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else (None, None)

    def take(self):
        tok = self.peek()
        self.i += 1
        return tok

    def parse(self):
        if not self.t:
            raise ParseError("Empty expression")
        node = self.expr()
        if self.i != len(self.t):
            raise ParseError(f"Unexpected token: '{self.peek()[1]}'")
        return node

    def expr(self):
        node = self.term()
        while self.peek() in (("OP", "+"), ("OP", "-")):
            op = self.take()[1]
            node = ("add" if op == "+" else "sub", node, self.term())
        return node

    def term(self):
        node = self.unary()
        while self.peek() in (("OP", "*"), ("OP", "/")):
            op = self.take()[1]
            node = ("mul" if op == "*" else "div", node, self.unary())
        return node

    def unary(self):
        if self.peek() in (("OP", "+"), ("OP", "-")):
            op = self.take()[1]
            a = self.unary()
            return a if op == "+" else ("neg", a)
        return self.power()

    def power(self):
        base = self.atom()
        if self.peek() == ("OP", "^"):
            self.take()
            exp_node = self.unary()
            if ast_vars(exp_node):
                raise ParseError("A variable appears in an exponent; exponential equations "
                                 "are outside the supported (level 1 and 2) families")
            k = ast_eval(exp_node, {})
            if k.denominator != 1:
                raise ParseError("Exponent must be an integer")
            return ("pow", base, int(k))
        return base

    def atom(self):
        kind, val = self.take()
        if kind == "NUM":
            return ("num", val)
        if kind == "VAR":
            return ("var", val)
        if (kind, val) == ("OP", "("):
            node = self.expr()
            if self.take() != ("OP", ")"):
                raise ParseError("Unclosed parenthesis")
            return node
        raise ParseError("Incomplete expression" if kind is None else f"Unexpected token: '{val}'")


def parse_expr(text):
    return _Parser(tokenize(normalize_text(text))).parse()


def ast_vars(node):
    k = node[0]
    if k == "num":
        return set()
    if k == "var":
        return {node[1]}
    if k in ("neg", "pow"):
        return ast_vars(node[1])
    return ast_vars(node[1]) | ast_vars(node[2])


def ast_eval(node, values):
    k = node[0]
    if k == "num":
        return node[1]
    if k == "var":
        return Fraction(values[node[1]])
    if k == "neg":
        return -ast_eval(node[1], values)
    if k == "pow":
        b, e = ast_eval(node[1], values), node[2]
        if e < 0:
            if b == 0:
                raise ZeroDivisionError
            return 1 / b ** (-e)
        return b ** e
    a, b = ast_eval(node[1], values), ast_eval(node[2], values)
    if k == "add":
        return a + b
    if k == "sub":
        return a - b
    if k == "mul":
        return a * b
    if b == 0:
        raise ZeroDivisionError
    return a / b


def ast_to_rat(node):
    k = node[0]
    if k == "num":
        return Rat(Poly.const(node[1]))
    if k == "var":
        return Rat(Poly.var(node[1]))
    if k == "neg":
        return -ast_to_rat(node[1])
    if k == "pow":
        return ast_to_rat(node[1]) ** node[2]
    a, b = ast_to_rat(node[1]), ast_to_rat(node[2])
    return {"add": a.__add__, "sub": a.__sub__, "mul": a.__mul__, "div": a.__truediv__}[k](b)


def integerize(p):
    from .numtheory import lcm, gcd_list
    if p.is_zero():
        return p
    L = 1
    for c in p.terms.values():
        L = lcm(L, c.denominator)
    ints = {m: int(c * L) for m, c in p.terms.items()}
    g = gcd_list(ints.values())
    return Poly({m: c // g for m, c in ints.items()})


class Equation:
    def __init__(self, text, lhs, rhs):
        self.text, self.lhs, self.rhs = text, lhs, rhs

    def variables(self):
        return sorted(ast_vars(self.lhs) | ast_vars(self.rhs))

    def check(self, values):
        try:
            return ast_eval(self.lhs, values) == ast_eval(self.rhs, values)
        except ZeroDivisionError:
            return None

    def denominators(self):
        out = []

        def walk(node):
            k = node[0]
            if k in ("num", "var"):
                return
            if k == "div" and ast_vars(node[2]):
                out.append(ast_to_rat(node[2]).num)
            if k == "pow" and node[2] < 0 and ast_vars(node[1]):
                out.append(ast_to_rat(node[1]).num)
            for ch in node[1:]:
                if isinstance(ch, tuple):
                    walk(ch)

        walk(self.lhs)
        walk(self.rhs)
        return out

    def raw_poly(self):
        try:
            r = ast_to_rat(self.lhs) - ast_to_rat(self.rhs)
        except ZeroDivisionError:
            raise ParseError("The expression divides by zero")
        return r.num

    def integer_poly(self):
        return integerize(self.raw_poly())


def parse_equation(text):
    t = normalize_text(text)
    parts = re.split(r"==|=|≡", t)
    if len(parts) == 1:
        parts.append("0")
    if len(parts) != 2:
        raise ParseError("An equation may contain only one '=' sign")
    lhs, rhs = (_Parser(tokenize(p)).parse() for p in parts)
    return Equation(text.strip(), lhs, rhs)
