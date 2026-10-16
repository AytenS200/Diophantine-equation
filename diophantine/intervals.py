from fractions import Fraction

from .numtheory import isqrt


def _floor(q):
    q = Fraction(q)
    return q.numerator // q.denominator


def _val(cs, t):
    return sum(c * t ** i for i, c in enumerate(cs))


def poly_int_set(cs, strict):
    cs = [Fraction(c) for c in cs]
    while len(cs) > 1 and cs[-1] == 0:
        cs.pop()

    def ok(t):
        v = _val(cs, t)
        return v > 0 if strict else v >= 0

    deg = len(cs) - 1
    crit = []
    if deg == 1:
        crit.append(_floor(-cs[0] / cs[1]))
    elif deg == 2:
        disc = cs[1] ** 2 - 4 * cs[2] * cs[0]
        if disc >= 0:
            sq = Fraction(isqrt(disc.numerator * disc.denominator), disc.denominator)
            for sg in (1, -1):
                crit.append(_floor((-cs[1] + sg * sq) / (2 * cs[2])))
    elif deg > 2:
        raise ValueError("degree > 2")
    if not crit:
        return [(None, None)] if ok(0) else []
    pts = sorted({c + d for c in crit for d in range(-3, 4)})
    segs = [(None, pts[0] - 1, ok(pts[0] - 1))]
    for i, p in enumerate(pts):
        segs.append((p, p, ok(p)))
        if i + 1 < len(pts) and pts[i + 1] > p + 1:
            segs.append((p + 1, pts[i + 1] - 1, ok(p + 1)))
    segs.append((pts[-1] + 1, None, ok(pts[-1] + 1)))
    out = []
    for lo, hi, good in segs:
        if not good:
            continue
        if out and out[-1][1] is not None and lo is not None and out[-1][1] + 1 == lo:
            out[-1] = (out[-1][0], hi)
        else:
            out.append((lo, hi))
    return out


def _key(iv):
    return float("-inf") if iv[0] is None else iv[0]


def intersect(a, b):
    out = []
    for lo1, hi1 in a:
        for lo2, hi2 in b:
            lo = lo2 if lo1 is None else (lo1 if lo2 is None else max(lo1, lo2))
            hi = hi2 if hi1 is None else (hi1 if hi2 is None else min(hi1, hi2))
            if lo is None or hi is None or lo <= hi:
                out.append((lo, hi))
    return sorted(out, key=_key)


def intersect_all(sets):
    res = [(None, None)]
    for s in sets:
        res = intersect(res, s)
    return res


def is_bounded(ivs):
    return all(lo is not None and hi is not None for lo, hi in ivs)


def size(ivs):
    return sum(hi - lo + 1 for lo, hi in ivs)


def contains(ivs, t):
    return any((lo is None or t >= lo) and (hi is None or t <= hi) for lo, hi in ivs)


def iterate(ivs):
    for lo, hi in ivs:
        yield from range(lo, hi + 1)


def describe(ivs, name="t"):
    parts = []
    for lo, hi in ivs:
        if lo is None and hi is None:
            parts.append(f"{name} ∈ ℤ")
        elif lo is None:
            parts.append(f"{name} ≤ {hi}")
        elif hi is None:
            parts.append(f"{name} ≥ {lo}")
        elif lo == hi:
            parts.append(f"{name} = {lo}")
        else:
            parts.append(f"{lo} ≤ {name} ≤ {hi}")
    return " or ".join(parts) if parts else "empty"
