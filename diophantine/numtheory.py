import random


def gcd(a, b):
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def gcd_list(xs):
    g = 0
    for x in xs:
        g = gcd(g, x)
        if g == 1:
            return 1
    return g


def lcm(a, b):
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


def egcd(a, b):
    old_r, r = a, b
    old_u, u = 1, 0
    old_v, v = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_u, u = u, old_u - q * u
        old_v, v = v, old_v - q * v
    if old_r < 0:
        old_r, old_u, old_v = -old_r, -old_u, -old_v
    return old_r, old_u, old_v


def modinv(a, m):
    g, u, _ = egcd(a, m)
    if g != 1:
        raise ValueError(f"{a} has no inverse modulo {m}")
    return u % m


def isqrt(n):
    if n < 0:
        raise ValueError("square root of a negative number")
    if n < 2:
        return n
    x = 1 << ((n.bit_length() + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            return x
        x = y


def is_square(n):
    if n < 0:
        return False
    r = isqrt(n)
    return r * r == n


_SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n):
    if n < 2:
        return False
    for p in _SMALL_PRIMES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _SMALL_PRIMES:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _pollard_brent(n):
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    rng = random.Random(n)
    while True:
        y, c, m = rng.randrange(1, n), rng.randrange(1, n), 128
        g = r = q = 1
        x = ys = y
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = gcd(q, n)
                k += m
            r <<= 1
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = gcd(abs(x - ys), n)
        if g != n:
            return g


def _factor_into(n, f):
    if n == 1:
        return
    if is_prime(n):
        f[n] = f.get(n, 0) + 1
        return
    d = _pollard_brent(n)
    _factor_into(d, f)
    _factor_into(n // d, f)


def factorint(n):
    n = abs(n)
    f = {}
    if n < 2:
        return f
    for p in (2, 3, 5):
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
    p, i = 7, 0
    wheel = (4, 2, 4, 2, 4, 6, 2, 6)
    while p * p <= n and p < 20000:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += wheel[i]
        i = (i + 1) % 8
    if n > 1:
        _factor_into(n, f)
    return dict(sorted(f.items()))


def divisors(n):
    n = abs(n)
    if n == 0:
        raise ValueError("0 has infinitely many divisors")
    ds = [1]
    for p, e in factorint(n).items():
        ds = [d * p ** k for d in ds for k in range(e + 1)]
    return sorted(ds)


def signed_divisors(n):
    ds = divisors(n)
    return [-d for d in reversed(ds)] + ds


def squarefree_decomp(n):
    if n == 0:
        return 0, 1
    sign = -1 if n < 0 else 1
    sf, s = 1, 1
    for p, e in factorint(n).items():
        s *= p ** (e // 2)
        if e % 2:
            sf *= p
    return sign * sf, s


def factor_str(n):
    f = factorint(n)
    if not f:
        return str(n)
    sup = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    return " · ".join(str(p) + (str(e).translate(sup) if e > 1 else "") for p, e in f.items())


def is_qr_prime(a, p):
    a %= p
    if a == 0 or p == 2:
        return True
    return pow(a, (p - 1) // 2, p) == 1


def sqrt_mod_prime(a, p):
    a %= p
    if a == 0:
        return 0
    if p == 2:
        return a
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def _sqrt_unit_pp(a, p, e):
    if e == 0:
        return [0]
    pe = p ** e
    if p == 2:
        if e == 1:
            return [1]
        if e == 2:
            return [1, 3] if a % 4 == 1 else []
        if a % 8 != 1:
            return []
        r = 1
        for k in range(3, e):
            if (r * r - a) % (1 << (k + 1)):
                r += 1 << (k - 1)
        half = 1 << (e - 1)
        return sorted({r % pe, (-r) % pe, (r + half) % pe, (-r + half) % pe})
    r = sqrt_mod_prime(a, p)
    if r is None:
        return []
    mod = p
    for _ in range(1, e):
        mod *= p
        r = (r - (r * r - a) * modinv(2 * r, mod)) % mod
    return sorted({r % pe, (-r) % pe})


def _sqrt_pp(a, p, e):
    pe = p ** e
    a %= pe
    if pe <= 1 << 16:
        return [z for z in range(pe) if z * z % pe == a]
    if a == 0:
        step = p ** ((e + 1) // 2)
        if pe // step > 10 ** 6:
            raise OverflowError("too many roots")
        return list(range(0, pe, step))
    k, a2 = 0, a
    while a2 % p == 0:
        a2 //= p
        k += 1
    if k % 2:
        return []
    base = _sqrt_unit_pp(a2, p, e - k)
    h = p ** (k // 2)
    mod_w = p ** (e - k)
    return sorted({h * (w + j * mod_w) % pe for w in base for j in range(h)})


def crt_pair(r1, m1, r2, m2):
    g, u, _ = egcd(m1, m2)
    if (r2 - r1) % g:
        return None
    l = m1 // g * m2
    k = ((r2 - r1) // g * u) % (m2 // g)
    return (r1 + m1 * k) % l, l


def sqrt_mod_all(a, n):
    if n == 1:
        return [0]
    a %= n
    if n <= 1 << 16:
        return [z for z in range(n) if z * z % n == a]
    res = [(0, 1)]
    for p, e in factorint(n).items():
        roots = _sqrt_pp(a, p, e)
        if not roots:
            return []
        pe = p ** e
        res = [crt_pair(r, m, s, pe) for (r, m) in res for s in roots]
    return sorted(r for r, _ in res)


def sqrt_cf(D):
    a0 = isqrt(D)
    if a0 * a0 == D:
        return a0, []
    m, d, a = 0, 1, a0
    period = []
    while True:
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d
        period.append(a)
        if a == 2 * a0:
            return a0, period


def _convergent(a0, per, k):
    p_prev, p = 1, a0
    q_prev, q = 0, 1
    L = len(per)
    for i in range(1, k + 1):
        a = per[(i - 1) % L]
        p_prev, p = p, a * p + p_prev
        q_prev, q = q, a * q + q_prev
    return p, q


def pell_fundamental(D):
    a0, per = sqrt_cf(D)
    if not per:
        raise ValueError("D is a perfect square")
    L = len(per)
    x, y = _convergent(a0, per, L - 1 if L % 2 == 0 else 2 * L - 1)
    assert x * x - D * y * y == 1
    return x, y


def pell_neg_fundamental(D):
    a0, per = sqrt_cf(D)
    if not per or len(per) % 2 == 0:
        return None
    x, y = _convergent(a0, per, len(per) - 1)
    assert x * x - D * y * y == -1
    return x, y


def lll(basis, delta=None):
    from fractions import Fraction
    delta = delta or Fraction(3, 4)
    b = [list(v) for v in basis]
    n = len(b)
    if n <= 1:
        return b

    def dot(u, v):
        return sum(x * y for x, y in zip(u, v))

    def gram_schmidt():
        bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = [Fraction(x) for x in b[i]]
            for j in range(i):
                mu[i][j] = dot(b[i], bs[j]) / dot(bs[j], bs[j])
                v = [x - mu[i][j] * y for x, y in zip(v, bs[j])]
            bs.append(v)
        return bs, mu

    bs, mu = gram_schmidt()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                b[k] = [x - q * y for x, y in zip(b[k], b[j])]
                bs, mu = gram_schmidt()
        if dot(bs[k], bs[k]) >= (delta - mu[k][k - 1] ** 2) * dot(bs[k - 1], bs[k - 1]):
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            bs, mu = gram_schmidt()
            k = max(k - 1, 1)
    return b


def size_reduce(x, basis, rounds=4):
    from fractions import Fraction
    x = list(x)
    for _ in range(rounds):
        changed = False
        for v in basis:
            vv = sum(c * c for c in v)
            if vv:
                q = round(Fraction(sum(a * c for a, c in zip(x, v)), vv))
                if q:
                    x = [a - q * c for a, c in zip(x, v)]
                    changed = True
        if not changed:
            break
    return x
