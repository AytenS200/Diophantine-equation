import heapq

from ..numtheory import gcd_list
from ..result import Result


def frobenius_table(coins):
    a1 = coins[0]
    INF = float("inf")
    dist = [INF] * a1
    dist[0] = 0
    heap = [(0, 0)]
    while heap:
        d, r = heapq.heappop(heap)
        if d > dist[r]:
            continue
        for c in coins[1:]:
            nr, nd = (r + c) % a1, d + c
            if nd < dist[nr]:
                dist[nr] = nd
                heapq.heappush(heap, (nd, nr))
    return dist


def representable(n, coins, dist=None):
    dist = dist or frobenius_table(coins)
    return n >= 0 and n >= dist[n % coins[0]]


def solve(coins, limit=10):
    coins = sorted({abs(c) for c in coins if c})
    text = "frobenius(" + ", ".join(map(str, coins)) + ")"
    r = Result(text, [], "N", "A5", "Frobenius (coin) problem", limit)
    r.status = "info"
    step = r.trace.append
    if len(coins) < 2:
        r.extra_lines.append("At least two coins are needed.")
        return r
    g = gcd_list(coins)
    if g != 1:
        r.extra_lines.append(f"gcd = {g} > 1 ⇒ no number not divisible by {g} is representable; "
                             "the Frobenius number does not exist (it is infinite).")
        return r
    if coins[0] == 1:
        r.extra_lines.append("There is a coin of value 1 ⇒ every n ≥ 0 is representable; Frobenius number −1.")
        return r
    if coins[0] > 2_000_000:
        r.extra_lines.append("The smallest coin is too large.")
        return r
    if len(coins) == 2:
        a, b = coins
        step(f"Closed formula for two coins: F = ab − a − b = {a}·{b} − {a} − {b}")
    step(f"Residue graph: vertices 0..{coins[0] - 1} (mod {coins[0]}), edges r → (r + c) mod {coins[0]} "
         f"of weight c; Dijkstra finds the smallest representable number in every residue class")
    dist = frobenius_table(coins)
    F = max(dist) - coins[0]
    cnt = sum((dist[r0] - r0) // coins[0] for r0 in range(coins[0]))
    r.extra_lines.append(f"Frobenius number (largest non-representable value): {F}")
    r.extra_lines.append(f"Count of non-representable natural numbers: {cnt}")
    if cnt <= 150:
        bad = [n for n in range(1, F + 1) if n < dist[n % coins[0]]]
        r.extra_lines.append("They are: " + ", ".join(map(str, bad)))
    r.extra_lines.append(f"Every number greater than {F} is representable. For a specific n, solve the "
                         f"linear equation with domain N.")
    return r
