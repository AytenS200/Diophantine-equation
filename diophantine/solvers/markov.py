import heapq
from itertools import permutations

from ..poly import M
from ..result import GeneratorFamily, in_domain

CODE, TITLE = "B7", "Markov-type equation (Vieta jumping)"


def match(ctx):
    if len(ctx.vars) != 3:
        return False
    x, y, z = ctx.vars
    sq = [M((v, 2)) for v in ctx.vars]
    xyz = M((x, 1), (y, 1), (z, 1))
    if set(ctx.ip) != set(sq) | {xyz}:
        return False
    a = ctx.ip[sq[0]]
    if any(ctx.ip[m] != a for m in sq) or ctx.ip[xyz] % a:
        return False
    ctx.data["markov_k"] = -ctx.ip[xyz] // a
    return True


def markov_tree(k, root):
    heap, seen = [(max(root), root)], {root}
    while heap:
        _, t = heapq.heappop(heap)
        yield t
        for i in range(3):
            o = [t[j] for j in range(3) if j != i]
            nv = k * o[0] * o[1] - t[i]
            if nv > 0:
                nt = tuple(sorted(o + [nv]))
                if nt not in seen:
                    seen.add(nt)
                    heapq.heappush(heap, (max(nt), nt))


def solve(ctx):
    k = ctx.data["markov_k"]
    names = ctx.vars
    r = ctx.result(CODE, f"{TITLE}: {names[0]}² + {names[1]}² + {names[2]}² = {k}·{''.join(names)}")
    kk = abs(k)
    ctx.step(f"The equation is quadratic in {names[2]}: {names[2]}² − {k}{names[0]}{names[1]}·{names[2]} "
             f"+ ({names[0]}² + {names[1]}²) = 0 ⇒ the other root is {k}{names[0]}{names[1]} − {names[2]} (Vieta)")
    if kk not in (1, 3):
        ctx.step(f"Hurwitz's theorem (proved by Vieta descent): no positive solutions for k = {kk}")
        r.solutions = [(0, 0, 0)]
        r.notes.append("Only the trivial solution (0, 0, 0).")
        return r
    root = (1, 1, 1) if kk == 3 else (3, 3, 3)
    ctx.step(f"Vieta descent takes any positive solution down to the fundamental solution {root}; "
             f"conversely, Vieta moves from {root} produce every positive solution")
    if k < 0:
        ctx.step("k < 0: solutions are the variants with an odd number of negative coordinates")

    def gen(n, domain):
        out = [] if domain == "P" else [(0, 0, 0)]
        for cnt, t in enumerate(markov_tree(kk, root)):
            if len(out) >= n or cnt > 300:
                break
            variants = set()
            for p in set(permutations(t)):
                for sx in (1, -1):
                    for sy in (1, -1):
                        for sz in (1, -1):
                            neg = (sx < 0) + (sy < 0) + (sz < 0)
                            if (neg % 2 == 0) == (k > 0):
                                v = (sx * p[0], sy * p[1], sz * p[2])
                                if in_domain(v, domain):
                                    variants.add(v)
            out += sorted(variants)
        return out[:n]

    lines = [f"Root: {root};  Vieta move: ({names[0]}, {names[1]}, {names[2]}) → "
             f"({names[0]}, {names[1]}, {kk}{names[0]}{names[1]} − {names[2]}) and its symmetric variants",
             "Every positive solution lies in this tree (Markov tree); plus (0, 0, 0) and sign changes of "
             + ("an even" if k > 0 else "an odd") + " number of coordinates"]
    r.families.append(GeneratorFamily(lines, gen))
    r.domain_applied = True
    return r
