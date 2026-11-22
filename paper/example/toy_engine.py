"""Small-scale Python versions of the searches of src/dbs.cpp, used to compute the facts of the 8x8 example of
Sec. 3 of the paper (toy_data.py). They follow the rules of dbs.cpp (candidate gates,
pre-emptive gates, MD filter, tie rules, beam search, backward splits, neighborhood search), with Python's
random generator instead of std::mt19937_64, so single runs differ from dbs runs; the distributions agree.
The search program (src/dbs.cpp) does not use them.

Values are bitmasks over the inputs (bit c = input x_c). A target with depth bound b has the budget
B = 2^b, and a signal of depth d costs 2^d (Li et al., ToSC 2019(1), Lemma 1). Budget None means no depth
bound (the Boyar-Peralta heuristic of SEA 2010)."""
import itertools
import random

import numpy as np

INF = 10 ** 6


def wt(v):
    return bin(v).count("1")


def dmin(v):
    w, d = wt(v), 0
    while (1 << d) < w:
        d += 1
    return d


class Cand:
    __slots__ = ("S", "N", "v", "red", "d", "a", "b")

    def __init__(self, v, d, a, b):
        self.S, self.N, self.v, self.red, self.d, self.a, self.b = 0, 0, v, 0, d, a, b


class Solver:
    """Forward search state, as in dbs.cpp: base signals (inputs first), gates, targets with budgets TB
    (None = no depth bound), exact distances delta, found flags."""

    def __init__(self, n, T, TB, D, tie="all", md=False):
        self.n, self.T, self.TB, self.D = n, list(T), list(TB), D
        self.BUD = None if D is None else 1 << D
        self.tie, self.md = tie, md
        self.val = [1 << c for c in range(n)]
        self.dep = [0] * n
        self.gates = []                                   # (v, a, b, d), base index n + k
        self.where = {1 << c: c for c in range(n)}        # value -> base index of minimum depth
        self.delta = [0] * len(T)
        self.found = [False] * len(T)
        self.log = None                                   # list of pre-emptive gates added by cascade

    def copy(self):
        s = Solver.__new__(Solver)
        s.__dict__.update(self.__dict__)
        s.T, s.TB = self.T, self.TB
        s.val, s.dep, s.gates = list(self.val), list(self.dep), list(self.gates)
        s.where, s.delta, s.found = dict(self.where), list(self.delta), list(self.found)
        s.log = None
        return s

    # ------------------------------------------------------------------ distances (exact, by DP)
    def order(self):
        """base indices usable inside a sum (depth < D), stable-sorted by depth"""
        idx = [i for i in range(len(self.val)) if self.BUD is None or (1 << self.dep[i]) < self.BUD]
        return sorted(idx, key=lambda i: self.dep[i])

    def table(self):
        """R[w][x]: fewest signals of order() with XOR x and sum of 2^depth <= w (one layer if unbounded)"""
        size = 1 << self.n
        perm = np.arange(size)
        if self.BUD is None:
            R = np.full((1, size), INF, dtype=np.int64)
            R[0, 0] = 0
            for i in self.order():
                R[0] = np.minimum(R[0], R[0][perm ^ self.val[i]] + 1)
            return R
        R = np.full((self.BUD + 1, size), INF, dtype=np.int64)
        R[:, 0] = 0
        for i in self.order():
            c, v = 1 << self.dep[i], self.val[i]
            for w in range(self.BUD, c - 1, -1):
                R[w] = np.minimum(R[w], R[w - c][perm ^ v] + 1)
        return R

    def need(self, R, x, budget):
        """fewest signals of order() with XOR x within the budget (None = unbounded)"""
        if budget is None:
            return int(R[0][x])
        if budget < 0:
            return INF
        return int(R[min(budget, R.shape[0] - 1)][x])

    def recompute(self):
        R = self.table()
        for j, t in enumerate(self.T):
            if self.found[j]:
                continue
            k = self.need(R, t, self.TB[j])
            self.delta[j] = k - 1 if k < INF else INF
        return R

    # ------------------------------------------------------------------ search steps
    def init(self):
        for j, t in enumerate(self.T):
            w = wt(t)
            if self.TB[j] is not None and w > self.TB[j]:
                return False
            self.delta[j] = w - 1
            self.found[j] = w == 1
        self.log = []
        self.cascade()
        return True

    def complete(self):
        return all(self.found)

    def add(self, v, a, b):
        d = max(self.dep[a], self.dep[b]) + 1
        i = len(self.val)
        self.val.append(v)
        self.dep.append(d)
        self.gates.append((v, a, b, d))
        if v not in self.where or self.dep[self.where[v]] > d:
            self.where[v] = i
        for j, t in enumerate(self.T):
            if not self.found[j] and t == v and (self.TB[j] is None or (1 << d) <= self.TB[j]):
                self.found[j], self.delta[j] = True, 0
        return i

    def cascade(self):
        """add every target at distance 1, at the smallest possible depth (pre-emptive gates)"""
        while True:
            best = None
            order = self.order()
            for j, t in enumerate(self.T):
                if self.found[j] or self.delta[j] != 1:
                    continue
                for x in order:
                    y = self.where.get(t ^ self.val[x])
                    if y is None or y == x:
                        continue
                    d = max(self.dep[x], self.dep[y]) + 1
                    if (self.TB[j] is None or (1 << d) <= self.TB[j]) and (best is None or d < best[0]):
                        best = (d, j, x, y)
            if best is None:
                return
            d, j, x, y = best
            self.add(self.T[j], x, y)
            if self.log is not None:
                self.log.append(len(self.gates) - 1)
            self.recompute()

    def score(self):
        """all useful candidate gates, with S = sum of the new distances and N = sum of their squares"""
        order = self.order()
        cv = {}
        for p, q in itertools.combinations(range(len(order)), 2):
            a, b = order[p], order[q]
            d = max(self.dep[a], self.dep[b]) + 1
            if self.D is not None and d > self.D - 1:
                continue
            v = self.val[a] ^ self.val[b]
            if v in self.where:
                continue
            c = cv.get(v)
            if c is None:
                cv[v] = Cand(v, d, a, b)
            elif d < c.d:
                c.d, c.a, c.b = d, a, b
        R = self.table()
        open_j = [j for j in range(len(self.T)) if not self.found[j]]
        base = sum(self.delta[j] for j in open_j)
        base2 = sum(self.delta[j] ** 2 for j in open_j)
        mind = min((self.delta[j] for j in open_j), default=0)
        out = []
        for c in cv.values():
            for j in open_j:
                B = self.TB[j]
                rest = None if B is None else B - (1 << c.d)
                if self.need(R, self.T[j] ^ c.v, rest) <= self.delta[j] - 1:
                    c.red |= 1 << j
            if not c.red:
                continue
            if self.md and not any(c.red >> j & 1 and self.delta[j] == mind for j in open_j):
                continue
            c.S, c.N = base, base2
            for j in open_j:
                if c.red >> j & 1:
                    c.S -= 1
                    c.N -= 2 * self.delta[j] - 1
            out.append(c)
        return out

    def ties(self, cands):
        if not cands:
            return []
        S = min(c.S for c in cands)
        best = [c for c in cands if c.S == S]
        if self.tie == "bp":
            N = max(c.N for c in best)
            best = [c for c in best if c.N == N]
        elif self.tie == "sfx":
            N = min(c.N for c in best)
            best = [c for c in best if c.N == N]
        return best

    def apply(self, c):
        self.add(c.v, c.a, c.b)
        for j in range(len(self.T)):
            if c.red >> j & 1 and not self.found[j]:
                self.delta[j] -= 1
        self.log = []
        self.cascade()

    def gate_set(self):
        return frozenset((g[0], g[3]) for g in self.gates)


# ---------------------------------------------------------------------- netlists
def build_net(s, splits, outs):
    """forward state + split gates -> (gates [(a, b)], outputs) over node ids (inputs 0..n-1, gate k = n + k)"""
    n = s.n
    idx, dep, vals = {1 << c: c for c in range(n)}, [0] * n, [1 << c for c in range(n)]
    gates = []
    for v, a, b, d in s.gates:
        i = n + len(gates)
        gates.append((a, b))
        dep.append(d)
        vals.append(v)
        if v not in idx or dep[idx[v]] > d:
            idx[v] = i
    for w, p1, p2 in splits:
        if p1 not in idx or p2 not in idx:
            return None
        a, b = idx[p1], idx[p2]
        d = max(dep[a], dep[b]) + 1
        i = n + len(gates)
        gates.append((a, b))
        dep.append(d)
        vals.append(w)
        if w not in idx or dep[idx[w]] > d:
            idx[w] = i
    if any(t not in idx for t in outs):
        return None
    return gates, [idx[t] for t in outs]


def prune(n, gates, out):
    need = set(out)
    for k in range(len(gates) - 1, -1, -1):
        if n + k in need:
            need.update(gates[k])
    newid = {i: i for i in range(n)}
    g2 = []
    for k, (a, b) in enumerate(gates):
        if n + k in need:
            newid[n + k] = n + len(g2)
            g2.append((newid[a], newid[b]))
    return g2, [newid[o] for o in out]


def evaluate(n, gates, out):
    """values and depths of a circuit with gates of 2 or 3 operands"""
    val, dep = [1 << c for c in range(n)], [0] * n
    for ops in gates:
        v = 0
        for o in ops:
            v ^= val[o]
        val.append(v)
        dep.append(max(dep[o] for o in ops) + 1)
    return val, dep


def merge3(n, gates, out):
    """maximum merge of the non-output gates of fan-out one into their user (XOR3), as merge3 of dbs.cpp:
    returns matched[k] = the user gate that absorbs gate k, or -1"""
    ng = len(gates)
    fan, cons = [0] * (n + ng), [-1] * (n + ng)
    isout = set(out)
    for k, (a, b) in enumerate(gates):
        for x in (a, b):
            fan[x] += 1
            cons[x] = k
    ch = [[] for _ in range(ng)]
    has_parent = [False] * ng
    for k in range(ng):
        if (n + k) not in isout and fan[n + k] == 1:
            ch[cons[n + k]].append(k)
            has_parent[k] = True
    f0, f1, pick = [0] * ng, [-1] * ng, [-1] * ng
    for v in range(ng):
        base = sum(max(f0[c], f1[c]) for c in ch[v])
        f0[v] = base
        for c in ch[v]:
            val = base - max(f0[c], f1[c]) + f0[c] + 1
            if val > f1[v]:
                f1[v], pick[v] = val, c
    matched = [-1] * ng
    st = [(v, True) for v in range(ng - 1, -1, -1) if not has_parent[v]]
    while st:
        v, may = st.pop()
        use = may and f1[v] > f0[v]
        for c in ch[v]:
            if use and c == pick[v]:
                matched[c] = v
                st.append((c, False))
            else:
                st.append((c, True))
    return matched


def merged_circuit(n, gates, out):
    """the circuit after the maximum merge: gates with 2 or 3 operands (an XOR3 gate lists the absorbed pair last)"""
    m = merge3(n, gates, out)
    ops = [list(g) for g in gates]
    for c, user in enumerate(m):
        if user >= 0:
            ops[user] = [o for o in ops[user] if o != n + c] + list(gates[c])
    newid, g2 = {i: i for i in range(n)}, []
    for k in range(len(gates)):
        if m[k] < 0:
            newid[n + k] = n + len(g2)
            g2.append(ops[k])
    return [[newid[o] for o in g] for g in g2], [newid[o] for o in out], m


UNIT = {"umc180": (3.00, 4.67), "stm130": (3.33, 4.66)}


def area(n2, n3, lib="umc180"):
    a2, a3 = UNIT[lib]
    return round(a2 * n2 + a3 * n3, 2)


# ---------------------------------------------------------------------- forward runs
def run_forward(n, T, TB, D, tie, md, choose, record=None):
    """one forward run; choose(ties, step) -> index into the ties; record collects every iteration"""
    s = Solver(n, T, TB, D, tie, md)
    if not s.init():
        return None
    if record is not None:
        record.append(dict(kind="init", pre=list(s.log), delta=list(s.delta)))
    step = 0
    while not s.complete():
        cands = s.score()
        T_ = s.ties(cands)
        if not T_:
            return None
        c = T_[choose(T_, step)]
        if record is not None:
            Smin = min(x.S for x in cands)
            record.append(dict(kind="step", cands=len(cands), near=[(x.v, x.d, x.S, x.N) for x in cands if x.S == Smin],
                               ties=[(x.v, x.d) for x in T_], pick=(c.v, c.a, c.b, c.d), S=c.S, N=c.N))
        s.apply(c)
        if record is not None:
            record[-1].update(pre=list(s.log), delta=list(s.delta))
        step += 1
    return s


def outcome_distribution(n, T, TB, D, tie, md, splits=(), outs=None):
    """exact distribution of the pruned XOR2 count of a forward run whose ties are broken uniformly at random:
    {count: probability}; memoized on the set of gates of the state"""
    outs = list(T) if outs is None else outs
    memo = {}

    def rec(s):
        key = s.gate_set()
        if key in memo:
            return memo[key]
        if s.complete():
            net = build_net(s, splits, outs)
            g, o = prune(n, *net)
            res = {len(g): 1.0}
        else:
            T_ = s.ties(s.score())
            res = {}
            for c in T_:
                t = s.copy()
                t.apply(c)
                for k, p in rec(t).items():
                    res[k] = res.get(k, 0.0) + p / len(T_)
        memo[key] = res
        return res

    s = Solver(n, T, TB, D, tie, md)
    if not s.init():
        return {}
    return rec(s)


# ---------------------------------------------------------------------- beam search over the near-ties
def beam_search(n, T, TB, D, W, slack, tie, md, rng, record=None):
    """the beam search of dbs.cpp (Sec. 3.2 of the paper): returns the completed solvers"""
    B = [Solver(n, T, TB, D, tie, md)]
    if not B[0].init():
        return []
    while not all(b.complete() for b in B):
        opts, done = [], [b for b in B if b.complete()]
        for bi, b in enumerate(B):
            if b.complete():
                continue
            for c in b.score():
                opts.append((bi, c))
        if not opts:
            break
        bestS = min(c.S for _, c in opts)
        keep = [(bi, c) for bi, c in opts if c.S <= bestS + slack]
        rng.shuffle(keep)
        if tie == "bp":
            keep.sort(key=lambda o: (o[1].S, -o[1].N))
        elif tie == "sfx":
            keep.sort(key=lambda o: (o[1].S, o[1].N))
        else:
            keep.sort(key=lambda o: o[1].S)
        NB, seen, parents, dup = list(done), set(), set(), 0
        for bi, c in keep:
            if len(NB) >= W:
                break
            nb = B[bi].copy()
            nb.apply(c)
            h = frozenset(g[0] for g in nb.gates)
            if h in seen:
                dup += 1
                continue
            seen.add(h)
            NB.append(nb)
            parents.add(bi)
        if record is not None:
            record.append(dict(seqs=len(B), complete=len(done), options=len(opts), kept=len(keep), Smin=bestS,
                               parents_kept=len({bi for bi, _ in keep}), new=len(NB) - len(done), dup=dup,
                               parents_new=len(parents),
                               beam=[dict(gates=[(g[0], g[3]) for g in b.gates], delta=list(b.delta),
                                          done=b.complete()) for b in NB]))
        if not NB:
            break
        B = NB
    return [b for b in B if b.complete()]


# ---------------------------------------------------------------------- backward splitting (BPBS)
class Prefix:
    def __init__(self):
        self.Wv, self.P, self.splits = [], [], []      # unsplit outputs; (value, is predecessor); (w, p1, p2)

    def copy(self):
        p = Prefix()
        p.Wv, p.P, p.splits = list(self.Wv), list(self.P), list(self.splits)
        return p

    def targets(self, s):
        T, TB, pos = [], [], {}

        def put(v, b):
            if v in pos:
                TB[pos[v]] = min(TB[pos[v]], 1 << b)
            else:
                pos[v] = len(T)
                T.append(v)
                TB.append(1 << b)
        for w in self.Wv:
            put(w, s)
        for v, pred in self.P:
            put(v, s - 1 if pred else s)
        return T, TB

    def key(self):
        return (frozenset(self.Wv), frozenset(self.P))

    def inP(self, v):
        for i, (x, _) in enumerate(self.P):
            if x == v:
                return i
        return -1


def prefix_init(outs, s):                                     # Rule 1
    pf = Prefix()
    for v in outs:
        if dmin(v) < s:
            if pf.inP(v) < 0:
                pf.P.append((v, False))
        elif v not in pf.Wv:
            pf.Wv.append(v)
    return pf


def random_split(w, s, rng):
    bits = [i for i in range(64) if w >> i & 1]
    n, cap = len(bits), 1 << (s - 1)
    lo, hi = max(1, n - cap), min(cap, n - 1)
    if lo > hi:
        return 0
    a = rng.randint(lo, hi)
    rng.shuffle(bits)
    return sum(1 << b for b in bits[:a])


def draw_split(pf, s, q, rng):
    """the next split by the priority (3) of SFX23: Rule 2 > Rule 3 (with probability q), Rule 5"""
    if not pf.Wv:
        return None
    r2 = [(wi, pi) for wi in range(len(pf.Wv)) for pi in range(len(pf.P)) if pf.inP(pf.Wv[wi] ^ pf.P[pi][0]) > pi]
    if r2:
        wi, pi = rng.choice(r2)
        w, p = pf.Wv[wi], pf.P[pi][0]
        return (w, p, w ^ p), 2
    r3 = [(wi, pi) for wi in range(len(pf.Wv)) for pi in range(len(pf.P))
          if pf.Wv[wi] ^ pf.P[pi][0] and dmin(pf.Wv[wi] ^ pf.P[pi][0]) < s]
    if r3 and rng.random() < q:
        wi, pi = rng.choice(r3)
        w, p = pf.Wv[wi], pf.P[pi][0]
        return (w, p, w ^ p), 3
    w = rng.choice(pf.Wv)
    p1 = random_split(w, s, rng)
    if not p1:
        return None
    return (w, p1, w ^ p1), 5


def apply_split(pf, sp):
    w, p1, p2 = sp
    pf.Wv.remove(w)
    for p in (p1, p2):
        i = pf.inP(p)
        if i < 0:
            pf.P.append((p, True))
        else:
            pf.P[i] = (p, True)
    pf.splits.append(sp)


def evaluate_prefix(n, pf, s, outs, R, tie, md, rng):
    """R randomized forward runs on the targets of the prefix; returns the best pruned circuit (gates, out)"""
    T, TB = pf.targets(s)
    best = None
    for _ in range(R):
        st = run_forward(n, T, TB, s, tie, md, lambda T_, k: rng.randrange(len(T_)))
        if st is None:
            continue
        net = build_net(st, pf.splits, outs)
        if net is None:
            continue
        g, o = prune(n, *net)
        if best is None or len(g) < len(best[0]):
            best = (g, o)
    return best


def merged_umc(n, net):
    """(XOR2 count, XOR2 after merging, XOR3 after merging, UMC area) of a pruned circuit"""
    g, o = net
    m = merge3(n, g, o)
    n3 = sum(x >= 0 for x in m)
    return len(g), len(g) - 2 * n3, n3, area(len(g) - 2 * n3, n3)


def key_xor(n, net):
    """the objective of dbs.cpp with --obj xor: the XOR2 count, then the UMC area"""
    e = merged_umc(n, net)
    return (e[0], e[3])


def bpbs(n, outs, s, N, maxsplits, R, q, tie, md, rng, record=None, pool=None):
    """stage 2: N random trajectories of splits; every prefix is evaluated by R forward runs"""
    best = None
    for tr in range(N):
        pf = prefix_init(outs, s)
        steps = []
        while pf.Wv and len(pf.splits) < maxsplits:
            d = draw_split(pf, s, q, rng)
            if d is None:
                break
            sp, rule = d
            apply_split(pf, sp)
            net = evaluate_prefix(n, pf, s, outs, R, tie, md, rng)
            steps.append(dict(rule=rule, split=sp, result=None if net is None else len(net[0])))
            if net is not None:
                if pool is not None:
                    pool.append(net)
                if best is None or key_xor(n, net) < key_xor(n, best):
                    best = net
        if record is not None:
            record.append(steps)
    return best


def split_beam(n, outs, s, reps, maxsplits, Bw, K, R, q, tie, md, rng, record=None, pool=None):
    """stage 3: beam search over the splits (Bw prefixes kept, K split options drawn per prefix)"""
    best = None
    for rep in range(reps):
        cur = [prefix_init(outs, s)]
        for step in range(maxsplits):
            nxt, seen = [], set()
            for pf in cur:
                for _ in range(K):
                    d = draw_split(pf, s, q, rng)
                    if d is None:
                        continue
                    ch = pf.copy()
                    apply_split(ch, d[0])
                    if ch.key() in seen:
                        continue
                    seen.add(ch.key())
                    net = evaluate_prefix(n, ch, s, outs, R, tie, md, rng)
                    if net is None:
                        continue
                    if pool is not None:
                        pool.append(net)
                    nxt.append((key_xor(n, net), ch, net))
            if not nxt:
                break
            nxt.sort(key=lambda x: x[0])
            nxt = nxt[:Bw]
            if record is not None:
                record.append(dict(rep=rep, step=step + 1, prefixes=len(nxt), best=nxt[0][0][0],
                                   splits=[len(ch.splits) for _, ch, _ in nxt]))
            if best is None or nxt[0][0] < key_xor(n, best):
                best = nxt[0][2]
            cur = [ch for _, ch, _ in nxt]
    return best


# ---------------------------------------------------------------------- large neighborhood search
def lns(n, outs, D, gates, out, iters, kmin, kmax, W, slack, tie, rng, record=None):
    """stage 4: keep the last gate of k random outputs, re-solve all the other gates by beam search"""
    inc = (list(gates), list(out))
    for it in range(iters):
        val, dep = evaluate(n, inc[0], inc[1])
        js, seen = [], set()
        for j, t in enumerate(outs):
            if inc[1][j] >= n and t not in seen:
                seen.add(t)
                js.append(j)
        rng.shuffle(js)
        k = min(rng.randint(kmin, max(kmin, kmax)), len(js))
        pf = Prefix()
        fixed = set()
        for j in js[:k]:
            a, b = inc[0][inc[1][j] - n]
            pf.splits.append((outs[j], val[a], val[b]))
            fixed.add(outs[j])
            for p in (val[a], val[b]):
                i = pf.inP(p)
                if i < 0:
                    pf.P.append((p, True))
                else:
                    pf.P[i] = (p, True)
        for t in outs:
            if t in fixed or t in pf.Wv:
                continue
            if dmin(t) < D:
                if pf.inP(t) < 0:
                    pf.P.append((t, False))
            else:
                pf.Wv.append(t)
        T, TB = pf.targets(D)
        done = beam_search(n, T, TB, D, W, slack, tie, False, rng)
        best = None
        for st in done:
            net = build_net(st, pf.splits, outs)
            if net is None:
                continue
            g, o = prune(n, *net)
            if best is None or len(g) < len(best[0]):
                best = (g, o)
        accepted = best is not None and len(best[0]) <= len(inc[0])
        if record is not None:
            record.append(dict(it=it, k=k, fixed=[outs[j] for j in js[:k]], splits=list(pf.splits),
                               result=None if best is None else len(best[0]), before=len(inc[0]), accepted=accepted,
                               net=best if accepted and len(best[0]) < len(inc[0]) else None))
        if accepted:
            inc = best
    return inc
