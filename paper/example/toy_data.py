"""Compute and check the facts of the 8x8 example of Sec. 3 of the paper; paper/tables/worked_example.py reads them.
usage: python3 paper/example/toy_data.py > paper/example/toy_data.json     (needs numpy and ortools; about 3 min)
A re-run gives the same file, except that CP-SAT (8 workers) may return another 16-gate optimum in "opt_circuit",
which paper/tables/worked_example.py does not use.

The example is the 8x8 circulant matrix with y_j = x_j + x_{j+1} + x_{j+2} + x_{j+3} + x_{j+4} (indices mod 8),
at its minimum depth D = 3. Every circuit below is checked against the matrix and the depth bound (slpd), and
its areas are computed as in the paper: UMC by the maximum merge of the gates of fan-out one, STM by the exact
merge with duplication of the gates of fan-out two (paper/tables/stm_merge.py)."""
import collections
import itertools
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "paper", "tables"))
from toy_engine import (Solver, run_forward, beam_search, build_net, prune, evaluate, merge3, merged_umc,  # noqa: E402
                        prefix_init, draw_split, apply_split, evaluate_prefix, split_beam, lns, key_xor, area, wt)
from slpd.circuit import Circuit                     # noqa: E402
from slpd.matrix import Matrix                       # noqa: E402
import stm_merge                                     # noqa: E402

N, D = 8, 3
ROWS = [sum(1 << ((j + i) % N) for i in range(5)) for j in range(N)]
M = Matrix(ROWS, N, "window5")
RUNS = 2000                    # forward runs per configuration (random restarts)
POOL = {}                      # every distinct circuit found by the toy searches: key -> (gates, out)


def check(gates, out, depth_bound=D):
    C = Circuit(N, [list(g) for g in gates], out)
    problems = C.check(M, depth_bound)
    assert not problems, problems
    return C


def add_pool(net):
    g, o = net
    val, _ = evaluate(N, g, o)
    key = frozenset((val[N + k], frozenset((val[a], val[b]))) for k, (a, b) in enumerate(g))
    POOL.setdefault(key, (list(g), list(o)))


def stm_exact(gates, out):
    C = check(gates, out)
    se, st = stm_merge.best_merge(C)
    K = stm_merge.apply_merge(C, se, st)
    assert not K.check(M, D)
    n2, n3 = K.counts()
    return dict(n2=n2, n3=n3, stm=round(K.area("stm130"), 2), dup=[t[0] for t in st], merges=[e[0] for e in se])


def circ_json(net, depth_bound=D):
    """a pruned XOR2 circuit with its values, depths, outputs, merges and areas"""
    g, o = net
    C = check(g, o, depth_bound)
    val, dep = evaluate(N, g, o)
    m = merge3(N, g, o)
    n3 = sum(x >= 0 for x in m)
    res = dict(gates=[dict(a=a, b=b, val=val[N + k], depth=dep[N + k], out=[j for j, x in enumerate(o) if x == N + k],
                           absorbed_by=m[k]) for k, (a, b) in enumerate(g)],
               out=list(o), count=len(g), depth=max(dep[x] for x in o), n2=len(g) - 2 * n3, n3=n3,
               umc=area(len(g) - 2 * n3, n3), stm_plain=area(len(g) - 2 * n3, n3, "stm130"))
    if depth_bound == D:
        res["stm_exact"] = stm_exact(g, o)
    return res


# ------------------------------------------------------------------------------------------------ optimum
def optimum():
    """exact minimum at depth 3 (CP-SAT): depth-1 gates are pairs of inputs, depth-2 gates have weight 3 (pair +
    input) or 4 (two disjoint pairs), and the 8 outputs are distinct gates of depth 3 (weight 5)"""
    from ortools.sat.python import cp_model
    inputs = [1 << c for c in range(N)]
    pairs = [a | b for a, b in itertools.combinations(inputs, 2)]
    m = cp_model.CpModel()
    P = {p: m.NewBoolVar("") for p in pairs}
    Q, dec = {}, {}
    for S in itertools.combinations(range(N), 3):
        v = sum(1 << c for c in S)
        Q[v] = m.NewBoolVar("")
        dec[v] = [(v & ~(1 << c), 1 << c) for c in S]
    for S in itertools.combinations(range(N), 4):
        v = sum(1 << c for c in S)
        Q[v] = m.NewBoolVar("")
        dec[v] = [((1 << S[0]) | (1 << b), v & ~((1 << S[0]) | (1 << b))) for b in S[1:]]
    avail = {v: None for v in inputs}
    avail.update(P)
    avail.update(Q)
    qsel = {}
    for v, ds in dec.items():
        opts = []
        for A, B in ds:
            x = m.NewBoolVar("")
            for p in (A, B):
                if avail[p] is not None:
                    m.AddImplication(x, avail[p])
            opts.append((x, A, B))
        m.AddBoolOr([x for x, _, _ in opts]).OnlyEnforceIf(Q[v])
        qsel[v] = opts
    osel = []
    for r in ROWS:
        opts = []
        for A in avail:
            B = r ^ A
            if B in avail and A < B:
                x = m.NewBoolVar("")
                for p in (A, B):
                    if avail[p] is not None:
                        m.AddImplication(x, avail[p])
                opts.append((x, A, B))
        m.AddBoolOr([x for x, _, _ in opts])
        osel.append(opts)
    m.Minimize(sum(P.values()) + sum(Q.values()))
    s = cp_model.CpSolver()
    s.parameters.num_workers = 8
    s.parameters.max_time_in_seconds = 600
    assert s.Solve(m) == cp_model.OPTIMAL
    k = int(round(s.ObjectiveValue()))
    # build the optimal circuit from the solution
    node, gates = {1 << c: c for c in range(N)}, []

    def gate(A, B):
        gates.append((node[A], node[B]))
        node[A ^ B] = N + len(gates) - 1
    for p in pairs:
        if s.Value(P[p]):
            gate(p & -p, p & (p - 1))
    for v in sorted(Q, key=wt):
        if s.Value(Q[v]):
            x, A, B = next(t for t in qsel[v] if s.Value(t[0]))
            gate(A, B)
    out = []
    for j, r in enumerate(ROWS):
        x, A, B = next(t for t in osel[j] if s.Value(t[0]))
        gate(A, B)
        out.append(N + len(gates) - 1)
    g, o = prune(N, gates, out)
    assert len(g) == len(ROWS) + k
    return len(g), (g, o)


# ------------------------------------------------------------------------------------------------ forward runs
def traced_run(T, TB, Dd, tie, md, rng):
    """a forward run that records, per iteration, the candidates with and without the MD filter, the near-ties,
    the ties after the norm, the chosen gate and the pre-emptive gates"""
    s = Solver(N, T, TB, Dd, tie, md)
    s.init()
    steps = [dict(kind="init", pre=[dict(val=s.gates[i][0], depth=s.gates[i][3]) for i in s.log], delta=list(s.delta))]
    while not s.complete():
        allc = Solver.score(_with(s, md=False))
        cands = s.score()
        ties = s.ties(cands)
        Smin = min(c.S for c in cands)
        c = ties[rng.randrange(len(ties))]
        st = dict(kind="step", cands=len(allc), filtered=len(cands), near=sum(x.S == Smin for x in cands),
                  ties=len(ties), S=c.S, N=c.N, gate=dict(val=c.v, a=s.val[c.a], b=s.val[c.b], depth=c.d),
                  lowered=[j for j in range(len(T)) if c.red >> j & 1], before=list(s.delta))
        s.apply(c)
        st.update(pre=[dict(val=s.gates[i][0], depth=s.gates[i][3]) for i in s.log], delta=list(s.delta))
        steps.append(st)
    return s, steps


def _with(s, md):
    t = s.copy()
    t.md = md
    return t


def distribution(TB, Dd, tie, md, seed, runs=RUNS, depth=False):
    rng = random.Random(seed)
    h = collections.Counter()
    for _ in range(runs):
        st = run_forward(N, ROWS, TB, Dd, tie, md, lambda T_, k: rng.randrange(len(T_)))
        net = prune(N, *build_net(st, [], ROWS))
        if Dd == D:
            add_pool(net)
        if depth:
            val, dep = evaluate(N, *net)
            h[(len(net[0]), max(dep[x] for x in net[1]))] += 1
        else:
            h[len(net[0])] += 1
    return h


def main():
    out = dict(n=N, D=D, rows=ROWS, naive=sum(wt(r) - 1 for r in ROWS))
    assert M.rank() == N and M.min_depth() == D

    # the exact optimum at depth 3
    opt, opt_net = optimum()
    out["optimum"] = opt
    out["opt_circuit"] = circ_json(opt_net)

    # Boyar-Peralta without a depth bound: distribution over random tie choices, and one traced run
    h = distribution([None] * N, None, "bp", False, seed=11, depth=True)
    out["bp_dist"] = sorted([k[0], k[1], v] for k, v in h.items())
    for seed in itertools.count(1):
        rng = random.Random(1000 + seed)
        s, steps = traced_run(ROWS, [None] * N, None, "bp", False, rng)
        net = prune(N, *build_net(s, [], ROWS))
        val, dep = evaluate(N, *net)
        if len(net[0]) == 15 and max(dep[x] for x in net[1]) == 7:
            break
    out["bp_run"] = dict(seed=1000 + seed, steps=steps, circuit=circ_json(net, depth_bound=None))

    # depth-bounded forward runs (D = 3): the tie rules, with and without the MD filter
    TB = [1 << D] * N
    dists = {}
    for name, tie, md in (("sfx_md", "sfx", True), ("all", "all", False), ("bp", "bp", False), ("sfx", "sfx", False),
                          ("all_md", "all", True), ("bp_md", "bp", True)):
        dists[name] = dict(sorted(distribution(TB, D, tie, md, seed=21).items()))
    out["fwd_dist"] = dists
    assert dists["sfx_md"] == {18: RUNS}
    s, steps = traced_run(ROWS, TB, D, "sfx", True, random.Random(31))
    sfx18 = prune(N, *build_net(s, [], ROWS))
    out["sfx_run"] = dict(steps=steps, circuit=circ_json(sfx18))
    for seed in itertools.count(1):                   # a run with random ties and no norm that ends at 16
        s, steps = traced_run(ROWS, TB, D, "all", False, random.Random(40 + seed))
        net = prune(N, *build_net(s, [], ROWS))
        if len(net[0]) == 16:
            break
    out["all_run"] = dict(seed=40 + seed, steps=steps, circuit=circ_json(net))

    # stage 2 / Shi-Feng-Xu: backward splitting (random trajectories, each prefix evaluated by 4 forward runs)
    rng = random.Random(51)
    trajs, first16, hits, best = [], None, 0, None
    for tr in range(400):
        pf = prefix_init(ROWS, D)
        steps = []
        while pf.Wv and len(pf.splits) < 8:
            dsp = draw_split(pf, D, 0.5, rng)
            if dsp is None:
                break
            sp, rule = dsp
            apply_split(pf, sp)
            net = evaluate_prefix(N, pf, D, ROWS, 4, "sfx", True, rng)
            steps.append(dict(rule=rule, split=list(sp), result=None if net is None else len(net[0])))
            if net is not None:
                add_pool(net)
                if best is None or key_xor(N, net) < key_xor(N, best):
                    best = net
                if len(net[0]) == 16 and first16 is None:
                    first16 = dict(traj=tr, steps=list(steps), targets=[list(x) for x in zip(*pf.targets(D))],
                                   circuit=circ_json(net))
                    net16 = net
        trajs.append(min((x["result"] for x in steps if x["result"] is not None), default=None))
        hits += any(x["result"] == 16 for x in steps)
    out["bpbs"] = dict(trajectories=len(trajs), reach16=hits, best=key_xor(N, best)[0],
                       best_dist=dict(sorted(collections.Counter(trajs).items())), first16=first16)
    assert first16 is not None

    # HILL: maximum merge (Rule 1) of the XOR2 circuits of Shi-Feng-Xu (the 16-gate and the 18-gate circuit)
    out["hill"] = dict(of16=merged_umc(N, net16)[1:], of18=merged_umc(N, sfx18)[1:])

    # stage 1: beam search over the near-ties
    for seed in itertools.count(1):
        rec = []
        done = beam_search(N, ROWS, TB, D, 4, 0, "all", False, random.Random(60 + seed), rec)
        nets = [prune(N, *build_net(b, [], ROWS)) for b in done]
        if min(len(x[0]) for x in nets) == 16:
            break
    for x in nets:
        add_pool(x)
    out["beam"] = dict(seed=60 + seed, W=4, s=0, iters=rec, results=sorted(len(x[0]) for x in nets),
                       circuit=circ_json(min(nets, key=lambda x: key_xor(N, x))))
    stats = collections.Counter()
    for sd in range(1, 21):
        done = beam_search(N, ROWS, TB, D, 4, 0, "all", False, random.Random(sd))
        stats[min(len(prune(N, *build_net(b, [], ROWS))[0]) for b in done)] += 1
    out["beam_seeds"] = dict(sorted(stats.items()))
    rec = []
    done = beam_search(N, ROWS, TB, D, 4, 0, "sfx", True, random.Random(61), rec)
    out["beam_sfxmd"] = dict(iters=rec, results=sorted(len(prune(N, *build_net(b, [], ROWS))[0]) for b in done))
    res = collections.Counter()
    for sd in range(1, 21):                       # the rule of Shi-Feng-Xu with and without their filter
        for md in (True, False):
            done = beam_search(N, ROWS, TB, D, 4, 0, "sfx", md, random.Random(sd))
            res[("md" if md else "nomd", min(len(prune(N, *build_net(b, [], ROWS))[0]) for b in done))] += 1
    out["beam_sfx_seeds"] = sorted([k[0], k[1], v] for k, v in res.items())
    rec = []
    done = beam_search(N, ROWS, TB, D, 4, 1, "all", False, random.Random(62), rec)
    nets = [prune(N, *build_net(b, [], ROWS)) for b in done]
    for x in nets:
        add_pool(x)
    out["beam_s1"] = dict(iters=[dict((k, v) for k, v in r.items() if k != "beam") for r in rec],
                          results=sorted(len(x[0]) for x in nets),
                          best_umc=circ_json(min(nets, key=lambda x: merged_umc(N, x)[3])))

    # stage 3: beam search over the splits
    pool3 = []
    rec = []
    b3 = split_beam(N, ROWS, D, 1, 8, 4, 8, 2, 0.5, "sfx", True, random.Random(71), rec, pool3)
    for x in pool3:
        add_pool(x)
    out["split_beam"] = dict(steps=rec, best=key_xor(N, b3)[0])

    # stage 4: neighborhood search from the 18-gate circuit of the Shi-Feng-Xu forward search
    for seed in itertools.count(1):
        rec = []
        inc = lns(N, ROWS, D, sfx18[0], sfx18[1], 30, 1, 4, 4, 0, "all", random.Random(80 + seed), rec)
        if len(inc[0]) == 16:
            break
    add_pool(inc)
    out["lns"] = dict(seed=80 + seed, iters=rec, circuit=circ_json(inc), start=circ_json(sfx18))
    # the circuit after each improvement
    imp = []
    for r in rec:
        if r["net"] is not None:
            imp.append(dict(it=r["it"], k=r["k"], fixed=r["fixed"], splits=[list(x) for x in r["splits"]],
                            before=r["before"], after=r["result"], circuit=circ_json(r["net"])))
        del r["net"]
    out["lns"]["improvements"] = imp

    # selection by area over every distinct circuit found by the toy searches
    rows = []
    for key, net in POOL.items():
        c, n2, n3, umc = merged_umc(N, net)
        rows.append((c, umc, n2, n3, net))
    by_count = collections.defaultdict(list)
    for r in rows:
        by_count[r[0]].append(r)
    table = []
    for c in sorted(by_count):
        rs = by_count[c]
        b = min(rs, key=lambda r: r[1])
        table.append(dict(count=c, circuits=len(rs), best_umc=b[1], n2=b[2], n3=b[3]))
    out["pool"] = dict(circuits=len(rows), by_count=table)
    best_umc = min(rows, key=lambda r: (r[1], r[0]))
    out["best_umc"] = circ_json(best_umc[4])
    stm_rows = []
    for r in rows:                                    # the exact STM merge of every circuit
        e = stm_exact(*r[4])
        stm_rows.append((e["stm"], r[0], r[4], e))
    best_stm = min(stm_rows, key=lambda x: (x[0], x[1]))
    out["best_stm"] = circ_json(best_stm[2])
    out["opt_stm"] = out["opt_circuit"]["stm_exact"]
    json.dump(out, sys.stdout, indent=1)


if __name__ == "__main__":
    main()
