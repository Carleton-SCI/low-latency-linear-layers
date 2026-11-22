"""STM 130 nm areas with an exact XOR3 merge that may duplicate gates of fan-out two (Table 3, "ours STM").

With the STM unit areas (XOR2 3.33 GE, XOR3 4.66 GE), duplicating a non-output XOR2 gate of fan-out two
and absorbing one copy into each of the two gates that use it saves 2 x (2 x 3.33 - 4.66) - 3.33 =
0.67 GE. With the UMC unit areas (3.00, 4.67) the same change costs 0.34 GE, and duplicating a gate of
fan-out three or more never pays with either library; so this applies to the STM areas only.

A merge absorbs a non-output gate c into a gate p that uses it (p becomes an XOR3 gate). Each gate absorbs
at most one gate, and an absorbed gate absorbs none. The candidates are the edges (c, p) of gates c of
fan-out one (saving 2.00 GE) and the triples (c, p1, p2) of gates c of fan-out two (saving 0.67 GE); a
selection is valid when no gate is in two selected candidates. This script finds a selection of maximum
saving exactly, as a 0-1 program solved by CP-SAT (OR-tools), for every recorded circuit of a run, and
keeps the smallest area per matrix. It verifies each input circuit and each merged circuit (outputs and
depth); duplication does not change the XOR2 depth, and merging never adds gate levels.

usage: python3 stm_merge.py RUN_PARENT BENCH_DIR OUT_DIR
       RUN_PARENT contains run/{forward,bpbs,bbeam,lns}/ with the recorded circuits (*.xor2.txt,
       *.umc180.txt, *.stm130.txt); writes OUT_DIR/stm_merge.json and OUT_DIR/<matrix>.stm130.circ.txt
requires: pip install ortools"""
import csv
import glob
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from slpd.circuit import Circuit, read_any, circ_text          # noqa: E402
from slpd.matrix import read_matrices                         # noqa: E402

C2, C3 = 3.33, 4.66
W_EDGE = round(100 * (2 * C2 - C3))                 # 200: one merge
W_DUP = round(100 * (2 * (2 * C2 - C3) - C2))       # 67: duplication of a fan-out-two gate, two merges


def candidates(C):
    n, G = C.n, len(C.gates)
    users = [[] for _ in range(n + G)]
    for k, ops in enumerate(C.gates):
        for o in ops:
            users[o].append(k)
    is_out = set(C.out)
    edges, triples = [], []
    for k in range(G):
        if n + k in is_out:
            continue
        u = users[n + k]
        if len(u) == 1:
            edges.append((k, u[0]))
        elif len(u) == 2 and u[0] != u[1]:
            triples.append((k, u[0], u[1]))
    return edges, triples


def best_merge(C):
    """(selected edges, selected triples) of maximum saving, exactly (CP-SAT)"""
    from ortools.sat.python import cp_model
    edges, triples = candidates(C)
    if not edges and not triples:
        return [], []
    m = cp_model.CpModel()
    xe = [m.NewBoolVar(f"e{i}") for i in range(len(edges))]
    xt = [m.NewBoolVar(f"t{i}") for i in range(len(triples))]
    inc = {}
    for i, (c, p) in enumerate(edges):
        for v in (c, p):
            inc.setdefault(v, []).append(xe[i])
    for i, (c, p1, p2) in enumerate(triples):
        for v in (c, p1, p2):
            inc.setdefault(v, []).append(xt[i])
    for xs in inc.values():
        if len(xs) > 1:
            m.Add(sum(xs) <= 1)
    m.Maximize(W_EDGE * sum(xe) + W_DUP * sum(xt))
    s = cp_model.CpSolver()
    s.parameters.num_workers = 1
    s.parameters.max_time_in_seconds = 600
    if s.Solve(m) != cp_model.OPTIMAL:
        raise RuntimeError("CP-SAT did not prove optimality")
    return ([e for e, x in zip(edges, xe) if s.Value(x)], [t for t, x in zip(triples, xt) if s.Value(x)])


def apply_merge(C, sel_edges, sel_triples):
    n, G = C.n, len(C.gates)
    ops = [list(g) for g in C.gates]
    absorbed = set()
    for c, p in sel_edges:
        ops[p] = [o for o in ops[p] if o != n + c] + C.gates[c]
        absorbed.add(c)
    for c, p1, p2 in sel_triples:                 # one copy of c into each user
        for p in (p1, p2):
            ops[p] = [o for o in ops[p] if o != n + c] + C.gates[c]
        absorbed.add(c)
    new_id, gates = list(range(n)), []
    for k in range(G):
        if k in absorbed:
            new_id.append(None)
        else:
            new_id.append(n + len(gates))
            gates.append(ops[k])
    return Circuit(n, [[new_id[o] for o in g] for g in gates], [new_id[o] for o in C.out])


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    run, bench, out = sys.argv[1:]
    os.makedirs(out, exist_ok=True)
    rows = list(csv.DictReader(open(os.path.join(bench, "benchmark.tsv")), delimiter="\t"))
    summary = {}
    solve = []                                    # time of each exact merge (model and CP-SAT solve), in s
    for r in rows:
        name, D = r["file"][:-4], int(r["min_depth"])
        M = read_matrices(os.path.join(bench, r["file"]))[0]
        files = sorted(glob.glob(os.path.join(run, "run", "forward", name + "__*.txt"))
                       + glob.glob(os.path.join(run, "run", "*", name + ".seed*.txt")))
        files = [f for f in files if f.endswith((".xor2.txt", ".umc180.txt", ".stm130.txt"))]
        best = None
        for f in files:
            C = read_any(f, M.n).prune()
            if any(len(g) != 2 for g in C.gates) or C.check(M, D):
                raise SystemExit(f"{f}: not a correct XOR2 circuit of depth <= {D}")
            t0 = time.perf_counter()
            se, st = best_merge(C)
            solve.append(time.perf_counter() - t0)
            Mg = apply_merge(C, se, st)
            if Mg.check(M, D):
                raise SystemExit(f"{f}: the merged circuit is wrong: {Mg.check(M, D)}")
            area = round(Mg.area("stm130"), 2)
            if best is None or area < best[0] - 1e-9:
                best = (area, f, C, Mg, len(st))
        area, f, C, Mg, dups = best
        n2, n3 = Mg.counts()
        assert n2 + 2 * n3 == len(C.gates) + dups
        summary[name] = dict(stm=area, n2=n2, n3=n3, xor2_gates=len(C.gates), duplicated=dups,
                             xor2_depth=C.depth(), gate_levels=Mg.depth(), circuits=len(files),
                             source=os.path.relpath(f, run))
        extra = [f"XOR2 circuit before merging: {len(C.gates)} gates, depth {C.depth()}"
                 + (f"; {dups} gate(s) of fan-out two duplicated before merging" if dups else ""),
                 "Merge: exact for the STM unit areas, with duplication of gates of fan-out two (paper/tables/stm_merge.py)"]
        with open(os.path.join(out, name + ".stm130.circ.txt"), "w") as fh:
            fh.write(circ_text(Mg, M, name, source=os.path.relpath(f, run), extra=extra))
        print(f"{name:34s} {len(files):4d} circuits  STM {area:7.2f} GE  [{n2} XOR2 + {n3} XOR3]  "
              f"from {len(C.gates)} gates, {dups} duplicated", flush=True)
    with open(os.path.join(out, "stm_merge.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(f"exact merges: {len(solve)} circuits, optimal; time per circuit (one thread): "
          f"median {sorted(solve)[len(solve) // 2]:.4f} s, max {max(solve):.2f} s")


if __name__ == "__main__":
    main()
