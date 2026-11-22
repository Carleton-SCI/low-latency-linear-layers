"""Starting circuits for stage 4 from the published circuits of Shi, Feng and Xu (SFX23).

For every benchmark matrix of minimum depth 3, takes the SFX23 circuits whose matrix equals it up to row and
column permutations (as an isomorphism of bipartite graphs) and whose depth is at most 3, maps the one with
the fewest XOR2 gates (then the smallest UMC area after a maximum merge) onto the benchmark matrix by the
permutation, verifies it, and writes it in the native format of src/dbs.cpp (read by its --lns option).

usage: python3 paper/compare/sfx_to_dbs.py SFX_RESULTS_DIR BENCH_DIR OUT_DIR
       SFX_RESULTS_DIR: the folder of result files of the SFX23 repository (our_results/)
       writes OUT_DIR/<matrix>.sfx.txt for paper/extra_runs.sh stage4fromsfx
requires: networkx"""
import csv
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
from slpd.circuit import Circuit, read_sfx_text, read_text, dbs_text   # noqa: E402
from slpd.matrix import read_matrices                                 # noqa: E402
import networkx as nx                                                  # noqa: E402
from networkx.algorithms import isomorphism                           # noqa: E402


def rows_of(C):
    val, _ = C.evaluate()
    return [val[o] for o in C.out]


def graph(rows, n):
    G = nx.Graph()
    for j in range(len(rows)):
        G.add_node(("r", j), side=0)
    for c in range(n):
        G.add_node(("c", c), side=1)
    for j, r in enumerate(rows):
        for c in range(n):
            if (r >> c) & 1:
                G.add_edge(("r", j), ("c", c))
    return G


def mapped(C, M):
    """C relabeled to compute M exactly (row and column permutations), or None"""
    GM, GC = graph(M.rows, M.n), graph(rows_of(C), C.n)
    gm = isomorphism.GraphMatcher(GC, GM, node_match=lambda a, b: a["side"] == b["side"])
    if not gm.is_isomorphic():
        return None
    mp = gm.mapping                                          # node of GC -> node of GM
    col = [mp[("c", c)][1] for c in range(C.n)]              # input c of C is input col[c] of M
    row = [mp[("r", j)][1] for j in range(len(C.out))]       # output j of C is output row[j] of M
    gates = [[col[o] if o < C.n else o for o in g] for g in C.gates]
    out = [None] * len(C.out)
    for j, o in enumerate(C.out):
        out[row[j]] = col[o] if o < C.n else o
    P = Circuit(C.n, gates, out)
    return P if not P.check(M, 3) else None


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sfx_dir, bench, out_dir = sys.argv[1:]
    os.makedirs(out_dir, exist_ok=True)
    circuits = []
    for f in sorted(glob.glob(os.path.join(sfx_dir, "*.txt"))):
        for k, C in enumerate(read_sfx_text(read_text(f), f)):
            C = C.prune()
            if C.depth() <= 3:
                circuits.append((os.path.basename(f), k, C))
    for r in csv.DictReader(open(os.path.join(bench, "benchmark.tsv")), delimiter="\t"):
        name = r["file"][:-4]
        if r["min_depth"] != "3":
            continue
        M = read_matrices(os.path.join(bench, r["file"]))[0]
        best = None
        for fname, k, C in sorted(circuits, key=lambda x: (len(x[2].gates), x[2].merge_xor3().area("umc180"))):
            if C.n != M.n or len(C.out) != M.m or len(C.gates) > (best[1] if best else 10 ** 9):
                continue
            P = mapped(C, M)
            if P is not None:
                best = (f"{fname}#{k}", len(P.gates), P)
                break
        if best is None:
            print(f"{name:32s} no SFX23 circuit")
            continue
        open(os.path.join(out_dir, name + ".sfx.txt"), "w").write(dbs_text(best[2]))
        print(f"{name:32s} {best[1]:4d} XOR2 gates from {best[0]}")


if __name__ == "__main__":
    main()
