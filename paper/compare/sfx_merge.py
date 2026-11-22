"""Areas of the published circuits of Shi, Feng and Xu (SFX23) after merging, as a baseline for Table 3.

For every benchmark matrix, takes the SFX23 circuits whose matrix equals it up to row and column
permutations (as an isomorphism of bipartite graphs) and whose depth is the minimum depth, and merges
each of them as in the paper: a maximum merge for the UMC areas, and the exact merge with duplication of
gates of fan-out two for the STM areas (paper/tables/stm_merge.py). Reports the smallest areas, with the
source circuit and whether it computes the benchmark matrix itself (same_matrix) or a permutation of it.

usage: python3 paper/compare/sfx_merge.py SFX_RESULTS_DIR BENCH_DIR OUT.json
       SFX_RESULTS_DIR: the folder of result files of the SFX23 repository (our_results/)
requires: networkx, ortools"""
import csv
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "tables"))
from slpd.circuit import read_sfx_text, read_text          # noqa: E402
from slpd.matrix import read_matrices                      # noqa: E402
from stm_merge import best_merge, apply_merge              # noqa: E402
import networkx as nx                                      # noqa: E402


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


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sfx_dir, bench, out = sys.argv[1:]
    circuits = []
    for f in sorted(glob.glob(os.path.join(sfx_dir, "*.txt"))):
        for k, C in enumerate(read_sfx_text(read_text(f), f)):
            C = C.prune()
            circuits.append((os.path.basename(f), k, C, rows_of(C), C.depth()))
    res = {}
    for r in csv.DictReader(open(os.path.join(bench, "benchmark.tsv")), delimiter="\t"):
        name, D = r["file"][:-4], int(r["min_depth"])
        M = read_matrices(os.path.join(bench, r["file"]))[0]
        GM = graph(M.rows, M.n)
        best = None
        for fname, k, C, rows, dep in circuits:
            if C.n != M.n or len(rows) != M.m or dep > D:
                continue
            if sorted(bin(x).count("1") for x in rows) != sorted(bin(x).count("1") for x in M.rows):
                continue
            if not nx.is_isomorphic(GM, graph(rows, C.n), node_match=lambda a, b: a["side"] == b["side"]):
                continue
            umc = C.merge_xor3()
            se, st = best_merge(C)
            stm = apply_merge(C, se, st)
            cand = dict(source=f"{fname}#{k}", xor2=len(C.gates), depth=dep,
                        umc180=round(umc.area("umc180"), 2), stm130=round(stm.area("stm130"), 2),
                        same_matrix=rows == M.rows)      # False: a circuit of a permutation-equivalent matrix
            if best is None:
                best = dict(xor2=cand, umc180=cand, stm130=cand)
            else:
                if cand["xor2"] < best["xor2"]["xor2"]:
                    best["xor2"] = cand
                if cand["umc180"] < best["umc180"]["umc180"]:
                    best["umc180"] = cand
                if cand["stm130"] < best["stm130"]["stm130"]:
                    best["stm130"] = cand
        res[name] = best
        if best:
            print(f"{name:32s} XOR2 {best['xor2']['xor2']:4d} ({best['xor2']['source']})  UMC {best['umc180']['umc180']:7.2f} "
                  f"({best['umc180']['xor2']} gates, {best['umc180']['source']})  STM {best['stm130']['stm130']:7.2f}")
        else:
            print(f"{name:32s} no SFX23 circuit at depth {D}")
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
