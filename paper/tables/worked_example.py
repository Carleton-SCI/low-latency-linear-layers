"""The worked example of Sec. 3 of the paper: one round of each stage of the search on the 8x8 matrix
y_j = x_j + x_{j+1} + ... + x_{j+4} (indices mod 8), at depth D = 3.

The distances of Tables 1 and 2 are computed with the forward search of paper/example/toy_engine.py, which follows the
rules of src/dbs.cpp (candidate gates, exact distances under the depth bound, norm, MD filter, backward splits);
the facts of the runs and the two listed circuits come from paper/example/toy_data.json (written by toy_data.py), and
the first two steps of its stage-3 run (B = 4, K = 8) are replayed here. The script checks every fact that the
example of Sec. 3 states (the asserts below), the two circuits against the matrix and the depth bound, and their
gate counts and areas.

usage: python3 paper/tables/worked_example.py OUT_DIR
writes OUT_DIR/example_matrix.tex (the matrix display), example_iterations.tex (Table 1: iterations 1 and 2 of the
forward search, with the tie rules and the MD filter), example_splits.tex (Table 2: two backward splits) and
example_circuits.tex (the 16-gate circuit of stage 1 and the 17-gate circuit of stage 1 with s = 1, merged), as in
the paper. Notation of the listings as in Appendix A (paper/tables/appendix_circuits.py)."""
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "paper", "example"))
from slpd.circuit import Circuit           # noqa: E402
from slpd.matrix import Matrix             # noqa: E402
from appendix_circuits import listing      # noqa: E402
from toy_engine import Solver, prefix_init, apply_split, dmin, split_beam   # noqa: E402

N, D = 8, 3
ROWS = [sum(1 << ((j + i) % N) for i in range(5)) for j in range(N)]
HEAD = "% Written by paper/tables/worked_example.py of the repository."


def v(*cs):
    return sum(1 << c for c in cs)


def tex(x):
    """a value as a sum of inputs, e.g. x_0 \\oplus x_1"""
    return "\\oplus ".join(f"x_{c}" for c in range(N) if x >> c & 1)


def after(s, c):
    """the distances of the targets after candidate c, before the pre-emptive gates"""
    return [s.delta[j] - (c.red >> j & 1) for j in range(len(s.T))]


def iterations_table(it1, it2):
    """Table 1: one row per situation (the base, then a candidate gate), one column per output, then the sum, the
    norm and, for iteration 2, whether the MD filter keeps the gate; it1, it2: lists of (label, distances, MD)"""
    out = [r"\begin{table}[tp]", r"\centering",
           r"\caption{Distances $\delta$ of the outputs of the example in iterations~1 and~2 of the forward search, "
           r"with the base alone and with four or three candidate gates, and their sum $\Delta$ and norm "
           r"$\sum\delta^2$. In iteration~1, the base holds the inputs, and the other 24 candidate gates are "
           r"rotations of the four shown. In iteration~2, the base also holds $x_0\oplus x_1$; $x_5\oplus x_6$ and "
           r"$x_6\oplus x_7$, the other two ties, are the mirror images of $x_3\oplus x_4$ and $x_2\oplus x_3$. MD: "
           r"whether the MD filter keeps the gate.}",
           r"\label{tab:exiter}", r"\small", r"\setlength{\tabcolsep}{4pt}",
           r"\begin{tabular}{@{}l l cccccccc rr l@{}}", r"\toprule",
           " & & " + " & ".join(f"$y_{j}$" for j in range(N)) + r" & $\Delta$ & $\sum\delta^2$ & MD \\",
           r"\midrule"]
    for k, (name, rows) in enumerate((("Iteration~1", it1), ("Iteration~2", it2))):
        if k:
            out.append(r"\midrule")
        for i, (label, dist, md) in enumerate(rows):
            vals = [str(x) for x in dist] + [str(sum(dist)), str(sum(x * x for x in dist))]
            out.append((name + " " if i == 0 else "") + f"& {label} & " + " & ".join(vals) + f" & {md} \\\\")
    return out + [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]


def splits_table(rows):
    """Table 2: one row per prefix of splits, one column per target (the outputs, then the part p), then the sum"""
    out = [r"\begin{table}[tp]", r"\centering",
           r"\caption{Targets of the forward search in the example and their distances $\delta$: without splits, after "
           r"Rule~5 splits $y_1$ into $x_5$ and $p=x_1\oplus x_2\oplus x_3\oplus x_4$, and after Rule~3 splits $y_0$ "
           r"into $p$ and $x_0$. The outputs have depth bound~3 and $p$ has depth bound~2; --: not a target. The parts "
           r"$x_5$ and $x_0$ are inputs ($\delta=0$).}",
           r"\label{tab:exsplits}", r"\small", r"\setlength{\tabcolsep}{4pt}",
           r"\begin{tabular}{@{}l cccccccc c r@{}}", r"\toprule",
           " & " + " & ".join(f"$y_{j}$" for j in range(N)) + r" & $p$ & $\Delta$ \\", r"\midrule"]
    for label, vals in rows:
        out.append(f"{label} & " + " & ".join(str(x) for x in vals) + r" \\")
    return out + [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    outdir = sys.argv[1]
    M = Matrix(ROWS, N, "window5")
    assert M.min_depth() == D and M.rank() == N
    d = json.load(open(os.path.join(ROOT, "paper", "example", "toy_data.json")))
    assert d["rows"] == ROWS and (d["n"], d["D"]) == (N, D)
    files = {}

    # ------------------------------------------------------------ the matrix
    body = " \\\\\n".join(" & ".join(str(r >> c & 1) for c in range(N)) for r in ROWS)
    files["matrix"] = [r"\[", r"\begin{bmatrix} y_0\\ y_1\\ \vdots\\ y_7 \end{bmatrix}=",
                       r"\begin{bmatrix}", body, r"\end{bmatrix}",
                       r"\begin{bmatrix} x_0\\ x_1\\ \vdots\\ x_7 \end{bmatrix}.", r"\]", ""]

    # ------------------------------------------------------------ iteration 1 (Table 1, tab:exiter)
    s = Solver(N, ROWS, [1 << D] * N, D, "all", False)
    s.init()
    assert s.delta == [4] * N and not s.log
    c1 = {c.v: c for c in s.score()}
    assert len(c1) == 28
    show = [v(0, 1), v(0, 2), v(0, 3), v(0, 4)]
    assert [c1[x].S for x in show] == [28, 29, 30, 30] and [c1[x].N for x in show] == [100, 107, 114, 114]
    assert sorted(c.S for c in c1.values()).count(28) == 8 and all(c1[v(i, (i + 1) % N)].S == 28 for i in range(N))
    assert {(c.S, c.N) for c in c1.values() if c.S == 28} == {(28, 100)}
    assert sorted(c.S for c in c1.values()).count(29) == 8 and sorted(c.S for c in c1.values()).count(30) == 12
    assert [j for j in range(N) if c1[v(0, 1)].red >> j & 1] == [0, 5, 6, 7]
    cols = [s.delta] + [after(s, c1[x]) for x in show]
    it1 = [("base", list(s.delta), "")] + [(f"${tex(x)}$", after(s, c1[x]), "") for x in show]   # copies
    assert [sum(x) for _, x, _ in it1] == [32, 28, 29, 30, 30]

    # ------------------------------------------------------------ iteration 2, after x0 + x1 (Table 1)
    s.apply(c1[v(0, 1)])
    assert s.delta == [3, 4, 4, 4, 4, 3, 3, 3] and not s.log
    c2 = {c.v: c for c in s.score()}
    m = Solver(N, ROWS, [1 << D] * N, D, "sfx", True)
    m.init()
    m.apply({c.v: c for c in m.score()}[v(0, 1)])
    c2md = {c.v: c for c in m.score()}
    ties = sorted((c.N, c.v) for c in c2.values() if c.S == 24)
    assert len(c2) == 29 and len(c2md) == 15 and min(c.S for c in c2.values()) == 24
    assert ties == sorted([(72, v(4, 5)), (74, v(3, 4)), (74, v(5, 6)), (76, v(2, 3)), (76, v(6, 7))])
    assert sum(1 for c in c2.values() if c.S <= 25) == 11
    s_bp = {x for n_, x in ties if n_ == max(t[0] for t in ties)}
    s_sfx = {x for n_, x in ties if n_ == min(t[0] for t in ties)}
    tmd = [(c.N, c.v) for c in c2md.values() if c.S == min(c.S for c in c2md.values())]
    s_sfxmd = {x for n_, x in tmd if n_ == min(t[0] for t in tmd)}
    assert s_bp == {v(2, 3), v(6, 7)} and s_sfx == {v(4, 5)} and s_sfxmd == {v(3, 4), v(5, 6)}
    assert v(4, 5) not in c2md and all(x in c2md for x in (v(2, 3), v(3, 4), v(5, 6), v(6, 7)))
    show = [v(2, 3), v(3, 4), v(4, 5)]
    cols = [s.delta] + [after(s, c2[x]) for x in show]
    it2 = [("base", list(s.delta), "")] + [(f"${tex(x)}$", after(s, c2[x]), "kept" if x in c2md else "removed")
                                     for x in show]
    assert [sum(x) for _, x, _ in it2] == [28, 24, 24, 24] and [m_ for _, _, m_ in it2[1:]] == ["kept", "kept", "removed"]
    files["iterations"] = iterations_table(it1, it2)

    # ------------------------------------------------------------ two backward splits (Table 2, tab:exsplits)
    p = v(1, 2, 3, 4)
    sp1, sp2 = (ROWS[1], v(5), p), (ROWS[0], p, v(0))
    st = d["bpbs"]["first16"]["steps"]
    assert [(x["rule"], tuple(x["split"])) for x in st[:2]] == [(5, sp1), (3, sp2)]
    assert (st[2]["rule"], tuple(st[2]["split"]), st[2]["result"]) == (5, (ROWS[5], v(1), v(0, 5, 6, 7)), 16)
    pf = prefix_init(ROWS, D)
    assert sorted(pf.Wv) == sorted(ROWS) and not pf.P              # Rule 1 applies to no output

    def applicable(pf):
        """the pairs (output, target within depth D-1) to which Rules 2 and 3 apply (draw_split of src/dbs.cpp)"""
        r2 = [(w, x) for w in pf.Wv for i, (x, _) in enumerate(pf.P) if pf.inP(w ^ x) > i]
        r3 = [(w, x) for w in pf.Wv for x, _ in pf.P if w ^ x and dmin(w ^ x) < D]
        return r2, r3
    # Sec. 3.3, Rule 5: a part has 1 to 4 of the 5 terms (the sizes of random_split, with the budget 2^(D-1) = 4)
    assert all(r.bit_count() == 5 for r in ROWS) and (max(1, 5 - (1 << (D - 1))), min(1 << (D - 1), 5 - 1)) == (1, 4)
    # Sec. 3.3: in the beginning, no target can be built within depth 2, so Rules 2 and 3 do not apply
    assert all(dmin(r) == D for r in ROWS) and applicable(pf) == ([], [])

    def dist(pf):
        T, TB = pf.targets(D)
        t = Solver(N, T, TB, D, "all", False)
        assert t.init()
        return {x: (b.bit_length() - 1, t.delta[i]) for i, (x, b) in enumerate(zip(T, TB))}
    d0 = dist(pf)
    apply_split(pf, sp1)
    d1 = dist(pf)
    # Sec. 3.3: the parts x5 and p are the only targets within depth 2; Rule 2 applies to no output (x5 + p = y1 is
    # split), and Rule 3 applies to 7 pairs, among them y0 and p (y0 + p = x0)
    r2, r3 = applicable(pf)
    assert sorted(x for x, _ in pf.P) == sorted([v(5), p]) and ROWS[1] == v(5) ^ p and ROWS[1] not in pf.Wv
    assert not r2 and len(r3) == 7 and (ROWS[0], p) in r3 and ROWS[0] ^ p == v(0)
    assert dmin(ROWS[0] ^ p) < D                                   # Rule 3 applies to y0 with p
    apply_split(pf, sp2)
    d2 = dist(pf)
    # Sec. 3.3, stage 2: after the two splits, Rule 2 still applies to no output and Rule 3 applies again, so that the
    # third split of the trajectory (Rule 5 on y5) was drawn with probability 1 - q; results 17, 17, 16
    r2, r3 = applicable(pf)
    assert not r2 and r3 and sorted(x for x, _ in pf.P) == sorted([v(5), p, v(0)])
    assert [x["result"] for x in st] == [17, 17, 16]
    # Sec. 3.3, stage 3: with B = 4 and K = 8 (paper/example/toy_data.py), the initial prefix draws its splits by Rule 5 (see
    # above), and 4 prefixes are kept after steps 1 and 2; the replay of these two steps reproduces the record
    sb = d["split_beam"]["steps"]
    assert [(x["step"], x["prefixes"], x["splits"]) for x in sb[:2]] == [(1, 4, [1] * 4), (2, 4, [2] * 4)]
    rec = []
    split_beam(N, ROWS, D, 1, 2, 4, 8, 2, 0.5, "sfx", True, random.Random(71), rec)
    assert rec == sb[:2]
    assert [sum(x[1] for x in dd.values()) for dd in (d0, d1, d2)] == [32, 31, 27]
    assert d1[p] == (2, 3) and d1[v(5)] == (2, 0) and d2[v(0)] == (2, 0) and ROWS[1] not in d1 and ROWS[0] not in d2
    targets = list(ROWS) + [p]
    files["splits"] = splits_table(
        [(label, [dd[x][1] if x in dd else "--" for x in targets] + [sum(x[1] for x in dd.values())])
         for label, dd in (("No split", d0), (r"Rule~5, $y_1$", d1), (r"Rule~3, $y_0$", d2))])

    # ------------------------------------------------------------ stage 4: the facts of the text
    L = d["lns"]
    assert L["start"]["count"] == 18
    i2, i6 = L["iters"][1], L["iters"][5]
    assert i2["k"] == 2 and sorted(map(tuple, i2["splits"])) == sorted([(ROWS[0], v(0, 1), v(2, 3, 4)),
                                                                        (ROWS[3], v(3, 4), v(5, 6, 7))])
    assert (i2["before"], i2["result"], i2["accepted"]) == (18, 17, True)
    assert i6["k"] == 1 and i6["splits"] == [[ROWS[7], v(3), v(0, 1, 2, 7)]] and i6["result"] == 16
    tg = {v(0, 1): 1, v(2, 3, 4): 2, v(3, 4): 1, v(5, 6, 7): 2}
    for x, dd in tg.items():
        assert x.bit_count() - 1 == dd and x.bit_count() <= 1 << (D - 1)

    # ------------------------------------------------------------ the listings, merging and selection
    def circuit(c):
        return Circuit(N, [[g["a"], g["b"]] for g in c["gates"]], c["out"])
    b = d["beam"]
    assert (b["W"], b["s"], b["results"]) == (4, 0, [16])
    assert [[x["gates"][0][0]] for x in b["iters"][0]["beam"]] == [[v(0, 1)], [v(2, 3)], [v(4, 5)], [v(1, 2)]]
    assert (b["iters"][1]["options"], b["iters"][1]["kept"], b["iters"][1]["Smin"]) == (116, 20, 24)
    C16 = circuit(b["circuit"])
    assert not C16.check(M, D) and len(C16.gates) == 16 and C16.merge_xor3().counts() == (16, 0)
    fan = [0] * (N + 16)
    for ops in C16.gates:
        for o in ops:
            fan[o] += 1
    inner = [N + k for k in range(16) if N + k not in C16.out]
    assert len(inner) == 8 and all(fan[x] == 2 for x in inner)
    e16 = b["circuit"]["stm_exact"]
    assert (e16["n2"], e16["n3"], e16["stm"], len(e16["dup"]), len(e16["merges"])) == (4, 8, 50.6, 4, 0)
    s1 = d["beam_s1"]
    assert s1["results"] == [17] and s1["iters"][0]["kept"] == 16
    c17 = s1["best_umc"]
    C17 = circuit(c17)
    K17 = C17.merge_xor3()
    assert not C17.check(M, D) and len(C17.gates) == 17 and not K17.check(M, D) and K17.counts() == (11, 3)
    assert round(K17.area("umc180"), 2) == c17["umc"] == 47.01
    absorbed = sorted(g["val"] for g in c17["gates"] if g["absorbed_by"] >= 0)
    assert absorbed == sorted([v(1, 6), v(0, 4), v(3, 7)])
    e17 = c17["stm_exact"]
    assert (e17["n2"], e17["n3"], e17["stm"], len(e17["dup"]), len(e17["merges"])) == (5, 7, 49.27, 2, 3)

    def block(head, C):
        """the heading and the listing in one unbreakable box, so that they stay on one page"""
        return [r"\par\medskip\noindent\begin{minipage}{\linewidth}", head,
                r"\begin{multicols}{4}\scriptsize\ttfamily\raggedright",
                "\\\\\n".join("\\mbox{" + x + "}" for x in listing(C)), r"\end{multicols}", r"\end{minipage}"]
    files["circuits"] = (["% Line numbers (submission version) do not work inside multicols: none for the listings.",
                          r"\ifdefined\nolinenumbers\nolinenumbers\fi"] +
                         block(r"\paragraph{The 16-gate circuit of stage~1 (Sec.~\ref{sec:beam}).} Each of its 8 "
                               r"non-output gates, $g_1,\dots,g_8$, is used by two gates.", C16) +
                         block(r"\paragraph{The 17-gate circuit of stage~1 with $s=1$, merged: [11 XOR2 + 3 XOR3] gates, "
                               r"47.01\,GE.} Its XOR3 gates are written $a + (b + c)$, where $b + c$ is the absorbed gate.",
                               K17) + [r"\par\ifdefined\linenumbers\linenumbers\fi", ""])
    for name, lines in files.items():
        with open(os.path.join(outdir, f"example_{name}.tex"), "w") as f:
            f.write("\n".join([HEAD] + lines))
    print("Sec. 3: matrix, Tables 1 and 2, stage 4 and the two circuits; all checks passed")


if __name__ == "__main__":
    main()
