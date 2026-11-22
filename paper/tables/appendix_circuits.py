"""Appendix A of the paper (supplementary material): the record circuits of Table 3 and the circuits of its
"ours" UMC areas, with their matrices.

A.1: for each problem with a new XOR2 record in Table 3, the XOR2 circuit, with the previous record of the problem.
A.2: for each problem whose UMC area
in Table 3 is below the baseline, our circuit of that area (XOR2/XOR3), and, for comparison, the area after merging
of our circuit with the fewest XOR2 gates (the circuit of A.1, or results/<matrix>.xor2.circ.txt), with the
baseline of the problem (Sec. 5.2). The script reads the matrices (bench/)
and the circuits (results/), checks that each circuit computes its matrix, its gate counts and area, and its
depth: for an XOR2 circuit, at most D; for an XOR2/XOR3 circuit, at most D gate levels, and at most D for the
XOR2 circuit obtained by splitting each XOR3 gate a + b + c into a + (b + c), where (b, c) is the pair of the
absorbed gate (the last two operands, as written by slpd and paper/tables/stm_merge.py).

usage: python3 paper/tables/appendix_circuits.py OUT.tex (the paper needs \\usepackage{needspace})
Notation of the listing (as in the paper): inputs x0..x(n-1) (x_c = column c of the matrix), outputs
y0..y(m-1) (y_j = row j), other gates g1, g2, ... in the order of the circuit; row j of the matrix is printed
as its 0/1 entries, the coefficients of x0, x1, ..., x(n-1) from left to right (as matrices are written in
papers), in groups of 4 (16 inputs) or 8 (32 inputs). Each printed row is compared with the tokens of the
bench/ file. (The NIST .circ files in results/ store each row as the integer sum_c M[j][c] 2^c in their
headers, as the code of Shi, Feng and Xu and NIST's linear-map files do.)"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
from slpd.circuit import read_any          # noqa: E402
from slpd.matrix import read_matrices      # noqa: E402

# A.1: (label, size, matrix file, circuit file, depth bound, XOR2 gates, note)
S4 = "From stage~4 started from the circuit of SFX23 (Sec.~\\ref{sec:setup})."
AR = "From stage~4 with the UMC area as the objective (Sec.~\\ref{sec:setup})."
RECORDS = [
    ("SmallScale AES", 16, "SmallAES", "SmallAES.xor2.circ.txt", 3, 45, S4),
    ("Joltik", 16, "JOLTIK", "JOLTIK.xor2.circ.txt", 3, 46,
     "[SKOP15] involutory $16\\times 16$ is the same matrix, and [LW16] involutory $16\\times 16$ is a row "
     "and column permutation of it (Sec.~\\ref{sec:equiv}). As noted in Table~\\ref{tab:results}, the GitHub "
     "repository~\\cite{Joe26} also lists a 46-gate circuit of depth~3."),
    ("[SS16] involutory", 16, "ToSC_SarSye16_16_Involutory", "ToSC_SarSye16_16_Involutory.xor2.circ.txt", 3, 40, ""),
    ("[SKOP15] Hadamard", 32, "FSE_SKOP15_32_Hadamard", "FSE_SKOP15_32_Hadamard.xor2.circ.txt", 3, 94, S4),
    ("[LS16] circulant", 32, "FSE_LiuSim16_32_Circulant_SFX", "FSE_LiuSim16_32_Circulant_SFX.xor2.circ.txt", 3, 106,
     "The matrix of SFX23 (Sec.~\\ref{sec:ls16}). " + S4),
    ("[LW16]", 32, "FSE_LiWang16_32", "FSE_LiWang16_32.xor2.circ.txt", 3, 98, ""),
    ("[BKL16] circulant", 32, "C_BeiKraLea16_32_Circulant", "C_BeiKraLea16_32_Circulant.xor2.circ.txt", 3, 104, ""),
    ("[SS16] Toeplitz", 32, "ToSC_SarSye16_32_Toeplitz", "ToSC_SarSye16_32_Toeplitz.xor2.circ.txt", 3, 101, S4),
    ("[JPST17] subfield", 32, "ePrint_JeaPeySim_32_Subfield", "ePrint_JeaPeySim_32_Subfield.xor2.circ.txt", 3, 86, ""),
    ("[SKOP15] involutory", 32, "FSE_SKOP15_32_Involutory", "FSE_SKOP15_32_Involutory.xor2.circ.txt", 3, 94, ""),
    ("[LW16] involutory", 32, "FSE_LiWang16_32_Involutory", "FSE_LiWang16_32_Involutory.xor2.circ.txt", 3, 87, ""),
    ("[SS16] involutory", 32, "ToSC_SarSye16_32_Involutory", "ToSC_SarSye16_32_Involutory.xor2.circ.txt", 4, 94, ""),
    ("[KLSW17]", 32, "KLSW", "KLSW.xor2.circ.txt", 3, 80,
     "The previous record comes from the circuit of SFX23 for [LSL+19] HW160~\\#93 (Sec.~\\ref{sec:equiv})."),
    ("[LSL+19] HW168 \\#130", 32, "168_130", "168_130.xor2.circ.txt", 3, 82, ""),
]
# A.2: (label, size, matrix file, circuit file, depth bound, UMC area of Table 3, XOR2 gates before merging, note)
AREAS = [
    ("SmallScale AES", 16, "SmallAES", "SmallAES.umc180.circ.txt", 3, 118.39, 47, ""),
    ("Joltik", 16, "JOLTIK", "JOLTIK.umc180.circ.txt", 3, 122.72, 48,
     "The same area for [SKOP15] and [LW16] involutory $16\\times 16$ (Table~\\ref{tab:results}). " + AR),
    ("[SKOP15] Hadamard", 16, "FSE_SKOP15_16_Hadamard", "FSE_SKOP15_16_Hadamard.umc180.circ.txt", 3, 122.38, 47, ""),
    ("[SS16] Toeplitz", 16, "ToSC_SarSye16_16_Toeplitz", "ToSC_SarSye16_16_Toeplitz.umc180.circ.txt", 3, 111.37, 42, ""),
    ("[JPST17]", 16, "ePrint_JeaPeySim_16", "ePrint_JeaPeySim_16.umc180.circ.txt", 3, 111.71, 43, ""),
    ("[SS16] involutory", 16, "ToSC_SarSye16_16_Involutory", "ToSC_SarSye16_16_Involutory.umc180.circ.txt", 3, 109.05, 43, ""),
    ("[JPST17] involutory", 16, "ePrint_JeaPeySim_16_Involutory", "ePrint_JeaPeySim_16_Involutory.umc180.circ.txt", 3, 121.05, 47, ""),
    ("AES MixColumns", 32, "AES", "AES.umc180.circ.txt", 3, 258.77, 100,
     "Column $8i+k$ (row $8i+k$) is bit~$k$ of input (output) byte~$i$, bit~0 being the least significant "
     "one. " + AR),
    ("[SKOP15] Hadamard", 32, "FSE_SKOP15_32_Hadamard", "FSE_SKOP15_32_Hadamard.umc180.circ.txt", 3, 244.76, 94, S4),
    ("[LS16] circulant", 32, "FSE_LiuSim16_32_Circulant_SFX", "FSE_LiuSim16_32_Circulant_SFX.umc180.circ.txt", 3, 274.79, 108, ""),
    ("[LW16]", 32, "FSE_LiWang16_32", "FSE_LiWang16_32.umc180.circ.txt", 3, 247.45, 98, S4),
    ("[BKL16] circulant", 32, "C_BeiKraLea16_32_Circulant", "C_BeiKraLea16_32_Circulant.umc180.circ.txt", 3, 264.12, 104, ""),
    ("[SS16] Toeplitz", 32, "ToSC_SarSye16_32_Toeplitz", "ToSC_SarSye16_32_Toeplitz.umc180.circ.txt", 3, 259.45, 102, S4),
    ("[JPST17] subfield", 32, "ePrint_JeaPeySim_32_Subfield", "ePrint_JeaPeySim_32_Subfield.umc180.circ.txt", 3, 223.42, 86, AR),
    ("[SKOP15] involutory", 32, "FSE_SKOP15_32_Involutory", "FSE_SKOP15_32_Involutory.umc180.circ.txt", 3, 250.42, 95, ""),
    ("[LW16] involutory", 32, "FSE_LiWang16_32_Involutory", "FSE_LiWang16_32_Involutory.umc180.circ.txt", 3, 243.06, 89, S4),
    ("[SS16] involutory", 32, "ToSC_SarSye16_32_Involutory", "ToSC_SarSye16_32_Involutory.umc180.circ.txt", 4, 244.11, 96, AR),
    ("[KLSW17]", 32, "KLSW", "KLSW.umc180.circ.txt", 3, 219.40, 82, S4),
    ("[LSL+19] HW168 \\#130", 32, "168_130", "168_130.umc180.circ.txt", 3, 221.07, 83, S4),
]


# Our fewest XOR2 count for the problems of A.2 without a new record (Table 3, column "ours"); their circuit of
# fewest XOR2 gates is results/<matrix>.xor2.circ.txt, as for the problems of A.1.
FEWEST = {"FSE_SKOP15_16_Hadamard": 47, "ToSC_SarSye16_16_Toeplitz": 42, "ePrint_JeaPeySim_16": 43,
          "ePrint_JeaPeySim_16_Involutory": 46, "AES": 97}

# The previous record of each problem of A.1, as in Table 3 (paper/tables/make_tables.py: ref_x): the count of Table 5
# of SFX23 (bench/benchmark.tsv, column sfx23_xor2), except for [KLSW17]: the 84-gate circuit of SFX23 for
# [LSL+19] HW160 #93, a row and column permutation of [KLSW17].
PREV = {"SmallAES": 46, "JOLTIK": 47, "ToSC_SarSye16_16_Involutory": 42, "FSE_SKOP15_32_Hadamard": 96,
        "FSE_LiuSim16_32_Circulant_SFX": 110, "FSE_LiWang16_32": 101, "C_BeiKraLea16_32_Circulant": 107,
        "ToSC_SarSye16_32_Toeplitz": 105, "ePrint_JeaPeySim_32_Subfield": 89, "FSE_SKOP15_32_Involutory": 98,
        "FSE_LiWang16_32_Involutory": 89, "ToSC_SarSye16_32_Involutory": 98, "KLSW": 84, "168_130": 84}
PREV_FROM_TSV_EXCEPT = {"KLSW"}
# The baseline of each problem of A.2, as in Table 3 (make_tables.py: ref_a): the smallest UMC area of the circuits of
# the HILL repository and of the published circuits of SFX23 after a maximum merge, with its source.
BASE = {"SmallAES": (127.70, "HILL"), "JOLTIK": (125.38, "HILL"), "FSE_SKOP15_16_Hadamard": (125.04, "SFX23, HILL"),
        "ToSC_SarSye16_16_Toeplitz": (114.03, "SFX23, HILL"), "ePrint_JeaPeySim_16": (113.04, "SFX23, HILL"),
        "ToSC_SarSye16_16_Involutory": (110.04, "SFX23, HILL"), "ePrint_JeaPeySim_16_Involutory": (123.37, "SFX23, HILL"),
        "AES": (264.09, "SFX23"), "FSE_SKOP15_32_Hadamard": (245.44, "SFX23, HILL"),
        "FSE_LiuSim16_32_Circulant_SFX": (283.45, "SFX23"), "FSE_LiWang16_32": (261.77, "SFX23, HILL"),
        "C_BeiKraLea16_32_Circulant": (278.44, "SFX23, HILL"), "ToSC_SarSye16_32_Toeplitz": (269.78, "SFX23"),
        "ePrint_JeaPeySim_32_Subfield": (240.40, "SFX23, HILL"), "FSE_SKOP15_32_Involutory": (256.76, "SFX23, HILL"),
        "FSE_LiWang16_32_Involutory": (247.42, "HILL"), "ToSC_SarSye16_32_Involutory": (263.41, "SFX23"),
        "KLSW": (230.72, "SFX23, HILL, for HW160~\\#93"), "168_130": (229.39, "SFX23, HILL")}


def xor2_values(C):
    """the set of gate values of the XOR2 circuit in which each XOR3 gate [a, b, c] is a + (b + c)"""
    val = [1 << c for c in range(C.n)]
    out = set()
    for ops in C.gates:
        if len(ops) == 3:
            h = val[ops[1]] ^ val[ops[2]]
            out.add(h)
            v = val[ops[0]] ^ h
        else:
            v = val[ops[0]] ^ val[ops[1]]
        out.add(v)
        val.append(v)
    return out


def comparison(mat, D, M, C, in_a1):
    """the sentence that compares an A.2 circuit with our circuit of fewest XOR2 gates (A.1, or results/)"""
    F = read_any(os.path.join(ROOT, "results", mat + ".xor2.circ.txt"))
    g = len(F.gates)
    if F.check(M, D) or F.counts()[1] or F.depth() > D or (not in_a1 and g != FEWEST[mat]):
        sys.exit(f"{mat}: the circuit of fewest XOR2 gates does not check")
    K = F.merge_xor3()
    n2, n3 = K.counts()
    a = round(3.00 * n2 + 4.67 * n3, 2)
    same = xor2_values(C) == xor2_values(F)
    if same:
        assert K.counts() == C.counts()
        txt = ("It is the circuit of Appendix~\\ref{app:records}, merged." if in_a1
               else "It is our circuit with the fewest XOR2 gates, merged.")
    elif in_a1:
        txt = (f"The circuit of Appendix~\\ref{{app:records}}, with {g} XOR2 gates, gives [{n2} XOR2 + {n3} XOR3] "
               f"gates, i.e., {a:.2f}\\,GE.")
    else:
        txt = (f"Our circuit with the fewest XOR2 gates, {g}, gives [{n2} XOR2 + {n3} XOR3] gates, i.e., "
               f"{a:.2f}\\,GE.")
    return txt, (g, n2, n3, a, same)


def names(C):
    """display names: inputs x0.., outputs y0.. (first output that a node computes), others g1.."""
    out_of = {}
    for j, o in enumerate(C.out):
        out_of.setdefault(o, j)
    nm, t = {}, 0
    for c in range(C.n):
        nm[c] = f"x{c}"
    for k in range(len(C.gates)):
        node = C.n + k
        if node in out_of:
            nm[node] = f"y{out_of[node]}"
        else:
            t += 1
            nm[node] = f"g{t}"
    return nm


def xor2_depth(C):
    """the depth of the XOR2 circuit in which each XOR3 gate [a, b, c] is a + (b + c)"""
    dep = [0] * C.n
    for ops in C.gates:
        d = [dep[o] for o in ops]
        dep.append(max(d) + 1 if len(ops) == 2 else max(d[0], max(d[1], d[2]) + 1) + 1)
    return max(dep[o] for o in C.out)


def listing(C):
    nm = names(C)
    lines = []
    for k, ops in enumerate(C.gates):
        rhs = f"{nm[ops[0]]} + {nm[ops[1]]}" if len(ops) == 2 else f"{nm[ops[0]]} + ({nm[ops[1]]} + {nm[ops[2]]})"
        lines.append(f"{nm[C.n + k]} = {rhs}")
    first = {}
    for j, o in enumerate(C.out):          # outputs that are inputs, or equal to an earlier output
        if o < C.n or first.setdefault(o, j) != j:
            lines.append(f"y{j} = {nm[o]}")
    return lines


def plain(s):
    """the text of a LaTeX paragraph, roughly, to estimate its number of lines"""
    s = re.sub(r"\\(paragraph|ref)\{", "{", s)
    s = s.replace("\\times", "x").replace("\\,", "").replace("\\dots", "...").replace("\\#", "#").replace("~", " ")
    return re.sub(r"[{}$\\]", "", s)


def block(label, size, M, C, D, head, note, matrix_from=None, cols=5, raw=None, section=None):
    """the lines of one entry; it starts with \\needspace (package needspace), so that its heading paragraph is not
    left at the bottom of a page without its listing: the lines of the paragraph (about 85 characters per line),
    plus 6 (multicols starts a new page when less than 50pt are left), plus 3 for the heading of a subsection"""
    out = [f"\\paragraph{{{label} (${size}\\times {size}$), $D={D}$: {head}.}}" + (f" {note}" if note else "")]
    if matrix_from:
        out[0] += f" Matrix: as in Appendix~\\ref{{{matrix_from}}}."      # one line per paragraph, as in the paper
    else:
        grp = 8 if M.n > 16 else 4
        bits = ["".join(str((r >> c) & 1) for c in range(M.n)) for r in M.rows]     # x0 first
        for j, b in enumerate(bits):            # read back against the tokens of the bench/ file
            assert b == "".join(raw[j]), f"row {j}"
        width = len(str(M.m - 1))
        rows = [f"y{j}" + "~" * (width - len(str(j))) + "~~" + " ".join(b[i:i + grp] for i in range(0, M.n, grp))
                for j, b in enumerate(bits)]
        out[0] += " Matrix rows $y_0,\\dots,y_{{{}}}$:".format(M.m - 1)
        out.append(f"\\begin{{multicols}}{{{4 if M.n <= 16 else 2}}}\\scriptsize\\ttfamily\\raggedright")
        out.append("\\\\\n".join("\\mbox{" + x + "}" for x in rows))
        out.append("\\end{multicols}")
    out.append(f"\\begin{{multicols}}{{{cols}}}\\scriptsize\\ttfamily\\raggedright")
    out.append("\\\\\n".join("\\mbox{" + x + "}" for x in listing(C)))
    out.append("\\end{multicols}")
    out.append("")
    k = next(i for i, x in enumerate(out) if x.startswith("\\begin{multicols}"))
    n = math.ceil(len(plain(" ".join(out[:k]))) / 85) + 6 + (3 if section else 0)
    return [f"\\needspace{{{n}\\baselineskip}}"] + ([section] if section else []) + out


def sfx23_counts():
    """the XOR2 counts of Table 5 of SFX23, from bench/benchmark.tsv"""
    rows = [line.rstrip("\n").split("\t") for line in open(os.path.join(ROOT, "bench", "benchmark.tsv"))]
    head = rows[0]
    i, j = head.index("file"), head.index("sfx23_xor2")
    return {r[i][:-4]: int(r[j]) for r in rows[1:] if len(r) > j and r[j]}


def load(mat, circ):
    path = os.path.join(ROOT, "bench", mat + ".txt")
    M = read_matrices(path)[0]
    C = read_any(os.path.join(ROOT, "results", circ))
    raw = [line.split() for line in open(path) if line.strip() and not line.startswith("#")]
    return M, C, raw


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    out = ["% Written by paper/tables/appendix_circuits.py from bench/ and results/ of the repository."]
    section = "\\subsection{XOR2 Circuits of the New Records}\\label{app:records}"
    listed = set()
    for label, size, mat, circ, D, g, note in RECORDS:
        M, C, raw = load(mat, circ)
        problems = C.check(M, D)
        n2, n3 = C.counts()
        if problems or n3 or len(C.gates) != g or C.depth() > D:
            sys.exit(f"{circ}: {problems[:3]}, {len(C.gates)} gates ({n3} XOR3), depth {C.depth()}; expected {g}, D = {D}")
        prev = PREV[mat]
        tsv = sfx23_counts()
        if g >= prev or (mat not in PREV_FROM_TSV_EXCEPT and prev != tsv[mat]):
            sys.exit(f"{mat}: previous record {prev} (Table 5 of SFX23: {tsv.get(mat)}), ours {g}")
        out += block(label, size, M, C, D, f"{g} XOR2 gates, previous record {prev}", note, raw=raw, section=section)
        section = None
        listed.add(mat)
        print(f"A.1 {label}: {g} XOR2 gates, depth {C.depth()}, checked")
    section = "\\subsection{Circuits of Our UMC Areas}\\label{app:areas}"
    for label, size, mat, circ, D, area, g, note in AREAS:
        M, C, raw = load(mat, circ)
        problems = C.check(M, D)
        n2, n3 = C.counts()
        a = round(3.00 * n2 + 4.67 * n3, 2)
        if problems or abs(a - area) > 1e-9 or n2 + 2 * n3 != g or C.depth() > D or xor2_depth(C) > D:
            sys.exit(f"{circ}: {problems[:3]}, [{n2} + {n3}] = {a} GE from {n2 + 2 * n3} gates, gate levels {C.depth()}, "
                     f"XOR2 depth {xor2_depth(C)}; expected {area} GE from {g} gates, D = {D}")
        base, src = BASE[mat]
        if not area < base:
            sys.exit(f"{mat}: area {area} not below the baseline {base}")
        head = f"[{n2} XOR2 + {n3} XOR3] gates, {area:.2f}\\,GE, from {g} XOR2 gates; baseline {base:.2f}\\,GE ({src})"
        cmp_txt, cmp = comparison(mat, D, M, C, mat in listed)
        out += block(label, size, M, C, D, head, (note + " " if note else "") + cmp_txt,
                     matrix_from="app:records" if mat in listed else None, cols=4, raw=raw, section=section)
        section = None
        print(f"A.2 {label}: [{n2} + {n3}] = {a:.2f} GE from {g} gates, XOR2 depth {xor2_depth(C)}, checked; "
              f"fewest XOR2: {cmp[0]} -> [{cmp[1]} + {cmp[2]}] = {cmp[3]:.2f} GE{' (same circuit)' if cmp[4] else ''}")
    with open(sys.argv[1], "w") as f:
        f.write("\n".join(out).rstrip("\n") + "\n")


if __name__ == "__main__":
    main()
