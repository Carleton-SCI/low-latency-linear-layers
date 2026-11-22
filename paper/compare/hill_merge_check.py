"""For every matrix in the public HILL repository, compare HILL's published XOR2/XOR3 circuit
(<name>_result.txt) with a plain maximum fan-out-1 XOR3 merge (xor3_merge.merge) of the XOR2 input
circuit that the repository stores next to it (seq.txt, or <name>.txt). Handles the three input
formats found there: 't[a] = t[b] + t[c] y[k]', the same with 'yk', and in-place 't[a] = t[a] ^ t[b]'.
Every circuit is checked against the matrix derived from HILL's result before it is compared.

usage: python3 hill_merge_check.py HILL_REPO"""
import re, os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from xor3_merge import merge

UMC = (3.00, 4.67)
LINE = re.compile(r"\s*t\[(\d+)\]\s*=\s*(.+?)\s*(?:\b(y\[?\d+\]?(?:\s+y\[?\d+\]?)*))?\s*$")

def parse(path):
    """Returns (gates as (name, [operand names]), outs {k: name}, n_inputs). In-place lines
    't[a] = t[a] ^ t[b]' are converted to fresh nodes (single assignment)."""
    raw = []
    for line in open(path):
        if not line.strip() or line.startswith("Optimization"): continue
        m = LINE.match(line)
        if not m: raise ValueError(f"{path}: {line!r}")
        g = int(m.group(1)); ops = [int(x) for x in re.findall(r"t\[(\d+)\]", m.group(2))]
        ys = [int(y) for y in re.findall(r"y\[?(\d+)\]?", m.group(3) or "")]
        raw.append((g, ops, ys, "^" in m.group(2)))
    inplace = any(r[3] for r in raw)
    if inplace:
        n = 1 + max(max(r[0], *r[1]) for r in raw)
        cur = {i: f"x{i}" for i in range(n)}
        gates, outs = [], {}
        for k, (g, ops, ys, _) in enumerate(raw):
            name = f"g{k}"
            gates.append((name, [cur[o] for o in ops])); cur[g] = name
            for y in ys: outs[y] = name
        return gates, outs, n
    assigned = {g for g, *_ in raw}
    ins = {o for _, ops, *_ in raw for o in ops} - assigned
    n = max(ins) + 1
    nm = lambda i: f"x{i}" if i < n else f"g{i}"
    gates = [(nm(g), [nm(o) for o in ops]) for g, ops, _, _ in raw]
    outs = {y: nm(g) for g, _, ys, _ in raw for y in ys}
    return gates, outs, n

def evaluate(gates, outs, n):
    v = {f"x{i}": 1 << i for i in range(n)}; d = {f"x{i}": 0 for i in range(n)}
    for g, ops in gates:
        x = 0
        for o in ops: x ^= v[o]
        v[g] = x; d[g] = max(d[o] for o in ops) + 1
    return {k: v[outs[k]] for k in outs}, max(d[outs[k]] for k in outs)

def area(gates):
    n2 = sum(len(o) == 2 for _, o in gates); n3 = sum(len(o) == 3 for _, o in gates)
    assert all(len(o) <= 3 for _, o in gates)
    return n2, n3, UMC[0] * n2 + UMC[1] * n3

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    repo = sys.argv[1]
    same = better = worse = 0
    for d in sorted(glob.glob(f"{repo}/results/*/")):
        name = os.path.basename(d.rstrip("/"))
        res = glob.glob(f"{d}*_result.txt")
        inp = [p for p in glob.glob(f"{d}*.txt") if not p.endswith("_result.txt")]
        if len(res) != 1 or len(inp) != 1:
            print(f"{name:32s} result/input files: {len(res)}/{len(inp)}, skipped"); continue
        hg, ho, hn = parse(res[0]); hv, hd = evaluate(hg, ho, hn)
        sg, so, sn = parse(inp[0]); sv, sd = evaluate(sg, so, max(sn, hn))
        if sv != hv:
            print(f"{name:32s} input circuit computes a different matrix, skipped"); continue
        mg = merge(sg, so); mv, md = evaluate(mg, so, max(sn, hn))
        assert mv == hv and md <= sd
        h2, h3, ha = area(hg); m2, m3, ma = area(mg)
        tag = "same" if abs(ha - ma) < 0.005 else ("merge better" if ma < ha else "HILL better")
        same += tag == "same"; better += tag == "merge better"; worse += tag == "HILL better"
        print(f"{name:32s} input {len(sg):3d} XOR2 (depth {sd}) | HILL {h2:3d}+{h3:2d} = {ha:7.2f} (depth {hd}) | "
              f"merge {m2:3d}+{m3:2d} = {ma:7.2f} (depth {md}) | {tag}")
    print(f"\nsame area: {same}, merge smaller: {better}, HILL smaller: {worse}")
