"""Derive the exact binary matrices of the HILL benchmark from the circuits in the public HILL
repository (results/<name>/<name>_result.txt): every output y[k] is evaluated symbolically over the
inputs t[0..n-1]. Also re-checks HILL's reported XOR2/XOR3 counts, depth and area, and cross-checks
the matrix against seq.txt (an in-place XOR program for the same matrix) when it parses.

usage: python3 hill_to_matrices.py HILL_REPO OUTDIR"""
import re, os, sys, glob, json

def parse_result(path):
    head, gates, outs = None, [], {}
    for line in open(path):
        line = line.strip()
        if not line: continue
        if line.startswith("Optimization"):
            m = re.search(r"Area:\s*([\d.]+)\).*XOR2:\s*(\d+),\s*XOR3:\s*(\d+).*Depth:\s*(\d+)", line)
            head = (float(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4))) if m else line
            continue
        m = re.match(r"t\[(\d+)\]\s*=\s*(.+?)(?:\s+((?:y\[\d+\]\s*)+))?$", line)
        if not m: raise ValueError(f"{path}: cannot parse {line!r}")
        g = int(m.group(1)); ops = [int(x) for x in re.findall(r"t\[(\d+)\]", m.group(2))]
        gates.append((g, ops))
        if m.group(3):
            for k in re.findall(r"y\[(\d+)\]", m.group(3)): outs[int(k)] = g
    return head, gates, outs

def evaluate(gates, n):
    v = {i: 1 << i for i in range(n)}; dep = {i: 0 for i in range(n)}
    for g, ops in gates:
        x = 0
        for o in ops: x ^= v[o]
        v[g] = x; dep[g] = max(dep[o] for o in ops) + 1
    return v, dep

def n_inputs(gates):
    assigned = {g for g, _ in gates}
    ins = {o for _, ops in gates for o in ops} - assigned
    return max(ins) + 1, ins

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    repo, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    summary = {}
    for d in sorted(glob.glob(f"{repo}/results/*/")):
        name = os.path.basename(d.rstrip("/"))
        res = [p for p in glob.glob(f"{d}*_result.txt")]
        if len(res) != 1: print(f"{name}: {len(res)} result files, skipped"); continue
        head, gates, outs = parse_result(res[0])
        n, ins = n_inputs(gates)                     # inputs: operands that are never assigned
        assert ins == set(range(n)) or ins <= set(range(n)), name
        assert all(g >= n for g, _ in gates), f"{name}: gate overwrites an input"
        v, dep = evaluate(gates, n)
        m = max(outs) + 1 if outs else 0
        missing = [k for k in range(m) if k not in outs]
        rows = [v[outs[k]] if k in outs else None for k in range(m)]
        n2 = sum(len(o) == 2 for _, o in gates); n3 = sum(len(o) == 3 for _, o in gates)
        nw = sum(len(o) > 3 for _, o in gates)
        area = 3.00 * n2 + 4.67 * n3
        odep = max(dep[outs[k]] for k in outs)
        # cross-check with seq.txt (same matrix up to an output permutation?)
        seqok, seq2, seqd = None, None, None
        if os.path.exists(f"{d}seq.txt") and not missing:
            try:
                _, sg, so = parse_result(f"{d}seq.txt")
                sn, _ = n_inputs(sg)
                sv, sd = evaluate(sg, max(sn, n))
                seqok = sorted(so) == list(range(m)) and all(sv[so[k]] == rows[k] for k in range(m))
                seq2 = (sum(len(o) == 2 for _, o in sg), sum(len(o) == 3 for _, o in sg))
                seqd = max(sd[so[k]] for k in so)
            except Exception as e:
                seqok = f"error {e}"
        # matrix: row k = output y[k]; column c = input t[c]
        if not missing:
            with open(f"{out}/{name}.txt", "w") as f:
                for r in rows:
                    f.write(" ".join(str((r >> c) & 1) for c in range(n)) + "\n")
        wts = [bin(r).count("1") for r in rows if r is not None]
        mind = max((w - 1).bit_length() for w in wts) if wts else None
        summary[name] = dict(n=n, m=m, missing=missing, hill_header=head, xor2=n2, xor3=n3, wider=nw,
                             area=round(area, 2), depth=odep, min_depth=mind, d_xor=sum(wts) - m,
                             seq_match=seqok, seq_xor2_xor3=seq2, seq_depth=seqd)
        print(f"{name:32s} {m:3d}x{n:<3d} rows-missing {len(missing)}  HILL {head}  recount XOR2 {n2} XOR3 {n3} "
              f"(>3: {nw}) area {area:.2f} depth {odep} (min {mind}) | seq.txt same matrix {seqok}, XOR2/XOR3 {seq2}, depth {seqd}")
    json.dump(summary, open(f"{out}/summary.json", "w"), indent=1)
