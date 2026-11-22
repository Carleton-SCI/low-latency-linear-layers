"""Post-process a 2-input XOR circuit into XOR2/XOR3 gates (HILL Rule 1): a non-output gate c with
fan-out 1 is absorbed into its consumer p (p = c ^ z, c = a ^ b  ->  p = a ^ b ^ z). An XOR3 gate can
neither absorb another gate nor be absorbed (that would need XOR4), so the merges form a matching in
the forest of (c -> p) edges; a maximum matching is computed exactly (tree DP). Depth never grows.
Reports XOR2/XOR3 counts and GE with the cell areas of HILL (ToSC 2026) Table 1 (UMC 180 nm:
XOR2 3.00, XOR3 4.67; STM 130 nm: XOR2 3.33, XOR3 4.66), and verifies the merged circuit against the
matrix (outputs and depth).

usage: python3 xor3_merge.py matrix.txt circuit.txt [more circuits ...]
       circuits in dbeam format (gK = A + B [d] ... yJ = gK); add --write to save <circuit>.xor3.txt"""
import re, sys
from collections import defaultdict

G2 = {"UMC180": (3.00, 4.67), "STM130": (3.33, 4.66)}

def load_dbeam(path):
    gates, outs = [], {}
    for l in open(path):
        m = re.match(r"(g\d+) = (\w+) \+ (\w+)", l)
        if m: gates.append((m.group(1), [m.group(2), m.group(3)]))
        m = re.match(r"y(\d+) = (\w+)\s*$", l)
        if m: outs[int(m.group(1))] = m.group(2)
    return gates, outs

def merge(gates, outs):
    use = defaultdict(list)
    for g, ops in gates:
        for o in ops: use[o].append(g)
    outset = set(outs.values())
    names = [g for g, _ in gates]
    # forest edges child -> parent
    parent = {g: use[g][0] for g in names if g not in outset and len(use[g]) == 1}
    children = defaultdict(list)
    for c, p in parent.items(): children[p].append(c)
    # tree DP for maximum matching: f0[v] = best in subtree with v unmatched, f1[v] = v matched to a child
    f0, f1, pick = {}, {}, {}
    for v in names:                             # topological: children come first in gate order
        ch = children[v]
        base = sum(max(f0[c], f1[c]) for c in ch)
        f0[v] = base; f1[v] = -1; pick[v] = None
        for c in ch:                            # match v with child c (c must be unmatched below)
            val = base - max(f0[c], f1[c]) + f0[c] + 1
            if val > f1[v]: f1[v], pick[v] = val, c
    # reconstruct from roots (nodes without a forest parent)
    matched = {}
    def walk(v, may_match):
        use_match = may_match and f1[v] > f0[v]
        if use_match:
            c = pick[v]; matched[c] = v
            for cc in children[v]:
                walk(cc, cc != c)
        else:
            for cc in children[v]: walk(cc, True)
    for v in names:
        if v not in parent: walk(v, True)
    # build merged circuit
    ops = {g: list(o) for g, o in gates}
    for c, p in matched.items():
        ops[p] = [o for o in ops[p] if o != c] + ops[c]
    merged = [(g, ops[g]) for g in names if g not in matched]
    return merged

def check(gates, outs, mat):
    rows = [[int(x) for x in l.split()] for l in open(mat) if l.strip()]
    n = len(rows[0])
    v = {f"x{c}": 1 << c for c in range(n)}; d = {f"x{c}": 0 for c in range(n)}
    for g, ops in gates:
        x = 0
        for o in ops: x ^= v[o]
        v[g] = x; d[g] = max(d[o] for o in ops) + 1
    ok = all(v[outs[j]] == sum(1 << c for c, b in enumerate(r) if b) for j, r in enumerate(rows))
    return ok, max(d[outs[j]] for j in outs)

if __name__ == "__main__":
    write = "--write" in sys.argv
    args = [a for a in sys.argv[1:] if a != "--write"]
    if len(args) < 2:
        sys.exit(__doc__)
    mat = args[0]
    for path in args[1:]:
        gates, outs = load_dbeam(path)
        merged = merge(gates, outs)
        ok, dep = check(merged, outs, mat)
        n2 = sum(len(o) == 2 for _, o in merged); n3 = sum(len(o) == 3 for _, o in merged)
        ge = ", ".join(f"{k} {a * n2 + b * n3:.2f} GE" for k, (a, b) in G2.items())
        print(f"{path}: {len(gates)} XOR2 -> {len(merged)} gates ({n2} XOR2 + {n3} XOR3), correct {ok}, depth {dep}; {ge}")
        if write:
            with open(path + ".xor3.txt", "w") as f:
                for g, o in merged: f.write(f"{g} = {' + '.join(o)}\n")
                for j in sorted(outs): f.write(f"y{j} = {outs[j]}\n")
