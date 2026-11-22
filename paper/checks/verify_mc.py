"""Independent check of a dbp/dbeam circuit: every output equals its matrix row, every output has
XOR depth <= D, and the gate count. Inputs are named x<c> with c the 0-based matrix column.

usage: python3 verify_mc.py MATRIX CIRCUIT D      (bare rows; native format of dbs)"""
import re, sys
if len(sys.argv) != 4:
    sys.exit(__doc__)
mat, circ, D = sys.argv[1], sys.argv[2], int(sys.argv[3])
rows = [[int(x) for x in l.split()] for l in open(mat) if l.strip()]
val, dep, ngates = {}, {}, 0
for c in range(len(rows[0])): val[f"x{c}"] = 1 << c; dep[f"x{c}"] = 0
outs = {}
for l in open(circ):
    m = re.match(r"(g\d+) = (\w+) \+ (\w+)\s+\[(\d+)\]", l)
    if m:
        g, a, b = m.group(1), m.group(2), m.group(3)
        val[g] = val[a] ^ val[b]; dep[g] = max(dep[a], dep[b]) + 1; ngates += 1
        assert dep[g] == int(m.group(4)), f"depth mismatch at {g}"
    m = re.match(r"y(\d+) = (\w+)", l)
    if m: outs[int(m.group(1))] = m.group(2)
ok = all(val[outs[j]] == sum(1 << c for c, x in enumerate(r) if x) for j, r in enumerate(rows))
print(f"{circ}: {ngates} XOR gates, outputs correct: {ok}, max output depth: {max(dep[outs[j]] for j in outs)} (limit {D}), "
      f"outputs: {len(outs)}/{len(rows)}")
