"""Summarize a benchmark run (run_benchmark.sh step 5): for every matrix, the best circuit by XOR2
count and by XOR2/XOR3 area (UMC 180 nm and STM 130 nm cell areas, HILL ToSC 2026 Table 1),
re-verified independently from the saved netlists (outputs against the matrix, depth bound, area
recomputed with xor3_merge.py), next to the published references in bench/benchmark.tsv.

usage: python3 summarize.py BENCH_DIR RUN_DIR [--hill HILL_REPO] [--sfx SFX_DIR [--export OUT_DIR]]

  --hill HILL_REPO  also compare every best circuit with HILL's published circuit for the same matrix
                    (results/<name>/*_result.txt of github.com/13est-Wishes/Linear_Layers_by_YongchaoLi):
                    both are evaluated on all unit inputs, so this check does not use bench/ at all
  --sfx SFX_DIR     find each bench matrix among Shi-Feng-Xu's files: the result circuits of their
                    repository (our_results/*.txt; every circuit is re-evaluated, its gate count and depth
                    recomputed, and the matrix it computes rebuilt) and any plain 0/1 matrix files (such
                    as testmatrix.txt). Reports whether each bench matrix appears identical, with rows and
                    columns permuted, or transposed (needs networkx), with the best SFX circuit for it; a
                    matrix not found is compared entrywise with the nearest SFX matrix. With --export
                    OUT_DIR, every SFX matrix is written in bench format and every SFX circuit in the
                    netlist format of dbs (usable as a --lns starting point)"""
import os, sys, re, csv, glob, json, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xor3_merge import load_dbeam, merge, check

def best_lines(log):
    out = {}
    for line in open(log):
        m = re.match(r"BEST-(\w+): (\d+) XOR2 gates -> (\d+) XOR2 \+ (\d+) XOR3", line)
        if m: out[m.group(1)] = (int(m.group(2)), int(m.group(3)), int(m.group(4)))
    return out

def area(n2, n3, lib): return 3.00 * n2 + 4.67 * n3 if lib == "umc180" else 3.33 * n2 + 4.66 * n3

def num(x, f): return f(x) if x not in (None, "") else None          # empty benchmark.tsv cell -> None
def fmt(x, spec): return format(x, spec) if x is not None else format("-", ">" + spec.split(".")[0].rstrip("df"))

# ---- HILL cross-check: evaluate both circuits on every unit input ------------------------------
def hill_outputs(path):
    """HILL result format 't[g] = t[a] + t[b] [+ t[c]]  y[k] ...'; returns {k: value over inputs t[i]}"""
    gates, outs = [], {}
    for line in open(path):
        m = re.match(r"\s*t\[(\d+)\]\s*=\s*(.+?)\s*((?:y\[\d+\]\s*)*)$", line)
        if not m or line.startswith("Optimization"): continue
        gates.append((int(m.group(1)), [int(x) for x in re.findall(r"t\[(\d+)\]", m.group(2))]))
        for k in re.findall(r"y\[(\d+)\]", m.group(3)): outs[int(k)] = int(m.group(1))
    assigned = {g for g, _ in gates}
    n = max({o for _, ops in gates for o in ops} - assigned) + 1
    v = {i: 1 << i for i in range(n)}
    for g, ops in gates:
        x = 0
        for o in ops: x ^= v[o]
        v[g] = x
    return {k: v[g] for k, g in outs.items()}, n

def our_outputs(path, n):
    gates, outs = load_dbeam(path)
    v = {f"x{c}": 1 << c for c in range(n)}
    for g, ops in gates:
        x = 0
        for o in ops: x ^= v[o]
        v[g] = x
    return {j: v[o] for j, o in outs.items()}

# ---- Shi-Feng-Xu matrix files: flexible reader and permutation equivalence ----------------------
def read_blocks(path):
    """all maximal runs of lines that are rows of 0/1 (with optional separators) of equal width"""
    blocks, cur = [], []
    try:
        lines = open(path, errors="ignore").read().splitlines()
    except OSError:
        return []
    for line in lines + [""]:
        s = line.strip()
        bits = re.sub(r"[\s,;\[\](){}]", "", s)
        if len(bits) >= 4 and set(bits) <= {"0", "1"} and (not cur or len(bits) == len(cur[0])):
            cur.append([int(b) for b in bits]); continue
        if cur: blocks.append(cur)
        cur = []
        if len(bits) >= 4 and set(bits) <= {"0", "1"}: cur = [[int(b) for b in bits]]
    out = []
    for b in blocks:                       # split stacked square matrices
        w = len(b[0])
        if len(b) > w and len(b) % w == 0: out += [b[i:i + w] for i in range(0, len(b), w)]
        else: out.append(b)
    return out

def bigraph(M):
    import networkx as nx
    G = nx.Graph()
    for i in range(len(M)): G.add_node(("r", i), side=0)
    for j in range(len(M[0])): G.add_node(("c", j), side=1)
    for i, row in enumerate(M):
        for j, b in enumerate(row):
            if b: G.add_edge(("r", i), ("c", j))
    return G

def same_up_to_perm(A, B):
    import networkx as nx
    if (len(A), len(A[0])) != (len(B), len(B[0])): return False
    if sorted(map(sum, A)) != sorted(map(sum, B)) or sorted(map(sum, zip(*A))) != sorted(map(sum, zip(*B))): return False
    GA, GB = bigraph(A), bigraph(B)
    if nx.weisfeiler_lehman_graph_hash(GA, node_attr="side") != nx.weisfeiler_lehman_graph_hash(GB, node_attr="side"):
        return False
    return nx.is_isomorphic(GA, GB, node_match=lambda a, b: a["side"] == b["side"])

def relation(A, B):
    """how matrix B (an SFX matrix) relates to A (a bench matrix), or None if unrelated"""
    At = [list(c) for c in zip(*A)]
    for X, tag in ((A, ""), (At, "TRANSPOSED")):
        if (len(X), len(X[0])) != (len(B), len(B[0])): continue
        if X == B: return tag or "identical"
        if same_up_to_perm(X, B): return (tag + ", " if tag else "") + "rows/columns permuted"
    return None

SFX_GATE = re.compile(r"x(\d+) x(\d+) x(\d+)((?: y\d+)*)(?: depth (\d+))?")
SFX_HEAD = re.compile(r"use\s*(\d+) xors")

def read_sfx_circuits(path):
    """implementations in a Shi-Feng-Xu result file: header '... use N xors, split size s', a line with
    the n unit vectors 1 2 4 ... (inputs x0..x(n-1)) followed by the node values, the gates
    'xK xA xB [yJ ...] [depth d]' (the first one sometimes on the same line as the numbers), and a
    trailer 'xorsN'; a file may hold several implementations. Returns one dict per implementation with
    the matrix it computes (row j = output yj, column c = input xc, the bench convention), the gate
    count and depth recomputed from the gates, the netlist in dbs format, and any problems found"""
    try:
        text = open(path, errors="ignore").read().replace("\r", "")
    except OSError:
        return []
    starts = [m.start() for m in SFX_HEAD.finditer(text)]
    impls = []
    for k, s in enumerate(starts):
        chunk = text[s:starts[k + 1] if k + 1 < len(starts) else len(text)]
        lines = chunk.split("\n")
        nums = re.findall(r"\d+", re.split(r"x\d", lines[1])[0]) if len(lines) > 1 else []
        n = 0
        while n < len(nums) and int(nums[n]) == 1 << n: n += 1
        if n == 0: continue
        name = lambda i: f"x{i}" if i < n else f"g{i - n}"
        v = {i: 1 << i for i in range(n)}; d = {i: 0 for i in range(n)}
        outs, problems, net = {}, [], []
        for m in SFX_GATE.finditer(chunk):
            g, a, b = int(m.group(1)), int(m.group(2)), int(m.group(3))
            if a not in v or b not in v or g in v:
                problems.append(f"bad gate '{m.group(0)}'"); continue
            v[g] = v[a] ^ v[b]; d[g] = max(d[a], d[b]) + 1
            net.append(f"{name(g)} = {name(a)} + {name(b)}")
            if m.group(5) and int(m.group(5)) != d[g]: problems.append(f"'{m.group(0)}' has depth {d[g]}")
            for y in re.findall(r"y(\d+)", m.group(4)): outs.setdefault(int(y), g)
        if not outs: continue
        m_out = max(outs) + 1
        missing = [j for j in range(m_out) if j not in outs]
        if missing: problems.append(f"outputs {missing} not computed by a gate")
        claimed = int(SFX_HEAD.match(chunk).group(1))
        trailer = re.search(r"xors(\d+)", chunk)
        if claimed != len(net) or (trailer and int(trailer.group(1)) != len(net)):
            problems.append(f"{len(net)} gates, header says {claimed}" + (f", trailer {trailer.group(1)}" if trailer else ""))
        M = None if missing else [[(v[outs[j]] >> c) & 1 for c in range(n)] for j in range(m_out)]
        impls.append(dict(idx=k, n=n, xor2=len(net), depth=max(d[g] for g in outs.values()), M=M,
                          problems=problems, net=net + [f"y{j} = {name(g)}" for j, g in sorted(outs.items())]))
    return impls

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("bench"); ap.add_argument("run")
    ap.add_argument("--hill"); ap.add_argument("--sfx"); ap.add_argument("--export")
    a = ap.parse_args()
    bench, run = a.bench, a.run
    rows = list(csv.DictReader(open(os.path.join(bench, "benchmark.tsv")), delimiter="\t"))
    summ = {}
    hdr = f"{'matrix':36s} {'depth':>5s} | {'SFX23':>5s} {'ours':>4s} {'method':8s} | {'HILL UMC':>8s} {'ours UMC':>8s} | {'ours STM':>8s}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        name, mat, D = r["file"][:-4], os.path.join(bench, r["file"]), int(r["min_depth"])
        logs = [l for pat in (f"{run}/*/{name}__*.log", f"{run}/*/{name}.seed*.log", f"{run}/*/{name}.log")
                for l in glob.glob(pat)]
        best = {}
        for log in logs:
            method = log.split("/")[-2]
            for key, (g, n2, n3) in best_lines(log).items():
                score = g if key == "xor2" else area(n2, n3, key)
                if key not in best or score < best[key][0] - 1e-9:
                    best[key] = (score, log[:-4] + f".{key}.txt", method)
        if not best: continue
        res = {}
        for key, (score, path, method) in best.items():
            gates, outs = load_dbeam(path)
            ok2, dep2 = check(gates, outs, mat)
            merged = merge(gates, outs); ok3, dep3 = check(merged, outs, mat)
            n2 = sum(len(o) == 2 for _, o in merged); n3 = sum(len(o) == 3 for _, o in merged)
            ok = ok2 and ok3 and dep2 <= D and dep3 <= D
            res[key] = dict(xor2=len(gates), n2=n2, n3=n3, umc180=round(area(n2, n3, "umc180"), 2),
                            stm130=round(area(n2, n3, "stm130"), 2), depth=dep2, verified=ok, method=method, file=path)
            if not ok: print(f"  !! {path} failed verification")
            if a.hill:
                hp = glob.glob(os.path.join(a.hill, "results", name, "*_result.txt"))
                if len(hp) != 1:
                    res[key]["same_map_as_hill"] = None
                else:
                    ho, n = hill_outputs(hp[0])
                    res[key]["same_map_as_hill"] = ho == our_outputs(path, n)
                    if not res[key]["same_map_as_hill"]: print(f"  !! {path} computes a different map than {hp[0]}")
        sx, hu, hd = num(r["sfx23_xor2"], int), num(r["hill_umc180"], float), num(r["hill_depth"], int)
        summ[name] = dict(label=r["label"], min_depth=D, sfx23_xor2=sx, hill_umc180=hu, hill_depth=hd, best=res)
        x, u, s = res["xor2"], res["umc180"], res["stm130"]
        print(f"{r['label'][:36]:36s} {D:5d} | {fmt(sx, '5d')} {x['xor2']:4d} {x['method']:8s} | "
              f"{fmt(hu, '8.2f')} {u['umc180']:8.2f} | {s['stm130']:8.2f}"
              + ("" if hd in (None, D) else f"  (HILL at depth {hd})"))
    json.dump(summ, open(os.path.join(run, "summary.json"), "w"), indent=1)    # written again at the end
    if a.hill:
        checked = [b.get("same_map_as_hill") for m in summ.values() for b in m["best"].values()]
        print(f"\nHILL cross-check: {sum(c is True for c in checked)} circuits compute the same map as HILL's "
              f"published circuit, {sum(c is False for c in checked)} differ, {sum(c is None for c in checked)} not found")
    if a.sfx:
        try:
            import networkx  # noqa: F401
        except ImportError:
            print("\n--sfx needs networkx: python3 -m pip install --user networkx")
            a.sfx = None
    if a.sfx:
        files = sorted(p for p in glob.glob(os.path.join(a.sfx, "**", "*"), recursive=True)
                       if os.path.isfile(p) and os.path.getsize(p) < 50_000_000 and "/.git/" not in p)
        cands = []     # (source, matrix, circuit or None)
        for p in files:
            rel = os.path.relpath(p, a.sfx)
            for I in read_sfx_circuits(p):
                for msg in I["problems"]: print(f"  !! {rel} #{I['idx']}: {msg}")
                if I["M"]: cands.append((f"{rel} #{I['idx']}", I["M"], I))
            for i, B in enumerate(read_blocks(p)): cands.append((f"{rel} [0/1 matrix {i}]", B, None))
        nc = sum(c[2] is not None for c in cands)
        print(f"\nShi-Feng-Xu files: {len(files)} files, {nc} circuits, {len(cands) - nc} plain 0/1 matrices")
        if a.export:
            os.makedirs(a.export, exist_ok=True)
            for src, B, I in cands:
                fn = re.sub(r"[^\w.-]+", "_", src).strip("_")
                open(os.path.join(a.export, fn + ".matrix.txt"), "w").write("".join(" ".join(map(str, row)) + "\n" for row in B))
                if I: open(os.path.join(a.export, fn + ".dbs.txt"), "w").write("\n".join(I["net"]) + "\n")
        print(f"{'matrix':36s} SFX circuit for it (XOR2, depth; table value) / where it appears")
        for r in rows:
            name, D = r["file"][:-4], int(r["min_depth"])
            M = [[int(x) for x in l.split()] for l in open(os.path.join(bench, r["file"])) if l.strip()]
            hits = [(src, rel_, I) for src, B, I in cands for rel_ in [relation(M, B)] if rel_]
            circ = sorted((I["xor2"], I["depth"], rel_ != "identical", src, rel_) for src, rel_, I in hits
                          if I and "TRANSPOSED" not in rel_)
            info = dict(found=[(src, rel_) for src, rel_, _ in hits])
            if circ:
                g, dep, _, src, rel_ = circ[0]
                info.update(sfx_xor2=g, sfx_depth=dep, sfx_circuit=src, sfx_relation=rel_)
                tab = num(r.get("sfx23_xor2"), int)
                line = f"{g} XOR2, depth {dep} (table {fmt(tab, 'd')}); {src}, {rel_}"
                ident = [c for c in circ if c[4] == "identical"]
                if rel_ != "identical" and ident: line += f"; identical: {ident[0][3]}, {ident[0][0]} XOR2"
                if tab is not None and g < tab: line += "  <-- below the table value"
                if tab is not None and g > tab: line += "  <-- table value not reproduced"
                if dep > D: line += "  <-- depth above the minimum"
            elif hits:
                line = "; ".join(f"{src} ({rel_})" for src, rel_, _ in hits[:3])
            else:
                near = sorted((sum(x != y for ra, rb in zip(M, B) for x, y in zip(ra, rb)), src, B) for src, B, _ in cands
                              if (len(B), len(B[0])) == (len(M), len(M[0])))
                line = "NOT FOUND"
                if near and near[0][0] <= 16:
                    k, src, B = near[0]
                    diff = [(i, j) for i in range(len(M)) for j in range(len(M[0])) if M[i][j] != B[i][j]]
                    line += f"; nearest {src} differs in {k} entries (row, column): {diff}"
                    info["nearest"] = dict(source=src, entries=diff)
            print(f"{name[:36]:36s} {line}")
            summ.setdefault(name, {})["sfx"] = info
    json.dump(summ, open(os.path.join(run, "summary.json"), "w"), indent=1)
