"""Run dbs over a benchmark list with a grid of search configurations, in parallel, and summarize
the best result per matrix in three measures (XOR2 count, XOR2/XOR3 area with UMC 180 nm and
STM 130 nm cell areas) next to the published references. Resumable: finished runs are skipped.

usage: python3 bench_run.py BENCH_DIR OUT_DIR [--workers N] [--size 16|32|all] [--grid quick|full]
                            [--only name1,name2] [--dbs PATH]
BENCH_DIR holds the matrices and benchmark.tsv (file, label, size, min_depth, sfx23_xor2,
hill_umc180, hill_depth, hill_input_xor2)."""
import os, sys, csv, re, json, time, argparse, subprocess, concurrent.futures as cf

HERE = os.path.dirname(os.path.abspath(__file__))
GRIDS = {   # (tie, width, slack, seeds[, extra options of dbs and the name suffix of their runs])
    "quick": [("all", 2000, 1, range(1, 5)), ("bp", 2000, 1, range(1, 3)), ("sfx", 2000, 0, range(1, 3))],
    "full": [("all", 2000, 1, range(1, 33)), ("bp", 2000, 1, range(1, 9)), ("sfx", 2000, 0, range(1, 9)),
             ("all", 8000, 1, range(1, 5)), ("all", 2000, 2, range(1, 5))],
    # rank the options of the beam search by gates so far + Sigma (seeds 1-8 of all/2000/1)
    "rankgates": [("all", 2000, 1, range(1, 9), ["--rank", "gates"], "_rgates")],
}

def job(args):
    dbs, mat, depth, tie, width, slack, seed, out, extra = args
    if os.path.exists(out + ".log") and "BEST-stm130" in open(out + ".log").read(): return out, 0.0
    t0 = time.time()
    r = subprocess.run([dbs, mat, str(depth), "--tie", tie, "--beam", str(width), "--slack", str(slack),
                        "--seed", str(seed), "--save", out] + list(extra), capture_output=True, text=True)
    open(out + ".log", "w").write(r.stdout + f"# {time.time() - t0:.1f} s\n")
    return out, time.time() - t0

def parse(log):
    best = {}
    for line in open(log):
        m = re.match(r"BEST-(\w+): (\d+) XOR2 gates -> (\d+) XOR2 \+ (\d+) XOR3, umc180 ([\d.]+) GE, stm130 ([\d.]+) GE", line)
        if m: best[m.group(1)] = dict(gates=int(m.group(2)), n2=int(m.group(3)), n3=int(m.group(4)),
                                      umc=float(m.group(5)), stm=float(m.group(6)))
    return best

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("bench"); ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=os.cpu_count()); ap.add_argument("--size", default="all")
    ap.add_argument("--grid", default="quick", choices=GRIDS); ap.add_argument("--only", default="")
    ap.add_argument("--dbs", default=os.path.join(HERE, "..", "bin", "dbs"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rows = list(csv.DictReader(open(os.path.join(a.bench, "benchmark.tsv")), delimiter="\t"))
    if a.size != "all": rows = [r for r in rows if r["size"].startswith(a.size + "x")]
    if a.only: rows = [r for r in rows if r["file"][:-4] in a.only.split(",")]
    jobs = []
    for r in rows:
        name = r["file"][:-4]
        for g in GRIDS[a.grid]:
            tie, width, slack, seeds = g[:4]
            extra, suffix = (g[4], g[5]) if len(g) > 4 else ([], "")
            for sd in seeds:
                jobs.append((a.dbs, os.path.join(a.bench, r["file"]), int(r["min_depth"]), tie, width, slack, sd,
                             os.path.join(a.out, f"{name}__{tie}_w{width}_s{slack}{suffix}_seed{sd}"), extra))
    # larger matrices first
    jobs.sort(key=lambda j: -os.path.getsize(j[1]))
    print(f"{len(rows)} matrices, {len(jobs)} runs, {a.workers} workers", flush=True)
    t0 = time.time(); done = 0
    with cf.ProcessPoolExecutor(a.workers) as ex:
        for out, dt in ex.map(job, jobs):
            done += 1
            if done % 10 == 0 or done == len(jobs):
                print(f"  {done}/{len(jobs)} runs, {(time.time() - t0) / 60:.1f} min", flush=True)
    # summary
    summ = {}
    print(f"\n{'matrix':38s} {'SFX23':>5s} {'ours':>5s} | {'HILL UMC':>8s} {'ours UMC':>9s} | {'ours STM':>8s}  (XOR2 at the min depth; areas after XOR3 merging)")
    for r in rows:
        name = r["file"][:-4]
        bx, bu, bs = None, None, None
        for fn in os.listdir(a.out):
            if not (fn.startswith(name + "__") and fn.endswith(".log")): continue
            b = parse(os.path.join(a.out, fn))
            if not b: continue
            cfg = fn[len(name) + 2:-4]
            if bx is None or b["xor2"]["gates"] < bx[0]: bx = (b["xor2"]["gates"], cfg)
            if bu is None or b["umc180"]["umc"] < bu[0]: bu = (b["umc180"]["umc"], cfg)
            if bs is None or b["stm130"]["stm"] < bs[0]: bs = (b["stm130"]["stm"], cfg)
        if bx is None: continue
        sx = int(r["sfx23_xor2"]) if r["sfx23_xor2"] else None
        hu = float(r["hill_umc180"]) if r["hill_umc180"] else None
        summ[name] = dict(label=r["label"], sfx23=sx, hill_umc=hu, hill_depth=r["hill_depth"] or None, xor2=bx, umc=bu, stm=bs)
        print(f"{r['label']:38s} {sx if sx is not None else '-':>5} {bx[0]:5d} | "
              f"{format(hu, '8.2f') if hu is not None else '-':>8} {bu[0]:9.2f} | {bs[0]:8.2f}"
              + ("   (HILL at depth %s)" % r["hill_depth"] if r["hill_depth"] not in ("", r["min_depth"]) else ""))
    json.dump(summ, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
