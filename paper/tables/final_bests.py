"""Best values and best circuits over all the runs of the paper (Table 3, Appendix A, results/).

Every recorded circuit of the runs at the minimum depth is checked against its matrix and the depth bound.
Its UMC area is that of the maximum merge of the gates of fan-out one (Circuit.merge_xor3), and its STM area
that of the exact merge with duplication of gates of fan-out two (stm_merge.py, CP-SAT). The runs form four
groups:

  main      the four stages: run/forward, run/bpbs, run/bbeam and run/lns of results/run.tar.xz, and, for
            [SS16] involutory 32x32 at its minimum depth 4, stages 2-4 in runs/depth4 of results/server_runs.tar.xz
  s4sfx     stage 4 started from the published circuits of SFX23 (runs/stage4fromsfx), marked by a double
            dagger in Table 3
  area      stage 4 with the UMC area as the objective (runs/area), marked by a section sign in Table 3
  analysis  the other runs at the minimum depth (runs/sfx1, runs/rank and runs/beamstats)

The best circuit of a matrix for a measure is the one with the smallest value; ties go to the group in the
order above, then to fewer XOR2 gates, then to the earliest stage.

usage: python3 paper/tables/final_bests.py RUN_PARENT SERVER_DIR OUT_DIR [RESULTS_DIR]
  RUN_PARENT  the folder that contains run/ (results/run.tar.xz extracted)
  SERVER_DIR  the folder that contains runs/ (results/server_runs.tar.xz extracted)
  writes OUT_DIR/all_circuits.jsonl (one record per circuit) and OUT_DIR/bests.json (per matrix: the best
  values of each group, the best circuit of each measure, and the best XOR2 count of each sfx1 run, for
  Table 5); with RESULTS_DIR, also writes there the best circuit of each matrix and measure (NIST .circ
  format), summary.json and README.md
requires: pip install ortools"""
import csv
import json
import os
import re
import statistics
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
from slpd.circuit import circ_text, read_any        # noqa: E402
from slpd.matrix import read_matrices               # noqa: E402
import stm_merge                                    # noqa: E402

GROUPS = ("main", "s4sfx", "area", "analysis")
RANK = {g: i for i, g in enumerate(GROUPS)}
SPLIT = re.compile(r"__|\.seed|\.s\d|\.d\d|\.bpbs|\.bbeam|\.lns|\.stage1|\.(?:sfxmd|sfx|all)\d")
DESC = {"main": "the four stages", "s4sfx": "stage 4 started from the published circuit of SFX23",
        "area": "stage 4 with the UMC area as the objective", "analysis": "a run of Secs. 5.3, 5.4 or 6"}


def sources(run_parent, server):
    run = os.path.join(run_parent, "run")
    for sub, stage in (("forward", 1), ("bpbs", 2), ("bbeam", 3), ("lns", 4)):
        yield "main", stage, os.path.join(run, sub), "paper run"
    runs = os.path.join(server, "runs")
    yield "main", None, os.path.join(runs, "depth4"), "server run"
    yield "s4sfx", 4, os.path.join(runs, "stage4fromsfx"), "server run"
    yield "area", 4, os.path.join(runs, "area"), "server run"
    for e in ("sfx1", "rank", "beamstats"):
        yield "analysis", 1, os.path.join(runs, e), "server run"


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    run_parent, server, out = sys.argv[1:4]
    results = sys.argv[4] if len(sys.argv) == 5 else None
    os.makedirs(out, exist_ok=True)
    rows = list(csv.DictReader(open(os.path.join(ROOT, "bench", "benchmark.tsv")), delimiter="\t"))
    tsv = {r["file"][:-4]: r for r in rows}
    mats, recs, solve = {}, [], []
    for group, stage0, d, where in sources(run_parent, server):
        if not os.path.isdir(d):
            sys.exit(f"missing folder {d}")
        for f in sorted(os.listdir(d)):
            m = re.match(r"(.+)\.(xor2|umc180|stm130)\.txt$", f)
            if not m:
                continue
            run, kind = m.group(1), m.group(2)
            name = SPLIT.split(run)[0]
            if name not in tsv:
                sys.exit(f"{f}: unknown matrix {name}")
            if name not in mats:
                M = read_matrices(os.path.join(ROOT, "bench", name + ".txt"))[0]
                mats[name] = (M, int(tsv[name]["min_depth"]))
            M, D = mats[name]
            path = os.path.join(d, f)
            C = read_any(path, M.n).prune()
            if any(len(g) != 2 for g in C.gates) or C.check(M, D):
                sys.exit(f"{path}: not a correct XOR2 circuit of depth <= {D}")
            U = C.merge_xor3()
            x2, x3 = U.counts()
            t0 = time.perf_counter()
            se, st = stm_merge.best_merge(C)
            solve.append(time.perf_counter() - t0)
            Mg = stm_merge.apply_merge(C, se, st)
            if Mg.check(M, D) or U.check(M, D):
                sys.exit(f"{path}: a merged circuit is wrong")
            s2, s3 = Mg.counts()
            stage = stage0 if stage0 else {"bpbs": 2, "bbeam": 3, "lns": 4}[re.search(r"\.(bpbs|bbeam|lns)\.", run).group(1)]
            recs.append(dict(group=group, exp=os.path.basename(d), stage=stage, matrix=name, run=run, kind=kind,
                             gates=len(C.gates), depth=C.depth(), x2=x2, x3=x3, umc=round(U.area("umc180"), 2),
                             stm=round(Mg.area("stm130"), 2), s2=s2, s3=s3, dup=len(st),
                             stm_nodup=round(3.33 * x2 + 4.66 * x3, 2),
                             path=os.path.relpath(path, run_parent if group == "main" and stage0 else server),
                             where=where))
    with open(os.path.join(out, "all_circuits.jsonl"), "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    print(f"{len(recs)} circuits checked; exact STM merges: median {statistics.median(solve):.4f} s, "
          f"max {max(solve):.2f} s (one thread), all optimal")

    keys = {"xor2": lambda r: (r["gates"], RANK[r["group"]], r["umc"], r["stage"], r["path"]),
            "umc180": lambda r: (r["umc"], RANK[r["group"]], r["gates"], r["stage"], r["path"]),
            "stm130": lambda r: (r["stm"], RANK[r["group"]], r["gates"], r["stage"], r["path"])}
    bests = {}
    for name in tsv:
        rs = [r for r in recs if r["matrix"] == name]
        if not rs:
            continue
        b = dict(groups={}, best={}, sfx1={}, nodup_stm=min(r["stm_nodup"] for r in rs))
        for g in GROUPS:
            gs = [r for r in rs if r["group"] == g]
            if gs:
                b["groups"][g] = dict(xor2=min(r["gates"] for r in gs),
                                      umc180=list(min((r["umc"], r["gates"]) for r in gs)),
                                      stm130=list(min((r["stm"], r["gates"]) for r in gs)))
        for k, key in keys.items():
            r = min(rs, key=key)
            v = r["gates"] if k == "xor2" else r["umc" if k == "umc180" else "stm"]
            reach = sorted({x["group"] for x in rs if abs((x["gates"] if k == "xor2" else x["umc" if k == "umc180" else "stm"]) - v) < 1e-9},
                           key=RANK.get)
            b["best"][k] = dict(value=v, gates=r["gates"], group=r["group"], groups_reaching=reach, exp=r["exp"],
                                stage=r["stage"], path=r["path"], where=r["where"])
        for r in rs:
            if r["exp"] == "sfx1" and r["kind"] == "xor2":
                b["sfx1"][re.search(r"seed(\d+)$", r["run"]).group(1)] = r["gates"]
        # the smallest UMC area among the circuits with the fewest XOR2 gates
        g0 = b["best"]["xor2"]["value"]
        b["umc_of_fewest"] = min(r["umc"] for r in rs if r["gates"] == g0)
        bests[name] = b
    json.dump(bests, open(os.path.join(out, "bests.json"), "w"), indent=1)
    if results:
        write_results(run_parent, server, results, tsv, mats, bests)


def write_results(run_parent, server, results, tsv, mats, bests):
    os.makedirs(results, exist_ok=True)
    summary, lines = {}, []
    lines += ["| Matrix | Size | D | XOR2 (SFX23) | XOR2 (ours) | UMC area (HILL) | UMC area (ours) | STM area (ours) |",
              "|" + "---|" * 8]
    for name, b in bests.items():
        M, D = mats[name]
        M.name = name
        res = {}
        for k in ("xor2", "umc180", "stm130"):
            e = b["best"][k]
            base = run_parent if e["where"] == "paper run" else server
            C = read_any(os.path.join(base, e["path"]), M.n).prune()
            extra = None
            if k == "xor2":
                K = C
            elif k == "umc180":
                K = C.merge_xor3()
                extra = [f"XOR2 circuit before merging: {len(C.gates)} gates, depth {C.depth()}"]
            else:
                se, st = stm_merge.best_merge(C)
                K = stm_merge.apply_merge(C, se, st)
                extra = [f"XOR2 circuit before merging: {len(C.gates)} gates, depth {C.depth()}"
                         + (f"; {len(st)} gate(s) of fan-out two duplicated before merging" if st else ""),
                         "Merge: exact for the STM unit areas, with duplication of gates of fan-out two "
                         "(paper/tables/stm_merge.py)"]
            if K.check(M, D):
                sys.exit(f"{name} {k}: the written circuit is wrong")
            src = f"{e['where']} ({DESC[e['group']]}, stage {e['stage']}), {e['path']}"
            fn = f"{name}.{k}.circ.txt"
            with open(os.path.join(results, fn), "w") as f:
                f.write(circ_text(K, M, f"{name}--{k}", src, extra))
            n2, n3 = K.counts()
            res[k] = dict(file=fn, group=e["group"], stage=e["stage"], source=e["path"], xor2_before_merge=len(C.gates),
                          xor2=n2, xor3=n3, depth=K.depth(), umc180=round(K.area("umc180"), 2),
                          stm130=round(K.area("stm130"), 2))
        r = tsv[name]
        summary[name] = dict(label=r["label"], size=r["size"], min_depth=D, sfx23_xor2=r["sfx23_xor2"] or None,
                             hill_umc180=r["hill_umc180"] or None, best=res)
        f2 = lambda v: f"{float(v):.2f}" if v not in (None, "", "-") else "-"
        def mark(k):          # the extra stage-4 runs that reach the value, if the four stages do not
            reach = b["best"][k]["groups_reaching"]
            return "" if "main" in reach else " " + "".join(c for g, c in (("s4sfx", "‡"), ("area", "§")) if g in reach)
        lines.append(f"| {r['label']} | {r['size']} | {D} | {r['sfx23_xor2'] or '-'} | {res['xor2']['xor2']}{mark('xor2')} | "
                     f"{f2(r['hill_umc180'])} | {f2(res['umc180']['umc180'])}{mark('umc180')} | "
                     f"{f2(res['stm130']['stm130'])}{mark('stm130')} |")
        print(f"{name}: XOR2 {res['xor2']['xor2']}, UMC {res['umc180']['umc180']}, STM {res['stm130']['stm130']}")
    with open(os.path.join(results, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    with open(os.path.join(results, "README.md"), "w") as f:
        f.write(open(os.path.join(HERE, "results_readme.md")).read().rstrip("\n") + "\n\n" + "\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
