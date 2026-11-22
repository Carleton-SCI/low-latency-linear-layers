"""The four-stage search of the paper, run with src/dbs.cpp on one matrix.

Stage 1  beam search over the ties of the depth-bounded forward search
Stage 2  backward splitting with random trajectories (our reimplementation of BPBS, Shi-Feng-Xu)
Stage 3  beam search over the backward splits
Stage 4  large neighborhood search from the best circuit of stages 1-3

Runs inside a stage are independent and run in parallel; the stages run in order.
"""
import concurrent.futures as cf
import json
import os
import re
import subprocess
import time

from .circuit import read_dbs_text, circ_text, read_text
from .matrix import write_rows

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DBS = os.path.join(ROOT, "bin", "dbs")


def ensure_dbs():
    """build bin/dbs with make if it is missing or older than src/dbs.cpp"""
    src = os.path.join(ROOT, "src", "dbs.cpp")
    if os.path.exists(DBS) and os.path.getmtime(DBS) >= os.path.getmtime(src):
        return DBS
    r = subprocess.run(["make", "-C", ROOT, "bin/dbs"], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(DBS):
        raise RuntimeError("building bin/dbs failed:\n" + r.stdout + r.stderr)
    return DBS


# ---------------------------------------------------------------------- presets
def preset_jobs(preset, n, m, stage):
    """argument lists for one stage. 'paper' is the exact configuration of the paper (its 16x16
    settings for n <= 16 inputs, its 32x32 settings otherwise); 'quick' and 'standard' are smaller
    budgets with the same structure."""
    small = n <= 16
    if preset == "paper":
        if stage == 1:
            grid = [("all", 2000, 1, range(1, 33)), ("bp", 2000, 1, range(1, 9)), ("sfx", 2000, 0, range(1, 9)),
                    ("all", 8000, 1, range(1, 5)), ("all", 2000, 2, range(1, 5))]
            return [(f"s1_{t}_w{w}_s{s}_seed{sd}", ["--tie", t, "--beam", str(w), "--slack", str(s), "--seed", str(sd)])
                    for t, w, s, seeds in grid for sd in seeds]
        seeds = range(1, 9)
        if stage == 2:
            T, MS = (2000, 16) if small else (300, 12)
            return [(f"s2_seed{sd}", ["--bpbs", str(T), "--maxsplits", str(MS), "--evals", "4", "--tie", "sfx", "--md",
                                      "--seed", str(sd)]) for sd in seeds]
        if stage == 3:
            R = 40 if small else 6
            return [(f"s3_seed{sd}", ["--bpbs", str(R), "--bbeam", "8", "--bopts", "16", "--evals", "2", "--maxsplits", "14",
                                      "--tie", "sfx", "--md", "--seed", str(sd)]) for sd in seeds]
        W, IT = (1000, 300) if small else (200, 150)
        return [(f"s4_seed{sd}", ["--iters", str(IT), "--fix", "4:28", "--beam", str(W), "--slack", "1", "--tie", "all",
                                  "--seed", str(sd)]) for sd in seeds]
    # quick and standard: smaller budgets, for trying the tool (README.md)
    q = preset == "quick"
    kmin = max(1, round(m / 8))
    kmax = max(kmin, round(7 * m / 8))
    if stage == 1:
        if q:
            W = 500 if small else 200
            grid = [("all", W, 1, range(1, 5 if small else 3)), ("bp", W, 1, range(1, 3 if small else 2)),
                    ("sfx", W, 0, range(1, 3 if small else 2))]
        else:
            grid = [("all", 2000, 1, range(1, 9 if small else 5)), ("bp", 2000, 1, range(1, 3)),
                    ("sfx", 2000, 0, range(1, 3)), ("all", 2000, 2, range(1, 3))]
        return [(f"s1_{t}_w{w}_s{s}_seed{sd}", ["--tie", t, "--beam", str(w), "--slack", str(s), "--seed", str(sd)])
                for t, w, s, seeds in grid for sd in seeds]
    seeds = range(1, 3) if q else range(1, 5)
    if stage == 2:
        T = (100 if small else 2) if q else (500 if small else 20)
        MS = 16 if small else 12
        return [(f"s2_seed{sd}", ["--bpbs", str(T), "--maxsplits", str(MS), "--evals", "4", "--tie", "sfx", "--md",
                                  "--seed", str(sd)]) for sd in seeds]
    if stage == 3:
        if small:
            R, B, K, MS = (4 if q else 10), 8, 16, 14
        else:
            R, B, K, MS = (1, 2, 4, 6) if q else (1, 4, 8, 10)
        return [(f"s3_seed{sd}", ["--bpbs", str(R), "--bbeam", str(B), "--bopts", str(K), "--evals", "2",
                                  "--maxsplits", str(MS), "--tie", "sfx", "--md", "--seed", str(sd)]) for sd in seeds]
    W, IT = ((300, 60) if small else (50, 10)) if q else ((1000, 150) if small else (200, 75))
    return [(f"s4_seed{sd}", ["--iters", str(IT), "--fix", f"{kmin}:{kmax}", "--beam", str(W), "--slack", "1",
                              "--tie", "all", "--seed", str(sd)]) for sd in seeds]


# ---------------------------------------------------------------------- running dbs
BEST_RE = re.compile(r"BEST-(\w+): (\d+) XOR2 gates -> (\d+) XOR2 \+ (\d+) XOR3, umc180 ([\d.]+) GE, stm130 ([\d.]+) GE")


def run_one(dbs, rows_file, D, prefix, args, timeout):
    t0 = time.time()
    cmd = [dbs, rows_file, str(D)] + args + ["--save", prefix]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        out, err, status = r.stdout, r.stderr, ("ok" if r.returncode == 0 else f"exit {r.returncode}")
    except subprocess.TimeoutExpired as e:
        out, err, status = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""), "", "timeout"
    with open(prefix + ".log", "w") as f:
        f.write("# " + " ".join(cmd) + "\n" + out + ("\n# stderr:\n" + err if err else "") +
                f"\n# {time.time() - t0:.1f} s, {status}\n")
    best = {}
    for line in out.splitlines():
        b = BEST_RE.match(line)
        if b:
            best[b.group(1)] = dict(gates=int(b.group(2)), n2=int(b.group(3)), n3=int(b.group(4)),
                                   umc180=float(b.group(5)), stm130=float(b.group(6)), file=f"{prefix}.{b.group(1)}.txt")
    return dict(prefix=prefix, seconds=round(time.time() - t0, 1), status=status, best=best)


def key_of(measure, b):
    if measure == "xor2":
        return (b["gates"], b["umc180"])
    return (b[measure], b["gates"])


def run_search(M, out_dir, D=None, preset="quick", stages=(1, 2, 3, 4), jobs=None, timeout=None, log=print):
    """Run the stages on matrix M; write the best circuits (NIST .circ format) and a summary to
    out_dir. Returns the summary dict."""
    dbs = ensure_dbs()
    s = M.min_depth()
    D = s if D is None else D
    if D < s:
        raise ValueError(f"depth bound {D} is below the minimum depth {s}")
    stages = sorted(set(stages))
    if D > s and any(st > 1 for st in stages):
        log(f"note: stages 2-4 assume the depth bound equals the minimum depth ({s}); with D = {D} "
            "only stage 1 runs")
        stages = [1]
    work = os.path.join(out_dir, "work", M.name)
    os.makedirs(work, exist_ok=True)
    rows_file = os.path.join(work, M.name + ".rows.txt")
    write_rows(M, rows_file)
    jobs = jobs or os.cpu_count() or 1
    results, per_stage = [], {}
    t_start = time.time()
    for st in stages:
        todo = preset_jobs(preset, M.n, M.m, st)
        if st == 4:
            done = [r for r in results if "xor2" in r["best"]]
            if not done:
                log("stage 4 skipped: no circuit from stages 1-3")
                continue
            inc = min(done, key=lambda r: key_of("xor2", r["best"]["xor2"]))["best"]["xor2"]["file"]
            todo = [(lab, ["--lns", inc] + args) for lab, args in todo]
        log(f"stage {st}: {len(todo)} runs on {min(jobs, len(todo))} workers")
        t0 = time.time()
        with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = [ex.submit(run_one, dbs, rows_file, D, os.path.join(work, lab), args, timeout) for lab, args in todo]
            stage_res = [f.result() for f in futs]
        for r in stage_res:
            r["stage"] = st
        results += stage_res
        bad = [r for r in stage_res if r["status"] != "ok"]
        got = [r["best"]["xor2"]["gates"] for r in stage_res if "xor2" in r["best"]]
        per_stage[st] = dict(runs=len(stage_res), failed=len(bad), seconds=round(time.time() - t0, 1),
                             best_xor2=min(got) if got else None)
        log(f"  done in {time.time() - t0:.0f} s; best XOR2 count {per_stage[st]['best_xor2']}"
            + (f"; {len(bad)} runs failed or timed out" if bad else ""))
    summary = dict(matrix=M.name, rows=M.m, columns=M.n, min_depth=s, depth_bound=D, preset=preset,
                   stages=per_stage, seconds=round(time.time() - t_start, 1), best={})
    done = [r for r in results if r["best"]]
    if not done:
        raise RuntimeError("no run produced a circuit; see the logs in " + work)
    for measure in ("xor2", "umc180", "stm130"):
        cand = [r for r in done if measure in r["best"]]
        r = min(cand, key=lambda r: key_of(measure, r["best"][measure]))
        b = r["best"][measure]
        C = read_dbs_text(read_text(b["file"]), b["file"], M.n).prune()
        problems = C.check(M, D)
        if problems:
            raise RuntimeError(f"{b['file']} does not verify: {problems}")
        target = C if measure == "xor2" else C.merge_xor3()
        problems = target.check(M, D)
        if problems:
            raise RuntimeError(f"merged {b['file']} does not verify: {problems}")
        n2, n3 = target.counts()
        path = os.path.join(out_dir, f"{M.name}.{measure}.circ.txt")
        src = f"slpd {preset} preset, stage {r['stage']}, run {os.path.basename(r['prefix'])}"
        extra = [f"XOR2 circuit before merging: {len(C.gates)} gates, depth {C.depth()}"] if measure != "xor2" else None
        with open(path, "w") as f:
            f.write(circ_text(target, M, f"{M.name}--{measure}", src, extra))
        summary["best"][measure] = dict(file=os.path.basename(path), stage=r["stage"], run=os.path.basename(r["prefix"]),
                                        xor2_gates_before_merge=len(C.gates), xor2=n2, xor3=n3, depth=target.depth(),
                                        umc180=round(target.area("umc180"), 2), stm130=round(target.area("stm130"), 2),
                                        verified=True)
    with open(os.path.join(out_dir, f"{M.name}.summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    return summary
