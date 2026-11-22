"""Tie exploration against random restarts, at equal CPU time.

The question: does the beam search over the ties (stage 1 of the paper, which extends
the tie exploration of Improved-BP, TCHES 2018) find smaller circuits than independent greedy runs
with random tie-breaking (the approach of SFX23) when both get the same CPU time? Everything else is
the same: the same program (src/dbs.cpp), the same exact depth-aware distances, the same tie rule,
the same depth bound (the minimum depth), and forward search only.

For every matrix, configuration ("pair") and seed, one task runs:
  1. the beam search:  dbs MATRIX D <rule> --beam W --slack s --seed k
     and measures its CPU time T (user + system time of the process);
  2. random restarts:   dbs MATRIX D <rule> --runs R --seed ...
     in chunks, until their total CPU time reaches T (each chunk is sized from the measured time per
     run, to use 85% of the remaining budget, so the restarts get T up to the error of one estimate).
     Every restart has its own seed.

Pairs (the rule is the same on both sides, except for xall1):
  all0    --tie all, slack 0: the tie sets of the BFS and of the greedy runs are the same. This is
          the clean test of tie exploration against random tie-breaking.
  all1    --tie all, slack 1: the main stage-1 configuration of the paper (seeds 1-4 of this pair
          reproduce the paper's stage-1 runs of the same seeds).
  sfxmd0  --tie sfx --md, slack 0: the forward rule of SFX23 (minimum norm, MD filter) on both sides.
  sfx0, bp1: also available (sfx0 separates the norm ordering from the MD filter).
  sfxmd1  --tie sfx --md, slack 1: the forward rule of SFX23 on both sides, with the slack of stage 1.
  xall1   beam search with --tie all and slack 1 (the main stage-1 configuration) against random restarts
          with the forward rule of SFX23 (--tie sfx --md): the comparison of the best beam configuration
          with the restarts of SFX23.
For sfx and bp, the BFS keeps the options of minimum sum of distances ordered by the norm (the first W),
while a greedy run picks at random among the options of minimum sum and best norm; for all, both sides
use the same set of options. Backward splitting is not part of this comparison.

usage:
  python3 paper/ties_vs_restarts.py run [-o OUT] [--dbs PATH] [--bench DIR] [-j N] [--pin] [--seeds 1-4]
                                        [--pairs all0,all1,sfxmd0] [--width 2000] [--only A,B] [--skip A,B]
                                        [--with-hill-ls16] [--timeout S] [--status-every MIN]
                                        [--tex FILE] [--paper-run RUN_DIR] [--dry-run]
  python3 paper/ties_vs_restarts.py summary OUT [--tex FILE] [--paper-run RUN_DIR] [--keep-duplicates]

'run' runs in the foreground (use screen or tmux for long runs). It prints a line for every finished
task and a status line every 10 minutes, keeps OUT/progress.txt with the running tasks, and prints
the summary at the end. Ctrl-C stops it cleanly; it is resumable: finished tasks are recorded in
OUT/results.jsonl and skipped when it is started again. Only one 'run' at a time can use OUT (a
second one stops with a message). 'summary' also works while 'run' is going, on the tasks finished
so far; it uses one record per (matrix, pair, seed, width), the first one, and warns about
duplicates (from older versions, where two runs could share OUT), unless --keep-duplicates. Use as many
jobs as physical cores (the default), so that the two sides of a task run under the same conditions;
with hyper-threads, the CPU time per unit of work depends on the load. --pin also runs each task (its
beam search and all its restarts) on one logical CPU of its own physical core, so that no two tasks
share a core.
Requires only Python 3.6+ and the compiled search program."""
import argparse
import concurrent.futures as cf
import queue
import csv
import json
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PAIRS = {   # name: (rule of the beam search, its slack, rule of the random restarts)
    "all0": (["--tie", "all"], 0, ["--tie", "all"]),
    "all1": (["--tie", "all"], 1, ["--tie", "all"]),
    "sfxmd0": (["--tie", "sfx", "--md"], 0, ["--tie", "sfx", "--md"]),
    "sfx0": (["--tie", "sfx"], 0, ["--tie", "sfx"]),
    "sfxmd1": (["--tie", "sfx", "--md"], 1, ["--tie", "sfx", "--md"]),
    "bp1": (["--tie", "bp"], 1, ["--tie", "bp"]),
    "xall1": (["--tie", "all"], 1, ["--tie", "sfx", "--md"]),
}
HILL_LS16 = "FSE_LiuSim16_32_Circulant"           # the [LS16] 32x32 matrix of the HILL circuits (Sec. 4)
BEST = re.compile(r"BEST-(\w+): (\d+) XOR2 gates -> (\d+) XOR2 \+ (\d+) XOR3, umc180 ([\d.]+) GE, stm130 ([\d.]+) GE")
MEASURES = ("xor2", "umc180", "stm130")
LOCK = threading.Lock()
STOP = threading.Event()            # set on Ctrl-C: no new process is started
DONE = threading.Event()            # set at the end: the status thread stops
RUNNING = {}                        # pid -> Popen of the running search processes
STATE = {}                          # (matrix, pair, seed) -> [phase, start time, detail] of the running tasks
TL = threading.local()              # TL.cpu: the logical CPU of the task of this thread (--pin), or None
CPUQ = queue.Queue()                # free logical CPUs, one per physical core (--pin)


# ------------------------------------------------------------------ running one process
def run_cpu(cmd, log, timeout=None):
    """run cmd with stdout and stderr to log; return (exit code, CPU seconds, wall seconds, timed out)"""
    t0 = time.time()
    fired = []
    with open(log, "w") as f:
        p = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
        cpu = getattr(TL, "cpu", None)
        if cpu is not None:                                    # --pin: this task runs on its own core
            try:
                os.sched_setaffinity(p.pid, {cpu})
            except OSError:
                pass
        with LOCK:
            RUNNING[p.pid] = p
        timer = None
        if timeout:
            def kill():
                fired.append(1)
                try:
                    p.kill()
                except OSError:
                    pass
            timer = threading.Timer(timeout, kill)
            timer.start()
        _, status, ru = os.wait4(p.pid, 0)
        if timer:
            timer.cancel()
    p.returncode = os.WEXITSTATUS(status) if os.WIFEXITED(status) else -os.WTERMSIG(status)
    with LOCK:
        RUNNING.pop(p.pid, None)
    return p.returncode, ru.ru_utime + ru.ru_stime, time.time() - t0, bool(fired)


def parse_log(log):
    best, hist = {}, {}
    with open(log) as f:
        for line in f:
            m = BEST.match(line)
            if m:
                best[m.group(1)] = dict(gates=int(m.group(2)), n2=int(m.group(3)), n3=int(m.group(4)),
                                        umc180=float(m.group(5)), stm130=float(m.group(6)))
            elif line.startswith("histogram:"):
                for tok in line.split()[1:]:
                    g, c = tok.split(":")
                    hist[int(g)] = hist.get(int(g), 0) + int(c)
    return best, hist


def key_of(measure, b):
    return (b["gates"], b["umc180"]) if measure == "xor2" else (b[measure], b["gates"])


# ------------------------------------------------------------------ one task: BFS, then restarts
def task(a, t):
    """returns the result record, an error record, or None if the run was interrupted"""
    name, D, pair, seed = t["matrix"], t["D"], t["pair"], t["seed"]
    key = (name, pair, seed)
    rule, slack, rule_rr = PAIRS[pair]
    mat = os.path.join(a.bench, name + ".txt")
    ldir = os.path.join(a.out, "logs", name, pair)
    cdir = os.path.join(a.out, "circuits", name, pair)
    os.makedirs(ldir, exist_ok=True)
    os.makedirs(cdir, exist_ok=True)

    def stopped():
        with LOCK:
            STATE.pop(key, None)
        return None

    # 1. beam search over the ties
    if STOP.is_set():
        return None
    with LOCK:
        STATE[key] = ["BFS", time.time(), ""]
    blog = os.path.join(ldir, f"seed{seed}.bfs.log")
    cmd = [a.dbs, mat, str(D)] + rule + ["--beam", str(a.width), "--slack", str(slack), "--seed", str(seed),
                                         "--save", os.path.join(cdir, f"seed{seed}.bfs")]
    rc, cpu, wall, timed_out = run_cpu(cmd, blog, a.timeout)
    if timed_out:
        stopped()
        return dict(t, error=f"BFS over the time limit of {a.timeout:.0f} s")
    if rc < 0 or STOP.is_set():                                # killed: interrupted run
        return stopped()
    bbest, bhist = parse_log(blog)
    if rc != 0 or "xor2" not in bbest:
        stopped()
        return dict(t, error=f"BFS failed (exit {rc}), see {blog}")
    bfs = dict(cpu=round(cpu, 2), wall=round(wall, 2), hist=bhist,
               **{m: bbest[m] for m in MEASURES if m in bbest})

    # 2. random restarts (same rule, except for the cross pair xall1), until their CPU time reaches the BFS time
    used, runs, chunks, rbest, rhist, rwall = 0.0, 0, 0, {}, {}, 0.0
    size = 10
    with LOCK:
        STATE[key] = ["restarts", time.time(), f"BFS {bfs['xor2']['gates']}, budget {cpu:.0f} s"]
    while True:
        if STOP.is_set():
            return stopped()
        chunk_seed = 1000000 + 1000 * seed + chunks          # every restart gets its own seed
        rlog = os.path.join(ldir, f"seed{seed}.rr{chunks}.log")
        pre = os.path.join(cdir, f"seed{seed}.rr{chunks}")
        cmd = [a.dbs, mat, str(D)] + rule_rr + ["--runs", str(size), "--seed", str(chunk_seed), "--save", pre]
        rc, c_cpu, c_wall, _ = run_cpu(cmd, rlog)
        if rc < 0 or STOP.is_set():
            return stopped()
        if rc != 0:
            stopped()
            return dict(t, error=f"restarts failed (exit {rc}), see {rlog}")
        used += c_cpu
        rwall += c_wall
        runs += size
        chunks += 1
        cb, ch = parse_log(rlog)
        for g, c in ch.items():
            rhist[g] = rhist.get(g, 0) + c
        for m in MEASURES:                                     # keep the best circuit of each measure
            if m in cb and (m not in rbest or key_of(m, cb[m]) < key_of(m, rbest[m])):
                rbest[m] = cb[m]
                shutil.copyfile(f"{pre}.{m}.txt", os.path.join(cdir, f"seed{seed}.rr.{m}.txt"))
        for m in MEASURES:
            if os.path.exists(f"{pre}.{m}.txt"):
                os.remove(f"{pre}.{m}.txt")
        with LOCK:
            if key in STATE:
                STATE[key][2] = (f"BFS {bfs['xor2']['gates']}, restarts {rbest['xor2']['gates'] if 'xor2' in rbest else '-'}"
                                 f" after {runs} runs, {used:.0f} of {cpu:.0f} s")
        per_run = used / runs
        left = cpu - used
        if left < per_run or chunks >= 50:
            break
        size = max(1, min(200000, int(0.85 * left / per_run)))
    with LOCK:
        STATE.pop(key, None)
    rr = dict(cpu=round(used, 2), wall=round(rwall, 2), runs=runs, chunks=chunks, hist=rhist,
              **{m: rbest[m] for m in MEASURES if m in rbest})
    return dict(t, bfs=bfs, rr=rr)


def task_pinned(a, t):
    """task() with all its processes (the beam search and its restarts) on one logical CPU of a free core"""
    cpu = CPUQ.get()
    TL.cpu = cpu
    try:
        return task(a, t)
    finally:
        TL.cpu = None
        CPUQ.put(cpu)


# ------------------------------------------------------------------ the plan
def core_cpus():
    """one logical CPU per physical core (the first sibling), from /sys; [] if unknown"""
    seen, cpus = set(), []
    base = "/sys/devices/system/cpu"
    try:
        names = sorted((d for d in os.listdir(base) if re.fullmatch(r"cpu\d+", d)), key=lambda d: int(d[3:]))
    except OSError:
        return []
    allowed = os.sched_getaffinity(0) if hasattr(os, "sched_getaffinity") else None
    for d in names:
        c = int(d[3:])
        if allowed is not None and c not in allowed:
            continue
        try:
            pkg = open(os.path.join(base, d, "topology", "physical_package_id")).read().strip()
            core = open(os.path.join(base, d, "topology", "core_id")).read().strip()
        except OSError:
            return []
        if (pkg, core) not in seen:
            seen.add((pkg, core))
            cpus.append(c)
    return cpus


def physical_cores():
    try:
        cores, cur = set(), {}
        with open("/proc/cpuinfo") as f:
            for line in f:
                if ":" in line:
                    k, v = [x.strip() for x in line.split(":", 1)]
                    cur[k] = v
                elif cur:
                    cores.add((cur.get("physical id"), cur.get("core id")))
                    cur = {}
        if cur:
            cores.add((cur.get("physical id"), cur.get("core id")))
        n = len([c for c in cores if c != (None, None)])
        return n if n > 0 else (os.cpu_count() or 1)
    except OSError:
        return os.cpu_count() or 1


def parse_seeds(s):
    out = []
    for part in s.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            out += list(range(int(lo), int(hi) + 1))
        elif part:
            out.append(int(part))
    return out


def plan(a):
    rows = list(csv.DictReader(open(os.path.join(a.bench, "benchmark.tsv")), delimiter="\t"))
    only = set(a.only.split(",")) if a.only else None
    skip = set(a.skip.split(",")) if a.skip else set()
    tasks = []
    for r in rows:
        name = r["file"][:-4]
        if only is not None and name not in only:
            continue
        if name in skip or (name == HILL_LS16 and not a.with_hill_ls16 and only is None):
            continue
        if not os.path.exists(os.path.join(a.bench, r["file"])):
            print(f"warning: {r['file']} is not in {a.bench}; skipped")
            continue
        D, size = int(r["min_depth"]), r["size"]
        weight = (100 if D >= 4 else 10) if size.startswith("32") else 1      # for the order only
        for pair in a.pairs:
            for seed in a.seeds:
                tasks.append(dict(matrix=name, label=r["label"], size=size, D=D, pair=pair, seed=seed,
                                  width=a.width, weight=weight))
    tasks.sort(key=lambda t: -t["weight"])                    # longest first
    return tasks


# ------------------------------------------------------------------ the run
def hms(sec):
    sec = int(sec)
    return f"{sec // 3600}:{sec % 3600 // 60:02d}"


def status(a, total, done, finished, t0):
    """one status line, and OUT/progress.txt with the running tasks"""
    now = time.time()
    with LOCK:
        st = sorted(STATE.items(), key=lambda kv: kv[1][1])
    nb = sum(1 for _, v in st if v[0] == "BFS")
    line = (f"[{hms(now - t0)} elapsed, {time.strftime('%H:%M')}] {done}/{total} tasks done, {len(st)} running "
            f"({nb} BFS, {len(st) - nb} restarts)")
    tally = {}
    for r in finished:
        b, q = r["bfs"]["xor2"]["gates"], r["rr"]["xor2"]["gates"]
        tally.setdefault(r["pair"], [0, 0, 0])[0 if b < q else (1 if b == q else 2)] += 1
    if tally:
        line += "; BFS lower/equal/higher than restarts: " + ", ".join(
            f"{p} {v[0]}/{v[1]}/{v[2]}" for p, v in sorted(tally.items()))
    with open(os.path.join(a.out, "progress.txt"), "w") as f:
        f.write(line + "\n\nrunning (oldest first):\n")
        for (m, p, sd), (phase, t1, detail) in st:
            f.write(f"  {m:38s} {p:7s} seed {sd}: {phase:8s} {hms(now - t1):>6s}  {detail}\n")
    return line


def cmd_run(a):
    a.pairs = [p for p in a.pairs.split(",") if p]
    bad = [p for p in a.pairs if p not in PAIRS]
    if bad:
        sys.exit(f"unknown pair(s): {', '.join(bad)}; available: {', '.join(PAIRS)}")
    a.seeds = parse_seeds(a.seeds)
    if not os.path.isfile(a.dbs) or not os.access(a.dbs, os.X_OK):
        sys.exit(f"no executable search program at {a.dbs} (build it with 'make', or give --dbs)")
    tasks = plan(a)
    os.makedirs(a.out, exist_ok=True)
    lock_out(a.out)
    res_path = os.path.join(a.out, "results.jsonl")
    done, finished = set(), []
    if os.path.exists(res_path):
        with open(res_path) as f:
            for line in f:
                r = json.loads(line)
                if "error" not in r and r.get("width") == a.width:
                    done.add((r["matrix"], r["pair"], r["seed"]))
                    finished.append(r)
    todo = [t for t in tasks if (t["matrix"], t["pair"], t["seed"]) not in done]
    print(f"{len(tasks)} tasks ({len(done)} already done), {len(todo)} to run on {a.jobs} workers; "
          f"width {a.width}, pairs {','.join(a.pairs)}, seeds {a.seeds[0]}-{a.seeds[-1]}", flush=True)
    pinned = []
    if a.pin:
        pinned = core_cpus()
        if not pinned:
            sys.exit("--pin: the core topology is not available (/sys/devices/system/cpu); run without --pin")
        if a.jobs > len(pinned):
            sys.exit(f"--pin: {a.jobs} jobs but only {len(pinned)} physical cores; use -j {len(pinned)} or fewer")
        for c in pinned[:a.jobs]:
            CPUQ.put(c)
        print(f"--pin: tasks run on the logical CPUs {pinned[:a.jobs]}, one per physical core", flush=True)
    if a.dry_run:
        for t in todo:
            print(f"  {t['matrix']:38s} {t['pair']:7s} seed {t['seed']}")
        return
    print(f"The longest tasks run first, so the first results can take an hour or more. A status line is "
          f"printed every {a.status_every:g} min; {os.path.join(a.out, 'progress.txt')} lists the running tasks. "
          f"Ctrl-C stops the run; the finished tasks are kept.", flush=True)
    with open(os.path.join(a.out, "config.json"), "w") as f:
        json.dump(dict(dbs=os.path.abspath(a.dbs), bench=os.path.abspath(a.bench), width=a.width, pairs=a.pairs,
                       seeds=a.seeds, jobs=a.jobs, host=os.uname().nodename, cores=physical_cores(),
                       pinned_cpus=pinned[:a.jobs] if a.pin else None,
                       started=time.strftime("%Y-%m-%d %H:%M:%S")), f, indent=1)
    t0, n = time.time(), 0
    total = len(tasks)

    def beat():
        wait = min(60.0, 60.0 * a.status_every)
        while not DONE.wait(wait):
            print(status(a, total, len(done) + n, finished, t0), flush=True)
            wait = 60.0 * a.status_every

    threading.Thread(target=beat, daemon=True).start()
    ex = cf.ThreadPoolExecutor(max_workers=a.jobs)
    futs = [ex.submit(task_pinned if a.pin else task, a, t) for t in todo]
    try:
        with open(res_path, "a") as out:
            for fu in cf.as_completed(futs):
                r = fu.result()
                if r is None:
                    continue
                n += 1
                out.write(json.dumps(r) + "\n")
                out.flush()
                if "error" in r:
                    print(f"[{len(done) + n}/{total}] {r['matrix']} {r['pair']} seed {r['seed']}: {r['error']}", flush=True)
                else:
                    finished.append(r)
                    b, q = r["bfs"], r["rr"]
                    print(f"[{len(done) + n}/{total}] {r['matrix']} {r['pair']} seed {r['seed']}: BFS {b['xor2']['gates']} "
                          f"({b['cpu'] / 60:.1f} min) | {q['runs']} restarts {q['xor2']['gates']} ({q['cpu'] / 60:.1f} min)"
                          f"  [{hms(time.time() - t0)} elapsed]", flush=True)
    except KeyboardInterrupt:
        STOP.set()
        DONE.set()
        for fu in futs:
            fu.cancel()
        with LOCK:
            procs = list(RUNNING.values())
        for p in procs:
            try:
                p.kill()
            except OSError:
                pass
        ex.shutdown(wait=True)
        print(f"\nstopped after {hms(time.time() - t0)}; {len(done) + n} of {total} tasks are in {res_path}. "
              f"Run the same command again to continue.", flush=True)
        sys.exit(130)
    DONE.set()
    ex.shutdown(wait=True)
    status(a, total, len(done) + n, finished, t0)
    print(f"\nall tasks done in {hms(time.time() - t0)} (h:mm)\n", flush=True)
    print_summary(a.out, a.tex if a.tex else os.path.join(a.out, "table.tex"), a.paper_run)


_LOCK_FILE = None                   # kept open while 'run' is going


def lock_out(out):
    """one 'run' per output directory: two runs on the same OUT would run every task twice and write
    the same log files"""
    global _LOCK_FILE
    try:
        import fcntl
    except ImportError:                                            # not on Unix: no lock
        return
    path = os.path.join(out, ".run.lock")
    f = open(path, "a+")
    try:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        f.seek(0)
        other = f.read().strip()
        sys.exit(f"another 'run' is using {out} ({other or 'unknown process'}); wait for it, or use another -o")
    f.seek(0)
    f.truncate()
    f.write(f"pid {os.getpid()} on {os.uname().nodename}, started {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.flush()
    _LOCK_FILE = f


# ------------------------------------------------------------------ summary
def load(out, keep_duplicates=False):
    res, seen, dup = [], set(), 0
    with open(os.path.join(out, "results.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            if "error" not in r:
                key = (r["matrix"], r["pair"], r["seed"], r.get("width"))
                if key in seen:
                    dup += 1
                    if not keep_duplicates:
                        continue
                seen.add(key)
                for side in ("bfs", "rr"):
                    r[side]["hist"] = {int(k): v for k, v in r[side]["hist"].items()}
                res.append(r)
    if dup:
        how = "kept (--keep-duplicates)" if keep_duplicates else "skipped: the summary uses the first record of each"
        print(f"warning: {dup} duplicate records (same matrix, pair, seed and width) in {out}/results.jsonl, "
              f"probably from two 'run' processes on the same directory; {how}", flush=True)
    return res


def paper_value(run_dir, name, seed):
    p = os.path.join(run_dir, "forward", f"{name}__all_w2000_s1_seed{seed}.log")
    if not os.path.exists(p):
        return None
    b, _ = parse_log(p)
    return b.get("xor2", {}).get("gates")


def summarize(res):
    groups = {}
    for r in res:
        groups.setdefault((r["matrix"], r["pair"], r["width"]), []).append(r)
    rows = []
    for (name, pair, _), rs in groups.items():
        rs.sort(key=lambda r: r["seed"])
        bx = [r["bfs"]["xor2"]["gates"] for r in rs]
        rx = [r["rr"]["xor2"]["gates"] for r in rs]
        bu = [r["bfs"]["umc180"]["umc180"] for r in rs]
        ru = [r["rr"]["umc180"]["umc180"] for r in rs]
        hist = {}
        for r in rs:
            for g, c in r["rr"]["hist"].items():
                hist[g] = hist.get(g, 0) + c
        nruns = sum(hist.values())
        rr_cpu = sum(r["rr"]["cpu"] for r in rs)
        t_run = rr_cpu / max(1, sum(r["rr"]["runs"] for r in rs))
        t_bfs = statistics.median(r["bfs"]["cpu"] for r in rs)
        target = math.ceil(statistics.median(bx))                      # median of the BFS results (rounded up)
        hits = sum(c for g, c in hist.items() if g <= target)
        p = hits / nruns if nruns else 0.0
        n_budget = statistics.median(r["rr"]["runs"] for r in rs)
        if p > 0:
            p_hit = 1 - (1 - p) ** n_budget
            ratio = (t_run / p) / t_bfs                                # expected CPU time to the target
            ratio_txt = f"{ratio:.2g}"
        else:
            p_hit = 0.0
            ratio = None
            ratio_txt = f">{nruns * t_run / t_bfs:.0f}"
        sd = lambda v: statistics.stdev(v) if len(v) > 1 else 0.0
        rows.append(dict(matrix=name, label=rs[0]["label"], size=rs[0]["size"], D=rs[0]["D"], pair=pair,
                         bfs_sd=sd(bx), rr_sd=sd(rx), bfs_cpu=[r["bfs"]["cpu"] for r in rs], rr_cpu=[r["rr"]["cpu"] for r in rs],
                         width=rs[0]["width"],
                         seeds=[r["seed"] for r in rs], bfs_xor2=bx, rr_xor2=rx,
                         bfs_best=min(bx), bfs_mean=statistics.mean(bx), rr_best=min(rx), rr_mean=statistics.mean(rx),
                         wins=sum(b < q for b, q in zip(bx, rx)), ties=sum(b == q for b, q in zip(bx, rx)),
                         losses=sum(b > q for b, q in zip(bx, rx)),
                         bfs_umc_best=min(bu), rr_umc_best=min(ru), bfs_cpu_median=round(t_bfs, 1),
                         rr_runs_per_budget=n_budget, rr_runs_total=nruns, rr_cpu_per_run=round(t_run, 4),
                         rr_hist=dict(sorted(hist.items())), target=target, p_single=p, p_hit_budget=p_hit,
                         time_to_target_ratio=ratio, time_to_target_txt=ratio_txt))
    rows.sort(key=lambda x: (x["pair"], x["size"], x["matrix"]))
    return rows


def cmd_summary(a):
    print_summary(a.out, a.tex, a.paper_run, a.keep_duplicates)


def print_summary(out, tex, paper_run, keep_duplicates=False):
    if not os.path.exists(os.path.join(out, "results.jsonl")):
        sys.exit(f"no results in {out}")
    rows = summarize(load(out, keep_duplicates))
    if paper_run:                                  # seeds of 'all1' reproduce the paper's stage-1 runs
        checked = differ = missing = 0
        for x in rows:
            if x["pair"] != "all1":
                continue
            if x["width"] != 2000:
                print(f"paper-run check skipped for {x['matrix']}: width {x['width']}, the paper used 2000")
                continue
            for s, g in zip(x["seeds"], x["bfs_xor2"]):
                pv = paper_value(paper_run, x["matrix"], s)
                if pv is None:
                    missing += 1
                    continue
                checked += 1
                if pv != g:
                    differ += 1
                    print(f"warning: {x['matrix']} all1 seed {s}: BFS {g}, paper run {pv}")
        print(f"paper-run check: {checked} BFS runs of pair all1 compared with the stage-1 logs in {paper_run}, "
              f"{differ} differ; {missing} without a log (the paper ran seeds 1-32 of all/2000/1)")
    for pair in sorted({x["pair"] for x in rows}):
        sub = [x for x in rows if x["pair"] == pair]
        print(f"\n## {pair}: beam search against random restarts at equal CPU time\n")
        print("| Matrix | Size | D | Beam XOR2 (best / mean) | Restarts XOR2 (best / mean) | Beam lower / equal / higher "
              "(seeds) | Restarts per budget | Restarts: expected CPU time to reach the beam median, in beam times "
              "| Beam per seed | Restarts per seed | Beam sd | Restarts sd | Beam CPU s (median) |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for x in sub:
            print(f"| {x['label']} | {x['size']} | {x['D']} | {x['bfs_best']} / {x['bfs_mean']:.2f} | "
                  f"{x['rr_best']} / {x['rr_mean']:.2f} | {x['wins']} / {x['ties']} / {x['losses']} | "
                  f"{x['rr_runs_per_budget']:.0f} | {x['time_to_target_txt']} | "
                  f"{' '.join(map(str, x['bfs_xor2']))} | {' '.join(map(str, x['rr_xor2']))} | "
                  f"{x['bfs_sd']:.2f} | {x['rr_sd']:.2f} | {x['bfs_cpu_median']:.0f} |")
        better = sum(x["bfs_mean"] < x["rr_mean"] for x in sub)
        equal = sum(x["bfs_mean"] == x["rr_mean"] for x in sub)
        worse = sum(x["bfs_mean"] > x["rr_mean"] for x in sub)
        bb = sum(x["bfs_best"] < x["rr_best"] for x in sub)
        be = sum(x["bfs_best"] == x["rr_best"] for x in sub)
        bw = sum(x["bfs_best"] > x["rr_best"] for x in sub)
        print(f"\nMean over the seeds: beam lower for {better}, equal for {equal}, higher for {worse} of {len(sub)} matrices.")
        print(f"Best over the seeds: beam lower for {bb}, equal for {be}, higher for {bw} of {len(sub)} matrices.")
        print("(The 23 matrices are 20 distinct problems: Joltik = [SKOP15] inv. 16 = [LW16] inv. 16 and "
              "[LS16] circ. 16 = [LW16] circ. 16, up to permutations.)")
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(rows, f, indent=1)
    if tex:
        write_tex(rows, tex)
        print(f"\nLaTeX table: {tex}")


def write_tex(rows, path):
    """one row per matrix; for every pair, the mean XOR2 count over the seeds of the beam search and of the
    restarts at equal CPU time (the lower one in bold); last row: lower / equal / higher counts"""
    pairs = [p for p in PAIRS if any(x["pair"] == p for x in rows)]
    names = {"all0": "\\emph{all}, $s=0$", "all1": "\\emph{all}, $s=1$", "sfxmd0": "\\emph{sfx} + MD",
             "sfx0": "\\emph{sfx}", "bp1": "\\emph{bp}, $s=1$", "xall1": "\\emph{all}, $s=1$ vs \\emph{sfx} + MD"}
    by = {(x["matrix"], x["pair"]): x for x in rows}
    mats = []
    for x in rows:
        if x["matrix"] not in [m for m, _ in mats]:
            mats.append((x["matrix"], x["label"]))
    cols = "l" + " r r" * len(pairs)
    lines = ["% Written by paper/ties_vs_restarts.py summary: mean XOR2 count over the seeds, beam search (BFS)",
             "% against random restarts with the same rule and the same CPU time (RR); the lower one in bold.",
             f"\\begin{{tabular}}{{@{{}}{cols}@{{}}}}", "\\toprule",
             " & " + " & ".join(f"\\multicolumn{{2}}{{c}}{{{names.get(p, p)}}}" for p in pairs) + " \\\\",
             "".join(f"\\cmidrule(lr){{{2 + 2 * i}-{3 + 2 * i}}}" for i in range(len(pairs))),
             "Matrix & " + " & ".join("BFS & RR" for _ in pairs) + " \\\\", "\\midrule"]
    tally = {p: [0, 0, 0] for p in pairs}
    for m, label in mats:
        cells = []
        for p in pairs:
            x = by.get((m, p))
            if not x:
                cells += ["--", "--"]
                continue
            b, q = x["bfs_mean"], x["rr_mean"]
            tally[p][0 if b < q else (1 if b == q else 2)] += 1
            fb, fq = f"{b:.2f}", f"{q:.2f}"
            if b < q:
                fb = f"\\textbf{{{fb}}}"
            elif q < b:
                fq = f"\\textbf{{{fq}}}"
            cells += [fb, fq]
        lines.append(f"{label} & " + " & ".join(cells) + " \\\\")
    lines.append("\\midrule")
    lines.append("BFS lower / equal / higher & " + " & ".join(
        f"\\multicolumn{{2}}{{c}}{{{t[0]} / {t[1]} / {t[2]}}}" for t in (tally[p] for p in pairs)) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabular}"]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    r = sub.add_parser("run", help="run the comparison")
    r.add_argument("-o", "--out", default="ties_vs_restarts")
    r.add_argument("--dbs", default=os.path.join(ROOT, "bin", "dbs"))
    r.add_argument("--bench", default=os.path.join(ROOT, "bench"))
    r.add_argument("-j", "--jobs", type=int, default=physical_cores(), help="parallel tasks (default: physical cores)")
    r.add_argument("--pin", action="store_true",
                   help="run each task (the beam search and its restarts) on its own physical core, one logical CPU per "
                        "core, so that the hardware threads of a core are never shared (-j at most the number of cores)")
    r.add_argument("--seeds", default="1-4")
    r.add_argument("--pairs", default="all0,all1,sfxmd0")
    r.add_argument("--width", type=int, default=2000)
    r.add_argument("--only", default="", help="comma-separated matrix names (file names without .txt)")
    r.add_argument("--skip", default="")
    r.add_argument("--with-hill-ls16", action="store_true", help=f"also run the [LS16] 32x32 matrix of the HILL circuits, {HILL_LS16}")
    r.add_argument("--timeout", type=float, default=None, help="time limit of a BFS run in seconds")
    r.add_argument("--dry-run", action="store_true")
    r.add_argument("--status-every", type=float, default=10, help="minutes between status lines (default 10)")
    r.add_argument("--tex", default=None, help="LaTeX table written at the end (default OUT/table.tex)")
    r.add_argument("--paper-run", default=None, help="the run/ folder of the paper, to check the all1 seeds at the end")
    s = sub.add_parser("summary", help="tables from the results")
    s.add_argument("out")
    s.add_argument("--tex", default=None, help="write a LaTeX table")
    s.add_argument("--paper-run", default=None, help="the run/ folder of the paper, to check the all1 seeds")
    s.add_argument("--keep-duplicates", action="store_true",
                   help="use every record, also repeated (matrix, pair, seed) records of two runs on one OUT")
    a = ap.parse_args()
    if a.cmd == "run":
        cmd_run(a)
    elif a.cmd == "summary":
        cmd_summary(a)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
