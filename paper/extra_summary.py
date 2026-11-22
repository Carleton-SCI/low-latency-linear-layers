"""Summarize the runs of paper/extra_runs.sh against the paper's run.

usage: python3 paper/extra_summary.py OUT PAPER_RUN
       OUT: the folder of extra_runs.sh (with rank/, area/, stage4/, stage4from2/, stage4fromsfx/, stage4long/, bpbs/, depth4/, sfx1/, aesdepth/, fwdtime/, invmc/ for the experiments
       that were run; beamstats/ is summarized by paper/tables/beam_stats.py); PAPER_RUN: the run/ folder of the paper (or of
       results/run.tar.xz)."""
import csv
import glob
import os
import re
import statistics
import sys

if len(sys.argv) != 3:
    sys.exit(__doc__)
OUT, RUN = sys.argv[1], sys.argv[2]
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BEST = re.compile(r"BEST-(\w+): (\d+) XOR2 gates -> (\d+) XOR2 \+ (\d+) XOR3, umc180 ([\d.]+) GE, stm130 ([\d.]+) GE")
rows = list(csv.DictReader(open(os.path.join(ROOT, "bench", "benchmark.tsv")), delimiter="\t"))
label = {r["file"][:-4]: r["label"] for r in rows}
sfx = {r["file"][:-4]: r["sfx23_xor2"] for r in rows}


def best(log):
    b = {}
    for line in open(log):
        m = BEST.match(line)
        if m:
            b[m.group(1)] = dict(g=int(m.group(2)), umc=float(m.group(5)), stm=float(m.group(6)))
    return b


def logs(pattern):
    return [f for f in glob.glob(pattern) if f.endswith(".log")]


def cpu_hours(folder):
    t = 0.0
    for f in glob.glob(os.path.join(folder, "*.time")):
        m = re.search(r"cpu ([\d.]+) user ([\d.]+) sys", open(f).read())
        if m:
            t += float(m.group(1)) + float(m.group(2))
    return t / 3600


def paper_best(name, measure):
    vals = [best(f).get(measure) for f in logs(os.path.join(RUN, "*", name + "__*")) + logs(os.path.join(RUN, "*", name + ".seed*"))]
    vals = [v for v in vals if v]
    return min(vals, key=lambda v: (v["g"], v["umc"]) if measure == "xor2" else (v[measure[:3]], v["g"])) if vals else None


def section(title):
    print(f"\n## {title}\n")


if os.path.isdir(os.path.join(OUT, "rank")):
    section("rank: beam search ranked by gates so far + Sigma (all/2000/1, seeds 1-8) against the paper's ranking by Sigma")
    print("| Matrix | Sigma: best / mean | gates + Sigma: best / mean | per seed: lower / equal / higher |")
    print("|---|---|---|---|")
    tot = [0, 0, 0]
    for name in label:
        new = {}
        for f in logs(os.path.join(OUT, "rank", name + "__all_w2000_s1_rgates_seed*")):
            sd = int(re.search(r"seed(\d+)\.log$", f).group(1))
            b = best(f)
            if "xor2" in b:
                new[sd] = b["xor2"]["g"]
        if not new:
            continue
        old = {}
        for sd in new:
            p = os.path.join(RUN, "forward", f"{name}__all_w2000_s1_seed{sd}.log")
            if os.path.exists(p) and "xor2" in best(p):
                old[sd] = best(p)["xor2"]["g"]
        sds = sorted(set(new) & set(old))
        w = [sum(new[s] < old[s] for s in sds), sum(new[s] == old[s] for s in sds), sum(new[s] > old[s] for s in sds)]
        tot = [a + b for a, b in zip(tot, w)]
        print(f"| {label[name]} | {min(old[s] for s in sds)} / {statistics.mean(old[s] for s in sds):.2f} | "
              f"{min(new[s] for s in sds)} / {statistics.mean(new[s] for s in sds):.2f} | {w[0]} / {w[1]} / {w[2]} |")
    print(f"\nPer seed, gates + Sigma against Sigma: lower {tot[0]}, equal {tot[1]}, higher {tot[2]}. "
          f"CPU time of these runs: {cpu_hours(os.path.join(OUT, 'rank')):.1f} h.")

for exp, start_dir, what in (("area", "forward", "large neighborhood search with the UMC area as the objective"),
                             ("stage4", "forward", "stage 4 from the best XOR2 circuit of stage 1 only"),
                             ("stage4from2", "bpbs", "stage 4 from the best XOR2 circuit of stage 2 only"),
                             ("stage4fromsfx", "bpbs", "stage 4 from the published circuits of SFX23 (last column: our stage 2)"),
                             ("stage4long", "forward", "stage 4 from the best XOR2 circuit of stage 1, 10 times the iterations")):
    d = os.path.join(OUT, exp)
    if not os.path.isdir(d):
        continue
    section(f"{exp}: {what} (8 seeds)")
    print(f"| Matrix | paper: best XOR2 / best UMC | this run: best XOR2 / best UMC | {start_dir} alone (paper): best XOR2 |")
    print("|---|---|---|---|")
    for name in label:
        fs = logs(os.path.join(d, name + ".seed*"))
        if not fs:
            continue
        bs = [best(f) for f in fs]
        bx = min(b["xor2"]["g"] for b in bs if "xor2" in b)
        bu = min(b["umc180"]["umc"] for b in bs if "umc180" in b)
        px, pu = paper_best(name, "xor2"), paper_best(name, "umc180")
        s1 = min(best(f)["xor2"]["g"] for f in logs(os.path.join(RUN, start_dir, name + "__*")) + logs(os.path.join(RUN, start_dir, name + ".seed*")))
        print(f"| {label[name]} | {px['g']} / {pu['umc']:.2f} | {bx} / {bu:.2f} | {s1} |")
    print(f"\nCPU time of these runs: {cpu_hours(d):.1f} h.")

if os.path.isdir(os.path.join(OUT, "bpbs")):
    section("bpbs: stage 2 with seeds 1-23 (1-8 from the paper's run, 9-23 here)")
    print("| Matrix | SFX23 | stage 2, seeds 1-8 | stage 2, seeds 1-23 | our final (paper) |")
    print("|---|---|---|---|---|")
    for name in label:
        new = logs(os.path.join(OUT, "bpbs", name + ".seed*"))
        if not new:
            continue
        old = logs(os.path.join(RUN, "bpbs", name + ".seed*"))
        o = min(best(f)["xor2"]["g"] for f in old)
        a = min([o] + [best(f)["xor2"]["g"] for f in new if "xor2" in best(f)])
        print(f"| {label[name]} | {sfx[name]} | {o} | {a} | {paper_best(name, 'xor2')['g']} |")
    print(f"\nCPU time of the new runs: {cpu_hours(os.path.join(OUT, 'bpbs')):.1f} h.")

if os.path.isdir(os.path.join(OUT, "depth4")):
    section("depth4: stages 2-4 for [SS16] involutory 32x32 at depth 4")
    name = "ToSC_SarSye16_32_Involutory"
    for st in ("bpbs", "bbeam", "lns"):
        fs = logs(os.path.join(OUT, "depth4", f"{name}.{st}.seed*"))
        if fs:
            bs = [best(f) for f in fs]
            print(f"- {st}: best XOR2 {min(b['xor2']['g'] for b in bs if 'xor2' in b)}, "
                  f"best UMC {min(b['umc180']['umc'] for b in bs if 'umc180' in b):.2f} ({len(fs)} runs)")
    print(f"- stage 1 (paper): best XOR2 {min(best(f)['xor2']['g'] for f in logs(os.path.join(RUN, 'forward', name + '__*')))}")
    print(f"\nCPU time of these runs: {cpu_hours(os.path.join(OUT, 'depth4')):.1f} h.")

if os.path.isdir(os.path.join(OUT, "aesdepth")):
    section("aesdepth: AES MixColumns under the depth bounds 4, 5 and 8 (NIST lists 90 XOR2 gates at depth 4, 88 at depth 5)")
    for D in (4, 5, 8):
        s1 = logs(os.path.join(OUT, "aesdepth", f"AES.d{D}.stage1.log"))
        s4 = logs(os.path.join(OUT, "aesdepth", f"AES.d{D}.lns.seed*"))
        if s1:
            b = best(s1[0])
            line = f"- depth {D}: stage 1 (W = 500, all, s = 1, seed 1): {b['xor2']['g']} XOR2 gates"
            if s4:
                bs = [best(f) for f in s4]
                line += (f"; stage 4 from it: best {min(x['xor2']['g'] for x in bs if 'xor2' in x)} XOR2 gates, "
                         f"best UMC {min(x['umc180']['umc'] for x in bs if 'umc180' in x):.2f} GE ({len(s4)} runs)")
            print(line)
    print(f"\nCPU time of these runs: {cpu_hours(os.path.join(OUT, 'aesdepth')):.1f} h.")

if os.path.isdir(os.path.join(OUT, "sfx1")):
    section("sfx1: stage 1 with the sfx rule at slack 1 (8 seeds) against the configurations of Table 3")
    print("| Matrix | sfx/2000/1 (new): best / mean | all/2000/1 (seeds 1-8) | bp/2000/1 | sfx/2000/0 |")
    print("|---|---|---|---|---|")
    for name in label:
        fs = logs(os.path.join(OUT, "sfx1", name + "__sfx_w2000_s1_seed*"))
        if not fs:
            continue
        new = [best(f)["xor2"]["g"] for f in fs if "xor2" in best(f)]
        cols = []
        for cfg, seeds in (("all_w2000_s1", range(1, 9)), ("bp_w2000_s1", range(1, 9)), ("sfx_w2000_s0", range(1, 9))):
            v = []
            for sd in seeds:
                p = os.path.join(RUN, "forward", f"{name}__{cfg}_seed{sd}.log")
                if os.path.exists(p) and "xor2" in best(p):
                    v.append(best(p)["xor2"]["g"])
            cols.append(f"{min(v)} / {statistics.mean(v):.2f}" if v else "-")
        print(f"| {label[name]} | {min(new)} / {statistics.mean(new):.2f} | " + " | ".join(cols) + " |")
    print(f"\nCPU time of these runs: {cpu_hours(os.path.join(OUT, 'sfx1')):.1f} h.")


def time_of(prefix):
    """(cpu seconds, wall seconds) from the .time file of /usr/bin/time, or None"""
    try:
        t = open(prefix + ".time").read()
    except OSError:
        return None
    m = re.search(r"cpu ([\d.]+) user ([\d.]+) sys, wall ([\d.]+) s", t)
    return (float(m.group(1)) + float(m.group(2)), float(m.group(3))) if m else None


if os.path.isdir(os.path.join(OUT, "fwdtime")):
    section("fwdtime: one randomized forward run (bp rule) for AES MixColumns, seeds 1-3 (Sec. 6)")
    for D in (3, 4, 5):
        cells = []
        for sd in (1, 2, 3):
            f = os.path.join(OUT, "fwdtime", f"AES.d{D}.bp.seed{sd}")
            t, b = time_of(f), (best(f + ".log") if os.path.exists(f + ".log") else {})
            tr = re.findall(r"TRACE cpu=([\d.]+)", open(f + ".err").read()) if os.path.exists(f + ".err") else []
            cells.append(f"seed {sd}: " + (f"{t[0]:.2f} s CPU, {t[1]:.2f} s wall" if t else
                                           (f"{float(tr[-1]):.2f} s CPU (TRACE)" if tr else "no time")) +
                         (f", {b['xor2']['g']} XOR2" if "xor2" in b else ", no result"))
        print(f"- depth {D}: " + "; ".join(cells))

if os.path.isdir(os.path.join(OUT, "invmc")):
    section("invmc: one randomized forward run for AES InvMixColumns at depth 5 (Sec. 6; NIST lists 119 XOR2 gates)")
    for t in ("bp", "sfx", "all"):
        f = os.path.join(OUT, "invmc", f"AESinv.d5.{t}.seed1")
        if not os.path.exists(f + ".err"):
            continue
        b = best(f + ".log") if os.path.exists(f + ".log") else {}
        last = added = None
        for line in open(f + ".err"):
            if line.startswith("PROGRESS"):
                last = line.split(None, 1)[1].strip()
                if " added " in line:
                    added = last
        tm = time_of(f)
        res = f"finished: {b['xor2']['g']} XOR2 gates" if "xor2" in b else "not finished (time limit)"
        print(f"- rule {t}: {res}; " + (f"{tm[0] / 3600:.2f} h CPU, {tm[1] / 3600:.2f} h wall; " if tm else "") +
              f"last gate added: {added or 'none'}; last line: {last or 'none'}")
