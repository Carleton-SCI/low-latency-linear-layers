"""Diagnostics of the beam search (Sec. 5.4 of the paper): runs of dbs with --stats.

usage: python3 paper/tables/beam_stats.py RUN_DIR
       RUN_DIR: the folder of paper/extra_runs.sh beamstats (or runs/beamstats of results/server_runs.tar.xz), with
       <matrix>.<config>.w<W>.s<seed>.{log|out} (the output of dbs) and .err (the STAT lines of --stats)
configs: sfxmd0 = --tie sfx --md --slack 0, sfx0 = --tie sfx --slack 0, all0 = --tie all --slack 0,
         sfxmd1 = --tie sfx --md --slack 1, sfxmd0rg = --tie sfx --md --slack 0 --rank gates
prints, per matrix, configuration and width: the best XOR2 count of each seed and their mean; then, for the
beam searches with W = 2000 of the 16x16 matrices: the share of the iterations with fewer than W near-ties
(options kept within the slack), the median number of sequences that the new sequences of an iteration
descend from, and whether the best circuit is among the circuits completed first (STATEND line)."""
import glob
import os
import re
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def size(mat):
    with open(os.path.join(ROOT, "bench", mat + ".txt")) as f:
        return len(f.readline().split())


NAME = re.compile(r"(.+)\.(sfxmd0rg|sfxmd0|sfxmd1|sfx0|all0)\.w(\d+)\.s(\d+)$")
BEST = re.compile(r"BEST-xor2: (\d+) XOR2")


def load(run):
    runs = {}
    for err in sorted(glob.glob(os.path.join(run, "*.err"))):
        base = err[:-4]
        m = NAME.match(os.path.basename(base))
        if not m:
            continue
        outf = base + ".log" if os.path.exists(base + ".log") else base + ".out"
        if not os.path.exists(outf):
            continue
        b = BEST.search(open(outf).read())
        if not b:
            continue
        it, end = [], None
        for line in open(err):
            if line.startswith("STAT "):
                it.append({k: int(v) for k, v in (kv.split("=") for kv in line.split()[1:])})
            elif line.startswith("STATEND"):
                end = {k: int(v) for k, v in (kv.split("=") for kv in line.split()[1:])}
        runs[(m.group(1), m.group(2), int(m.group(3)), int(m.group(4)))] = dict(best=int(b.group(1)), it=it, end=end)
    return runs


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    runs = load(sys.argv[1])
    groups = {}
    for (mat, cfg, w, seed), r in runs.items():
        groups.setdefault((mat, cfg, w), {})[seed] = r
    print("matrix | config | W | best XOR2 per seed | mean")
    for (mat, cfg, w), rs in sorted(groups.items()):
        best = [rs[s]["best"] for s in sorted(rs)]
        print(f"{mat} | {cfg} | {w} | {' '.join(map(str, best))} | {statistics.mean(best):.2f}")
    print()
    for cfg in ("sfxmd0", "sfx0", "all0"):
        share, parents, n, first = [], [], 0, 0
        for (mat, c, w), rs in sorted(groups.items()):
            if c != cfg or w != 2000 or size(mat) != 16:
                continue
            kept = [x["kept"] for s in rs for x in rs[s]["it"]]
            pn = [x["parents_new"] for s in rs for x in rs[s]["it"]]
            share.append(sum(k < w for k in kept) / len(kept))
            parents.append(statistics.median(pn))
            for s in rs:
                e = rs[s]["end"]
                n += 1
                first += e["best_pruned"] == e["best_pruned_at_unpruned_min"]
        if share:
            print(f"{cfg}, W = 2000, 16x16 ({len(share)} matrices): share of iterations with fewer than W near-ties "
                  f"{min(share):.0%} to {max(share):.0%}; median parents of the new sequences {min(parents):g} to "
                  f"{max(parents):g}; best circuit among the first completed in {first} of {n} runs")
    sm = {}
    for (mat, cfg, w), rs in groups.items():
        if cfg == "sfxmd0" and size(mat) == 16 and len(rs) == 4:
            sm.setdefault(mat, {})[w] = statistics.mean(r["best"] for r in rs.values())
    same = sum(1 for v in sm.values() if 100 in v and 2000 in v and v[100] == v[2000])
    tot = sum(1 for v in sm.values() if 100 in v and 2000 in v)
    print(f"sfxmd0, 16x16: W = 100 and W = 2000 give the same mean for {same} of {tot} matrices")
    for other, what in (("sfxmd1", "s = 1 against s = 0 (W = 2000), mean per matrix"),
                        ("sfxmd0rg", "ranking by gates so far + Delta against Delta (W = 2000, s = 0), per run")):
        lo = eq = hi = 0
        for (mat, cfg, w), rs in sorted(groups.items()):
            if cfg != other or w != 2000 or (mat, "sfxmd0", 2000) not in groups:
                continue
            base = groups[(mat, "sfxmd0", 2000)]
            sds = sorted(set(rs) & set(base))
            if not sds:
                continue
            if other == "sfxmd1":
                a, b = statistics.mean(rs[s]["best"] for s in sds), statistics.mean(base[s]["best"] for s in sds)
                lo += a < b; eq += a == b; hi += a > b
            else:
                for s in sds:
                    lo += rs[s]["best"] < base[s]["best"]; eq += rs[s]["best"] == base[s]["best"]; hi += rs[s]["best"] > base[s]["best"]
        if lo + eq + hi:
            print(f"{other}: {what}: lower {lo}, equal {eq}, higher {hi}")


if __name__ == "__main__":
    main()
