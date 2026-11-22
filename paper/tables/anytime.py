"""Anytime curves of the runs of paper/extra_runs.sh (dbs --trace): the best XOR2 count (or UMC area) of each run
as a function of its CPU time.

usage: python3 paper/tables/anytime.py RUN_DIR [xor2|umc180]
       RUN_DIR: the folder of one experiment of paper/extra_runs.sh, with <matrix>.<...>.err files that hold the
       TRACE lines of dbs --trace, and the .time files of /usr/bin/time when it was available
prints, per matrix: the number of runs, the best value, the median CPU time at which a run reached its own best,
and the mean over the runs of the best value reached within 1/8, 1/4, 1/2 and all of the median CPU time T of
the runs (T from the .time files, or from the last TRACE line when there is none)."""
import glob
import os
import re
import statistics
import sys

TRACE = re.compile(r"TRACE cpu=([\d.]+) measure=(\w+) xor2=(\d+) umc180=([\d.]+)")
TIME = re.compile(r"cpu ([\d.]+) user ([\d.]+) sys")


def main():
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    run = sys.argv[1]
    measure = sys.argv[2] if len(sys.argv) == 3 else "xor2"
    curves = {}
    for err in sorted(glob.glob(os.path.join(run, "*.err"))):
        base = os.path.basename(err)[:-4]
        mat = re.split(r"\.seed\d+|__|\.lns\.|\.bpbs\.|\.bbeam\.|\.stage1", base)[0]
        pts = []
        for line in open(err):
            m = TRACE.match(line)
            if m and m.group(2) == measure:
                v = int(m.group(3)) if measure == "xor2" else float(m.group(4))
                pts.append((float(m.group(1)), v))
        if not pts:
            continue
        t_end = pts[-1][0]
        tf = err[:-4] + ".time"
        if os.path.exists(tf):
            m = TIME.search(open(tf).read())
            if m:
                t_end = float(m.group(1)) + float(m.group(2))
        curves.setdefault(mat, []).append((pts, t_end))
    print(f"measure: {measure}")
    print("| Matrix | runs | best | median CPU s to the run's best | mean best within T/8, T/4, T/2, T | T (s) |")
    print("|---|---|---|---|---|---|")
    for mat, cs in sorted(curves.items()):
        best = min(p[-1][1] for p, _ in cs)
        t_best = statistics.median(next(t for t, v in p if v == p[-1][1]) for p, _ in cs)
        T = statistics.median(t for _, t in cs)
        def at(p, t):
            vals = [v for tt, v in p if tt <= t]
            return vals[-1] if vals else None
        cols = []
        for frac in (1 / 8, 1 / 4, 1 / 2, 1):
            vals = [at(p, frac * T) for p, _ in cs]
            vals = [v for v in vals if v is not None]
            cols.append(f"{statistics.mean(vals):.2f} ({len(vals)})" if vals else "-")
        print(f"| {mat} | {len(cs)} | {best} | {t_best:.0f} | {' / '.join(cols)} | {T:.0f} |")
    print("(n): number of runs with a circuit by then")


if __name__ == "__main__":
    main()
