"""Results of the extra stage-4 runs that start from the published circuits of SFX23 (Sec. 5.3 of the paper;
paper/extra_runs.sh stage4fromsfx), for Table 3 and the text.

usage: python3 paper/tables/extra_stage4.py RUN_DIR OUT.json
       RUN_DIR: the folder of the experiment (extra_runs/stage4fromsfx), with <matrix>.seed<k>.log files
writes, per matrix: the best XOR2 count of each seed, the best XOR2 count, and the best UMC area with the
XOR2 count of its circuit (over the circuits recorded by these runs)"""
import glob
import json
import os
import re
import sys

BEST = re.compile(r"BEST-(\w+): (\d+) XOR2 gates -> (\d+) XOR2 \+ (\d+) XOR3, umc180 ([\d.]+) GE, stm130 ([\d.]+) GE")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    run, out = sys.argv[1:]
    res = {}
    for f in sorted(glob.glob(os.path.join(run, "*.seed*.log"))):
        m = re.match(r"(.+)\.seed(\d+)\.log$", os.path.basename(f))
        name, seed = m.group(1), int(m.group(2))
        b = {}
        for line in open(f):
            x = BEST.match(line)
            if x:
                b[x.group(1)] = dict(g=int(x.group(2)), umc=float(x.group(5)))
        if "xor2" not in b or "umc180" not in b:
            continue
        r = res.setdefault(name, dict(seeds={}, xor2=10 ** 9, umc180=1e9, umc_g=None))
        r["seeds"][str(seed)] = b["xor2"]["g"]
        r["xor2"] = min(r["xor2"], b["xor2"]["g"])
        if (b["umc180"]["umc"], b["umc180"]["g"]) < (r["umc180"], r["umc_g"] or 10 ** 9):
            r["umc180"], r["umc_g"] = b["umc180"]["umc"], b["umc180"]["g"]
    json.dump(res, open(out, "w"), indent=1, sort_keys=True)
    for name, r in sorted(res.items()):
        seeds = " ".join(f"{r['seeds'][k]}" for k in sorted(r["seeds"], key=int))
        print(f"{name:32s} best XOR2 {r['xor2']:4d} (seeds: {seeds})  best UMC {r['umc180']:.2f} ({r['umc_g']})")


if __name__ == "__main__":
    main()
