"""Run times of stages 2-4 (Table 7 of the paper) from the files of a benchmark run: a run starts
when its .log file is created (birth time) and ends at its last write. Needs a file system that
records birth times (stat %W); stage 1 times are in the logs themselves ('# <seconds> s').

usage: python3 paper/tables/stage_times.py RUN_DIR [--csv FILE]   (--csv: also write the time of every run)"""
import glob
import os
import statistics
import subprocess
import sys

if len(sys.argv) not in (2, 4) or (len(sys.argv) == 4 and sys.argv[2] != "--csv"):
    sys.exit(__doc__)
run = sys.argv[1]
rows = []
for stage in ("bpbs", "bbeam", "lns"):
    by = {16: [], 32: []}
    for log in glob.glob(os.path.join(run, stage, "*.seed*.log")):
        born, mod = subprocess.run(["stat", "-c", "%W %Y", log], capture_output=True, text=True).stdout.split()
        if born == "0":
            sys.exit("this file system does not record birth times")
        name = os.path.basename(log).split(".seed")[0]
        size = 32 if ("32" in name or name in ("AES", "KLSW", "168_130")) else 16
        by[size].append((int(mod) - int(born)) / 60)
        rows.append((stage, name, os.path.basename(log).split(".seed")[1].split(".")[0], int(born), int(mod),
                     round((int(mod) - int(born)) / 60, 2)))
    for size, t in by.items():
        if t:
            print(f"{stage:6s} {size}x{size}: {len(t):3d} runs, median {statistics.median(t):6.0f} min, "
                  f"total {sum(t) / 60:7.1f} h")
if len(sys.argv) == 4:
    with open(sys.argv[3], "w") as f:
        f.write("stage,matrix,seed,birth_unix,last_write_unix,minutes\n")
        for r in sorted(rows):
            f.write(",".join(map(str, r)) + "\n")
    print(f"per-run times: {sys.argv[3]} ({len(rows)} runs)")
