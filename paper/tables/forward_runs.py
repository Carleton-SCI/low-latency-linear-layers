"""Number of forward runs done by stages 2 and 3 (Sec. 5.3 of the paper), a measure of work that does not depend on
the load of the server.

In a stage-2 or stage-3 run (dbs --bpbs, with or without --bbeam), the histogram line of the log counts the R
forward runs of the initial prefix one by one, and then one entry per evaluated prefix (the best of its R forward
runs). Hence, the number of forward runs of a run is R + (histogram total - R) * R, with R = 4 for stage 2 and
R = 2 for stage 3 (the --evals of paper/run_benchmark.sh).

usage: python3 paper/tables/forward_runs.py RUN_DIR [MATRIX ...]
       RUN_DIR: the run/ folder of the paper (or of results/run.tar.xz); default matrices: the 11 benchmark
       matrices of size 32x32 and depth 3"""
import glob
import os
import re
import sys

H = re.compile(r"histogram: (.*)")
DEFAULT = ["AES", "FSE_SKOP15_32_Hadamard", "FSE_LiuSim16_32_Circulant_SFX", "FSE_LiWang16_32",
           "C_BeiKraLea16_32_Circulant", "ToSC_SarSye16_32_Toeplitz", "ePrint_JeaPeySim_32_Subfield",
           "FSE_SKOP15_32_Involutory", "FSE_LiWang16_32_Involutory", "KLSW", "168_130"]


def histogram_total(path):
    with open(path) as f:
        for line in f:
            m = H.match(line)
            if m:
                return sum(int(x.split(":")[1]) for x in m.group(1).split())
    return None


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    run, mats = sys.argv[1], sys.argv[2:] or DEFAULT
    for stage, folder, R in (("stage 2", "bpbs", 4), ("stage 3", "bbeam", 2)):
        runs = total = 0
        for m in mats:
            for log in sorted(glob.glob(os.path.join(run, folder, m + ".seed*.log"))):
                t = histogram_total(log)
                if t is None:
                    sys.exit(f"{log}: no histogram line")
                runs += 1
                total += R + (t - R) * R
        print(f"{stage}: {runs} runs, {total} forward runs ({total / 1e6:.2f} million)")


if __name__ == "__main__":
    main()
