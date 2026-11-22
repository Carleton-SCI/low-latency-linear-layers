#!/bin/bash
# Full benchmark run of dbs on a multi-core server. Every step is resumable (finished runs are
# skipped), so the same command can be re-run after an interruption.
#
#   ./paper/run_benchmark.sh [WORKERS]                 all steps, results in ./run (override with OUT=dir)
#   STEPS="3 5" ./paper/run_benchmark.sh 48            only some steps
#   ONLY=AES,KLSW ./paper/run_benchmark.sh 48          only some matrices (file names in bench/ without .txt)
#   SKIP=AES ./paper/run_benchmark.sh 40               all matrices except some
#   HILL=path/to/hill_repo SFX=path/to/sfx_repo STEPS=5 ./paper/run_benchmark.sh
#                                             summary plus the two cross-checks of summarize.py
#
# Steps (matrices from bench/benchmark.tsv; steps 2-4 need minimum depth 3)
#   1  forward search grid, beam search over the ties: 56 runs per matrix (bench_run.py --grid full)
#   2  Shi-Feng-Xu BPBS baseline: random backward-split trajectories, IBPD-MD evaluations, SEEDS seeds
#   3  beam search over the backward splits (--bbeam), BB_SEEDS seeds
#   4  large neighborhood search (LNS) from the best circuit of steps 1-3, LNS_SEEDS seeds
#   5  summary: best XOR2 count and best XOR2/XOR3 areas, every netlist re-verified independently
# Steps 2-4 write their job list to OUT/stepN.jobs, print the number of jobs, and print one line
# per finished job (time, job, best XOR2 count, minutes).
#
# Run times measured for the paper (one core per run, 2x Xeon E5-2670 v3, up to 48 runs at a time;
# Table 7 of the paper), median per run for 16x16 / 32x32 matrices of depth 3:
#   step 1: 26 s / 20.6 min    step 2: 38 / 388 min    step 3: 26 / 347 min    step 4: 2 / 78 min
# All steps for the 24 matrices: about 2,040 core-hours.
set -euo pipefail
J=${1:-$(nproc)}
PKG=$(cd "$(dirname "$0")" && pwd)          # paper/ (scripts)
ROOT=$(cd "$PKG/.." && pwd)                 # repository root (bin/dbs, bench/)
OUT=${OUT:-$PWD/run}
STEPS=${STEPS:-"1 2 3 4 5"}; GRID=${GRID:-full}
SEEDS=${SEEDS:-8}; TRAJ16=${TRAJ16:-2000}; TRAJ32=${TRAJ32:-300}
BB_SEEDS=${BB_SEEDS:-8}; BB_REPS16=${BB_REPS16:-40}; BB_REPS32=${BB_REPS32:-6}
LNS_SEEDS=${LNS_SEEDS:-8}
ONLY=${ONLY:-}; SKIP=${SKIP:-}
[ -x "$ROOT/bin/dbs" ] || make -C "$ROOT" bin/dbs
mkdir -p "$OUT"
names() {
  cut -f1 "$ROOT/bench/benchmark.tsv" | tail -n +2 | sed 's/\.txt$//' |
    { if [ -n "$ONLY" ]; then grep -xE "$(echo "$ONLY" | tr ',' '|')" || true; else cat; fi; } |
    { if [ -n "$SKIP" ]; then grep -vxE "$(echo "$SKIP" | tr ',' '|')" || true; else cat; fi; }
}
depth_of() { awk -F'\t' -v f="$1.txt" '$1 == f {print $4}' "$ROOT/bench/benchmark.tsv"; }
size_of() { awk 'NR == 1 {print NF; exit}' "$ROOT/bench/$1.txt"; }
best_of() {   # best XOR2 circuit of matrix $1 found by steps 1-3 (empty if none)
  local files=()
  shopt -s nullglob
  files=("$OUT"/forward/"$1"__*.log "$OUT"/bpbs/"$1".seed*.log "$OUT"/bbeam/"$1".seed*.log)
  shopt -u nullglob
  [ ${#files[@]} -gt 0 ] || return 0
  awk '/^BEST-xor2:/ {print $2, FILENAME}' "${files[@]}" | sort -n -k1,1 |
    awk 'NR == 1 {f = $2; sub(/\.log$/, ".xor2.txt", f); print f}'
}
run_one() {   # $1 = output prefix, rest = dbs arguments; skipped when already finished
  local f=$1; shift
  [ -f "$f.log" ] && grep -q BEST-stm130 "$f.log" && return 0
  local t0=$SECONDS r
  "$ROOT/bin/dbs" "$@" --save "$f" > "$f.log" 2> "$f.err" || true
  r=$(grep -m1 -oE 'BEST-xor2: [0-9]+' "$f.log" || true)
  echo "$(date '+%m-%d %H:%M') $(basename "$f"): ${r:-NO RESULT} ($(( (SECONDS - t0) / 60 )) min)"
}
export -f run_one; export PKG ROOT
run_jobs() {   # $1 = step number; runs every line of OUT/step$1.jobs with J workers
  local jf="$OUT/step$1.jobs" n
  n=$(wc -l < "$jf")
  echo "step $1: $n jobs, $J workers, started $(date '+%m-%d %H:%M')"
  [ "$n" -gt 0 ] || return 0
  xargs -P "$J" -L 1 bash -c 'run_one "$@"' _ < "$jf"
  echo "step $1: finished $(date '+%m-%d %H:%M')"
}

for s in $STEPS; do case $s in
1)  python3 "$PKG/bench_run.py" "$ROOT/bench" "$OUT/forward" --workers "$J" --grid "$GRID" \
        ${ONLY:+--only "$ONLY"} --dbs "$ROOT/bin/dbs" | tee "$OUT/forward.txt" ;;
2)  mkdir -p "$OUT/bpbs"; : > "$OUT/step2.jobs"
    for m in $(names); do
      [ "$(depth_of "$m")" = 3 ] || continue
      if [ "$(size_of "$m")" -le 16 ]; then T=$TRAJ16; MS=16; else T=$TRAJ32; MS=12; fi
      for sd in $(seq 1 "$SEEDS"); do
        echo "$OUT/bpbs/$m.seed$sd $ROOT/bench/$m.txt 3 --bpbs $T --maxsplits $MS --evals 4 --tie sfx --md --seed $sd" >> "$OUT/step2.jobs"
      done
    done
    run_jobs 2 ;;
3)  mkdir -p "$OUT/bbeam"; : > "$OUT/step3.jobs"
    for m in $(names); do
      [ "$(depth_of "$m")" = 3 ] || continue
      if [ "$(size_of "$m")" -le 16 ]; then R=$BB_REPS16; else R=$BB_REPS32; fi
      for sd in $(seq 1 "$BB_SEEDS"); do
        echo "$OUT/bbeam/$m.seed$sd $ROOT/bench/$m.txt 3 --bpbs $R --bbeam 8 --bopts 16 --evals 2 --maxsplits 14 --tie sfx --md --seed $sd" >> "$OUT/step3.jobs"
      done
    done
    run_jobs 3 ;;
4)  mkdir -p "$OUT/lns"; : > "$OUT/step4.jobs"
    for m in $(names); do
      [ "$(depth_of "$m")" = 3 ] || continue
      best=$(best_of "$m") || best=""
      if [ -z "$best" ] || [ ! -f "$best" ]; then echo "step 4: no circuit for $m, skipped"; continue; fi
      if [ "$(size_of "$m")" -le 16 ]; then W=1000; IT=300; else W=200; IT=150; fi
      for sd in $(seq 1 "$LNS_SEEDS"); do
        echo "$OUT/lns/$m.seed$sd $ROOT/bench/$m.txt 3 --lns $best --iters $IT --fix 4:28 --beam $W --slack 1 --tie all --seed $sd" >> "$OUT/step4.jobs"
      done
    done
    run_jobs 4 ;;
5)  python3 "$PKG/summarize.py" "$ROOT/bench" "$OUT" ${HILL:+--hill "$HILL"} ${SFX:+--sfx "$SFX"} | tee "$OUT/summary.txt" ;;
*)  echo "unknown step $s"; exit 1 ;;
esac; done
