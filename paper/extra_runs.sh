#!/usr/bin/env bash
# Additional experiments. Each experiment writes its runs (logs and circuits, as in run_benchmark.sh)
# into $OUT/<experiment>/ and can be restarted: finished runs are kept.
#
#   PAPER_RUN=path/to/run ./paper/extra_runs.sh EXPERIMENT [JOBS]
#
# PAPER_RUN is the run/ folder of the paper's benchmark run (or results/run.tar.xz extracted): the
# starting circuits are taken from it (not needed by beamstats, sfx1, aesdepth, fwdtime and invmc). JOBS
# defaults to nproc.
#
#   rank      stage 1 with the options ranked by gates so far + Sigma (dbs --rank gates), all/2000/1,
#             seeds 1-8, all 23 matrices (Sec. 6 of the paper: 88 runs on the 16x16 matrices).
#             About 55 core-hours; [SS16] involutory 32x32 (depth 4) dominates the wall time (~2.5 h).
#   area      stage 4 with the UMC area as the objective (--obj umc), 8 seeds, from the circuit with the
#             smallest UMC area of the paper's run; same widths and iterations as stage 4. All 23 matrices.
#             About 150-200 core-hours.
#   stage4    stage 4 as in the paper (--obj xor), 8 seeds, but from the best XOR2 circuit of stage 1 only
#             (to measure stage 4 without stages 2-3). The 22 matrices of depth 3. About 150 core-hours.
#   stage4from2  the same, from the best XOR2 circuit of stage 2 only. About 150 core-hours.
#   stage4fromsfx  stage 4 from the published circuits of SFX23 (needs SFX_START, the folder written by
#             python3 paper/audit/sfx_to_dbs.py SFX_RESULTS_DIR bench SFX_START). The 22 matrices of depth 3.
#             About 150 core-hours.
#   stage4long   stage 4 from the best XOR2 circuit of stage 1, with 10 times the iterations (about the
#             per-matrix time of stages 2-3), for 4 matrices of size 32x32 whose record came from stage 4
#             ([BKL16], [JPST17] subfield, [SKOP15] involutory, [KLSW17]). About 500 core-hours.
#   bpbs      stage 2 (our BPBS) with seeds 9-23 for the 12 matrices of size 32x32 and depth 3, so that
#             stage 2 gets about the per-matrix budget of the whole search. About 1,200 core-hours.
#   depth4    stages 2-4 for [SS16] involutory 32x32 at its minimum depth 4 (8 seeds each). Several
#             hundred core-hours (stage-1 runs at depth 4 took about 7 times longer than at depth 3).
#   beamstats the beam searches of Sec. 5.4 with diagnostics (dbs --stats; summarized by
#             paper/tables/beam_stats.py): for the 11 matrices of size 16x16, seeds 1-4, the sfx rule with
#             the MD filter and s = 0 at W = 1, 10, 100 and 2000, and the sfx and all rules at s = 0 with
#             W = 2000; for [LW16] involutory 32x32, the sfx rule with the MD filter and s = 0 at
#             W = 1, 10, 100 (seeds 1-4) and 2000 (seeds 1-2); and, with W = 2000, the sfx rule with the MD
#             filter at s = 1 (sfxmd1: 16x16 and [LW16] involutory 32x32, seeds 1-4) and ranked by gates so far
#             + Sigma at s = 0 (sfxmd0rg: 16x16 seeds 1-4, [LW16] involutory 32x32 seeds 1-2). About 4
#             core-hours.
#   sfx1      stage 1 with the sfx rule at slack 1 (--tie sfx --beam 2000 --slack 1), seeds 1-8, all 23
#             matrices, to compare the tie rules of Table 3 at the same slack. About 55 core-hours ([SS16]
#             involutory 32x32, at depth 4, dominates).
#   aesdepth  AES MixColumns under the depth bounds 4, 5 and 8: stage 1 (rule all, W = 500, s = 1,
#             seed 1; Sec. 6 of the paper), then stage 4 from its best circuit (8 seeds, as for the
#             32x32 matrices). Stage 1 takes minutes; stage 4 about 3 x 8 runs of 1-2 hours.
#   fwdtime   the time of one randomized forward run (bp rule, --runs 1) for AES MixColumns at the depth bounds
#             3, 4 and 5, seeds 1-3 (Sec. 6 of the paper). Seconds; run it with JOBS = 1 on an idle server.
#   invmc     AES InvMixColumns (examples/aes_invmixcolumns.mat.txt) at its minimum depth 5: one randomized
#             forward run (seed 1) with each of the rules bp, sfx and all, with --progress (CPU time after each
#             target scanned and each gate added) and a time limit TIMEOUT (default 14h) per run (Sec. 6).
#             3 cores, overnight.
#
# Every run but those of beamstats has --trace: the .err file of a run holds the CPU time of each improvement of
# its best XOR2 count and UMC area (anytime curves: python3 paper/tables/anytime.py $OUT/EXPERIMENT).
# Options: ONLY=A,B (matrix names) restricts the matrices; DRY=1 prints the runs without starting them.
# Afterwards: python3 paper/extra_summary.py $OUT $PAPER_RUN
set -euo pipefail
EXP=${1:?experiment: rank, area, stage4, stage4from2, stage4fromsfx, stage4long, bpbs, depth4, beamstats, sfx1, aesdepth, fwdtime or invmc}
PKG=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$PKG/.." && pwd)
J=${2:-$(nproc)}
OUT=${OUT:-$PWD/extra_runs}
case $EXP in beamstats|aesdepth|sfx1|fwdtime|invmc) ;; *)
  PAPER_RUN=${PAPER_RUN:?set PAPER_RUN to the run/ folder of the paper (or of results/run.tar.xz)} ;; esac
[ -x "$ROOT/bin/dbs" ] || make -C "$ROOT" bin/dbs
mkdir -p "$OUT/$EXP"
names() { cut -f1 "$ROOT/bench/benchmark.tsv" | tail -n +2 | sed 's/\.txt$//' | grep -vx 'FSE_LiuSim16_32_Circulant' |
  { if [ -n "${ONLY:-}" ]; then grep -xE "$(echo "$ONLY" | tr ',' '|')" || true; else cat; fi; }; }
depth_of() { awk -F'\t' -v f="$1.txt" '$1 == f {print $4}' "$ROOT/bench/benchmark.tsv"; }
size_of() { awk 'NR == 1 {print NF; exit}' "$ROOT/bench/$1.txt"; }
best_in() {   # best_in MEASURE MATRIX DIRS...: recorded circuit with the smallest MEASURE among the logs in DIRS
  local measure=$1 m=$2; shift 2
  local files=() d
  shopt -s nullglob
  for d in "$@"; do files+=("$d/$m"__*.log "$d/$m".seed*.log); done
  shopt -u nullglob
  [ ${#files[@]} -gt 0 ] || return 0
  case $measure in
    xor2)   awk '/^BEST-xor2:/ {print $2, FILENAME}' "${files[@]}" | sort -n -k1,1 | awk 'NR == 1 {f = $2; sub(/\.log$/, ".xor2.txt", f); print f}' ;;
    umc180) awk '/^BEST-umc180:/ {print $12, $2, FILENAME}' "${files[@]}" | sort -n -k1,1 -k2,2 | awk 'NR == 1 {f = $3; sub(/\.log$/, ".umc180.txt", f); print f}' ;;
  esac
}
run_one() {   # $1 = output prefix, rest = dbs arguments; skipped when already finished
  local f=$1; shift
  [ -f "$f.log" ] && grep -q BEST-stm130 "$f.log" && return 0
  local t0=$SECONDS r pre=()
  [ -n "${TIMEOUT:-}" ] && pre=(timeout "$TIMEOUT")   # time limit per run (invmc)
  if [ -x /usr/bin/time ]; then
    /usr/bin/time -f "# cpu %U user %S sys, wall %e s" -o "$f.time" ${pre[@]+"${pre[@]}"} "$ROOT/bin/dbs" "$@" --save "$f" > "$f.log" 2> "$f.err" || true
  else
    ${pre[@]+"${pre[@]}"} "$ROOT/bin/dbs" "$@" --save "$f" > "$f.log" 2> "$f.err" || true
  fi
  r=$(grep -m1 -oE 'BEST-xor2: [0-9]+' "$f.log" || true)
  echo "$(date '+%m-%d %H:%M') $(basename "$f"): ${r:-NO RESULT} ($(( (SECONDS - t0) / 60 )) min)"
}
export -f run_one; export ROOT
JOBS="$OUT/$EXP/jobs.txt"; : > "$JOBS"
case $EXP in
rank)
  for m in $(names); do
    for sd in 1 2 3 4 5 6 7 8; do
      echo "$OUT/rank/${m}__all_w2000_s1_rgates_seed$sd $ROOT/bench/$m.txt $(depth_of "$m") --tie all --beam 2000 --slack 1 --rank gates --seed $sd --trace" >> "$JOBS"
    done
  done ;;
area|stage4|stage4from2|stage4fromsfx|stage4long)
  for m in $(names); do
    D=$(depth_of "$m")
    [ "$EXP" = area ] || [ "$D" = 3 ] || continue
    if [ "$EXP" = stage4long ]; then
      case $m in C_BeiKraLea16_32_Circulant|ePrint_JeaPeySim_32_Subfield|FSE_SKOP15_32_Involutory|KLSW) ;; *) continue ;; esac
    fi
    case $EXP in
      area)        start=$(best_in umc180 "$m" "$PAPER_RUN/forward" "$PAPER_RUN/bpbs" "$PAPER_RUN/bbeam" "$PAPER_RUN/lns"); obj=umc ;;
      stage4from2) start=$(best_in xor2 "$m" "$PAPER_RUN/bpbs"); obj=xor ;;
      stage4fromsfx) start="${SFX_START:?set SFX_START (see the header)}/$m.sfx.txt"; obj=xor ;;
      *)           start=$(best_in xor2 "$m" "$PAPER_RUN/forward"); obj=xor ;;
    esac
    [ -n "$start" ] && [ -f "$start" ] || { echo "$EXP: no starting circuit for $m"; continue; }
    if [ "$(size_of "$m")" -le 16 ]; then W=1000; IT=300; else W=200; IT=150; fi
    [ "$EXP" = stage4long ] && IT=$((IT * 10))
    for sd in 1 2 3 4 5 6 7 8; do
      echo "$OUT/$EXP/$m.seed$sd $ROOT/bench/$m.txt $D --lns $start --obj $obj --iters $IT --fix 4:28 --beam $W --slack 1 --tie all --seed $sd --trace" >> "$JOBS"
    done
  done ;;
bpbs)
  for m in $(names); do
    [ "$(depth_of "$m")" = 3 ] && [ "$(size_of "$m")" -gt 16 ] || continue
    for sd in $(seq 9 23); do
      echo "$OUT/bpbs/$m.seed$sd $ROOT/bench/$m.txt 3 --bpbs 300 --maxsplits 12 --evals 4 --tie sfx --md --seed $sd --trace" >> "$JOBS"
    done
  done ;;
depth4)
  m=ToSC_SarSye16_32_Involutory
  for sd in 1 2 3 4 5 6 7 8; do
    echo "$OUT/depth4/$m.bpbs.seed$sd $ROOT/bench/$m.txt 4 --bpbs 300 --maxsplits 12 --evals 4 --tie sfx --md --seed $sd --trace" >> "$JOBS"
    echo "$OUT/depth4/$m.bbeam.seed$sd $ROOT/bench/$m.txt 4 --bpbs 6 --bbeam 8 --bopts 16 --evals 2 --maxsplits 14 --tie sfx --md --seed $sd --trace" >> "$JOBS"
  done ;;
beamstats)
  for m in $(names); do
    [ "$(size_of "$m")" -eq 16 ] || [ "$m" = FSE_LiWang16_32_Involutory ] || continue
    for sd in 1 2 3 4; do
      for W in 1 10 100 2000; do
        [ "$(size_of "$m")" -eq 16 ] || [ "$W" -lt 2000 ] || [ "$sd" -le 2 ] || continue
        echo "$OUT/beamstats/$m.sfxmd0.w$W.s$sd $ROOT/bench/$m.txt 3 --tie sfx --md --beam $W --slack 0 --seed $sd --stats" >> "$JOBS"
      done
      if [ "$(size_of "$m")" -eq 16 ]; then
        echo "$OUT/beamstats/$m.sfx0.w2000.s$sd $ROOT/bench/$m.txt 3 --tie sfx --beam 2000 --slack 0 --seed $sd --stats" >> "$JOBS"
        echo "$OUT/beamstats/$m.all0.w2000.s$sd $ROOT/bench/$m.txt 3 --tie all --beam 2000 --slack 0 --seed $sd --stats" >> "$JOBS"
      fi
      echo "$OUT/beamstats/$m.sfxmd1.w2000.s$sd $ROOT/bench/$m.txt 3 --tie sfx --md --beam 2000 --slack 1 --seed $sd --stats" >> "$JOBS"
      [ "$(size_of "$m")" -eq 16 ] || [ "$sd" -le 2 ] || continue
      echo "$OUT/beamstats/$m.sfxmd0rg.w2000.s$sd $ROOT/bench/$m.txt 3 --tie sfx --md --beam 2000 --slack 0 --rank gates --seed $sd --stats" >> "$JOBS"
    done
  done ;;
sfx1)
  for m in $(names); do
    for sd in 1 2 3 4 5 6 7 8; do
      echo "$OUT/sfx1/${m}__sfx_w2000_s1_seed$sd $ROOT/bench/$m.txt $(depth_of "$m") --tie sfx --beam 2000 --slack 1 --seed $sd --trace" >> "$JOBS"
    done
  done ;;
fwdtime)
  for D in 3 4 5; do
    for sd in 1 2 3; do
      echo "$OUT/fwdtime/AES.d$D.bp.seed$sd $ROOT/bench/AES.txt $D --tie bp --runs 1 --seed $sd --trace" >> "$JOBS"
    done
  done ;;
invmc)
  tail -n +3 "$ROOT/examples/aes_invmixcolumns.mat.txt" > "$OUT/invmc/AESinv.txt"   # bare rows for dbs
  export TIMEOUT=${TIMEOUT:-14h}
  for t in bp sfx all; do
    echo "$OUT/invmc/AESinv.d5.$t.seed1 $OUT/invmc/AESinv.txt 5 --tie $t --runs 1 --seed 1 --progress --trace" >> "$JOBS"
  done ;;
aesdepth)
  for D in 4 5 8; do
    echo "$OUT/aesdepth/AES.d$D.stage1 $ROOT/bench/AES.txt $D --tie all --beam 500 --slack 1 --seed 1 --trace" >> "$JOBS"
  done ;;
*) echo "unknown experiment $EXP"; exit 1 ;;
esac
echo "$EXP: $(wc -l < "$JOBS") runs, $J at a time, started $(date '+%m-%d %H:%M')"
if [ -n "${DRY:-}" ]; then cat "$JOBS"; exit 0; fi
xargs -P "$J" -L 1 bash -c 'run_one "$@"' _ < "$JOBS"
if [ "$EXP" = aesdepth ]; then   # stage 4 from the stage-1 circuit of each depth bound
  : > "$JOBS"
  for D in 4 5 8; do
    for sd in 1 2 3 4 5 6 7 8; do
      echo "$OUT/aesdepth/AES.d$D.lns.seed$sd $ROOT/bench/AES.txt $D --lns $OUT/aesdepth/AES.d$D.stage1.xor2.txt --iters 150 --fix 4:28 --beam 200 --slack 1 --tie all --seed $sd --trace" >> "$JOBS"
    done
  done
  xargs -P "$J" -L 1 bash -c 'run_one "$@"' _ < "$JOBS"
fi
if [ "$EXP" = depth4 ]; then   # stage 4 from the best circuit of stage 1 (paper's run) and stages 2-3 (above)
  m=ToSC_SarSye16_32_Involutory
  shopt -s nullglob
  cand=$(awk '/^BEST-xor2:/ {print $2, FILENAME}' "$PAPER_RUN/forward/$m"__*.log "$OUT/depth4/$m".*.log | sort -n -k1,1 | awk 'NR == 1 {f = $2; sub(/\.log$/, ".xor2.txt", f); print f}')
  shopt -u nullglob
  : > "$JOBS"
  for sd in 1 2 3 4 5 6 7 8; do
    echo "$OUT/depth4/$m.lns.seed$sd $ROOT/bench/$m.txt 4 --lns $cand --iters 150 --fix 4:28 --beam 200 --slack 1 --tie all --seed $sd --trace" >> "$JOBS"
  done
  echo "depth4: stage 4 from $cand"
  xargs -P "$J" -L 1 bash -c 'run_one "$@"' _ < "$JOBS"
fi
echo "$EXP: finished $(date '+%m-%d %H:%M')"
