#!/usr/bin/env bash
# The server runs for the paper, back to back, into one folder, with one archive to send back.
#
#   ./paper/server_runs.sh                 (from the root of the repository; nothing else to set)
#
# It runs, with paper/extra_runs.sh, the experiments that the text and the tables of the paper use:
#   beamstats      Sec. 5.4: the diagnostics of the beam search (near-ties, widths, s = 1 and the ranking
#                  with the sfx rule and the MD filter); about 4 core-hours
#   stage4fromsfx  Table 1 and Sec. 5.3: stage 4 from the published circuits of SFX23 (paper/sfx_start/),
#                  all 22 matrices of depth 3; about 150 core-hours
#   sfx1           Sec. 5.3 (Table 3): the sfx rule at slack 1, all 23 matrices; about 55 core-hours
#   rank           Sec. 6: the ranking of the beam search by gates so far + Delta, all 23 matrices; about
#                  55 core-hours
#   aesdepth       Sec. 6: AES MixColumns under the depth bounds 4, 5 and 8 (stage 1, then stage 4);
#                  about 40 core-hours
#   area           Table 1: stage 4 with the UMC area as the objective, all 23 matrices; about 175
#                  core-hours
#   depth4         Table 1: stages 2-4 for [SS16] involutory 32x32 at its minimum depth 4; the longest
#                  (1-3 days), so it runs in the background from the start, on DEPTH4_JOBS threads
# The paper's benchmark run, which several experiments start from or compare with, is taken from
# results/run.tar.xz (or from PAPER_RUN, if set).
#
# Options (environment variables):
#   OUT=DIR            output folder (default ./server_results); the archive is DIR.tar.gz
#   JOBS=N             threads in total (default: nproc)
#   DEPTH4_JOBS=N      threads of depth4 (default 16; 0 skips it)
#   INVMC=1            also AES InvMixColumns at depth 5, each run limited to INVMC_TIMEOUT (default 14h),
#                      on 3 threads in the background (not needed by the paper)
#   PAPER_RUN=DIR      the run/ folder of the paper's benchmark run, instead of results/run.tar.xz
#   ONLY=A,B           restrict the matrices (for a quick test), as in extra_runs.sh
# Arguments: a subset of the experiments above, e.g. ./paper/server_runs.sh beamstats rank
#
# Run it in screen or tmux; it can be stopped and started again with the same command: finished runs are
# kept, and a second copy on the same folder stops at once. The archive is refreshed after every
# experiment, so a partial archive can be sent at any time. Send back DIR.tar.gz.
set -euo pipefail
unset TIMEOUT                                   # no time limit for the runs of the paper (see INVMC_TIMEOUT)
ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT=${OUT:-$PWD/server_results}
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)
JOBS=${JOBS:-$(nproc)}
DEPTH4_JOBS=${DEPTH4_JOBS:-16}
ALL="beamstats stage4fromsfx sfx1 rank aesdepth area depth4"
EXPS=${*:-$ALL}
for e in $EXPS; do case " $ALL " in *" $e "*) ;; *) echo "unknown experiment $e (known: $ALL)"; exit 1 ;; esac; done

if command -v flock > /dev/null; then          # one copy per output folder
  exec 9> "$OUT/.lock"
  flock -n 9 || { echo "another server_runs.sh is using $OUT; stop it first"; exit 1; }
fi
LOG=$OUT/server_runs.log
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$LOG"; }
log "start: experiments: $EXPS; JOBS=$JOBS DEPTH4_JOBS=$DEPTH4_JOBS INVMC=${INVMC:-0} OUT=$OUT"

# the program, its tests, and the machine
make -s -C "$ROOT" bin/dbs
if ! make -s -C "$ROOT" test > "$OUT/make_test.txt" 2>&1; then log "make test failed (see $OUT/make_test.txt)"; exit 1; fi
{
  echo "date: $(date)"; echo "host: $(hostname)"; uname -a; echo "nproc: $(nproc)"
  command -v lscpu > /dev/null && lscpu | grep -E "^(Model name|Socket|Core|Thread|CPU\(s\))"
  command -v free > /dev/null && free -g | head -2
  g++ --version | head -1; python3 --version
  sha256sum "$ROOT/src/dbs.cpp" "$ROOT/paper/extra_runs.sh" "$0"
} > "$OUT/system.txt" 2>&1 || true

# the paper's benchmark run
if [ -z "${PAPER_RUN:-}" ]; then
  if [ ! -d "$OUT/paper_run/run/forward" ]; then
    mkdir -p "$OUT/paper_run"
    tar -xJf "$ROOT/results/run.tar.xz" -C "$OUT/paper_run"
  fi
  PAPER_RUN=$OUT/paper_run/run
fi
for d in forward bpbs bbeam lns; do
  [ -d "$PAPER_RUN/$d" ] || { log "PAPER_RUN=$PAPER_RUN has no $d/"; exit 1; }
done
export PAPER_RUN SFX_START=${SFX_START:-$ROOT/paper/sfx_start}
RUNS=$OUT/runs
mkdir -p "$RUNS"

experiment() {   # experiment NAME THREADS: one experiment of extra_runs.sh, its output in $OUT/runs_NAME.txt
  log "start $1 ($2 threads)"
  if OUT="$RUNS" "$ROOT/paper/extra_runs.sh" "$1" "$2" >> "$OUT/runs_$1.txt" 2>&1; then
    log "end $1: $(grep -c ' (' "$OUT/runs_$1.txt" || true) runs listed in runs_$1.txt"
  else
    log "$1: extra_runs.sh stopped with an error (see runs_$1.txt)"
  fi
}
summarize() {    # the summaries the paper uses; each one covers whatever has finished so far
  python3 "$ROOT/paper/extra_summary.py" "$RUNS" "$PAPER_RUN" > "$OUT/extra_summary.md" 2> "$OUT/extra_summary.err" || true
  if [ -d "$RUNS/beamstats" ]; then
    python3 "$ROOT/paper/tables/beam_stats.py" "$RUNS/beamstats" > "$OUT/beam_stats.txt" 2>&1 || true
  fi
  if [ -d "$RUNS/stage4fromsfx" ]; then
    python3 "$ROOT/paper/tables/extra_stage4.py" "$RUNS/stage4fromsfx" "$OUT/extra_stage4.json" > "$OUT/extra_stage4.txt" 2>&1 || true
  fi
  for e in stage4fromsfx sfx1 rank aesdepth area depth4 invmc; do
    if [ -d "$RUNS/$e" ]; then python3 "$ROOT/paper/tables/anytime.py" "$RUNS/$e" > "$OUT/anytime_$e.txt" 2>&1 || true; fi
  done
}
pack() {         # the archive to send back (everything but the extracted paper run)
  local tmp="$OUT.tar.gz.part"
  tar -czf "$tmp" --warning=no-file-changed --exclude="$(basename "$OUT")/paper_run" --exclude=.lock \
      -C "$(dirname "$OUT")" "$(basename "$OUT")" 2> /dev/null || true
  if [ -s "$tmp" ]; then
    mv -f "$tmp" "$OUT.tar.gz"
    log "archive: $OUT.tar.gz ($(du -h "$OUT.tar.gz" | cut -f1))"
  else
    log "archive: tar failed"
  fi
}

# background: depth4 (and InvMixColumns if asked); the foreground uses the other threads
bg=()
FG=$JOBS
if [[ " $EXPS " == *" depth4 "* ]] && [ "$DEPTH4_JOBS" -gt 0 ]; then
  experiment depth4 "$DEPTH4_JOBS" &
  bg+=($!)
  FG=$((FG - DEPTH4_JOBS))
fi
if [ "${INVMC:-0}" = 1 ]; then
  log "start invmc (3 threads, ${INVMC_TIMEOUT:-14h} per run)"
  ( TIMEOUT=${INVMC_TIMEOUT:-14h} OUT="$RUNS" "$ROOT/paper/extra_runs.sh" invmc 3 >> "$OUT/runs_invmc.txt" 2>&1; log "end invmc" ) &
  bg+=($!)
  FG=$((FG - 3))
fi
[ "$FG" -ge 1 ] || FG=1

for e in $EXPS; do
  [ "$e" = depth4 ] && continue
  experiment "$e" "$FG"
  summarize
  pack
done
if [ ${#bg[@]} -gt 0 ]; then
  log "waiting for the background experiments (depth4: 1-3 days); the archive is refreshed every hour"
  n=0
  while :; do
    alive=0
    for p in "${bg[@]}"; do
      if kill -0 "$p" 2> /dev/null; then alive=1; fi
    done
    [ "$alive" = 1 ] || break
    sleep 60
    n=$((n + 1))
    if [ $((n % 60)) -eq 0 ]; then summarize; pack; fi
  done
  wait || true
fi
summarize
pack
log "done. Send $OUT.tar.gz"
