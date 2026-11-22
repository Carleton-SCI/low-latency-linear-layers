# Best circuits of the paper

One file per matrix and measure, in the circuit format of the NIST Circuit Complexity Project: one gate per line,
`XOR out a b`, with the inputs `x1, x2, ...` and the outputs `y1, y2, ...` numbered from 1 (`x1` is column 0 of the
matrix in `bench/`), and one extension, `XOR3 out a b c`, a 3-input XOR gate whose last two operands form the XOR2
gate that it absorbed. The comment lines give the matrix as integers (as in the NIST files), the gate counts, the
depth (gate levels; an XOR3 gate counts as one level) and the areas:

- `NAME.xor2.circ.txt`: the fewest XOR2 gates at the minimum depth.
- `NAME.umc180.circ.txt`: the smallest area after XOR3 merging, with the UMC 180 nm unit areas (maximum
  merge of the gates of fan-out one).
- `NAME.stm130.circ.txt`: the smallest area with the STM 130 nm unit areas, where the merge may also duplicate a
  gate of fan-out two and absorb one copy into each of its two users (saves 0.67 GE with these unit areas). This
  merge is exact (0-1 program, CP-SAT). The header of each file gives the XOR2 circuit before merging and the
  duplications.

The values are the best over all the runs of the paper at the minimum depth, i.e., over the 9,060 circuits
recorded by its 3,020 runs: the four stages (`run.tar.xz`, and `server_runs.tar.xz`, `runs/depth4`, for stages
2-4 of [SS16] involutory 32x32 at depth 4), stage 4 started from the published circuits of SFX23 (marked ‡ below,
as in Table 3 of the paper), stage 4 with the UMC area as the objective (marked §), and the runs of Secs. 5.3,
5.4 and 6 (`server_runs.tar.xz`). The `Source` line of each file names the run. Ties go to the four stages.

All the files are written and re-verified by `paper/tables/final_bests.py` (needs `pip install ortools`), which
also writes `summary.json`. `make verify` (`paper/checks/verify_results.py`) checks all of them against the
matrices of `bench/` and against `summary.json`; for AES MixColumns, `python3 paper/checks/verify_gf_aes.py
results/AES.xor2.circ.txt` also checks the circuit against MixColumns computed in GF(2^8) (bit order:
`bench/README.md`).

Other files:

- `run.tar.xz`: the logs of all the runs of the four stages and the circuits that they recorded.
- `server_runs.tar.xz`: the server runs of Oct 30 - Nov 1, 2026 (`paper/server_runs.sh`): stages 2-4 for [SS16]
  involutory 32x32 (`depth4`), stage 4 started from the circuits of SFX23 (`stage4fromsfx`), stage 4 with the UMC
  objective (`area`), the `sfx` rule at slack 1 (`sfx1`), the ranking by gates so far + Delta (`rank`), the
  diagnostics of the beam search (`beamstats`, `dbs --stats`), AES MixColumns under the depth bounds 4, 5 and 8
  (`aesdepth`) and AES InvMixColumns at depth 5 (`invmc`), with their logs, circuits, summaries
  (`extra_summary.md`, `beam_stats.txt`, `anytime_*.txt`) and the description of the server (`system.txt`, with the
  SHA-256 hashes of `src/dbs.cpp`, `paper/extra_runs.sh` and `paper/server_runs.sh`, which are unchanged since).
- `restarts/`: the raw data of Table 6 (`results.jsonl`, `config.json`); `paper/ties_vs_restarts.py
  summary` on it gives Table 6 exactly (see its README).

Areas are sums of unit areas in GE. The published values are from Table 5 of Shi–Feng–Xu and from the HILL
repository (see `bench/README.md`). "XOR2 (ours)" is the best of our runs on each matrix; Table 3 of the paper
also transfers circuits between permutation-equivalent matrices (Sec. 4 of the paper): for example, LW16
involutory 16 gets 46 gates through Joltik.
