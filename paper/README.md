# Reproducing the paper

All commands run from the repository root.

## Runs

```bash
./paper/run_benchmark.sh 48                   # stages 1-4 for the 24 matrices of bench/, then the summary
STEPS="2 3" ONLY=AES ./paper/run_benchmark.sh 8   # some stages, some matrices
HILL=path/to/hill_repo SFX=path/to/sfx_repo STEPS=5 ./paper/run_benchmark.sh   # summary with the cross-checks
python3 paper/tables/final_bests.py RUN_PARENT server_runs OUT results   # results/ of the paper: all runs
```

| Step | What runs | Script |
|---|---|---|
| 1 | 56 forward runs per matrix: tie rules `all`, `bp` and `sfx`, widths 2000 and 8000, slack 0–2 | `bench_run.py` (grid `full`) |
| 2 | backward splitting (BPBS), 8 seeds | `run_benchmark.sh` |
| 3 | beam search over the backward splits, 8 seeds | `run_benchmark.sh` |
| 4 | large neighborhood search from the best circuit of steps 1–3, 8 seeds | `run_benchmark.sh` |
| 5 | summary: best circuits, independent re-verification, the HILL and SFX23 cross-checks | `summarize.py` |

The server used two Intel Xeon E5-2670 v3 processors (48 hardware threads) and 128 GB of RAM. The
software was g++ 11.4.0 (`-O3 -march=native -std=c++17`) on Ubuntu 20.04.6 LTS. Up to 48 runs ran at
a time.

## Where each number of the paper comes from

| Paper | Script |
|---|---|
| Algorithm 1 | `src/dbs.cpp`: `beam_search` (stage 1, and in stage 4), `draw_split` (priority (3) of SFX23), and the `--bpbs`, `--bbeam` and `--lns` branches of `main` (stages 2, 3 and 4) |
| Sec. 2.2, the binary tree method gives the minimum depth | `tables/check_tree.py` |
| Sec. 2.4, the exact merge is optimal for every circuit | printed by `tables/final_bests.py` (and by `tables/stm_merge.py`, which also prints the time per circuit) |
| Sec. 3, the worked example on the 8 × 8 matrix (the matrix, Tables 1 and 2, the two circuit listings) | `tables/worked_example.py OUT_DIR` (computes the distances of Tables 1 and 2 with the forward search of `example/toy_engine.py`, reads the facts of the runs and the two circuits from `example/toy_data.json`, which `example/toy_data.py` writes (needs `numpy` and `ortools`; about 3 min), and asserts every fact of the example) |
| Sec. 4: the circuits of SFX23, equivalences, orientation | `summarize.py ... --sfx SFX_DIR` (needs `networkx`) |
| Sec. 4: the matrices of the circuits of the HILL repository (including its [LS16] 32x32 variant), their counts, areas and depths | `compare/hill_to_matrices.py HILL_REPO OUT_DIR` |
| Sec. 4, HILL areas: for 54 of the 56 matrices, the area equals a maximum merge of the stored XOR2 circuit | `compare/hill_merge_check.py HILL_REPO` |
| Sec. 4: every benchmark matrix is invertible and MDS (69 square block submatrices) | `compare/mds_check.py bench` |
| Table 3 (XOR2 counts and areas), Table 4 (each stage), Table 5 (forward beam search by configuration) | Extract `results/run.tar.xz` into `RUN_PARENT` (the four stages) and `results/server_runs.tar.xz` into a folder (`server_runs/`, the other runs). Then `tables/ablation.py RUN_PARENT bench OUT_DIR server_runs/runs/depth4` (the four stages, with stages 2-4 of [SS16] involutory 32x32 at depth 4; writes `OUT_DIR/ablation.json`), `tables/final_bests.py RUN_PARENT server_runs OUT_DIR [results]` (all the 9,060 circuits at the minimum depth, their UMC and exact STM areas, the best values of each group of runs; needs `pip install ortools`; writes `OUT_DIR/bests.json`, and with the last argument the files of `results/`), `compare/sfx_merge.py SFX_RESULTS_DIR bench OUT_DIR/sfx_merge.json` (the SFX23 columns: the published circuits of SFX23 merged as ours; needs `networkx` and `ortools`), then `tables/make_tables.py OUT_DIR/ablation.json bench/benchmark.tsv OUT_DIR/bests.json OUT_DIR/sfx_merge.json OUT_DIR`. |
| Table 3 (marks ‡ and §) and Sec. 5.3, the additional runs of stage 4 | `SFX_START=paper/sfx_start ./paper/extra_runs.sh stage4fromsfx` (stage 4 from the published circuits of SFX23; `paper/sfx_start/` was written by `compare/sfx_to_dbs.py SFX_RESULTS_DIR bench paper/sfx_start`) and `PAPER_RUN=RUN_PARENT/run ./paper/extra_runs.sh area` (stage 4 with the UMC objective), or `./paper/server_runs.sh`; the runs are in `results/server_runs.tar.xz` (`runs/stage4fromsfx`, `runs/area`), summarized by `tables/extra_stage4.py` and `extra_summary.py` |
| Sec. 5.2: the gains (1 to 4 gates: 1 for 2 problems, 2 for 4, 3 for 3, 4 for 5), the area range (0.3–7.3 %), the area comparison by problem (smaller for 19 of 20; fewer XOR2 gates for 11, as many or more for 8; with the circuits of fewest XOR2 gates only, below / above for 16 / 4), count against area (16 of 23 matrices, 0.31–9.96 GE), the STM duplication (15 matrices, 0.67–5.36 GE), and the counts of Sec. 5.3 | printed by `tables/make_tables.py` and saved in its `facts.json` (over all the runs at the minimum depth) |
| Sec. 5.3, forward runs of stages 2 and 3 (1.27 and 1.33 million on the 32x32 matrices of depth 3) | `tables/forward_runs.py RUN_DIR` |
| Sec. 5.3 and Table 5, the `sfx` rule at slack 1; Table 3, the area columns of SFX23 and HILL per problem | `extra_runs.sh sfx1` (all 23 matrices; `runs/sfx1` of `results/server_runs.tar.xz`, read by `tables/final_bests.py` into `bests.json`); `compare/sfx_merge.py` (its `same_matrix` field gives the daggers) and `tables/make_tables.py` |
| Table 6 (beam search against random restarts) and Sec. 5.4 | `ties_vs_restarts.py run -o tvr`, then `ties_vs_restarts.py summary tvr > tvr/summary.md` and `tables/restarts_table.py tvr/summary.md restarts.tex` (see below); the raw data of the paper's run is `results/restarts/` (its summary gives Table 6 exactly) |
| Sec. 5.4, the parent sequences of the near-ties per iteration and the widths W = 1, 10, 100, 2000 with `sfx` + MD | `src/dbs.cpp --stats` (read-only diagnostics on stderr), `extra_runs.sh beamstats`, then `tables/beam_stats.py server_runs/runs/beamstats` |
| Table 7 (run times) and Sec. 5.5 | stage 1: the `# <s> s` line of each forward log (`tables/ablation.py`); stages 2–4: `tables/stage_times.py RUN_DIR`; stages 2-4 of [SS16] involutory 32x32: the `.time` files of `runs/depth4` in `results/server_runs.tar.xz`; Sec. 5.5, the CPU time of the additional runs of stage 4: `extra_summary.py` |
| Sec. 6, ranking by gates so far + Delta in stage 1 (all 23 matrices, seeds 1-8) | `extra_runs.sh rank` (`runs/rank` of `results/server_runs.tar.xz`), summarized by `extra_summary.py` |
| Sec. 6, AES MixColumns under the depth bounds 4, 5 and 8 | `extra_runs.sh aesdepth` (stage 1 under each depth bound, then stage 4; `runs/aesdepth` of `results/server_runs.tar.xz`) |
| Sec. 6, AES InvMixColumns at depth 5 | `extra_runs.sh invmc` (`dbs --progress`; time limit `TIMEOUT`, default 14 h; `runs/invmc` of `results/server_runs.tar.xz`), summarized by `extra_summary.py`; the matrix is `examples/aes_invmixcolumns.mat.txt`, the matrix of NIST's 119-gate circuit |
| Appendix A (the record circuits and the circuits of our UMC areas, with their matrices) | `tables/appendix_circuits.py OUT.tex` (reads `bench/` and `results/`, checks every circuit, and the XOR2 depth of each merged circuit) |
| Sec. 5.1, the circuits of `results/` | `checks/verify_results.py` (`make verify`): each of the 72 circuits computes its matrix within the minimum depth, with the gate counts and areas of `results/summary.json` |
| AES MixColumns circuits checked over GF(2^8) | `checks/verify_gf_aes.py results/AES.xor2.circ.txt` (also the native format of `dbs`), `checks/verify_mc.py` |

The tables that these scripts write are those of the paper's LaTeX source, between its `% begin tables/X.tex`
and `% end tables/X.tex` markers.

`RUN_PARENT` is the folder that contains `run/`.

**Run scripts as they ran.** `extra_runs.sh` and `server_runs.sh` are the versions of the server runs of Oct 30 -
Nov 1, 2026, whose SHA-256 hashes are in `system.txt` of `results/server_runs.tar.xz`, as is `src/dbs.cpp`. Their
comments, and the section titles of `extra_summary.py` and of its output in that archive, number the tables as
the draft of that time: its Tables 1 and 3 are Tables 3 and 5 of the paper. `paper/audit/` in these comments is
now `paper/compare/`. The paper does not use `fwdtime`, nor the `beamstats` runs at s = 1 and with the ranking by
gates (`sfxmd1`, `sfxmd0rg`). Delta, the sum of the distances in the paper, is `Sigma` in `src/dbs.cpp` (`--rank
sigma|gates`) and in these summaries.

## Tie exploration against random restarts, at equal CPU time

`ties_vs_restarts.py` compares the beam search over the ties (stage 1) with independent greedy runs
that break the ties at random, the approach of SFX23. Both sides get the same CPU time and use the
same program, distances, tie rule and depth bound, with the forward search only. The header of the
script describes the design and the three default pairs (`all0`, `all1`, `sfxmd0`).

```bash
make
python3 paper/ties_vs_restarts.py run --only JOLTIK,AES --width 200 --seeds 1 -o tvr_smoke   # smoke test, about 10 min
python3 paper/ties_vs_restarts.py run -o tvr                                                # 23 matrices x 3 pairs x 4 seeds
python3 paper/ties_vs_restarts.py summary tvr --paper-run RUN_DIR > tvr/summary.md
python3 paper/tables/restarts_table.py tvr/summary.md restarts.tex                          # Table 6 of the paper
```

`run` stays in the foreground and prints a line for every finished task and a status line every 10
minutes: start it in `screen` or `tmux`, not with `nohup ... &`. Only one `run` can use an output
folder at a time; a second one stops with a message.

- **Cost** (rough estimate from the stage-1 times of the paper's run-time table): about 4 h on 24 cores without
  [SS16] involutory 32×32 (`--skip ToSC_SarSye16_32_Involutory`), and about 8 h with it, because
  each of its `sfxmd0` tasks takes about that long. `--pairs all0,all1` drops the slowest pair.
- **Jobs.** The default is one job per physical core, so that both sides of a task run under the
  same load.
- **Resuming.** The run can be restarted at any time: finished tasks are kept in
  `tvr/results.jsonl`.
- **Summary.** It uses one result per (matrix, pair, seed), the first one, and warns if the results
  file has repeated records (`--keep-duplicates` uses them all). For each pair and matrix, it prints:
  - the XOR2 counts of both sides for every seed, and how often the BFS is lower, equal or higher;
  - the number of restarts that fit in the BFS time;
  - from the distribution of all the restarts, the expected CPU time they need to reach the median
    BFS result, in multiples of the BFS time.
- **Check.** `--paper-run` checks that the BFS runs of `all1` reproduce the paper's stage-1 runs of
  the same seeds.

## Additional experiments

`server_runs.sh` runs, back to back and with no setting, the experiments that the paper uses: `beamstats`,
`stage4fromsfx` (from `paper/sfx_start/`), `sfx1`, `rank`, `aesdepth`, `area` and, in the background,
`depth4`. It reads the paper's run from `results/run.tar.xz`, writes everything into one folder
(`server_results/`), and refreshes one archive (`server_results.tar.gz`) after each experiment, with the
summaries (`extra_summary.md`, `beam_stats.txt`, `extra_stage4.json`, `anytime_*.txt`). It can be stopped
and started again. The archive of the paper's server runs (Oct 30 - Nov 1, 2026, on the server of Sec. 5.1;
the server path replaced by `<repo>`) is `results/server_runs.tar.xz`.

`extra_runs.sh` runs them (ranking of the beam search by gates + Σ, the area as the objective of the
neighborhood search, stage 4 from stage 1 only, from stage 2 only or from the circuits of SFX23, stage 4
with 10 times the iterations, stage 2 with a larger budget, stages 2-4 at depth 4, AES MixColumns at the
depth bounds 4, 5 and 8, the time of one forward run for AES MixColumns, AES InvMixColumns at depth 5, and
the diagnostics of the beam search), and `extra_summary.py` compares their
results with a run of the paper. Some of these runs also give numbers that the final version of the paper
does not report (for example, `beamstats` with `sfxmd1` and `sfxmd0rg`); their summaries are in
`results/server_runs.tar.xz`. Their runs use `dbs --trace`, and `tables/anytime.py RUN_DIR` gives their
anytime curves (best XOR2 count or UMC area against CPU time). `ties_vs_restarts.py` also has the pairs
`xall1` (beam with `all` against restarts with `sfx` + MD), `sfxmd1` (`sfx` + MD on both sides, slack 1) and `sfx0`.

## External repositories

The comparisons read two published repositories:

- the HILL repository: <https://github.com/13est-Wishes/Linear_Layers_by_YongchaoLi>;
- the result files of Shi, Feng and Xu, from the code repository cited in their ToSC 2023(4)
  paper.

Neither is redistributed here.
