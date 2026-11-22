# Starting circuits from the published circuits of SFX23

These files are the starting circuits of `paper/extra_runs.sh stage4fromsfx` (Sec. 5.3 of the paper).
Each one is the published circuit of Shi, Feng and Xu (ToSC 2023(4)) with the fewest XOR2 gates at
depth 3 for one benchmark matrix (then the smallest UMC area after a maximum merge). It is mapped onto
the benchmark matrix by the row and column permutations, verified against the matrix, and written in
the native format of `src/dbs.cpp` (read by `--lns`).

- **Source:** the `our_results/` folder of the SFX23 code repository,
  <https://gitee.com/Haotian-Shi/A_framework_with_improved_heuristics_to_optimize_low-latency_implementations_of_linear_layers>.
- **Regenerate:** `python3 paper/compare/sfx_to_dbs.py SFX_RESULTS_DIR bench OUT_DIR` (needs `networkx`).
- **Files:** 22, one per benchmark matrix of minimum depth 3. The [LS16] 32×32 matrix of the HILL
  repository has none (Sec. 4 of the paper).
- **Identity:** they are the starting circuits of the runs of stage 4 from SFX23 in
  `results/server_runs.tar.xz` (`runs/stage4fromsfx`), and byte for byte those of the earlier runs on the
  16×16 matrices, whose results the server reproduced.
