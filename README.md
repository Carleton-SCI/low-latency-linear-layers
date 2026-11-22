# Ties, splits and neighborhoods: new area records for low-latency linear layers

This repository holds the code, data and results of the paper *Ties, Splits and Neighborhoods: New
Area Records for Low-Latency Linear Layers* (submitted to IACR Transactions on Symmetric Cryptology).

**The problem.** A linear layer is a binary matrix `M`, and it is implemented as a circuit of XOR
gates that computes `y = M x`. We look for the fewest 2-input XOR (XOR2) gates at the minimum
depth, which is the low-latency setting. We also look for the smallest area when 3-input XOR
(XOR3) gates are available.

**The search** combines known techniques in four stages, all in `src/dbs.cpp`:

1. Beam search (width `W`) over the near-ties of a depth-bounded Boyar–Peralta forward search with exact
   distances, i.e., the candidate gates within a slack `s` of the best one in the beam.
2. The backward splitting of Shi, Feng and Xu (our reimplementation of their BPBS).
3. Beam search over these backward splits.
4. A large neighborhood search: keep the last gate of `k` random outputs of the best circuit and
   re-solve the rest. It can also start from another circuit (for example, a published one) and
   minimize the area after merging instead of the XOR2 count.

Every circuit that is found is also merged into XOR2/XOR3 gates, and the circuit with the smallest
area is kept. For the STM 130 nm unit areas, an exact merge that may also duplicate gates of fan-out
two is applied afterwards to all the recorded circuits (`paper/tables/final_bests.py`, with
`paper/tables/stm_merge.py`).

**Results.** On the 23 benchmark matrices of Shi, Feng and Xu (20 distinct problems), we obtain:

- new best XOR2 counts at the minimum depth for 14 problems, by 1 to 4 gates (12 from the four stages;
  for SmallScale AES and [SKOP15] Hadamard 32x32, from stage 4 started from the published circuits of
  Shi, Feng and Xu);
- equal counts for the other 6 problems;
- smaller XOR2/XOR3 areas (UMC 180 nm unit areas) than both the circuits of the HILL repository and the
  circuits of Shi, Feng and Xu merged in the same way, for 19 of the 20 problems, and an equal area for the
  other one.

See `results/`, which also holds the logs and all the circuits recorded by the runs (`results/run.tar.xz`
for the four stages, `results/server_runs.tar.xz` for the other runs of the paper).

## Contents

| Path | Content |
|---|---|
| `slpd.py`, `slpd/` | Command-line tool: checks a matrix file, runs the search, verifies, merges and converts circuits; its library also serves the scripts of `paper/` and shares no code with `src/dbs.cpp` |
| `web/` | The search in the browser: `web/index.html` runs the presets of the tool on a matrix of your choice (see `web/README.md`) |
| `src/dbs.cpp` | The search program (C++17, single-threaded, deterministic for a given seed) |
| `bench/` | The 24 matrices of the paper (bare rows) and `benchmark.tsv` with the published values |
| `results/` | The best circuits of the paper for every matrix, a summary, and the logs and circuits of all the runs (see `results/README.md`) |
| `paper/` | Scripts that run the experiments, reproduce the tables and check the circuits; `paper/README.md` maps the paper to them |
| `examples/` | A 4 × 6 toy matrix, AES MixColumns, AES InvMixColumns and Joltik, in the NIST format |
| `docs/FORMATS.md` | Matrix and circuit formats |
| `tests/` | Tests of `slpd/` and of the search program |

## Requirements

- A C++17 compiler (for example g++ 9 or later) and `make`. The default flags are
  `-O3 -march=native -std=c++17`; with a compiler that does not accept `-march=native`, build with
  `make CXXFLAGS="-O3 -std=c++17"`.
- Python 3.8 or later. The tool and the checks need only the standard library. Some scripts of `paper/` also need
  `ortools` (the exact merge for the STM areas), `networkx` (the comparison with the circuits of Shi–Feng–Xu)
  or `numpy` (the data of the example of Sec. 3); `paper/README.md` says which.

```bash
make            # builds bin/dbs (the tool also builds it on first use)
make test       # tests of slpd/, of the search program and of web/
make verify     # checks the 72 circuits of results/ (paper/checks/verify_results.py)
```

`make verify` checks that every circuit of `results/` computes its matrix of `bench/` within the minimum depth,
with the gate counts and areas of `results/summary.json`. For AES MixColumns,
`python3 paper/checks/verify_gf_aes.py results/AES.xor2.circ.txt` also checks the circuit against MixColumns
computed in GF(2^8).

## Quick start

```bash
python3 slpd.py check examples/toy.mat.txt
python3 slpd.py run   examples/toy.mat.txt -o out
python3 slpd.py verify examples/toy.mat.txt out/toy.umc180.circ.txt -D 2
```

`run` writes three files per matrix:

- `out/NAME.xor2.circ.txt`: the fewest XOR2 gates;
- `out/NAME.umc180.circ.txt` and `out/NAME.stm130.circ.txt`: the smallest area after merging, for
  the UMC 180 nm and STM 130 nm unit areas;
- `out/NAME.summary.json`: the run summary.

Every circuit is re-verified in Python, independently of the search program, before it is written.
The logs and native circuits of every run stay in `out/work/NAME/`.

### Your own matrix

1. Write the matrix in the NIST / Boyar–Peralta format: a line with the number of matrices, a line
   with the numbers of rows and columns, then the rows of `0` and `1`. Bare rows also work. See
   `docs/FORMATS.md`.
2. Check the file:

   ```bash
   python3 slpd.py check my_matrix.mat.txt
   ```

   The check reports, with line and column numbers:
   - entries that are not `0` or `1`, rows of different lengths, and headers that do not match
     the data;
   - all-zero rows;
   - sizes over the limits (64 inputs, 128 outputs).

   It also warns about duplicate rows, rows of weight 1 and unused inputs, and prints the minimum
   depth.
3. Run the search:

   ```bash
   python3 slpd.py run my_matrix.mat.txt --preset standard -j 8 -o out
   ```

**Options of `run`.**

- `-D d`: a depth bound above the minimum. Stages 2–4 assume the minimum depth, as in Shi–Feng–Xu,
  so with a larger bound only stage 1 runs.
- `--stages 1,4`: run only some of the stages.
- `--timeout S`: drop any single run that takes longer than `S` seconds.

**Presets.** `paper` is the exact configuration of the paper: its 16×16 settings for at most 16
inputs, and its 32×32 settings otherwise.

| Preset | Stage 1 | Stages 2–4 |
|---|---|---|
| `quick` | 8 runs of width 500 (32×32: 4 runs of width 200) | 2 seeds each, small budgets |
| `standard` | 14 runs of width 2000 (32×32: 10 runs) | 4 seeds each, larger budgets |
| `paper` | 56 runs, widths 2000 and 8000 | 8 seeds each |

`quick` and `standard` are smaller budgets for trying the tool (`slpd/pipeline.py`); they do not reproduce
the paper. With `paper`, the server of the paper spent about 8.5 core-hours per 16×16 matrix and about 150
core-hours per 32×32 matrix of depth 3 (Table 7 of the paper). For Joltik, these four stages reached
46 XOR2 gates and 124.05 GE, and for AES MixColumns, 97 XOR2 gates and 260.10 GE (UMC 180 nm areas after
merging; Table 4 of the paper).

### In the browser

`web/index.html` runs the same presets in the browser, on a pasted or opened matrix or on an example, and gives the
same files. Download it and open it: it needs no installation and sends nothing over the network (see
`web/README.md`).

## Reproducing the paper

The 1,896 runs of the four stages used `paper/run_benchmark.sh` on a server with two Intel Xeon
E5-2670 v3 processors and 48 hardware threads. The software was g++ 11.4.0 with the flags
`-O3 -march=native -std=c++17`, on Ubuntu 20.04.6 LTS.

```bash
./paper/run_benchmark.sh 48                       # stages 1-4 and the summary, results in ./run
HILL=path/to/hill_repo SFX=path/to/sfx_repo STEPS=5 ./paper/run_benchmark.sh
```

The second command adds the cross-checks against the published circuits of HILL and Shi–Feng–Xu.
`paper/README.md` maps the tables and the values of the paper to the scripts that produce them.

## License

MIT; see `LICENSE`. The matrices in `bench/` come from published papers and from the public
repositories of Shi–Feng–Xu and HILL; see `bench/README.md`.

## Use of AI tools

See the section "Use of AI Tools" of the paper.

## Citation

To be added after publication.

## Authors

Mostafa Taha (Carleton University) and Arash Reyhani-Masoleh (Western University).
