# The search in the browser

`index.html` runs the four-stage search on a matrix of your choice, in the browser. Download the file (or open it
from a clone of the repository) and open it in a recent browser; it uses only standard features (WebAssembly, Web
Workers, BigInt) and was tested with Chromium. It needs no server and no installation, and it sends nothing over the
network: the search program, compiled to WebAssembly, runs on the page, one Web Worker per parallel run.

The page:

- reads a matrix (pasted, from a file, or an example: the toy matrix, the example of Sec. 3 of the paper, the 24
  matrices of `bench/` and AES InvMixColumns) and checks it with the checks of `slpd.py check`;
- runs a preset of `slpd.py run` (`quick`, `standard` or `paper`), with a depth bound, a choice of stages and a
  number of parallel runs, and shows the progress (the improvements that `dbs --trace` reports);
- re-verifies the best circuits, independently of the search program, and gives the files of `slpd.py run`: the
  three circuits in the NIST format (with `XOR3` for the merged gates) and the summary (`.summary.json`);
- draws each circuit (View): the inputs, then one column per gate level, with the gate names of the text; hovering a
  gate highlights the gates and inputs that it uses, and the drawing can be saved as SVG.

For the matrices of `bench/`, it also gives the best XOR2 count and UMC area of `results/`.

**Differences from `slpd.py run`.**

- The random choices of a seed also depend on the C++ standard library (libc++ here, usually libstdc++ on Linux),
  so a run can give other circuits than the same run of `slpd.py`.
- The WebAssembly build is slower than the native one; for the `paper` preset, use `slpd.py`.
- There is no `--timeout`. The Stop button ends the search and keeps the best circuits of the runs that finished.

## Files

| File | Content |
|---|---|
| `index.html` | The page, built by `build.py`; edit `src/` and rebuild instead |
| `src/core.js` | Ports of the matrix checks (`slpd/matrix.py`), of the evaluation, merge and NIST writer of `slpd/circuit.py`, and of the presets and the selection of `slpd/pipeline.py`; the drawing of a circuit (SVG) |
| `src/worker.js` | Runs the module on one job: a minimal WASI layer with the arguments, an in-memory file system, a clock, stdout and stderr |
| `src/app.js`, `src/index.html`, `src/style.css` | The page |
| `patch_dbs.py` | Writes the copy of `src/dbs.cpp` that is compiled, with three changes for wasm32 (below) |
| `build.py` | Builds `index.html` |
| `THIRD_PARTY.txt` | The licenses of the libraries in the module (also on the page) |

## Building

```bash
pip install ziglang==0.14.1
python3 web/build.py
```

`build.py` writes the copy of `src/dbs.cpp` with `patch_dbs.py`, compiles it with Zig 0.14.1 (clang, libc++ and
wasi-libc for the target `wasm32-wasi`), and inlines the module, the examples and `results/summary.json` into
`index.html`. The page gives the SHA-256 of `src/dbs.cpp` and of the module. `src/dbs.cpp` itself is not changed.
The three changes:

1. `long` (32 bits on wasm32, 64 bits on x86-64 Linux) becomes `long long`, with its literals and `printf` formats;
2. the output of `--save` goes through a `FILE*` instead of an assignment to `stdout` (a constant in wasi-libc);
3. `std::clock()` becomes a monotonic clock, since the browser gives no CPU time; it times only the diagnostics
   (`--progress`, `--trace`).

`make test` checks that `index.html` is the build of `web/src/` and of the data of the repository, that the changed
copy of `src/dbs.cpp`, compiled with g++, gives the same output as `bin/dbs`, and, with Node.js, that `src/core.js`
gives the presets, checks, merges and circuit files of `slpd/` (`tests/test_web.py`).

## Licenses

The page is under the MIT license of the repository. The module also contains parts of wasi-libc, of libc++ and
libc++abi, and of the compiler-rt of Zig; `THIRD_PARTY.txt` gives their licenses.
