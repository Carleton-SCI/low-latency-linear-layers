# Matrix and circuit formats

The tool reads the matrix formats used by earlier work, and writes circuits in the format of the
NIST Circuit Complexity Project. Throughout, row `j` of a matrix is output `y_j`, column `c` is
input `x_c`, and `y_j` is the XOR of the inputs whose column has a 1 in row `j`.

## Input matrices

### 1. NIST / Boyar–Peralta format (recommended)

This is the format of the matrix files of the NIST Circuit Complexity Project (for example
`Joltik.mat.txt` in its list of linear maps). It is also the input format of the code of Shi, Feng
and Xu (`testmatrix.txt`).

```
1
4 6
1 0 1 1 0 1
0 1 1 1 0 1
1 0 1 0 0 1
1 1 1 0 1 0
```

- Line 1: the number of matrices `k` in the file.
- Line 2: the number of rows `m` (outputs) and of columns `n` (inputs).
- Then `m` lines of `n` entries, each `0` or `1`, separated by spaces.
- The last two parts repeat `k` times.
- The example is a 4 × 6 toy matrix (`examples/toy.mat.txt`).

### 2. Bare rows

The same rows with no header lines. This is the input of `src/dbs.cpp` and the format of the
benchmark files in `bench/`.

### Common rules

- Entries are separated by spaces or tabs. Rows written without separators (`101101`) and
  comma-separated rows are rejected, with a message that says so.
- Lines that start with `#` are comments. This is an extension; NIST files have none.
- Windows line ends and a UTF-8 byte-order mark are accepted.
- Limits of the search program: at most 64 inputs and 128 outputs.

### Checks made by `python3 slpd.py check FILE`, before every run

| Check | Result |
|---|---|
| file missing, empty or not text | error |
| header inconsistent with the data (count, rows, columns) | error, with the line number |
| an entry other than `0` or `1`, or rows of different lengths | error, with line and column |
| a row of zeros (an output that is the constant 0) | error |
| more than 64 inputs or 128 outputs | error |
| every row of weight 1 (nothing to compute) | error |
| duplicate rows | warning: computed once and shared |
| rows of weight 1 | warning: the output is an input, no gate |
| all-zero columns | warning: an unused input |

The checker also reports the size, the number of ones, the row weights, the rank, whether the
matrix is invertible, the minimum depth max<sub>j</sub> ⌈log2(HW(row j))⌉, and the XOR2 count
without sharing.

## Output circuits: the NIST `.circ` format

`slpd.py run` writes three circuits per matrix:

- `NAME.xor2.circ.txt`: the fewest XOR2 gates;
- `NAME.umc180.circ.txt`: the smallest UMC 180 nm area after XOR3 merging;
- `NAME.stm130.circ.txt`: the smallest STM 130 nm area after XOR3 merging.

The format is that of the circuit files of the NIST Circuit Complexity Project:

```
# File timestamp (UTC): 2026-10-28T02:00:00
# Written by slpd (low-latency linear layers); format of the NIST Circuit Complexity Project
# https://csrc.nist.gov/projects/circuit-complexity

begin circuit toy--umc180

# Boolean Circuit for a linear system y = A.x defined by a 4x6 bit-matrix A
# Matrix A weights (W): totalW=15, adjW=11, minWInRow=3, maxWInRow=4, minWInCol=1, maxWICol=4
# Matrix represented as 4 rows: vecRows=UInt8[0x2d, 0x2e, 0x25, 0x17]
# Matrix represented as 6 columns: vecCols=UInt8[0x0d, 0x0a, 0x0f, 0x03, 0x08, 0x07]
# The element in the i-th row and j-th col of A is A[i,j] = 1 & (vecRows[i]>>(j-1)) = 1 & (vecCols[j]>>(i-1))

# Tally: 6 inputs, 4 outputs, 6 gates (4 XOR, 2 XOR3)
# Depth: 2 (gate levels; an XOR3 gate counts as one level)
# Area (unit-area sums, HILL ToSC 2026(1) Table 1): UMC 180 nm 21.34 GE, STM 130 nm 22.64 GE

Inputs: x1:x6
Outputs: y1:y4
Internal: t1:t2
GateSyntax: GateName Output Inputs

begin SLP
XOR t1 x3 x1
XOR y3 x6 t1
XOR3 y4 t1 x5 x2
XOR t2 x6 x4
XOR y1 t1 t2
XOR3 y2 t2 x3 x2
end SLP
end circuit
```

This example is the output of `python3 slpd.py run examples/toy.mat.txt`.

- **Numbering.** As in the NIST files, inputs and outputs are numbered from 1: `x1` is the first
  column and `y1` the first row. In this repository's other files (the paper, `bench/`, the native
  format of `src/dbs.cpp`), they are numbered from 0.
- **Gates.** One gate per line: `XOR out a b`. Operands are inputs `x`, internal signals `t` or
  outputs `y` defined on an earlier line.
- **Extensions** (not in the NIST files):
  - `XOR3 out a b c`: a 3-input XOR gate.
  - `BUF out s`: an output equal to an input or to another output. This happens only for rows of
    weight 1 and for duplicate rows; it costs nothing.
- **Header.** The comment lines repeat the matrix as integers, as the NIST files do (bit `j-1` of
  `vecRows[i]` is the entry in row `i`, column `j`). They also give the tally, the depth and the
  areas, with the unit areas of HILL (ToSC 2026(1), Table 1): XOR2 = 3.00 GE and XOR3 = 4.67 GE
  for UMC 180 nm, XOR2 = 3.33 GE and XOR3 = 4.66 GE for STM 130 nm.
- **Depth.** The depth counts gate levels; an XOR3 gate is one level. Merging never increases it.

## Other circuit formats read by the tool

`slpd.py verify`, `merge` and `convert` detect these formats from the file's content:

| Format | Gate lines | Numbering |
|---|---|---|
| native format of `src/dbs.cpp` | `gK = A + B  [d]` and `yJ = S` | `x0` = column 0 |
| result files of the Shi–Feng–Xu code | `xK xA xB [yJ ...] depth d`, after a header `... use N xors` | `x0` = column 0 |
| result files of the HILL repository | `t[k] = t[a] ^ t[b] [^ t[c]] [y[j] ...]` | `t[0]` = column 0 |
| NIST `.circ` | as above | `x1` = column 0 |

`python3 slpd.py convert MATRIX CIRCUIT -o OUT.circ.txt` turns any of them into the NIST format,
after checking that the circuit computes the matrix.

## Sources

- NIST Circuit Complexity Project, list of circuits and linear maps:
  <https://csrc.nist.gov/projects/circuit-complexity>.
- Shi, Feng and Xu, ToSC 2023(4): the README of their code describes `testmatrix.txt`.
- Li, Wei, Li and Guo, "HILL", ToSC 2026(1): the result files of its repository, and the unit areas
  (their Table 1).
