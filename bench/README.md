# Benchmark matrices

These are the 23 matrices of Table 5 of Shi, Feng and Xu (ToSC 2023(4)), plus the variant of the
32×32 [LS16] matrix that the circuits of the HILL repository compute (Sec. 4 of the paper).

**Format.** One row per output, entries `0`/`1` separated by spaces; row `j` is output `y_j` and
column `c` is input `x_c`. This is the input format of `src/dbs.cpp`; the scripts read it with
`slpd/matrix.py`.

**Bit order of `AES.txt`.** As in the HILL repository, the bits of each byte are numbered from the
least significant one: column `8i+k` is bit `k` of input byte `i`, and row `8i+k` is bit `k` of
output byte `i`. `paper/checks/verify_gf_aes.py` checks a circuit against MixColumns computed in
GF(2^8), and reports this bit order.

**`benchmark.tsv`** lists, for each file:

- the label;
- the size;
- the minimum depth;
- the XOR2 count of Shi–Feng–Xu (their Table 5);
- the UMC 180 nm area and depth of the circuit in the HILL repository;
- the XOR2 count of the circuit that the HILL repository stores as its starting point.

**Where the matrices come from.** We obtained each matrix by evaluating the published circuits on
the unit vectors:

- the result circuits of the HILL repository (<https://github.com/13est-Wishes/Linear_Layers_by_YongchaoLi>,
  MIT license) give every file except `FSE_LiuSim16_32_Circulant_SFX.txt`;
- `FSE_LiuSim16_32_Circulant_SFX.txt` comes from the 110-gate circuit in the repository of
  Shi–Feng–Xu.

For 22 of the 23 matrices, the HILL repository and the circuits of Shi–Feng–Xu give the same matrix,
entry by entry. The exception is the 32×32 [LS16] matrix (`FSE_LiuSim16_32_Circulant.txt` is the
HILL version; it differs from the matrix of Shi–Feng–Xu in two entries). `paper/compare/` holds the
scripts of these checks.

| File | Matrix | Original design |
|---|---|---|
| `AES.txt` | AES MixColumns | [DR02] |
| `SmallAES.txt` | SmallScale AES | [CMR05] |
| `JOLTIK.txt` | Joltik | [JNP15] |
| `FSE_SKOP15_16_Hadamard.txt`, `FSE_SKOP15_32_Hadamard.txt` | Hadamard MDS | [SKOP15] |
| `FSE_SKOP15_16_Involutory.txt`, `FSE_SKOP15_32_Involutory.txt` | involutory MDS | [SKOP15] |
| `FSE_LiuSim16_16_Circulant.txt`, `FSE_LiuSim16_32_Circulant_SFX.txt`, `FSE_LiuSim16_32_Circulant.txt` | circulant MDS | [LS16] |
| `FSE_LiWang16_16_Circulant.txt`, `FSE_LiWang16_32.txt` | circulant / MDS | [LW16] |
| `FSE_LiWang16_16_Involutory.txt`, `FSE_LiWang16_32_Involutory.txt` | involutory MDS | [LW16] |
| `ToSC_SarSye16_16_Toeplitz.txt`, `ToSC_SarSye16_32_Toeplitz.txt` | Toeplitz MDS | [SS16] |
| `ToSC_SarSye16_16_Involutory.txt`, `ToSC_SarSye16_32_Involutory.txt` | involutory MDS | [SS16] |
| `ePrint_JeaPeySim_16.txt`, `ePrint_JeaPeySim_16_Involutory.txt`, `ePrint_JeaPeySim_32_Subfield.txt` | MDS | [JPST17] |
| `C_BeiKraLea16_32_Circulant.txt` | circulant MDS | [BKL16] |
| `KLSW.txt` | MDS | [KLSW17] |
| `168_130.txt` | involutory MDS, HW 168 #130 (matrix M of Shi–Feng–Xu) | [LSL+19] |

The full references are in the paper.

**Equivalences** (Sec. 4 of the paper):

- Joltik and [SKOP15] involutory 16×16 are the same matrix, and [LW16] involutory 16×16 is a row
  and column permutation of it.
- [LS16] and [LW16] circulant 16×16 are row and column permutations of each other.
- [KLSW17] is a row and column permutation of [LSL+19] HW160 #93.
