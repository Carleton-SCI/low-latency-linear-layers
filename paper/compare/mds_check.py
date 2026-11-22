"""Invertibility and MDS property of the benchmark matrices (Sec. 4 of the paper).

Each n x n binary matrix (n = 16 or 32) is read as a 4 x 4 matrix of (n/4) x (n/4) blocks, i.e., blocks of
4 x 4 bits for the 16 x 16 matrices and of 8 x 8 bits for the 32 x 32 matrices. The matrix is MDS when all
its 69 square block submatrices (16 + 36 + 16 + 1) are nonsingular over F_2. Prints the rank of each matrix
and the number of singular block submatrices.

usage: python3 paper/compare/mds_check.py BENCH_DIR      (BENCH_DIR contains benchmark.tsv and the matrices)"""
import csv
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from slpd.matrix import read_matrices  # noqa: E402


def rank(rows, ncols):
    rows, r = list(rows), 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(rows)) if (rows[i] >> c) & 1), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(len(rows)):
            if i != r and (rows[i] >> c) & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    bench = sys.argv[1]
    for r in csv.DictReader(open(os.path.join(bench, "benchmark.tsv")), delimiter="\t"):
        M = read_matrices(os.path.join(bench, r["file"]))[0]
        n, b = M.n, M.n // 4
        singular = total = 0
        for k in range(1, 5):
            for R in itertools.combinations(range(4), k):
                for C in itertools.combinations(range(4), k):
                    cols = [c for cb in C for c in range(cb * b, cb * b + b)]
                    sub = []
                    for rb in R:
                        for i in range(rb * b, rb * b + b):
                            sub.append(sum(1 << j for j, c in enumerate(cols) if (M.rows[i] >> c) & 1))
                    total += 1
                    singular += rank(sub, len(cols)) < len(cols)
        print(f"{r['file'][:-4]:32s} {n}x{n}, blocks of {b}x{b} bits: rank {rank(M.rows, n)}, "
              f"singular block submatrices {singular} of {total}")


if __name__ == "__main__":
    main()
