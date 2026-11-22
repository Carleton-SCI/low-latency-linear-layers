"""Binary matrices: reading, sanity checks and writing.

A matrix M (m x n over F2) defines the linear map y = M x: row j is output y_j and column c is
input x_c (0-based here; the NIST files number both from 1). A row is stored as an integer whose
bit c is the entry in column c.

Accepted input formats (see docs/FORMATS.md):

1. NIST / Boyar-Peralta matrix file (".mat.txt"), also the input format of the Shi-Feng-Xu code:
       <number of matrices k>
       <rows m> <columns n>
       m lines of n entries 0/1 separated by spaces
       (repeated k times)
2. Bare rows: m lines of n entries 0/1 separated by spaces, nothing else.

Lines starting with '#' are comments in both formats (an extension; NIST files have none).
"""
import math
import os

MAX_INPUTS = 64     # limits of the search program (src/dbs.cpp)
MAX_OUTPUTS = 128


class MatrixError(ValueError):
    """The file is not a valid matrix file; the message says where and why."""


class Matrix:
    def __init__(self, rows, n, name="matrix"):
        self.rows = list(rows)          # row j as an int, bit c = column c
        self.n = n                      # number of inputs (columns)
        self.m = len(self.rows)         # number of outputs (rows)
        self.name = name

    def entry(self, j, c):
        return (self.rows[j] >> c) & 1

    def as_lists(self):
        return [[self.entry(j, c) for c in range(self.n)] for j in range(self.m)]

    def row_weights(self):
        return [bin(r).count("1") for r in self.rows]

    def col_weights(self):
        return [sum(self.entry(j, c) for j in range(self.m)) for c in range(self.n)]

    def min_depth(self):
        """theoretically minimum depth: max over the rows of ceil(log2(HW(row)))"""
        return max((math.ceil(math.log2(w)) if w > 1 else 0) for w in self.row_weights())

    def rank(self):
        rows, r = list(self.rows), 0
        for c in range(self.n):
            piv = next((i for i in range(r, len(rows)) if rows[i] >> c & 1), None)
            if piv is None:
                continue
            rows[r], rows[piv] = rows[piv], rows[r]
            for i in range(len(rows)):
                if i != r and rows[i] >> c & 1:
                    rows[i] ^= rows[r]
            r += 1
        return r

    def transpose(self):
        cols = [sum(self.entry(j, c) << j for j in range(self.m)) for c in range(self.n)]
        return Matrix(cols, self.m, self.name + "_T")


def _tokens(line):
    return line.split()


def read_matrices(path):
    """Read every matrix of a file in either accepted format. Raises MatrixError with the line
    number of the first problem."""
    if not os.path.isfile(path):
        raise MatrixError(f"{path}: no such file")
    try:
        with open(path, encoding="utf-8-sig") as f:      # utf-8-sig: a byte-order mark is dropped
            raw = f.read()
    except UnicodeDecodeError:
        raise MatrixError(f"{path}: not a text file")
    lines = []                                    # (line number, tokens) of the non-comment lines
    for k, line in enumerate(raw.splitlines(), 1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        lines.append((k, _tokens(s)))
    if not lines:
        raise MatrixError(f"{path}: the file holds no matrix")
    base = os.path.basename(path)
    for suf in (".mat.txt", ".txt", ".mat"):
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    header = (len(lines) >= 2 and len(lines[0][1]) == 1 and len(lines[1][1]) == 2
              and all(t.isdigit() for t in lines[0][1] + lines[1][1]))
    if not header:
        return [_parse_rows(path, lines, base)]
    count = int(lines[0][1][0])
    if count < 1:
        raise MatrixError(f"{path}, line {lines[0][0]}: the number of matrices is {count}")
    out, pos = [], 1
    for k in range(count):
        if pos >= len(lines):
            raise MatrixError(f"{path}: the header announces {count} matrices, the file holds {k}")
        ln, t = lines[pos]
        if len(t) != 2 or not all(x.isdigit() for x in t):
            raise MatrixError(f"{path}, line {ln}: expected '<rows> <columns>', found {' '.join(t)!r}")
        m, n = int(t[0]), int(t[1])
        if m < 1 or n < 1:
            raise MatrixError(f"{path}, line {ln}: a matrix needs at least one row and one column")
        body = lines[pos + 1: pos + 1 + m]
        if len(body) < m:
            raise MatrixError(f"{path}, line {ln}: the header announces {m} rows, the file holds {len(body)}")
        name = base if count == 1 else f"{base}_{k + 1}"
        M = _parse_rows(path, body, name)
        if M.n != n:
            raise MatrixError(f"{path}, line {body[0][0]}: the header announces {n} columns, the rows have {M.n}")
        out.append(M)
        pos += 1 + m
    if pos < len(lines):
        raise MatrixError(f"{path}, line {lines[pos][0]}: unexpected data after the last matrix")
    return out


def _parse_rows(path, lines, name):
    n = len(lines[0][1])
    rows = []
    for ln, t in lines:
        if len(t) != n:
            raise MatrixError(f"{path}, line {ln}: {len(t)} entries, the first row has {n}")
        v = 0
        for c, x in enumerate(t):
            if x not in ("0", "1"):
                hint = ""
                if len(x) > 1 and set(x) <= {"0", "1"}:
                    hint = " (separate the entries by spaces)"
                elif "," in x or ";" in x:
                    hint = " (separate the entries by spaces or tabs, not commas)"
                raise MatrixError(f"{path}, line {ln}, column {c + 1}: entry {x!r} is not 0 or 1{hint}")
            if x == "1":
                v |= 1 << c
        rows.append(v)
    return Matrix(rows, n, name)


def _rows(word, xs):
    return f"{word}{'s' if len(xs) > 1 else ''} {', '.join(map(str, xs))}"


def _verb(xs, one, many):
    return many if len(xs) > 1 else one


def check_matrix(M):
    """Sanity checks for the search. Returns (errors, warnings, info); the search runs only
    without errors. (web/src/core.js has the same checks and messages.)"""
    errors, warnings = [], []
    if M.n > MAX_INPUTS:
        errors.append(f"{M.n} inputs: the search program supports at most {MAX_INPUTS}")
    if M.m > MAX_OUTPUTS:
        errors.append(f"{M.m} outputs: the search program supports at most {MAX_OUTPUTS}")
    zero = [j + 1 for j, r in enumerate(M.rows) if r == 0]
    if zero:
        errors.append(f"{_rows('row', zero)} " + _verb(zero, "is all zero (an output that is the constant 0)",
                                                       "are all zero (outputs that are the constant 0)"))
    w = M.row_weights()
    if not errors and max(w) < 2:
        errors.append("every row has weight 1: every output is an input, there is nothing to compute")
    seen, dup = {}, []
    for j, r in enumerate(M.rows):
        if r in seen and r:
            dup.append((j + 1, seen[r] + 1))
        seen.setdefault(r, j)
    if dup:
        warnings.append("duplicate rows (computed once, the other outputs reuse the signal): "
                        + ", ".join(f"row {a} = row {b}" for a, b in dup))
    ones = [j + 1 for j, x in enumerate(w) if x == 1]
    if ones:
        warnings.append(f"{_rows('row', ones)} {_verb(ones, 'has', 'have')} weight 1: "
                        f"{_verb(ones, 'this output is an input', 'these outputs are inputs')} and need{_verb(ones, 's', '')} no gate")
    unused = [c + 1 for c, x in enumerate(M.col_weights()) if x == 0]
    if unused:
        warnings.append(f"{_rows('column', unused)} {_verb(unused, 'is', 'are')} all zero: "
                        f"{_verb(unused, 'this input is', 'these inputs are')} not used")
    rk = M.rank()
    info = dict(name=M.name, rows=M.m, columns=M.n, ones=sum(w),
                row_weights=(min(w), max(w)), min_depth=M.min_depth() if not zero else None,
                rank=rk, invertible=(M.m == M.n and rk == M.n),
                naive_xor2=sum(x - 1 for x in w if x > 0))
    return errors, warnings, info


def write_rows(M, path):
    """bare-row format, the input of src/dbs.cpp"""
    with open(path, "w") as f:
        for row in M.as_lists():
            f.write(" ".join(map(str, row)) + "\n")


def write_mat(Ms, path):
    """NIST / Boyar-Peralta format, one or more matrices"""
    if isinstance(Ms, Matrix):
        Ms = [Ms]
    with open(path, "w") as f:
        f.write(f"{len(Ms)}\n")
        for M in Ms:
            f.write(f"{M.m} {M.n}\n")
            for row in M.as_lists():
                f.write(" ".join(map(str, row)) + "\n")
