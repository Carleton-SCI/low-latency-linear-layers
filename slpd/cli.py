"""Command line: python3 slpd.py <command> ...  (python3 slpd.py <command> -h for details)

  check   MATRIX                    sanity checks and properties of the matrices in a file
  run     MATRIX [options]          the four-stage search; writes the best circuits and a summary
  verify  MATRIX CIRCUIT [-D d]     re-evaluate a circuit (any supported format) against a matrix
  merge   MATRIX CIRCUIT [-o FILE]  merge an XOR2 circuit into XOR2/XOR3 gates (maximum merge)
  convert MATRIX CIRCUIT -o FILE    convert a circuit to the NIST .circ format (or --to dbs)
"""
import argparse
import os
import sys

from .matrix import read_matrices, check_matrix, MatrixError
from .circuit import read_any, circ_text, dbs_text, CircuitError
from . import pipeline


def _load(path):
    try:
        return read_matrices(path)
    except MatrixError as e:
        sys.exit(f"error: {e}")


def _print_check(M):
    errors, warnings, info = check_matrix(M)
    print(f"{info['name']}: {info['rows']} x {info['columns']} matrix, {info['ones']} ones, "
          f"row weights {info['row_weights'][0]}..{info['row_weights'][1]}, rank {info['rank']}"
          + (", invertible" if info["invertible"] else ""))
    if info["min_depth"] is not None:
        print(f"  minimum depth {info['min_depth']} (max over the rows of ceil(log2(weight))); "
              f"{info['naive_xor2']} XOR2 gates without sharing")
    for w in warnings:
        print(f"  warning: {w}")
    for e in errors:
        print(f"  ERROR: {e}")
    return not errors


def cmd_check(a):
    ok = all([_print_check(M) for M in _load(a.matrix)])
    if ok:
        print("OK: the file is a valid input for 'run'")
    sys.exit(0 if ok else 1)


def cmd_run(a):
    Ms = _load(a.matrix)
    if not all([_print_check(M) for M in Ms]):
        sys.exit("error: fix the matrix file first (see the ERROR lines above)")
    stages = [int(x) for x in a.stages.split(",")]
    if any(s not in (1, 2, 3, 4) for s in stages):
        sys.exit("error: --stages takes a comma-separated subset of 1,2,3,4")
    os.makedirs(a.out, exist_ok=True)
    for M in Ms:
        print(f"\n== {M.name}: preset {a.preset}, depth bound {a.depth or M.min_depth()} ==")
        try:
            s = pipeline.run_search(M, a.out, a.depth, a.preset, stages, a.jobs, a.timeout)
        except (RuntimeError, ValueError) as e:
            sys.exit(f"error: {e}")
        for k, b in s["best"].items():
            what = "fewest XOR2 gates" if k == "xor2" else f"smallest {k.upper()} area"
            print(f"  {what:26s}: {b['xor2']} XOR2 + {b['xor3']} XOR3, depth {b['depth']}, "
                  f"UMC {b['umc180']:.2f} GE, STM {b['stm130']:.2f} GE  -> {os.path.join(a.out, b['file'])}")
        print(f"  {s['seconds']:.0f} s in total; summary: {os.path.join(a.out, M.name + '.summary.json')}")


def _pick(Ms, a):
    if len(Ms) > 1 and a.index is None:
        sys.exit(f"error: {a.matrix} holds {len(Ms)} matrices; choose one with --index")
    return Ms[(a.index or 1) - 1]


def _read_circuit(path, M):
    try:
        return read_any(path, M.n, M.m)
    except (CircuitError, OSError) as e:
        sys.exit(f"error: {e}")


def cmd_verify(a):
    M = _pick(_load(a.matrix), a)
    C = _read_circuit(a.circuit, M)
    problems = C.check(M, a.depth)
    n2, n3 = C.counts()
    print(f"{a.circuit}: {n2} XOR2 + {n3} XOR3 gates, depth {C.depth() if not problems else '?'}, "
          f"UMC {C.area('umc180'):.2f} GE, STM {C.area('stm130'):.2f} GE")
    if n3 == 0 and not problems:
        Mg = C.prune().merge_xor3()
        g2, g3 = Mg.counts()
        print(f"  after a maximum XOR3 merge: {g2} XOR2 + {g3} XOR3, UMC {Mg.area('umc180'):.2f} GE, "
              f"STM {Mg.area('stm130'):.2f} GE")
    if problems:
        print("FAILED: " + "; ".join(problems[:10]) + (" ..." if len(problems) > 10 else ""))
        sys.exit(1)
    print("OK: the circuit computes the matrix" + (f" within depth {a.depth}" if a.depth is not None else ""))


def cmd_merge(a):
    M = _pick(_load(a.matrix), a)
    C = _read_circuit(a.circuit, M).prune()
    if C.check(M):
        sys.exit("error: the circuit does not compute the matrix: " + "; ".join(C.check(M)[:5]))
    try:
        Mg = C.merge_xor3()
    except CircuitError as e:
        sys.exit(f"error: {e}")
    out = a.out or os.path.splitext(a.circuit)[0] + ".merged.circ.txt"
    with open(out, "w") as f:
        f.write(circ_text(Mg, M, M.name + "--merged", f"maximum XOR3 merge of {os.path.basename(a.circuit)}",
                          [f"XOR2 circuit before merging: {len(C.gates)} gates, depth {C.depth()}"]))
    g2, g3 = Mg.counts()
    print(f"{len(C.gates)} XOR2 -> {g2} XOR2 + {g3} XOR3, UMC {Mg.area('umc180'):.2f} GE, "
          f"STM {Mg.area('stm130'):.2f} GE -> {out}")


def cmd_convert(a):
    M = _pick(_load(a.matrix), a)
    C = _read_circuit(a.circuit, M)
    if C.check(M):
        sys.exit("error: the circuit does not compute the matrix: " + "; ".join(C.check(M)[:5]))
    text = dbs_text(C) if a.to == "dbs" else circ_text(C, M, M.name, f"converted from {os.path.basename(a.circuit)}")
    with open(a.out, "w") as f:
        f.write(text)
    print(f"written {a.out}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="slpd.py", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("check", help="sanity checks of a matrix file")
    p.add_argument("matrix")
    p.set_defaults(fn=cmd_check)
    p = sub.add_parser("run", help="run the four-stage search on every matrix of a file")
    p.add_argument("matrix")
    p.add_argument("-D", "--depth", type=int, help="depth bound (default: the minimum depth)")
    p.add_argument("--preset", choices=("quick", "standard", "paper"), default="quick",
                   help="search budget; 'paper' is the configuration of the paper (default: quick)")
    p.add_argument("--stages", default="1,2,3,4", help="comma-separated subset of 1,2,3,4 (default: all)")
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count(), help="parallel runs (default: all CPUs)")
    p.add_argument("--timeout", type=float, help="time limit per run in seconds (a run over it is dropped)")
    p.add_argument("-o", "--out", default="slpd_out", help="output directory (default: slpd_out)")
    p.set_defaults(fn=cmd_run)
    for name, fn, hlp in (("verify", cmd_verify, "re-evaluate a circuit against a matrix"),
                          ("merge", cmd_merge, "merge an XOR2 circuit into XOR2/XOR3 gates"),
                          ("convert", cmd_convert, "convert a circuit to the NIST .circ format or to dbs")):
        p = sub.add_parser(name, help=hlp)
        p.add_argument("matrix")
        p.add_argument("circuit")
        p.add_argument("--index", type=int, help="which matrix of a file with several (1-based)")
        if name == "verify":
            p.add_argument("-D", "--depth", type=int, help="also check this depth bound")
        if name == "merge":
            p.add_argument("-o", "--out", help="output file (default: <circuit>.merged.circ.txt)")
        if name == "convert":
            p.add_argument("-o", "--out", required=True)
            p.add_argument("--to", choices=("circ", "dbs"), default="circ")
        p.set_defaults(fn=fn)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
