"""Check the circuits of results/ against the matrices of bench/: each circuit computes its matrix, every output is
within the minimum depth D of the matrix (also with each XOR3 gate a + (b + c) read as two XOR2 gates), its gate
counts and areas are those of results/summary.json, and its XOR2 + 2 XOR3 equals the XOR2 count before merging
plus the duplicated gates that its header states. Standard library only; less than a second.

usage: python3 paper/checks/verify_results.py [NAME ...]      (NAME: a file of bench/ without .txt; default: all)"""
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "paper", "tables"))
from slpd.circuit import read_any, read_text   # noqa: E402
from slpd.matrix import read_matrices      # noqa: E402
from appendix_circuits import xor2_depth   # noqa: E402

summary = json.load(open(os.path.join(ROOT, "results", "summary.json")))
names = sys.argv[1:] or list(summary)
failed = 0
for name in names:
    if name not in summary:
        sys.exit(f"{name}: not in results/summary.json")
    M = read_matrices(os.path.join(ROOT, "bench", name + ".txt"))[0]
    D = summary[name]["min_depth"]
    for measure, b in summary[name]["best"].items():
        path = os.path.join(ROOT, "results", b["file"])
        C = read_any(path, M.n, M.m)
        text = read_text(path)
        n2, n3 = C.counts()
        problems = C.check(M, D)
        if not problems and xor2_depth(C) > D:
            problems.append(f"depth {xor2_depth(C)} > {D} with the XOR3 gates read as two XOR2 gates")
        got = (n2, n3, C.area("umc180"), C.area("stm130"))
        want = (b["xor2"], b["xor3"], b["umc180"], b["stm130"])
        if got[:2] != want[:2] or any(abs(x - y) > 1e-6 for x, y in zip(got[2:], want[2:])):
            problems.append(f"[{n2} XOR2 + {n3} XOR3], {got[2]:.2f} / {got[3]:.2f} GE; summary.json: "
                            f"[{want[0]} + {want[1]}], {want[2]:.2f} / {want[3]:.2f} GE")
        before = b["xor2_before_merge"]
        dup = re.search(r"(\d+) gate\(s\) of fan-out two duplicated", text)
        stated = re.search(r"XOR2 circuit before merging: (\d+) gates", text)
        if n2 + 2 * n3 != before + (int(dup.group(1)) if dup else 0) or (stated and int(stated.group(1)) != before):
            problems.append(f"XOR2 + 2 XOR3 = {n2 + 2 * n3}, not {before} gates before merging plus the duplicated gates")
        status = "OK" if not problems else "FAILED: " + "; ".join(problems[:3])
        print(f"{b['file']:<48} [{n2:>3} XOR2 + {n3:>2} XOR3] from {before:>3} XOR2, D = {D}, "
              f"UMC {got[2]:7.2f} GE, STM {got[3]:7.2f} GE  {status}")
        failed += bool(problems)
total = sum(len(summary[n]["best"]) for n in names)
print(f"{total - failed} of {total} circuits verified")
sys.exit(1 if failed else 0)
