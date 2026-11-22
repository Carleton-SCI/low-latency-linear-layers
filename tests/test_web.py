"""Tests of the browser version (web/): python3 -m unittest discover -s tests -v
(the second test compiles web/patch_dbs.py's copy of src/dbs.cpp with g++; the third needs Node.js and is
skipped without it)"""
import base64
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from slpd.matrix import read_matrices, check_matrix, Matrix
from slpd.circuit import Circuit, circ_text
from slpd.pipeline import preset_jobs
from test_slpd import random_circuit, tmpfile

WEB = os.path.join(ROOT, "web")


def text(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def ensure_dbs():
    if not os.path.exists(os.path.join(ROOT, "bin", "dbs")):
        subprocess.run(["make", "-s", "-C", ROOT], check=True)
    return os.path.join(ROOT, "bin", "dbs")


class WebTests(unittest.TestCase):
    def test_page_is_current(self):
        """web/index.html is what web/build.py makes from web/src/, examples/, bench/ and results/"""
        page = text(os.path.join(WEB, "index.html"))
        wasm = base64.b64decode(re.search(r'const WASM_BASE64 = "([A-Za-z0-9+/=]+)"', page).group(1))
        d = tempfile.mkdtemp()
        with open(os.path.join(d, "dbs.wasm"), "wb") as f:
            f.write(wasm)
        subprocess.run([sys.executable, os.path.join(WEB, "build.py"), "--wasm", os.path.join(d, "dbs.wasm"),
                        "-o", os.path.join(d, "index.html")], check=True, capture_output=True)
        self.assertEqual(text(os.path.join(d, "index.html")), page)

    def test_patch_keeps_results(self):
        """the changes of web/patch_dbs.py leave the output of the search program unchanged (x86-64 Linux)"""
        dbs, d = ensure_dbs(), tempfile.mkdtemp()
        src = os.path.join(d, "dbs_web.cpp")
        subprocess.run([sys.executable, os.path.join(WEB, "patch_dbs.py"), os.path.join(ROOT, "src", "dbs.cpp"), src],
                       check=True, capture_output=True)
        web = os.path.join(d, "dbs_web")
        subprocess.run([os.environ.get("CXX", "g++"), "-O2", "-std=c++17", "-o", web, src], check=True)
        joltik = os.path.join(ROOT, "bench", "JOLTIK.txt")
        runs = [["--tie", "all", "--beam", "100", "--slack", "1", "--seed", "3"],
                ["--bpbs", "5", "--maxsplits", "16", "--evals", "4", "--tie", "sfx", "--md", "--seed", "2"],
                ["--bpbs", "1", "--bbeam", "2", "--bopts", "4", "--evals", "2", "--maxsplits", "6", "--tie", "sfx", "--md", "--seed", "1"]]
        for k, args in enumerate(runs):
            outs = []
            for prog in (dbs, web):
                pre = os.path.join(d, f"{os.path.basename(prog)}_{k}")
                r = subprocess.run([prog, joltik, "3"] + args + ["--save", pre], check=True, capture_output=True, text=True)
                outs.append([r.stdout] + [text(f"{pre}.{m}.txt") for m in ("xor2", "umc180", "stm130")])
            self.assertEqual(outs[0], outs[1], args)
        start = os.path.join(d, "dbs_0.xor2.txt")
        outs = [subprocess.run([prog, joltik, "3", "--lns", start, "--iters", "5", "--fix", "2:14", "--beam", "50", "--slack", "1",
                                "--tie", "all", "--seed", "1"], check=True, capture_output=True, text=True).stdout for prog in (dbs, web)]
        self.assertEqual(outs[0], outs[1])

    @unittest.skipUnless(shutil.which("node"), "needs Node.js")
    def test_core_js_matches_slpd(self):
        """web/src/core.js: the presets of slpd/pipeline.py, its matrix checks, and its merge and .circ writer"""
        rng = random.Random(7)
        circuits = []
        for n, G, m in ((6, 9, 4), (8, 20, 8), (16, 40, 16), (32, 90, 32)):
            C = random_circuit(rng, n, G, m).prune()
            circuits.append(dict(n=n, gates=C.gates, out=C.out))
        texts = [text(os.path.join(ROOT, "examples", "toy.mat.txt")), "1 1 0 1\n1 1 0 1\n0 0 1 0\n1 0 0 0\n",
                 "2\n2 3\n1 1 0\n0 1 1\n2 2\n1 1\n1 0\n", "1 1 0\n0 0 0\n", "1 1 0 0 0\n0 1 1 0 0\n"]
        presets = [[p, n, m, st] for p in ("quick", "standard", "paper") for n, m in ((6, 4), (16, 16), (32, 32), (64, 128))
                   for st in (1, 2, 3, 4)]
        script = r"""
const core = require(process.argv[1]);
const inp = JSON.parse(require("fs").readFileSync(0, "utf8"));
const skip = (t) => t.split("\n").filter((l) => !l.startsWith("# File timestamp") && !l.startsWith("# Written by"));
const out = { presets: inp.presets.map(([p, n, m, st]) => core.presetJobs(p, n, m, st)), checks: [], circuits: [] };
for (const t of inp.texts) for (const M of core.parseMatrices(t, "m")) {
  const c = core.checkMatrix(M);
  out.checks.push([M.name, c.errors, c.warnings, c.info.minDepth, c.info.naiveXor2, c.info.rank]);
}
for (const c of inp.circuits) {
  const C = new core.Circuit(c.n, c.gates, c.out), { val } = C.evaluate();
  const M = new core.Matrix(c.out.map((o) => val[o]), c.n, "r");
  const G = C.mergeXor3();
  out.circuits.push([G.gates, G.out, G.check(M, C.depth()), skip(core.circText(G, M, "r--umc180", "test", ["extra"]))]);
}
console.log(JSON.stringify(out));
"""
        r = subprocess.run(["node", "-e", script, os.path.join(WEB, "src", "core.js")], check=True, capture_output=True, text=True,
                           input=json.dumps(dict(presets=presets, texts=texts, circuits=circuits)))
        js = json.loads(r.stdout)
        self.assertEqual(js["presets"], [[[a, b] for a, b in preset_jobs(p, n, m, st)] for p, n, m, st in presets])
        checks = []
        for t in texts:
            for M in read_matrices(tmpfile(t)):
                e, w, i = check_matrix(M)
                checks.append([e, w, i["min_depth"], i["naive_xor2"], i["rank"]])
        self.assertEqual([c[1:] for c in js["checks"]], checks)
        skip = lambda t: [l for l in t.split("\n") if not l.startswith(("# File timestamp", "# Written by"))]
        for c, (gates, out, problems, circ) in zip(circuits, js["circuits"]):
            C = Circuit(c["n"], c["gates"], c["out"])
            val, _ = C.evaluate()
            M = Matrix([val[o] for o in C.out], c["n"], "r")
            G = C.merge_xor3()
            self.assertEqual((gates, out, problems), ([list(g) for g in G.gates], list(G.out), []))
            self.assertEqual(circ, skip(circ_text(G, M, "r--umc180", "test", ["extra"])))


if __name__ == "__main__":
    unittest.main()
