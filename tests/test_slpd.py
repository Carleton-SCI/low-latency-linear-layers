"""Tests of slpd/ (library and command-line tool) and of the search program: python3 -m unittest discover -s tests -v
(the last three tests build and run bin/dbs)"""
import itertools
import os
import random
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from slpd.matrix import read_matrices, check_matrix, write_mat, write_rows, MatrixError, Matrix
from slpd.circuit import Circuit, circ_text, read_circ_text, read_dbs_text, dbs_text, read_sfx_text, read_hill_text, read_text

TOY = ["1 0 1 1 0 1", "0 1 1 1 0 1", "1 0 1 0 0 1", "1 1 1 0 1 0"]


def tmpfile(text):
    f = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
    f.write(text)
    f.close()
    return f.name


def random_circuit(rng, n, G, m):
    gates = []
    for k in range(G):
        a, b = rng.sample(range(n + k), 2)
        gates.append([a, b])
    out = [rng.randrange(n, n + G) for _ in range(m)]
    return Circuit(n, gates, out)


class MatrixTests(unittest.TestCase):
    def test_formats_agree(self):
        a = read_matrices(tmpfile("1\n4 6\n" + "\n".join(TOY) + "\n"))[0]
        b = read_matrices(tmpfile("# bare rows\n" + "\n".join(TOY) + "\n"))[0]
        self.assertEqual(a.rows, b.rows)
        self.assertEqual((a.m, a.n), (4, 6))
        self.assertEqual(a.min_depth(), 2)

    def test_several_matrices(self):
        Ms = read_matrices(tmpfile("2\n2 3\n1 1 0\n0 1 1\n2 2\n1 1\n1 0\n"))
        self.assertEqual([(M.m, M.n) for M in Ms], [(2, 3), (2, 2)])

    def test_errors(self):
        bad = {"1\n3 4\n1 1 0 0\n0 1 1 0\n": "announces 3 rows",
               "1 1 0\n0 2 1\n": "is not 0 or 1",
               "1 1 0\n0 1\n": "the first row has 3",
               "1\n2 3\n1 1 0 1\n0 1 1 0\n": "announces 3 columns",
               "": "holds no matrix",
               "2\n2 2\n1 1\n0 1\n": "holds 1"}
        for text, msg in bad.items():
            with self.assertRaises(MatrixError) as e:
                read_matrices(tmpfile(text))
            self.assertIn(msg, str(e.exception))

    def test_checks(self):
        M = read_matrices(tmpfile("1 1 0\n0 0 0\n"))[0]
        errors, _, _ = check_matrix(M)
        self.assertTrue(any("all zero" in e for e in errors))
        M = read_matrices(tmpfile("1 1 0 1\n1 1 0 1\n0 0 1 0\n"))[0]
        errors, warnings, _ = check_matrix(M)
        self.assertEqual(errors, [])
        self.assertTrue(any("duplicate" in w for w in warnings))
        self.assertTrue(any("weight 1" in w for w in warnings))

    def test_write_read(self):
        M = read_matrices(tmpfile("\n".join(TOY)))[0]
        p = tmpfile("")
        write_mat(M, p)
        self.assertEqual(read_matrices(p)[0].rows, M.rows)
        write_rows(M, p)
        self.assertEqual(read_matrices(p)[0].rows, M.rows)


class CircuitTests(unittest.TestCase):
    def test_merge_is_correct_and_maximum(self):
        rng = random.Random(1)
        for _ in range(300):
            n, G, m = rng.randint(3, 6), rng.randint(2, 9), rng.randint(1, 4)
            C = random_circuit(rng, n, G, m).prune()
            if not C.gates:
                continue
            val, _ = C.evaluate()
            M = Matrix([val[o] for o in C.out], n)
            if any(r == 0 for r in M.rows):
                continue
            Mg = C.merge_xor3()
            self.assertEqual(Mg.check(M), [])
            self.assertLessEqual(Mg.depth(), C.depth())
            self.assertEqual(Mg.xor2_equivalent(), len(C.gates))
            # brute force: the largest set of (child -> consumer) edges, each node in one edge at most
            fan = {}
            for k, ops in enumerate(C.gates):
                for o in ops:
                    fan.setdefault(o, []).append(k)
            outs = set(C.out)
            edges = [(k, fan[n + k][0]) for k in range(len(C.gates))
                     if n + k not in outs and len(fan.get(n + k, [])) == 1]
            best = 0
            for r in range(len(edges), 0, -1):
                for S in itertools.combinations(edges, r):
                    nodes = [x for e in S for x in e]
                    if len(nodes) == len(set(nodes)):
                        best = r
                        break
                if best:
                    break
            self.assertEqual(Mg.counts()[1], best)

    def test_circ_round_trip(self):
        rng = random.Random(2)
        M = read_matrices(tmpfile("\n".join(TOY)))[0]
        for _ in range(50):
            C = random_circuit(rng, 6, 8, 4).prune()
            val, _ = C.evaluate()
            M2 = Matrix([val[o] for o in C.out], 6)
            if any(r == 0 for r in M2.rows):
                continue
            for K in (C, C.merge_xor3()):
                back = read_circ_text(circ_text(K, M2, "t"))
                self.assertEqual(back.check(M2), [])
                self.assertEqual(back.counts(), K.counts())
            self.assertEqual(read_dbs_text(dbs_text(C), n=6).check(M2), [])
        self.assertIsNotNone(M)

    def test_buf_outputs(self):
        M = read_matrices(tmpfile("1 1 0 1\n1 1 0 1\n0 0 1 0\n"))[0]
        C = Circuit(4, [[0, 1], [4, 3]], [5, 5, 2])
        text = circ_text(C, M, "dup")
        self.assertIn("BUF y2 y1", text)
        self.assertIn("BUF y3 x3", text)
        self.assertEqual(read_circ_text(text).check(M), [])

    def test_other_readers(self):
        sfx = "21use3 xors, split size 1\n1 2 4 3 7 5\nx3 x0 x1 y0 depth 1\nx4 x3 x2 y1 depth 2\nx5 x0 x2 y2 depth 1\n\nxors3\n"
        C = read_sfx_text(sfx)[0]
        M = Matrix([3, 7, 5], 3)
        self.assertEqual(C.check(M), [])
        hill = "Optimization (Area: 7.67) XOR2: 1, XOR3: 1, Depth: 1\nt[3] = t[0] ^ t[1] y[0]\nt[4] = t[0] ^ t[1] ^ t[2] y[1]\n"
        H = read_hill_text(hill)
        self.assertEqual(H.check(Matrix([3, 7], 3)), [])


class SearchProgramTests(unittest.TestCase):
    def test_toy_quick(self):
        from slpd.pipeline import run_search
        M = read_matrices(os.path.join(ROOT, "examples", "toy.mat.txt"))[0]
        out = tempfile.mkdtemp()
        s = run_search(M, out, preset="quick", jobs=2, log=lambda *a: None)
        self.assertEqual(s["best"]["xor2"]["xor2"], 8)            # the optimum at depth 2
        self.assertAlmostEqual(s["best"]["umc180"]["umc180"], 21.34)
        for k in ("xor2", "umc180", "stm130"):
            C = read_circ_text(read_text(os.path.join(out, s["best"][k]["file"])))
            self.assertEqual(C.check(M, 2), [])

    def test_stats_and_trace_leave_the_search_unchanged(self):
        import subprocess
        subprocess.run(["make", "-s", "-C", ROOT, "bin/dbs"], check=True)
        dbs = os.path.join(ROOT, "bin", "dbs")
        mat = os.path.join(ROOT, "bench", "SmallAES.txt")
        for opts in (["--tie", "sfx", "--md", "--beam", "200", "--slack", "0"], ["--tie", "all", "--beam", "100", "--slack", "1"]):
            a = subprocess.run([dbs, mat, "3"] + opts + ["--seed", "3"], capture_output=True, text=True, check=True)
            b = subprocess.run([dbs, mat, "3"] + opts + ["--seed", "3", "--stats", "--trace"], capture_output=True, text=True, check=True)
            self.assertEqual(a.stdout, b.stdout)
            self.assertIn("STATEND", b.stderr)
            self.assertIn("TRACE cpu=", b.stderr)

    def test_progress_leaves_the_forward_runs_unchanged(self):
        import subprocess
        subprocess.run(["make", "-s", "-C", ROOT, "bin/dbs"], check=True)
        dbs = os.path.join(ROOT, "bin", "dbs")
        mat = os.path.join(ROOT, "bench", "SmallAES.txt")
        for tie in ("bp", "sfx", "all"):
            opts = ["--tie", tie, "--runs", "3", "--seed", "7"]
            a = subprocess.run([dbs, mat, "3"] + opts, capture_output=True, text=True, check=True)
            b = subprocess.run([dbs, mat, "3"] + opts + ["--progress"], capture_output=True, text=True, check=True)
            self.assertEqual(a.stdout, b.stdout)
            self.assertIn("PROGRESS cpu=", b.stderr)
            self.assertIn(" added gates=", b.stderr)


if __name__ == "__main__":
    unittest.main()
