"""XOR circuits: evaluation and verification, merging into XOR3 gates, areas, and the circuit
formats read and written by the tool (see docs/FORMATS.md).

Internal model: nodes 0..n-1 are the inputs (node c = input x_c = matrix column c); node n+k is
gate k, whose operands are 2 (XOR2) or 3 (XOR3) earlier nodes; out[j] is the node that gives
output y_j (several outputs may share a node, an output may be an input).
"""
import datetime
import re

# unit areas in GE (Li, Wei, Li and Guo, "HILL", ToSC 2026(1), Table 1)
AREAS = {"umc180": (3.00, 4.67), "stm130": (3.33, 4.66)}


class CircuitError(ValueError):
    pass


def read_text(path):
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        return f.read()


class Circuit:
    def __init__(self, n, gates, out):
        self.n = n
        self.gates = [list(g) for g in gates]   # operand node ids, each earlier than the gate
        self.out = list(out)

    # ------------------------------------------------------------------ evaluation
    def evaluate(self):
        """values (bit c = input c) and depths (gate levels; XOR3 counts as one level)"""
        val = [1 << c for c in range(self.n)]
        dep = [0] * self.n
        for k, ops in enumerate(self.gates):
            me = self.n + k
            if not 2 <= len(ops) <= 3:
                raise CircuitError(f"gate {k} has {len(ops)} operands")
            v = 0
            for o in ops:
                if not 0 <= o < me:
                    raise CircuitError(f"gate {k} uses node {o}, which is not defined before it")
                v ^= val[o]
            val.append(v)
            dep.append(max(dep[o] for o in ops) + 1)
        return val, dep

    def counts(self):
        n2 = sum(len(g) == 2 for g in self.gates)
        return n2, len(self.gates) - n2

    def xor2_equivalent(self):
        """number of XOR2 gates if every XOR3 gate were two XOR2 gates"""
        n2, n3 = self.counts()
        return n2 + 2 * n3

    def area(self, lib):
        a2, a3 = AREAS[lib]
        n2, n3 = self.counts()
        return a2 * n2 + a3 * n3

    def depth(self):
        _, dep = self.evaluate()
        return max((dep[o] for o in self.out), default=0)

    def check(self, M, D=None):
        """problems found when the circuit is compared with matrix M (and depth bound D)"""
        problems = []
        if self.n != M.n:
            return [f"the circuit has {self.n} inputs, the matrix {M.n}"]
        if len(self.out) != M.m:
            return [f"the circuit has {len(self.out)} outputs, the matrix {M.m}"]
        try:
            val, dep = self.evaluate()
        except CircuitError as e:
            return [str(e)]
        for j, o in enumerate(self.out):
            if val[o] != M.rows[j]:
                problems.append(f"output y{j + 1} is wrong")
            elif D is not None and dep[o] > D:
                problems.append(f"output y{j + 1} has depth {dep[o]} > {D}")
        return problems

    def prune(self):
        """remove the gates that no output uses"""
        need = [False] * (self.n + len(self.gates))
        for o in self.out:
            need[o] = True
        for k in range(len(self.gates) - 1, -1, -1):
            if need[self.n + k]:
                for o in self.gates[k]:
                    need[o] = True
        new_id, gates = list(range(self.n)), []
        for k, ops in enumerate(self.gates):
            if need[self.n + k]:
                new_id.append(self.n + len(gates))
                gates.append([new_id[o] for o in ops])
            else:
                new_id.append(None)
        return Circuit(self.n, gates, [new_id[o] for o in self.out])

    # ------------------------------------------------------------------ XOR3 merging
    def merge_xor3(self):
        """Absorb XOR2 gates of fan-out one that are not outputs into their consumer, which becomes
        an XOR3 gate. A gate absorbs at most one gate and an absorbed gate absorbs none, so the
        merges are a matching on the forest of fan-out-one gates; this takes a maximum matching
        (dynamic programming over the trees), as src/dbs.cpp does. Merging never increases the
        number of gate levels on a path."""
        if any(len(g) != 2 for g in self.gates):
            raise CircuitError("merge_xor3 expects an XOR2 circuit")
        n, G = self.n, len(self.gates)
        fan, cons = [0] * (n + G), [-1] * (n + G)
        is_out = [False] * (n + G)
        for o in self.out:
            is_out[o] = True
        for k, ops in enumerate(self.gates):
            for o in ops:
                fan[o] += 1
                cons[o] = k
        children = [[] for _ in range(G)]
        has_parent = [False] * G
        for k in range(G):
            if not is_out[n + k] and fan[n + k] == 1:
                children[cons[n + k]].append(k)
                has_parent[k] = True
        f0, f1, pick = [0] * G, [-1] * G, [-1] * G
        for v in range(G):                        # children have smaller indices
            base = sum(max(f0[c], f1[c]) for c in children[v])
            f0[v], f1[v] = base, -1
            for c in children[v]:
                val = base - max(f0[c], f1[c]) + f0[c] + 1
                if val > f1[v]:
                    f1[v], pick[v] = val, c
        matched_to = [-1] * G
        stack = [(v, True) for v in range(G - 1, -1, -1) if not has_parent[v]]
        while stack:
            v, may = stack.pop()
            use = may and f1[v] > f0[v]
            for c in children[v]:
                if use and c == pick[v]:
                    matched_to[c] = v
                    stack.append((c, False))
                else:
                    stack.append((c, True))
        ops = [list(g) for g in self.gates]
        for c in range(G):
            if matched_to[c] >= 0:
                p = matched_to[c]
                ops[p] = [o for o in ops[p] if o != n + c] + self.gates[c]
        new_id, gates = list(range(n)), []
        for k in range(G):
            if matched_to[k] >= 0:
                new_id.append(None)
            else:
                new_id.append(n + len(gates))
                gates.append(ops[k])
        gates = [[new_id[o] for o in g] for g in gates]
        return Circuit(n, gates, [new_id[o] for o in self.out])


# ====================================================================== readers
def read_any(path, n=None, m=None):
    """read a circuit in any supported format (detected from the content)"""
    text = read_text(path)
    if "begin SLP" in text:
        return read_circ_text(text, path)
    if re.search(r"^\s*g\d+\s*=", text, re.M):
        return read_dbs_text(text, path, n)
    if re.search(r"use\s*\d+\s*xors", text):
        impls = read_sfx_text(text, path)
        if not impls:
            raise CircuitError(f"{path}: no complete implementation in the Shi-Feng-Xu file")
        return impls[0]
    if re.search(r"^\s*t\[\d+\]\s*=", text, re.M):
        return read_hill_text(text, path, n)
    raise CircuitError(f"{path}: unknown circuit format")


def read_circ_text(text, path="circuit"):
    """NIST Circuit Complexity Project format (gates 'XOR out a b'); also XOR3 and BUF gates"""
    n = m = None
    mi = re.search(r"^Inputs:\s*x1:x(\d+)", text, re.M)
    mo = re.search(r"^Outputs:\s*y1:y(\d+)", text, re.M)
    if not mi or not mo:
        raise CircuitError(f"{path}: missing 'Inputs: x1:xN' or 'Outputs: y1:yM'")
    n, m = int(mi.group(1)), int(mo.group(1))
    body = re.search(r"begin SLP\s*\n(.*?)\nend SLP", text, re.S)
    if not body:
        raise CircuitError(f"{path}: missing 'begin SLP' ... 'end SLP'")
    node = {f"x{c + 1}": c for c in range(n)}
    gates, out = [], [None] * m
    for ln, line in enumerate(body.group(1).splitlines(), 1):
        line = line.split("#")[0].strip()
        if not line:
            continue
        t = line.split()
        op, dst, src = t[0].upper(), t[1] if len(t) > 1 else None, t[2:]
        for s in src:
            if s not in node:
                raise CircuitError(f"{path}: SLP line {ln}: operand {s} is used before it is defined")
        if op == "BUF" and len(src) == 1:
            ident = node[src[0]]
        elif (op == "XOR" and len(src) == 2) or (op == "XOR3" and len(src) == 3):
            gates.append([node[s] for s in src])
            ident = n + len(gates) - 1
        else:
            raise CircuitError(f"{path}: SLP line {ln}: unsupported gate {line!r} (XOR with 2 inputs, "
                               "XOR3 with 3 inputs or BUF with 1 input)")
        if dst in node:
            raise CircuitError(f"{path}: SLP line {ln}: {dst} is assigned twice")
        node[dst] = ident
        my = re.fullmatch(r"y(\d+)", dst)
        if my:
            j = int(my.group(1)) - 1
            if not 0 <= j < m:
                raise CircuitError(f"{path}: SLP line {ln}: output {dst} outside y1:y{m}")
            out[j] = ident
    missing = [f"y{j + 1}" for j, o in enumerate(out) if o is None]
    if missing:
        raise CircuitError(f"{path}: outputs {', '.join(missing)} are not assigned")
    return Circuit(n, gates, out)


def read_dbs_text(text, path="circuit", n=None):
    """native format of src/dbs.cpp: 'gK = A + B [d]' and 'yJ = S', with x<c> = column c"""
    gate_re = re.compile(r"^\s*(g\d+)\s*=\s*(\w+)\s*\+\s*(\w+)(?:\s*\+\s*(\w+))?")
    out_re = re.compile(r"^\s*y(\d+)\s*=\s*(\w+)\s*$")
    lines = text.splitlines()
    xs = [int(x) for x in re.findall(r"\bx(\d+)\b", text)]
    if n is None:
        n = max(xs) + 1 if xs else 0
    node = {f"x{c}": c for c in range(n)}
    gates, outs = [], {}
    for ln, line in enumerate(lines, 1):
        g = gate_re.match(line)
        if g:
            ops = [s for s in g.groups()[1:] if s]
            for s in ops:
                if s not in node:
                    raise CircuitError(f"{path}, line {ln}: operand {s} is used before it is defined")
            gates.append([node[s] for s in ops])
            node[g.group(1)] = n + len(gates) - 1
            continue
        o = out_re.match(line)
        if o:
            if o.group(2) not in node:
                raise CircuitError(f"{path}, line {ln}: output uses undefined {o.group(2)}")
            outs[int(o.group(1))] = node[o.group(2)]
    if not outs:
        raise CircuitError(f"{path}: no output lines 'yJ = S'")
    m = max(outs) + 1
    missing = [j for j in range(m) if j not in outs]
    if missing:
        raise CircuitError(f"{path}: outputs {missing} are not assigned")
    return Circuit(n, gates, [outs[j] for j in range(m)])


SFX_GATE = re.compile(r"x(\d+) x(\d+) x(\d+)((?: y\d+)*)(?: depth (\d+))?")
SFX_HEAD = re.compile(r"use\s*(\d+) xors")


def read_sfx_text(text, path="circuit"):
    """result files of the Shi-Feng-Xu code: gates 'xK xA xB [yJ ...] [depth d]' after a header
    'use N xors' and a line of node values; a file may hold several implementations"""
    text = text.replace("\r", "")
    starts = [mm.start() for mm in SFX_HEAD.finditer(text)]
    impls = []
    for k, s in enumerate(starts):
        chunk = text[s: starts[k + 1] if k + 1 < len(starts) else len(text)]
        lines = chunk.split("\n")
        nums = re.findall(r"\d+", re.split(r"x\d", lines[1])[0]) if len(lines) > 1 else []
        n = 0
        while n < len(nums) and int(nums[n]) == 1 << n:
            n += 1
        if n == 0:
            continue
        node = {i: i for i in range(n)}
        gates, outs = [], {}
        for g in SFX_GATE.finditer(chunk):
            a, b, z = int(g.group(2)), int(g.group(3)), int(g.group(1))
            if a not in node or b not in node or z in node:
                raise CircuitError(f"{path}: bad gate {g.group(0)!r}")
            gates.append([node[a], node[b]])
            node[z] = n + len(gates) - 1
            for y in re.findall(r"y(\d+)", g.group(4)):
                outs.setdefault(int(y), node[z])
        if not outs:
            continue
        m = max(outs) + 1
        if any(j not in outs for j in range(m)):
            continue
        impls.append(Circuit(n, gates, [outs[j] for j in range(m)]))
    return impls


def read_hill_text(text, path="circuit", n=None):
    """result files of the HILL repository: 't[k] = t[a] ^ t[b] [^ t[c]] [y[j] ...]', inputs t[0..n-1]"""
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("Optimization"):
            continue
        g = re.match(r"t\[(\d+)\]\s*=\s*(.+?)(?:\s+((?:y\[\d+\]\s*)+))?$", line)
        if not g:
            raise CircuitError(f"{path}: cannot parse {line!r}")
        rows.append((int(g.group(1)), [int(x) for x in re.findall(r"t\[(\d+)\]", g.group(2))],
                     [int(x) for x in re.findall(r"y\[(\d+)\]", g.group(3) or "")]))
    assigned = {z for z, _, _ in rows}
    used = {o for _, ops, _ in rows for o in ops} - assigned
    if n is None:
        n = max(used) + 1
    node = {i: i for i in range(n)}
    gates, outs = [], {}
    for z, ops, ys in rows:
        for o in ops:
            if o not in node:
                raise CircuitError(f"{path}: t[{o}] is used before it is defined")
        gates.append([node[o] for o in ops])
        node[z] = n + len(gates) - 1
        for y in ys:
            outs[y] = node[z]
    m = max(outs) + 1
    return Circuit(n, gates, [outs.get(j, -1) for j in range(m)])


# ====================================================================== writers
def _uint(width_bits):
    for w in (8, 16, 32, 64, 128):
        if width_bits <= w:
            return w
    return width_bits


def circ_text(C, M, name, source=None, extra=None):
    """the circuit in the NIST Circuit Complexity Project format, with the header lines of that
    project (matrix weights, rows and columns as integers, tally, depth). XOR3 gates are written
    'XOR3 out a b c' and an output that equals an input or another output 'BUF out s'."""
    n, m = C.n, len(C.out)
    val, dep = C.evaluate()
    uses_y = {}                                   # node -> name of the first output it gives
    for j, o in enumerate(C.out):
        uses_y.setdefault(o, f"y{j + 1}")
    names, lines, t = {c: f"x{c + 1}" for c in range(n)}, [], 0
    for k, ops in enumerate(C.gates):
        me = n + k
        if me in uses_y:
            names[me] = uses_y[me]
        else:
            t += 1
            names[me] = f"t{t}"
        op = "XOR" if len(ops) == 2 else "XOR3"
        lines.append(f"{op} {names[me]} " + " ".join(names[o] for o in ops))
    bufs = []
    for j, o in enumerate(C.out):
        if names[o] != f"y{j + 1}":
            bufs.append(f"BUF y{j + 1} {names[o]}")
    n2, n3 = C.counts()
    rows = M.rows
    cols = [sum(((rows[j] >> c) & 1) << j for j in range(m)) for c in range(n)]
    rw, cw = [bin(r).count("1") for r in rows], [bin(c).count("1") for c in cols]
    ur, uc = _uint(n), _uint(m)
    hexr = ", ".join(f"0x{r:0{ur // 4}x}" for r in rows)
    hexc = ", ".join(f"0x{c:0{uc // 4}x}" for c in cols)
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    gate_tally = f"{n2} XOR" + (f", {n3} XOR3" if n3 else "")
    out = [f"# File timestamp (UTC): {now}",
           "# Written by slpd (low-latency linear layers); format of the NIST Circuit Complexity Project",
           "# https://csrc.nist.gov/projects/circuit-complexity", "",
           f"begin circuit {name}", "",
           f"# Boolean Circuit for a linear system y = A.x defined by a {m}x{n} bit-matrix A",
           f"# Matrix A weights (W): totalW={sum(rw)}, adjW={sum(rw) - sum(1 for x in rw if x)}, "
           f"minWInRow={min(rw)}, maxWInRow={max(rw)}, minWInCol={min(cw)}, maxWICol={max(cw)}",
           f"# Matrix represented as {m} rows: vecRows=UInt{ur}[{hexr}]",
           f"# Matrix represented as {n} columns: vecCols=UInt{uc}[{hexc}]",
           "# The element in the i-th row and j-th col of A is A[i,j] = 1 & (vecRows[i]>>(j-1)) "
           "= 1 & (vecCols[j]>>(i-1))", ""]
    if source:
        out.append(f"# Source: {source}")
    out += [f"# Tally: {n} inputs, {m} outputs, {len(C.gates)} gates ({gate_tally})",
            f"# Depth: {max(dep[o] for o in C.out)}" + (" (gate levels; an XOR3 gate counts as one level)" if n3 else ""),
            f"# Area (unit-area sums, HILL ToSC 2026(1) Table 1): UMC 180 nm {C.area('umc180'):.2f} GE, "
            f"STM 130 nm {C.area('stm130'):.2f} GE"]
    for e in extra or []:
        out.append(f"# {e}")
    out += ["", f"Inputs: x1:x{n}", f"Outputs: y1:y{m}", f"Internal: t1:t{t}" if t else "Internal: none",
            "GateSyntax: GateName Output Inputs", "", "begin SLP"] + lines + bufs + ["end SLP", "end circuit", ""]
    return "\n".join(out)


def dbs_text(C):
    """native format of src/dbs.cpp (XOR2 circuits; read back by its --lns option)"""
    if any(len(g) != 2 for g in C.gates):
        raise CircuitError("the dbs format holds XOR2 circuits only")
    _, dep = C.evaluate()
    nm = lambda i: f"x{i}" if i < C.n else f"g{i - C.n}"
    lines = [f"{nm(C.n + k)} = {nm(a)} + {nm(b)}  [{dep[C.n + k]}]" for k, (a, b) in enumerate(C.gates)]
    lines += [f"y{j} = {nm(o)}" for j, o in enumerate(C.out)]
    return "\n".join(lines) + "\n"
