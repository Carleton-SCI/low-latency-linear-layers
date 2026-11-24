// Matrices, circuits and the four-stage search of slpd/pipeline.py, for the browser.
// Ports of slpd/matrix.py (reading and checks), slpd/circuit.py (evaluation, pruning, maximum XOR3 merge,
// the NIST .circ writer, the reader of the native format of src/dbs.cpp) and slpd/pipeline.py (presets,
// stages, selection of the best circuits). Rows and values are BigInts: bit c is input x_c.
"use strict";

const MAX_INPUTS = 64, MAX_OUTPUTS = 128;            // limits of src/dbs.cpp
const AREAS = { umc180: [3.00, 4.67], stm130: [3.33, 4.66] };   // GE (HILL, ToSC 2026(1), Table 1)

class MatrixError extends Error {}
class CircuitError extends Error {}

const popcount = (x) => { let c = 0; while (x) { x &= x - 1n; c++; } return c; };
const listOf = (word, xs) => `${word}${xs.length > 1 ? "s" : ""} ${xs.join(", ")}`;
const verb = (xs, one, many) => (xs.length > 1 ? many : one);
const bit = (x, c) => Number((x >> BigInt(c)) & 1n);

class Matrix {
  constructor(rows, n, name = "matrix") { this.rows = rows; this.n = n; this.m = rows.length; this.name = name; }
  rowWeights() { return this.rows.map(popcount); }
  colWeights() { return Array.from({ length: this.n }, (_, c) => this.rows.reduce((s, r) => s + bit(r, c), 0)); }
  minDepth() { return Math.max(...this.rowWeights().map((w) => (w > 1 ? Math.ceil(Math.log2(w)) : 0))); }
  rank() {
    const rows = [...this.rows];
    let r = 0;
    for (let c = 0; c < this.n; c++) {
      const piv = rows.findIndex((x, i) => i >= r && bit(x, c));
      if (piv < 0) continue;
      [rows[r], rows[piv]] = [rows[piv], rows[r]];
      for (let i = 0; i < rows.length; i++) if (i !== r && bit(rows[i], c)) rows[i] ^= rows[r];
      r++;
    }
    return r;
  }
  rowsText() {
    return this.rows.map((r) => Array.from({ length: this.n }, (_, c) => bit(r, c)).join(" ")).join("\n") + "\n";
  }
}

function parseRows(lines, name) {
  const n = lines[0][1].length, rows = [];
  for (const [ln, t] of lines) {
    if (t.length !== n) throw new MatrixError(`line ${ln}: ${t.length} entries, the first row has ${n}`);
    let v = 0n;
    t.forEach((x, c) => {
      if (x !== "0" && x !== "1") {
        let hint = "";
        if (x.length > 1 && /^[01]+$/.test(x)) hint = " (separate the entries by spaces)";
        else if (/[,;]/.test(x)) hint = " (separate the entries by spaces or tabs, not commas)";
        throw new MatrixError(`line ${ln}, column ${c + 1}: entry '${x}' is not 0 or 1${hint}`);
      }
      if (x === "1") v |= 1n << BigInt(c);
    });
    rows.push(v);
  }
  return new Matrix(rows, n, name);
}

// Every matrix of a text in the NIST / Boyar-Peralta format (k, then "m n" and m rows per matrix) or in bare rows;
// lines starting with '#' are comments. Throws MatrixError with the line of the first problem.
function parseMatrices(text, base = "matrix") {
  const lines = [];
  text.replace(/^﻿/, "").split(/\r\n|\r|\n/).forEach((line, i) => {
    const s = line.trim();
    if (s && !s.startsWith("#")) lines.push([i + 1, s.split(/\s+/)]);
  });
  if (!lines.length) throw new MatrixError("the input holds no matrix");
  const digits = (t) => /^\d+$/.test(t);
  const header = lines.length >= 2 && lines[0][1].length === 1 && lines[1][1].length === 2 &&
    [...lines[0][1], ...lines[1][1]].every(digits);
  if (!header) return [parseRows(lines, base)];
  const count = Number(lines[0][1][0]);
  if (count < 1) throw new MatrixError(`line ${lines[0][0]}: the number of matrices is ${count}`);
  const out = [];
  let pos = 1;
  for (let k = 0; k < count; k++) {
    if (pos >= lines.length) throw new MatrixError(`the header announces ${count} matrices, the input holds ${k}`);
    const [ln, t] = lines[pos];
    if (t.length !== 2 || !t.every(digits)) throw new MatrixError(`line ${ln}: expected '<rows> <columns>', found '${t.join(" ")}'`);
    const m = Number(t[0]), n = Number(t[1]);
    if (m < 1 || n < 1) throw new MatrixError(`line ${ln}: a matrix needs at least one row and one column`);
    const body = lines.slice(pos + 1, pos + 1 + m);
    if (body.length < m) throw new MatrixError(`line ${ln}: the header announces ${m} rows, the input holds ${body.length}`);
    const M = parseRows(body, count === 1 ? base : `${base}_${k + 1}`);
    if (M.n !== n) throw new MatrixError(`line ${body[0][0]}: the header announces ${n} columns, the rows have ${M.n}`);
    out.push(M);
    pos += 1 + m;
  }
  if (pos < lines.length) throw new MatrixError(`line ${lines[pos][0]}: unexpected data after the last matrix`);
  return out;
}

// Sanity checks for the search (slpd/matrix.py, check_matrix): the search runs only without errors.
function checkMatrix(M) {
  const errors = [], warnings = [];
  if (M.n > MAX_INPUTS) errors.push(`${M.n} inputs: the search program supports at most ${MAX_INPUTS}`);
  if (M.m > MAX_OUTPUTS) errors.push(`${M.m} outputs: the search program supports at most ${MAX_OUTPUTS}`);
  const zero = M.rows.flatMap((r, j) => (r === 0n ? [j + 1] : []));
  if (zero.length) errors.push(`${listOf("row", zero)} ${verb(zero, "is all zero (an output that is the constant 0)", "are all zero (outputs that are the constant 0)")}`);
  const w = M.rowWeights();
  if (!errors.length && Math.max(...w) < 2) errors.push("every row has weight 1: every output is an input, there is nothing to compute");
  const seen = new Map(), dup = [];
  M.rows.forEach((r, j) => { if (r && seen.has(r)) dup.push(`row ${j + 1} = row ${seen.get(r) + 1}`); if (!seen.has(r)) seen.set(r, j); });
  if (dup.length) warnings.push(`duplicate rows (computed once, the other outputs reuse the signal): ${dup.join(", ")}`);
  const ones = w.flatMap((x, j) => (x === 1 ? [j + 1] : []));
  if (ones.length) warnings.push(`${listOf("row", ones)} ${verb(ones, "has", "have")} weight 1: ${verb(ones, "this output is an input", "these outputs are inputs")} and need${verb(ones, "s", "")} no gate`);
  const unused = M.colWeights().flatMap((x, c) => (x === 0 ? [c + 1] : []));
  if (unused.length) warnings.push(`${listOf("column", unused)} ${verb(unused, "is", "are")} all zero: ${verb(unused, "this input is", "these inputs are")} not used`);
  const rank = M.rank();
  const info = {
    rows: M.m, columns: M.n, ones: w.reduce((a, b) => a + b, 0), minWeight: Math.min(...w), maxWeight: Math.max(...w),
    minDepth: zero.length ? null : M.minDepth(), rank, invertible: M.m === M.n && rank === M.n,
    naiveXor2: w.reduce((s, x) => s + (x > 0 ? x - 1 : 0), 0),
  };
  return { errors, warnings, info };
}

class Circuit {
  constructor(n, gates, out) { this.n = n; this.gates = gates.map((g) => [...g]); this.out = [...out]; }
  evaluate() {
    const val = Array.from({ length: this.n }, (_, c) => 1n << BigInt(c)), dep = new Array(this.n).fill(0);
    this.gates.forEach((ops, k) => {
      const me = this.n + k;
      if (ops.length < 2 || ops.length > 3) throw new CircuitError(`gate ${k} has ${ops.length} operands`);
      let v = 0n, d = 0;
      for (const o of ops) {
        if (!(o >= 0 && o < me)) throw new CircuitError(`gate ${k} uses node ${o}, which is not defined before it`);
        v ^= val[o]; d = Math.max(d, dep[o]);
      }
      val.push(v); dep.push(d + 1);
    });
    return { val, dep };
  }
  counts() { const n2 = this.gates.filter((g) => g.length === 2).length; return [n2, this.gates.length - n2]; }
  area(lib) { const [a2, a3] = AREAS[lib], [n2, n3] = this.counts(); return a2 * n2 + a3 * n3; }
  depth() { const { dep } = this.evaluate(); return Math.max(0, ...this.out.map((o) => dep[o])); }
  check(M, D = null) {
    if (this.n !== M.n) return [`the circuit has ${this.n} inputs, the matrix ${M.n}`];
    if (this.out.length !== M.m) return [`the circuit has ${this.out.length} outputs, the matrix ${M.m}`];
    let e;
    try { e = this.evaluate(); } catch (err) { return [err.message]; }
    const problems = [];
    this.out.forEach((o, j) => {
      if (e.val[o] !== M.rows[j]) problems.push(`output y${j + 1} is wrong`);
      else if (D !== null && e.dep[o] > D) problems.push(`output y${j + 1} has depth ${e.dep[o]} > ${D}`);
    });
    return problems;
  }
  prune() {
    const need = new Array(this.n + this.gates.length).fill(false);
    for (const o of this.out) need[o] = true;
    for (let k = this.gates.length - 1; k >= 0; k--) if (need[this.n + k]) for (const o of this.gates[k]) need[o] = true;
    const id = Array.from({ length: this.n }, (_, c) => c), gates = [];
    this.gates.forEach((ops, k) => {
      if (need[this.n + k]) { id.push(this.n + gates.length); gates.push(ops.map((o) => id[o])); } else id.push(null);
    });
    return new Circuit(this.n, gates, this.out.map((o) => id[o]));
  }
  // Absorb the XOR2 gates of fan-out one that are not outputs into the gate that uses them (a maximum matching
  // on the forest of these gates, by dynamic programming), as slpd/circuit.py and src/dbs.cpp do.
  mergeXor3() {
    if (this.gates.some((g) => g.length !== 2)) throw new CircuitError("the merge expects an XOR2 circuit");
    const n = this.n, G = this.gates.length;
    const fan = new Array(n + G).fill(0), cons = new Array(n + G).fill(-1), isOut = new Array(n + G).fill(false);
    for (const o of this.out) isOut[o] = true;
    this.gates.forEach((ops, k) => { for (const o of ops) { fan[o]++; cons[o] = k; } });
    const children = Array.from({ length: G }, () => []), hasParent = new Array(G).fill(false);
    for (let k = 0; k < G; k++) if (!isOut[n + k] && fan[n + k] === 1) { children[cons[n + k]].push(k); hasParent[k] = true; }
    const f0 = new Array(G).fill(0), f1 = new Array(G).fill(-1), pick = new Array(G).fill(-1);
    for (let v = 0; v < G; v++) {
      const base = children[v].reduce((s, c) => s + Math.max(f0[c], f1[c]), 0);
      f0[v] = base; f1[v] = -1;
      for (const c of children[v]) {
        const val = base - Math.max(f0[c], f1[c]) + f0[c] + 1;
        if (val > f1[v]) { f1[v] = val; pick[v] = c; }
      }
    }
    const matchedTo = new Array(G).fill(-1), stack = [];
    for (let v = G - 1; v >= 0; v--) if (!hasParent[v]) stack.push([v, true]);
    while (stack.length) {
      const [v, may] = stack.pop();
      const use = may && f1[v] > f0[v];
      for (const c of children[v]) {
        if (use && c === pick[v]) { matchedTo[c] = v; stack.push([c, false]); } else stack.push([c, true]);
      }
    }
    const ops = this.gates.map((g) => [...g]);
    for (let c = 0; c < G; c++) if (matchedTo[c] >= 0) { const p = matchedTo[c]; ops[p] = ops[p].filter((o) => o !== n + c).concat(this.gates[c]); }
    const id = Array.from({ length: n }, (_, c) => c), gates = [];
    for (let k = 0; k < G; k++) { if (matchedTo[k] >= 0) id.push(null); else { id.push(n + gates.length); gates.push(ops[k]); } }
    return new Circuit(n, gates.map((g) => g.map((o) => id[o])), this.out.map((o) => id[o]));
  }
}

// The native format of src/dbs.cpp: 'gK = A + B  [d]' and 'yJ = S', with x<c> = column c.
function readDbsText(text, n) {
  const node = new Map(Array.from({ length: n }, (_, c) => [`x${c}`, c]));
  const gates = [], outs = new Map();
  text.split(/\r?\n/).forEach((line, i) => {
    const g = line.match(/^\s*(g\d+)\s*=\s*(\w+)\s*\+\s*(\w+)(?:\s*\+\s*(\w+))?/);
    if (g) {
      const ops = g.slice(2).filter(Boolean);
      for (const s of ops) if (!node.has(s)) throw new CircuitError(`line ${i + 1}: operand ${s} is used before it is defined`);
      gates.push(ops.map((s) => node.get(s)));
      node.set(g[1], n + gates.length - 1);
      return;
    }
    const o = line.match(/^\s*y(\d+)\s*=\s*(\w+)\s*$/);
    if (o) {
      if (!node.has(o[2])) throw new CircuitError(`line ${i + 1}: output uses undefined ${o[2]}`);
      outs.set(Number(o[1]), node.get(o[2]));
    }
  });
  if (!outs.size) throw new CircuitError("no output lines 'yJ = S'");
  const m = Math.max(...outs.keys()) + 1;
  for (let j = 0; j < m; j++) if (!outs.has(j)) throw new CircuitError(`output ${j} is not assigned`);
  return new Circuit(n, gates, Array.from({ length: m }, (_, j) => outs.get(j)));
}

const uintWidth = (bits) => [8, 16, 32, 64, 128].find((w) => bits <= w) || bits;

// The circuit in the format of the NIST Circuit Complexity Project, as slpd/circuit.py writes it.
function circText(C, M, name, source, extra = []) {
  const n = C.n, m = C.out.length, { dep } = C.evaluate();
  const usesY = new Map();
  C.out.forEach((o, j) => { if (!usesY.has(o)) usesY.set(o, `y${j + 1}`); });
  const names = new Map(Array.from({ length: n }, (_, c) => [c, `x${c + 1}`])), lines = [];
  let t = 0;
  C.gates.forEach((ops, k) => {
    const me = n + k;
    names.set(me, usesY.has(me) ? usesY.get(me) : `t${++t}`);
    lines.push(`${ops.length === 2 ? "XOR" : "XOR3"} ${names.get(me)} ` + ops.map((o) => names.get(o)).join(" "));
  });
  const bufs = C.out.flatMap((o, j) => (names.get(o) !== `y${j + 1}` ? [`BUF y${j + 1} ${names.get(o)}`] : []));
  const [n2, n3] = C.counts();
  const rows = M.rows;
  const cols = Array.from({ length: n }, (_, c) => rows.reduce((s, r, j) => s | (BigInt(bit(r, c)) << BigInt(j)), 0n));
  const rw = rows.map(popcount), cw = cols.map(popcount);
  const ur = uintWidth(n), uc = uintWidth(m);
  const hex = (x, w) => "0x" + x.toString(16).padStart(w / 4, "0");
  const sum = rw.reduce((a, b) => a + b, 0);
  const now = new Date().toISOString().slice(0, 19);
  const out = [
    `# File timestamp (UTC): ${now}`,
    "# Written by slpd (low-latency linear layers), in the browser; format of the NIST Circuit Complexity Project",
    "# https://csrc.nist.gov/projects/circuit-complexity", "",
    `begin circuit ${name}`, "",
    `# Boolean Circuit for a linear system y = A.x defined by a ${m}x${n} bit-matrix A`,
    `# Matrix A weights (W): totalW=${sum}, adjW=${sum - rw.filter((x) => x).length}, minWInRow=${Math.min(...rw)}, ` +
      `maxWInRow=${Math.max(...rw)}, minWInCol=${Math.min(...cw)}, maxWICol=${Math.max(...cw)}`,
    `# Matrix represented as ${m} rows: vecRows=UInt${ur}[${rows.map((r) => hex(r, ur)).join(", ")}]`,
    `# Matrix represented as ${n} columns: vecCols=UInt${uc}[${cols.map((c) => hex(c, uc)).join(", ")}]`,
    "# The element in the i-th row and j-th col of A is A[i,j] = 1 & (vecRows[i]>>(j-1)) = 1 & (vecCols[j]>>(i-1))", "",
  ];
  if (source) out.push(`# Source: ${source}`);
  out.push(`# Tally: ${n} inputs, ${m} outputs, ${C.gates.length} gates (${n2} XOR${n3 ? `, ${n3} XOR3` : ""})`,
    `# Depth: ${Math.max(...C.out.map((o) => dep[o]))}` + (n3 ? " (gate levels; an XOR3 gate counts as one level)" : ""),
    `# Area (unit-area sums, HILL ToSC 2026(1) Table 1): UMC 180 nm ${C.area("umc180").toFixed(2)} GE, ` +
      `STM 130 nm ${C.area("stm130").toFixed(2)} GE`);
  for (const e of extra) out.push(`# ${e}`);
  out.push("", `Inputs: x1:x${n}`, `Outputs: y1:y${m}`, t ? `Internal: t1:t${t}` : "Internal: none",
    "GateSyntax: GateName Output Inputs", "", "begin SLP", ...lines, ...bufs, "end SLP", "end circuit", "");
  return out.join("\n");
}

// ------------------------------------------------------------------ the search (slpd/pipeline.py)
const range = (a, b) => Array.from({ length: b - a }, (_, i) => a + i);

// Argument lists of one stage. 'paper' is the configuration of the paper (its 16x16 settings for n <= 16 inputs,
// its 32x32 settings otherwise); 'quick' and 'standard' are smaller budgets with the same structure.
function presetJobs(preset, n, m, stage) {
  const small = n <= 16;
  const s1 = (grid) => grid.flatMap(([t, w, s, seeds]) => seeds.map((sd) =>
    [`s1_${t}_w${w}_s${s}_seed${sd}`, ["--tie", t, "--beam", String(w), "--slack", String(s), "--seed", String(sd)]]));
  if (preset === "paper") {
    if (stage === 1) return s1([["all", 2000, 1, range(1, 33)], ["bp", 2000, 1, range(1, 9)], ["sfx", 2000, 0, range(1, 9)],
      ["all", 8000, 1, range(1, 5)], ["all", 2000, 2, range(1, 5)]]);
    const seeds = range(1, 9);
    if (stage === 2) {
      const [T, MS] = small ? [2000, 16] : [300, 12];
      return seeds.map((sd) => [`s2_seed${sd}`, ["--bpbs", String(T), "--maxsplits", String(MS), "--evals", "4", "--tie", "sfx", "--md", "--seed", String(sd)]]);
    }
    if (stage === 3) {
      const R = small ? 40 : 6;
      return seeds.map((sd) => [`s3_seed${sd}`, ["--bpbs", String(R), "--bbeam", "8", "--bopts", "16", "--evals", "2", "--maxsplits", "14",
        "--tie", "sfx", "--md", "--seed", String(sd)]]);
    }
    const [W, IT] = small ? [1000, 300] : [200, 150];
    return seeds.map((sd) => [`s4_seed${sd}`, ["--iters", String(IT), "--fix", "4:28", "--beam", String(W), "--slack", "1", "--tie", "all", "--seed", String(sd)]]);
  }
  const q = preset === "quick";
  const kmin = Math.max(1, Math.round(m / 8)), kmax = Math.max(kmin, Math.round((7 * m) / 8));
  if (stage === 1) {
    if (q) {
      const W = small ? 500 : 200;
      return s1([["all", W, 1, range(1, small ? 5 : 3)], ["bp", W, 1, range(1, small ? 3 : 2)], ["sfx", W, 0, range(1, small ? 3 : 2)]]);
    }
    return s1([["all", 2000, 1, range(1, small ? 9 : 5)], ["bp", 2000, 1, range(1, 3)], ["sfx", 2000, 0, range(1, 3)], ["all", 2000, 2, range(1, 3)]]);
  }
  const seeds = q ? range(1, 3) : range(1, 5);
  if (stage === 2) {
    const T = q ? (small ? 100 : 2) : (small ? 500 : 20), MS = small ? 16 : 12;
    return seeds.map((sd) => [`s2_seed${sd}`, ["--bpbs", String(T), "--maxsplits", String(MS), "--evals", "4", "--tie", "sfx", "--md", "--seed", String(sd)]]);
  }
  if (stage === 3) {
    const [R, B, K, MS] = small ? [q ? 4 : 10, 8, 16, 14] : q ? [1, 2, 4, 6] : [1, 4, 8, 10];
    return seeds.map((sd) => [`s3_seed${sd}`, ["--bpbs", String(R), "--bbeam", String(B), "--bopts", String(K), "--evals", "2",
      "--maxsplits", String(MS), "--tie", "sfx", "--md", "--seed", String(sd)]]);
  }
  const [W, IT] = q ? (small ? [300, 60] : [50, 10]) : (small ? [1000, 150] : [200, 75]);
  return seeds.map((sd) => [`s4_seed${sd}`, ["--iters", String(IT), "--fix", `${kmin}:${kmax}`, "--beam", String(W), "--slack", "1",
    "--tie", "all", "--seed", String(sd)]]);
}

const BEST_RE = /^BEST-(\w+): (\d+) XOR2 gates -> (\d+) XOR2 \+ (\d+) XOR3, umc180 ([\d.]+) GE, stm130 ([\d.]+) GE/;

function parseBest(stdout, saved) {
  const best = {};
  for (const line of stdout.split("\n")) {
    const b = line.match(BEST_RE);
    if (b && saved[b[1]]) best[b[1]] = { gates: +b[2], n2: +b[3], n3: +b[4], umc180: +b[5], stm130: +b[6], text: saved[b[1]] };
  }
  return best;
}

const keyOf = (measure, b) => (measure === "xor2" ? [b.gates, b.umc180] : [b[measure], b.gates]);
const lessKey = (a, b) => a[0] < b[0] || (a[0] === b[0] && a[1] < b[1]);

// The best circuit of each measure over the runs, re-verified (independently of the search program) and written
// in the NIST format; the circuits of the area measures are the maximum XOR3 merge of their XOR2 circuit.
function selectBest(M, D, results, preset) {
  const best = {};
  for (const measure of ["xor2", "umc180", "stm130"]) {
    let r = null;
    for (const x of results) if (x.best[measure] && (!r || lessKey(keyOf(measure, x.best[measure]), keyOf(measure, r.best[measure])))) r = x;
    if (!r) continue;
    const C = readDbsText(r.best[measure].text, M.n).prune();
    let problems = C.check(M, D);
    if (problems.length) throw new CircuitError(`the ${measure} circuit of run ${r.label} does not verify: ${problems.join("; ")}`);
    const target = measure === "xor2" ? C : C.mergeXor3();
    problems = target.check(M, D);
    if (problems.length) throw new CircuitError(`the merged ${measure} circuit of run ${r.label} does not verify: ${problems.join("; ")}`);
    const [n2, n3] = target.counts();
    const extra = measure === "xor2" ? [] : [`XOR2 circuit before merging: ${C.gates.length} gates, depth ${C.depth()}`];
    best[measure] = {
      stage: r.stage, run: r.label, xor2Before: C.gates.length, xor2: n2, xor3: n3, depth: target.depth(),
      umc180: +target.area("umc180").toFixed(2), stm130: +target.area("stm130").toFixed(2), circuit: target, native: C,
      circ: circText(target, M, `${M.name}--${measure}`, `slpd in the browser, ${preset} preset, stage ${r.stage}, run ${r.label}`, extra),
    };
  }
  return best;
}

// ------------------------------------------------------------------ drawing
// The names of circText: inputs x1..xn; a gate that computes an output, the name of its first output; other gates t1, t2, ...
function circuitNames(C) {
  const n = C.n, first = new Map(), outs = new Map();
  C.out.forEach((o, j) => {
    if (!first.has(o)) first.set(o, `y${j + 1}`);
    if (!outs.has(o)) outs.set(o, []);
    outs.get(o).push(`y${j + 1}`);
  });
  const names = Array.from({ length: n }, (_, c) => `x${c + 1}`);
  let t = 0;
  C.gates.forEach((_, k) => names.push(first.has(n + k) ? first.get(n + k) : `t${++t}`));
  return { names, outs };
}

const LEVEL_COLORS = ["#2563eb", "#7c3aed", "#db2777", "#ea580c", "#059669", "#0891b2", "#4f46e5", "#be123c"];
const escXml = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));

// The circuit as an SVG drawing (a string): the inputs on the left, the outputs on the right (in their order), and one
// column per gate level in between (an XOR3 gate, drawn with a double ring, is one level). A wire that skips levels runs
// through a free slot of each column in between (one slot per signal and column, shared by its wires), and the gates of
// each column are ordered to reduce the crossings (barycenter sweeps, as in the layered drawings of Sugiyama et al.).
// Nodes carry data-id (the node index, or p<j> for output j) and wires data-a/data-b, for the page's highlighting.
function circuitSvg(C) {
  const n = C.n, m = C.out.length, N = n + C.gates.length;
  const { val, dep } = C.evaluate();
  const { names } = circuitNames(C);
  const L = Math.max(1, ...dep), P = L + 1;
  const ops = (v) => (v < n ? [] : C.gates[v - n]);
  // layered graph: the N nodes, the slots of the wires that skip levels, then the m output ports (level P)
  const level = [...dep], ins = Array.from({ length: N }, () => []), nxt = Array.from({ length: N }, () => []);
  const slot = new Map(), wires = [];
  const route = (o, to, toLevel) => {
    const chain = [o];
    for (let l = dep[o] + 1; l < toLevel; l++) {
      const key = `${o}:${l}`;
      if (!slot.has(key)) {
        const d = level.length, prev = chain[chain.length - 1];
        slot.set(key, d); level.push(l); ins.push([prev]); nxt.push([]); nxt[prev].push(d);
      }
      chain.push(slot.get(key));
    }
    const last = chain[chain.length - 1];
    ins[to].push(last); nxt[last].push(to); chain.push(to);
    return chain;
  };
  for (let v = n; v < N; v++) for (const o of ops(v)) wires.push([o, v, route(o, v, dep[v])]);
  const port = [];
  C.out.forEach((o, j) => {
    const p = level.length;
    level.push(P); ins.push([]); nxt.push([]); port.push(p);
    wires.push([o, `p${j}`, route(o, p, P)]);
  });
  const K = level.length;
  const cols = Array.from({ length: P + 1 }, () => []);
  for (let u = 0; u < K; u++) cols[level[u]].push(u);
  const maxLen = Math.max(...cols.map((c) => c.length));
  const rowH = Math.max(22, Math.min(40, Math.floor(620 / maxLen)));
  const small = N <= 80;
  const colW = Math.max(140, Math.min(small ? 200 : 240, Math.floor(720 / L)));
  const y = new Array(K);
  const place = (l) => cols[l].forEach((u, i) => { y[u] = (i - (cols[l].length - 1) / 2) * rowH; });
  cols.forEach((_, l) => place(l));
  const mean = (a, self) => (a.length ? a.reduce((s, u) => s + y[u], 0) / a.length : y[self]);
  const sortBy = (l, key) => {
    const k = new Map(cols[l].map((u) => [u, key(u)]));
    cols[l].sort((a, b) => k.get(a) - k.get(b));
    place(l);
  };
  for (let it = 0; it < 8; it++) {
    for (let l = 1; l <= L; l++) sortBy(l, (u) => mean(ins[u], u));
    for (let l = L; l >= 1; l--) sortBy(l, (u) => mean(nxt[u], u));
  }
  for (let l = 1; l <= L; l++) sortBy(l, (u) => mean(ins[u], u));

  const padL = 40, padT = 54, padB = 24, r = small ? 11 : 9, boxW = 38, boxH = small ? 22 : 18;
  const portGap = Math.max(120, Math.round(colW * 0.75));
  const x0 = padL + boxW / 2;
  const W = x0 + L * colW + portGap + boxW / 2 + 24, H = padT + maxLen * rowH + padB;
  const X = (u) => (level[u] <= L ? x0 + level[u] * colW : x0 + L * colW + portGap);
  const Y = (u) => padT + ((maxLen - 1) * rowH) / 2 + y[u] + rowH / 2;
  const color = (l) => LEVEL_COLORS[(l - 1) % LEVEL_COLORS.length];
  const f = (x) => (Math.round(x * 10) / 10).toString();
  const terms = (v) => Array.from({ length: n }, (_, c) => c).filter((c) => (val[v] >> BigInt(c)) & 1n).map((c) => `x${c + 1}`).join(" + ");
  const font = 'font-family="system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif"';
  const s = [`<svg xmlns="http://www.w3.org/2000/svg" width="${f(W)}" height="${f(H)}" viewBox="0 0 ${f(W)} ${f(H)}" role="img" ` +
    `aria-label="Circuit: ${n} inputs, ${m} outputs, ${C.gates.length} gates, depth ${L}">`, `<rect width="${f(W)}" height="${f(H)}" fill="#ffffff"/>`];
  const head = (x, text, fill) => s.push(`<text x="${f(x)}" y="26" text-anchor="middle" font-size="13" font-weight="600" ${font} fill="${fill}">${text}</text>`);
  head(x0, "inputs", "#4b5563");
  for (let l = 1; l <= L; l++) head(x0 + l * colW, `depth ${l}`, color(l));
  head(x0 + L * colW + portGap, "outputs", "#4b5563");
  for (const [o, b, chain] of wires) {
    const toPort = typeof b === "string";
    const pts = chain.map((u, i) => [i === 0 ? X(u) + (u < n ? boxW / 2 : r) : i === chain.length - 1 ? X(u) - (toPort ? boxW / 2 : r) : X(u), Y(u)]);
    let d = `M${f(pts[0][0])},${f(pts[0][1])}`;
    for (let i = 1; i < pts.length; i++) {
      const [x1, y1] = pts[i - 1], [x2, y2] = pts[i], mx = (x1 + x2) / 2;
      d += ` C${f(mx)},${f(y1)} ${f(mx)},${f(y2)} ${f(x2)},${f(y2)}`;
    }
    const stroke = toPort ? "#94a3b8" : color(dep[b]);
    s.push(`<path class="edge" data-a="${o}" data-b="${b}" d="${d}" fill="none" stroke="${stroke}" ` +
      `stroke-width="${small ? 1.8 : 1.3}" stroke-opacity="${toPort ? ".8" : ".55"}"/>`);
  }
  const box = (x, yy, text) => {
    s.push(`<rect x="${f(x - boxW / 2)}" y="${f(yy - boxH / 2)}" width="${boxW}" height="${boxH}" rx="6" fill="#f3f4f6" stroke="#cbd5e1"/>`,
      `<text x="${f(x)}" y="${f(yy + 4.5)}" text-anchor="middle" font-size="${small ? 13 : 11}" ${font} fill="#111827">${text}</text>`);
  };
  for (let v = 0; v < N; v++) {
    const x = X(v), yy = Y(v);
    const tip = v < n ? `${names[v]}: input`
      : `${names[v]} = ${terms(v)}; depth ${dep[v]}` + (ops(v).length === 3 ? "; XOR3" : "") +
        `; inputs of the gate: ${ops(v).map((o) => names[o]).join(", ")}`;
    s.push(`<g class="node" data-id="${v}"><title>${escXml(tip)}</title>`);
    if (v < n) box(x, yy, names[v]);
    else {
      const c = color(dep[v]), p = r - 4;
      s.push(`<circle cx="${f(x)}" cy="${f(yy)}" r="${r}" fill="#ffffff" stroke="${c}" stroke-width="2.2"/>`,
        `<path d="M${f(x - p)},${f(yy)} H${f(x + p)} M${f(x)},${f(yy - p)} V${f(yy + p)}" stroke="${c}" stroke-width="2"/>`);
      if (ops(v).length === 3) s.push(`<circle cx="${f(x)}" cy="${f(yy)}" r="${r + 3.5}" fill="none" stroke="${c}" stroke-width="1.3"/>`);
      if (small) s.push(`<text x="${f(x + r + 2)}" y="${f(yy - r + 1)}" font-size="10.5" ${font} fill="#6b7280">${names[v]}</text>`);
    }
    s.push("</g>");
  }
  C.out.forEach((o, j) => {
    const u = port[j];
    const tip = `y${j + 1} = ${terms(o)}; ` + (o < n ? "an input" : `computed by gate ${names[o]} at depth ${dep[o]}`);
    s.push(`<g class="node port" data-id="p${j}"><title>${escXml(tip)}</title>`);
    box(X(u), Y(u), `y${j + 1}`);
    s.push("</g>");
  });
  s.push("</svg>");
  return s.join("");
}

if (typeof module !== "undefined") {
  module.exports = { Matrix, MatrixError, CircuitError, Circuit, parseMatrices, checkMatrix, readDbsText, circText, presetJobs,
    parseBest, keyOf, selectBest, popcount, circuitNames, circuitSvg, MAX_INPUTS, MAX_OUTPUTS };
}
