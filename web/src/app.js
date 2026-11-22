// The page: the matrix and its checks, the stages of slpd/pipeline.py on a pool of Web Workers that run
// src/dbs.cpp compiled to WebAssembly (worker.js), the progress (the TRACE lines of --trace), and the circuits,
// re-verified by core.js before they are shown.
"use strict";
(function () {
  const $ = (id) => document.getElementById(id);
  const STAGES = { 1: "beam search over the near-ties", 2: "backward splitting", 3: "beam search over the splits",
    4: "large neighborhood search" };
  const MEASURES = [["xor2", "Fewest XOR2 gates"], ["umc180", "Smallest area, UMC 180 nm"], ["stm130", "Smallest area, STM 130 nm"]];
  const TRACE_RE = /^TRACE cpu=[\d.]+ measure=(\w+) xor2=(\d+) umc180=([\d.]+)/;
  const wasmBytes = Uint8Array.from(atob(WASM_BASE64), (c) => c.charCodeAt(0));
  const cores = Math.max(1, Math.min(32, navigator.hardwareConcurrency || 1));

  const state = { base: "matrix", example: null, matrices: [], M: null, info: null, key: "", search: null, last: null, status: "" };

  // ---------------------------------------------------------------- small helpers
  function el(tag, props = {}, ...kids) {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(props)) {
      if (k === "class") e.className = v;
      else if (k === "text") e.textContent = v;
      else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
      else e.setAttribute(k, v);
    }
    for (const k of kids) if (k !== null && k !== undefined) e.append(k);
    return e;
  }
  const msg = (kind, text) => el("div", { class: `msg ${kind}`, text });
  const plural = (k, one, many = one + "s") => `${k} ${k === 1 ? one : many}`;
  function clock(seconds) {
    const s = Math.floor(seconds), h = Math.floor(s / 3600), m = Math.floor((s % 3600) / 60), r = s % 60;
    return h ? `${h}:${String(m).padStart(2, "0")}:${String(r).padStart(2, "0")}` : `${m}:${String(r).padStart(2, "0")}`;
  }
  function download(name, data, type = "text/plain") {
    const url = URL.createObjectURL(data instanceof Blob ? data : new Blob([data], { type }));
    const a = el("a", { href: url, download: name });
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 10000);
  }
  const safeName = (s) => s.replace(/\.(mat\.txt|txt|mat|csv)$/i, "").replace(/[^A-Za-z0-9_.-]+/g, "_").replace(/^[._]+/, "") || "matrix";

  // A ZIP archive without compression (method 0), enough for a few text files.
  const CRC = Array.from({ length: 256 }, (_, n) => { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; return c >>> 0; });
  function zip(files) {
    const enc = new TextEncoder(), parts = [], central = [], d = new Date();
    const time = (d.getHours() << 11) | (d.getMinutes() << 5) | (d.getSeconds() >> 1);
    const date = ((d.getFullYear() - 1980) << 9) | ((d.getMonth() + 1) << 5) | d.getDate();
    let offset = 0;
    for (const f of files) {
      const name = enc.encode(f.name), data = enc.encode(f.text);
      let crc = 0xffffffff;
      for (const b of data) crc = CRC[(crc ^ b) & 0xff] ^ (crc >>> 8);
      crc = (crc ^ 0xffffffff) >>> 0;
      const h = new DataView(new ArrayBuffer(30)), c = new DataView(new ArrayBuffer(46));
      [[0, 0x04034b50, 4], [4, 20, 2], [6, 0x0800, 2], [8, 0, 2], [10, time, 2], [12, date, 2], [14, crc, 4], [18, data.length, 4],
        [22, data.length, 4], [26, name.length, 2], [28, 0, 2]].forEach(([o, v, w]) => (w === 4 ? h.setUint32(o, v, true) : h.setUint16(o, v, true)));
      [[0, 0x02014b50, 4], [4, 20, 2], [6, 20, 2], [8, 0x0800, 2], [10, 0, 2], [12, time, 2], [14, date, 2], [16, crc, 4],
        [20, data.length, 4], [24, data.length, 4], [28, name.length, 2], [30, 0, 2], [32, 0, 2], [34, 0, 2], [36, 0, 2], [38, 0, 4],
        [42, offset, 4]].forEach(([o, v, w]) => (w === 4 ? c.setUint32(o, v, true) : c.setUint16(o, v, true)));
      parts.push(new Uint8Array(h.buffer), name, data);
      central.push(new Uint8Array(c.buffer), name);
      offset += 30 + name.length + data.length;
    }
    const size = central.reduce((s, p) => s + p.length, 0), e = new DataView(new ArrayBuffer(22));
    [[0, 0x06054b50, 4], [4, 0, 2], [6, 0, 2], [8, files.length, 2], [10, files.length, 2], [12, size, 4], [16, offset, 4], [20, 0, 2]]
      .forEach(([o, v, w]) => (w === 4 ? e.setUint32(o, v, true) : e.setUint16(o, v, true)));
    return new Blob([...parts, ...central, new Uint8Array(e.buffer)], { type: "application/zip" });
  }

  // ---------------------------------------------------------------- the matrix
  function setMatrixText(text, base, example = null) {
    $("matrix").value = text;
    state.base = base;
    state.example = example ? { ...example, text } : null;
    $("pick").replaceChildren();
    onMatrix();
  }

  function onMatrix() {
    const box = $("check"), text = $("matrix").value;
    box.replaceChildren();
    state.M = null;
    state.info = null;
    if (state.example && state.example.text !== text) state.example = null;
    if (!text.trim()) {
      $("pick-wrap").hidden = true;
      box.append(el("div", { class: "muted small", text: "Paste a matrix, open a file or choose an example." }));
      return update();
    }
    try {
      state.matrices = parseMatrices(text, state.base);
    } catch (e) {
      if (!(e instanceof MatrixError)) throw e;
      $("pick-wrap").hidden = true;
      box.append(msg("err", e.message));
      return update();
    }
    const pick = $("pick"), k = state.matrices.length;
    if (k > 1) {
      if (pick.options.length !== k) {
        pick.replaceChildren(...state.matrices.map((M, i) => el("option", { value: String(i), text: `${i + 1} of ${k}: ${M.m} × ${M.n}` })));
      }
      $("pick-wrap").hidden = false;
    } else $("pick-wrap").hidden = true;
    const M = state.matrices[k > 1 ? Number(pick.value || 0) : 0];
    const { errors, warnings, info } = checkMatrix(M);
    const facts = el("div", { class: "facts" },
      el("span", {}, el("b", { text: `${info.rows} × ${info.columns}` }), " (outputs × inputs)"),
      el("span", { text: `${info.ones} ones` }),
      el("span", { text: `row weights ${info.minWeight}–${info.maxWeight}` }),
      el("span", { text: `rank ${info.rank}${info.invertible ? ", invertible" : ""}` }));
    if (info.minDepth !== null) {
      facts.append(el("span", {}, "minimum depth ", el("b", { text: String(info.minDepth) })),
        el("span", { text: `${info.naiveXor2} XOR2 gates without sharing` }));
    }
    box.append(facts);
    for (const e of errors) box.append(msg("err", e));
    for (const w of warnings) box.append(msg("warn", w));
    if (state.example && state.example.ref) box.append(el("div", { class: "muted small", text: state.example.ref }));
    if (!errors.length) {
      state.M = M;
      state.info = info;
      const key = `${M.n}|${M.rowsText()}`;
      if (key !== state.key) {
        state.key = key;
        state.status = "";
        $("depth").value = String(info.minDepth);
      }
      $("depth").min = String(info.minDepth);
    }
    update();
  }

  // ---------------------------------------------------------------- settings
  const preset = () => document.querySelector("input[name=preset]:checked").value;
  const chosenStages = () => [...document.querySelectorAll("input[name=stage]:checked")].map((x) => Number(x.value));

  // What the run would do, and why it cannot start: the same rules as slpd/pipeline.py (run_search).
  function plan() {
    const p = { problems: [], notes: [], stages: [], runs: 0 };
    const M = state.M;
    if (!M) { p.problems.push(state.matrices.length || $("matrix").value.trim() ? "fix the matrix first" : "no matrix"); return p; }
    const s = state.info.minDepth, D = Number($("depth").value);
    if (!Number.isInteger(D) || D < 1) p.problems.push("the depth bound must be a positive integer");
    else if (D < s) p.problems.push(`the depth bound ${D} is below the minimum depth ${s}`);
    let stages = chosenStages();
    if (!stages.length) p.problems.push("choose at least one stage");
    if (Number.isInteger(D) && D > s && stages.some((st) => st > 1)) {
      p.notes.push(`stages 2–4 assume the minimum depth (${s}); with D = ${D} only stage 1 runs`);
      stages = [1];
    }
    p.D = D;
    p.stages = stages.map((st) => ({ st, jobs: presetJobs(preset(), M.n, M.m, st) }));
    p.runs = p.stages.reduce((a, x) => a + x.jobs.length, 0);
    return p;
  }

  function update() {
    const running = !!state.search;
    const p = plan();
    const s = state.info ? state.info.minDepth : null, D = Number($("depth").value);
    $("depth-note").textContent = s === null ? "" : D < s ? `below the minimum depth, ${s}` : D > s ? `minimum depth: ${s}; only stage 1 runs` : "the minimum depth";
    $("depth-note").className = s !== null && D < s ? "err-text" : "";
    const w = Number($("workers").value);
    $("workers-note").textContent = `this computer reports ${plural(cores, "logical processor")}` + (w > cores ? "; more runs than processors only share them" : "");
    const pr = preset();
    let note = { quick: "Smallest budget, to try the search.", standard: "Larger budgets, same stages.",
      paper: "The configuration of the paper (its 16×16 settings for at most 16 inputs, its 32×32 settings otherwise). " +
        "On the server of the paper, it took about 8.5 core-hours per 16×16 matrix and about 150 core-hours per 32×32 matrix of " +
        "depth 3 (Table 7 of the paper)." }[pr];
    if (pr !== "paper") note += " It does not reproduce the paper.";
    if (state.M && p.stages.length) {
      note += ` ${plural(p.runs, "run")}: ` + p.stages.map((x) => `stage ${x.st}: ${x.jobs.length}`).join(", ") + ".";
    }
    $("preset-note").textContent = note;
    $("run").disabled = running || p.problems.length > 0;
    $("stop").disabled = !running;
    for (const id of ["matrix", "example", "open-file", "clear", "pick", "depth", "workers"]) $(id).disabled = running;
    for (const x of document.querySelectorAll("input[name=preset], input[name=stage]")) x.disabled = running;
    if (!running) {
      $("run-note").textContent = p.problems.length ? `Cannot run: ${p.problems[0]}.` : p.notes.length ? `Note: ${p.notes[0]}.` : state.status;
    }
  }

  // ---------------------------------------------------------------- the workers
  // A pool of Web Workers; run(job) resolves with the result of the worker (see worker.js), or with
  // {stopped: true} when the pool is terminated first.
  function makePool(size, onLine) {
    const url = URL.createObjectURL(new Blob([WORKER_SOURCE], { type: "text/javascript" }));
    const workers = [], idle = [], queue = [], pending = new Map();
    let nextId = 1, closed = false;
    const pump = () => {
      while (!closed && idle.length && queue.length) {
        const w = idle.pop(), [job, resolve] = queue.shift();
        w.job = job.id;
        pending.set(job.id, resolve);
        w.postMessage({ type: "run", job });
      }
    };
    const ready = Promise.all(Array.from({ length: size }, () => new Promise((resolve, reject) => {
      const w = new Worker(url);
      w.onmessage = (e) => {
        const m = e.data;
        if (m.type === "ready") { idle.push(w); resolve(); pump(); }
        else if (m.type === "line") onLine(m.id, m.text);
        else if (m.type === "done") {
          const done = pending.get(m.id);
          pending.delete(m.id);
          w.job = null;
          idle.push(w);
          if (done) done(m.result);
          pump();
        }
      };
      w.onerror = (e) => {
        e.preventDefault();
        const text = e.message || "the worker failed";
        if (w.job === undefined) return reject(new Error(text));
        const done = pending.get(w.job);
        pending.delete(w.job);
        w.job = null;
        idle.push(w);
        if (done) done({ code: -1, error: text, stdout: "", stderr: "", saved: {}, seconds: 0 });
        pump();
      };
      w.postMessage({ type: "init", wasm: wasmBytes });
      workers.push(w);
    })));
    return {
      ready,
      run(args, files) {
        const id = nextId++;
        const job = { id, args, files };
        return { id, result: new Promise((resolve) => { if (closed) resolve({ stopped: true }); else { queue.push([job, resolve]); pump(); } }) };
      },
      terminate() {
        closed = true;
        for (const w of workers) w.terminate();
        for (const [, resolve] of queue) resolve({ stopped: true });
        for (const resolve of pending.values()) resolve({ stopped: true });
        queue.length = 0;
        pending.clear();
        URL.revokeObjectURL(url);
      },
    };
  }

  // ---------------------------------------------------------------- progress
  const logLines = [];
  let paintQueued = false;
  function paint() {
    if (paintQueued) return;
    paintQueued = true;
    requestAnimationFrame(() => {
      paintQueued = false;
      const S = state.search || state.last;
      if (!S) return;
      $("live-xor2").textContent = S.liveXor2 === null ? "–" : String(S.liveXor2);
      $("live-umc").textContent = S.liveUmc === null ? "–" : `${S.liveUmc.toFixed(2)} GE`;
      $("live-runs").textContent = `${S.done} / ${S.total}`;
      const pct = S.total ? Math.round((100 * S.done) / S.total) : 0;
      $("bar").style.width = `${pct}%`;
      $("bar").parentElement.setAttribute("aria-valuenow", String(pct));
      for (const row of S.rows) {
        const x = S.stage[row.st];
        row.tr.className = x.status;
        row.runs.textContent = x.status === "skipped" ? "skipped" : `${x.done} / ${x.total}` + (x.failed ? ` (${x.failed} failed)` : "");
        row.best.textContent = x.best === null ? "–" : String(x.best);
        row.time.textContent = x.t0 === null ? "–" : clock(((x.t1 ?? performance.now()) - x.t0) / 1000);
      }
      const pre = $("log"), atEnd = pre.scrollTop + pre.clientHeight >= pre.scrollHeight - 4;
      pre.textContent = logLines.join("\n");
      if (atEnd) pre.scrollTop = pre.scrollHeight;
    });
  }
  function log(text) {
    logLines.push(text);
    if (logLines.length > 4000) logLines.splice(0, logLines.length - 4000);
    paint();
  }

  // ---------------------------------------------------------------- the search
  async function runSearch() {
    const p = plan();
    if (p.problems.length || state.search) return;
    const M = state.M, D = p.D, pr = preset(), example = state.example;
    const workers = Math.max(1, Math.min(32, Math.round(Number($("workers").value)) || cores));
    const S = {
      M, D, preset: pr, example, stopped: false, results: [], perStage: {}, t0: performance.now(), done: 0, total: p.runs,
      liveXor2: null, liveUmc: null, jobs: new Map(), stage: {}, rows: [],
    };
    state.search = S;
    logLines.length = 0;
    $("results-card").hidden = true;
    $("progress-card").hidden = false;
    const tbody = $("stage-rows");
    tbody.replaceChildren();
    for (const { st, jobs } of p.stages) {
      S.stage[st] = { status: "", done: 0, total: jobs.length, failed: 0, best: null, t0: null, t1: null };
      const row = { st, tr: el("tr"), runs: el("td"), best: el("td"), time: el("td") };
      row.tr.append(el("td", { text: `${st} · ${STAGES[st]}` }), row.runs, row.best, row.time);
      tbody.append(row.tr);
      S.rows.push(row);
    }
    log(`${M.name}: ${M.m} × ${M.n}, depth bound ${D}, ${pr} preset, ${plural(p.runs, "run")} on ${plural(workers, "worker")}`);
    for (const n of p.notes) log(`note: ${n}`);
    $("run-note").textContent = "Running… Closing or reloading the page stops the search.";
    update();
    const timer = setInterval(() => { $("elapsed").textContent = clock((performance.now() - S.t0) / 1000); paint(); }, 500);
    $("elapsed").textContent = "0:00";
    paint();

    const onLine = (id, text) => {
      const job = S.jobs.get(id);
      if (!job) return;
      const t = text.match(TRACE_RE);
      if (!t) return log(`[${job.label}] ${text}`);
      const x2 = Number(t[2]), umc = Number(t[3]), st = S.stage[job.stage];
      if (t[1] === "xor2") {
        if (S.liveXor2 === null || x2 < S.liveXor2) S.liveXor2 = x2;
        if (st.best === null || x2 < st.best) st.best = x2;
      } else if (S.liveUmc === null || umc < S.liveUmc) S.liveUmc = umc;
      paint();
    };
    let pool = null;
    try {
      const most = Math.max(...p.stages.map((x) => x.jobs.length));
      pool = makePool(Math.min(workers, most), onLine);
      S.pool = pool;
      await pool.ready;
      const rowsFile = M.rowsText();
      for (const { st, jobs } of p.stages) {
        if (S.stopped) break;
        const files = { "m.txt": rowsFile };
        let todo = jobs;
        if (st === 4) {
          const done = S.results.filter((r) => r.best.xor2);
          if (!done.length) {
            log("stage 4 skipped: no circuit from stages 1-3");
            S.stage[4].status = "skipped";
            S.total -= jobs.length;
            continue;
          }
          const inc = done.reduce((a, b) => (lessKey(keyOf("xor2", b.best.xor2), keyOf("xor2", a.best.xor2)) ? b : a));
          files["start.txt"] = inc.best.xor2.text;
          log(`stage 4 starts from the ${inc.best.xor2.gates}-gate circuit of run ${inc.label}`);
          todo = jobs.map(([label, args]) => [label, ["--lns", "/start.txt", ...args]]);
        }
        const X = S.stage[st];
        X.status = "running";
        X.t0 = performance.now();
        log(`stage ${st}: ${plural(todo.length, "run")} on ${plural(Math.min(workers, todo.length), "worker")}`);
        const runs = todo.map(([label, args]) => {
          const r = pool.run(["/m.txt", String(D), ...args, "--save", "/run", "--trace"], files);
          S.jobs.set(r.id, { label, stage: st });
          return r.result.then((res) => {
            if (res.stopped) return null;
            const ok = res.code === 0 && !res.error;
            const best = ok ? parseBest(res.stdout, res.saved) : {};
            const out = { label, stage: st, args, ok, seconds: res.seconds, best };
            X.done++;
            S.done++;
            if (!ok) {
              X.failed++;
              log(`[${label}] failed (${res.error || `exit ${res.code}`})${res.stderr ? ": " + res.stderr.trim().split("\n").pop() : ""}`);
            } else {
              const b = best.xor2;
              log(`[${label}] ${res.seconds.toFixed(1)} s` + (b ? `: fewest XOR2 gates ${b.gates}` : ": no circuit") +
                (best.umc180 ? `, smallest UMC area ${best.umc180.umc180.toFixed(2)} GE` : ""));
            }
            paint();
            return out;
          });
        });
        const res = (await Promise.all(runs)).filter(Boolean);
        S.results.push(...res);
        X.t1 = performance.now();
        const got = res.filter((r) => r.best.xor2).map((r) => r.best.xor2.gates);
        X.best = got.length ? Math.min(...got) : null;      // the runs that finished
        X.status = S.stopped ? "stopped" : "done";
        S.perStage[st] = { runs: res.length, failed: res.filter((r) => !r.ok).length, seconds: Math.round((X.t1 - X.t0) / 100) / 10,
          best_xor2: got.length ? Math.min(...got) : null };
        if (!S.stopped) log(`  done in ${clock((X.t1 - X.t0) / 1000)}; best XOR2 count ${S.perStage[st].best_xor2}`);
        paint();
      }
    } catch (e) {
      log(`error: ${e.message || e}`);
      S.error = e.message || String(e);
    } finally {
      if (pool) pool.terminate();
      clearInterval(timer);
    }
    S.seconds = (performance.now() - S.t0) / 1000;
    $("elapsed").textContent = clock(S.seconds);
    finish(S);
  }

  function finish(S) {
    state.search = null;
    state.last = S;
    const done = S.results.filter((r) => Object.keys(r.best).length);
    let best = null, error = S.error || null;
    if (!error && done.length) {
      try { best = selectBest(S.M, S.D, done, S.preset); } catch (e) { error = e.message || String(e); }
    }
    const runs = `${S.done} of ${plural(S.total, "run")}`;
    if (error) state.status = `Error: ${error}`;
    else if (!best) state.status = S.stopped ? `Stopped after ${runs}; no circuit yet.` : "No run produced a circuit; see the log.";
    else state.status = S.stopped ? `Stopped after ${runs}, ${clock(S.seconds)}.` : `Done: ${plural(S.done, "run")}, ${clock(S.seconds)}.`;
    if (error) log(`error: ${error}`);
    S.liveXor2 = best && best.xor2 ? best.xor2.xor2 : null;     // the runs that finished (all of them, unless stopped)
    S.liveUmc = best && best.umc180 ? best.umc180.umc180 : null;
    if (best) showResults(S, best);
    update();
    paint();
  }

  // ---------------------------------------------------------------- results
  function showResults(S, best) {
    const M = S.M, files = [], cards = [];
    const summary = {
      matrix: M.name, rows: M.m, columns: M.n, min_depth: M.minDepth(), depth_bound: S.D, preset: S.preset, stages: {},
      seconds: Math.round(S.seconds * 10) / 10, best: {}, program: "src/dbs.cpp, WebAssembly build (web/index.html)",
    };
    if (S.stopped) summary.stopped = `after ${S.done} of ${S.total} runs`;
    for (const [st, x] of Object.entries(S.perStage)) summary.stages[st] = x;
    const ref = S.example && S.example.refs;
    for (const [measure, title] of MEASURES) {
      const b = best[measure];
      if (!b) continue;
      const file = `${M.name}.${measure}.circ.txt`;
      files.push({ name: file, text: b.circ });
      summary.best[measure] = { file, stage: b.stage, run: b.run, xor2_gates_before_merge: b.xor2Before, xor2: b.xor2, xor3: b.xor3,
        depth: b.depth, umc180: b.umc180, stm130: b.stm130, verified: true };
      const big = measure === "xor2" ? `${b.xor2} XOR2 gates` : `${b[measure].toFixed(2)} GE`;
      const lines = measure === "xor2"
        ? [`depth ${b.depth}`]
        : [`${b.xor2} XOR2 + ${b.xor3} XOR3, depth ${b.depth}`, `merged from ${b.xor2Before} XOR2 gates`];
      if (measure === "umc180") lines.push(`STM 130 nm: ${b.stm130.toFixed(2)} GE`);
      if (measure === "stm130") lines.push(`UMC 180 nm: ${b.umc180.toFixed(2)} GE`);
      lines.push(`stage ${b.stage}, run ${b.run}`);
      if (ref && ref[measure]) lines.push(`best of the paper: ${ref[measure]}`);
      const pre = el("pre", { hidden: "" });
      pre.textContent = b.circ;
      const show = el("button", { type: "button", text: "Show", onclick: () => { pre.hidden = !pre.hidden; show.textContent = pre.hidden ? "Show" : "Hide"; } });
      const copy = el("button", { type: "button", text: "Copy", onclick: async () => {
        try { await navigator.clipboard.writeText(b.circ); copy.textContent = "Copied"; } catch (e) { copy.textContent = "Copy failed"; }
        setTimeout(() => { copy.textContent = "Copy"; }, 1500);
      } });
      cards.push(el("div", { class: "result" },
        el("h3", { text: title }),
        el("div", { class: "big", text: big }),
        ...lines.map((t) => el("div", { class: "sub", text: t })),
        el("div", { class: "verified", text: `✓ re-verified: computes the matrix, depth ≤ ${S.D}` }),
        el("div", { class: "btns" },
          el("button", { type: "button", text: "Download", title: file, onclick: () => download(file, b.circ) }), copy, show),
        pre));
    }
    const sumFile = `${M.name}.summary.json`;
    const sumText = JSON.stringify(summary, null, 1) + "\n";
    files.push({ name: sumFile, text: sumText });
    $("results").replaceChildren(...cards);
    $("results-note").textContent = "Circuit format of the NIST Circuit Complexity Project (XOR3 out a b c: a 3-input gate whose " +
      "last two operands form the absorbed XOR2 gate). The STM areas use the same merge as the UMC areas (no duplication of gates)." +
      (S.stopped ? " The search was stopped: these are the best circuits of the runs that finished." : "");
    $("download-all").onclick = () => download(`${M.name}.circuits.zip`, zip(files));
    $("results-card").hidden = false;
    state.files = files;
  }

  // ---------------------------------------------------------------- wiring
  function init() {
    const sel = $("example"), groups = new Map();
    EXAMPLES.forEach((x, i) => {
      if (!groups.has(x.group)) { const g = el("optgroup", { label: x.group }); groups.set(x.group, g); sel.append(g); }
      groups.get(x.group).append(el("option", { value: String(i), text: x.label }));
    });
    sel.addEventListener("change", () => {
      const x = EXAMPLES[Number(sel.value)];
      sel.value = "";
      if (x) setMatrixText(x.text, x.name, x);
    });
    $("open-file").addEventListener("click", () => $("file").click());
    const readFile = (f) => {
      if (!f) return;
      if (f.size > 1 << 20) { $("check").replaceChildren(msg("err", `${f.name}: ${f.size} bytes; a matrix file of at most 1 MB is expected`)); return; }
      f.text().then((t) => setMatrixText(t, safeName(f.name)));
    };
    $("file").addEventListener("change", (e) => { readFile(e.target.files[0]); e.target.value = ""; });
    const ta = $("matrix");
    ta.addEventListener("dragover", (e) => { e.preventDefault(); ta.classList.add("drag"); });
    ta.addEventListener("dragleave", () => ta.classList.remove("drag"));
    ta.addEventListener("drop", (e) => { e.preventDefault(); ta.classList.remove("drag"); readFile(e.dataTransfer.files[0]); });
    let timer = null;
    ta.addEventListener("input", () => { clearTimeout(timer); timer = setTimeout(onMatrix, 250); });
    ta.addEventListener("keydown", (e) => { if (e.key === "Enter" && (e.ctrlKey || e.metaKey)) { e.preventDefault(); onMatrix(); runSearch(); } });
    $("clear").addEventListener("click", () => { setMatrixText("", "matrix"); ta.focus(); });
    $("pick").addEventListener("change", onMatrix);
    $("workers").value = String(cores);
    for (const id of ["depth", "workers"]) $(id).addEventListener("input", update);
    for (const x of document.querySelectorAll("input[name=preset], input[name=stage]")) x.addEventListener("change", update);
    $("run").addEventListener("click", runSearch);
    $("stop").addEventListener("click", () => {
      const S = state.search;
      if (!S || S.stopped) return;
      S.stopped = true;
      log("stopped by the user");
      if (S.pool) S.pool.terminate();
    });
    window.addEventListener("beforeunload", (e) => { if (state.search) { e.preventDefault(); e.returnValue = ""; } });
    $("build-info").textContent = BUILD_INFO;
    onMatrix();
  }

  init();
})();
