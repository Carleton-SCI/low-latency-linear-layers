// Web Worker: runs src/dbs.cpp, compiled to WebAssembly (wasm32-wasi), on one job.
// A minimal WASI shim gives the program its arguments, an in-memory file system under "/" (the matrix, the
// starting circuit of stage 4, the circuits that --save writes), a monotonic clock, stdout and stderr.
// stderr is streamed line by line to the page (the TRACE lines of --trace give the progress of the run).
"use strict";

let module = null;

class Exit {
  constructor(code) { this.code = code; }
}

const E = { SUCCESS: 0, BADF: 8, INVAL: 28, NOENT: 44, NOTDIR: 54 };

function run(job) {
  const enc = new TextEncoder();
  const dec = new TextDecoder();
  const files = new Map(Object.entries(job.files).map(([k, v]) => [k.replace(/^\/+/, ""), enc.encode(v)]));
  const fds = new Map();                     // fd -> {path, pos} for open files; 3 is the preopened "/"
  let nextFd = 4;
  let stdout = "", stderr = "", line = "";
  let mem = null;
  const view = () => new DataView(mem.buffer);
  const bytes = () => new Uint8Array(mem.buffer);
  const args = ["dbs", ...job.args].map((a) => enc.encode(a + "\0"));
  const t0 = performance.now();

  const iovs = (ptr, len) => {
    const v = view(), out = [];
    for (let i = 0; i < len; i++) out.push([v.getUint32(ptr + 8 * i, true), v.getUint32(ptr + 8 * i + 4, true)]);
    return out;
  };
  const errText = (s) => {
    stderr += s;
    line += s;
    let k;
    while ((k = line.indexOf("\n")) >= 0) {
      postMessage({ type: "line", id: job.id, text: line.slice(0, k) });
      line = line.slice(k + 1);
    }
  };

  const wasi = {
    args_sizes_get(argc, size) {
      const v = view();
      v.setUint32(argc, args.length, true);
      v.setUint32(size, args.reduce((s, a) => s + a.length, 0), true);
      return E.SUCCESS;
    },
    args_get(argv, buf) {
      const v = view(), b = bytes();
      for (const a of args) { v.setUint32(argv, buf, true); argv += 4; b.set(a, buf); buf += a.length; }
      return E.SUCCESS;
    },
    environ_sizes_get(count, size) { const v = view(); v.setUint32(count, 0, true); v.setUint32(size, 0, true); return E.SUCCESS; },
    environ_get() { return E.SUCCESS; },
    clock_time_get(id, precision, ptr) {         // the monotonic clock starts with the run
      const ns = id === 0 ? BigInt(Date.now()) * 1000000n : BigInt(Math.round((performance.now() - t0) * 1e6));
      view().setBigUint64(ptr, ns, true);
      return E.SUCCESS;
    },
    fd_prestat_get(fd, ptr) {
      if (fd !== 3) return E.BADF;
      const v = view();
      v.setUint8(ptr, 0);                      // a directory
      v.setUint32(ptr + 4, 1, true);           // its name, "/", has one byte
      return E.SUCCESS;
    },
    fd_prestat_dir_name(fd, ptr, len) {
      if (fd !== 3) return E.BADF;
      bytes().set(enc.encode("/").subarray(0, len), ptr);
      return E.SUCCESS;
    },
    fd_fdstat_get(fd, ptr) {
      const v = view();
      const type = fd <= 2 ? 2 : fd === 3 ? 3 : fds.has(fd) ? 4 : -1;   // character device, directory, file
      if (type < 0) return E.BADF;
      v.setUint8(ptr, type);
      v.setUint16(ptr + 2, 0, true);
      v.setBigUint64(ptr + 8, 0xffffffffffffffffn, true);
      v.setBigUint64(ptr + 16, 0xffffffffffffffffn, true);
      return E.SUCCESS;
    },
    fd_fdstat_set_flags() { return E.SUCCESS; },
    path_open(dirfd, dirflags, pathPtr, pathLen, oflags, rightsBase, rightsInh, fdflags, fdPtr) {
      if (dirfd !== 3) return E.BADF;
      const path = dec.decode(bytes().slice(pathPtr, pathPtr + pathLen)).replace(/^\/+|^\.\//g, "");
      if (oflags & 2) return E.NOTDIR;         // O_DIRECTORY: there are no subdirectories
      if (!files.has(path)) {
        if (!(oflags & 1)) return E.NOENT;     // no O_CREAT
        files.set(path, new Uint8Array(0));
      }
      if (oflags & 8) files.set(path, new Uint8Array(0));   // O_TRUNC
      const fd = nextFd++;
      fds.set(fd, { path, pos: fdflags & 1 ? files.get(path).length : 0 });
      view().setUint32(fdPtr, fd, true);
      return E.SUCCESS;
    },
    fd_read(fd, iovPtr, iovLen, nPtr) {
      let n = 0;
      const f = fds.get(fd);
      if (f) {
        const data = files.get(f.path), b = bytes();
        for (const [p, len] of iovs(iovPtr, iovLen)) {
          const k = Math.max(0, Math.min(len, data.length - f.pos));
          b.set(data.subarray(f.pos, f.pos + k), p);
          f.pos += k; n += k;
          if (k < len) break;
        }
      } else if (fd !== 0) return E.BADF;
      view().setUint32(nPtr, n, true);
      return E.SUCCESS;
    },
    fd_write(fd, iovPtr, iovLen, nPtr) {
      let n = 0;
      const b = bytes(), parts = [];
      for (const [p, len] of iovs(iovPtr, iovLen)) { parts.push(b.slice(p, p + len)); n += len; }
      if (fd === 1) for (const p of parts) stdout += dec.decode(p);
      else if (fd === 2) for (const p of parts) errText(dec.decode(p));
      else {
        const f = fds.get(fd);
        if (!f) return E.BADF;
        for (const p of parts) {
          let data = files.get(f.path);
          if (f.pos + p.length > data.length) {
            const grown = new Uint8Array(Math.max(f.pos + p.length, 2 * data.length));
            grown.set(data);
            data = grown.subarray(0, f.pos + p.length);
            files.set(f.path, data);
          }
          data.set(p, f.pos);
          f.pos += p.length;
        }
      }
      view().setUint32(nPtr, n, true);
      return E.SUCCESS;
    },
    fd_seek(fd, offset, whence, ptr) {
      const f = fds.get(fd);
      if (!f) return fd <= 2 ? E.SUCCESS : E.BADF;
      const size = files.get(f.path).length;
      const base = whence === 0 ? 0 : whence === 1 ? f.pos : size;
      const pos = base + Number(offset);
      if (pos < 0) return E.INVAL;
      f.pos = pos;
      view().setBigUint64(ptr, BigInt(pos), true);
      return E.SUCCESS;
    },
    fd_close(fd) { return fds.delete(fd) || fd <= 3 ? E.SUCCESS : E.BADF; },
    proc_exit(code) { throw new Exit(code); },
  };

  let code = 0, error = null;
  try {
    const inst = new WebAssembly.Instance(module, { wasi_snapshot_preview1: wasi });
    mem = inst.exports.memory;
    inst.exports._start();
  } catch (e) {
    if (e instanceof Exit) code = e.code;
    else { code = -1; error = String(e && e.message || e); }
  }
  if (line) errText("\n");
  const saved = {};
  for (const m of ["xor2", "umc180", "stm130"]) {
    const data = files.get(`run.${m}.txt`);
    if (data) saved[m] = dec.decode(data);
  }
  return { code, error, stdout, stderr, saved, seconds: (performance.now() - t0) / 1000 };
}

onmessage = async (e) => {
  const msg = e.data;
  if (msg.type === "init") {
    module = await WebAssembly.compile(msg.wasm);
    postMessage({ type: "ready" });
  } else if (msg.type === "run") {
    postMessage({ type: "done", id: msg.job.id, result: run(msg.job) });
  }
};
