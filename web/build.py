"""Build web/index.html, the search in the browser, as one self-contained file.

1. web/patch_dbs.py writes a copy of src/dbs.cpp with three changes for wasm32 (src/dbs.cpp is not changed);
2. Zig 0.14.1 (pip install ziglang==0.14.1) compiles it to WebAssembly (wasm32-wasi), and the notices of the
   libraries that it links are written to web/THIRD_PARTY.txt;
3. the page of web/src/ gets the module (base64), the examples (examples/, bench/), the best values of the paper
   (results/summary.json) and the notices inlined.

usage: python3 web/build.py               # compile with Zig, write web/index.html and web/THIRD_PARTY.txt
       python3 web/build.py --wasm FILE   # use a compiled module (and the existing web/THIRD_PARTY.txt)
"""
import argparse
import base64
import csv
import hashlib
import html
import json
import os
import re
import subprocess
import sys
import tempfile

WEB = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WEB)
ZIG_VERSION = "0.14.1"
FLAGS = ["-target", "wasm32-wasi", "-O2", "-std=c++17", "-fno-exceptions", "-s", "-Wl,-z,stack-size=8388608"]


def read(*p):
    with open(os.path.join(*p), encoding="utf-8") as f:
        return f.read()


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def compile_wasm(zig, tmp):
    src = os.path.join(tmp, "dbs_web.cpp")
    subprocess.run([sys.executable, os.path.join(WEB, "patch_dbs.py"), os.path.join(ROOT, "src", "dbs.cpp"), src], check=True)
    env = json.loads(subprocess.run(zig + ["env"], check=True, capture_output=True, text=True).stdout)
    if env["version"] != ZIG_VERSION:
        sys.exit(f"error: Zig {env['version']}; the build expects Zig {ZIG_VERSION} (pip install ziglang=={ZIG_VERSION})")
    out = os.path.join(tmp, "dbs.wasm")
    subprocess.run(zig + ["c++"] + FLAGS + ["-o", out, src], check=True)
    with open(out, "rb") as f:
        return f.read(), env["lib_dir"]


def notices(lib):
    """The notices of the libraries that the module links (see web/README.md)."""
    def part(title, path):
        return f"==== {title}\n\n{read(lib, *path.split('/')).strip()}\n"
    llvm = read(lib, "libcxx", "LICENSE.TXT")
    i = llvm.index("---- LLVM Exceptions to the Apache 2.0 License ----")
    exc = llvm[i:llvm.index("=====", i)].strip()
    emm = read(lib, "libc", "wasi", "emmalloc", "emmalloc.c")
    emm = emm[emm.index(" * Copyright"):emm.index("found in the LICENSE file.") + 26].replace(" * ", "").replace(" *", "")
    mit = read(lib, "libc", "wasi", "LICENSE-MIT").strip()
    return "\n".join([
        "Third-party code in the WebAssembly module of web/index.html",
        "",
        f"web/build.py compiles src/dbs.cpp to WebAssembly (wasm32-wasi) with Zig {ZIG_VERSION}, which links it with",
        "libraries that come with Zig. The module contains parts of:",
        "",
        "- wasi-libc, the C library (Apache License 2.0 with LLVM Exceptions, Apache License 2.0, or MIT), with",
        "  code from musl (MIT), cloudlibc (BSD-2-Clause) and emmalloc of Emscripten (MIT);",
        "- libc++ and libc++abi of the LLVM Project (Apache License 2.0 with LLVM Exceptions);",
        "- compiler-rt of Zig (MIT; the 128-bit shifts and multiplication).",
        "",
        "Their notices follow.",
        "",
        part("wasi-libc: LICENSE", "libc/wasi/LICENSE"),
        part("wasi-libc: LICENSE-MIT", "libc/wasi/LICENSE-MIT"),
        part("musl: COPYRIGHT", "libc/wasi/libc-top-half/musl/COPYRIGHT"),
        part("cloudlibc: LICENSE", "libc/wasi/libc-bottom-half/cloudlibc/LICENSE"),
        f"==== emmalloc (Emscripten): the notice of emmalloc.c, and the MIT license\n\n{emm}\n\n{mit}\n",
        part("Zig (compiler-rt): LICENSE", "../LICENSE"),
        f"==== libc++ and libc++abi (LLVM Project): Apache License 2.0 with LLVM Exceptions; the exceptions\n\n{exc}\n",
    ])


def examples():
    """The matrices offered by the page, with the best values of the paper for the benchmark matrices."""
    out = []
    toy = read(ROOT, "examples", "toy.mat.txt")
    out.append(dict(group="Small examples", label="Toy example, 4 × 6", name="toy",
                    text="# examples/toy.mat.txt\n" + toy))
    rows = [" ".join("1" if (c - j) % 8 < 5 else "0" for c in range(8)) for j in range(8)]
    out.append(dict(group="Small examples", label="Example of Sec. 3 of the paper, 8 × 8", name="example8x8",
                    text="# The 8 x 8 example of Sec. 3 of the paper: y_j = x_j + x_(j+1) + ... + x_(j+4), indices mod 8\n"
                    + "\n".join(rows) + "\n"))
    best = json.loads(read(ROOT, "results", "summary.json"))
    with open(os.path.join(ROOT, "bench", "benchmark.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            name = r["file"][:-4]
            size = r["size"].replace("x", " × ")
            b = best[name]["best"]
            x2, umc = b["xor2"]["xor2"] + b["xor2"]["xor3"], b["umc180"]["umc180"]
            out.append(dict(group="Benchmark of the paper (bench/)", label=f"{r['label']}, {size}", name=name,
                            text=f"# {r['label']}, {r['size']} (bench/{r['file']})\n" + read(ROOT, "bench", r["file"]),
                            ref=f"Best circuits of the paper (results/): {x2} XOR2 gates; {umc:.2f} GE with the UMC 180 nm "
                                "unit areas, after merging.",
                            refs=dict(xor2=f"{x2} XOR2 gates", umc180=f"{umc:.2f} GE")))
    out.append(dict(group="Other", label="AES InvMixColumns, 32 × 32", name="aes_invmixcolumns",
                    text="# examples/aes_invmixcolumns.mat.txt\n" + read(ROOT, "examples", "aes_invmixcolumns.mat.txt")))
    return out


def js(value):
    return json.dumps(value, ensure_ascii=False).replace("</", "<\\/")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--wasm", help="use this compiled module instead of compiling src/dbs.cpp")
    ap.add_argument("--zig", default=f"{sys.executable} -m ziglang", help="command that runs Zig (default: %(default)s)")
    ap.add_argument("-o", "--out", default=os.path.join(WEB, "index.html"))
    a = ap.parse_args()
    third = os.path.join(WEB, "THIRD_PARTY.txt")
    if a.wasm:
        with open(a.wasm, "rb") as f:
            wasm = f.read()
    else:
        with tempfile.TemporaryDirectory() as tmp:
            wasm, lib = compile_wasm(a.zig.split(), tmp)
        with open(third, "w", encoding="utf-8") as f:
            f.write(notices(lib))
    with open(os.path.join(ROOT, "src", "dbs.cpp"), "rb") as f:
        dbs_hash = sha256(f.read())
    info = (f"Search program: src/dbs.cpp (SHA-256 {dbs_hash}) with the three changes of web/patch_dbs.py, compiled to "
            f"WebAssembly by web/build.py with Zig {ZIG_VERSION} (clang, libc++, wasi-libc); module of {len(wasm):,} bytes, "
            f"SHA-256 {sha256(wasm)}.")
    parts = {
        "STYLE": read(WEB, "src", "style.css"),
        "CORE": read(WEB, "src", "core.js"),
        "WASM": base64.b64encode(wasm).decode(),
        "WORKER": js(read(WEB, "src", "worker.js")),
        "EXAMPLES": js(examples()),
        "BUILD_INFO": js(info),
        "APP": read(WEB, "src", "app.js"),
        "NOTICES": html.escape(read(third)),
    }
    for k in ("CORE", "APP", "STYLE"):
        assert "</script" not in parts[k].lower() and "</style" not in parts[k].lower(), k
    page = read(WEB, "src", "index.html")
    found = re.findall(r"/\*([A-Z_]+)\*/|<!--([A-Z_]+)-->", page)
    assert sorted(x or y for x, y in found) == sorted(parts), found
    page = re.sub(r"/\*([A-Z_]+)\*/|<!--([A-Z_]+)-->", lambda m: parts[m.group(1) or m.group(2)], page)
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(page)
    print(f"{os.path.relpath(a.out)}: {len(page.encode()):,} bytes; module {len(wasm):,} bytes, SHA-256 {sha256(wasm)}")


if __name__ == "__main__":
    main()
