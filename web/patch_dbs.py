"""Write the copy of src/dbs.cpp that web/build.py compiles to WebAssembly (wasm32-wasi), and print what changes.
src/dbs.cpp itself is not changed (its SHA-256 is recorded with the server runs of the paper).

Three changes, needed on wasm32; on x86-64 Linux, the changed copy gives the same output as src/dbs.cpp
(tests/test_web.py):
  1. `long` (32 bits on wasm32, 64 bits on x86-64 Linux) becomes `long long`, with the literals and printf formats;
  2. print_net() writes to a FILE* instead of assigning to `stdout` (a constant in the C library of WASI);
  3. std::clock() becomes a monotonic clock (the browser gives no CPU time); it times only the diagnostics on
     stderr (--progress, --trace), not the search.
usage: python3 web/patch_dbs.py src/dbs.cpp OUT.cpp"""
import re
import sys

if len(sys.argv) != 3:
    sys.exit(__doc__)
src = open(sys.argv[1]).read()
t = re.sub(r"(?<!long )(?<!unsigned )\blong\b(?! long)", "long long", src)
t = re.sub(r"\b(\d+)L\b", r"\1LL", t)
t = t.replace("%ld", "%lld")
i = t.index("static void print_net(const Net& n, bool xor3) {")
j = t.index("\n}\n", i) + 3
body = t[i:j].replace("printf(", "fprintf(NET_OUT ? NET_OUT : stdout, ")
t = t[:i] + "static FILE* NET_OUT = nullptr;   // web build: where print_net() writes (stdout if null)\n" + body + t[j:]
old = "if (o) { FILE* keep = stdout; stdout = o; print_net(BEST3[i].net, false); fflush(o); stdout = keep; fclose(o); }"
assert t.count(old) == 1, "the --save code of src/dbs.cpp changed: update web/patch_dbs.py"
t = t.replace(old, "if (o) { NET_OUT = o; print_net(BEST3[i].net, false); fflush(o); NET_OUT = nullptr; fclose(o); }")
assert t.count("std::clock()") == 2, "the uses of std::clock() in src/dbs.cpp changed: update web/patch_dbs.py"
t = t.replace("std::clock()", "dbs_clock()")
i = t.index("static double cpu_seconds()")
t = (t[:i] + "static clock_t dbs_clock() {   // web build: monotonic clock in units of CLOCKS_PER_SEC\n"
     "  struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);\n"
     "  return (clock_t)((double)ts.tv_sec * CLOCKS_PER_SEC + (double)ts.tv_nsec * CLOCKS_PER_SEC / 1e9);\n}\n" + t[i:])
open(sys.argv[2], "w").write(t)
import difflib
n = sum(1 for x in difflib.unified_diff(src.splitlines(), t.splitlines(), lineterm="", n=0) if x[:1] == "+" and x[:3] != "+++")
print(f"{sys.argv[2]}: {n} lines added or changed against {sys.argv[1]}")
