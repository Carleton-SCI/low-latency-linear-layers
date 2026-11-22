"""Independent check of an XOR circuit for AES MixColumns, written without reusing the matrix code:
- computes MixColumns in GF(2^8) (xtime), not from the matrix file;
- checks the 32 unit inputs (sufficient for a linear map) plus 5000 random inputs;
- recomputes the depths from the netlist (an XOR3 gate is one level) and ignores any annotation;
- reports the gate counts, the depth histogram, dead gates, the maximum fan-out and duplicate values.

Reads the native format of dbs ('gK = A + B [d]', 'yJ = S'; x<c> = matrix column c, 0-based) and
the NIST .circ format of results/ ('XOR out a b', 'XOR3 out a b c', 'BUF out s'; x1 = column 0).

Bit order: the matrix of bench/AES.txt (from the HILL repository) numbers the bits of each byte from
the least significant one: column 8i+k is bit k of input byte i, row 8i+k is bit k of output byte
i. --order msb checks the other convention (bit 7-k); the default, auto, accepts either and says
which one holds.

usage: python3 verify_gf_aes.py CIRCUIT [D] [--order auto|lsb|msb]      (D defaults to 3)"""
import random
import re
import sys
from collections import Counter


def xt(b):
    return ((b << 1) ^ (0x1B if b & 0x80 else 0)) & 0xFF


def mixcol(a):
    a0, a1, a2, a3 = a
    return [xt(a0) ^ xt(a1) ^ a1 ^ a2 ^ a3,
            a0 ^ xt(a1) ^ xt(a2) ^ a2 ^ a3,
            a0 ^ a1 ^ xt(a2) ^ xt(a3) ^ a3,
            xt(a0) ^ a0 ^ a1 ^ a2 ^ xt(a3)]


def load(path):
    """returns (gates [(name, [operands])], outputs {j: signal}) with inputs named x0..x31"""
    text = open(path).read()
    gates, outs = [], {}
    if "begin SLP" in text:                                          # NIST .circ
        ren = {}

        def nm(s):
            m = re.fullmatch(r"x(\d+)", s)
            if m:
                return f"x{int(m.group(1)) - 1}"
            return ren.get(s, s)
        for l in text.splitlines():
            t = l.split()
            if not t or t[0] not in ("XOR", "XOR3", "BUF"):
                continue
            out, ops = t[1], [nm(s) for s in t[2:]]
            m = re.fullmatch(r"y(\d+)", out)
            if t[0] == "BUF":
                ren[out] = ops[0]
                if m:
                    outs[int(m.group(1)) - 1] = ops[0]
                continue
            g = f"n{len(gates)}"
            gates.append((g, ops))
            ren[out] = g
            if m:
                outs[int(m.group(1)) - 1] = g
    else:                                                            # native dbs format
        for l in text.splitlines():
            m = re.match(r"\s*(g\d+)\s*=\s*(\w+)\s*\+\s*(\w+)", l)
            if m:
                gates.append((m.group(1), [m.group(2), m.group(3)]))
                continue
            m = re.match(r"\s*y(\d+)\s*=\s*(\w+)\s*$", l)
            if m:
                outs[int(m.group(1))] = m.group(2)
    return gates, outs


def evaluate(gates, outs, bits):
    v = dict(bits)
    for g, ops in gates:
        x = 0
        for o in ops:
            x ^= v[o]
        v[g] = x
    return [v[outs[j]] for j in range(32)]


def check(gates, outs, pos):
    def to_bits(col):
        return {f"x{c}": (col[c // 8] >> pos(c % 8)) & 1 for c in range(32)}

    def from_bits(y):
        return [sum(y[8 * i + k] << pos(k) for k in range(8)) for i in range(4)]
    units = [[(1 << pos(c % 8)) if c // 8 == i else 0 for i in range(4)] for c in range(32)]
    ok_unit = all(from_bits(evaluate(gates, outs, to_bits(u))) == mixcol(u) for u in units)
    rng = random.Random(2026)
    ok_rand = True
    for _ in range(5000):
        a = [rng.randrange(256) for _ in range(4)]
        ok_rand &= from_bits(evaluate(gates, outs, to_bits(a))) == mixcol(a)
    return ok_unit, ok_rand


def main():
    args = [x for x in sys.argv[1:] if not x.startswith("--order")]
    order = "auto"
    for i, x in enumerate(sys.argv):
        if x.startswith("--order="):
            order = x.split("=", 1)[1]
        elif x == "--order" and i + 1 < len(sys.argv):
            order = sys.argv[i + 1]
            args.remove(order)
    if len(args) not in (1, 2) or order not in ("auto", "lsb", "msb"):
        sys.exit(__doc__)
    path, D = args[0], int(args[1]) if len(args) > 1 else 3
    gates, outs = load(path)
    names = {g for g, _ in gates}
    assert len(names) == len(gates), "gate redefined"
    assert sorted(outs) == list(range(32)), "missing outputs"
    dep = {f"x{c}": 0 for c in range(32)}
    for g, ops in gates:
        assert all(o in dep for o in ops), f"{g} uses an undefined signal"
        dep[g] = max(dep[o] for o in ops) + 1
    fan = Counter(o for _, ops in gates for o in ops)
    used = set(outs.values()) | set(fan)
    dead = [g for g in names if g not in used]
    val = {f"x{c}": 1 << c for c in range(32)}
    for g, ops in gates:
        v = 0
        for o in ops:
            v ^= val[o]
        val[g] = v
    dups = len(gates) - len({val[g] for g in names})
    orders = {"lsb": ("LSB first", lambda k: k), "msb": ("MSB first", lambda k: 7 - k)}
    tried = ["lsb", "msb"] if order == "auto" else [order]
    found = None
    for o in tried:
        ok_unit, ok_rand = check(gates, outs, orders[o][1])
        if ok_unit and ok_rand:
            found = orders[o][0]
            break
    n3 = sum(len(ops) == 3 for _, ops in gates)
    odep = max(dep[outs[j]] for j in range(32))
    print(f"{path}: {len(gates) - n3} XOR2 + {n3} XOR3; MixColumns over GF(2^8) on the 32 unit inputs and 5000 random "
          f"inputs: {'OK, bit order ' + found if found else 'FAILED (' + ', '.join(orders[o][0] for o in tried) + ')'}; "
          f"max output depth {odep} (limit {D}); gates per level {dict(sorted(Counter(dep[g] for g in names).items()))}; "
          f"dead {len(dead)}; duplicate values {dups}; max fan-out {max(fan.values())}")
    sys.exit(0 if found and odep <= D and not dead else 1)


if __name__ == "__main__":
    main()
