"""Data of each stage from the server logs (data only): best XOR2 / UMC / STM per matrix and stage, and
umc_x, the best UMC area among the recorded circuits with the fewest XOR2 gates.

usage: python3 ablation.py RUN_PARENT BENCH_DIR OUT_DIR [DEPTH4_DIR]
       (RUN_PARENT contains run/; writes OUT_DIR/ablation.json). DEPTH4_DIR, e.g. runs/depth4 of
       results/server_runs.tar.xz, adds the runs NAME.bpbs.seedN, NAME.bbeam.seedN and NAME.lns.seedN that
       paper/extra_runs.sh depth4 made for stages 2, 3 and 4 of [SS16] involutory 32x32 at depth 4."""
import os, re, glob, json, csv, sys
if len(sys.argv) not in (4, 5):
    sys.exit(__doc__)
LOGS, BENCH, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
DEPTH4 = sys.argv[4] if len(sys.argv) == 5 else None
rows = list(csv.DictReader(open(os.path.join(BENCH, 'benchmark.tsv')), delimiter='\t'))
pat = re.compile(r"BEST-(\w+): (\d+) XOR2 gates -> (\d+) XOR2 \+ (\d+) XOR3, umc180 ([\d.]+) GE, stm130 ([\d.]+) GE")
def upd(s, b):
    """update the bests of a stage with the BEST-* records of one run"""
    gx, ux = b['xor2']['g'], b['xor2']['umc']
    if gx < s['xor2']: s['xor2'], s['umc_x'] = gx, ux
    elif gx == s['xor2']: s['umc_x'] = min(s['umc_x'], ux)
    s['stm'] = min(s['stm'], b['stm130']['stm'])
    if (b['umc180']['umc'], b['umc180']['g']) < (s['umc'], s.get('umc_g', 10**9)): s['umc'], s['umc_g'] = b['umc180']['umc'], b['umc180']['g']
def parse(log):
    b, t = {}, None
    for line in open(log):
        m = pat.match(line)
        if m: b[m.group(1)] = dict(g=int(m.group(2)), n2=int(m.group(3)), n3=int(m.group(4)), umc=float(m.group(5)), stm=float(m.group(6)))
        m = re.match(r"# ([\d.]+) s", line)
        if m: t = float(m.group(1))
    return b, t
res = {}
for r in rows:
    name = r['file'][:-4]
    d = dict(label=r['label'], D=int(r['min_depth']))
    steps = {}
    cfg = {}
    tsum = 0.0; tmax = 0.0; nrun = 0
    for f in glob.glob(os.path.join(LOGS, 'run', 'forward', name + '__*.log')):
        b, t = parse(f)
        if not b: continue
        c = os.path.basename(f)[len(name) + 2:-4]
        conf = re.sub(r'_seed\d+$', '', c)
        sd = int(re.search(r'_seed(\d+)$', c).group(1))
        cfg.setdefault(conf, {})[sd] = b['xor2']['g']
        if t: tsum += t; tmax = max(tmax, t); nrun += 1
        upd(steps.setdefault('1', dict(xor2=10**9, umc=1e9, stm=1e9, umc_x=1e9)), b)
    for st, sub in (('2', 'bpbs'), ('3', 'bbeam'), ('4', 'lns')):
        per_seed = {}
        files = glob.glob(os.path.join(LOGS, 'run', sub, name + '.seed*.log'))
        if DEPTH4:
            files += glob.glob(os.path.join(DEPTH4, name + '.' + sub + '.seed*.log'))
        for f in files:
            b, _ = parse(f)
            if not b: continue
            sd = int(re.search(r'seed(\d+)', f).group(1))
            per_seed[sd] = b['xor2']['g']
            upd(steps.setdefault(st, dict(xor2=10**9, umc=1e9, stm=1e9, umc_x=1e9)), b)
        if per_seed: steps[st]['seeds'] = per_seed
    d['steps'] = steps
    d['forward_cfg'] = {k: dict(best=min(v.values()), n=len(v), mean=round(sum(v.values()) / len(v), 2),
                               seeds={str(sd): g for sd, g in sorted(v.items())}) for k, v in sorted(cfg.items())}
    d['forward_time'] = dict(total_s=round(tsum, 1), max_s=round(tmax, 1), runs=nrun)
    cum, best = {}, dict(xor2=10**9, umc=1e9, stm=1e9, umc_x=1e9)
    for st in ('1', '2', '3', '4'):
        if st in steps:
            if steps[st]['xor2'] < best['xor2']: best['xor2'], best['umc_x'] = steps[st]['xor2'], steps[st]['umc_x']
            elif steps[st]['xor2'] == best['xor2']: best['umc_x'] = min(best['umc_x'], steps[st]['umc_x'])
            best['stm'] = min(best['stm'], steps[st]['stm'])
            if (steps[st]['umc'], steps[st]['umc_g']) < (best['umc'], best.get('umc_g', 10**9)):
                best['umc'], best['umc_g'] = steps[st]['umc'], steps[st]['umc_g']
        cum[st] = dict(best)
    d['cumulative'] = cum
    res[name] = d
json.dump(res, open(os.path.join(OUT, 'ablation.json'), 'w'), indent=1)
hdr = f"{'matrix':34s} {'D':>1s} | {'s1':>4s} {'s2':>4s} {'s3':>4s} {'s4':>4s} | {'c1':>4s} {'c12':>4s} {'c123':>4s} {'all':>4s} | umc s1 s2 s3 s4"
print(hdr)
for name, d in res.items():
    s, c = d['steps'], d['cumulative']
    g = lambda st, k='xor2': (s[st][k] if st in s else None)
    print(f"{d['label'][:34]:34s} {d['D']} | " + ' '.join(f"{g(st) if g(st) is not None else '-':>4}" for st in '1234') + ' | '
          + ' '.join(f"{c[st]['xor2']:>4}" for st in '1234') + ' | ' + ' '.join(f"{g(st,'umc') if g(st,'umc') is not None else '-'}" for st in '1234'))
