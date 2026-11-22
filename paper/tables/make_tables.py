"""Generate Tables 3, 4 and 5 of the paper (XOR2 counts and areas, each stage, forward beam search by
configuration) and the facts used in the prose, from:
  ablation.json   the four stages, per matrix and stage (ablation.py, with the depth-4 runs of [SS16]
                  involutory 32x32 as its DEPTH4_DIR)
  benchmark.tsv   the benchmark (bench/)
  bests.json      the best values of each group of runs and the best circuits over all runs (final_bests.py):
                  Table 3, our columns; Table 5, the sfx/2000/1 column
  sfx_merge.json  the published circuits of SFX23 merged in the same way (paper/compare/sfx_merge.py)
Every number in these tables comes from these files.

usage: python3 make_tables.py ablation.json benchmark.tsv bests.json sfx_merge.json OUT_DIR
       (writes OUT_DIR/{results,ablation,forward}.tex and OUT_DIR/facts.json)"""
import json, csv, statistics, sys
from collections import Counter
if len(sys.argv) != 6:
    sys.exit(__doc__)
d = json.load(open(sys.argv[1]))
rows = list(csv.DictReader(open(sys.argv[2]), delimiter='\t'))
bests = json.load(open(sys.argv[3]))
sfxm = json.load(open(sys.argv[4]))
OUT = sys.argv[5]
disp = {
 'AES': (r'AES MixColumns', 32), 'SmallAES': (r'SmallScale AES', 16), 'JOLTIK': (r'Joltik', 16),
 'FSE_SKOP15_16_Hadamard': (r'[SKOP15] Hadamard', 16), 'FSE_LiuSim16_16_Circulant': (r'[LS16] circulant', 16),
 'FSE_LiWang16_16_Circulant': (r'[LW16] circulant', 16), 'ToSC_SarSye16_16_Toeplitz': (r'[SS16] Toeplitz', 16),
 'ePrint_JeaPeySim_16': (r'[JPST17]', 16), 'FSE_SKOP15_16_Involutory': (r'[SKOP15] involutory', 16),
 'FSE_LiWang16_16_Involutory': (r'[LW16] involutory', 16), 'ToSC_SarSye16_16_Involutory': (r'[SS16] involutory', 16),
 'ePrint_JeaPeySim_16_Involutory': (r'[JPST17] involutory', 16), 'FSE_SKOP15_32_Hadamard': (r'[SKOP15] Hadamard', 32),
 'FSE_LiuSim16_32_Circulant_SFX': (r'[LS16] circulant', 32), 'FSE_LiWang16_32': (r'[LW16]', 32),
 'C_BeiKraLea16_32_Circulant': (r'[BKL16] circulant', 32), 'ToSC_SarSye16_32_Toeplitz': (r'[SS16] Toeplitz', 32),
 'ePrint_JeaPeySim_32_Subfield': (r'[JPST17] subfield', 32), 'FSE_SKOP15_32_Involutory': (r'[SKOP15] involutory', 32),
 'FSE_LiWang16_32_Involutory': (r'[LW16] involutory', 32), 'ToSC_SarSye16_32_Involutory': (r'[SS16] involutory', 32),
 'KLSW': (r'[KLSW17]', 32), '168_130': (r'[LSL+19] HW168 \#130', 32)}
tsv = {r['file'][:-4]: r for r in rows}
bench_order = [r['file'][:-4] for r in rows if r['file'][:-4] in disp]
# table order: the 16x16 matrices, then the 32x32 matrices with AES MixColumns first
order = [k for k in bench_order if disp[k][1] == 16] + ['AES'] + [k for k in bench_order if disp[k][1] == 32 and k != 'AES']
# distinct problems (Sec. 4): matrices equal up to row and column permutations
problems = [['AES'], ['SmallAES'], ['JOLTIK', 'FSE_SKOP15_16_Involutory', 'FSE_LiWang16_16_Involutory'],
            ['FSE_SKOP15_16_Hadamard'], ['FSE_LiuSim16_16_Circulant', 'FSE_LiWang16_16_Circulant'],
            ['ToSC_SarSye16_16_Toeplitz'], ['ePrint_JeaPeySim_16'], ['ToSC_SarSye16_16_Involutory'],
            ['ePrint_JeaPeySim_16_Involutory'], ['FSE_SKOP15_32_Hadamard'], ['FSE_LiuSim16_32_Circulant_SFX'],
            ['FSE_LiWang16_32'], ['C_BeiKraLea16_32_Circulant'], ['ToSC_SarSye16_32_Toeplitz'],
            ['ePrint_JeaPeySim_32_Subfield'], ['FSE_SKOP15_32_Involutory'], ['FSE_LiWang16_32_Involutory'],
            ['ToSC_SarSye16_32_Involutory'], ['KLSW'], ['168_130']]
assert sorted(sum(problems, [])) == sorted(order)
cls = {k: pr for pr in problems for k in pr}
final = {k: d[k]['cumulative']['4'] for k in order}                  # the four stages (Table 4)
RANK = {'main': 0, 's4sfx': 1, 'area': 2, 'analysis': 3}
MARK = {'s4sfx': r'\ddagger', 'area': r'\S'}
def ours(k, m):
    """our value of the problem of matrix k for measure m ('xor2', 'umc180', 'stm130') over all runs, with the
    XOR2 count of its circuit, the matrix that gives it, and the groups of runs that reach it"""
    v = min(bests[j]['best'][m]['value'] for j in cls[k])
    js = [j for j in cls[k] if abs(bests[j]['best'][m]['value'] - v) < 1e-9]
    j = min(js, key=lambda j: (RANK[bests[j]['best'][m]['group']], bests[j]['best'][m]['gates'], j != k))
    groups = set().union(*[bests[i]['best'][m]['groups_reaching'] for i in js])
    return v, bests[j]['best'][m]['gates'], j, groups
own_x = {k: final[k]['xor2'] for k in order}                        # the four stages, each matrix alone
main_x = {k: min(own_x[j] for j in cls[k]) for k in order}            # the four stages, with transfers
ours_x = {k: ours(k, 'xor2')[0] for k in order}                       # all runs, with transfers
# XOR2 records before this work, with transfers within a problem and from outside the benchmark
ref_x = {k: min(int(tsv[j]['sfx23_xor2']) for j in cls[k]) for k in order}
ref_x['AES'] = 97                       # NIST AES circuit collection (2026)
ref_x['KLSW'] = 84                      # SFX23's circuit for LSL+19 HW160 #93, a row/column permutation of [KLSW17]
# UMC areas before this work: the HILL repository and the SFX23 circuits after a maximum merge (sfx_merge.json),
# with transfers within a problem; [KLSW17] via #93 (HILL 230.72); [SS16] involutory 32x32: maximum merge of
# the XOR2 circuit stored in the HILL repository (263.41); [LS16] 32x32: not in the HILL repository
hill_a = {k: (float(tsv[k]['hill_umc180']) if tsv[k]['hill_umc180'] else None) for k in order}
hill_g = {k: (int(tsv[k]['hill_input_xor2']) if tsv[k]['hill_input_xor2'] else None) for k in order}
hill_g['SmallAES'] = 47                 # published HILL circuit: 27 XOR2 + 10 XOR3 (its stored circuit has depth 5)
# per problem, as the SFX23 columns: the smallest HILL area over the matrices of the problem, and for [KLSW17]
# also HILL's circuit for #93 (230.72 GE, 52 XOR2 + 16 XOR3, i.e., 84 XOR2 gates); dagger when it comes from
# another matrix
def hill_problem(k):
    c = [(hill_a[j], hill_g[j], j != k) for j in cls[k] if hill_a[j] is not None]
    if k == 'KLSW':
        c.append((230.72, 84, True))
    return min(c) if c else (None, None, False)
hill_p = {k: hill_problem(k) for k in order}
sfx_a = {k: (sfxm[k]['umc180']['umc180'] if sfxm.get(k) else None) for k in order}
sfx_ag = {k: (sfxm[k]['umc180']['xor2'] if sfxm.get(k) else None) for k in order}
sfx_s = {k: (sfxm[k]['stm130']['stm130'] if sfxm.get(k) else None) for k in order}
sfx_sg = {k: (sfxm[k]['stm130']['xor2'] if sfxm.get(k) else None) for k in order}
def mn(*v):
    v = [x for x in v if x is not None]
    return min(v) if v else None
ref_a = {k: mn(*[hill_a[j] for j in cls[k]], *[sfx_a[j] for j in cls[k]]) for k in order}
ref_a['KLSW'] = mn(ref_a['KLSW'], 230.72)
ref_a['ToSC_SarSye16_32_Involutory'] = mn(ref_a['ToSC_SarSye16_32_Involutory'], 263.41)
ref_s = {k: mn(*[sfx_s[j] for j in cls[k]]) for k in order}
def fx(v, ref): return (r'\textbf{%d}' % v) if v < ref else str(v)
def fa(v, ref): return (r'\textbf{%.2f}' % v) if (ref is not None and v < ref - 1e-9) else '%.2f' % v
def rowlab(k): return '{' + disp[k][0] + '}'
def blocks(lines, emit, ncol):
    prev = None
    for k in order:
        n = disp[k][1]
        if n != prev:
            if prev is not None: lines.append(r'\midrule')
            lines.append(r'\multicolumn{%d}{@{}l}{\emph{$%d\times %d$ matrices}} \\' % (ncol, n, n))
        prev = n
        lines.append(emit(k))

HEAD = '% Written by paper/tables/make_tables.py of the repository.'   # first line of each table, as in the paper

# ---------------------------------------------------------------- Table 3: XOR2 counts and areas
notes = {'SmallAES': 'a', 'JOLTIK': 'b', 'AES': 'c', 'FSE_LiuSim16_32_Circulant_SFX': 'd',
         'FSE_SKOP15_32_Involutory': 'e', 'KLSW': 'f'}
DAG = r'$^{\dagger}$'
DDAG = r'$^{\ddagger}$'
def cell(a, g, other=False): return '--' if a is None else '%.2f\\,(%d)' % (a, g) + (DAG if other else '')
def sfx_other(k, m): return bool(sfxm.get(k)) and not sfxm[k][m].get('same_matrix', True)
def marks(k, m):
    v, g, j, groups = ours(k, m)
    sup = '' if 'main' in groups else ''.join(MARK[x] for x in ('s4sfx', 'area') if x in groups)
    assert 'main' in groups or sup, (k, m, groups)
    if bests[k]['best'][m]['value'] > v + 1e-9:
        sup += r'\dagger'
    return (r'$^{%s}$' % sup) if sup else ''
def ours_x_cell(k):        # dagger: from another matrix of the problem; double dagger / section: from the
    return fx(ours_x[k], ref_x[k]) + marks(k, 'xor2')   # extra stage-4 runs only
def ours_umc(k):
    a, g, _, _ = ours(k, 'umc180')
    return fa(a, ref_a[k]) + '\\,(%d)' % g + marks(k, 'umc180')
def ours_stm(k):
    a, g, _, _ = ours(k, 'stm130')
    return fa(a, ref_s[k]) + '\\,(%d)' % g + marks(k, 'stm130')
t1 = [r"""\begin{table}[tp]
\centering
\caption{XOR2 count and area (GE) at the minimum depth $D$. SFX23: Table~5 of~\cite{SFX23}, and their published circuits after our merge; HILL: the circuits of its repository~\cite{HILLcode}, which has no STM areas. The areas and our counts are the best for each problem over all our runs; in parentheses, the XOR2 count before merging. $\dagger$: from another matrix of the problem (Sec.~\ref{sec:equiv}); $\ddagger$: from stage~4 started from the circuits of SFX23; $\S$: from stage~4 with the UMC area as the objective (Sec.~\ref{sec:setup}). Bold: below the previous record (Sec.~\ref{sec:results}), below SFX23 and HILL (UMC), or below SFX23 (STM).}
\label{tab:results}
\scriptsize
\setlength{\tabcolsep}{2.2pt}
\begin{threeparttable}
\begin{tabular}{@{}l c r r r r r r r@{}}
\toprule
& & \multicolumn{2}{c}{XOR2 count} & \multicolumn{3}{c}{UMC 180\,nm area} & \multicolumn{2}{c}{STM 130\,nm area} \\
\cmidrule(lr){3-4}\cmidrule(lr){5-7}\cmidrule(l){8-9}
Matrix & $D$ & SFX23 & ours & SFX23 & HILL & ours & SFX23 & ours \\
\midrule"""]
def nt(k): return (r'\tnote{%s}' % notes[k]) if k in notes else ''
blocks(t1, lambda k: f"{rowlab(k)}{nt(k)} & {d[k]['D']} & {tsv[k]['sfx23_xor2']} & "
                     f"{ours_x_cell(k)} & {cell(sfx_a[k], sfx_ag[k], sfx_other(k, 'umc180'))} & {cell(*hill_p[k])} & "
                     f"{ours_umc(k)} & {cell(sfx_s[k], sfx_sg[k], sfx_other(k, 'stm130'))} & {ours_stm(k)} \\\\", 9)
t1.append(r"""\bottomrule
\end{tabular}
\begin{tablenotes}\footnotesize
\item[a] The HILL circuit has [27 XOR2 + 10 XOR3] gates; in parentheses, $27+2\times 10$.
\item[b] The GitHub repository~\cite{Joe26} lists a 46-gate circuit of depth~3 for Joltik, dated October~30, 2026, which we verified.
\item[c] The NIST collection of AES circuits~\cite{NISTcircuits} lists a 97-gate circuit of depth~3, added in 2026.
\item[d] The circuits of the HILL repository compute a variant of this matrix (Sec.~\ref{sec:ls16}); we use the matrix of SFX23.
\item[e] Table~3 of HILL lists 240.09~GE; we follow its repository (Sec.~\ref{sec:values}).
\item[f] SFX23 published 84 gates for [LSL+19] HW160 \#93, a row and column permutation of [KLSW17]. Daggers: the circuits of SFX23 and HILL for \#93; HILL lists 231.40~GE for [KLSW17].
\end{tablenotes}
\end{threeparttable}
\end{table}""")
open(f'{OUT}/results.tex', 'w').write('\n'.join([HEAD] + t1) + '\n')

# ---------------------------------------------------------------- Table 4: each stage
ab = [r'''\begin{table}[tp]
\centering
\caption{Best XOR2 count and UMC area (GE) of each stage (Sec.~\ref{sec:setup}); stage~4 starts from the best XOR2 circuit of stages 1--3, which its values include. Underlined: the best of the four stages; bold: as in Table~\ref{tab:results}.}
\label{tab:ablation}
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{@{}l r r r r r r r r@{}}
\toprule
& \multicolumn{4}{c}{XOR2 count, stage} & \multicolumn{4}{c}{UMC area (GE), stage} \\
\cmidrule(lr){2-5}\cmidrule(l){6-9}
Matrix & 1 & 2 & 3 & 4 & 1 & 2 & 3 & 4 \\
\midrule''']
below = {st: [] for st in '1234'}
first_x, first_a = {}, {}
def stage_row(k):
    s = d[k]['steps']
    bx = min(s[st]['xor2'] for st in s); ba = min(s[st]['umc'] for st in s)
    c = d[k]['cumulative']
    first_x[k] = next(st for st in '1234' if c[st]['xor2'] == final[k]['xor2'])
    first_a[k] = next(st for st in '1234' if abs(c[st]['umc'] - final[k]['umc']) < 1e-9)
    for st in s:
        if s[st]['xor2'] < ref_x[k]: below[st].append(k)
    def cl(st, key):
        if st not in s: return '--'
        v = s[st][key]
        txt = ('%d' % v) if key == 'xor2' else ('%.2f' % v)
        if (key == 'xor2' and v < ref_x[k]) or (key == 'umc' and ref_a[k] is not None and v < ref_a[k] - 1e-9):
            txt = r'\textbf{%s}' % txt
        if abs(v - (bx if key == 'xor2' else ba)) < 1e-9:
            txt = r'\underline{%s}' % txt
        return txt
    return f"{rowlab(k)} & " + ' & '.join(cl(st, 'xor2') for st in '1234') + ' & ' + ' & '.join(cl(st, 'umc') for st in '1234') + r' \\'
blocks(ab, stage_row, 9)
ab.append(r'''\bottomrule
\end{tabular}
\end{table}''')
open(f'{OUT}/ablation.tex', 'w').write('\n'.join([HEAD] + ab) + '\n')

# ---------------------------------------------------------------- Table 5: stage-1 configurations
cfgs = [('all_w2000_s1', 'all/2000/1', 8), ('bp_w2000_s1', 'bp/2000/1', 8), ('sfx_w2000_s0', 'sfx/2000/0', 8),
        ('sfx1', 'sfx/2000/1', 8), ('all_w8000_s1', 'all/8000/1', 4), ('all_w2000_s2', 'all/2000/2', 4),
        ('all_w2000_s1', 'all/2000/1', 32)]
NCMP = 4                                    # the configurations with 8 seeds, compared in bold
fw = [r'''\begin{table}[tp]
\centering
\caption{Forward beam search by configuration rule/$W$/$s$ (tie rule, width, slack; \emph{sfx} without the MD filter), with the number of seeds below: best XOR2 count (mean). The first column uses seeds 1--8 of the \emph{all}/2000/1 runs of stage~1, the last one all 32; \emph{sfx}/2000/1 is not part of stage~1 (Sec.~\ref{sec:setup}). Bold: the best of the first four columns.}
\label{tab:forward}
\scriptsize
\setlength{\tabcolsep}{2.5pt}
\resizebox{\textwidth}{!}{%
\begin{tabular}{@{}l ''' + ' '.join('c' for _ in cfgs) + r'''@{}}
\toprule
Matrix & ''' + ' & '.join(c[1] for c in cfgs) + r''' \\
 & ''' + ' & '.join(f'({c[2]})' for c in cfgs) + r''' \\
\midrule''']
def stats(k, cfg, nseeds):
    seeds = bests[k]['sfx1'] if cfg == 'sfx1' else d[k]['forward_cfg'][cfg]['seeds']
    vals = [seeds[str(i)] for i in range(1, nseeds + 1)] if nseeds < len(seeds) else list(seeds.values())
    assert len(vals) == nseeds, (k, cfg, nseeds, len(vals))
    return min(vals), statistics.mean(vals)
win8 = Counter(); lowmean8 = Counter()
def fw_row(k):
    st = [stats(k, c, ns) for c, _, ns in cfgs]
    best8 = min(b for b, _ in st[:NCMP]); mean8 = min(m for _, m in st[:NCMP])
    for i, (b, m) in enumerate(st[:NCMP]):
        if b == best8: win8[cfgs[i][1]] += 1
        if abs(m - mean8) < 1e-9: lowmean8[cfgs[i][1]] += 1
    return f"{rowlab(k)} & " + ' & '.join((r'\textbf{%d}' % b if (i < NCMP and b == best8) else '%d' % b) + ' (%.2f)' % m for i, (b, m) in enumerate(st)) + r' \\'
blocks(fw, fw_row, len(cfgs) + 1)
fw.append(r'''\bottomrule
\end{tabular}}
\end{table}''')
open(f'{OUT}/forward.tex', 'w').write('\n'.join([HEAD] + fw) + '\n')

# ---------------------------------------------------------------- facts used in the prose
improved = [k for k in order if ours_x[k] < ref_x[k]]                 # all runs, with transfers
improved_own = [k for k in order if own_x[k] < ref_x[k]]              # the four stages, each matrix alone
pc, pa, eq_x, eq_a, pc_main, gain = {}, {}, 0, 0, {}, {}
for pr in problems:
    rx = min(ref_x[k] for k in pr); ox = ours_x[pr[0]]; om = min(own_x[k] for k in pr)
    if ox < rx: pc[pr[0]] = round(100 * (rx - ox) / rx, 1); gain[pr[0]] = rx - ox
    elif ox == rx: eq_x += 1
    if om < rx: pc_main[pr[0]] = round(100 * (rx - om) / rx, 1)
    ra = min(ref_a[k] for k in pr); oa = ours(pr[0], 'umc180')[0]
    if oa < ra - 1e-9: pa[pr[0]] = (round(100 * (ra - oa) / ra, 2), round(ra - oa, 2))
    elif abs(oa - ra) < 1e-9: eq_a += 1
# the XOR2 count of our smallest UMC circuit against that of the baseline circuit (the published circuit of the
# smallest UMC area for the same problem), for the problems with a smaller area
base_g = {}
for pr in problems:
    c = [(hill_a[j], hill_g[j]) for j in pr if hill_a[j] is not None] + [(sfx_a[j], sfx_ag[j]) for j in pr if sfx_a[j] is not None]
    if 'KLSW' in pr: c.append((230.72, 84))                         # HILL's circuit for #93
    if 'ToSC_SarSye16_32_Involutory' in pr: c.append((263.41, 98))  # maximum merge of HILL's stored circuit
    ra = min(a for a, _ in c); gs = {g for a, g in c if abs(a - ra) < 1e-9}
    assert abs(ra - min(ref_a[k] for k in pr)) < 1e-9 and len(gs) == 1, (pr, c)
    base_g[pr[0]] = gs.pop()
gates_vs_base = {'fewer': [], 'equal': [], 'more': []}
for pr in problems:
    if pr[0] not in pa: continue
    g, h = ours(pr[0], 'umc180')[1], base_g[pr[0]]
    gates_vs_base['fewer' if g < h else ('equal' if g == h else 'more')].append(pr[0])
# areas if only the recorded circuits with the fewest XOR2 gates were merged (all runs)
mincount_vs_base = {'below': [], 'equal': [], 'above': []}
for pr in problems:
    g = ours_x[pr[0]]
    a = min(bests[k]['umc_of_fewest'] for k in pr if bests[k]['best']['xor2']['value'] == g)
    ra = min(ref_a[k] for k in pr)
    mincount_vs_base['below' if a < ra - 1e-9 else ('equal' if abs(a - ra) < 1e-9 else 'above')].append(pr[0])
# per matrix (all runs): the smallest UMC circuit has more XOR2 gates than the fewest, and the area gap
count_vs_area = {k: round(bests[k]['umc_of_fewest'] - bests[k]['best']['umc180']['value'], 2) for k in order
                 if bests[k]['best']['umc180']['gates'] > bests[k]['best']['xor2']['value']}
# STM: the duplication of gates of fan-out two, per matrix, against the smallest STM area without it
stm_dup_gain = {k: round(bests[k]['nodup_stm'] - bests[k]['best']['stm130']['value'], 2) for k in order
                if bests[k]['best']['stm130']['value'] < bests[k]['nodup_stm'] - 1e-9}
# the extra stage-4 runs against the four stages, per matrix
g4 = lambda k, grp, m: bests[k]['groups'][grp][m] if grp in bests[k]['groups'] else None
s4_vs_main = {'below': [], 'equal': [], 'above': []}
area_vs_main = {}
for k in order:
    x = g4(k, 's4sfx', 'xor2')
    if x is not None:
        s4_vs_main['below' if x < own_x[k] else ('equal' if x == own_x[k] else 'above')].append(k)
    u = g4(k, 'area', 'umc180')
    if u is not None and u[0] < final[k]['umc'] - 1e-9:
        area_vs_main[k] = round(final[k]['umc'] - u[0], 2)
facts = dict(improved=improved, improved_own_runs=improved_own, first_x=first_x, first_a=first_a,
             stage_alone_below={st: [k for k in v if k in improved_own] for st, v in below.items()},
             forward_best8=dict(win8), forward_lowest_mean8=dict(lowmean8),
             stm_dup_gain=stm_dup_gain,
             sfx_merge_below_hill=[k for k in order if sfx_a[k] is not None and hill_a[k] is not None and sfx_a[k] < hill_a[k] - 1e-9],
             hill_vs_sfx_merged_by_problem={w: [pr[0] for pr in problems if hill_p[pr[0]][0] is not None and sfx_a[pr[0]] is not None
                                                and (w == 'smaller') == (hill_p[pr[0]][0] < sfx_a[pr[0]] - 1e-9)
                                                and (w == 'larger') == (hill_p[pr[0]][0] > sfx_a[pr[0]] + 1e-9)]
                                            for w in ('smaller', 'same', 'larger')},
             gates_of_umc_circuit_vs_baseline=gates_vs_base, mincount_circuits_vs_baseline=mincount_vs_base,
             count_vs_area_gap=count_vs_area, s4sfx_vs_four_stages=s4_vs_main, area_runs_below_four_stages=area_vs_main,
             improved_matrices_by_size={n: len([k for k in improved if disp[k][1] == n]) for n in (16, 32)},
             problems=dict(n=len(problems), xor2_improved=len(pc), xor2_improved_four_stages=len(pc_main), xor2_equal=eq_x,
                           xor2_gain_counts=dict(sorted(Counter(gain.values()).items())), xor2_gain=gain,
                           area_below=len(pa), area_equal=eq_a,
                           area_pct_range=[min(v[0] for v in pa.values()), max(v[0] for v in pa.values())],
                           area_ge_range=[min(v[1] for v in pa.values()), max(v[1] for v in pa.values())], area_pct=pa),
             step1_hours=round(sum(d[k]['forward_time']['total_s'] for k in d) / 3600, 1))
json.dump(facts, open(f'{OUT}/facts.json', 'w'), indent=1)
print('improved XOR2 entries (all runs, with transfers):', len(improved), facts['improved_matrices_by_size'],
      '; by the four stages, each matrix alone:', len(improved_own))
for st in '1234':
    print(f'stage {st} alone below the previous record: {len(facts["stage_alone_below"][st])} of the {len(improved_own)}')
print('final XOR2 first reached (cumulative) at stage:', Counter(first_x[k] for k in improved_own))
print('final area first reached (cumulative) at stage:', Counter(first_a.values()))
print('8 seeds: best of the four configurations:', dict(win8), '; lowest mean:', dict(lowmean8))
print('STM duplication lowers the smallest STM area for', len(stm_dup_gain), 'matrices, by',
      min(stm_dup_gain.values()), '-', max(stm_dup_gain.values()), 'GE')
print('HILL against SFX23 merged, by problem:', {w: len(v) for w, v in facts['hill_vs_sfx_merged_by_problem'].items()})
print('gates of our smallest-UMC circuit vs the baseline circuit (problems with a smaller area):', {k: len(v) for k, v in gates_vs_base.items()})
print('only the circuits with the fewest XOR2 gates, against the baseline (problems):', {k: len(v) for k, v in mincount_vs_base.items()},
      mincount_vs_base['above'], mincount_vs_base['equal'])
print('matrices whose smallest UMC circuit has more XOR2 gates than the fewest:', len(count_vs_area),
      'gap', min(count_vs_area.values()), '-', max(count_vs_area.values()), 'GE')
print('stage 4 from SFX23 against the four stages (XOR2, per matrix):', {k: len(v) for k, v in s4_vs_main.items()}, s4_vs_main['below'], s4_vs_main['above'])
print('stage 4 with the UMC objective below the four stages (UMC, per matrix):', area_vs_main)
p = facts['problems']
print(f"problems: {p['n']}; XOR2 improved {p['xor2_improved']} (gains {p['xor2_gain_counts']}), equal {p['xor2_equal']}, "
      f"by the four stages {p['xor2_improved_four_stages']}; area below {p['area_below']} ({p['area_pct_range']} %, "
      f"{p['area_ge_range']} GE), equal {p['area_equal']}")
print('step-1 hours:', facts['step1_hours'])
