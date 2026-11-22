# Check: greedy "add the two signals of least depth" gives depth ceil(log2(sum 2^d_i)),
# and this equals the optimum over all binary trees (exhaustive for small k).
import itertools, math, random, functools
def greedy(ds):
    ds = sorted(ds)
    while len(ds) > 1:
        a, b = ds[0], ds[1]; ds = sorted(ds[2:] + [max(a, b) + 1])
    return ds[0]
@functools.lru_cache(None)
def opt(ds):            # ds: sorted tuple; exhaustive over all pairings
    if len(ds) == 1: return ds[0]
    best = 99
    for i in range(len(ds)):
        for j in range(i + 1, len(ds)):
            rest = list(ds[:i] + ds[i+1:j] + ds[j+1:]) + [max(ds[i], ds[j]) + 1]
            best = min(best, opt(tuple(sorted(rest))))
    return best
bad = 0
for k in range(1, 8):
    for ds in itertools.combinations_with_replacement(range(0, 4), k):
        kraft = math.ceil(math.log2(sum(2 ** d for d in ds)))
        g, o = greedy(list(ds)), opt(tuple(ds))
        if not (g == o == max(kraft, max(ds))): bad += 1; print(ds, g, o, kraft)
print("mismatches:", bad)
