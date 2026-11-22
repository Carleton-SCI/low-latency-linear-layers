# Raw data of Table 6 (beam search against random restarts at equal CPU time)

These files come from the server run of Oct 28-29, 2026 (`paper/ties_vs_restarts.py run`, seeds 1-4,
W = 2000, pairs all0, all1, sfxmd0).
- `results.jsonl`: one record per task: matrix, pair, seed; the beam search (its CPU and wall time,
  the histogram of its completed circuits, and its best circuits by XOR2 count, UMC and STM area); and
  the random restarts (CPU and wall time, number of runs and chunks, histogram, best circuits).
- `config.json`: the settings of the run. The server paths of `dbs` and `bench` are replaced by the
  repository paths.

Two copies of the run wrote into the same folder, so every task appears twice (552 records, 276 tasks).
The two records of a task have the same beam search, which is deterministic; the restarts differ only
in how many runs fit in the CPU time.

To reproduce the summary that Table 6 and Sec. 5.4 use:

```
mkdir tvr && cp results/restarts/*.json* tvr/
python3 paper/ties_vs_restarts.py summary tvr
```

The summary keeps the first record of each task and prints a warning about the duplicates. Its 69 rows
(23 matrices × 3 pairs) equal those of the paper's summary. Then
`python3 paper/tables/restarts_table.py SUMMARY.md restarts.tex` writes the table.
