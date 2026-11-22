// dbs: depth-bounded search for low-latency implementations of binary linear layers: the
// shortest linear program (XOR2 gates) under a depth bound, and the XOR2/XOR3 area after merging.
// This is the search program of the paper; slpd.py runs it on user matrices (see README.md).
//
//   * forward search with exact depth-aware distances (Boyar-Peralta type; the depth bound of
//     Li et al., ToSC 2019(1), Lemma 1), per-target depth bounds, the IBPD-MD filter of
//     Shi-Feng-Xu (ToSC 2023(4)) (--md) and the tie rules bp, sfx and all (--tie);
//   * beam search over the ties (--beam W, --slack s);
//   * backward splitting of the outputs around the forward search (Shi-Feng-Xu's BPBS framework,
//     their priority (3): Rule 1 > Rule 2 > Rule 3, Rule 5), with random trajectories (--bpbs N)
//     or a beam search over the splits (--bbeam B --bopts K);
//   * large neighborhood search around a circuit (--lns): keep the last gate of k random outputs and
//     re-solve everything else;
//   * selection by XOR2 count or by XOR2/XOR3 area after absorbing fan-out-1 gates into XOR3
//     gates (maximum matching; --obj xor|umc|stm, cell areas from HILL, ToSC 2026(1), Table 1).
//
// Distances: target t with budget B = 2^b has distance
//   delta(t) = min |S| - 1 over basis subsets S with XOR(S) = t and sum_{s in S} 2^depth(s) <= B.
// (The bound ceil(log2(sum 2^d_i)) is Lemma 1 of Li et al., ToSC 2019(1).)
//
// usage: dbs matrix.txt D [options]
//   forward search:             [--tie bp|sfx|all|sdf|fs] [--md] [--runs R] [--seed S]
//                               [--beam W [--slack s] [--rank sigma|gates] [--stats]]   (--rank gates: rank the
//                               options of the beam search by gates so far + 1 + Sigma; default: Sigma;
//                               --stats: per-iteration diagnostics of the beam search on stderr)
//   backward splits + forward:  --bpbs N [--evals R] [--p3 q] [--maxsplits K] [--topk M --beam W [--slack s]]
//                               [--btie RULE] [--bmd] (tie rule / MD filter of the beam refinement) [--verbose]
//                               --bbeam B --bopts K: beam search over the backward splits (B prefixes kept,
//                               K split options drawn per prefix); --bpbs N then repeats it N times
//   large neighborhood search:  --lns circuit.txt --iters N --fix kmin[:kmax] [--plateau q]
//                               (+ --beam W [--slack s] or --evals R for the forward search)
//   selection / output:         [--obj xor|umc|stm] [--print] [--print3] [--save PREFIX] [--trace]
//                               (--trace: CPU time of every improvement of the best XOR2 count and UMC area, stderr)
//   (--save writes the best circuit by each measure to PREFIX.{xor2,umc180,stm130}.txt)
//   [--progress]: for randomized forward runs (--runs), PROGRESS lines on stderr with the CPU time of the
//                 process after every target scanned for candidate gates and after every gate added
// matrix.txt: one row per output (0/1 separated by spaces); inputs x<c> = column c.
// Tie rules sdf and fs are the Shortest-Dist-First and Focused-Search selections of
// Reyhani-Masoleh, Taha and Ashmawy (TCHES 2018(2)) on exact distances (not used in the paper).
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <random>
#include <fstream>
#include <sstream>
#include <map>
#include <ctime>
typedef uint64_t u64;
typedef unsigned __int128 mask_t;

static int NIN = 0, DMAX = 3, BUDMAX = 8;
static inline int W(int d) { return 1 << d; }
static inline int popc(u64 x) { return __builtin_popcountll(x); }
static inline int dmin_of(u64 x) { int w = popc(x), d = 0; while ((1 << d) < w) d++; return d; }
static inline mask_t bit(int j) { return ((mask_t)1) << j; }
static inline int ctz128(mask_t m) { u64 lo = (u64)m; return lo ? __builtin_ctzll(lo) : 64 + __builtin_ctzll((u64)(m >> 64)); }

struct Gate { u64 v; int a, b, d; };
struct State {
  std::vector<u64> val; std::vector<int> dep;    // basis (inputs first)
  std::vector<Gate> gates;                       // gates in creation order (basis index = NIN + k)
  std::unordered_map<u64, int> where;            // value -> basis index of minimum depth
  std::vector<u64> T; std::vector<int> TB;       // targets and their Kraft budgets 2^b
  std::vector<int> delta; std::vector<char> found;
};

template <class F>
static void subsets(const std::vector<int>& idx, const std::vector<u64>& V, const std::vector<int>& Wt,
                    int k, int wmax, F&& f) {
  int n = (int)idx.size();
  std::vector<int> st(k); std::vector<u64> acc(k + 1, 0); std::vector<int> wac(k + 1, 0);
  int lvl = 0; if (k == 0) { f(0ULL, 0); return; }
  st[0] = 0;
  while (lvl >= 0) {
    if (st[lvl] > n - (k - lvl)) { lvl--; if (lvl >= 0) st[lvl]++; continue; }
    int e = idx[st[lvl]];
    int w = wac[lvl] + Wt[e];
    if (w > wmax) { lvl--; if (lvl >= 0) st[lvl]++; continue; }   // idx sorted by weight
    acc[lvl + 1] = acc[lvl] ^ V[e]; wac[lvl + 1] = w;
    if (lvl + 1 == k) { f(acc[k], w); st[lvl]++; continue; }
    st[lvl + 1] = st[lvl] + 1; lvl++;
  }
}

struct Cand { long S, N; u64 v; mask_t red; int d, a, b; };

// CPU time of the process in seconds (std::clock; the program is single-threaded)
static double cpu_seconds() { return (double)std::clock() / CLOCKS_PER_SEC; }

struct Solver {
  State s; std::mt19937_64 rng; int tie = 2; bool firsttie = false, md = false;
  bool progress = false;                         // --progress (randomized forward runs only)
  void progress_line(const char* what) {         // state of the run on stderr; reads the state only
    long sum = 0; int mx = 0, left = 0;
    for (size_t j = 0; j < s.T.size(); ++j) if (!s.found[j]) { sum += s.delta[j]; mx = std::max(mx, s.delta[j]); left++; }
    fprintf(stderr, "PROGRESS cpu=%.1f %s gates=%zu targets_left=%d sum_delta=%ld max_delta=%d\n",
            cpu_seconds(), what, s.gates.size(), left, sum, mx);
  }
  std::vector<int> order;
  void rebuild_order() {
    order.clear();
    for (int i = 0; i < (int)s.val.size(); ++i) if (W(s.dep[i]) < BUDMAX) order.push_back(i);
    std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return s.dep[x] < s.dep[y]; });
  }
  std::vector<int> wt() { std::vector<int> w(s.val.size()); for (size_t i = 0; i < w.size(); ++i) w[i] = W(s.dep[i]); return w; }
  int add(u64 v, int a, int b) {
    int d = std::max(s.dep[a], s.dep[b]) + 1;
    int id = (int)s.val.size();
    s.val.push_back(v); s.dep.push_back(d); s.gates.push_back({v, a, b, d});
    auto it = s.where.find(v);
    if (it == s.where.end() || s.dep[it->second] > d) s.where[v] = id;
    return id;
  }
  bool exists(u64 x, int k, int wmax, int skip) {
    if (k == 0) return x == 0;
    std::vector<int> idx; idx.reserve(order.size());
    for (int i : order) if (i != skip) idx.push_back(i);
    bool ok = false; auto w = wt();
    subsets(idx, s.val, w, k, wmax, [&](u64 acc, int) { if (acc == x) ok = true; });
    return ok;
  }
  void update_after(int e) {
    u64 ve = s.val[e]; int de = s.dep[e];
    for (size_t j = 0; j < s.T.size(); ++j) {
      if (s.found[j]) continue;
      if (s.T[j] == ve && W(de) <= s.TB[j]) { s.found[j] = 1; s.delta[j] = 0; continue; }
      int k = s.delta[j] - 1;
      if (k < 0) continue;
      if (s.TB[j] - W(de) < 0) continue;
      if (exists(s.T[j] ^ ve, k, s.TB[j] - W(de), e)) s.delta[j]--;
    }
  }
  bool cascade() {                                // add every target at distance 1 (best depth)
    bool any = false;
    while (true) {
      int bj = -1, ba = -1, bb = -1, bd = 99;
      for (size_t j = 0; j < s.T.size(); ++j) {
        if (s.found[j] || s.delta[j] != 1) continue;
        for (int x : order) {
          auto it = s.where.find(s.T[j] ^ s.val[x]);
          if (it == s.where.end()) continue;
          int y = it->second;
          if (y == x) continue;
          int d = std::max(s.dep[x], s.dep[y]) + 1;
          if (W(d) <= s.TB[j] && d < bd) { bd = d; bj = (int)j; ba = x; bb = y; }
        }
      }
      if (bj < 0) break;
      int e = add(s.T[bj], ba, bb);
      rebuild_order();
      update_after(e);
      any = true;
    }
    return any;
  }
  bool init(const std::vector<u64>& T, const std::vector<int>& TB) {
    s = State(); s.T = T; s.TB = TB;
    for (int i = 0; i < NIN; ++i) { s.val.push_back(1ULL << i); s.dep.push_back(0); s.where[1ULL << i] = i; }
    s.delta.resize(T.size()); s.found.assign(T.size(), 0);
    for (size_t j = 0; j < T.size(); ++j) {
      int w = popc(T[j]);
      if (w > TB[j]) return false;               // not reachable within its depth bound
      s.delta[j] = w - 1; if (w == 1) s.found[j] = 1;
    }
    rebuild_order(); cascade();
    return true;
  }
  bool complete() const { for (char f : s.found) if (!f) return false; return true; }
  // all useful candidates with their scores (Sum delta', Sum delta'^2); IBPD-MD filter if md
  void score(std::vector<Cand>& out) {
    out.clear();
    const std::vector<u64>& T = s.T;
    std::unordered_map<u64, int> cid; std::vector<Cand> cv;
    for (size_t x = 0; x < order.size(); ++x)
      for (size_t y = x + 1; y < order.size(); ++y) {
        int a = order[x], b = order[y];
        int d = std::max(s.dep[a], s.dep[b]) + 1;
        if (d > DMAX - 1) continue;
        u64 v = s.val[a] ^ s.val[b];
        if (s.where.count(v)) continue;
        auto it = cid.find(v);
        if (it == cid.end()) { cid[v] = (int)cv.size(); cv.push_back({0, 0, v, 0, d, a, b}); }
        else if (d < cv[it->second].d) { cv[it->second].d = d; cv[it->second].a = a; cv[it->second].b = b; }
      }
    auto w = wt();
    for (size_t j = 0; j < T.size(); ++j) {
      if (s.found[j]) continue;
      int k = s.delta[j] - 1, tb = s.TB[j];
      u64 tj = T[j]; mask_t bj = bit((int)j);
      subsets(order, s.val, w, k, tb - 2, [&](u64 acc, int ws) {
        auto it = cid.find(tj ^ acc);
        if (it != cid.end() && ws + W(cv[it->second].d) <= tb) cv[it->second].red |= bj;
      });
      if (progress) fprintf(stderr, "PROGRESS cpu=%.1f scanned gates=%zu target=%zu delta=%d\n", cpu_seconds(), s.gates.size(), j, s.delta[j]);
    }
    long base = 0, base2 = 0; int mind = 1 << 30;
    for (size_t j = 0; j < T.size(); ++j) if (!s.found[j]) { base += s.delta[j]; base2 += (long)s.delta[j] * s.delta[j]; mind = std::min(mind, s.delta[j]); }
    mask_t near = 0;
    if (md) for (size_t j = 0; j < T.size(); ++j) if (!s.found[j] && s.delta[j] == mind) near |= bit((int)j);
    for (auto& c : cv) {
      if (!c.red) continue;
      if (md && !(c.red & near)) continue;
      c.S = base; c.N = base2;
      for (mask_t m = c.red; m; m &= m - 1) { int j = ctz128(m); c.S -= 1; c.N -= 2L * s.delta[j] - 1; }
      if (tie >= 3) {                              // sdf / fs rules (TCHES 2018) on exact distances
        int z = 99, cnt = 0;
        for (size_t j = 0; j < T.size(); ++j) {
          if (s.found[j]) continue;
          int dj = s.delta[j] - (int)((c.red >> j) & 1);
          if (dj <= 1 && tie == 3) continue;
          if (dj < z) { z = dj; cnt = 1; } else if (dj == z) cnt++;
        }
        if (tie == 3) { c.S = -((10L - std::min(z, 10)) * 100 + cnt); c.N = 0; }
        else { c.S = z; c.N = 0; }
      }
      out.push_back(c);
    }
  }
  void apply(const Cand& c) {
    add(c.v, c.a, c.b);
    for (mask_t m = c.red; m; m &= m - 1) s.delta[ctz128(m)]--;
    rebuild_order(); cascade();
  }
  // one greedy run with random tie-breaking (per-target budgets, optional --md filter)
  int run_greedy(int maxg) {
    std::vector<Cand> cs;
    if (progress) progress_line("start");
    while ((int)s.gates.size() < maxg) {
      if (complete()) return (int)s.gates.size();
      score(cs);
      if (cs.empty()) return -1;
      long bestS = 1L << 60, bestN = 0; std::vector<int> ties;
      for (size_t c = 0; c < cs.size(); ++c) {
        long S = cs[c].S, N = cs[c].N;
        bool better = S < bestS || (S == bestS && (tie == 0 ? N > bestN : (tie == 1 ? N < bestN : false)));
        bool equal = S == bestS && (tie >= 2 || N == bestN);
        if (better) { bestS = S; bestN = N; ties.assign(1, (int)c); }
        else if (equal) ties.push_back((int)c);
      }
      int c = firsttie ? ties[0] : ties[std::uniform_int_distribution<int>(0, (int)ties.size() - 1)(rng)];
      apply(cs[c]);
      if (progress) progress_line("added");
    }
    return -1;
  }
  u64 sethash() const {
    u64 h = 0; for (auto& g : s.gates) { u64 x = g.v * 0x9E3779B97F4A7C15ULL; x ^= x >> 29; h += x * 0xBF58476D1CE4E5B9ULL; } return h;
  }
};

// ---------------------------------------------------------------- netlists ----------------------
struct Net { std::vector<Gate> g; std::vector<int> out; };   // node ids: 0..NIN-1 inputs, NIN+k gate k
struct Split { u64 w, p1, p2; };

// forward state + split gates -> netlist computing the original outputs `outs`
static bool build_net(const State& s, const std::vector<Split>& splits, const std::vector<u64>& outs, Net& net) {
  net = Net();
  std::unordered_map<u64, int> idx; std::vector<int> dep(NIN, 0);
  for (int i = 0; i < NIN; ++i) idx[1ULL << i] = i;
  for (auto& G : s.gates) {
    int id = NIN + (int)net.g.size(); net.g.push_back(G); dep.push_back(G.d);
    auto it = idx.find(G.v); if (it == idx.end() || dep[it->second] > G.d) idx[G.v] = id;
  }
  for (auto& sp : splits) {
    auto a = idx.find(sp.p1), b = idx.find(sp.p2);
    if (a == idx.end() || b == idx.end()) return false;
    int d = std::max(dep[a->second], dep[b->second]) + 1;
    int id = NIN + (int)net.g.size(); net.g.push_back({sp.w, a->second, b->second, d}); dep.push_back(d);
    auto it = idx.find(sp.w); if (it == idx.end() || dep[it->second] > d) idx[sp.w] = id;
  }
  for (u64 t : outs) { auto it = idx.find(t); if (it == idx.end()) return false; net.out.push_back(it->second); }
  return true;
}

static Net prune(const Net& n) {
  int ng = (int)n.g.size();
  std::vector<char> need(NIN + ng, 0);
  for (int o : n.out) need[o] = 1;
  for (int k = ng - 1; k >= 0; --k) if (need[NIN + k]) { need[n.g[k].a] = 1; need[n.g[k].b] = 1; }
  std::vector<int> nid(NIN + ng, -1); for (int i = 0; i < NIN; ++i) nid[i] = i;
  Net r;
  for (int k = 0; k < ng; ++k) if (need[NIN + k]) {
    Gate G = n.g[k]; G.a = nid[G.a]; G.b = nid[G.b];
    nid[NIN + k] = NIN + (int)r.g.size(); r.g.push_back(G);
  }
  for (int o : n.out) r.out.push_back(nid[o]);
  return r;
}

// recompute values and depths; check outputs and depth bound
static bool verify(const Net& n, const std::vector<u64>& outs, int D) {
  std::vector<u64> v(NIN + n.g.size()); std::vector<int> d(NIN + n.g.size(), 0);
  for (int i = 0; i < NIN; ++i) v[i] = 1ULL << i;
  for (size_t k = 0; k < n.g.size(); ++k) {
    const Gate& G = n.g[k];
    if (G.a >= NIN + (int)k || G.b >= NIN + (int)k) return false;
    v[NIN + k] = v[G.a] ^ v[G.b]; d[NIN + k] = std::max(d[G.a], d[G.b]) + 1;
    if (v[NIN + k] != G.v || d[NIN + k] != G.d) return false;
  }
  for (size_t j = 0; j < outs.size(); ++j) if (v[n.out[j]] != outs[j] || d[n.out[j]] > D) return false;
  return true;
}

// absorb fan-out-1, non-output gates into their consumer as XOR3 (maximum matching on the forest)
static void merge3(const Net& n, std::vector<int>& matchedTo) {    // matchedTo[k] = consumer gate or -1
  int ng = (int)n.g.size();
  std::vector<int> fan(NIN + ng, 0), cons(NIN + ng, -1); std::vector<char> isout(NIN + ng, 0);
  for (int o : n.out) isout[o] = 1;
  for (int k = 0; k < ng; ++k) for (int x : {n.g[k].a, n.g[k].b}) { fan[x]++; cons[x] = k; }
  std::vector<std::vector<int>> ch(ng);
  std::vector<char> hasParent(ng, 0);
  for (int k = 0; k < ng; ++k)
    if (!isout[NIN + k] && fan[NIN + k] == 1) { ch[cons[NIN + k]].push_back(k); hasParent[k] = 1; }
  std::vector<int> f0(ng), f1(ng), pick(ng, -1);
  for (int v = 0; v < ng; ++v) {                  // children have smaller indices
    int base = 0; for (int c : ch[v]) base += std::max(f0[c], f1[c]);
    f0[v] = base; f1[v] = -1;
    for (int c : ch[v]) { int val = base - std::max(f0[c], f1[c]) + f0[c] + 1; if (val > f1[v]) { f1[v] = val; pick[v] = c; } }
  }
  matchedTo.assign(ng, -1);
  std::vector<std::pair<int, bool>> st;          // (node, may be matched with a child)
  for (int v = ng - 1; v >= 0; --v) if (!hasParent[v]) st.push_back({v, true});
  while (!st.empty()) {
    auto [v, may] = st.back(); st.pop_back();
    bool use = may && f1[v] > f0[v];
    for (int c : ch[v]) {
      if (use && c == pick[v]) { matchedTo[c] = v; st.push_back({c, false}); }
      else st.push_back({c, true});
    }
  }
}
static int OBJ = 0;                               // 0 = XOR2 count, 1 = UMC 180 nm area, 2 = STM 130 nm area
static double area_of(int n2, int n3, int lib) { return lib == 1 ? 3.00 * n2 + 4.67 * n3 : 3.33 * n2 + 4.66 * n3; }
struct Eval { int gates = 1 << 30, n2 = 0, n3 = 0; double key = 1e18; };
static Eval evaluate_net(const Net& n) {
  Eval e; e.gates = (int)n.g.size();
  std::vector<int> m; merge3(n, m);
  int merged = 0; for (int x : m) merged += x >= 0;
  e.n3 = merged; e.n2 = e.gates - 2 * merged;
  e.key = OBJ == 0 ? e.gates + 1e-3 * area_of(e.n2, e.n3, 1) : area_of(e.n2, e.n3, OBJ) + 1e-6 * e.gates;
  return e;
}
static std::string nm(int i) { char b[16]; if (i < NIN) snprintf(b, 16, "x%d", NIN - 1 - i); else snprintf(b, 16, "g%d", i - NIN); return b; }
static void print_net(const Net& n, bool xor3) {
  if (!xor3) {
    for (size_t k = 0; k < n.g.size(); ++k)
      printf("%s = %s + %s  [%d]\n", nm(NIN + (int)k).c_str(), nm(n.g[k].a).c_str(), nm(n.g[k].b).c_str(), n.g[k].d);
  } else {
    std::vector<int> m; merge3(n, m);
    std::vector<std::vector<int>> ops(n.g.size());
    for (size_t k = 0; k < n.g.size(); ++k) ops[k] = {n.g[k].a, n.g[k].b};
    for (size_t c = 0; c < n.g.size(); ++c) if (m[c] >= 0) {
      auto& o = ops[m[c]]; o.erase(std::find(o.begin(), o.end(), NIN + (int)c));
      o.push_back(n.g[c].a); o.push_back(n.g[c].b);
    }
    for (size_t k = 0; k < n.g.size(); ++k) {
      if (m[k] >= 0) continue;
      printf("%s =", nm(NIN + (int)k).c_str());
      for (size_t t = 0; t < ops[k].size(); ++t) printf("%s %s", t ? " +" : "", nm(ops[k][t]).c_str());
      printf("\n");
    }
  }
  for (size_t j = 0; j < n.out.size(); ++j) printf("y%zu = %s\n", j, nm(n.out[j]).c_str());
}

// ---------------------------------------------------------------- forward searches --------------
struct Result { Eval e; Net net; bool ok = false; };
static Result BEST3[3];                           // best circuit seen so far by XOR2 count, UMC area, STM area
// --trace: one TRACE line on stderr (with the CPU time of the process) whenever the best circuit of the run by
// XOR2 count or by UMC area improves, for anytime curves; the search and its output are unchanged
static bool TRACE = false;
static void track3(const Eval& e, const Net& p) {
  double k[3] = {e.gates + 1e-3 * area_of(e.n2, e.n3, 1), area_of(e.n2, e.n3, 1) + 1e-6 * e.gates, area_of(e.n2, e.n3, 2) + 1e-6 * e.gates};
  for (int i = 0; i < 3; ++i)
    if (!BEST3[i].ok || k[i] < BEST3[i].e.key) {
      if (TRACE && i < 2)
        fprintf(stderr, "TRACE cpu=%.2f measure=%s xor2=%d umc180=%.2f\n", (double)std::clock() / CLOCKS_PER_SEC,
                i == 0 ? "xor2" : "umc180", e.gates, area_of(e.n2, e.n3, 1));
      BEST3[i].e = e; BEST3[i].e.key = k[i]; BEST3[i].net = p; BEST3[i].ok = true;
    }
}
static void consider(Result& best, const State& s, const std::vector<Split>& sp, const std::vector<u64>& outs, std::map<int, int>* hist) {
  Net n; if (!build_net(s, sp, outs, n)) return;
  Net p = prune(n);
  if (!verify(p, outs, DMAX)) { fprintf(stderr, "internal error: circuit does not verify\n"); exit(2); }
  Eval e = evaluate_net(p);
  if (hist) (*hist)[e.gates]++;
  track3(e, p);
  if (!best.ok || e.key < best.e.key) { best.e = e; best.net = p; best.ok = true; }
}
static bool PROGRESS = false;                     // --progress (see the header)
static void greedy_runs(const std::vector<u64>& T, const std::vector<int>& TB, int runs, u64 seed, int tie, bool md,
                        const std::vector<Split>& sp, const std::vector<u64>& outs, Result& best, std::map<int, int>* hist,
                        bool firsttie0 = true) {
  Solver S; S.tie = tie; S.md = md; S.progress = PROGRESS;
  for (int r = 0; r < runs; ++r) {
    S.rng.seed(seed * 1000003ULL + r); S.firsttie = firsttie0 && (seed == 0 && r == 0);
    if (!S.init(T, TB)) return;
    if (S.run_greedy(400) < 0) continue;
    consider(best, S.s, sp, outs, hist);
  }
}
// --rank gates: rank the options of the beam search by (gates so far + 1 + Sigma) instead of Sigma alone
static bool RANK_GATES = false;
// --stats: diagnostics of the beam search on stderr (one STAT line per iteration, one STATEND line per
// beam search); they only read the state, so the search and its output are unchanged
static bool STATS = false;

static void beam_search(const std::vector<u64>& T, const std::vector<int>& TB, int width, long slack, u64 seed, int tie, bool md,
                        const std::vector<Split>& sp, const std::vector<u64>& outs, Result& best, std::map<int, int>* hist) {
  std::mt19937_64 rng(seed);
  std::vector<Solver> B(1); B[0].tie = tie; B[0].md = md;
  if (!B[0].init(T, TB)) return;
  std::vector<Cand> cs; int st_iter = 0;
  while (true) {
    bool all = true; for (auto& b : B) if (!b.complete()) all = false;
    if (all) break;
    struct Opt { int bi; Cand c; };
    std::vector<Opt> opts; long bestS = 1L << 60;
    std::vector<Solver> nextDone;
    for (int bi = 0; bi < (int)B.size(); ++bi) {
      if (B[bi].complete()) { nextDone.push_back(B[bi]); continue; }
      B[bi].score(cs);
      if (RANK_GATES) for (auto& c : cs) c.S += (long)B[bi].s.gates.size() + 1;
      for (auto& c : cs) { opts.push_back({bi, c}); bestS = std::min(bestS, c.S); }
    }
    std::vector<Opt> keep;
    for (auto& o : opts) if (o.c.S <= bestS + slack) keep.push_back(o);
    int st_parents_keep = 0, st_dup = 0;
    if (STATS) { std::unordered_set<int> pk; for (auto& o : keep) pk.insert(o.bi); st_parents_keep = (int)pk.size(); }
    std::shuffle(keep.begin(), keep.end(), rng);
    if (tie <= 1) std::stable_sort(keep.begin(), keep.end(), [&](const Opt& x, const Opt& y) {
      if (x.c.S != y.c.S) return x.c.S < y.c.S;
      return tie == 0 ? x.c.N > y.c.N : x.c.N < y.c.N; });
    else std::stable_sort(keep.begin(), keep.end(), [&](const Opt& x, const Opt& y) { return x.c.S < y.c.S; });
    std::vector<Solver> NB = nextDone; std::unordered_map<u64, int> seen;
    std::unordered_set<int> st_parents_new;
    for (auto& o : keep) {
      if ((int)NB.size() >= width) break;
      Solver nb = B[o.bi]; nb.apply(o.c);
      u64 h = nb.sethash();
      if (seen.count(h)) { st_dup++; continue; }
      seen[h] = 1; NB.push_back(std::move(nb));
      if (STATS) st_parents_new.insert(o.bi);
    }
    if (STATS) {
      fprintf(stderr, "STAT iter=%d seqs=%zu complete=%zu options=%zu kept=%zu dmin=%ld parents_kept=%d new=%zu dup=%d parents_new=%zu\n",
              ++st_iter, B.size(), nextDone.size(), opts.size(), keep.size(), bestS, st_parents_keep,
              NB.size() - nextDone.size(), st_dup, st_parents_new.size());
    }
    if (NB.empty()) break;
    B.swap(NB);
  }
  if (STATS) {                                   // gate counts before pruning (unpruned) and after pruning
    int umin = 1 << 30, best_at_umin = 1 << 30, best_all = 1 << 30, umin_best = 1 << 30, n = 0;
    for (auto& b : B) {
      Net nt; if (!build_net(b.s, sp, outs, nt)) continue;
      int u = (int)nt.g.size(), p = (int)prune(nt).g.size(); n++;
      if (u < umin) { umin = u; best_at_umin = p; } else if (u == umin) best_at_umin = std::min(best_at_umin, p);
      if (p < best_all || (p == best_all && u < umin_best)) { best_all = p; umin_best = u; }
    }
    fprintf(stderr, "STATEND circuits=%d unpruned_min=%d best_pruned_at_unpruned_min=%d best_pruned=%d its_unpruned=%d\n",
            n, umin, best_at_umin, best_all, umin_best);
  }
  for (auto& b : B) consider(best, b.s, sp, outs, hist);
}

// ---------------------------------------------------------------- backward splitting -----------
struct Prefix {
  std::vector<u64> Wv;                           // unsplit outputs (bound s)
  std::vector<std::pair<u64, char>> P;           // (node, is predecessor -> bound s-1, else s)
  std::vector<Split> splits;
  void targets(int s, std::vector<u64>& T, std::vector<int>& TB) const {
    T.clear(); TB.clear();
    std::unordered_map<u64, int> pos;
    auto put = [&](u64 v, int b) {
      auto it = pos.find(v);
      if (it == pos.end()) { pos[v] = (int)T.size(); T.push_back(v); TB.push_back(W(b)); }
      else TB[it->second] = std::min(TB[it->second], W(b));
    };
    for (u64 w : Wv) put(w, s);
    for (auto& p : P) put(p.first, p.second ? s - 1 : s);
  }
  u64 hash() const {
    std::vector<u64> a(Wv); for (auto& p : P) a.push_back(p.first * 2 + (u64)p.second * 0x5bd1e995ULL);
    std::sort(a.begin(), a.end()); u64 h = 1469598103934665603ULL;
    for (u64 x : a) { h ^= x; h *= 1099511628211ULL; h ^= h >> 31; }
    return h;
  }
};

static u64 random_split(u64 w, int s, std::mt19937_64& rng) {    // cancellation-free: w = p1 ^ p2
  std::vector<int> bits; for (int i = 0; i < 64; ++i) if ((w >> i) & 1) bits.push_back(i);
  int n = (int)bits.size(), cap = 1 << (s - 1);
  int lo = std::max(1, n - cap), hi = std::min(cap, n - 1);
  if (lo > hi) return 0;
  int a = std::uniform_int_distribution<int>(lo, hi)(rng);
  std::shuffle(bits.begin(), bits.end(), rng);
  u64 p1 = 0; for (int i = 0; i < a; ++i) p1 |= 1ULL << bits[i];
  return p1;
}

// ---- backward splitting helpers (priority (3) of Shi-Feng-Xu: Rule 1 > Rule 2 > Rule 3, Rule 5)
static int inP(const Prefix& pf, u64 v) { for (size_t i = 0; i < pf.P.size(); ++i) if (pf.P[i].first == v) return (int)i; return -1; }
static void add_pred(Prefix& pf, u64 v) { int i = inP(pf, v); if (i < 0) pf.P.push_back({v, 1}); else pf.P[i].second = 1; }
static void prefix_init(Prefix& pf, const std::vector<u64>& outs, int s) {             // Rule 1
  std::unordered_set<u64> wv;
  for (u64 v : outs) {
    if (dmin_of(v) < s) { if (inP(pf, v) < 0) pf.P.push_back({v, 0}); }
    else if (!wv.count(v)) { wv.insert(v); pf.Wv.push_back(v); }
  }
}
static bool draw_split(const Prefix& pf, int s, double p3, std::mt19937_64& rng, Split& sp, int& rule) {
  if (pf.Wv.empty()) return false;
  std::vector<std::pair<int, int>> r2;                       // Rule 2: w = p1 ^ p2 with p1, p2 in P
  for (int wi = 0; wi < (int)pf.Wv.size(); ++wi)
    for (int pi = 0; pi < (int)pf.P.size(); ++pi) { int q = inP(pf, pf.Wv[wi] ^ pf.P[pi].first); if (q > pi) r2.push_back({wi, pi}); }
  if (!r2.empty()) {
    auto c = r2[std::uniform_int_distribution<int>(0, (int)r2.size() - 1)(rng)];
    sp = {pf.Wv[c.first], pf.P[c.second].first, pf.Wv[c.first] ^ pf.P[c.second].first}; rule = 2; return true;
  }
  std::vector<std::pair<int, int>> r3;                       // Rule 3: w ^ p1 has d_min < s
  for (int wi = 0; wi < (int)pf.Wv.size(); ++wi)
    for (int pi = 0; pi < (int)pf.P.size(); ++pi) { u64 q = pf.Wv[wi] ^ pf.P[pi].first; if (q && dmin_of(q) < s) r3.push_back({wi, pi}); }
  if (!r3.empty() && std::uniform_real_distribution<double>(0, 1)(rng) < p3) {
    auto c = r3[std::uniform_int_distribution<int>(0, (int)r3.size() - 1)(rng)];
    sp = {pf.Wv[c.first], pf.P[c.second].first, pf.Wv[c.first] ^ pf.P[c.second].first}; rule = 3; return true;
  }
  int wi = std::uniform_int_distribution<int>(0, (int)pf.Wv.size() - 1)(rng);   // Rule 5: random split
  u64 p1 = random_split(pf.Wv[wi], s, rng);
  if (!p1) return false;
  sp = {pf.Wv[wi], p1, pf.Wv[wi] ^ p1}; rule = 5; return true;
}
static void apply_split(Prefix& pf, const Split& sp) {
  pf.Wv.erase(std::find(pf.Wv.begin(), pf.Wv.end(), sp.w));
  add_pred(pf, sp.p1); add_pred(pf, sp.p2); pf.splits.push_back(sp);
}

// parse a 2-input circuit in this tool's output format ("gK = A + B [d]" lines, "yJ = name" lines)
static bool load_net(const char* path, const std::vector<u64>& outs, Net& net) {
  std::ifstream f(path); std::string line;
  if (!f) return false;
  std::unordered_map<std::string, int> id;
  for (int c = 0; c < NIN; ++c) id["x" + std::to_string(c)] = NIN - 1 - c;
  std::vector<u64> v(NIN); std::vector<int> d(NIN, 0); for (int i = 0; i < NIN; ++i) v[i] = 1ULL << i;
  net = Net(); std::map<int, std::string> outn;
  while (std::getline(f, line)) {
    char g[64], a[64], b[64]; int j;
    if (sscanf(line.c_str(), " %63[g0-9] = %63[gx0-9] + %63[gx0-9]", g, a, b) == 3 && g[0] == 'g') {
      auto ia = id.find(a), ib = id.find(b);
      if (ia == id.end() || ib == id.end()) return false;
      int k = NIN + (int)net.g.size();
      u64 val = v[ia->second] ^ v[ib->second]; int dd = std::max(d[ia->second], d[ib->second]) + 1;
      net.g.push_back({val, ia->second, ib->second, dd}); v.push_back(val); d.push_back(dd); id[g] = k;
    } else if (sscanf(line.c_str(), " y%d = %63s", &j, a) == 2) outn[j] = a;
  }
  if (outn.size() != outs.size()) return false;
  for (size_t j = 0; j < outs.size(); ++j) { auto it = id.find(outn[(int)j]); if (it == id.end()) return false; net.out.push_back(it->second); }
  return verify(net, outs, DMAX);
}

int main(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: dbs matrix.txt D [--tie bp|sfx|all|sdf|fs] [--md] [--runs R] [--seed S] [--beam W] [--slack s] [--rank sigma|gates] [--stats] [--trace] [--progress]\n"
                    "                        [--bpbs N --evals R --p3 q --topk M] [--obj xor|umc|stm] [--print] [--print3]\n");
    return 1;
  }
  DMAX = atoi(argv[2]); BUDMAX = 1 << DMAX;
  int tie = 0, runs = 1, beam = 0, bpbs = 0, evals = 4, topk = 0, maxsplits = 1 << 30; long slack = 0; double p3 = 0.5;
  const char* lnsfile = nullptr; int iters = 100, kmin = 8, kmax = 8; double plateau = 1.0;
  const char* savepre = nullptr; int btie = -1; bool bmd = false; int bbeam = 0, bopts = 8;
  unsigned long long seed = 0; bool print = false, print3 = false, md = false, verbose = false;
  for (int a = 3; a < argc; ++a) {
    std::string o = argv[a];
    if (o == "--tie") { std::string t = argv[++a]; tie = t == "bp" ? 0 : t == "sfx" ? 1 : t == "all" ? 2 : t == "sdf" ? 3 : 4; }
    else if (o == "--runs") runs = atoi(argv[++a]);
    else if (o == "--seed") seed = strtoull(argv[++a], 0, 10);
    else if (o == "--print") print = true;
    else if (o == "--print3") print3 = true;
    else if (o == "--beam") beam = atoi(argv[++a]);
    else if (o == "--slack") slack = atol(argv[++a]);
    else if (o == "--md") md = true;
    else if (o == "--stats") STATS = true;
    else if (o == "--trace") TRACE = true;
    else if (o == "--progress") PROGRESS = true;
    else if (o == "--rank") { std::string t = argv[++a]; if (t != "sigma" && t != "gates") { fprintf(stderr, "--rank sigma|gates\n"); return 1; } RANK_GATES = t == "gates"; }
    else if (o == "--bpbs") bpbs = atoi(argv[++a]);
    else if (o == "--evals") evals = atoi(argv[++a]);
    else if (o == "--p3") p3 = atof(argv[++a]);
    else if (o == "--topk") topk = atoi(argv[++a]);
    else if (o == "--maxsplits") maxsplits = atoi(argv[++a]);
    else if (o == "--verbose") verbose = true;
    else if (o == "--lns") lnsfile = argv[++a];
    else if (o == "--iters") iters = atoi(argv[++a]);
    else if (o == "--fix") { std::string t = argv[++a]; auto c = t.find(':'); kmin = atoi(t.substr(0, c).c_str()); kmax = c == std::string::npos ? kmin : atoi(t.substr(c + 1).c_str()); }
    else if (o == "--plateau") plateau = atof(argv[++a]);
    else if (o == "--save") savepre = argv[++a];
    else if (o == "--btie") { std::string t = argv[++a]; btie = t == "bp" ? 0 : t == "sfx" ? 1 : t == "all" ? 2 : t == "sdf" ? 3 : 4; }
    else if (o == "--bmd") bmd = true;
    else if (o == "--bbeam") bbeam = atoi(argv[++a]);
    else if (o == "--bopts") bopts = atoi(argv[++a]);
    else if (o == "--obj") { std::string t = argv[++a]; OBJ = t == "xor" ? 0 : t == "umc" ? 1 : 2; }
    else { fprintf(stderr, "unknown option %s\n", o.c_str()); return 1; }
  }
  std::ifstream f(argv[1]); std::string line; std::vector<std::vector<int>> rows;
  while (std::getline(f, line)) { std::istringstream ss(line); std::vector<int> r; int x; while (ss >> x) r.push_back(x); if (!r.empty()) rows.push_back(r); }
  NIN = (int)rows[0].size();
  if (NIN > 64 || rows.size() > 128) { fprintf(stderr, "at most 64 inputs and 128 outputs\n"); return 1; }
  std::vector<u64> outs;
  for (auto& r : rows) { u64 v = 0; for (int c = 0; c < NIN; ++c) if (r[c]) v |= 1ULL << (NIN - 1 - c); if (!v) { fprintf(stderr, "zero row\n"); return 1; } outs.push_back(v); }
  int s = 0; for (u64 v : outs) s = std::max(s, dmin_of(v));
  if (s > DMAX) { fprintf(stderr, "depth bound %d below the minimum depth %d\n", DMAX, s); return 1; }
  std::vector<int> TB0(outs.size(), BUDMAX);
  const char* tn[] = {"bp", "sfx", "all", "sdf", "fs"};
  const char* on[] = {"xor", "umc", "stm"};
  Result best; std::map<int, int> hist;

  if (lnsfile) {
    // large-neighborhood search around an incumbent circuit: keep the top-level splits
    // (output = p1 ^ p2) of k random outputs, re-solve everything else with the forward search
    // (beam or randomized greedy), accept if not worse (equal: with probability --plateau)
    Net inc;
    if (!load_net(lnsfile, outs, inc)) { fprintf(stderr, "cannot load or verify %s\n", lnsfile); return 1; }
    inc = prune(inc); best.net = inc; best.e = evaluate_net(inc); best.ok = true; track3(best.e, inc);
    fprintf(stderr, "incumbent: %d gates, %d XOR2 + %d XOR3\n", best.e.gates, best.e.n2, best.e.n3);
    std::mt19937_64 rng(seed ^ 0xA5A5A5A5ULL);
    for (int it = 0; it < iters; ++it) {
      std::vector<u64> v(NIN + inc.g.size());
      for (int i = 0; i < NIN; ++i) v[i] = 1ULL << i;
      for (size_t k = 0; k < inc.g.size(); ++k) v[NIN + k] = inc.g[k].v;
      std::vector<int> js;                       // outputs computed by a gate (first occurrence of each value)
      std::unordered_set<u64> seenv;
      for (size_t j = 0; j < outs.size(); ++j) if (inc.out[j] >= NIN && !seenv.count(outs[j])) { seenv.insert(outs[j]); js.push_back((int)j); }
      std::shuffle(js.begin(), js.end(), rng);
      int k = std::uniform_int_distribution<int>(kmin, std::max(kmin, kmax))(rng);
      k = std::min(k, (int)js.size());
      Prefix pf; std::unordered_set<u64> fixedv;
      auto inP = [&](u64 x) { for (size_t i = 0; i < pf.P.size(); ++i) if (pf.P[i].first == x) return (int)i; return -1; };
      auto addPred = [&](u64 x) { int i = inP(x); if (i < 0) pf.P.push_back({x, 1}); else pf.P[i].second = 1; };
      for (int t = 0; t < k; ++t) {
        const Gate& G = inc.g[inc.out[js[t]] - NIN];
        pf.splits.push_back({G.v, v[G.a], v[G.b]}); fixedv.insert(G.v);
        addPred(v[G.a]); addPred(v[G.b]);
      }
      std::unordered_set<u64> wv;
      for (u64 u : outs) {
        if (fixedv.count(u) || wv.count(u)) continue;
        wv.insert(u);
        if (dmin_of(u) < s) { if (inP(u) < 0) pf.P.push_back({u, 0}); } else pf.Wv.push_back(u);
      }
      std::vector<u64> T; std::vector<int> TB; pf.targets(s, T, TB);
      Result r;
      if (T.size() <= 128) {
        if (beam > 0) beam_search(T, TB, beam, slack, seed * 1000003ULL + it, tie, md, pf.splits, outs, r, nullptr);
        else greedy_runs(T, TB, evals, seed * 1000003ULL + it, tie, md, pf.splits, outs, r, nullptr, false);
      }
      const char* what = "rejected";
      if (r.ok && r.e.key < best.e.key - 1e-9) { what = "IMPROVED"; best = r; inc = r.net; }
      else if (r.ok && r.e.key < best.e.key + 1e-9 && std::uniform_real_distribution<double>(0, 1)(rng) < plateau) { what = "plateau move"; best = r; inc = r.net; }
      if (r.ok) hist[r.e.gates]++;
      if (verbose || what[0] == 'I')
        fprintf(stderr, "iter %d: fixed %d outputs, %zu targets -> %d gates (%d+%d) %s; best %d gates (%d XOR2 + %d XOR3)\n",
                it, k, T.size(), r.ok ? r.e.gates : -1, r.ok ? r.e.n2 : 0, r.ok ? r.e.n3 : 0, what, best.e.gates, best.e.n2, best.e.n3);
    }
  } else if (!bpbs) {
    std::vector<Split> none;
    if (beam > 0) for (int r = 0; r < runs; ++r) beam_search(outs, TB0, beam, slack, seed, tie, md, none, outs, best, &hist);
    else greedy_runs(outs, TB0, runs, seed, tie, md, none, outs, best, &hist);
  } else {
    // backward splitting of the outputs at the top level (Shi-Feng-Xu Algorithm 4, priority (3))
    if (DMAX != s) { fprintf(stderr, "--bpbs expects D = minimum depth (%d)\n", s); return 1; }
    std::mt19937_64 rng(seed ^ 0x9E3779B97F4A7C15ULL);
    std::vector<Split> none; greedy_runs(outs, TB0, evals, seed, tie, md, none, outs, best, &hist);
    std::vector<std::pair<double, Prefix>> pool; std::unordered_set<u64> seenp;
    long evalsDone = 0;
    // evaluate a prefix with randomized forward search; update best, histogram, prefix pool
    auto evaluate = [&](const Prefix& pf, u64 sd, Result& r) {
      std::vector<u64> T; std::vector<int> TB; pf.targets(s, T, TB);
      if (T.size() > 128) return false;
      greedy_runs(T, TB, evals, sd, tie, md, pf.splits, outs, r, nullptr, false);
      evalsDone += evals;
      if (!r.ok) return false;
      if (!best.ok || r.e.key < best.e.key) { best = r; fprintf(stderr, "%zu splits: %d gates, %d XOR2 + %d XOR3\n", pf.splits.size(), r.e.gates, r.e.n2, r.e.n3); }
      hist[r.e.gates]++;
      if (topk > 0) { u64 h = pf.hash(); if (!seenp.count(h)) { seenp.insert(h); pool.push_back({r.e.key, pf}); } }
      if (topk > 0 && pool.size() > (size_t)(4 * topk)) {
        std::sort(pool.begin(), pool.end(), [](auto& x, auto& y) { return x.first < y.first; });
        pool.resize(topk);
      }
      return true;
    };
    if (bbeam <= 0) {
      for (int tr = 0; tr < bpbs; ++tr) {                       // independent random trajectories
        Prefix pf; prefix_init(pf, outs, s);
        while (!pf.Wv.empty() && (int)pf.splits.size() < maxsplits) {
          Split sp; int rule;
          if (!draw_split(pf, s, p3, rng, sp, rule)) break;
          apply_split(pf, sp);
          Result r; bool ok = evaluate(pf, seed * 7919 + tr * 131 + pf.splits.size(), r);
          if (verbose) fprintf(stderr, "  traj %d split %zu (rule %d): %d gates\n", tr, pf.splits.size(), rule, ok ? r.e.gates : -1);
        }
      }
    } else {
      // beam search over backward splits: every prefix draws --bopts split options,
      // each child is evaluated by --evals greedy runs, the --bbeam best children survive
      for (int rep = 0; rep < bpbs; ++rep) {
        Prefix p0; prefix_init(p0, outs, s);
        std::vector<std::pair<double, Prefix>> cur(1, {1e18, p0});
        for (int step = 0; step < maxsplits; ++step) {
          std::vector<std::pair<double, Prefix>> nxt; std::unordered_set<u64> seen;
          for (auto& cp : cur) {
            for (int o = 0; o < bopts; ++o) {
              Split sp; int rule;
              if (!draw_split(cp.second, s, p3, rng, sp, rule)) continue;
              Prefix ch = cp.second; apply_split(ch, sp);
              u64 h = ch.hash(); if (seen.count(h)) continue; seen.insert(h);
              Result r;
              if (evaluate(ch, seed * 7919 + rep * 1000003ULL + step * 131 + nxt.size(), r)) nxt.push_back({r.e.key, ch});
            }
          }
          if (nxt.empty()) break;
          std::sort(nxt.begin(), nxt.end(), [](auto& x, auto& y) { return x.first < y.first; });
          if ((int)nxt.size() > bbeam) nxt.resize(bbeam);
          if (verbose) fprintf(stderr, "  rep %d step %d: %zu prefixes, best %.2f\n", rep, step + 1, nxt.size(), nxt[0].first);
          cur.swap(nxt);
        }
      }
    }
    fprintf(stderr, "bpbs: %d %s, %ld greedy evaluations\n", bpbs, bbeam > 0 ? "backward beams" : "trajectories", evalsDone);
    if (topk > 0 && beam > 0) {                                      // beam refinement of the best prefixes
      std::sort(pool.begin(), pool.end(), [](auto& x, auto& y) { return x.first < y.first; });
      if ((int)pool.size() > topk) pool.resize(topk);
      for (size_t i = 0; i < pool.size(); ++i) {
        std::vector<u64> T; std::vector<int> TB; pool[i].second.targets(s, T, TB);
        Result r; beam_search(T, TB, beam, slack, seed + i, btie >= 0 ? btie : tie, btie >= 0 ? bmd : md, pool[i].second.splits, outs, r, nullptr);
        fprintf(stderr, "beam on prefix %zu (%zu splits, greedy %.2f): %s\n", i, pool[i].second.splits.size(), pool[i].first,
                r.ok ? (std::to_string(r.e.gates) + " gates").c_str() : "none");
        if (r.ok && r.e.key < best.e.key) best = r;
      }
    }
  }
  if (!best.ok) { printf("RESULT none\n"); return 0; }
  printf("RESULT best=%d depth<=%d tie=%s%s runs=%d seed=%llu obj=%s xor2/xor3 after merging: %d + %d, area umc180 %.2f stm130 %.2f\n",
         best.e.gates, DMAX, tn[tie], md ? "+md" : "", runs, seed, on[OBJ], best.e.n2, best.e.n3,
         area_of(best.e.n2, best.e.n3, 1), area_of(best.e.n2, best.e.n3, 2));
  printf("histogram:"); for (auto& kv : hist) printf(" %d:%d", kv.first, kv.second); printf("\n");
  if (BEST3[0].ok) {
    const char* tag[3] = {"xor2", "umc180", "stm130"};
    for (int i = 0; i < 3; ++i) {
      const Eval& e = BEST3[i].e;
      printf("BEST-%s: %d XOR2 gates -> %d XOR2 + %d XOR3, umc180 %.2f GE, stm130 %.2f GE\n", tag[i], e.gates, e.n2, e.n3,
             area_of(e.n2, e.n3, 1), area_of(e.n2, e.n3, 2));
      if (savepre) {
        std::string fn = std::string(savepre) + "." + tag[i] + ".txt";
        FILE* o = fopen(fn.c_str(), "w");
        if (o) { FILE* keep = stdout; stdout = o; print_net(BEST3[i].net, false); fflush(o); stdout = keep; fclose(o); }
      }
    }
  }
  if (print) print_net(best.net, false);
  if (print3) print_net(best.net, true);
  return 0;
}
