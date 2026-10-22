# Top-3 candidate screening and path validation (Table 5.3, Section 5.3).
import csv
import heapq
import io
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "networks"
SEEDS = ["s117", "s120", "s175", "s25", "s179", "s187"]
TOL = 1e-6
RTOL = 1e-12


def dijkstra(adj, src):
    dist = {src: 0.0}
    pq = [(0.0, src)]
    while pq:
        d, v = heapq.heappop(pq)
        if d > dist[v]:
            continue
        for w, wt in adj.get(v, []):
            nd = d + wt
            if nd < dist.get(w, float("inf")):
                dist[w] = nd
                heapq.heappush(pq, (nd, w))
    return dist


edges_by_cyc = defaultdict(list)
for r in csv.DictReader(io.open(BASE / "tip_edge_weights.csv")):
    edges_by_cyc[(r["seed"], int(r["cycle"]))].append(
        (abs(float(r["w"])), float(r["w"]), int(r["u"]), int(r["v"])))

candidates = []
for (seed, cyc), es in edges_by_cyc.items():
    es.sort(key=lambda x: -x[0])
    cutoff = es[2][0] if len(es) >= 3 else 0.0
    order = 0
    for a, w, u, v in es:
        if a < cutoff * (1.0 - RTOL):
            break
        order += 1
        candidates.append((seed, cyc, u, v, w, order))

rows = []
hits = defaultdict(set)
net_cache = {}
for seed in SEEDS:
    d = BASE / seed
    tree_keys = set()
    tree_adj, net_adj = {}, {}
    for r in csv.DictReader(io.open(d / "net_tree_edges.csv")):
        a, b, L = int(r["from"]), int(r["to"]), float(r["len"])
        tree_keys.add(frozenset((a, b)))
        tree_adj.setdefault(a, []).append((b, L))
        tree_adj.setdefault(b, []).append((a, L))
        net_adj.setdefault(a, []).append((b, L))
        net_adj.setdefault(b, []).append((a, L))
    retic = []
    retic_keys = set()
    for r in csv.DictReader(io.open(d / "net_retic_edges.csv")):
        if r["from"] == "":
            continue
        a, b, L = int(r["from"]), int(r["to"]), float(r["len"])
        retic_keys.add(frozenset((a, b)))
        net_adj.setdefault(a, []).append((b, L))
        net_adj.setdefault(b, []).append((a, L))
        retic.append((a, b, L))
    assert not (tree_keys & retic_keys), \
        "tree and reticulation edges share endpoints in %s" % seed
    net_cache[seed] = (net_adj, retic)

    for (s2, cyc, i, j, w, order) in [c for c in candidates if c[0] == seed]:
        di = dijkstra(net_adj, i)
        dj = dijkstra(net_adj, j)
        dn = di[j]
        hit_edges = []
        for k, (a, b, L) in enumerate(retic, 1):
            if (abs(di[a] + L + dj[b] - dn) < 1e-6 or
                    abs(di[b] + L + dj[a] - dn) < 1e-6):
                hit_edges.append(k)
                hits[(seed, k)].add((i, j))
        rows.append(dict(seed=seed, cycle=cyc, tip_i=i, tip_j=j,
                         harmonic_weight=w, selection_order=order,
                         n_retic_edges=len(hit_edges),
                         retic_edges=";".join("R%d" % k for k in hit_edges)
                         or "none"))

with open(BASE / "candidate_paths_before.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["seed", "cycle", "tip_i", "tip_j", "harmonic_weight",
                "selection_order", "n_retic_edges", "retic_edges"])
    for r in rows:
        w.writerow([r["seed"], r["cycle"], r["tip_i"], r["tip_j"],
                    r["harmonic_weight"], r["selection_order"],
                    r["n_retic_edges"], r["retic_edges"]])

uniq_pairs = {(r["seed"], r["tip_i"], r["tip_j"]) for r in rows}
pair_hit = {p for p in uniq_pairs if any(
    r["n_retic_edges"] > 0 and (r["seed"], r["tip_i"], r["tip_j"]) == p
    for r in rows)}
n_ev = sum(len([r for r in csv.DictReader(io.open(
    BASE / s / "net_retic_edges.csv")) if r["from"] != ""]) for s in SEEDS)
ev_hit = set(hits.keys())
print("candidate incidences: %d | unique candidate pairs: %d"
      % (len(rows), len(uniq_pairs)))
print("candidate pairs whose shortest path contains >=1 retic edge: "
      "%d/%d (%.1f%%)" % (len(pair_hit), len(uniq_pairs),
                          100.0 * len(pair_hit) / len(uniq_pairs)))
print("true reticulation events hit by >=1 candidate path: %d/%d (%.1f%%)"
      % (len(ev_hit), n_ev, 100.0 * len(ev_hit) / n_ev))
print("per seed (events hit / events):")
for s in SEEDS:
    n_s = len([r for r in csv.DictReader(io.open(
        BASE / s / "net_retic_edges.csv")) if r["from"] != ""])
    h_s = sum(1 for (s2, k) in ev_hit if s2 == s)
    print("  %s: %d/%d" % (s, h_s, n_s))
print("per seed (unique candidate pairs / hit pairs):")
for s in SEEDS:
    ps = [p for p in uniq_pairs if p[0] == s]
    hs = [p for p in ps if p in pair_hit]
    print("  %s: %d pairs, %d hit" % (s, len(ps), len(hs)))
