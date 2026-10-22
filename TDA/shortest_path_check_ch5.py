# Shortest-path checks: 11/12 and 22.8 (Section 5.4).
import csv
import heapq
import io
from pathlib import Path

BASE = Path(__file__).resolve().parents[1] / "data" / "networks"
SEEDS = ["s117", "s120", "s175", "s25", "s179", "s187"]
TOL = 1e-6


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


def on_path(net_adj, i, j, a, b, L, dn):
    di = dijkstra(net_adj, i)
    dj = dijkstra(net_adj, j)
    return (abs(di[a] + L + dj[b] - dn) < 1e-6 or
            abs(di[b] + L + dj[a] - dn) < 1e-6)


pairs = list(csv.DictReader(io.open(BASE / "interp_v2_pairs.csv")))
out = []
for seed in SEEDS:
    d = BASE / seed
    nt = {}
    for r in csv.DictReader(io.open(d / "net_name_table.csv")):
        if r["species"] != "NA":
            nt[int(r["species"])] = int(r["phy_node"])
    node2idx = {v: k for k, v in nt.items()}
    tree_adj, net_adj = {}, {}
    for r in csv.DictReader(io.open(d / "net_tree_edges.csv")):
        a, b, L = int(r["from"]), int(r["to"]), float(r["len"])
        tree_adj.setdefault(a, []).append((b, L))
        tree_adj.setdefault(b, []).append((a, L))
        net_adj.setdefault(a, []).append((b, L))
        net_adj.setdefault(b, []).append((a, L))
    retic = []
    for r in csv.DictReader(io.open(d / "net_retic_edges.csv")):
        if r["from"] == "":
            continue
        a, b, L = int(r["from"]), int(r["to"]), float(r["len"])
        net_adj.setdefault(a, []).append((b, L))
        net_adj.setdefault(b, []).append((a, L))
        retic.append((a, b, L))
    reg = {}
    for r in csv.DictReader(io.open(d / "pair_registry_v2.csv")):
        reg[(int(r["tip_i"]), int(r["tip_j"]))] = (int(r["r"]), r["events"])

    dists2 = {}
    with io.open(d / "dists2.csv") as f:
        rd = csv.reader(f)
        hdr = next(rd)[1:]
        keys = [nt[int(h)] for h in hdr]
        for r in rd:
            u = nt[int(r[0])]
            for k, v in zip(keys, r[1:]):
                dists2[(u, k)] = float(v)
    for p in pairs:
        if p["seed"] != seed:
            continue
        i, j = int(p["tip_i"]), int(p["tip_j"])
        dt = dijkstra(tree_adj, i)[j]
        dn = dijkstra(net_adj, i)[j]
        hits = [(a, b) for a, b, L in retic
                if on_path(net_adj, i, j, a, b, L, dn)]
        rv, evs = reg.get((i, j), (0, ""))
        target_hit = "na"
        if p["type"] == "reticulate" and evs:
            eid = int(evs.split(";")[0]) - 1
            target_hit = ("yes" if (retic[eid][0], retic[eid][1]) in hits
                          else "no")
        d2 = dists2.get((i, j), float("nan"))
        out.append(dict(
            seed=seed, type=p["type"], tip_i=i, tip_j=j, r=rv,
            d_tree=round(dt, 4), d_net=round(dn, 4),
            shortening=round(dt - dn, 4),
            dists2=round(d2, 4),
            d2_minus_dnet=round(d2 - dn, 4),
            n_retic_on_path=len(hits),
            retic_on_path=";".join("%d-%d" % h for h in hits) or "none",
            target_retic_on_path=target_hit))

with open(BASE / "shortest_path_results.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(list(out[0].keys()))
    for r in out:
        w.writerow([r[k] for k in out[0]])

print("%-5s %-10s %-8s %9s %9s %9s %10s %4s %s" % (
    "seed", "type", "pair", "d_tree", "d_net", "shorten", "d2-d_net",
    "#ret", "target_hit"))
for r in out:
    print("%-5s %-10s %-8s %9.3f %9.3f %9.3f %10.4f %4d %s" % (
        r["seed"], r["type"], "%d,%d" % (r["tip_i"], r["tip_j"]),
        r["d_tree"], r["d_net"], r["shortening"], r["d2_minus_dnet"],
        r["n_retic_on_path"], r["target_retic_on_path"]))

ret = [r for r in out if r["type"] == "reticulate"]
con = [r for r in out if r["type"] == "control"]
print("\ntarget retic edge on shortest path (r=1 pairs): %d/%d"
      % (sum(1 for r in ret if r["target_retic_on_path"] == "yes"), len(ret)))
print("any retic edge on shortest path (r=1 pairs): %d/%d"
      % (sum(1 for r in ret if r["n_retic_on_path"] > 0), len(ret)))
print("any retic edge on shortest path (controls):    %d/%d"
      % (sum(1 for r in con if r["n_retic_on_path"] > 0), len(con)))
print("median shortening: r=1 %.4f  controls %.4f"
      % (sorted(r["shortening"] for r in ret)[len(ret) // 2],
         sorted(r["shortening"] for r in con)[len(con) // 2]))
print("max |dists2 - d_net| over all pairs: %.4f"
      % max(abs(r["d2_minus_dnet"]) for r in out))
