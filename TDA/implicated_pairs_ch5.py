# r-values of implicated tip pairs (Section 5.3).
import csv
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "networks"
SEEDS = ["s117", "s120", "s175", "s25", "s179", "s187"]
TOL = 1e-6


def load_network(d):
    tree = []
    with open(d / "net_tree_edges.csv") as f:
        for r in csv.DictReader(f):
            tree.append((int(r["from"]), int(r["to"])))
    retic = []
    with open(d / "net_retic_edges.csv") as f:
        for r in csv.DictReader(f):
            if r["from"] == "":
                continue
            retic.append((int(r["from"]), int(r["to"])))
    tips = []
    with open(d / "net_extant_tips.csv") as f:
        for r in csv.DictReader(f):
            tips.append(int(r["phy_node"]))
    return tree, retic, tips


def desc_tips(tree, tips):

    children = defaultdict(list)
    nodes = set(tips)
    for a, b in tree:
        children[a].append(b)
        nodes.update((a, b))
    tip_set = set(tips)
    out = {}
    for x in nodes:
        seen = []
        stack = list(children[x])
        while stack:
            v = stack.pop()
            if v in tip_set:
                seen.append(v)
            stack.extend(children[v])
        out[x] = seen
    return out


rows = []
match_rows = []
print("%-5s %-8s %-8s %-12s %-10s" % ("seed", "tips", "retic", "pairs_r1",
                                        "pairs_r>=2"))
for seed in SEEDS:
    d = BASE / seed
    tree, retic, tips = load_network(d)
    DT = desc_tips(tree, tips)

    r_count = defaultdict(int)
    ev_of_pair = defaultdict(list)
    for e, (p_e, h_e) in enumerate(retic, 1):
        A, B = DT.get(p_e, []), DT.get(h_e, [])
        for i in A:
            for j in B:
                if i == j:
                    continue
                p = (min(i, j), max(i, j))
                r_count[p] += 1
                ev_of_pair[p].append(e)
    n_r1 = sum(1 for v in r_count.values() if v == 1)
    n_r2 = sum(1 for v in r_count.values() if v >= 2)
    print("%-5s %-8d %-8d %-12d %-10d" % (seed, len(tips), len(retic),
                                          n_r1, n_r2))
    with open(d / "pair_registry_v2.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tip_i", "tip_j", "r", "events"])
        for (i, j), r in sorted(r_count.items()):
            w.writerow([i, j, r, ";".join(map(str, ev_of_pair[(i, j)]))])


    with open(d / "eve_full_rip.json") as f:
        log = json.load(f)
    with open(d / "dists2.csv", newline="") as f:
        rd = csv.reader(f)
        hdr = next(rd)[1:]
    nt = {}
    with open(d / "net_name_table.csv") as f:
        for r in csv.DictReader(f):
            if r["species"] != "NA":
                nt[int(r["species"])] = int(r["phy_node"])
    idx2node = [nt[int(h)] for h in hdr]
    tip_set = set(tips)
    n_events_hit = set()
    n_pairs = 0
    for cyc in log["harmonic_cycles"]:
        for edge in cyc["edges"]:
            w = float(edge["weight"])
            if abs(w) <= TOL:
                continue
            s = edge["simplex"]
            if len(s) != 2:
                continue
            ni, nj = idx2node[int(s[0])], idx2node[int(s[1])]
            if ni not in tip_set or nj not in tip_set:
                continue
            p = (min(ni, nj), max(ni, nj))
            if p in r_count and r_count[p] >= 1:
                n_pairs += 1
                for e in ev_of_pair[p]:
                    n_events_hit.add(e)
                match_rows.append([seed, cyc["cycle_index"], p[0], p[1],
                                   w, r_count[p],
                                   ";".join(map(str, ev_of_pair[p]))])
    print("  T_e tip-pair matching: %d events hit (of %d), %d pair-incidences"
          % (len(n_events_hit), len(retic), n_pairs))

with open(BASE / "event_cycle_matches_genealogical.csv", "w",
          newline="") as f:
    w = csv.writer(f)
    w.writerow(["seed", "cycle", "tip_i", "tip_j", "harmonic_weight", "r",
                "events"])
    w.writerows(match_rows)
print("\nsaved: pair_registry_v2.csv per seed; "
      "event_cycle_matches_genealogical.csv")
