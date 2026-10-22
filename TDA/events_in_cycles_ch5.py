# Events with non-zero harmonic weight (Table 5.1).
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "networks"
SEEDS = ["s117", "s120", "s175", "s25", "s179", "s187"]
TOL = 1e-6


def load_maps(d):
    with open(d / "dists2.csv", newline="") as f:
        rd = csv.reader(f)
        hdr = next(rd)[1:]
        M = [[float(x) for x in r[1:]] for r in rd]
    keys = [int(h) for h in hdr]
    nt = {}
    with open(d / "net_name_table.csv") as f:
        for r in csv.DictReader(f):
            if r["species"] != "NA":
                nt[int(r["species"])] = int(r["phy_node"])
    idx2node = [nt[k] for k in keys]
    retic = []
    with open(d / "net_retic_edges.csv") as f:
        rows = [r for r in csv.DictReader(f) if r["from"] != ""]
    inher = {}
    with open(d / "inher.csv") as f:
        for r in csv.DictReader(f):
            inher[int(r["k"])] = float(r["inher"])
    for k, r in enumerate(rows, 1):
        retic.append((k, frozenset((int(r["from"]), int(r["to"]))),
                      inher.get(k, float("nan"))))
    return M, idx2node, retic


all_rows = []
print("%-5s %-9s %-12s %-22s %-10s %-7s %s" % (
    "seed", "cycle", "birth-death", "retic edge (event)", "h_weight",
    "rank", "closing"))
for seed in SEEDS:
    d = BASE / seed
    M, idx2node, retic = load_maps(d)
    with open(d / "eve_full_rip.json") as f:
        log = json.load(f)
    seed_rows = []
    for cyc in log["harmonic_cycles"]:
        edges = []
        for e in cyc["edges"]:
            w = float(e["weight"])
            if abs(w) <= TOL:
                continue
            s = e["simplex"]
            if len(s) != 2:
                continue
            i, j = int(s[0]), int(s[1])
            ek = frozenset((idx2node[i], idx2node[j]))
            edges.append((ek, w, M[i][j]))
        distinct = sorted({a for _, _, a in edges})
        rank = {t: r for r, t in enumerate(distinct, start=1)}
        last = distinct[-1] if distinct else None
        for eid, ek, inher_w in retic:
            hit = next(((w, a) for e2, w, a in edges if e2 == ek), None)
            if hit is None:
                continue
            w, a = hit
            row = [seed, cyc["cycle_index"], cyc["birth"], cyc["death"],
                   eid, ek, w, inher_w, rank[a], len(distinct),
                   a == last]
            seed_rows.append(row)
            all_rows.append(row)
            print("%-5s %-9s %-12s %-22s %-10.5f %-7s %s"
                  % (seed, cyc["cycle_index"],
                     "%.3f-%.3f" % (cyc["birth"], cyc["death"]),
                     "%s (R%d)" % (tuple(sorted(ek)), eid), w,
                     "%d/%d" % (rank[a], len(distinct)),
                     "CLOSING" if a == last else ""))


    ev_cyc = {}
    for row in seed_rows:
        ev_cyc.setdefault(row[4], set()).add(row[1])
    shared = {e: cs for e, cs in ev_cyc.items() if len(cs) > 1}
    print("  %s: %d/%d events in cycles, shared events: %s"
          % (seed, len(ev_cyc), len(retic),
             {e: sorted(cs) for e, cs in shared.items()} if shared else "-"))

with open(BASE / "event_cycle_direct.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["seed", "cycle_index", "birth", "death", "event_id",
                "retic_edge", "harmonic_weight", "inher_weight",
                "eve_rank", "n_distinct", "is_closing"])
    for r in all_rows:
        w.writerow([r[0], r[1], r[2], r[3], r[4], "%d-%d" % tuple(sorted(r[5])),
                    r[6], r[7], r[8], r[9], r[10]])


print("\n=== SUMMARY ===")
n_ev = n_hit = 0
shared_count = 0
max_w = 0.0
for seed in SEEDS:
    d = BASE / seed
    M, idx2node, retic = load_maps(d)
    with open(d / "eve_full_rip.json") as f:
        log = json.load(f)
    ev_cyc = {}
    wmax = 0.0
    for row in all_rows:
        if row[0] != seed:
            continue
        ev_cyc.setdefault(row[4], set()).add(row[1])
        wmax = max(wmax, abs(row[6]))
    n_ev += len(retic)
    n_hit += len(ev_cyc)
    shared_count += sum(1 for cs in ev_cyc.values() if len(cs) > 1)
    max_w = max(max_w, wmax)
print("events total: %d" % n_ev)
print("events directly present in >=1 harmonic cycle: %d (%.1f%%)"
      % (n_hit, 100.0 * n_hit / max(n_ev, 1)))
print("events shared across >=2 harmonic cycles: %d" % shared_count)
print("max |harmonic weight| on a retic edge: %.4f" % max_w)
