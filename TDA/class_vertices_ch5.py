# Class vertex sizes and tip/node split (Table 5.1, Section 5.2).
import csv
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "networks"
SEEDS = ["s117", "s120", "s175", "s25", "s179", "s187"]
TOL = 1e-6


def load(d):
    nt = {}
    with io.open(d / "net_name_table.csv", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["species"] != "NA":
                nt[int(r["species"])] = int(r["phy_node"])
    with io.open(d / "dists2.csv", encoding="utf-8") as f:
        hdr = next(csv.reader(f))
        idx2node = [nt[int(h)] for h in hdr[1:]]
    tips = {int(r["phy_node"])
            for r in csv.DictReader(io.open(d / "net_extant_tips.csv",
                                            encoding="utf-8"))}
    return idx2node, tips


cyc_rows = []
for seed in SEEDS:
    d = BASE / seed
    idx2node, tips = load(d)
    log = json.load(io.open(d / "eve_full_rip.json", encoding="utf-8"))
    for cyc in log["harmonic_cycles"]:
        V = set()
        for e in cyc["edges"]:
            w = float(e["weight"])
            if abs(w) <= TOL:
                continue
            s = e["simplex"]
            if len(s) != 2:
                continue
            V.add(idx2node[int(s[0])])
            V.add(idx2node[int(s[1])])
        t = sorted(V & tips)
        n = sorted(V - tips)
        cyc_rows.append([seed, cyc["cycle_index"], cyc["birth"], cyc["death"],
                         len(V), len(t), len(n),
                         "{%s}" % ",".join(map(str, t)),
                         "{%s}" % ",".join(map(str, n))])

with open(BASE / "event_cycle_vertices.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["seed", "cycle_index", "birth", "death", "n_vertices",
                "n_tips", "n_nodes", "tip_labels", "node_labels"])
    w.writerows(cyc_rows)

print("cycles:", len(cyc_rows))
print("per-seed class vertex sizes (min/max):")
for seed in SEEDS:
    sizes = sorted(r[4] for r in cyc_rows if r[0] == seed)
    print("  %s: %s" % (seed, sizes))
