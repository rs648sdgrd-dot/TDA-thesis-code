# 7/8 linkage check (Section 5.4); needs Eve's networks-lab-tda (github.com/EveWCheng/networks-lab-tda).
import csv
import json
import sys
from pathlib import Path

import numpy as np

import os
EVE_SRC = os.environ.get("EVE_SRC")
if not EVE_SRC:
    raise SystemExit("Set EVE_SRC to the src directory of Eve's networks-lab-tda.")
sys.path.insert(0, EVE_SRC)
from network_lab_tda.tda_analysis import harmonic_cycle  # noqa: E402

BASE = Path(__file__).resolve().parents[1] / "data" / "networks"
TMP = Path(__file__).resolve().parents[1] / "tmp"
TMP.mkdir(exist_ok=True)
TOL = 1e-6


def load_seed(seed):
    d = BASE / seed
    with open(d / "dists2.csv", newline="") as f:
        rd = csv.reader(f)
        hdr = next(rd)[1:]
        M = np.array([[float(x) for x in r[1:]] for r in rd])
    nt = {}
    with open(d / "net_name_table.csv") as f:
        for r in csv.DictReader(f):
            if r["species"] != "NA":
                nt[int(r["species"])] = int(r["phy_node"])
    node2idx = {nt[int(h)]: i for i, h in enumerate(hdr)}
    return M, node2idx


def run_full(M):
    hc = harmonic_cycle(M, cycle_dim=1, sim_log=True,
                        log_path=str(TMP / "mr3_rip.json"))
    result = hc.run_harmonics(threshold=float("inf"), save=True)
    int_cycles = result.homology.persistent_cycles.get(1, {})
    bars = result.homology.bars.get(1, {})
    simplices = hc.log["simplicies"]
    int_edge = {}
    for cid, chain in int_cycles.items():
        for sid, coeff in chain.items():
            s = simplices[int(sid)]
            if len(s) == 2:
                key = (min(s[0], s[1]), max(s[0], s[1]))
                int_edge.setdefault(key, []).append(cid)
    # harmonic representative edges per (birth,death) class
    harm_by_bd = {}
    for cyc in hc.log.get("harmonic_cycles", []):
        b, dd = cyc["birth"], cyc["death"]
        edges = []
        for e in cyc["edges"]:
            w = float(e["weight"])
            if abs(w) <= TOL:
                continue
            s = e["simplex"]
            if len(s) == 2:
                edges.append(((min(int(s[0]), int(s[1])),
                               max(int(s[0]), int(s[1]))), w))
        harm_by_bd[(round(b, 6), round(dd, 6))] = edges
    return int_cycles, bars, int_edge, harm_by_bd


pairs = list(csv.DictReader(open(BASE / "interp_v2_pairs.csv")))
fp_ctrl = {("s187", "14", "20")}
sel = [p for p in pairs
       if p["type"] == "reticulate"
       or (p["seed"], p["tip_i"], p["tip_j"]) in fp_ctrl]

rows = []
for p in sel:
    seed = p["seed"]
    M, node2idx = load_seed(seed)
    i, j = int(p["tip_i"]), int(p["tip_j"])
    ti, tj = node2idx[i], node2idx[j]
    d_ij = M[ti, tj]
    n = M.shape[0]
    Mm = np.zeros((n + 1, n + 1))
    Mm[:-1, :-1] = M
    row_m = 0.5 * (M[ti] + M[tj])
    row_m[ti] = row_m[tj] = d_ij / 2.0
    Mm[-1, :-1] = row_m
    Mm[:-1, -1] = row_m
    Mm[-1, -1] = 0.0
    int_cycles, bars, int_edge, harm_by_bd = run_full(Mm)
    e_im = (min(ti, n), max(ti, n))
    e_mj = (min(n, tj), max(n, tj))
    resp_cids = sorted(set(int_edge.get(e_im, []) + int_edge.get(e_mj, [])))
    pers = [(round(bars[c][0], 3), round(bars[c][1], 3),
             round(bars[c][1] - bars[c][0], 3)) for c in resp_cids]


    evs = []
    if p["type"] == "reticulate":
        with open(BASE / seed / "pair_registry_v2.csv") as f:
            for r in csv.DictReader(f):
                if (int(r["tip_i"]), int(r["tip_j"])) == (i, j):
                    evs = [int(x) for x in r["events"].split(";")]
    with open(BASE / seed / "net_retic_edges.csv") as f:
        r_edges = [r for r in csv.DictReader(f) if r["from"] != ""]
    ret_edges = []
    for eid in evs:
        a, b = int(r_edges[eid - 1]["from"]), int(r_edges[eid - 1]["to"])
        ret_edges.append((min(node2idx[a], node2idx[b]),
                          max(node2idx[a], node2idx[b])))


    link = []
    for cid in resp_cids:
        b, dd = bars[cid]
        hrep = harm_by_bd.get((round(b, 6), round(dd, 6)), [])
        found = [(ek, round(w, 4)) for ek, w in hrep if ek in ret_edges]
        link.append("c%d(%.3f,%.3f):%s" % (cid, b, dd,
                                           ";".join("%s@%s" % x
                                                    for x in found)
                                           if found else "none"))

    # cycle counts for the target reticulation edges
    before = {}
    with open(BASE / "event_cycle_direct.csv") as f:
        for r in csv.DictReader(f):
            if r["seed"] == seed and int(r["event_id"]) in evs:
                before.setdefault(int(r["event_id"]), set()).add(
                    r["cycle_index"])
    after = {}
    for bd, hedges in harm_by_bd.items():
        for (ek, w) in hedges:
            for k, re_ in enumerate(ret_edges):
                if ek == re_:
                    after.setdefault(evs[k], set()).add(bd)
    n_link_hits = sum(1 for l in link if l.split(":")[1] != "none")
    rows.append([seed, p["type"], i, j, evs, len(resp_cids),
                 ";".join("%s" % (p_,) for p_ in pers), n_link_hits,
                 "|".join(link),
                 ";".join("%d:%d->%d" % (e, len(before.get(e, ())),
                                         len(after.get(e, ())))
                          for e in evs)])
    print("%s %-10s (%d,%d) ev=%s n_resp_cyc=%d pers=%s n_link=%d"
          % (seed, p["type"], i, j, evs, len(resp_cids), pers, n_link_hits))
    for l in link:
        print("    %s" % l)

out = BASE / "mr3_linkage.csv"
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["seed", "type", "tip_i", "tip_j", "events", "n_resp_cyc",
                "persistence", "n_link_hits", "linkage", "before_after"])
    w.writerows(rows)
print("saved", out)
