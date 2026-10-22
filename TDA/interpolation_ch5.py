# Midpoint interpolation experiment (Section 5.4).
import csv
import os
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
# Eve's networks-lab-tda lives next to this folder;
# override with EVE_SRC if it lives elsewhere.
EVE_SRC = os.environ.get("EVE_SRC", str(ROOT.parent / "eve_code" / "networks-lab-tda" / "src"))
sys.path.insert(0, EVE_SRC)
from network_lab_tda.tda_analysis import harmonic_cycle

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "networks"
TOL = 1e-6


def load_seed(seed):
    d = BASE / seed
    with open(d / "dists2.csv", newline="") as f:
        rd = csv.reader(f)
        hdr = next(rd)[1:]
        M = np.array([[float(x) for x in r[1:]] for r in rd])
    keys = [int(h) for h in hdr]
    nt = {}
    with open(d / "net_name_table.csv") as f:
        for r in csv.DictReader(f):
            if r["species"] != "NA":
                nt[int(r["species"])] = int(r["phy_node"])
    idx2node = [nt[k] for k in keys]
    node2idx = {n: i for i, n in enumerate(idx2node)}
    return M, node2idx, nt, keys


def run_full(M, log_path):
    hc = harmonic_cycle(M, cycle_dim=1, sim_log=True,
                        log_path=log_path)
    result = hc.run_harmonics(threshold=float("inf"), save=True)
    int_cycles = result.homology.persistent_cycles.get(1, {})
    simplices = hc.log["simplicies"]
    int_edge = {}
    for cid, chain in int_cycles.items():
        for sid, coeff in chain.items():
            s = simplices[int(sid)]
            if len(s) == 2:
                int_edge.setdefault((min(s[0], s[1]), max(s[0], s[1])),
                                    set()).add(cid)
    harm_edge = {}
    for cyc in hc.log.get("harmonic_cycles", []):
        for e in cyc["edges"]:
            w = float(e["weight"])
            if abs(w) <= TOL:
                continue
            s = e["simplex"]
            if len(s) != 2:
                continue
            harm_edge.setdefault((min(int(s[0]), int(s[1])),
                                  max(int(s[0]), int(s[1]))), []).append(
                (cyc["cycle_index"], w))
    return int_edge, harm_edge


def rec(edge, int_edge, harm_edge):
    hs = harm_edge.get(edge, [])
    return (len(int_edge.get(edge, set())), len(hs),
            max((abs(w) for _, w in hs), default=0.0))


pairs = list(csv.DictReader(open(BASE / "interp_v2_pairs.csv")))
batch = int(sys.argv[1])
chunk = pairs[batch * 10:(batch + 1) * 10]

rows = []
for p in chunk:
    seed = p["type"] and p["seed"]
    M, node2idx, nt, keys = load_seed(seed)
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
    int_edge, harm_edge = run_full(
        Mm, str(BASE / ("_tmp_interp_%s_%d_%d.json" % (seed, i, j))))
    e_ij = (min(ti, tj), max(ti, tj))
    e_im = (min(ti, n), max(ti, n))
    e_mj = (min(n, tj), max(n, tj))
    r_ij = rec(e_ij, int_edge, harm_edge)
    r_im = rec(e_im, int_edge, harm_edge)
    r_mj = rec(e_mj, int_edge, harm_edge)
    ev_post = ""
    if p["type"] == "reticulate":
        evs = []
        with open(BASE / seed / "pair_registry_v2.csv") as f:
            for r in csv.DictReader(f):
                if (int(r["tip_i"]), int(r["tip_j"])) == (i, j):
                    evs = [int(x) for x in r["events"].split(";")]
        with open(BASE / seed / "net_retic_edges.csv") as f:
            r_edges = [r for r in csv.DictReader(f) if r["from"] != ""]
        parts = []
        for eid in evs:
            a = int(r_edges[eid - 1]["from"])
            b = int(r_edges[eid - 1]["to"])
            ek = (min(node2idx[a], node2idx[b]), max(node2idx[a],
                                                     node2idx[b]))
            ri = rec(ek, int_edge, harm_edge)
            parts.append("R%d:%d/%d" % (eid, ri[0], ri[1]))
        ev_post = ";".join(parts)
    rows.append([seed, p["type"], i, j, round(d_ij, 3),
                 r_ij[0], r_ij[1], round(r_ij[2], 4),
                 r_im[0], r_im[1], round(r_im[2], 4),
                 r_mj[0], r_mj[1], round(r_mj[2], 4), ev_post])
    print("%s %-10s (%d,%d) d=%.3f | (i,j): %d/%d/%.4f | (i,m): %d/%d/%.4f | "
          "(m,j): %d/%d/%.4f | %s"
          % (seed, p["type"], i, j, d_ij, r_ij[0], r_ij[1], r_ij[2],
             r_im[0], r_im[1], r_im[2], r_mj[0], r_mj[1], r_mj[2], ev_post))

out = BASE / ("interp_experiment_v2_batch%d.csv" % batch)
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["seed", "type", "tip_i", "tip_j", "d",
                "ij_int_cyc", "ij_harm_cyc", "ij_harm_maxw",
                "im_int_cyc", "im_harm_cyc", "im_harm_maxw",
                "mj_int_cyc", "mj_harm_cyc", "mj_harm_maxw",
                "event_retic_post"])
    w.writerows(rows)
print("saved", out)
