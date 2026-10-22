# Candidate pair selection (Section 5.3).
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "networks"
SEEDS = ["s117", "s120", "s175", "s25", "s179", "s187"]

old_ret_ok = {
    "s175": [(27, 31), (27, 30)],
    "s25": [(7, 33), (7, 32)],
    "s187": [(35, 37), (33, 37)],
}
old_ctrl = {
    "s117": [(3, 14), (5, 6)],
    "s120": [(2, 8), (2, 9)],
    "s175": [(12, 24), (7, 21)],
    "s25": [(19, 30), (19, 31)],
    "s179": [(25, 28), (17, 33)],
    "s187": [(14, 20), (28, 36)],
}

event_wmax = {}
direct_events = {}
with open(BASE / "event_cycle_direct.csv") as f:
    for r in csv.DictReader(f):
        key = (r["seed"], int(r["event_id"]))
        direct_events.setdefault(r["seed"], set()).add(int(r["event_id"]))
        event_wmax[key] = max(event_wmax.get(key, 0.0),
                              abs(float(r["harmonic_weight"])))

final = {}
for seed in SEEDS:
    d = BASE / seed
    reg = {}
    with open(d / "pair_registry_v2.csv") as f:
        for r in csv.DictReader(f):
            reg[(int(r["tip_i"]), int(r["tip_j"]))] = (int(r["r"]),
                                                       r["events"])
    with open(d / "net_extant_tips.csv") as f:
        tips = [int(r["phy_node"]) for r in csv.DictReader(f)]
    with open(d / "dists2.csv", newline="") as f:
        rd = csv.reader(f)
        hdr = next(rd)[1:]
        M = [[float(x) for x in r[1:]] for r in rd]
    nt = {}
    with open(d / "net_name_table.csv") as f:
        for r in csv.DictReader(f):
            if r["species"] != "NA":
                nt[int(r["species"])] = int(r["phy_node"])
    node2idx = {nt[int(h)]: i for i, h in enumerate(hdr)}
    tips_involved = {t for p in reg for t in p}


    ret_pairs = [p for p in old_ret_ok.get(seed, []) if reg.get(p, (0,))[0] == 1]
    used = set(ret_pairs)
    cand = []
    for p, (r, evs) in reg.items():
        if r != 1 or p in used:
            continue
        ev = int(evs.split(";")[0])
        if ev in direct_events.get(seed, set()):
            cand.append((event_wmax.get((seed, ev), 0.0), p))
    cand.sort(reverse=True)
    for _, p in cand:
        if len(ret_pairs) >= 2:
            break
        ret_pairs.append(p)
    assert len(ret_pairs) == 2, (seed, ret_pairs)


    ret_tips = {t for p in ret_pairs for t in p}
    free = [t for t in tips if t not in tips_involved]
    ctrl = []
    for (i, j) in ret_pairs:
        dref = M[node2idx[i]][node2idx[j]]
        best = None
        for p in old_ctrl.get(seed, []):
            if p in ctrl or set(p) & ret_tips:
                continue
            if reg.get(p, (0,))[0] != 0:
                continue
            dd = M[node2idx[p[0]]][node2idx[p[1]]]
            rel = abs(dd - dref) / max(dref, 1e-9)
            if rel <= 0.20 and (best is None or rel < best[0]):
                best = (rel, p)
        if best is not None:
            ctrl.append(best[1])
        else:
            fb = None
            pools = [free, [t for t in tips]]
            for pool in pools:
                for ai, a in enumerate(pool):
                    for b in pool[ai + 1:]:
                        if {a, b} & ret_tips or (a, b) in ctrl:
                            continue
                        if (a, b) in reg:
                            continue
                        dd = M[node2idx[a]][node2idx[b]]
                        rel = abs(dd - dref) / max(dref, 1e-9)
                        if fb is None or rel < fb[0]:
                            fb = (rel, (a, b))
                if fb is not None:
                    break
            assert fb is not None, seed
            ctrl.append(fb[1])
    final[seed] = (ret_pairs, ctrl)
    print("%s: ret=%s ctrl=%s" % (seed, ret_pairs, ctrl))

with open(BASE / "interp_v2_pairs.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["seed", "type", "tip_i", "tip_j"])
    for seed in SEEDS:
        for p in final[seed][0]:
            w.writerow([seed, "reticulate", p[0], p[1]])
        for p in final[seed][1]:
            w.writerow([seed, "control", p[0], p[1]])
print("saved interp_v2_pairs.csv")
