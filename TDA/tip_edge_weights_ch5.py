# Tip-edge weights and correlations (Table 5.2).
import csv
import io
import json
import statistics
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[1] / "data" / "networks"
SEEDS = ["s117", "s120", "s175", "s25", "s179", "s187"]
TOL = 1e-6


def descendants(children, start):
    out = set()
    stack = [start]
    while stack:
        v = stack.pop()
        for ch in children.get(v, []):
            out.add(ch)
            stack.append(ch)
    return out


def spearman(x, y):
    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i])
        out = [0.0] * len(v)
        i = 0
        while i < len(s):
            j = i
            while j + 1 < len(s) and v[s[j + 1]] == v[s[i]]:
                j += 1
            for t in range(i, j + 1):
                out[s[t]] = (i + j) / 2.0 + 1
            i = j + 1
        return out
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) *
           sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den


def _betacf(a, b, x):
    MAXIT = 300
    EPS = 3e-14
    FPMIN = 1e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d
    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            break
    return h


def _betai(a, b, x):
    import math
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def spearman_p(rho, n):
    """Two-sided asymptotic p (Student-t approximation, scipy-compatible;
    deterministic and reproducible from the CSV)."""
    import math
    if n <= 2:
        return float("nan")
    t = rho * math.sqrt((n - 2) / max(1e-12, 1.0 - rho * rho))
    return _betai((n - 2) / 2.0, 0.5, (n - 2) / ((n - 2) + t * t))


rows = []
for seed in SEEDS:
    d = BASE / seed
    nt = {}
    for r in csv.DictReader(io.open(d / "net_name_table.csv")):
        if r["species"] != "NA":
            nt[int(r["species"])] = int(r["phy_node"])
    hdr = next(csv.reader(io.open(d / "dists2.csv")))
    idx2node = [nt[int(h)] for h in hdr[1:]]
    tips = {int(r["phy_node"])
            for r in csv.DictReader(io.open(d / "net_extant_tips.csv"))}

    children = defaultdict(list)
    for r in csv.DictReader(io.open(d / "net_tree_edges.csv")):
        children[int(r["from"])].append(int(r["to"]))

    simple = {}
    for r in csv.DictReader(io.open(d / "inher.csv")):
        a, b = int(r["from"]), int(r["to"])
        da = descendants(children, a) & tips
        db = descendants(children, b) & tips
        if len(da) == 1 and len(db) == 1:
            u, v = sorted((next(iter(da)), next(iter(db))))
            simple[(u, v)] = int(r["k"])

    reg = {}
    reg_ev = {}
    for r in csv.DictReader(io.open(d / "pair_registry_v2.csv")):
        u, v = sorted((int(r["tip_i"]), int(r["tip_j"])))
        reg[(u, v)] = int(r["r"])
        reg_ev[(u, v)] = r["events"]

    log = json.load(io.open(d / "eve_full_rip.json"))
    for cyc in log["harmonic_cycles"]:
        for e in cyc["edges"]:
            w = float(e["weight"])
            if abs(w) <= TOL:
                continue
            s = e["simplex"]
            if len(s) != 2:
                continue
            u, v = sorted((idx2node[int(s[0])], idx2node[int(s[1])]))
            if u not in tips or v not in tips:
                continue
            key = (u, v)
            rv = reg.get(key, 0)
            sid = simple.get(key, "")
            rows.append(dict(seed=seed, cycle=cyc["cycle_index"], u=u, v=v,
                             w=w, r=rv, events=reg_ev.get(key, ""),
                             retic_related=1 if rv >= 1 else 0,
                             simple_event_id=sid,
                             is_simple_event=1 if sid != "" else 0,
                             clean_simple=1 if (sid != "" and rv == 1)
                             else 0))

with open(BASE / "tip_edge_weights.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["seed", "cycle", "u", "v", "w", "r", "events",
                "retic_related", "simple_event_id", "is_simple_event",
                "clean_simple"])
    for r in rows:
        w.writerow([r["seed"], r["cycle"], r["u"], r["v"], r["w"], r["r"],
                    r["events"], r["retic_related"], r["simple_event_id"],
                    r["is_simple_event"], r["clean_simple"]])

# pair-level W_ij = max_C |w_ij(C)| (as in the thesis)
W = defaultdict(float)
rel = {}
sim = {}
clean = {}
for r in rows:
    k = (r["seed"], r["u"], r["v"])
    W[k] = max(W[k], abs(float(r["w"])))
    rel[k] = r["retic_related"]
    sim[k] = r["is_simple_event"]
    clean[k] = r["clean_simple"]

pairs = sorted(W)
wv = [W[k] for k in pairs]
relv = [rel[k] for k in pairs]
simv = [sim[k] for k in pairs]
cleanv = [clean[k] for k in pairs]

print("tip-tip edge records with |w|>1e-6: %d (unique pairs: %d)"
      % (len(rows), len(pairs)))
print("max |w| over all tip-tip edges: %.6f" % max(wv))
print("exact |w|=1: %d/%d" % (
    sum(1 for r in rows if abs(abs(float(r["w"])) - 1.0) < 1e-9), len(rows)))
print()
print("pair-level W_ij = max_C |w_ij(C)| "
      "(pooled across the six networks; descriptive only):")
for name, x in [("retic_related", relv),
                ("is_simple_event", simv),
                ("clean_simple", cleanv)]:
    a = [wv[i] for i in range(len(pairs)) if x[i] == 1]
    b = [wv[i] for i in range(len(pairs)) if x[i] == 0]
    med = lambda z: statistics.median(z) if z else float("nan")
    print("  %-16s n=%-4d median=%.5f max=%.5f | not: n=%-4d median=%.5f max=%.5f"
          % (name, len(a), med(a), max(a) if a else 0.0,
             len(b), med(b), max(b) if b else 0.0))
    rho = spearman(x, wv)
    p = spearman_p(rho, len(pairs))
    print("  %-16s Spearman rho=%.4f  p=%.4f (asymptotic t; n=%d pairs; "
          "pairs are nested in 6 networks)"
          % (name, rho, p, len(pairs)))
