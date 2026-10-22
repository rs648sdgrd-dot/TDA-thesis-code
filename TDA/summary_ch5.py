# Summary (Tables 5.1-5.4).
import csv
import glob
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data" / "networks"
SEEDS = ["s117", "s120", "s175", "s25", "s179", "s187"]


def load_bars(seed):
    with open(BASE / seed / "eve_full_rip.json") as f:
        log = json.load(f)
    return log["harmonic_cycles"]


per = {}

all_pairs = []
for f2 in [BASE / "interp_experiment_v2_batch0.csv",
           BASE / "interp_experiment_v2_batch1.csv",
           BASE / "interp_experiment_v2_batch2.csv"]:
    all_pairs += list(csv.DictReader(open(f2)))
assert len(all_pairs) == 24
assert sum(1 for r in all_pairs if r["type"] == "reticulate") == 12
assert sum(1 for r in all_pairs if r["type"] == "control") == 12
keys = [(r["seed"], r["type"], r["tip_i"], r["tip_j"]) for r in all_pairs]
assert len(keys) == len(set(keys))

for seed in SEEDS:
    d = BASE / seed
    with open(d / "net_retic_edges.csv") as f:
        nr = sum(1 for r in csv.DictReader(f) if r["from"] != "")
    with open(d / "net_extant_tips.csv") as f:
        nt = sum(1 for _ in csv.DictReader(f))
    bars = load_bars(seed)

    ev = set()
    shared = defaultdict(set)
    wmax = 0.0
    for r in csv.DictReader(open(BASE / "event_cycle_direct.csv")):
        if r["seed"] == seed:
            ev.add(int(r["event_id"]))
            shared[int(r["event_id"])].add(r["cycle_index"])
            wmax = max(wmax, abs(float(r["harmonic_weight"])))
    n_shared = sum(1 for cs in shared.values() if len(cs) > 1)

    ev_gene = set()
    for r in csv.DictReader(open(BASE / "event_cycle_matches_genealogical.csv")):
        if r["seed"] == seed:
            for e in r["events"].split(";"):
                ev_gene.add(int(e))

    m2 = []
    for f2 in [BASE / "interp_experiment_v2_batch0.csv",
               BASE / "interp_experiment_v2_batch1.csv",
               BASE / "interp_experiment_v2_batch2.csv"]:
        for r in csv.DictReader(open(f2)):
            if r["seed"] == seed:
                m2.append(r)
    ret_resp = sum(1 for r in m2 if r["type"] == "reticulate"
                   and (int(r["im_int_cyc"]) > 0 or int(r["mj_int_cyc"]) > 0))
    ctrl_resp = sum(1 for r in m2 if r["type"] == "control"
                    and (int(r["im_int_cyc"]) > 0 or int(r["mj_int_cyc"]) > 0))

    m3 = [r for r in csv.DictReader(open(BASE / "mr3_linkage.csv"))
          if r["seed"] == seed and r["type"] == "reticulate"]
    n_resp = sum(1 for r in m3 if int(r["n_resp_cyc"]) > 0)
    n_link = sum(1 for r in m3 if int(r["n_link_hits"]) > 0)
    ba = [r["before_after"] for r in m3]
    per[seed] = dict(nr=nr, nt=nt, nbars=len(bars), ev=len(ev),
                     n_shared=n_shared, wmax=wmax, ev_gene=len(ev_gene),
                     ret_resp=ret_resp, ctrl_resp=ctrl_resp,
                     n_resp=n_resp, n_link=n_link, ba=ba)
    print("%-5s tips=%2d retic=%2d bars=%d | MR1 direct=%2d/%-2d (%.0f%%) "
          "shared=%d maxw=%.3f gene_tip=%d/%d | MR2 ret=%d/2 ctrl=%d/2 | "
          "MR3 resp=%d/2 link=%d/2"
          % (seed, nt, nr, len(bars), len(ev), nr, 100 * len(ev) / nr,
             n_shared, wmax, len(ev_gene), nr, ret_resp, ctrl_resp,
             n_resp, n_link))


tot = defaultdict(int)
for s in SEEDS:
    for k in ["nr", "nt", "nbars", "ev", "n_shared", "ev_gene",
              "ret_resp", "ctrl_resp", "n_resp", "n_link"]:
        tot[k] += per[s][k]
print("\n===== AGGREGATE (final 6) =====")
print("extant tips: %d  retic events: %d  persistent H1 classes: %d"
      % (tot["nt"], tot["nr"], tot["nbars"]))
print("MR1 direct: %d/%d = %.1f%% events in >=1 harmonic cycle"
      % (tot["ev"], tot["nr"], 100 * tot["ev"] / tot["nr"]))
print("MR1 shared events: %d" % tot["n_shared"])
print("genealogical tip-pair matching: %d/%d events"
      % (tot["ev_gene"], tot["nr"]))
print("MR2: reticulate responders %d/12 = %.0f%%; control FPs %d/12"
      % (tot["ret_resp"], 100 * tot["ret_resp"] / 12, tot["ctrl_resp"]))
print("MR3: conditional linkage %d/%d = %.0f%%; end-to-end %d/12 = %.0f%%"
      % (tot["n_link"], tot["n_resp"],
         100 * tot["n_link"] / max(tot["n_resp"], 1),
         tot["n_link"], 100 * tot["n_link"] / 12))
from scipy.stats import fisher_exact
print("Fisher MR2 (final 6):",
      fisher_exact([[tot["ret_resp"], 12 - tot["ret_resp"]],
                    [tot["ctrl_resp"], 12 - tot["ctrl_resp"]]]))
