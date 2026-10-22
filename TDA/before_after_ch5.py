# Before/after weights and Fisher test (Section 5.4).
import csv
import io
from pathlib import Path

BASE = Path(__file__).resolve().parents[1] / "data" / "networks"
SEEDS = {"s117", "s120", "s175", "s25", "s179", "s187"}

  # before: pair-level max |w|
before = {}
for r in csv.DictReader(io.open(BASE / "tip_edge_weights.csv")):
    key = (r["seed"], int(r["u"]), int(r["v"]))
    w = abs(float(r["w"]))
    if key not in before or w > before[key]:
        before[key] = w

  # after: interpolation batch fields (batches 0..2)
after = {}
for b in range(3):
    for r in csv.DictReader(io.open(
            BASE / ("interp_experiment_v2_batch%d.csv" % b))):
        if r["seed"] not in SEEDS:
            continue
        key = (r["seed"], int(r["tip_i"]), int(r["tip_j"]))
        after[key] = dict(
            ij_harm_maxw=float(r["ij_harm_maxw"]),
            im_harm_maxw=float(r["im_harm_maxw"]),
            mj_harm_maxw=float(r["mj_harm_maxw"]),
            im_int=int(r["im_int_cyc"]), mj_int=int(r["mj_int_cyc"]),
            ij_int=int(r["ij_int_cyc"]))


path = {}
for r in csv.DictReader(io.open(BASE / "shortest_path_results.csv")):
    path[(r["seed"], int(r["tip_i"]), int(r["tip_j"]))] = r


reg = {}
for seed in sorted(SEEDS):
    for r in csv.DictReader(io.open(BASE / seed / "pair_registry_v2.csv")):
        reg[(seed, int(r["tip_i"]), int(r["tip_j"]))] = (int(r["r"]),
                                                         r["events"])

pairs = list(csv.DictReader(io.open(BASE / "interp_v2_pairs.csv")))
rows = []
for p in pairs:
    key = (p["seed"], int(p["tip_i"]), int(p["tip_j"]))
    rv, evs = reg.get(key, (0, ""))
    a = after.get(key)
    pr = path.get(key, {})
    if a is None:
        resp, h1 = "", ""
    else:
        resp = 1 if (a["im_int"] >= 1 or a["mj_int"] >= 1) else 0
        h1 = 1 if (abs(a["im_harm_maxw"] - 1.0) < 1e-9 or
                   abs(a["mj_harm_maxw"] - 1.0) < 1e-9) \
            else 0
    rows.append(dict(
        seed=p["seed"], type=p["type"], tip_i=int(p["tip_i"]),
        tip_j=int(p["tip_j"]), r=rv, events=evs,
        before_ij_max_abs_w=round(before.get(key, 0.0), 6),
        before_ij_supported=1 if key in before else 0,
        after_ij_harm_maxw=a["ij_harm_maxw"] if a else "",
        after_im_harm_maxw=a["im_harm_maxw"] if a else "",
        after_mj_harm_maxw=a["mj_harm_maxw"] if a else "",
        after_im_int=a["im_int"] if a else "",
        after_mj_int=a["mj_int"] if a else "",
        responder_int=resp,
        harmonic_exact1=h1,
        d_tree=pr.get("d_tree", ""), d_net=pr.get("d_net", ""),
        shortening=pr.get("shortening", ""),
        target_retic_on_path=pr.get("target_retic_on_path", ""),
        n_retic_on_path=pr.get("n_retic_on_path", "")))

cols = ["seed", "type", "tip_i", "tip_j", "r", "events",
        "before_ij_max_abs_w", "before_ij_supported",
        "after_ij_harm_maxw", "after_im_harm_maxw", "after_mj_harm_maxw",
        "after_im_int", "after_mj_int", "responder_int", "harmonic_exact1",
        "d_tree", "d_net", "shortening",
        "target_retic_on_path", "n_retic_on_path"]
with open(BASE / "tip_pair_before_after.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(cols)
    for r in rows:
        w.writerow([r[c] for c in cols])

ret = [r for r in rows if r["type"] == "reticulate"]
con = [r for r in rows if r["type"] == "control"]
n_resp = sum(r["responder_int"] == 1 for r in ret)
n_h1 = sum(r["harmonic_exact1"] == 1 for r in ret)
n_con_resp = sum(r["responder_int"] == 1 for r in con)
n_con_h1 = sum(r["harmonic_exact1"] == 1 for r in con)
print("rows: %d" % len(rows))
print()
print("--- interpolation (frozen 12+12) ---")
print("integer responders: reticulate %d/12, controls %d/12"
      % (n_resp, n_con_resp))
print("harmonic exact |w|=1 on an m-edge: reticulate %d/12, controls %d/12"
      % (n_h1, n_con_h1))
print()
print("--- shortest path (selected-pair validation) ---")
print("any retic edge on path: reticulate %d/12 vs controls %d/12"
      % (sum(int(r["n_retic_on_path"]) > 0 for r in ret),
         sum(int(r["n_retic_on_path"]) > 0 for r in con)))
print("target retic edge on path (r=1 pairs): %d/12"
      % sum(r["target_retic_on_path"] == "yes" for r in ret))
print()
print("--- candidate-conditioned path hits ---")
print("Student-literal harmonic ±1 candidates: %d/12 -> target on path %d/%d"
      % (n_h1, sum(r["target_retic_on_path"] == "yes"
                   for r in ret if r["harmonic_exact1"] == 1), n_h1))
print("integer responders: %d/12 -> target on path %d/%d"
      % (n_resp, sum(r["target_retic_on_path"] == "yes"
                     for r in ret if r["responder_int"] == 1), n_resp))
print()
print("before: pairs with any (i,j) support %d/12, max |w| %.4f"
      % (sum(r["before_ij_supported"] for r in ret),
         max(r["before_ij_max_abs_w"] for r in ret)))
