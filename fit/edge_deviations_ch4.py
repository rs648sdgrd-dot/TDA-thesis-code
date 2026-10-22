# Per-edge deviations of Table 4.1 (Section 4.4).
import io, os, csv

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "networks")
SEEDS = ['s117', 's120', 's175', 's25', 's179', 's187']

def read_csv(p):
    with io.open(p, encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))

print('%-6s %4s %6s | %10s %10s %8s | %10s %10s %8s | %6s %6s | %6s %6s | %8s' % (
    'seed', 'nE', 'nHyb',
    'meanDevH', 'meanAbsH', 'shrDevH',
    'meanDevO', 'meanAbsO', 'shrDevO',
    'negH', 'posH', 'negO', 'posO', 'relDevH'))
for s in SEEDS:
    d0 = os.path.join(BASE, s)
    ret = read_csv(os.path.join(d0, 'net_retic_edges.csv'))
    ret = [r for r in ret if r.get('from', '') != '']
    hyb_nodes = set(r['to'] for r in ret)
    dev = read_csv(os.path.join(d0, 'edge_dev_xhat_x0.csv'))

    is_hyb = []
    for r in dev:
        to = r['to'].strip()

        try:
            to = str(int(float(to)))
        except ValueError:
            pass
        is_hyb.append(to in hyb_nodes)
    nE = len(dev)
    nHyb = sum(is_hyb)
    H = [float(r['dev']) for r, h in zip(dev, is_hyb) if h]
    O = [float(r['dev']) for r, h in zip(dev, is_hyb) if not h]
    import statistics as st
    mH, mO = st.mean(H), st.mean(O)
    aH, aO = st.mean([abs(v) for v in H]), st.mean([abs(v) for v in O])
    sH = sum(abs(v) for v in H) / sum(abs(v) for v in H + O)
    sO = 1.0 - sH
    negH = sum(1 for v in H if v < 0)
    posH = sum(1 for v in H if v > 0)
    negO = sum(1 for v in O if v < 0)
    posO = sum(1 for v in O if v > 0)
    x0H = [float(r['x0']) for r, h in zip(dev, is_hyb) if h]
    relH = 100.0 * st.mean([d / x for d, x in zip(H, x0H)])
    print('%-6s %4d %6d | %10.4f %10.4f %7.1f%% | %10.4f %10.4f %7.1f%% | %6d %6d | %6d %6d | %6.1f%%' % (
        s, nE, nHyb, mH, aH, 100 * sH, mO, aO, 100 * sO, negH, posH, negO, posO, relH))
    print('        ratio mean|dev| hyb/other = %.2f ; mean signed: hyb=%.4f other=%.4f' % (
        aH / aO, mH, mO))
