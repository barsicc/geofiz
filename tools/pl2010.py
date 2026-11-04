#!/usr/bin/env python3
"""Точные отсчёты с планшета ПЛ2010 из векторов PDF (кривые — залитые контуры штрихов).

python3 tools/pl2010.py 942.2 944.0          → экстремумы БК/МБК/ИМ/СП/КВ в интервале + средние
python3 tools/pl2010.py --at 945.6           → значения всех кривых на глубине (по центру штриха)
Калибровка (pt PDF): глубина = 940 + (y − 213,9)/12,126; БК/МБК: lg ρ = (x − 222,4)/71,9;
ИМ: σ = (x − 400,2)/(578,8 − 400,2)·200; КВ: dс = 280 + (x − 400,2)/(178,6)·60;
СП: мВ от левого края = (x − 42,75)/(221,55 − 42,75)·40 (10 мВ на клетку 1,6 см; линия глин ≈ 34–35 мВ).
"""
import sys, math, argparse, warnings
warnings.filterwarnings('ignore')
import pymupdf as fitz
import numpy as np, os

PDF = os.path.join(os.path.dirname(__file__), '..', '02_Практика', 'ПЛ2010_Планшет_СП_КС_БК_МБК_ИК_КВ_940-980м.pdf')
IDS = {'SP': [44], 'BK': [48, 49], 'MBK': [50, 51], 'IM': [52], 'KV': [53]}
dep = lambda y: 940 + (y - 213.9) / 12.126
ydep = lambda d: 213.9 + (d - 940) * 12.126
CONV = {
    'BK': lambda x: 10 ** ((x - 222.4) / 71.9), 'MBK': lambda x: 10 ** ((x - 222.4) / 71.9),
    'IM': lambda x: (x - 400.2) / 178.6 * 200, 'KV': lambda x: 280 + (x - 400.2) / 178.6 * 60,
    'SP': lambda x: (x - 42.75) / (221.55 - 42.75) * 40,
}


def segments(ids):
    p = fitz.open(PDF)[0]; dr = p.get_drawings(); segs = []
    for i in ids:
        for it in dr[i]['items']:
            if it[0] == 'l': segs.append((it[1], it[2]))
            elif it[0] == 'c':
                pts = [it[1], it[2], it[3], it[4]]
                for a, b in zip(pts[:-1], pts[1:]): segs.append((a, b))
            elif it[0] == 're':
                r = it[1]; segs += [((r.x0, r.y0), (r.x1, r.y0)), ((r.x1, r.y0), (r.x1, r.y1))]
    return [((a[0], a[1]), (b[0], b[1])) for a, b in segs]


_cache = {}
def crossings(name, y):
    """x-координаты пересечения горизонтали y с контурами → центры штрихов (пары)"""
    if name not in _cache: _cache[name] = segments(IDS[name])
    xs = []
    for (x1, y1), (x2, y2) in _cache[name]:
        if (y1 - y) * (y2 - y) < 0 or (y1 == y and y2 != y):
            xs.append(x1 + (x2 - x1) * (y - y1) / (y2 - y1))
    xs.sort()
    # центр каждой пары (левый край – правый край штриха)
    return [(xs[i] + xs[i + 1]) / 2 for i in range(0, len(xs) - 1, 2)]


def profile(name, d0, d1, step=0.02):
    out = []
    for d in np.arange(d0, d1 + 1e-9, step):
        c = crossings(name, ydep(d))
        if c: out.append((d, [CONV[name](x) for x in c]))
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('d', type=float, nargs='*'); ap.add_argument('--at', type=float)
    a = ap.parse_args()
    if a.at:
        for n in IDS:
            print(n, [round(v, 2) for v in [CONV[n](x) for x in crossings(n, ydep(a.at))]])
        sys.exit()
    d0, d1 = a.d
    for n in ['BK', 'MBK', 'IM', 'SP', 'KV']:
        pr = profile(n, d0, d1)
        vals = [(d, v) for d, vs in pr for v in vs]
        if not vals: continue
        mx = max(vals, key=lambda t: t[1]); mn = min(vals, key=lambda t: t[1])
        mean = np.mean([np.mean(vs) for d, vs in pr])
        print(f"{n:4s} max {mx[1]:7.2f} @ {mx[0]:.2f}   min {mn[1]:7.2f} @ {mn[0]:.2f}   mean {mean:7.2f}")
