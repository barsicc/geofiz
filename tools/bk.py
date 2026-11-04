#!/usr/bin/env python3
"""ρп по БК (зонд АБК-3) — по оцифрованным номограммам ПАЛ стр. 2–5 (h = 4d, 8d, 16d, 32d).

Алгоритм как в методичке МБК (п. 5–7):
  1) h/dс → выбрать номограмму (ближайшая из 4d/8d/16d/32d);
  2) левое поле: точка (ρвм/ρр; ρк/ρр) → «параллельно двум ближайшим кривым» на правый край поля;
  3) по горизонтали в правое поле до вертикали dс (см) → шифр кривой = ρп/ρр (интерполяция в логарифмах);
  4) ρп = (ρп/ρр)·ρр.

Пример (методичка МетБК, пласт 1):
  python3 tools/bk.py --rk 60 --rvm 8.2 --rr 2 --d 0.287 --h 1.4     → ρп/ρр ≈ 52 (в методичке 50), ρп ≈ 100
Можно задать ρвм как два отсчёта: --rvm 10 6.4 (среднее арифметическое).
Точность оцифровки: ±1–2 px ≈ ±1–2 % по оси; итог по шифру ≈ ±5–10 % (как при ручном снятии).
"""
import argparse, json, os, math
import numpy as np

D = json.load(open(os.path.join(os.path.dirname(__file__), 'pal', 'bk_abk3.json')))
PAGES = {'4d': 4, '8d': 8, '16d': 16, '32d': 32}


def edge(page, rvm_r, rk_r):
    """значение ρк/ρр, перенесённое параллельно кривым на правый край левого поля"""
    lx, ly = math.log10(rvm_r), math.log10(rk_r)
    vals, ends = [], []
    for t in D[page]['left']:
        a = np.array(t)
        vals.append(np.interp(lx, a[:, 0], a[:, 1])); ends.append(a[-1, 1])
    o = np.argsort(vals); vals = np.array(vals)[o]; ends = np.array(ends)[o]
    i = int(np.clip(np.searchsorted(vals, ly) - 1, 0, len(vals) - 2))
    t = (ly - vals[i]) / (vals[i + 1] - vals[i])
    return 10 ** (ends[i] + t * (ends[i + 1] - ends[i])), (lx < D[page]['left'][0][0][0] or ly < vals[0] or ly > vals[-1])


def cipher(page, d_cm, y):
    ld, ly = math.log10(d_cm), math.log10(y)
    pts = []
    for c, tr in D[page]['right']:
        a = np.array(tr); pts.append((np.interp(ld, a[:, 0], a[:, 1]), math.log10(c)))
    pts.sort()
    ys = np.array([p[0] for p in pts]); cs = np.array([p[1] for p in pts])
    out = ly < ys[0] or ly > ys[-1]
    return 10 ** float(np.interp(ly, ys, cs)), out


def solve(rk, rvm, rr, d_m, h=None, page=None):
    hd = h / d_m if h else None
    if page is None:
        page = min(PAGES, key=lambda p: abs(math.log(hd / PAGES[p])))
    e, o1 = edge(page, rvm / rr, rk / rr)
    c, o2 = cipher(page, d_m * 100, e)
    return dict(page=page, h_d=hd, rk_rr=rk / rr, rvm_rr=rvm / rr, edge=e, rp_rr=c, rp=c * rr, warn=o1 or o2)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--rk', type=float, required=True, help='ρк БК, Ом·м')
    ap.add_argument('--rvm', type=float, nargs='+', required=True, help='ρвм (или два отсчёта ρвм\' ρвм")')
    ap.add_argument('--rr', type=float, default=2.0, help='ρр (ρс), Ом·м')
    ap.add_argument('--d', type=float, required=True, help='dс, м')
    ap.add_argument('--h', type=float, help='hпл, м')
    ap.add_argument('--page', choices=list(PAGES))
    a = ap.parse_args()
    rvm = sum(a.rvm) / len(a.rvm)
    r = solve(a.rk, rvm, a.rr, a.d, a.h, a.page)
    print(f"ρвм = {rvm:.2f}; h/dс = {r['h_d']:.1f} → номограмма h={r['page']}" if r['h_d'] else f"номограмма {r['page']}")
    print(f"ρк/ρр = {r['rk_rr']:.1f}; ρвм/ρр = {r['rvm_rr']:.1f}; на краю левого поля ≈ {r['edge']:.1f}")
    print(f"ρп/ρр ≈ {r['rp_rr']:.1f} → ρп ≈ {r['rp']:.0f} Ом·м" + ("   (!) вне поля палетки — экстраполяция" if r['warn'] else ''))
    # соседние номограммы для контроля
    for p in PAGES:
        if p != r['page']:
            q = solve(a.rk, rvm, a.rr, a.d, a.h, p)
            print(f"   [контроль: h={p}: ρп/ρр ≈ {q['rp_rr']:.1f}]")
