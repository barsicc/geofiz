#!/usr/bin/env python3
"""ρп по ИК (зонд 6Ф1) — по оцифрованным палеткам ПАЛ стр. 7, 8, 9 (= рис. 2, 3, 4 методички МетИК).

Алгоритм (методичка МетИК, п. 5–7):
  1) σк, σвм (мСим/м) → ρк, ρвм по кривой 6Ф1 «скин-эффект» (стр. 7);
  2) поправка за скважину (стр. 8): группа кривых с шифром ρр, шифр кривой — dс (м) → ρк^н.с;
  3) поправка за толщину/вмещающие (стр. 9): группа кривых по ρвм (0,6; 1,2; 6,0; 12,0 — ближайшая), шифр — hпл (м) → ρп^ИМ.

Пример (методичка МетИК, пласт 1): python3 tools/ik.py --sk 43 --svm 75 140 --rr 2 --d 0.287 --h 1.4
   → ρк ≈ 21 (в мет. 23), ρвм ≈ 7.7 (7.5), ρк^н.с ≈ 18–19 (20), ρп ≈ 20 (23)
Отдельно: --rho 23 (сразу ρк, без скин-эффекта).
Точность ±5–10 %. Поле ρр = 5 Ом·м не оцифровано (на палетке 6 линий на 5 подписей) — снимать вручную.
"""
import argparse, json, os, math
import numpy as np

D = json.load(open(os.path.join(os.path.dirname(__file__), 'pal', 'ik_6f1.json')))
SK = np.array(D['skin_6F1'])  # lg ρк, lg σ


def rho_from_sigma(s):
    a = SK[np.argsort(SK[:, 1])]
    return 10 ** float(np.interp(math.log10(s), a[:, 1], a[:, 0]))


def _family_value(curves, key, lx):
    """curves: [[param,[[lgx,lgy],...]],...] → interpolate in param"""
    ps, vs = [], []
    for p, tr in curves:
        a = np.array(tr)
        if a[0, 0] - 0.05 <= lx <= a[-1, 0] + 0.05:
            ps.append(p); vs.append(np.interp(lx, a[:, 0], a[:, 1]))
    if not ps:
        return None, True
    o = np.argsort(ps); ps = np.array(ps)[o]; vs = np.array(vs)[o]
    out = key < ps[0] or key > ps[-1]
    return 10 ** float(np.interp(key, ps, vs)), out


def borehole(rk, rr, d):
    fams = {float(k): v for k, v in D['borehole'].items()}
    ks = sorted(fams)
    lx = math.log10(rk)
    if rr in fams:
        return _family_value(fams[rr], d, lx)
    # интерполяция между полями по lg ρр
    lo = max([k for k in ks if k < rr], default=ks[0]); hi = min([k for k in ks if k > rr], default=ks[-1])
    v1, o1 = _family_value(fams[lo], d, lx); v2, o2 = _family_value(fams[hi], d, lx)
    if lo == hi or v1 is None or v2 is None:
        return (v1 or v2), True
    t = (math.log10(rr) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
    return 10 ** (math.log10(v1) + t * (math.log10(v2) - math.log10(v1))), o1 or o2 or rr > ks[-1]


def shoulder(rns, rvm, h, field=None):
    fields = {float(k): (k, v) for k, v in D['shoulder'].items()}
    if field is None:
        field = min(fields, key=lambda f: abs(math.log10(f) - math.log10(rvm)))
    k, curves = fields[float(field)]
    v, o = _family_value(curves, h, math.log10(rns))
    return v, o, k


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--sk', type=float, help='σк, мСим/м')
    ap.add_argument('--rho', type=float, help='ρк, Ом·м (если уже переведено)')
    ap.add_argument('--svm', type=float, nargs='+', help="σвм (или σвм' σвм\"), мСим/м")
    ap.add_argument('--rvm', type=float, help='ρвм, Ом·м (если уже переведено)')
    ap.add_argument('--rr', type=float, default=2.0)
    ap.add_argument('--d', type=float, required=True, help='dс, м')
    ap.add_argument('--h', type=float, required=True, help='hпл, м')
    ap.add_argument('--field', type=float, help='принудительно поле ρвм палетки стр. 9')
    a = ap.parse_args()
    rk = a.rho if a.rho else rho_from_sigma(a.sk)
    if a.rvm:
        rvm = a.rvm
    else:
        svm = sum(a.svm) / len(a.svm); rvm = rho_from_sigma(svm)
        print(f"σвм = {svm:.1f} мСим/м")
    print(f"ρк = {rk:.1f} Ом·м (σк={a.sk});  ρвм = {rvm:.1f} Ом·м   [скин-эффект 6Ф1, стр. 7]" if a.sk else f"ρк = {rk}; ρвм = {rvm:.1f}")
    rns, o1 = borehole(rk, a.rr, a.d)
    print(f"ρк^н.с ≈ {rns:.1f} Ом·м   [скважина, ρр={a.rr}, dс={a.d}; стр. 8]" + ('  (!) вне диапазона шифров' if o1 else ''))
    rp, o2, f = shoulder(rns, rvm, a.h, a.field)
    print(f"ρп^ИМ ≈ {rp:.1f} Ом·м   [вмещающие: поле ρвм={f}, h={a.h}; стр. 9]" + ('  (!) h вне шифров 0,75–2,0 — экстраполяция' if o2 else ''))
    for g in D['shoulder']:
        if g != f:
            v, _, _ = shoulder(rns, rvm, a.h, float(g))
            if v: print(f"   [контроль: поле ρвм={g}: {v:.1f}]")
