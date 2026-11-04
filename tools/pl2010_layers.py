#!/usr/bin/env python3
"""Отсчёты по 6 пластам методичек (ПЛ2010) — из векторов PDF (tools/pl2010.py)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__)); import pl2010 as P
CLAY = 35.1  # линия глин СП, мВ от левого края (среднее максимумов в глинах 945.6, 963.3, 968.4, 972.7, 978.0)
L = [  # (№, кровля, подошва, тип, вышележащий (a,b), нижележащий (a,b))
 (1, 942.2, 944.0, 'плот.', (941.2, 942.2), (944.4, 946.0)),
 (2, 944.0, 946.2, 'глина', (942.6, 943.8), (946.4, 948.6)),
 (3, 960.5, 962.8, 'кол.',  (957.2, 960.45), (963.0, 965.0)),
 (4, 965.4, 967.0, 'плот.', (963.0, 965.2), (967.2, 969.0)),
 (5, 969.4, 970.6, 'кол.',  (967.2, 969.2), (970.8, 973.0)),
 (6, 974.6, 977.7, 'кол.',  (971.0, 974.4), (977.8, 978.8)),
]
def vals(n, a, b): return [v for d, vs in P.profile(n, a, b, 0.02) for v in vs]
def hump(a, b, side):
    """СП: ближайший к пласту локальный максимум (к линии глин) в соседнем пласте → Eвм = CLAY − SP"""
    pr = [(d, max(vs)) for d, vs in P.profile('SP', a, b, 0.05)]
    if side == 'above': pr = pr[::-1]          # идём от границы пласта наружу
    best = pr[0]
    for i in range(1, len(pr) - 2):            # первый локальный максимум (к линии глин)
        if pr[i][1] >= pr[i - 1][1] and pr[i][1] >= max(p[1] for p in pr[i + 1:i + 9]) and pr[i][1] - min(p[1] for p in pr[i + 1:i + 9]) >= 1.0:
            best = pr[i]; break
    else:
        best = max(pr, key=lambda t: t[1])
    return CLAY - min(CLAY, best[1]), best
if __name__ == '__main__':
    for n, t, b, typ, up, dn in L:
        bk = vals('BK', t, b); im = vals('IM', t, b); sp = vals('SP', t, b); kv = vals('KV', t, b); mbk = vals('MBK', t, b)
        res = typ != 'глина'
        rk = max(bk) if res else min(bk); sk = min(im) if res else max(im)
        print(f"№{n} {t}–{b} ({typ}) h={b-t:.1f} м  dс={np.mean(kv):.0f} мм (мин {min(kv):.0f}, макс {max(kv):.0f})")
        print(f"   БК ρк={'max' if res else 'min'} {rk:.1f} (сред {np.mean(bk):.1f}); МБК max {max(mbk):.1f} сред {np.mean(mbk):.1f}")
        print(f"   ИМ σк={'min' if res else 'max'} {sk:.1f} мСим/м (сред {np.mean(im):.0f})")
        for nm, (a, c) in (('выше', up), ('ниже', dn)):
            print(f"   {nm} {a}–{c}: БК медиана {np.median(vals('BK', a, c)):.1f} [{min(vals('BK', a, c)):.1f}…{max(vals('BK', a, c)):.1f}]  "
                  f"ИМ медиана {np.median(vals('IM', a, c)):.0f} [{min(vals('IM', a, c)):.0f}…{max(vals('IM', a, c)):.0f}]")
        e1, h1 = hump(*up, 'above'); e2, h2 = hump(*dn, 'below')
        du = CLAY - min(sp)
        print(f"   СП: min {min(sp):.1f} → ΔUсп={du:.1f} мВ ({du/10*1.6:.2f} см);  Eвм'={e1:.1f} (SP {h1[1]:.1f} @ {h1[0]:.1f}), Eвм\"={e2:.1f} (SP {h2[1]:.1f} @ {h2[0]:.1f}) → Eвм={(e1+e2)/2:.1f} мВ")
