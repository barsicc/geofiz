#!/usr/bin/env python3
"""Полный расчёт таблиц БК, ИК (и входных данных СП) по 6 пластам ПЛ2010: python3 tools/pl2010_solve.py"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import pl2010_layers as LY, bk, ik
RR = 2.0
rows = []
for n, t, b, typ, up, dn in LY.L:
    v = lambda c, a, z: LY.vals(c, a, z)
    res = typ != 'глина'
    h = round(b - t, 1); dc = round(np.mean(v('KV', t, b))) / 1000
    rk = max(v('BK', t, b)) if res else min(v('BK', t, b))
    r1, r2 = np.median(v('BK', *up)), np.median(v('BK', *dn))
    rvm = (r1 + r2) / 2
    B = bk.solve(rk, rvm, RR, dc, h)
    sk = min(v('IM', t, b)) if res else max(v('IM', t, b))
    ext = max if res else min
    s1, s2 = ext(v('IM', *up)), ext(v('IM', *dn)); svm = (s1 + s2) / 2
    rho_k = ik.rho_from_sigma(sk); rho_vm = ik.rho_from_sigma(svm)
    rns, _ = ik.borehole(rho_k, RR, dc); rpi, o, f = ik.shoulder(rns, rho_vm, h)
    rows.append((n, t, b, typ, h, rk, r1, r2, rvm, dc, B, sk, s1, s2, svm, rho_k, rho_vm, rns, rpi))
print('БК (АБК-3): №  интервал  h  ρк  ρвм(\'/")  ρвм  dс  h/dс  ρк/ρр  ρвм/ρр  ρп/ρр  ρп  [номограмма]')
for r in rows:
    n, t, b, typ, h, rk, r1, r2, rvm, dc, B = r[:11]
    print(f"{n}. {t}–{b} {typ:5s} {h:.1f}  {rk:5.1f}  ({r1:.1f}/{r2:.1f}) {rvm:5.2f}  {dc:.3f}  {B['h_d']:.1f}  {B['rk_rr']:.1f}  {B['rvm_rr']:.2f}  {B['rp_rr']:.1f}  {B['rp']:.0f}  [{B['page']}]{' (!)' if B['warn'] else ''}")
print('\nИК (6Ф1): №  интервал  h  σк  σвм(\'/")  σвм  dс  ρк  ρвм  ρк^н.с  ρп^ИМ')
for r in rows:
    n, t, b, typ, h = r[:5]; dc = r[9]; sk, s1, s2, svm, rho_k, rho_vm, rns, rpi = r[11:]
    print(f"{n}. {t}–{b} {typ:5s} {h:.1f}  {sk:5.1f}  ({s1:.0f}/{s2:.0f}) {svm:5.1f}  {dc:.3f}  {rho_k:5.1f}  {rho_vm:4.1f}  {rns:5.1f}  {rpi:5.1f}")
