#!/usr/bin/env python3
"""Расчёт строки таблицы СП (методичка МСП, п. 9–11).

  python3 tools/sp.py --du 29.6 --evm 4.4 0 --nu 0.51            → ΔEсп, αсп, ηгл (ΔEмакс = 65 мВ по умолчанию)
  python3 tools/sp.py --cm 4.8 --arrow 1.6 ...                    → ΔUсп из сантиметров: X = см·10/стрелка
  python3 tools/sp.py --kmax --k18 69 --t 30 --rf 1.5 --rv 0.19   → Kсп,t и ΔEсп,макс
Формулы:  ΔEсп = (ΔUсп − Eвм)/ν + Eвм;  αсп = ΔEсп/ΔEсп,макс;  ηгл = 1 − αсп;
          Kсп,t = Kсп,18·(t + 273)/291;  ΔEсп,макс = Kсп,t·lg(ρф/ρв)  (знак «−» опускаем: отклонение влево от линии глин).
ρф по ρр — палетка 3-I (ПАЛ стр. 25); ρв по C и t — палетка 1-I (ПАЛ стр. 24). Для варианта методички/пробника:
ρр = 2,0 Ом·м при 30 °C → ρф = 1,5; C = 30 кг/м³, 30 °C → ρв = 0,19 → ΔEмакс = 72·lg(1,5/0,19) = 64,6 ≈ 65 мВ.
"""
import argparse, math

ap = argparse.ArgumentParser()
ap.add_argument('--du', type=float, help='ΔUсп, мВ'); ap.add_argument('--cm', type=float, help='ΔUсп в см (линейкой)')
ap.add_argument('--arrow', type=float, default=1.6, help='длина масштабной стрелки 10 мВ, см')
ap.add_argument('--evm', type=float, nargs='+', default=[0.0], help="Eвм' Eвм\" (мВ) или одно Eвм")
ap.add_argument('--nu', type=float); ap.add_argument('--emax', type=float, default=65.0)
ap.add_argument('--kmax', action='store_true'); ap.add_argument('--k18', type=float, default=69); ap.add_argument('--t', type=float, default=30)
ap.add_argument('--rf', type=float, default=1.5); ap.add_argument('--rv', type=float, default=0.19)
a = ap.parse_args()
if a.kmax:
    k = a.k18 * (a.t + 273) / 291; e = k * math.log10(a.rf / a.rv)
    print(f"Kсп,t = {a.k18}·({a.t}+273)/291 = {k:.1f} мВ;  ΔEсп,макс = {k:.1f}·lg({a.rf}/{a.rv}) = {e:.1f} мВ")
du = a.du if a.du is not None else (a.cm * 10 / a.arrow if a.cm else None)
if du is not None and a.nu:
    evm = sum(a.evm) / len(a.evm)
    de = (du - evm) / a.nu + evm; al = de / a.emax
    print(f"ΔUсп = {du:.1f} мВ; Eвм = {evm:.1f} мВ; ν = {a.nu}")
    print(f"ΔEсп = ({du:.1f} − {evm:.1f})/{a.nu} + {evm:.1f} = {de:.1f} мВ;  αсп = {de:.1f}/{a.emax:g} = {al:.2f};  ηгл = 1 − {al:.2f} = {1 - al:.2f}")
