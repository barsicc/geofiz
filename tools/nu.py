#!/usr/bin/env python3
"""νсп по палеткам ПАЛ стр. 11–21 (картинки: рисунки/ПАЛ/nu/стрNN{в,н}.png).

  python3 tools/nu.py 15в 8.0            → пересечения всех кривых с вертикалью H/dс = 8 (ν сверху вниз)
  python3 tools/nu.py 15в 8.0 --c 8      → + ν для шифра ρп/ρр = 8 (интерполяция между соседними шифрами в lg)
  python3 tools/nu.py 11в 50 --c 4.9     → палетка без проникновения: X = ρп/ρр, шифр = H/dс
  python3 tools/nu.py --list             → какие палетки откалиброваны, их параметры (ρвм/ρр; D/dс; ρзп/ρр)
Всегда смотреть контрольную картинку (печатается путь): вертикаль X, отмеченные пересечения и их номера.
Порядок пересечений сверху вниз = шифры по возрастанию ρп/ρр (стр. 12–21) / по убыванию H/dс (стр. 11),
НО кривые, которые на этом X ещё не начались (начинаются при ν≈0,2–0,27) или уже упёрлись в ν=1, выпадают —
поэтому соответствие «пересечение ↔ шифр» подтверждать по картинке (параметр --skip N сдвигает нумерацию).
Выбор палетки — ближайшие (линейно, как в методичке МетСП) ρвм/ρр ∈ {1; 5}, D/dс ∈ {2; 4; 6}, ρзп/ρр ∈ {10; 20; 50; 100}.
"""
import sys, os, math, argparse
import numpy as np
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, '..', 'рисунки', 'ПАЛ', 'nu')
# x0 — пиксель X = xv0; ppd — пикселей на декаду; y1/y0 — пиксели ν = 1 и ν = 0; ys — диапазон строк (для стр. 11 два поля)
CAL = {  # k = пиксели углов (X_min, ν=0), (X_max, ν=0), (X_min, ν=1), (X_max, ν=1); X = (X_min, X_max)
    '11в': dict(f='стр11_нет_проникн.png', kind='noinv', par='ρвм/ρр=5; D/dс=1 (верхняя; в мет. МетСП — 4-V)', X=(1, 150), k=[(53, 689), (1309, 698), (51, 85), (1307, 94)],
                c=[30, 20, 15, 10, 8, 6, 4, 2]),
    '11н': dict(f='стр11_нет_проникн.png', kind='noinv', par='ρвм/ρр=1; D/dс=1 (нижняя)', X=(1, 150), k=[(65, 1399), (1312, 1405), (66, 806), (1314, 815)],
                c=[30, 20, 15, 10, 8, 6, 4, 2]),
    '15в': dict(f='стр15в.png', kind='inv', par='ρвм/ρр=1; D/dс=6; ρзп/ρр=10', X=(4, 30), k=[(67, 985), (859, 985), (67, 77), (859, 77)],
                c=[0.1, 1, 5, 10, 20, 40, 60, 80, 100, 150]),
    '15н': dict(f='стр15н.png', kind='inv', par='ρвм/ρр=1; D/dс=6; ρзп/ρр=20', X=(4, 30), k=[(125, 995), (900, 970), (108, 90), (880, 70)],
                c=[0.1, 1, 10, 20, 40, 60, 100, 150]),
}
PARAMS = {  # стр.: (верхняя), (нижняя) — ρвм/ρр; D/dс; ρзп/ρр
    12: ((1, 2, 10), (1, 2, 20)), 13: ((1, 2, 50), (1, 4, 10)), 14: ((1, 4, 20), (1, 4, 50)), 15: ((1, 6, 10), (1, 6, 20)),
    16: ((1, 6, 50), (1, 6, 100)), 17: ((5, 2, 10), (5, 2, 20)), 18: ((5, 2, 100), (5, 4, 10)), 19: ((5, 4, 20), (5, 4, 100)),
    20: ((5, 6, 10), (5, 6, 20)), 21: ((5, 6, 50), (5, 6, 100)),
}


# остальные палетки: диапазон X и шифры (по картинкам); калибровка k — по запросу (--ruler ID, снять 4 угла)
for _id, _X, _c in [
    ('12в', (2, 30), [0.1, 5, 10, 20, 40, 60, 100, 150]), ('12н', (2, 30), [0.1, 5, 10, 20, 40, 60, 100, 150]),
    ('13в', (4, 30), [0.1, 1, 5, 10, 20, 40, 60, 100, 150]), ('13н', (4, 30), [0.1, 5, 10, 20, 40, 60, 100, 150]),
    ('14в', (4, 30), [0.1, 5, 10, 20, 40, 60, 100, 150]), ('14н', (6, 30), [0.1, 1, 10, 20, 40, 60, 100, 150]),
    ('16в', (6, 30), [1, 10, 20, 40, 60, 100, 150]), ('16н', (8, 30), [1, 10, 20, 40, 60, 100, 150]),
    ('17в', (2, 30), [0.1, 5, 10, 20, 40, 60, 100, 150]), ('17н', (2, 30), [0.1, 5, 10, 20, 40, 60, 100, 150]),
    ('18в', (4, 30), [0.1, 1, 5, 10, 20, 40, 60, 100, 150]), ('18н', (2, 30), [0.1, 5, 10, 20, 40, 60, 80, 100, 150]),
    ('19в', (4, 30), [0.1, 5, 10, 20, 40, 60, 80, 100, 150]), ('19н', (6, 30), [0.1, 1, 10, 20, 40, 60, 100, 150]),
    ('20в', (2.5, 30), [0.1, 5, 10, 20, 40, 60, 100, 150]), ('20н', (4, 30), [0.1, 1, 10, 20, 40, 60, 100, 150]),
    ('21в', (6, 30), [0.1, 1, 10, 20, 40, 60, 100, 150]), ('21н', (6, 30), [1, 10, 20, 40, 60, 100, 150])]:
    _p = PARAMS[int(_id[:2])][0 if _id[2] == 'в' else 1]
    CAL[_id] = dict(f=f'стр{_id}.png', kind='inv', par=f'ρвм/ρр={_p[0]}; D/dс={_p[1]}; ρзп/ρр={_p[2]}', X=_X, k=None, c=_c)
# ⚠ шифры/диапазоны 18в, 20в и нижние кривые 12в — проверить по картинке перед использованием


def ruler(id_):
    im = Image.open(os.path.join(IMG, CAL[id_]['f'])).convert('RGB'); W, H = im.size
    big = Image.new('RGB', (W + 70, H + 40), 'white'); big.paste(im, (60, 0)); d = ImageDraw.Draw(big)
    for y in range(0, H, 25):
        d.line((48 if y % 50 else 40, y, 60, y), fill=(255, 0, 0))
        if y % 50 == 0: d.text((0, y - 5), str(y), fill=(255, 0, 0))
    for x in range(0, W, 25):
        d.line((x + 60, H, x + 60, H + (10 if x % 50 == 0 else 5)), fill=(255, 0, 0))
        if x % 50 == 0: d.text((x + 52, H + 14), str(x), fill=(255, 0, 0))
    out = f'/tmp/ruler_{id_}.png'; big.save(out); return out


def choose(rvm_r, D_d, rzp_r):
    """как в методичке МетСП: каждый параметр — к ближайшему (линейно) значению; при равенстве — оба варианта"""
    def near(v, opts):
        d = sorted(opts, key=lambda o: abs(o - v)); return [o for o in d if abs(abs(o - v) - abs(d[0] - v)) < 1e-9]
    out = []
    for r in near(rvm_r, [1, 5]):
        for D in near(D_d, [2, 4, 6]):
            for z in near(rzp_r, [10, 20, 50, 100]):
                for p, (a, b) in PARAMS.items():
                    for tag, q in (('в', a), ('н', b)):
                        if q == (r, D, z): out.append((f'{p}{tag}', q))
    return out


def _pt(cal, lx, nu):
    """билинейное отображение (lg X, ν) → пиксель по 4 углам: k = [(x,y) для (Xmin,0),(Xmax,0),(Xmin,1),(Xmax,1)]"""
    (a, b, c, d) = cal['k']; u = (lx - math.log10(cal['X'][0])) / (math.log10(cal['X'][1]) - math.log10(cal['X'][0]))
    p0 = (a[0] + (b[0] - a[0]) * u, a[1] + (b[1] - a[1]) * u); p1 = (c[0] + (d[0] - c[0]) * u, c[1] + (d[1] - c[1]) * u)
    return p0[0] + (p1[0] - p0[0]) * nu, p0[1] + (p1[1] - p0[1]) * nu


def _line_cross(dark, cal, lx):
    n = 1400; ts = np.linspace(0.005, 0.995, n)
    pts = [_pt(cal, lx, t) for t in ts]
    prof = np.array([dark[int(round(y)), int(round(x))] for x, y in pts])
    runs, cur = [], None
    for i, v in enumerate(prof):
        if v and cur is None: cur = i
        if (not v) and cur is not None: runs.append((cur, i - 1)); cur = None
    if cur is not None: runs.append((cur, n - 1))
    pxl = abs(pts[-1][1] - pts[0][1]) / n  # пикселей на шаг
    out = []
    for a, b in runs:
        w = (b - a + 1) * pxl
        if 2 <= w <= 16:
            t = ts[(a + b) // 2]; x, y = _pt(cal, lx, t)
            if w <= 6 and abs(t * 10 - round(t * 10)) * 0.1 < 0.0045: continue   # горизонтальная линия сетки (ν кратно 0,1)
            out.append((t, (x, y)))
    return out, prof.mean()


def crossings(cal, X, ys=None, half=6, minhits=4, tol=0.012):
    """пересечения кривых с линией X = const: 13 соседних линий (±half px), кластеризация по ν —
    устойчиво к линиям сетки и прерывистым тонким кривым на сканах"""
    dark = np.array(Image.open(os.path.join(IMG, cal['f'])).convert('L')) < 140
    lx = math.log10(X); step = (math.log10(cal['X'][1]) - math.log10(cal['X'][0])) / abs(cal['k'][1][0] - cal['k'][0][0])
    hits = []
    for dx in range(-half, half + 1):
        c, cov = _line_cross(dark, cal, lx + dx * step)
        if cov > 0.25: continue                     # сидим на вертикальной линии сетки
        hits += [(t, dx, p) for t, p in c]
    hits.sort()
    cl = []
    for h in hits:
        if cl and h[0] - cl[-1][-1][0] <= tol: cl[-1].append(h)
        else: cl.append([h])
    out = []
    for g in cl:
        if len({dx for _, dx, _ in g}) < minhits: continue
        # линейная регрессия ν(dx) → значение в dx = 0 (кривые наклонные)
        t = np.array([h[0] for h in g]); d = np.array([h[1] for h in g])
        A = np.vstack([d, np.ones_like(d)]).T; k, b = np.linalg.lstsq(A, t, rcond=None)[0]
        out.append((float(b), _pt(cal, lx, float(b))))
    return out


def draw(cal, X, nus, out):
    im = Image.open(os.path.join(IMG, cal['f'])).convert('RGB'); d = ImageDraw.Draw(im)
    lx = math.log10(X); d.line((*_pt(cal, lx, 0), *_pt(cal, lx, 1)), fill=(255, 0, 0), width=1)
    for v in np.arange(0, 1.01, 0.1):  # контроль калибровки: метки ν на краях
        for e in (cal['X'][0], cal['X'][1]):
            x, y = _pt(cal, math.log10(e), v); d.line((x - 8, y, x + 8, y), fill=(0, 160, 0), width=2)
    for i, (nu, (x, y)) in enumerate(nus):
        d.ellipse((x - 5, y - 5, x + 5, y + 5), outline=(0, 0, 255), width=2); d.text((x + 8, y - 6), f'{i}:{nu:.2f}', fill=(0, 0, 255))
    im.save(out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('id', nargs='?'); ap.add_argument('X', type=float, nargs='?')
    ap.add_argument('--c', type=float, help='шифр, для которого нужно ν'); ap.add_argument('--skip', type=int, default=0)
    ap.add_argument('--list', action='store_true'); ap.add_argument('--ruler', help='картинка с линейкой (пиксели) для калибровки'); ap.add_argument('--choose', type=float, nargs=3, metavar=('RVM_R', 'D_D', 'RZP_R'))
    a = ap.parse_args()
    if a.ruler:
        print(ruler(a.ruler), '— снять пиксели углов (Xmin,0),(Xmax,0),(Xmin,1),(Xmax,1) и вписать в CAL[...]["k"] (координата x = отсчёт линейки, без сдвига 60)'); sys.exit()
    if a.list:
        for k, v in CAL.items(): print(k, v['par'], 'калибровка есть' if v.get('k') else 'НЕТ калибровки (см. tools/nu.py)')
        sys.exit()
    if a.choose:
        print('ближайшая палетка:', choose(*a.choose)); sys.exit()
    cal = CAL[a.id]
    if not cal.get('k'): sys.exit('палетка не откалибрована: задать k (4 угла) в CAL (по линейке, см. регламент)')
    nus = crossings(cal, a.X); nus = sorted(nus, key=lambda t: -t[0])
    out = f'/tmp/nu_{a.id}_{a.X}.png'; draw(cal, a.X, nus, out)
    print(f"{a.id} [{cal['par']}]  X={a.X}; пересечения сверху вниз:")
    cs = cal['c'][a.skip:]
    for i, (nu, y) in enumerate(nus):
        print(f"  {i}: ν={nu:.3f}   (шифр по порядку: {cs[i] if i < len(cs) else '?'})")
    print('картинка:', out)
    if a.c is not None:
        pairs = [(c, nu) for c, (nu, y) in zip(cs, nus)]
        f = (lambda v: math.log10(v)) if cal['kind'] == 'inv' else (lambda v: math.log10(v))
        pairs.sort(key=lambda t: f(t[0]))
        xs = [f(c) for c, _ in pairs]; vs = [nu for _, nu in pairs]
        print(f"ν(шифр {a.c}) ≈ {np.interp(f(a.c), xs, vs):.2f}   (по парам {pairs})")
