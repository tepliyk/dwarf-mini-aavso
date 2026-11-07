#!/usr/bin/env python3
"""
transform_ab.py — тест для Роя Акселсена (07.11.2026): сирий Bayer (A) проти інтерполяції (B).

A: зелений суперпіксель 960x540 (два G усереднені), синій піксель (1,1) без інтерполяції — як у aavso_phot.py.
B: білінійна дебаєризація каліброваного кадру в повні 1920x1080 площини G і B (як у типовому DSLR-конвеєрі),
   фотометрія на повній роздільності, апертура rap*2 px (та сама площа неба).
Однакові кадри, еталони, калібрування, відсів хмар/насичення. Порівнюються коефіцієнти, розкид залишків
і покадровий розкид кожної зорі. Нічого не пише в AAVSO і не чіпає out/tcoef_*.json.

Запуск:
  python -X utf8 transform_ab.py --session "<сесія>" --field NGC869 --cloud "<хмарний стек>" [--step 1 --maxf 60]
  python -X utf8 transform_ab.py --summary        # таблиця A/B по всіх out/tcoef_ab_*.json + перехресна перевірка
"""
import argparse, glob, json, os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import convolve
import astroalign as aa
import sep
from aavso_phot import ASTRO, HERE, SAT, R_AP, master_dark, master_flat, green, hms, recentroid, log
from transform_coef import robust_fit

OUT = os.path.join(HERE, 'out')
K3 = np.ones((3, 3), np.float32)


def demosaic_bilinear(cal):
    """RGGB -> повні площини G і B білінійною інтерполяцією (середнє відомих сусідів у вікні 3x3)."""
    H, W = cal.shape
    mg = np.zeros((H, W), np.float32); mg[0::2, 1::2] = 1; mg[1::2, 0::2] = 1
    mb = np.zeros((H, W), np.float32); mb[1::2, 1::2] = 1
    out = []
    for m in (mg, mb):
        num = convolve(cal * m, K3, mode='mirror'); den = convolve(m, K3, mode='mirror')
        interp = num / np.maximum(den, 1e-6)
        out.append(np.where(m > 0, cal, interp).astype(np.float32))
    return out


def measure(a):
    sdir = os.path.join(ASTRO, a.session)
    frames = sorted(f for f in glob.glob(os.path.join(sdir, '*.fits'))
                    if not os.path.basename(f).startswith('failed') and 'stacked' not in os.path.basename(f))
    frames = frames[::a.step][:a.maxf]
    h0 = fits.getheader(frames[0])
    exp, gain, temp = float(h0['EXPTIME']), int(h0['GAIN']), int(h0.get('DET-TEMP', 25))
    log('[ab] %s: %d кадрів (крок %d), %gs g%d' % (a.session, len(frames), a.step, exp, gain))
    dark, _ = master_dark(exp, gain, temp)
    flat = master_flat()

    cd = fits.getdata(a.cloud).astype(np.float32); ch = fits.getheader(a.cloud)
    cG = cd[1] if cd.ndim == 3 else cd
    hh, ww = cG.shape[0] // 2 * 2, cG.shape[1] // 2 * 2
    cG2 = np.ascontiguousarray(cG[:hh, :ww].reshape(hh // 2, 2, ww // 2, 2).mean(axis=(1, 3)))
    cG2 -= np.median(cG2)
    w = WCS(ch, naxis=2)
    offx, offy = (1920 - cG.shape[1]) / 2.0, (1080 - cG.shape[0]) / 2.0

    chart = json.load(open(os.path.join(HERE, 'comps_%s.json' % a.field), encoding='utf-8'))
    stars = []
    for s in chart['photometry']:
        bd = {b['band']: b['mag'] for b in s['bands']}
        if 'B' in bd and 'V' in bd and a.vmin <= bd['V'] <= a.vmax:
            stars.append({'label': s['label'], 'auid': s['auid'], 'V': bd['V'], 'B': bd['B'],
                          'ra': hms(s['ra'], 15), 'dec': hms(s['dec'], 1)})
    x, y = w.all_world2pix([s['ra'] for s in stars], [s['dec'] for s in stars], 0)
    hx, hy = (x - offx + 0.5) / 2 - 0.5, (y - offy + 0.5) / 2 - 0.5
    det = sep.extract(cG2, 8.0, err=sep.Background(cG2).globalrms)
    d = [(det['x'][np.argmin(np.hypot(det['x'] - u, det['y'] - v))] - u,
          det['y'][np.argmin(np.hypot(det['x'] - u, det['y'] - v))] - v)
         for u, v in zip(hx, hy) if np.min(np.hypot(det['x'] - u, det['y'] - v)) < 15]
    if len(d) >= 3:
        hx, hy = hx + np.median([q[0] for q in d]), hy + np.median([q[1] for q in d])
    H, W = cG2.shape
    keep = [(s, u, v) for s, u, v in zip(stars, hx, hy) if 15 <= u <= W - 15 and 15 <= v <= H - 15]
    stars = [k[0] for k in keep]
    src = np.array([[k[1], k[2]] for k in keep])
    log('  еталонів з B і V у кадрі: %d' % len(stars))

    GA, BA, GB, BB = [], [], [], []
    tf_prev = None
    for f in frames:
        raw = fits.getdata(f).astype(np.float32)
        cal = (raw - dark) / flat
        # ---- A: суперпіксель ----
        g = np.ascontiguousarray(green(cal)); b = np.ascontiguousarray(cal[1::2, 1::2])
        bg = sep.Background(g, bw=32, bh=32); bb = sep.Background(b, bw=32, bh=32)
        gs = np.ascontiguousarray(g - bg.back()); bs = np.ascontiguousarray(b - bb.back())
        try:
            tf, _ = aa.find_transform(gs, cG2, detection_sigma=5, max_control_points=50); tf_prev = tf
        except Exception:
            if tf_prev is None:
                continue
            tf = tf_prev
        pos = tf.inverse(src)
        xs, ys = zip(*[recentroid(gs, u, v) for u, v in pos]); xs, ys = np.array(xs), np.array(ys)
        fgA, _, _ = sep.sum_circle(gs, xs, ys, a.rap, err=bg.globalrms)
        fbA, _, _ = sep.sum_circle(bs, xs - 0.25, ys - 0.25, a.rap, err=bb.globalrms)
        sat = np.zeros(len(xs), bool)
        for i, (u, v) in enumerate(zip(xs, ys)):
            xi, yi = int(round(u)), int(round(v))
            sub = raw[2 * (yi - 2):2 * (yi + 3), 2 * (xi - 2):2 * (xi + 3)]
            sat[i] = sub.size == 0 or sub.max() > SAT
        # ---- B: білінійна дебаєризація, повна роздільність ----
        gF, bF = demosaic_bilinear(cal)
        bgF = sep.Background(gF, bw=64, bh=64); bbF = sep.Background(bF, bw=64, bh=64)
        gFs = np.ascontiguousarray(gF - bgF.back()); bFs = np.ascontiguousarray(bF - bbF.back())
        xF, yF = 2 * xs + 0.5, 2 * ys + 0.5
        xF, yF = zip(*[recentroid(gFs, u, v) for u, v in zip(xF, yF)]); xF, yF = np.array(xF), np.array(yF)
        fgB, _, _ = sep.sum_circle(gFs, xF, yF, 2 * a.rap, err=bgF.globalrms)
        fbB, _, _ = sep.sum_circle(bFs, xF, yF, 2 * a.rap, err=bbF.globalrms)
        bad = sat | (fgA <= 0) | (fbA <= 0) | (fgB <= 0) | (fbB <= 0)
        GA.append(np.where(bad, np.nan, fgA)); BA.append(np.where(bad, np.nan, fbA))
        GB.append(np.where(bad, np.nan, fgB)); BB.append(np.where(bad, np.nan, fbB))
    if len(GA) < 3:
        sys.exit('замало вирівняних кадрів (%d)' % len(GA))
    GA, BA, GB, BB = map(np.array, (GA, BA, GB, BB))
    rel = np.nanmedian(GA / np.nanmedian(GA, axis=0), axis=1)
    fok = rel > 0.6 * np.nanpercentile(rel, 90)
    V = np.array([s['V'] for s in stars]); Bc = np.array([s['B'] for s in stars])
    res = {'session': a.session, 'field': a.field, 'exp': exp, 'gain': gain, 'n_frames': int(fok.sum()),
           'rap_A': a.rap, 'rap_B': 2 * a.rap, 'modes': {}}
    for mode, Gm, Bm in (('A', GA, BA), ('B', GB, BB)):
        Gm, Bm = Gm[fok], Bm[fok]
        gm, bm = -2.5 * np.log10(Gm), -2.5 * np.log10(Bm)
        gm -= np.nanmedian(gm - V, axis=1, keepdims=True)
        bm -= np.nanmedian(bm - Bc, axis=1, keepdims=True)
        fin = (np.mean(np.isfinite(gm), axis=0) > 0.5) & (np.mean(np.isfinite(bm), axis=0) > 0.5)
        g = np.nanmedian(gm, axis=0); b = np.nanmedian(bm, axis=0)
        ng = np.maximum(np.sum(np.isfinite(gm), axis=0), 1); nb = np.maximum(np.sum(np.isfinite(bm), axis=0), 1)
        gfr = np.nanstd(gm, axis=0); bfr = np.nanstd(bm, axis=0)          # покадровий розкид зорі
        gsd = gfr / np.sqrt(ng); bsd = bfr / np.sqrt(nb)
        ok = fin & (gsd < 0.02) & (bsd < 0.04)
        rows = [{'label': s['label'], 'V': float(V[i]), 'B': float(Bc[i]), 'g': float(g[i]), 'b': float(b[i]),
                 'gsd': float(gsd[i]), 'bsd': float(bsd[i]), 'g_frame_sd': float(gfr[i]), 'b_frame_sd': float(bfr[i])}
                for i, s in enumerate(stars) if ok[i]]
        res['modes'][mode] = {'stars': rows}
        fit(res['modes'][mode])
        f = res['modes'][mode]['fit']
        log('  %s: зір %2d | Tv_bv %+.3f ± %.3f (σ %.3f) | Tbv %.3f ± %.3f (σ %.3f) | мед. покадр. sd g %.3f b %.3f' % (
            mode, len(rows), f['Tv_bv'], f['Tv_bv_err'], f['Tv_bv_sd'], f['Tbv'], f['Tbv_err'], f['Tbv_sd'],
            np.median([r['g_frame_sd'] for r in rows]), np.median([r['b_frame_sd'] for r in rows])))
    jf = os.path.join(OUT, 'tcoef_ab_%s_%s.json' % (a.field, a.session[-23:-13].replace('-', '')))
    json.dump(res, open(jf, 'w', encoding='utf-8'), indent=1)
    log('[ab] -> %s' % jf)


def fit(m):
    st = m['stars']
    if len(st) < 6:
        m['fit'] = None; return
    V = np.array([s['V'] for s in st]); B = np.array([s['B'] for s in st])
    g = np.array([s['g'] for s in st]); b = np.array([s['b'] for s in st]); BV = B - V
    s1, c1, e1, m1, d1 = robust_fit(BV, V - g)
    s2, c2, e2, m2, d2 = robust_fit(b - g, BV)
    m['fit'] = {'Tv_bv': s1, 'c1': c1, 'Tv_bv_err': e1, 'Tv_bv_sd': d1, 'n1': int(m1.sum()),
                'Tbv': s2, 'c2': c2, 'Tbv_err': e2, 'Tbv_sd': d2, 'n2': int(m2.sum())}


def apply_check(fit_from, st_to):
    """Коефіцієнти з одного поля на інше: кожна зоря по черзі ціль, решта ансамбль. Повертає sd(V_kat - V_tr), sd(V_kat - g)."""
    V = np.array([s['V'] for s in st_to]); B = np.array([s['B'] for s in st_to])
    g = np.array([s['g'] for s in st_to]); b = np.array([s['b'] for s in st_to])
    tv, tbv = fit_from['Tv_bv'], fit_from['Tbv']
    rt, ru = [], []
    for i in range(len(st_to)):
        o = np.arange(len(st_to)) != i
        zp = np.median(V[o] - g[o] - tv * (B[o] - V[o]))
        zc = np.median((B[o] - V[o]) - tbv * (b[o] - g[o]))
        bv_i = tbv * (b[i] - g[i]) + zc
        rt.append(V[i] - (g[i] + zp + tv * bv_i))
        ru.append(V[i] - (g[i] + np.median(V[o] - g[o])))
    return float(np.std(rt, ddof=1)), float(np.std(ru, ddof=1))


def summary():
    files = sorted(glob.glob(os.path.join(OUT, 'tcoef_ab_*.json')))
    if not files:
        sys.exit('нема out/tcoef_ab_*.json')
    R = {os.path.basename(f): json.load(open(f, encoding='utf-8')) for f in files}
    log('%-32s %s %5s %9s %9s %7s %7s %8s %8s' % ('файл', 'м', 'зір', 'Tv_bv', 'Tbv', 'σ(V-g)', 'σ(BV)', 'sd_g/кадр', 'sd_b/кадр'))
    for fn, r in R.items():
        for mode in ('A', 'B'):
            m = r['modes'][mode]; f = m['fit']
            if not f:
                log('%-32s %s  замало зір' % (fn, mode)); continue
            log('%-32s %s %5d %+.3f±%.3f %.2f±%.2f %7.3f %7.3f %8.3f %8.3f' % (
                fn, mode, len(m['stars']), f['Tv_bv'], f['Tv_bv_err'], f['Tbv'], f['Tbv_err'], f['Tv_bv_sd'], f['Tbv_sd'],
                np.median([s['g_frame_sd'] for s in m['stars']]), np.median([s['b_frame_sd'] for s in m['stars']])))
    # спільний фіт і перехресна перевірка між полями
    for mode in ('A', 'B'):
        pooled = []
        for r in R.values():
            st = r['modes'][mode]['stars']
            if len(st) < 6: continue
            V = np.array([s['V'] for s in st]); g = np.array([s['g'] for s in st]); B = np.array([s['B'] for s in st]); b = np.array([s['b'] for s in st])
            z1 = np.median(V - g); z3 = np.median(B - b)
            pooled += [dict(s, g=s['g'] + z1, b=s['b'] + z3) for s in st]
        P = {'stars': pooled}; fit(P); f = P['fit']
        log('\n[%s] спільний фіт %d зір: Tv_bv %+.3f ± %.3f (σ %.3f)  Tbv %.3f ± %.3f (σ %.3f)' % (
            mode, len(pooled), f['Tv_bv'], f['Tv_bv_err'], f['Tv_bv_sd'], f['Tbv'], f['Tbv_err'], f['Tbv_sd']))
        names = list(R.keys())
        for i in names:
            for j in names:
                if i == j or not R[i]['modes'][mode]['fit'] or len(R[j]['modes'][mode]['stars']) < 6: continue
                st, su = apply_check(R[i]['modes'][mode]['fit'], R[j]['modes'][mode]['stars'])
                log('[%s] %s -> %s: sd(V-Vtr) %.3f, без трансформації %.3f (n %d)' % (
                    mode, i[9:-5], j[9:-5], st, su, len(R[j]['modes'][mode]['stars'])))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--session'); ap.add_argument('--field'); ap.add_argument('--cloud')
    ap.add_argument('--vmin', type=float, default=8.5); ap.add_argument('--vmax', type=float, default=13.0)
    ap.add_argument('--rap', type=float, default=R_AP)
    ap.add_argument('--step', type=int, default=1); ap.add_argument('--maxf', type=int, default=60)
    ap.add_argument('--summary', action='store_true')
    a = ap.parse_args()
    summary() if a.summary else measure(a)
