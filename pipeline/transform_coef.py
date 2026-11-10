#!/usr/bin/env python3
"""
transform_coef.py — коефіцієнти трансформації TG/TB -> Johnson V, B-V для DWARF Mini (23.10.2026).

Ідея (AAVSO DSLR Observing Manual, розд. трансформація):
  V - g = Tv_bv * (B - V) + c1        -> Tv_bv  (поправка зеленого каналу на колір)
  B - V = Tbv   * (b - g) + c2        -> Tbv    (колір зорі з інструментального b-g)
  B - b = Tb_bv * (B - V) + c3        -> Tb_bv  (контроль)
де g, b — інструментальні величини G- і B-каналу Байєра, B, V — каталог VSP еталонів поля.

Вхід — ті самі сирі сесії й хмарні стеки, що й aavso_phot.py; нічого не пише в AAVSO.
Запуск:
  python -X utf8 transform_coef.py --session "<сесія>" --field NGC869 --cloud "<хмарний стек>" [--step 3 --maxf 40]
  python -X utf8 transform_coef.py --summary     # зведення всіх out/tcoef_*.json
"""
import argparse, glob, json, os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astroalign as aa
import sep
from aavso_phot import ASTRO, HERE, SAT, R_AP, master_dark, master_flat, green, hms, recentroid, log

OUT = os.path.join(HERE, 'out')


def robust_fit(x, y, nsig=3.0, iters=5):
    """Лінійна регресія з ітеративним відсіком за MAD. Повертає slope, icpt, slope_err, mask, sd."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    for _ in range(iters):
        if m.sum() < 5:
            break
        A = np.vstack([x[m], np.ones(m.sum())]).T
        (s, c), *_ = np.linalg.lstsq(A, y[m], rcond=None)
        r = y - (s * x + c)
        sd = 1.4826 * np.median(np.abs(r[m] - np.median(r[m])))
        new = np.isfinite(r) & (np.abs(r) < nsig * max(sd, 0.01))
        if (new == m).all():
            break
        m = new
    A = np.vstack([x[m], np.ones(m.sum())]).T
    (s, c), res, *_ = np.linalg.lstsq(A, y[m], rcond=None)
    r = y[m] - (s * x[m] + c)
    sd = float(np.std(r, ddof=2)) if m.sum() > 2 else np.nan
    sxx = np.sum((x[m] - x[m].mean()) ** 2)
    se = sd / np.sqrt(sxx) if sxx > 0 else np.nan
    return float(s), float(c), float(se), m, sd


def measure(a):
    sdir = os.path.join(ASTRO, a.session)
    frames = sorted(f for f in glob.glob(os.path.join(sdir, '*.fits'))
                    if not os.path.basename(f).startswith('failed') and 'stacked' not in os.path.basename(f))
    frames = frames[::a.step][:a.maxf]
    h0 = fits.getheader(frames[0])
    if str(h0.get('BAYERPAT', 'RGGB')).strip().upper() != 'RGGB':
        sys.exit('BAYERPAT %s — очікувався RGGB' % h0.get('BAYERPAT'))
    exp, gain, temp = float(h0['EXPTIME']), int(h0['GAIN']), int(h0.get('DET-TEMP', 25))
    log('[tcoef] %s: %d кадрів (крок %d), %gs g%d' % (a.session, len(frames), a.step, exp, gain))
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
        mdx, mdy = np.median([q[0] for q in d]), np.median([q[1] for q in d])
        hx, hy = hx + mdx, hy + mdy
        log('  зсув WCS: (%.2f, %.2f) px за %d зорями' % (mdx, mdy, len(d)))
    H, W = cG2.shape
    keep = [(s, u, v) for s, u, v in zip(stars, hx, hy) if 15 <= u <= W - 15 and 15 <= v <= H - 15]
    stars = [k[0] for k in keep]
    src = np.array([[k[1], k[2]] for k in keep])
    log('  еталонів з B і V у кадрі: %d' % len(stars))

    G_all, B_all = [], []
    tf_prev = None
    for f in frames:
        raw = fits.getdata(f).astype(np.float32)
        cal = (raw - dark) / flat
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
        xs, ys = zip(*[recentroid(gs, u, v) for u, v in pos])
        xs, ys = np.array(xs), np.array(ys)
        # центр G-суперпікселя — сирі (0,5; 0,5), B — (1; 1): B-координата = G-координата − 0,25
        fg, _, _ = sep.sum_circle(gs, xs, ys, a.rap, err=bg.globalrms)
        fb, _, _ = sep.sum_circle(bs, xs - 0.25, ys - 0.25, a.rap, err=bb.globalrms)
        sat = np.zeros(len(xs), bool)
        for i, (u, v) in enumerate(zip(xs, ys)):
            xi, yi = int(round(u)), int(round(v))
            sub = raw[2 * (yi - 2):2 * (yi + 3), 2 * (xi - 2):2 * (xi + 3)]   # лише апертура, без сусідів
            sat[i] = sub.size == 0 or sub.max() > SAT
        fg = np.where(sat | (fg <= 0), np.nan, fg); fb = np.where(sat | (fb <= 0), np.nan, fb)
        G_all.append(fg); B_all.append(fb)
    if len(G_all) < 3:
        sys.exit('замало вирівняних кадрів (%d)' % len(G_all))
    G_all, B_all = np.array(G_all), np.array(B_all)
    # хмарні кадри геть: медіанний потік зір < 60 % від найкращих кадрів
    rel = np.nanmedian(G_all / np.nanmedian(G_all, axis=0), axis=1)
    fok = rel > 0.6 * np.nanpercentile(rel, 90)
    G_all, B_all = G_all[fok], B_all[fok]
    gm, bm = -2.5 * np.log10(G_all), -2.5 * np.log10(B_all)
    V = np.array([s['V'] for s in stars]); Bc = np.array([s['B'] for s in stars])
    # нормування прозорості кадру: медіана (g − V) по зорях
    gm -= np.nanmedian(gm - V, axis=1, keepdims=True)
    bm -= np.nanmedian(bm - Bc, axis=1, keepdims=True)
    fin = (np.mean(np.isfinite(gm), axis=0) > 0.5) & (np.mean(np.isfinite(bm), axis=0) > 0.5)
    g = np.nanmedian(gm, axis=0); b = np.nanmedian(bm, axis=0)
    ng = np.maximum(np.sum(np.isfinite(gm), axis=0), 1); nb = np.maximum(np.sum(np.isfinite(bm), axis=0), 1)
    gsd = np.nanstd(gm, axis=0) / np.sqrt(ng); bsd = np.nanstd(bm, axis=0) / np.sqrt(nb)   # похибка середнього
    log('  кадрів після відсіву хмар: %d/%d; зір з даними %d; мед. похибка g %.3f, b %.3f' % (
        fok.sum(), len(fok), fin.sum(), np.nanmedian(gsd[fin]) if fin.any() else np.nan, np.nanmedian(bsd[fin]) if fin.any() else np.nan))
    ok = fin & (gsd < 0.02) & (bsd < 0.04)     # лише зорі з точним середнім
    log('  зір із надійним g і b: %d (кадрів %d)' % (ok.sum(), len(G_all)))
    rows = [{'label': s['label'], 'V': float(V[i]), 'B': float(Bc[i]), 'g': float(g[i]), 'b': float(b[i]),
             'gsd': float(gsd[i]), 'bsd': float(bsd[i])} for i, s in enumerate(stars) if ok[i]]
    res = {'session': a.session, 'field': a.field, 'exp': exp, 'gain': gain, 'n_frames': len(G_all), 'stars': rows}
    fit_and_report(res)
    jf = os.path.join(OUT, 'tcoef_%s_%s.json' % (a.field, a.session[-23:-13].replace('-', '')))
    json.dump(res, open(jf, 'w', encoding='utf-8'), indent=1)
    log('[tcoef] -> %s' % jf)


def fit_and_report(res, plot=None):
    st = res['stars']
    if len(st) < 6:
        log('  ! замало зір для фіту (%d)' % len(st)); res['fit'] = None; return
    V = np.array([s['V'] for s in st]); B = np.array([s['B'] for s in st])
    g = np.array([s['g'] for s in st]); b = np.array([s['b'] for s in st])
    BV = B - V
    s1, c1, e1, m1, d1 = robust_fit(BV, V - g)
    s2, c2, e2, m2, d2 = robust_fit(b - g, BV)
    s3, c3, e3, m3, d3 = robust_fit(BV, B - b)
    res['fit'] = {'Tv_bv': s1, 'Tv_bv_err': e1, 'Tv_bv_sd': d1, 'n1': int(m1.sum()),
                  'Tbv': s2, 'Tbv_err': e2, 'Tbv_sd': d2, 'n2': int(m2.sum()),
                  'Tb_bv': s3, 'Tb_bv_err': e3, 'Tb_bv_sd': d3, 'n3': int(m3.sum()),
                  'BV_range': [float(np.min(BV)), float(np.max(BV))]}
    f = res['fit']
    log('  Tv_bv = %+.3f ± %.3f (σ %.3f, n %d) | Tbv = %.3f ± %.3f (σ %.3f, n %d) | Tb_bv = %+.3f ± %.3f | B-V %.2f..%.2f' % (
        f['Tv_bv'], f['Tv_bv_err'], f['Tv_bv_sd'], f['n1'], f['Tbv'], f['Tbv_err'], f['Tbv_sd'], f['n2'],
        f['Tb_bv'], f['Tb_bv_err'], *f['BV_range']))
    if plot:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
        for axx, xx, yy, mm, ss, cc, xl, yl in ((ax[0], BV, V - g, m1, s1, c1, 'B-V (VSP)', 'V - g'),
                                                (ax[1], b - g, BV, m2, s2, c2, 'b - g (інстр.)', 'B-V (VSP)')):
            axx.plot(xx[mm], yy[mm], 'o', ms=4); axx.plot(xx[~mm], yy[~mm], 'rx', ms=5)
            q = np.linspace(np.nanmin(xx), np.nanmax(xx), 10); axx.plot(q, ss * q + cc, 'k-')
            axx.set_xlabel(xl); axx.set_ylabel(yl); axx.set_title('нахил %.3f' % ss)
        fig.suptitle(plot[1]); fig.tight_layout(); fig.savefig(plot[0], dpi=110); plt.close(fig)


def summary():
    files = sorted(glob.glob(os.path.join(OUT, 'tcoef_*.json')))
    allst = []
    for fn in files:
        r = json.load(open(fn, encoding='utf-8'))
        log('%-34s' % os.path.basename(fn)); fit_and_report(r)
        allst += r['stars']
    if not allst:
        sys.exit('нема out/tcoef_*.json')
    # спільний фіт: нуль-пункти полів різні, тому віднімаємо медіану (V-g) і (B-V)-(b-g) кожного поля
    byf = {}
    for fn in files:
        r = json.load(open(fn, encoding='utf-8'))
        for s in r['stars']:
            byf.setdefault(fn, []).append(s)
    pooled = []
    for fn, st in byf.items():
        V = np.array([s['V'] for s in st]); g = np.array([s['g'] for s in st])
        B = np.array([s['B'] for s in st]); b = np.array([s['b'] for s in st])
        z1 = np.median((V - g)); z2 = np.median((B - V) - (b - g)); z3 = np.median(B - b)
        for s in st:
            pooled.append(dict(s, g=s['g'] + z1, b=s['b'] + z3))
    res = {'stars': pooled}
    log('\nСПІЛЬНИЙ ФІТ усіх сесій (%d зоре-вимірів):' % len(pooled))
    fit_and_report(res, plot=(os.path.join(OUT, 'tcoef_summary.png'), 'DWARF Mini KTFA: трансформація TG/TB'))
    json.dump(res['fit'], open(os.path.join(OUT, 'tcoef_summary.json'), 'w', encoding='utf-8'), indent=1)
    log('-> out/tcoef_summary.json, out/tcoef_summary.png')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--session'); ap.add_argument('--field'); ap.add_argument('--cloud')
    ap.add_argument('--vmin', type=float, default=8.5); ap.add_argument('--vmax', type=float, default=13.0)
    ap.add_argument('--rap', type=float, default=R_AP)
    ap.add_argument('--step', type=int, default=3); ap.add_argument('--maxf', type=int, default=40)
    ap.add_argument('--summary', action='store_true')
    a = ap.parse_args()
    summary() if a.summary else measure(a)
