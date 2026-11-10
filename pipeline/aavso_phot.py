#!/usr/bin/env python3
"""
aavso_phot.py — диференційна фотометрія змінних в архівних сесіях DWARF Mini (беклог A2/A3, 22.10.2026).

Сирий кадр (RGGB, 12 біт) -> дарк (exp/gain/найближча T) і флет (по каналах Байєра) ->
зелений суперпіксель (G1+G2)/2, 960x540, ~8,03"/px -> вирівнювання до WCS хмарного стека
(astroalign, хмарний G біновано 2x2) -> апертура sep r=2,5 px, кільце 5–8 ->
ансамбль еталонів VSP (V), контрольна зоря -> бінінг часового ряду -> JSON + PNG.

Запуск (ПК):
  python -X utf8 aavso_phot.py --session "DWARF_RAW_TELE_M 33_EXP_30_GAIN_60_2026-08-26-00-13-50-252" ^
      --field M33 --cloud "C:\\astro\\Fits хмара\\a01f82e099a373b6e4f5a72090833be3.fits" --stars "V Tri" --bin 180
"""
import argparse, glob, json, os, re, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.utils import iers
iers.conf.auto_download = False      # без мережі: старі IERS-таблиці дають помилку для дат після їхнього кінця
iers.conf.auto_max_age = None
import astroalign as aa
import sep

ASTRO = r'C:\astro\DWARF Mini\Astronomy'
HERE = os.path.dirname(os.path.abspath(__file__))
SITE = EarthLocation(lat=48.45 * u.deg, lon=35.0 * u.deg, height=100 * u.m)
TZ = 'Europe/Kyiv'  # DATE-OBS у кадрах DWARF — місцевий час: UTC+3 влітку, UTC+2 узимку (перехід 25.10.2026)
SAT = 3500          # сирий пік у ADU (12 біт, стеля 4095)
R_AP, R_IN, R_OUT = 2.0, 4.5, 7.5   # тест 22.10: r=2 найкраще (шум фону міста)
# трансформація TG -> V (transform_coef.py 23.10.2026: NGC 869 + M35, 33 зорі; перехресна перевірка OK)
#   V = g + zp + TV_BV*(B-V);  (B-V) = TBV*(b-g) + zpc;  g, b — інструментальні G і B канали Байєра
TV_BV, TV_BV_ERR, TBV = 0.048, 0.014, 1.94

def log(*a):
    print(*a, flush=True)

# ---------- калібрування ----------
def master_dark(exp, gain, temp):
    d = os.path.join(ASTRO, 'CALI_FRAME', 'dark', 'cam_0')
    rx = re.compile(r'dark_exp_([0-9.]+)_gain_(\d+)_bin_1_(\d+)C')
    c = []
    for f in os.listdir(d):
        m = rx.search(f)
        if m and abs(float(m.group(1)) - exp) < 1e-3 and int(m.group(2)) == gain:
            c.append((abs(int(m.group(3)) - temp), f))
    if not c:
        # нема дарка цієї витримки: найближча витримка того ж gain, тепловий сигнал масштабовано за часом
        alt = []
        for f in os.listdir(d):
            m = rx.search(f)
            if m and int(m.group(2)) == gain:
                alt.append((abs(float(m.group(1)) - exp), abs(int(m.group(3)) - temp), float(m.group(1)), f))
        if alt:
            de_, dt_, e0, _ = min(alt)
            fs = sorted(f for a_, b_, e_, f in alt if e_ == e0 and b_ == dt_)[:10]
            bias = fits.getdata(os.path.join(ASTRO, 'CALI_FRAME', 'bias', 'cam_0', 'bias_gain_2_bin_1.fits')).astype(np.float32)
            dk = np.median([fits.getdata(os.path.join(d, f)).astype(np.float32) for f in fs], axis=0)
            log('  дарк: %gs масштабовано ×%.2f до %gs (%d шт., ΔT=%d°C)' % (e0, exp / e0, exp, len(fs), dt_))
            return bias + (dk - bias) * (exp / e0), dt_
        log('  ! дарків для %gs g%d нема — лише bias' % (exp, gain))
        return fits.getdata(os.path.join(ASTRO, 'CALI_FRAME', 'bias', 'cam_0', 'bias_gain_2_bin_1.fits')).astype(np.float32), None
    best = min(x[0] for x in c)
    fs = sorted(f for dt, f in c if dt == best)[:10]
    log('  дарк: %d шт., ΔT=%d°C (%s…)' % (len(fs), best, fs[0]))
    return np.median([fits.getdata(os.path.join(d, f)).astype(np.float32) for f in fs], axis=0), best

def master_flat():
    bias = fits.getdata(os.path.join(ASTRO, 'CALI_FRAME', 'bias', 'cam_0', 'bias_gain_2_bin_1.fits')).astype(np.float32)
    fl = fits.getdata(os.path.join(ASTRO, 'CALI_FRAME', 'flat', 'cam_0_ir_1', 'flat_gain_2_bin_1_ir_1.fits')).astype(np.float32) - bias
    for i in (0, 1):
        for j in (0, 1):
            fl[i::2, j::2] /= np.median(fl[i::2, j::2])
    return fl

def green(a):
    return 0.5 * (a[0::2, 1::2] + a[1::2, 0::2])

def raw_peak_g(raw, x, y, r=3):
    """Максимум сирих G-пікселів у вікні навколо (x,y) у координатах суперпікселя."""
    xi, yi = int(round(x)), int(round(y))
    sub = raw[2 * (yi - r):2 * (yi + r + 1), 2 * (xi - r):2 * (xi + r + 1)]
    if sub.size == 0:
        return np.nan
    return float(max(sub[0::2, 1::2].max(), sub[1::2, 0::2].max()))

# ---------- допоміжне ----------
def hms(s, mult):
    p = [float(x) for x in s.replace('+', '').split(':')]
    sg = -1 if s.strip().startswith('-') else 1
    v = abs(p[0]) + p[1] / 60 + p[2] / 3600
    return sg * v * mult

def recentroid(img, x, y, r=2):
    xi, yi = int(round(x)), int(round(y))
    b = img[yi - r:yi + r + 1, xi - r:xi + r + 1]
    if b.shape != (2 * r + 1, 2 * r + 1):
        return x, y
    b = np.clip(b - np.median(b), 0, None)
    s = b.sum()
    if s <= 0:
        return x, y
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    dx, dy = (b * xx).sum() / s, (b * yy).sum() / s
    if abs(dx) > 1.5 or abs(dy) > 1.5:
        return x, y
    return xi + dx, yi + dy

def frame_time(h):
    from datetime import datetime
    from zoneinfo import ZoneInfo
    loc = datetime.fromisoformat(h['DATE-OBS'][:26]).replace(tzinfo=ZoneInfo(TZ))
    t = Time(loc.astimezone(ZoneInfo('UTC')).replace(tzinfo=None), scale='utc')
    return t - TimeDelta(float(h['EXPTIME']) / 2, format='sec')   # DATE-OBS = кінець витримки

# ---------- головне ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--session', required=True)
    ap.add_argument('--field', required=True)
    ap.add_argument('--cloud', required=True)
    ap.add_argument('--stars', required=True, help='імена з targets_archive.json через кому')
    ap.add_argument('--bin', type=float, default=180, help='бін, с (0 = 1 точка на сесію)')
    ap.add_argument('--vmin', type=float, default=9.5)
    ap.add_argument('--vmax', type=float, default=13.8)
    ap.add_argument('--ncomp', type=int, default=8)
    ap.add_argument('--limit', type=int, default=0, help='лише перші N кадрів (тест)')
    ap.add_argument('--rap', type=float, default=R_AP, help='радіус апертури, px суперпікселя')
    ap.add_argument('--out-suffix', default='', help='суфікс імені вихідного файлу (напр. _E)')
    ap.add_argument('--force', action='store_true', help='обробити навіть сесію з міткою _FAILED.txt')
    ap.add_argument('--notrans', action='store_true', help='без трансформації: сирий TG, як у v1')
    ap.add_argument('--with-failed', action='store_true', help='брати й кадри failed_ (бортовий стекер відкинув, але кадр може бути добрий)')
    ap.add_argument('--bkgmax', type=float, default=2.0, help='кадр відкидається, якщо фон > bkgmax x медіана (типово 2)')
    a = ap.parse_args()

    sdir = os.path.join(ASTRO, a.session)
    fail_mark = os.path.join(sdir, '_FAILED.txt')
    if os.path.isfile(fail_mark) and not a.force:
        reason = next((l.split(':', 1)[1].strip() for l in open(fail_mark, encoding='utf-8') if l.startswith('reason:')), '?')
        log('[phot] %s: фейлова сесія (_FAILED.txt, reason: %s) — пропуск; --force, щоб обробити' % (a.session, reason))
        sys.exit(3)
    frames = sorted(f for f in glob.glob(os.path.join(sdir, '*.fits'))
                    if (a.with_failed or not os.path.basename(f).startswith('failed')) and 'stacked' not in os.path.basename(f))
    if a.limit:
        frames = frames[:a.limit]
    h0 = fits.getheader(frames[0])
    if str(h0.get('FILTER', '')).strip().lower() != 'astro':
        sys.exit('фільтр %s — не Astro, фотометрію не робимо' % h0.get('FILTER'))
    exp, gain, temp = float(h0['EXPTIME']), int(h0['GAIN']), int(h0.get('DET-TEMP', 25))
    log('[phot] %s: %d кадрів, %gs g%d, T=%d°C' % (a.session, len(frames), exp, gain, temp))
    dark, dT = master_dark(exp, gain, temp)
    flat = master_flat()
    # кадри іншого розміру (напр. 1280x720 на старті сесії) з дарком/флетом 1920x1080 не калібруються
    nf0 = len(frames)
    frames = [f for f in frames if (fits.getheader(f)['NAXIS2'], fits.getheader(f)['NAXIS1']) == dark.shape]
    if len(frames) < nf0:
        log('  відкинуто %d кадрів іншого розміру, ніж дарк %s' % (nf0 - len(frames), dark.shape))

    # хмарний стек: G-площина, біновано 2x2, WCS переведено в напівроздільність
    cd = fits.getdata(a.cloud).astype(np.float32)
    ch = fits.getheader(a.cloud)
    cG = cd[1] if cd.ndim == 3 else cd
    hh, ww = cG.shape[0] // 2 * 2, cG.shape[1] // 2 * 2
    cG2 = cG[:hh, :ww].reshape(hh // 2, 2, ww // 2, 2).mean(axis=(1, 3))
    cG2 = np.ascontiguousarray(cG2 - np.median(cG2))
    w = WCS(ch, naxis=2)

    # WCS хмарного стека задано для необрізаного кадру 1920x1080, а дані обрізані до NAXIS1xNAXIS2 по центру
    offx, offy = (1920 - cG.shape[1]) / 2.0, (1080 - cG.shape[0]) / 2.0

    def sky2half(ra, dec):
        x, y = w.all_world2pix(np.atleast_1d(ra), np.atleast_1d(dec), 0)
        x, y = x - offx, y - offy
        return (x + 0.5) / 2 - 0.5, (y + 0.5) / 2 - 0.5

    tg_all = {t['name']: t for t in json.load(open(os.path.join(HERE, 'targets_archive.json'), encoding='utf-8'))}
    targets = [tg_all[s.strip()] for s in a.stars.split(',')]
    chart = json.load(open(os.path.join(HERE, 'comps_%s.json' % a.field), encoding='utf-8'))
    comps = []
    for s in chart['photometry']:
        v = next((b for b in s['bands'] if b['band'] == 'V'), None)
        bb = next((b for b in s['bands'] if b['band'] == 'B'), None)
        if v and a.vmin <= v['mag'] <= a.vmax:
            comps.append({'label': s['label'], 'auid': s['auid'], 'V': v['mag'], 'Verr': v.get('error') or 0.0,
                          'B': bb['mag'] if bb and bb.get('mag') is not None else None,
                          'ra': hms(s['ra'], 15), 'dec': hms(s['dec'], 1)})
    stars = [{'name': t['name'], 'auid': t['auid'], 'ra': t['ra'], 'dec': t['dec'], 'kind': 'var'} for t in targets] + \
            [dict(c, name='C' + c['label'], kind='comp') for c in comps]
    hx, hy = sky2half([s['ra'] for s in stars], [s['dec'] for s in stars])
    # уточнення зсуву WCS→пікселі за самими еталонами (обрізка хмарного стека не завжди симетрична)
    det = sep.extract(cG2, 8.0, err=sep.Background(cG2).globalrms)
    dxs, dys = [], []
    for s_, x_, y_ in zip(stars, hx, hy):
        if s_['kind'] != 'comp':
            continue
        dd = np.hypot(det['x'] - x_, det['y'] - y_)
        if dd.size and dd.min() < 15:
            k = dd.argmin(); dxs.append(det['x'][k] - x_); dys.append(det['y'][k] - y_)
    if len(dxs) >= 3:
        mdx, mdy = np.median(dxs), np.median(dys)
        sc = 1.4826 * np.median(np.hypot(np.array(dxs) - mdx, np.array(dys) - mdy))
        log('  зсув WCS за %d еталонами: (%.2f, %.2f) px, розкид %.2f' % (len(dxs), mdx, mdy, sc))
        if sc < 1.5:
            hx, hy = hx + mdx, hy + mdy
        else:
            log('  ! зсув WCS неоднозначний — не застосовано')
    H, W = cG2.shape
    keep = []
    for s, x, y in zip(stars, hx, hy):
        if 12 <= x <= W - 12 and 12 <= y <= H - 12:
            s['cx'], s['cy'] = float(x), float(y)
            keep.append(s)
        elif s['kind'] == 'var':
            log('  ! %s поза кадром хмарного стека (x=%.0f y=%.0f) — пропуск' % (s['name'], x, y))
    stars = keep
    log('  зір у кадрі: змінних %d, еталонів %d' % (sum(s['kind'] == 'var' for s in stars), sum(s['kind'] == 'comp' for s in stars)))
    src = np.array([[s['cx'], s['cy']] for s in stars])

    rows = []
    tf_prev = None
    for k, f in enumerate(frames):
        raw = fits.getdata(f).astype(np.float32)
        hd = fits.getheader(f)
        cal = (raw - dark) / flat
        g = np.ascontiguousarray(green(cal))
        bk = sep.Background(g, bw=32, bh=32)
        gs = np.ascontiguousarray(g - bk.back())
        try:
            tf, _ = aa.find_transform(gs, cG2, detection_sigma=5, max_control_points=50)
            tf_prev = tf
        except Exception as e:
            if tf_prev is None:
                log('  %s: вирівнювання не вдалось (%s)' % (os.path.basename(f), e)); continue
            tf = tf_prev
        pos = tf.inverse(src)
        xs, ys, peaks = [], [], []
        for si, (x, y) in enumerate(pos):
            # змінні — у фіксованій позиції з трансформації (залишок ~0,2 px): у мінімумі вони слабкі й центроїд блукає
            x2, y2 = (x, y) if stars[si]['kind'] == 'var' else recentroid(gs, x, y)
            xs.append(x2); ys.append(y2); peaks.append(raw_peak_g(raw, x2, y2))
        # фон — з карти sep.Background (bw 32); кільце sep не робастне до сусідніх зір (HR Gem: сусід 7 px)
        flux, ferr, flag = sep.sum_circle(gs, xs, ys, a.rap, err=bk.globalrms)
        # зорі біля краю кадру (дрейф/поворот поля) — не міряти
        Hf, Wf = gs.shape
        edge = (np.array(xs) < 12) | (np.array(xs) > Wf - 13) | (np.array(ys) < 12) | (np.array(ys) > Hf - 13)
        flux = np.where(edge, np.nan, flux)
        # B-канал Байєра (сирі (1,1)): центр зсунутий на −0,25 суперпікселя відносно G — лише для кольору цілі
        bpl = np.ascontiguousarray(cal[1::2, 1::2])
        bkb = sep.Background(bpl, bw=32, bh=32)
        fluxb, _, _ = sep.sum_circle(np.ascontiguousarray(bpl - bkb.back()), np.array(xs) - 0.25, np.array(ys) - 0.25,
                                     a.rap, err=bkb.globalrms)
        fluxb = np.where(edge, np.nan, fluxb)
        t = frame_time(hd)
        rows.append({'file': os.path.basename(f), 'jd': t.jd, 'bkg': float(np.median(bk.back())), 'rms': float(bk.globalrms),
                     'flux': flux.tolist(), 'ferr': ferr.tolist(), 'peak': peaks, 'fluxb': fluxb.tolist()})
        if k % 25 == 0:
            log('  кадр %d/%d  фон %.0f' % (k + 1, len(frames), rows[-1]['bkg']))
    if not rows:
        sys.exit('нема кадрів')

    F = np.array([r['flux'] for r in rows]); P = np.array([r['peak'] for r in rows])
    jd = np.array([r['jd'] for r in rows]); bkg = np.array([r['bkg'] for r in rows])
    ci = [i for i, s in enumerate(stars) if s['kind'] == 'comp']
    # еталони: без насичення (>5 % кадрів), стабільні, з нормальним потоком
    good = []
    for i in ci:
        if np.nanmean(P[:, i] > SAT) > 0.05 or np.nanmedian(F[:, i]) <= 0 or np.mean(~np.isfinite(F[:, i])) > 0.02:
            continue
        good.append(i)
    # прозорість кадру: медіана відношення потоку еталонів до їхньої медіани
    rel = np.nanmedian(F[:, good] / np.nanmedian(F[:, good], axis=0), axis=1)
    okf = (rel > 0.6) & (bkg < a.bkgmax * np.median(bkg))
    log('  кадрів добрих: %d/%d (хмари/фон відкинуто %d)' % (okf.sum(), len(rows), (~okf).sum()))
    # локальне затінення (козирок, край поля, роса): зоря «гасне» разом із сусідніми еталонами
    R = F[:, good] / np.nanmedian(F[okf][:, good], axis=0)
    Rf = np.nanmedian(R, axis=1, keepdims=True)
    F[:, good] = np.where(R < 0.7 * Rf, np.nan, F[:, good])        # еталон, що просів відносно поля
    gxy = np.array([[stars[i]['cx'], stars[i]['cy']] for i in good])
    nmask = 0
    for vi_, s_ in enumerate(stars):
        if s_['kind'] != 'var' and vi_ not in good:
            continue
        dd = np.hypot(gxy[:, 0] - s_['cx'], gxy[:, 1] - s_['cy'])
        near = [k for k in np.argsort(dd) if dd[k] > 0.5][:3]
        loc = np.nanmedian(R[:, near], axis=1) / Rf[:, 0]
        bad = loc < 0.8
        nmask += int(bad.sum())
        F[bad, vi_] = np.nan
    log('  маска затінення: %d зоре-кадрів' % nmask)
    FB = np.array([r['fluxb'] for r in rows])
    FB = np.where(np.isfinite(F), FB, np.nan)                   # ті самі маски, що й для G
    tr = not a.notrans

    def colour(si, ii):
        """Інструментальний колір b-g зорі si за кадрами ii (сумарні потоки) і його похибка; None, якщо нема сигналу."""
        fg, fb = F[ii, si], FB[ii, si]
        ok_ = np.isfinite(fg) & np.isfinite(fb)
        n_ = int(ok_.sum())
        if n_ < 2 or np.sum(fb[ok_]) <= 0 or np.sum(fg[ok_]) <= 0:
            return None
        fg, fb = fg[ok_], fb[ok_]
        bg = -2.5 * np.log10(np.sum(fb) / np.sum(fg))
        sbg = 1.0857 * np.hypot(np.std(fb, ddof=1) / np.sqrt(n_) / np.mean(fb), np.std(fg, ddof=1) / np.sqrt(n_) / np.mean(fg))
        return float(bg), float(sbg)

    def bv_of(si, ii, ens_):
        """(B-V) зорі si: TBV*(b-g) + нуль-пункт кольору з ансамблю в тих самих кадрах."""
        zc_ = []
        for c in ens_:
            cc = colour(c, ii)
            if cc is not None and cc[1] < 0.1:
                zc_.append((stars[c]['B'] - stars[c]['V']) - TBV * cc[0])
        ct = colour(si, ii)
        if ct is None or len(zc_) < 3:
            return None
        return TBV * ct[0] + float(np.median(zc_)), TBV * ct[1]

    out = {'session': a.session, 'field': a.field, 'chart': chart.get('chartid'), 'exp': exp, 'gain': gain,
           'dark_dT': dT, 'n_frames': len(rows), 'n_good': int(okf.sum()), 'bin_s': a.bin, 'stars': [],
           'filt': 'V' if tr else 'TG',
           'tcoef': {'Tv_bv': TV_BV, 'Tv_bv_err': TV_BV_ERR, 'Tbv': TBV} if tr else None,
           'frame_diag': [{'jd': float(j_), 'bkg': float(b_), 'rel': float(r_), 'ok': bool(o_)} for j_, b_, r_, o_ in zip(jd, bkg, rel, okf)]}
    for vi, s in enumerate(stars):
        if s['kind'] != 'var':
            continue
        if np.nanmean(P[okf, vi] > SAT) > 0.02:
            log('  ! %s насичена в >2 %% кадрів — пропуск' % s['name']); continue
        vmag_inst = -2.5 * np.log10(np.nanmedian(F[okf, vi]))
        # ансамбль — ncomp еталонів, найближчих за блиском; контрольна — наступна
        # ансамбль — найяскравіші ненасичені еталони (у міському небі шум фону домінує, слабкі еталони шумлять);
        # контрольна — найближча за блиском до змінної серед решти
        # K-кандидати: еталони, чий нуль-пункт узгоджений з полем (≤0,2m) — без блендів і помилкових ідентифікацій
        zc = {i: stars[i]['V'] + 2.5 * np.log10(np.nanmedian(F[okf, i])) for i in good}
        zmed = np.median(list(zc.values()))
        kc = [i for i in good if abs(zc[i] - zmed) <= 0.2]
        chk = min(kc, key=lambda i: abs(-2.5 * np.log10(np.nanmedian(F[okf, i])) - vmag_inst)) if len(kc) > 3 else None
        cand = sorted([i for i in good if i != chk and abs(zc[i] - zmed) <= 0.2 and (not tr or stars[i]['B'] is not None)],
                      key=lambda i: -np.nanmedian(F[okf, i]))
        ens = cand[:a.ncomp]
        # з трансформацією нуль-пункт еталона — V - TV_BV*(B-V) (тоді V цілі = g + zp + TV_BV*(B-V)ціль)
        Vc = np.array([stars[i]['V'] - (TV_BV * (stars[i]['B'] - stars[i]['V']) if tr else 0.0) for i in ens])
        zpf = Vc[None, :] + 2.5 * np.log10(F[:, ens])         # нуль-пункт кожного еталона в кожному кадрі
        # вага еталона = 1/σ² його zp у часі (відносно медіани ансамблю); змінні/погані — геть
        for _ in range(2):
            dev = zpf - np.nanmedian(zpf, axis=1, keepdims=True)
            sdv = np.nanstd(dev[okf] - np.nanmedian(dev[okf], axis=0), axis=0)
            use = sdv < max(0.03, 2.5 * np.nanmedian(sdv))
            ens = [e for e, u_ in zip(ens, use) if u_]; zpf = zpf[:, use]; sdv = sdv[use]
        wts = 1.0 / np.maximum(sdv, 0.005) ** 2
        off = np.nanmedian(zpf[okf], axis=0)                    # сталий колірний зсув кожного еталона
        zp_rel = np.nansum((zpf - off) * wts, axis=1) / np.nansum(np.isfinite(zpf) * wts, axis=1)
        zp = zp_rel + np.median(off)                           # абсолютний нуль-пункт — медіана (стійка до бленду)
        zps = np.full_like(zp, 1.0 / np.sqrt(np.sum(wts)))       # похибка ансамблю в кадрі
        log('    діагн.: ансамбль %s; σ еталонів у часі: %s' % (len(ens),
            ' '.join('%s:%.3f' % (stars[e]['label'], d) for e, d in zip(ens, sdv))))
        # бінінг
        # бінінг у потоці (слабка змінна в мінімумі може давати від'ємні кадрові потоки — їх не викидаємо)
        fs_ = F[:, vi] * 10 ** (-0.4 * zp)
        fk_ = (F[:, chk] * 10 ** (-0.4 * zp)) if chk is not None else np.full_like(fs_, np.nan)
        sel = okf & np.isfinite(fs_)
        if tr:
            # колір за всю сесію — запасний варіант, коли в біні B-сигналу замало
            bv_sess = bv_of(vi, np.where(sel)[0], ens)
            bk_sess = bv_of(chk, np.where(sel)[0], ens) if chk is not None else None
            bvc_med = float(np.median([stars[c]['B'] - stars[c]['V'] for c in ens]))
            if bv_sess is None:
                log('  ! %s: колір не визначається (B-сигнал) — пропуск у V' % s['name']); continue
            log('    колір %s за сесію: B-V %.2f ± %.2f (ансамбль %.2f)' % (s['name'], bv_sess[0], bv_sess[1], bvc_med))
        if a.bin > 0:
            edges = np.arange(jd[sel].min(), jd[sel].max() + a.bin / 86400 + 1e-9, a.bin / 86400)
            idx = np.digitize(jd, edges)
        else:
            idx = np.where(sel, 1, 0)
        bins = []
        for b in np.unique(idx[sel]):
            j = sel & (idx == b)
            if j.sum() < (3 if a.bin > 0 else 1):
                continue
            ff = fs_[j]; med = np.median(ff); sd0 = np.std(ff)
            keepb = np.abs(ff - med) < max(3 * sd0, 1e-12)
            ff = ff[keepb]; n = len(ff)
            if n < 2 or np.mean(ff) <= 0:
                continue
            mf = np.mean(ff)
            err = float(np.sqrt((1.0857 * np.std(ff, ddof=1) / np.sqrt(n) / mf) ** 2 + np.mean(zps[j]) ** 2 / n + 0.005 ** 2))
            jm = float(np.mean(jd[j][keepb]))
            kf = np.nanmean(fk_[j])
            altaz = SkyCoord(s['ra'] * u.deg, s['dec'] * u.deg).transform_to(AltAz(obstime=Time(jm, format='jd'), location=SITE))
            mag = float(-2.5 * np.log10(mf)); kmag = float(-2.5 * np.log10(kf)) if kf > 0 else float('nan')
            b_ = {'jd': jm, 'mag': mag, 'err': err, 'n': n, 'kmag': kmag, 'airmass': float(altaz.secz)}
            if tr:
                ii = np.where(j)[0][keepb]
                bv = bv_of(vi, ii, ens)
                if bv is None or bv[1] > 0.10:
                    bv = bv_sess                                  # колір біну шумний (B-сигнал слабкий) — колір сесії
                kb = bv_of(chk, ii, ens) if chk is not None else None
                if kb is None or kb[1] > 0.10:
                    kb = bk_sess                                  # те саме правило для контрольної зорі
                b_.update(mag_tg=mag, kmag_tg=kmag, bv=bv[0], bv_err=bv[1],
                          mag=mag + TV_BV * bv[0],
                          err=float(np.sqrt(err ** 2 + (TV_BV * bv[1]) ** 2 + (TV_BV_ERR * (bv[0] - bvc_med)) ** 2)),
                          kmag=(kmag + TV_BV * kb[0]) if kb is not None else float('nan'))
            bins.append(b_)
        kname = stars[chk]['auid'] if chk is not None else None
        kcat = stars[chk]['V'] if chk is not None else None
        ks = np.nanstd([b['kmag'] for b in bins]) if bins else np.nan
        km_mean = np.nanmean([b['kmag'] for b in bins]) if bins else np.nan
        if not bins:
            log('  ! %s: жодного біну (потік ≤0?) — пропуск' % s['name']); continue
        log('  %s: бінів %d, mag %.2f..%.2f, медіанна похибка %.3f | ансамбль %d, K %s %.3f (кат. %.3f, Δ%+.3f, σ %.3f)' % (
            s['name'], len(bins), min(b['mag'] for b in bins), max(b['mag'] for b in bins), np.median([b['err'] for b in bins]),
            len(ens), stars[chk]['label'] if chk is not None else '-', km_mean, kcat or 0, (km_mean - kcat) if kcat else 0, ks))
        out['stars'].append({'name': s['name'], 'auid': s['auid'], 'ensemble': [stars[i]['label'] for i in ens],
                             'ensemble_auid': [stars[i]['auid'] for i in ens], 'check_auid': kname,
                             'check_label': stars[chk]['label'] if chk is not None else None, 'check_cat': kcat,
                             'check_mean': float(km_mean), 'check_sd': float(ks), 'bins': bins})
    tag = re.search(r'(\d{4}-\d{2}-\d{2})-(\d{2})-(\d{2})', a.session)
    stem = '%s_%s%s%s' % (a.field, tag.group(1).replace('-', ''), tag.group(2), tag.group(3))
    stem += a.out_suffix
    odir = os.path.join(HERE, 'out_v' if tr else 'out')         # V окремо: TG-результати v1 в out\ не перезаписуємо
    os.makedirs(odir, exist_ok=True)
    jf = os.path.join(odir, stem + '.json')
    json.dump(out, open(jf, 'w', encoding='utf-8'), indent=1)
    log('[phot] -> %s' % jf)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        for st in out['stars']:
            fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
            x = [(b['jd'] - 2460000) * 24 for b in st['bins']]
            ax[0].errorbar(x, [b['mag'] for b in st['bins']], [b['err'] for b in st['bins']], fmt='o', ms=3)
            ax[0].invert_yaxis(); ax[0].set_ylabel(out['filt']); ax[0].set_title('%s  %s  KTFA' % (st['name'], stem))
            ax[1].plot(x, [b['kmag'] for b in st['bins']], 'k.', ms=3); ax[1].invert_yaxis()
            ax[1].set_ylabel('K %s' % st['check_label']); ax[1].set_xlabel('год від JD 2460000')
            fig.tight_layout(); fig.savefig(os.path.join(odir, '%s_%s.png' % (stem, st['name'].replace(' ', '_'))), dpi=110)
            plt.close(fig)
    except Exception as e:
        log('  (графік не зроблено: %s)' % e)

if __name__ == '__main__':
    main()
