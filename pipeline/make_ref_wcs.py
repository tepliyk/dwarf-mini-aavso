#!/usr/bin/env python3
"""
make_ref_wcs.py — опорний кадр з WCS для полів без хмарного стека (22.10.2026).

Медіана перших N каліброваних кадрів сесії (зелений суперпіксель) -> зорі sep ->
зіставлення з еталонами VSP (comps_<field>.json, проекція на дотичну площину, 8,03"/px,
перебір парності) через astroalign -> fit_wcs_from_points -> 2D FITS 1920x1080 (kron 2x)
з WCS для повної роздільності. Далі цей файл іде в aavso_phot.py як --cloud.

  python -X utf8 make_ref_wcs.py --session "DWARF_RAW_TELE_M 77_EXP_30_GAIN_60_2026-10-10-23-40-49-979" --field M77
"""
import argparse, glob, json, os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs.utils import fit_wcs_from_points
from astropy.coordinates import SkyCoord
import astropy.units as u
import astroalign as aa
import sep
import aavso_phot as P

SCALE_HALF = 8.03  # "/px зеленого суперпікселя

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--session', required=True)
    ap.add_argument('--field', required=True)
    ap.add_argument('--n', type=int, default=15)
    ap.add_argument('--maglim', type=float, default=14.5)
    ap.add_argument('--ndet', type=int, default=80, help='скільки найяскравіших зір кадру брати')
    ap.add_argument('--frame', type=float, default=0.5, help='частка сесії, звідки брати кадр (0..1)')
    ap.add_argument('--cat', choices=['vsp', 'ucac4'], default='vsp', help='каталог для астрометрії')
    ap.add_argument('--ncat', type=int, default=100, help='скільки найяскравіших зір каталогу (ucac4)')
    ap.add_argument('--amax', type=float, default=360, help='перебір кута 0..amax і 360-amax..360 (EQ: ~10)')
    a = ap.parse_args()
    sdir = os.path.join(P.ASTRO, a.session)
    frames = sorted(f for f in glob.glob(os.path.join(sdir, '*.fits'))
                    if not os.path.basename(f).startswith('failed') and 'stacked' not in os.path.basename(f))
    # один кадр із середини сесії (медіана кількох дає змазування через дрейф поля між кадрами)
    frames = [frames[min(len(frames) - 1, int(len(frames) * a.frame))]]
    h0 = fits.getheader(frames[0])
    exp, gain, temp = float(h0['EXPTIME']), int(h0['GAIN']), int(h0.get('DET-TEMP', 25))
    dark, _ = P.master_dark(exp, gain, temp)
    flat = P.master_flat()
    G = np.median([P.green((fits.getdata(f).astype(np.float32) - dark) / flat) for f in frames], axis=0)
    G = np.ascontiguousarray(G, dtype=np.float32)
    bk = sep.Background(G, bw=32, bh=32)
    Gs = np.ascontiguousarray(G - bk.back())
    det = sep.extract(Gs, 5.0, err=bk.globalrms, minarea=4)
    det = det[(det['flag'] == 0)]
    det = det[np.argsort(-det['flux'])][:a.ndet]
    dpts = np.column_stack([det['x'], det['y']])
    print('[ref] %d кадрів, зір виявлено %d' % (len(frames), len(dpts)))

    cat = []
    if a.cat == 'ucac4':
        # UCAC4 з VizieR: поле 140'x80' навколо наведення, V<12.5, найяскравіші a.ncat
        import urllib.request
        url = ('https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=I/322A/out&-c=%.5f%+.5f&-c.bm=140x80'
               '&-out=RAJ2000,DEJ2000,Vmag&Vmag=<12.5&-out.max=2000' % (float(h0['RA']), float(h0['DEC'])))
        t = urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (KTFA)'}), timeout=90).read().decode()
        for l in t.splitlines():
            p = l.split('\t')
            if len(p) == 3 and p[0][:1].isdigit():
                try:
                    cat.append((float(p[0]), float(p[1]), float(p[2])))
                except ValueError:
                    pass
        cat = sorted(set(cat), key=lambda c: c[2])[:a.ncat]
        print('[ref] UCAC4: %d зір' % len(cat))
    else:
        chart = json.load(open(os.path.join(P.HERE, 'comps_%s.json' % a.field), encoding='utf-8'))
        for s in chart['photometry']:
            v = next((b['mag'] for b in s['bands'] if b['band'] == 'V'), None)
            if v and v <= a.maglim:
                cat.append((P.hms(s['ra'], 15), P.hms(s['dec'], 1), v))
    cat = np.array(cat)
    c0 = SkyCoord(float(h0.get('RA', np.mean(cat[:, 0]))) * u.deg, float(h0.get('DEC', np.mean(cat[:, 1]))) * u.deg)
    sc = SkyCoord(cat[:, 0] * u.deg, cat[:, 1] * u.deg)
    off = c0.spherical_offsets_to(sc)
    xi, eta = off[0].arcsec / SCALE_HALF, off[1].arcsec / SCALE_HALF
    # перебір: парність × кут (крок 0,5°) × масштаб; зсув — голосуванням різниць «зоря − каталог»
    best = None
    for par in (1, -1):
        for sc_ in (0.98, 1.0, 1.02):
            for ang in [x for x in np.arange(0, 360, 0.5) if x <= a.amax or x >= 360 - a.amax]:
                t = np.radians(ang)
                x0, y0 = par * xi * sc_, -eta * sc_
                cx = x0 * np.cos(t) - y0 * np.sin(t); cy = x0 * np.sin(t) + y0 * np.cos(t)
                dx = (dpts[None, :, 0] - cx[:, None]).ravel(); dy = (dpts[None, :, 1] - cy[:, None]).ravel()
                Hh, xe, ye = np.histogram2d(dx, dy, bins=[np.arange(-200, 1200, 3), np.arange(-200, 800, 3)])
                v = Hh.max()
                if best is None or v > best[0]:
                    i, j = np.unravel_index(Hh.argmax(), Hh.shape)
                    best = (v, par, sc_, ang, xe[i] + 1.5, ye[j] + 1.5, np.column_stack([cx, cy]))
    v, par, sc_, ang, tx, ty, crot = best
    print('  найкраще: парність %+d, масштаб %.2f, кут %.1f°, зсув (%.0f, %.0f), голосів %d' % (par, sc_, ang, tx, ty, v))
    if v < 6:
        sys.exit('замало збігів — WCS не побудовано')
    pc = crot + [tx, ty]
    for _ in range(3):   # уточнення афінним МНК за збігами
        dd = np.hypot(*(pc[:, None, :] - dpts[None]).transpose(2, 0, 1))
        k = np.argmin(dd, axis=1); okm = dd[np.arange(len(k)), k] < 4.0
        A = np.column_stack([crot[okm], np.ones(okm.sum())])
        M, *_ = np.linalg.lstsq(A, dpts[k[okm]], rcond=None)
        pc = np.column_stack([crot, np.ones(len(crot))]) @ M
    dd = np.hypot(*(pc[:, None, :] - dpts[None]).transpose(2, 0, 1))
    k = np.argmin(dd, axis=1); ok = dd[np.arange(len(k)), k] < 2.0
    xs = 2 * dpts[k[ok], 0] + 0.5; ys = 2 * dpts[k[ok], 1] + 0.5
    w = fit_wcs_from_points((xs, ys), sc[ok], projection='TAN')
    rx, ry = w.world_to_pixel(sc[ok])
    rms = float(np.sqrt(np.mean((rx - xs) ** 2 + (ry - ys) ** 2)))
    print('[ref] WCS за %d зорями, rms %.2f px повної роздільності, масштаб %.3f"/px' % (
        ok.sum(), rms, float(np.mean([q.to(u.arcsec).value for q in w.proj_plane_pixel_scales()]))))
    if rms > 1.5:
        sys.exit('rms завеликий — перевірити')
    full = np.kron(G, np.ones((2, 2), dtype=np.float32))
    hdr = w.to_header(); hdr['OBJECT'] = a.field; hdr['FILTER'] = 'Astro'; hdr['HISTORY'] = 'make_ref_wcs.py ' + a.session
    os.makedirs(os.path.join(P.HERE, 'ref'), exist_ok=True)
    outf = os.path.join(P.HERE, 'ref', '%s_ref.fits' % a.field)
    fits.writeto(outf, full, hdr, overwrite=True)
    print('[ref] ->', outf, full.shape)

if __name__ == '__main__':
    main()
