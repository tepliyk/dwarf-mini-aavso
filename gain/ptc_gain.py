# -*- coding: utf-8 -*-
"""Орієнтовний system gain DWARF Mini (e-/ADU) по каналах R, G, B методом photon transfer (26.10.2026).
Пари сусідніх кадрів однієї сесії: var(різниці)/2 = RN^2 + S/g  ->  нахил по рівню неба S = 1/g.
Зсув bias на нахил не впливає, тому окремі bias-кадри не потрібні. Зорі відкидаються робастною (MAD) оцінкою.
Тільки читає FITS; пише ptc_gain.json і ptc_gain.png поруч.
  python -X utf8 ptc_gain.py "<тека сесії>" ["<тека сесії>" ...]"""
import glob, json, os, sys
import numpy as np
from astropy.io import fits
sys.stdout.reconfigure(encoding='utf-8')

HERE = os.path.dirname(os.path.abspath(__file__))
CH = {'R': (0, 0), 'G1': (0, 1), 'G2': (1, 0), 'B': (1, 1)}   # RGGB з хедера BAYERPAT
PED = 200.0                                                    # пьєдестал, лише для масштабування пари
SAT = 4000                                                     # 12-біт АЦП, стеля 4095

def frames(d):
    return sorted(f for f in glob.glob(os.path.join(d, '*.fits'))
                  if not os.path.basename(f).startswith(('failed', 'stacked')))

def robust_var(x):
    x = x[np.isfinite(x)]
    mad = np.median(np.abs(x - np.median(x)))
    return (1.4826 * mad) ** 2

T = 25   # плитка 25x25 пікселів одного кольору: градієнт сутінків і віньєтка всередині ~ сталі

def tile_var(p, q):
    """Медіана по плитках: q масштабується до рівня p локально, далі MAD-дисперсія різниці / 2."""
    vs = []
    for y in range(0, p.shape[0] - T + 1, T):
        for x in range(0, p.shape[1] - T + 1, T):
            a, b = p[y:y + T, x:x + T], q[y:y + T, x:x + T]
            ma, mb = np.median(a), np.median(b)
            if mb - PED < 5:
                continue
            vs.append(robust_var(a - (PED + (b - PED) * (ma - PED) / (mb - PED))) / 2.0)
    return float(np.median(vs))

pts = {k: [] for k in CH}
meta = []
for d in sys.argv[1:]:
    fs = frames(d)
    h = fits.getheader(fs[0])
    assert h.get('BAYERPAT', 'RGGB').strip() == 'RGGB', h.get('BAYERPAT')
    meta.append({'dir': os.path.basename(d), 'n': len(fs), 'gain': h['GAIN'], 'exp': h['EXPTIME'],
                 'filter': h.get('FILTER'), 'fw': h.get('FIRMWARE')})
    prev = None
    for f in fs:
        a = fits.getdata(f).astype(np.float64)[140:940, 460:1460]   # центр 1000x800, без віньєтування країв
        if prev is not None:
            for k, (y, x) in CH.items():
                p, q = prev[y::2, x::2], a[y::2, x::2]
                mp, mq = np.median(p), np.median(q)
                if max(mp, mq) > SAT or (mp - PED) < 20:
                    continue
                pts[k].append((0.5 * (mp + mq), tile_var(p, q)))
        prev = a

res = {'sessions': meta, 'channels': {}}
for k, v in pts.items():
    if len(v) < 5:
        res['channels'][k] = {'n_pairs': len(v)}; print(k, 'мало пар:', len(v)); continue
    s, var = np.array(v).T
    A = np.vstack([s, np.ones_like(s)]).T
    keep = np.ones(len(s), bool)
    for _ in range(3):                                   # 3-sigma ітеративне відкидання
        (m, c), *_ = np.linalg.lstsq(A[keep], var[keep], rcond=None)
        r = var - (m * s + c); sd = np.std(r[keep])
        keep = np.abs(r) < 3 * sd
    cov = np.linalg.inv(A[keep].T @ A[keep]) * np.var(r[keep], ddof=2)
    g = 1 / m; ge = np.sqrt(cov[0, 0]) / m ** 2
    rn2 = c + m * PED                                     # дисперсія при S = пьєдестал ≈ шум читання^2
    res['channels'][k] = {'n_pairs': int(keep.sum()), 'gain_e_per_adu': round(float(g), 3),
                          'gain_err': round(float(ge), 3), 'S_range_adu': [int(s.min()), int(s.max())],
                          'read_noise_adu': round(float(np.sqrt(rn2)), 2) if rn2 > 0 else None}
    print('%-2s пар %3d  gain %.3f ± %.3f e-/ADU  рівень %d..%d ADU  RN %s ADU' % (
        k, keep.sum(), g, ge, s.min(), s.max(), res['channels'][k]['read_noise_adu']))

json.dump(res, open(os.path.join(HERE, 'ptc_gain.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    col = {'R': 'tab:red', 'G1': 'tab:green', 'G2': 'olive', 'B': 'tab:blue'}
    for k, v in pts.items():
        if v:
            s, var = np.array(v).T; plt.plot(s, var, '.', color=col[k], label=k, ms=4)
    plt.xlabel('mean level, ADU (incl. ~200 pedestal)'); plt.ylabel('var(diff)/2, ADU^2')
    plt.title('DWARF Mini photon transfer, dawn 2026-10-07, gain 60'); plt.legend(); plt.grid(alpha=.3)
    plt.savefig(os.path.join(HERE, 'ptc_gain.png'), dpi=110, bbox_inches='tight')
except Exception as e:
    print('графік не збудовано:', e)
print('записано ptc_gain.json / ptc_gain.png')
