# -*- coding: utf-8 -*-
"""Огляд сутінкової серії 26.10: по кожній теці число кадрів, DATE-OBS першого/останнього, медіани R/G/B, насичення."""
import glob, os, sys
import numpy as np
from astropy.io import fits
sys.stdout.reconfigure(encoding='utf-8')
ROOT = r'C:\astro\DWARF Mini\Astronomy'
dirs = sorted(glob.glob(os.path.join(ROOT, 'DWARF_RAW_TELE_EXP_0.25_GAIN_0_2026-10-26*')))
tot = 0
print('%-19s %4s %-8s %-8s %5s %5s %5s %5s' % ('тека', 'n', 'перший', 'останній', 'R', 'G', 'B', 'sat%'))
for d in dirs:
    fs = sorted(f for f in glob.glob(os.path.join(d, '*.fits')) if 'stacked' not in os.path.basename(f))
    if not fs:
        print(os.path.basename(d)[-23:], 'порожня'); continue
    h0, h1 = fits.getheader(fs[0]), fits.getheader(fs[-1])
    a = fits.getdata(fs[len(fs) // 2]).astype(np.float32)
    tot += len(fs)
    print('%-19s %4d %-8s %-8s %5d %5d %5d %5.1f' % (
        os.path.basename(d)[-23:-4], len(fs), h0['DATE-OBS'][11:19], h1['DATE-OBS'][11:19],
        np.median(a[0::2, 0::2]), np.median(a[0::2, 1::2]), np.median(a[1::2, 1::2]), 100 * np.mean(a >= 4090)))
print('разом кадрів', tot, 'у', len(dirs), 'теках')
