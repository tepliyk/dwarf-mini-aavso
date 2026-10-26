# -*- coding: utf-8 -*-
"""Швидка перевірка тестових кадрів сутінкової серії: рівні по каналах RGGB, частка насичених пікселів. Тільки читає."""
import glob, os, sys
import numpy as np
from astropy.io import fits
sys.stdout.reconfigure(encoding='utf-8')
CH = {'R': (0, 0), 'G': (0, 1), 'B': (1, 1)}
for d in sys.argv[1:]:
    fs = sorted(glob.glob(os.path.join(d, '*.fits')))
    print('==', os.path.basename(d), '| файлів', len(fs))
    if not fs:
        continue
    h = fits.getheader(fs[0])
    print('   CAMERA %s FILTER %s GAIN %s EXPTIME %s DATE-OBS %s OBJECT %s' % (
        h.get('CAMERA'), h.get('FILTER'), h.get('GAIN'), h.get('EXPTIME'), h.get('DATE-OBS'), h.get('OBJECT')))
    for f in fs[:3] + fs[-2:] if len(fs) > 5 else fs:
        a = fits.getdata(f).astype(np.float32)
        s = []
        for k, (y, x) in CH.items():
            c = a[y::2, x::2]
            s.append('%s med %4d p99 %4d sat %.1f%%' % (k, np.median(c), np.percentile(c, 99), 100 * np.mean(c >= 4090)))
        print('   %-55s %s | max %d' % (os.path.basename(f)[:55], ' | '.join(s), a.max()))
