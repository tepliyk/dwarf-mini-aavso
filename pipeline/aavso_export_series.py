#!/usr/bin/env python3
"""
aavso_export_series.py — out/*.json (aavso_phot.py) -> AAVSO Extended Format (беклог A4).
Правила якості: K-зоря |Δ| ≤ 0,10 і σ ≤ max(0,04; 1,5×мед. похибки); біни з K-викидом (>3σ роб.) і err>3×медіани — геть; точка йде, якщо err ≤ 0,10.
  python -X utf8 aavso_export_series.py out\\M33_202608260013.json --stars "V Tri" [--write]
"""
import argparse, json, os, sys
OBSCODE = 'KTFA'
K_DMAX, K_SDMAX, ERR_MAX = 0.10, 0.04, 0.10

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+')
    ap.add_argument('--stars', default='', help='імена через кому; порожньо — усі')
    ap.add_argument('--write', default='', help='шлях вихідного .txt')
    ap.add_argument('--force', action='store_true', help='ігнорувати ворота K-зорі')
    ap.add_argument('--live', action='store_true', help='нові (не архівні) спостереження: без слова archival у NOTES')
    a = ap.parse_args()
    want = {s.strip() for s in a.stars.split(',') if s.strip()}
    filts = {json.load(open(fn, encoding='utf-8')).get('filt', 'TG') for fn in a.files}
    if len(filts) != 1:
        sys.exit('змішані фільтри у вхідних файлах: %s — експортуй V і TG окремо' % filts)
    filt = filts.pop()
    sw = ('aavso_phot.py v2 (sep+astroalign; VSP ensemble; G+B channels; transformed to V)' if filt == 'V'
          else 'aavso_phot.py v1 (sep+astroalign; VSP ensemble; green channel)')
    hdr = ['#TYPE=EXTENDED', '#OBSCODE=%s' % OBSCODE,
           '#SOFTWARE=%s' % sw,
           '#DELIM=,', '#DATE=JD', '#OBSTYPE=DSLR',
           '#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,GROUP,CHART,NOTES']
    rows, group = [], 0
    for fn in a.files:
        d = json.load(open(fn, encoding='utf-8'))
        for st in d['stars']:
            if want and st['name'] not in want:
                continue
            import numpy as np
            kk = np.array([b['kmag'] for b in st['bins']], float)
            kmed = np.nanmedian(kk); ksd = 1.4826 * np.nanmedian(np.abs(kk - kmed))     # робастна σ K
            emed = np.median([b['err'] for b in st['bins']])
            good = [b for b in st['bins'] if b['err'] <= min(ERR_MAX, 3 * emed)
                    and np.isfinite(b['kmag']) and abs(b['kmag'] - kmed) <= max(3 * ksd, 0.10)]
            # викиди змінної: відхилення від ковзної медіани (5 бінів) > max(0,15; 5σ роб.) — край кадру, супутник, хмара
            if len(good) >= 5:
                mg = np.array([b['mag'] for b in good])
                rm = np.array([np.median(mg[max(0, i - 2):i + 3]) for i in range(len(mg))])
                rs = 1.4826 * np.median(np.abs(mg - rm))
                good = [b for b, d_ in zip(good, np.abs(mg - rm)) if d_ <= max(0.15, 5 * rs)]
            kg = np.array([b['kmag'] for b in good]) if good else np.array([np.nan])
            dk = float(np.mean(kg) - st['check_cat']) if st['check_cat'] is not None else 99
            st['check_sd'] = float(np.std(kg))
            # σ K має відповідати шуму (K підбирається близько до блиску змінної): ≤ max(0,04; 1,5 × медіанна похибка)
            ok = abs(dk) <= K_DMAX and st['check_sd'] <= max(K_SDMAX, 1.5 * emed)
            print('%-28s %-10s бінів %3d, придатних %3d | K %s Δ%+.3f σ%.3f %s' % (
                os.path.basename(fn), st['name'], len(st['bins']), len(good), st['check_label'], dk, st['check_sd'],
                'OK' if ok else 'ВІДХИЛЕНО (K)'))
            if not ok and not a.force:
                continue
            group += 1
            mode = 'bin %ds' % d['bin_s'] if d['bin_s'] else 'session stack'
            floor = 0.02 if not d['bin_s'] else 0.0       # 1 точка/ніч: систематика ансамблю/кольору TG ≥ 0,02
            for b in good:
                b = dict(b, err=float((b['err'] ** 2 + floor ** 2) ** 0.5))
                notes = ('' if a.live else 'archival ') + ('DWARF Mini 30/150 IMX662 OSC green channel; %d x %gs %s; JD mid-bin; '
                         'ensemble %d %s comps; dark+flat calibrated' % (b['n'], d['exp'], mode, len(st['ensemble']),
                                                                         'APASS DR9' if d['chart'] == 'na' else 'VSP'))
                if filt == 'V':
                    tc = d['tcoef']
                    notes += ('; transformed TG to V: Tv_bv=%+.3f Tbv=%.2f; target B-V %.2f from B channel' %
                              (tc['Tv_bv'], tc['Tbv'], b['bv']))
                assert ',' not in notes
                rows.append(','.join([st['name'], '%.5f' % b['jd'], '%.3f' % b['mag'], '%.3f' % b['err'], filt,
                                      'YES' if filt == 'V' else 'NO', 'STD',
                                      'ENSEMBLE', 'na', st['check_auid'], '%.3f' % b['kmag'], '%.3f' % b['airmass'],
                                      str(group), d['chart'], notes]))
    print('рядків: %d' % len(rows))
    if a.write:
        open(a.write, 'w', newline='\n').write('\n'.join(hdr + rows) + '\n')
        print('записано', a.write)

if __name__ == '__main__':
    main()
