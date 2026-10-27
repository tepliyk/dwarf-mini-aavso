# Test set: NGC 869 field, 2026-09-03 (DWARF Mini, 45 x 30 s)

One clean photometric session from the DWARF Mini tele camera, with everything needed to re-run
differential photometry from the raw Bayer frames and compare against my own result.
Intended as a regression / interoperability set for OSC smart-telescope photometry software
(CAst, SeeVar, ASTAP, your own scripts). Data CC BY 4.0, observer KTFA.

## Observation

| | |
|---|---|
| Field | NGC 869 (h Persei), centre RA 02h 19m, Dec +57.1 deg, FOV 2.1 x 1.2 deg |
| Date / time | 2026-09-03, 01:39 to 02:05 local (UTC+3), i.e. 22:39 to 23:05 UTC on 2026-09-02 |
| Camera | Sony IMX662 OSC, 1920 x 1080, 2.9 um, RGGB, 12-bit data in 16-bit container |
| Optics | 30 mm f/5 (150 mm), 4.0 arcsec/px |
| Exposure | 45 x 30 s, GAIN 60, filter "Astro" (clear), EQ mode, no dithering |
| Sensor temperature | 31 C (from file names and header) |
| Site | Dnipro, Ukraine, city balcony, roughly Bortle 6, airmass 1.08 to 1.09 (field at about 67 deg altitude) |
| Sky | clear and steady (the best of the three NGC 869 nights; 01.09 and 02.09 had haze); waning gibbous Moon up in the east |

## Files

### Release asset `DWARFmini_NGC869_2026-09-03_raw.zip` (137 MB)

```
NGC869_2026-09-03/
  raw/                    45 sub-frames, raw Bayer FITS, untouched headers
  onboard_stack/          the scope's own stack of the same 45 frames (FITS + PNG + JPG) and shotsInfo.json
  calibration/dark/       master darks shipped with the scope: 30 s, GAIN 60, 30 to 32 C (several stacks)
  calibration/bias/       bias for gain index 2
  calibration/flat/       flats for gain index 2, both IR filter positions (ir_1 = Astro)
```

### Release asset `DWARFmini_NGC869_2026-09-03_cloudstack_wcs.fits` (12 MB)

DWARFLAB cloud re-stack of the three NGC 869 nights (5310 s total) with a plate-solved WCS.
My pipeline uses it only as the astrometric reference: frames are aligned to its green channel,
then the WCS offset is refined on the VSP comparison stars. Not needed if your software solves plates itself.

### In this folder

| Path | What |
|---|---|
| `shotsInfo.json` | the scope's session summary (target, RA/Dec, exposure, gain, temperatures) |
| `manifest.json` | which calibration files were packed |
| `results/aavso_extended_NGC869_2026-09-03.txt` | 77 observations as submitted to the AAVSO AID (transformed V, OBSTYPE DSLR) |
| `results/NGC869_202609030137*.json` | full pipeline output: per-star ensemble, check star, per-bin mag/err/airmass, untransformed TG, B-V |
| `results/targets_NGC869.json` | the 21 catalogued variables in the field (name, AUID, position, type, period) |
| `results/comps_NGC869_vsp.json` | VSP comparison-star sequence used (chart X43330GU, V and B) |
| `results/tcoef_NGC869_20260903.json` | per-star instrumental g, b vs catalogue V, B used for the transformation coefficients |
| `plots/*.png` | one light-curve plot per star (mine vs check star) |

## Header notes that matter for timing

* `DATE-OBS` is the **end** of the exposure (the header comment says so), in the **local time** of the paired phone
  (UTC+3 on this date), with no time-zone keyword. Mid-exposure UTC = DATE-OBS - 3 h - 15 s.
* `EXPTIME` 30, `GAIN` 60, `FILTER` 'Astro', `CAMERA` 'TELE'. No `JD`, no `BAYERPAT` keyword: pattern is RGGB,
  first pixel (0,0) is R.
* Pixel values: pedestal about 200 ADU, saturation 4095 (12-bit). The on-board and cloud stacks are scaled to 65535.

## How my numbers were produced

1. Dark (same exposure/gain, nearest temperature) and flat (per Bayer channel) calibration.
2. Green channel only: the two G sites of each 2x2 cell averaged into a 960 x 540 image. Blue channel the same way,
   used only for colour.
3. Source detection and centroids with `sep`, registration with `astroalign` to the cloud stack, WCS offset refined on
   the comparison stars.
4. Aperture photometry, r = 2 superpixels, background from the `sep` background map.
5. Ensemble zero point from the 8 brightest consistent VSP comps (weights 1/sigma^2, median), nearest-brightness comp
   as the check star. Bins of 300 s (10 frames) for the short-period stars, one point per session for the slow ones.
6. Transformation TG to V: V = g + zp + Tv_bv (B-V), with (B-V) = Tbv (b-g) + c from the sensor's own G and B
   channels. Coefficients from 33 VSP stars in this field and M35: Tv_bv = +0.048 +/- 0.014, Tbv = 1.94 +/- 0.09.
7. Error = statistical error combined with the scatter of the ensemble zero point.

Check star in this session: label 137 (AUID 000-BBD-554), catalogue V 13.749, measured 13.734 +/- 0.066 over the session.

## Things to compare

* Zero point and check-star residual against mine (`check_mean`, `check_sd` in the result json).
* GSC 03693-01705 (HADS(B), P = 0.0911 d, 13.4 mag): 5 bins across half a pulsation cycle, a good test of
  precision at the faint end.
* Mis V1395 (EW, 11.4 mag) and NSVS 1908107 (EA, 10.9 mag): 9 bins each, flat parts of the light curve, a good test
  of systematics.
* Blue-channel scatter: compare your blue values with `b` in `tcoef_NGC869_20260903.json`.

Questions and results: the thread "DWARF Mini (30 mm) photometry" in the Smart Telescope Photometry category
on forums.aavso.org, or an issue in this repository.
