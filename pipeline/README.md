# Photometry pipeline (KTFA, DWARF Mini)

The scripts I use to turn raw DWARF Mini sub-frames into AAVSO observations. They are working scripts,
not a package: paths at the top of each file point to my own folders and have to be edited.
Comments and help strings are in Ukrainian.

I am a web developer (PHP), not a Python programmer. These scripts were generated with Claude (Anthropic)
and I check their output on check stars and against other observers' data in the AAVSO database.

| File | What it does |
|---|---|
| `aavso_phot.py` | Main script (v2). Raw Bayer frame (RGGB, 12 bit) -> dark (matched by exposure, gain, sensor temperature) and flat -> green superpixel (G1+G2)/2, 960x540, 8.03"/px -> alignment to a plate-solved reference (astroalign) -> aperture photometry with SEP -> ensemble of VSP comparison stars with a check star -> transformation TG to V -> per-star JSON and PNG. |
| `aavso_export_series.py` | Writes the AAVSO Extended Format file from the JSON results (OBSTYPE DSLR, transformed V, ensemble, check star in KNAME/KMAG). |
| `transform_coef.py` | Transformation coefficients (Tv_bv, Tbv) from standard or well-covered fields. |
| `make_ref_wcs.py` | Builds the astrometric reference (median of aligned frames + WCS fitted to catalogue positions) when no plate-solved stack is available. |
| `targets_archive.json` | Target list used by `--stars`: name, AUID, coordinates. |

## Typical run

```
python -X utf8 aavso_phot.py --session "<session folder>" --field SZLyn --cloud <reference_wcs.fits> ^
    --stars "SZ Lyn" --bin 180 --vmin 8.3 --with-failed --bkgmax 5
python -X utf8 aavso_export_series.py out_v\SZLyn_<date>.json --live --write KTFA_AAVSO_<date>.txt
```

Notes:

* `--with-failed` also takes frames the on-board stacker renamed to `failed_*`: the raw data in them are intact.
* Frames whose size differs from the master dark are skipped (the scope sometimes saves 1280x720 frames at the start of a session).
* A frame is rejected when fewer than 60 % of the reference stars are matched or the sky background exceeds `--bkgmax` times the session median.
* Saturation guard: a star is dropped when its brightest green pixel exceeds 3500 raw ADU on more than 5 % of the frames.
* Requirements: numpy, astropy, sep, astroalign, matplotlib.

A complete worked example with results is in [`datasets/NGC869_2026-09-03`](../datasets/NGC869_2026-09-03/).
