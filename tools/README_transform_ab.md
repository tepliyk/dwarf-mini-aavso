# Raw Bayer vs bilinear demosaic: does it matter for transformation? (2026-11-07)

Question raised by Roy Axelsen in the AAVSO thread "Transformation Coefficients: Important for New Observers":
is there any measurable gain from doing photometry on the raw Bayer channels instead of interpolated (demosaiced) frames?

`transform_ab.py` measures the same calibrated DWARF Mini sub-frames two ways:

* **A, raw Bayer**: the two G sites of each 2x2 cell averaged into a half-resolution image, B taken from the (1,1) site. Aperture r = 2 superpixels.
* **B, interpolated**: bilinear demosaic (mean of the known neighbours in a 3x3 window) to full 1920x1080 G and B planes, aperture r = 4 px, i.e. the same sky area.

Same darks, flats, stars, alignment, saturation and cloud rejection. Fits as in `transform_coef.py`:
V - g = Tv_bv (B - V) + c, B - V = Tbv (b - g) + c.

## Results

NGC 869 field, 2026-09-03, 45 x 30 s, GAIN 60, 23 VSP stars with B and V (`tcoef_ab_NGC869_20260903.json`):

| | A raw Bayer | B interpolated |
|---|---|---|
| Tv_bv | +0.062 +/- 0.019 | +0.060 +/- 0.019 |
| Tbv | 2.33 +/- 0.11 | 2.25 +/- 0.12 |
| scatter of V - g about the fit | 0.032 | 0.032 |
| median frame-to-frame scatter per star, G | 0.030 | 0.029 |
| same, B | 0.085 | 0.086 |

M35 field, 2026-09-03, 60 x 30 s, 11 stars, low altitude (`tcoef_ab_M35_20260903.json`): Tv_bv +0.000 +/- 0.083 (A) vs +0.051 +/- 0.051 (B), frame-to-frame scatter 0.076 vs 0.078.

Pooled fit, 34 stars: Tv_bv +0.063 +/- 0.019 (A), +0.058 +/- 0.019 (B).
Cross-check (coefficients from one field applied to the other, each star in turn as the target): 0.05 to 0.06 mag either way for both reductions.

Conclusion: on this scope (4 arcsec/px, FWHM about 2 px) interpolation costs nothing measurable. The raw-channel reduction is kept only because the pipeline already works that way.

Full log: `ab_run_2026-11-07.log`.
