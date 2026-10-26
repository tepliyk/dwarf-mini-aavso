DWARF Mini twilight series for system gain measurement (AAVSO STWG request, forum post #114 by M. Munkacsy)
Observer: Kostiantyn Tepliy (KTFA), Dnipro, Ukraine (48.45 N, 35.0 E)

Telescope: DWARF Mini, firmware 1.1.4.7, tele camera (30/150 mm, Sony IMX662 OSC, Bayer RGGB)
Date: 2026-10-26, evening twilight. Sunset 16:35 local (UTC+2), civil dusk 17:04.
Pointing: blank sky, east, altitude ~40-60 deg, no tracking target (planetary RAW mode, mount static).
Settings: GAIN 0 (lowest available), exposure 0.25 s, filter "Astro" (UV/IR-cut), RAW FITS, 12-bit data in 16-bit container.
Series: 13 folders, 16:30 to 17:02 local, 60 frames each (first folder 5 frames), ~3 min apart, 725 frames total.
No dark/flat/bias applied. Files named "failed_*" are flagged by the scope only because no GOTO target was set; they are ordinary raw frames.

Sky level (median ADU, R/G/B, pedestal ~200, saturation 4095):
16:30 1903/3222/2470 -> 16:37 993/1551/1209 -> 16:46 412/548/465 -> 16:55 233/258/247 -> 17:02 209/216/213.
No saturated pixels in any frame.

Header notes: DATE-OBS = END of exposure, local time (phone time zone, UTC+2 on this date), no TIMEZONE keyword.
EXPTIME in seconds. GAIN is the scope's gain setting (not e-/ADU).

My own photon-transfer estimate (ptc_gain.py, pairs of consecutive frames, 25x25 tiles, robust variance, fit below 1500 ADU):
GAIN 0: R 9.5, G 10.0, B 9.4 e-/ADU (+-0.05).  Read noise ~1.3-1.6 ADU (rough).
Dawn 2026-10-07, GAIN 60 (separate archive): R 1.16, G 1.23, B 1.20 e-/ADU.
Both consistent with IMX662 full well ~38 ke- at 12 bits. Channels differ by <7%.
