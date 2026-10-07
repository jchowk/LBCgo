# LBCB gain and read noise, published values (Giallongo et al. 2008)

## Summary

- **Product file:** `LBCgo/conf/lbc_detector.ecsv` (the four `LBCB` rows)
- **Made by / date:** seeded 2026-10-04 on the PI's decision
  (PR jchowk/LBCgo#6)
- **Supersedes:** nothing; before this, the header values
  (`GAIN = 1.75` e⁻/ADU, `RDNOISE = 12` e⁻) were used
- **Valid for:** LBCB, all filters, open-ended date range (to be replaced or
  date-limited once measured values exist)

## Purpose

The `GAIN` and `RDNOISE` header keywords are the same nominal values on every
chip of both LBC cameras. The weight maps (`LBCgo/masks.py`) use per-chip gain
and read noise, so the published LBCB measurements are used until LBCgo's own
measurements replace them. LBCR has no published per-chip values and keeps
using header values.

## Inputs

None run by LBCgo: the values are transcribed from Giallongo et al. 2008,
A&A 482, 349, Table 1 ("variance method" on flat-field sequences, LBCB
commissioning, 2006).

| | chip 1 | chip 2 | chip 3 | chip 4 |
|---|---|---|---|---|
| gain (e⁻/ADU) | 1.96 | 2.09 | 2.06 | 1.98 |
| read noise (e⁻) | 11.4 | 11.6 | 11.6 | 11.2 |

## Next version

Measure per chip and per channel with
`LBCgo.detector.measure_gain_rdnoise_files(flat1, flat2, bias1, bias2)`
(two raw flats at similar, unsaturated levels and two raw biases), at several
epochs. Record those runs in a new directory (e.g. `gain_rdnoise_2027a/`)
following [`../NEW_PRODUCT_TEMPLATE.md`](../NEW_PRODUCT_TEMPLATE.md), compare with the table
above, and replace or date-limit these rows.

**Status (2026-10-07):** measured rows from
[`../gain_rdnoise_lbc_202505/`](../gain_rdnoise_lbc_202505/README.md)
(`mjd_start` = 60822) supersede these from MJD 60822 on. These rows stay
in force for earlier dates and when the date is unknown. The 2025 LBCB
gains are 0.82–0.92 × these values, while read noise in ADU agrees for
chip 2; the difference is not yet explained.
