# Büntgen et al. (2020) chronology reproduction

A reproducible reconstruction of the 9 regional tree-ring-width (TRW)
ensemble-median chronologies of

> Büntgen, U. et al. (2020) *Prominent role of global volcanism in Common Era
> climate variability and human history.* **Dendrochronologia** 64, 125757.

The original study built, for each site, an ensemble of 16 different
Regional Curve Standardization (RCS) chronologies and took their median as
the site chronology. Neither the code nor the final chronologies were archived,
so this repository reconstructs them from the publicly archived ring-width
measurements using [**dplPy**](https://github.com/OpenDendro/dplPy) (the Python
Dendrochronology Program Library).

## Sites and data

Raw ring-width files are in [`itrdb/`](itrdb/), one per site (`<CODE>_all.rwl`,
Tucson format). Eight are the current ITRDB versions; **YAM** is the file used in
the related Büntgen et al. (2021) exercise (the ITRDB holding differs) and is
included here for completeness.

| code | region | series | span (CE) | ensemble-median level |
|------|--------|-------:|-----------|-----------------------|
| GTB | Great Basin, USA | 721 | −2649 to 2009 | 0.996 |
| SCO | S Colorado Plateau, USA | 224 | −568 to 2002 | 0.971 |
| QUL | Quebec–Labrador, Canada | 2725 | 5 to 2016 | 1.058 |
| SCA | N Scandinavia (paper code NSC) | 1237 | −181 to 2010 | 0.960 |
| ALP | Austrian Alps | 1886 | −717 to 2016 | 0.917 |
| YAM | Yamal, Russia | 665 | −585 to 2015 | 0.970 |
| TAI | Taimyr, Russia | 356 | −207 to 2011 | 1.013 |
| ALT | Altai, Russia | 738 | −87 to 2017 | 1.049 |
| NYA | N Yakutia, Russia | 1883 | −379 to 2016 | 1.207 |

The original ITRDB measurement files are available from the NOAA/WDS
International Tree-Ring Data Bank: https://www.ncei.noaa.gov/products/paleoclimatology/tree-ring

**Note on NYA.** The raw ITRDB N Yakutia file stores some series IDs with an
internal space in the 8-character ID field (e.g. `Y8d 024b`). The version in
`itrdb/` has been reformatted to the standard left-justified ID (`Y8d024b`),
leaving the year and measurement columns untouched; this does not change the
chronology, only the series labels. (dplPy reads either form equivalently. This
cleanup is for compatibility with stricter readers. The
`fix_nya_ids()` helper in the script regenerates it from a fresh download.)

## The 16-member ensemble

Per site, the ensemble crosses standardization choices described in the paper's
supplement, plus signal-free RCS, for 16 members — run for **all** series and for
a **living/relict** split (last ring > 1950 = living) merged at equal weight:

| member | variance stabilization | index | RC smoothing | chronology |
|--------|------------------------|-------|--------------|------------|
| ptRCS-10_STD / _ARS | power transform (Cook & Peters 1997) | residual | −10 | STD / ARSTAN |
| RCS_STD / _ARS | none | ratio | default | STD / ARSTAN |
| RCS-10_STD / _ARS | none | ratio | −10 | STD / ARSTAN |
| sfptRCS10 | power transform | residual | −10 | signal-free |
| sfRCS10 | none | ratio | −10 | signal-free |

The ensemble **median** across the (up to) 16 members is the site chronology.

Method notes:

- **`−10`** indicates a *length-adaptive* cubic spline at **10 % of the regional-curve
  length** (ARSTAN's negative-stiffness = percent-of-length convention), not a
  fixed 10-year spline.
- **Common scale for the median.** Büntgen's STD/ARSTAN chronologies came from
  ARSTAN, which places power-transform (residual) indices on the same ~1 scale as
  ratio indices. dplPy's residual detrending is centered on 0, so the
  power-transform members here are rescaled to the mean and standard deviation of
  their matching ratio member before the median is taken, so all members share one
  scale. (This is applied as a post-processing step; see the script header.)
- **Signal-free RCS** uses dplPy's `sfrcs()`, a port of CRUST single-curve
  signal-free RCS (Melvin & Briffa 2014).
- Robust (biweight) means throughout, matching the original.

## Caveats

- **Pith offsets are unavailable** in the archived measurements, so a pith offset
  of 1 (first measured ring = cambial age 1) is assumed for every series. This is
  a well-known limitation of RCS for relict/sub-fossil wood, where the true
  cambial age at the first ring can be large. It is most visible at **NYA**, whose
  early record (1–1000 CE) is almost entirely relict and whose
  ensemble median is correspondingly elevated (~1.2). Supplying real pith offsets
  (via the `pith_offset` column) would most improve that site.
- **Some series counts may differ from the paper's Table S1** presumably because the
  online ITRDB holdings have grown/changed since 2020. Using the current archive
  is the reproducible choice, but may not fully replicate Buntgen et al. 2020.
- Some archived files contain minor data-quality quirks (e.g. a few overlapping
  duplicate IDs at GTB). dplPy handles these on read.

## Reproduce

```bash
pip install dplpy==0.7.0
python3 build_buntgen_chronologies_dplpy.py          # all nine sites
python3 build_buntgen_chronologies_dplpy.py QUL NYA  # selected sites
```

The signal-free members iterate, so the largest sites (QUL, NYA) take a few
minutes each.

## Outputs (`output/`)

- `<CODE>_all_variants.csv` — the 16 member chronologies for a site (year index).
- `<CODE>_chrono.csv` — the final ensemble chronology for a site: `crn_median`,
  `crn_mean`, `crn_min`, `crn_max`, `n_variants`, `n_series`, filtered to years
  with ≥ 6 members and ≥ 3 series.
- `buntgen_all_chronologies.csv` — the nine ensemble medians, one column per site.
- `buntgen_site_metadata.csv` — site code, name, lat, lon, calibration season, type.

## Acknowledgements

Chronologies built with [dplPy](https://github.com/OpenDendro/dplPy). Original
measurements courtesy of the contributing investigators and the ITRDB. See the
paper's Supplementary Table S1 for the source study of each regional dataset.
