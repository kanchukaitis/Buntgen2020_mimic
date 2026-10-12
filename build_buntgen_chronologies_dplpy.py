#!/usr/bin/env python3
"""
build_buntgen_chronologies_dplpy.py
===================================================================================
dplPy reproduction of the Buntgen et al. (2020, Dendrochronologia 64:125757)
regional tree-ring-width (TRW) chronologies: the 16-member standardization
ensemble per site, and the ensemble median used as the site chronology.

Implements the ensemble described in the paper's Supplementary "Chronology
development" paragraph: a Regional Curve Standardization (RCS) ensemble crossing
  * variance stabilization : none (ratios) vs Cook & Peters (1997) power
    transformation (residuals, i.e. subtractive index)
  * regional-curve smoothing: default vs a length-adaptive "-10" spline
  * data pooling           : all series vs a living/relict split, merged
  * chronology type        : Standard (STD) vs ARSTAN (ARS)
plus signal-free RCS (sfptRCS10 and sfRCS10), for 16 members total:

  {all, liv/rel} x [ ptRCS-10_STD, ptRCS-10_ARS, RCS_STD, RCS_ARS,
                     RCS-10_STD, RCS-10_ARS, sfptRCS10, sfRCS10 ]

Per site the 16 members are formed and their median is taken as the chronology.

Method notes (see README for the full discussion):
  * "-10" is a length-adaptive 10%-of-curve-length cubic spline (ARSTAN's
    negative-stiffness = percent-of-length convention), not a fixed 10-year spline.
  * Buntgen built STD/ARS in ARSTAN, whose index calculation places residual
    (power-transform) indices on the same ~1 scale as ratio indices. dplPy's
    rcs(method="difference") returns residuals centered on 0, so the power-transform
    members are rescaled here to the mean & SD of their matching ratio member
    (rescale_to_ratio) -- a chronology-level stand-in for ARSTAN's per-series
    rescaling, so all 16 members share one scale and the median is coherent.
  * Pith offsets are unavailable in the archived measurements, so po = 0 (pith
    present; first measured ring = cambial age 1) for all series -- dplPy 0.8's
    pith offset is a TRUE offset (0 = pith present). This is a known RCS limitation
    for relict wood and most affects the thin early portion of NYA (see README).

Tools: dplr-equivalent routines in dplPy -- readers(), powt(), rcs(), sfrcs()
       (CRUST single-curve signal-free RCS), chron(), chron_ars().

Run:   python3 build_buntgen_chronologies_dplpy.py           # all 9 sites
       python3 build_buntgen_chronologies_dplpy.py QUL NYA   # selected sites
Inputs:  itrdb/<CODE>_all.rwl       Outputs: output/<CODE>_all_variants.csv,
         output/<CODE>_chrono.csv, output/buntgen_all_chronologies.csv,
         output/buntgen_site_metadata.csv
Requires: pip install "dplpy>=0.8.0"  (uses the 0.8 harmonized RCS API). The signal-free members
          iterate, so the largest sites (QUL, NYA) take several minutes each.
"""
import sys, warnings
import numpy as np
import pandas as pd
from pathlib import Path
import dplpy as dpl
import tempfile, os

def read_site(path):
    """Read a site's Tucson .rwl. dplPy reads all nine files directly as
    downloaded. The shipped NYA file has already had its ITRDB series-ID
    formatting cleaned (see fix_nya_ids / README); dplPy reads the raw ITRDB
    NYA equivalently either way, so no per-read normalization is applied here
    (a blanket ID rewrite would wrongly collide IDs in other files, e.g. GTB)."""
    return dpl.readers(path)


def fix_nya_ids(in_path, out_path):
    """Utility (not used by the pipeline): regenerate the cleaned NYA file from
    the raw ITRDB download. The raw NYA archive stores some IDs with an internal
    space in the 8-char ID field (e.g. 'Y8d 024b'); ARSTAN needs the standard
    left-justified form ('Y8d024b'). This removes internal spaces in cols 1-8
    and left-justifies, leaving the year (col 9+) and data untouched. It does
    not change the chronology -- it only cleans the series-ID labels."""
    raw = open(in_path, encoding="latin1").read().splitlines()
    fixed = [(f"{ln[:8].replace(' ', ''):<8}" + ln[8:]) if len(ln) >= 9 else ln
             for ln in raw]
    open(out_path, "w", encoding="latin1").write("\n".join(fixed) + "\n")

# ---- paths (edit if needed) -------------------------------------------------------
HERE     = Path(__file__).resolve().parent
ITRDB_DIR = HERE / "itrdb"                # the nine *_all.rwl inputs
OUT_DIR   = HERE / "output"               # per-site CSVs + combined medians (created at run time)

FILE_MAP = {c: f"{c}_all.rwl" for c in
            ("GTB","SCO","QUL","SCA","ALP","YAM","TAI","ALT","NYA")}
SITES_META = {  # code: (name, lat, lon, season)  -- NSC is the paper code for the SCA file
    "GTB":("Great Basin",39,-117,"7, 8"),"SCO":("S Colorado Plateau",35,-111,"6"),
    "QUL":("Quebec-Labrador",54,-71,"8"),"SCA":("N Scandinavia",68,20,"7, 8"),
    "ALP":("Austrian Alps",47,12,"7, 8"),"YAM":("Yamal",69,71,"7"),
    "TAI":("Taimyr",72,90,"7"),"ALT":("Altai",50,88,"6, 7"),"NYA":("N Yakutia",68,128,"6")}
# NSC is the paper's code for N Scandinavia; the combined medians CSV uses SCA as the
# column name (buildBuntgenMat renames SCA->NSC when assembling buntgen.mat).

# ---- parameters (from the R script / supplement) ----------------------------------
MIN_LEN       = 40     # drop series shorter than this (years)
LR_THRESH     = 1950   # living if last ring > 1950, else relict
MIN_N_MERGE   = 20     # equal-weight merge when both subsets have >= this many series
MIN_N_SINGLE  = 3      # min series to use a value at all
MIN_VARIANTS  = 6      # min valid members for a year to enter the final chronology
MIN_SERIES    = 3      # min total series for a year

# ---- helpers ----------------------------------------------------------------------
def trim_short(r, n=MIN_LEN):
    return r.loc[:, r.notna().sum() >= n]

def make_po(r):
    # .rwl files carry no pith offset -> pith_offset = 0 for all (pith present,
    # first measured ring = cambial age 1). dplPy 0.8's pith offset is a TRUE
    # offset (0 = pith present). KNOWN LIMITATION for relict wood; supply real
    # offsets here if available.
    return pd.DataFrame({"series": list(r.columns),
                         "pith_offset": np.zeros(len(r.columns), dtype=int)})

def split_lr(r, t=LR_THRESH):
    last = r.apply(lambda c: c.last_valid_index())
    return r.loc[:, last > t], r.loc[:, last <= t]

def _std(rwi):
    c = dpl.chron(rwi, biweight=True, plot=False)          # year-indexed; cols std, samp_depth
    return c["std"], c["samp_depth"]

def _ars(rwi):
    c = dpl.chron_ars(rwi, biweight=True, verbose=False)   # year-indexed; std,res,ars,samp_depth
    return c["ars"], c["samp_depth"]                       # "ars" = re-reddened ARSTAN

def _rcs_indices(rwl_in, ratios, smoothing):
    """RCS detrend. smoothing='ten' -> the supplement's '-10' spline (10%-of-
    curve-length; dplPy 0.8 period=-10); smoothing='default' -> the distinct 'RCS
    without the -10 secondary spline' level (fit='AgeDepSpline').
    *** 'default' RC smoothing is the one unconfirmed mapping; its effect on the
    median is tiny (see header). pos_slope=True keeps the pre-0.8 AgeDepSpline
    behaviour (dplPy 0.8 changed its default to pos_slope=False). ***"""
    po = make_po(rwl_in)
    method = "ratio" if ratios else "difference"
    if smoothing == "ten":
        return dpl.rcs(rwl_in, po=po, method=method, fit="Spline", period=-10,
                       plot=False)
    return dpl.rcs(rwl_in, po=po, method=method, fit="AgeDepSpline",
                   pos_slope=True, plot=False)

# the 16 members (name, use_pt, smoothing, chron_kind, is_signal_free)
SPEC = [
    ("ptRCS-10_STD", True,  "ten",     "STD", False),
    ("ptRCS-10_ARS", True,  "ten",     "ARS", False),
    ("RCS_STD",      False, "default", "STD", False),
    ("RCS_ARS",      False, "default", "ARS", False),
    ("RCS-10_STD",   False, "ten",     "STD", False),
    ("RCS-10_ARS",   False, "ten",     "ARS", False),
    ("sfptRCS10",    True,  "ten",     "STD", True),
    ("sfRCS10",      False, "ten",     "STD", True),
]
VAR_NAMES = [f"{p}_{nm}" for p in ("all", "lr") for (nm, _a, _b, _c, _d) in SPEC]
# residual-scale member -> its ratio counterpart (for the mean-offset rescaling)
RATIO_COUNTERPART = {"ptRCS-10_STD": "RCS-10_STD", "ptRCS-10_ARS": "RCS-10_ARS",
                     "sfptRCS10": "sfRCS10"}

def build_variant(rwl, use_pt, smoothing, chron_kind, is_sf, min_len=MIN_LEN):
    """One ensemble member -> (crn_series, n_series_series) or (None, None)."""
    rwl = trim_short(rwl, min_len)
    if rwl.shape[1] < 2:
        return None, None
    try:
        if use_pt:
            rwl_in = dpl.powt(rwl, method="cook")   # Cook & Peters variance stabilization
            ratios = False                           # residuals after power transform
        else:
            rwl_in, ratios = rwl, True               # ratios on raw RW
        if is_sf:
            rwi = dpl.sfrcs(rwl_in, po=make_po(rwl_in),
                            method="ratio" if ratios else "difference",
                            plot=False, verbose=False)   # CRUST single-curve SF-RCS
        else:
            rwi = _rcs_indices(rwl_in, ratios, smoothing)
        crn, nser = _ars(rwi) if chron_kind == "ARS" else _std(rwi)
        return crn.replace([np.inf, -np.inf], np.nan), nser
    except Exception as e:
        print(f"      ERROR [{'pt' if use_pt else 'rw'},{smoothing},{chron_kind},sf={is_sf}]: {e}")
        return None, None

def rescale_to_ratio(diff_crn, ratio_crn):
    """MEAN-OFFSET POST-STEP (see header #1): put a difference (~0) chronology onto
    the mean & SD of its matching ratio (~1) chronology over their common years --
    a chronology-level stand-in for ARSTAN calind's per-series index rescaling."""
    if diff_crn is None or ratio_crn is None:
        return diff_crn
    c = diff_crn.index.intersection(ratio_crn.index)
    d, r = diff_crn.reindex(c), ratio_crn.reindex(c)
    ok = d.notna() & r.notna()
    if ok.sum() < 10:
        return diff_crn
    sd = d[ok].std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return diff_crn
    return (diff_crn - d[ok].mean()) / sd * r[ok].std(ddof=1) + r[ok].mean()

def merge_chrons(cl, nl, cr, nr):
    """Equal-weight merge of living & relict subset chronologies (Buntgen 2011)."""
    if cl is None and cr is None: return None, None
    if cl is None: return cr, nr
    if cr is None: return cl, nl
    yrs = sorted(set(cl.index) | set(cr.index))
    m  = pd.Series(index=yrs, dtype=float)
    ns = pd.Series(0, index=yrs, dtype=float)
    for y in yrs:
        a, b = cl.get(y, np.nan), cr.get(y, np.nan)
        na = nl.get(y, 0) if nl is not None else 0
        nb = nr.get(y, 0) if nr is not None else 0
        if np.isfinite(a) and np.isfinite(b) and na >= MIN_N_MERGE and nb >= MIN_N_MERGE:
            m[y] = 0.5 * (a + b)
        elif np.isfinite(a) and na >= MIN_N_SINGLE: m[y] = a
        elif np.isfinite(b) and nb >= MIN_N_SINGLE: m[y] = b
        else: m[y] = np.nan
        ns[y] = na + nb
    return m, ns

def build_all_members(rwl_all):
    """All 16 members for one site -> dict name -> (crn, n_series)."""
    out = {}
    liv, rel = split_lr(rwl_all)
    has_lr = liv.shape[1] >= 5 and rel.shape[1] >= 5
    for prefix in ("all", "lr"):
        for (nm, pt, sm, ck, sf) in SPEC:
            name = f"{prefix}_{nm}"
            if prefix == "all":
                out[name] = build_variant(rwl_all, pt, sm, ck, sf)
            elif has_lr:
                cl, nl = build_variant(liv, pt, sm, ck, sf)
                cr, nr = build_variant(rel, pt, sm, ck, sf)
                out[name] = merge_chrons(cl, nl, cr, nr)
            else:                                   # no usable split -> reuse all-series
                out[name] = out.get(f"all_{nm}", (None, None))
    # MEAN-OFFSET post-step: rescale residual-scale members onto their ratio counterpart
    for prefix in ("all", "lr"):
        for res_nm, rat_nm in RATIO_COUNTERPART.items():
            k, kr = f"{prefix}_{res_nm}", f"{prefix}_{rat_nm}"
            if out.get(k, (None,))[0] is not None and out.get(kr, (None,))[0] is not None:
                out[k] = (rescale_to_ratio(out[k][0], out[kr][0]), out[k][1])
    return out

def process(code):
    rwl = read_site(str(ITRDB_DIR / FILE_MAP[code]))
    members = build_all_members(rwl)
    yrs = sorted(set().union(*[c.index for c, _ in members.values() if c is not None]))
    mat  = pd.DataFrame(np.nan, index=yrs, columns=VAR_NAMES)
    ndep = pd.DataFrame(0,      index=yrs, columns=VAR_NAMES, dtype=float)
    for nm, (c, n) in members.items():
        if c is None: continue
        cc = c.reindex(yrs)
        nn = n.reindex(yrs) if n is not None else pd.Series(0, index=yrs)
        keep = cc.notna() & (nn.fillna(0) >= MIN_N_SINGLE)
        mat.loc[keep, nm] = cc[keep]; ndep.loc[keep, nm] = nn[keep]
    med   = mat.median(axis=1)
    nvar  = mat.notna().sum(axis=1)
    nser  = ndep.max(axis=1)
    valid = (nvar >= MIN_VARIANTS) & (nser >= MIN_SERIES)
    res = pd.DataFrame({"year": yrs, "crn_median": med.values, "crn_mean": mat.mean(axis=1).values,
                        "crn_min": mat.min(axis=1).values, "crn_max": mat.max(axis=1).values,
                        "n_variants": nvar.values, "n_series": nser.values})[valid.values].reset_index(drop=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mat.to_csv(OUT_DIR / f"{code}_all_variants.csv", index_label="year")
    res.to_csv(OUT_DIR / f"{code}_chrono.csv", index=False)
    nb = sum(v[0] is not None for v in members.values())
    print(f"{code} ({SITES_META[code][0]}): {rwl.shape[1]} series | built {nb}/16 | "
          f"final {len(res)} yrs {int(res.year.min())}..{int(res.year.max())} | "
          f"median mean={np.nanmean(res.crn_median):.3f}")
    return res.set_index("year")["crn_median"].rename(code)

if __name__ == "__main__":
    warnings.simplefilter("ignore")
    codes = sys.argv[1:] or list(FILE_MAP)
    meds = {c: process(c) for c in codes}
    # combined summary (one column per site)
    comb = None
    for c, s in meds.items():
        comb = s.to_frame() if comb is None else comb.join(s, how="outer")
    comb.sort_index().to_csv(OUT_DIR / "buntgen_all_chronologies.csv", index_label="year")
    pd.DataFrame([{"code": c, "name": SITES_META[c][0], "lat": SITES_META[c][1],
                   "lon": SITES_META[c][2], "season": SITES_META[c][3], "type": "TRW"}
                  for c in codes]).to_csv(OUT_DIR / "buntgen_site_metadata.csv", index=False)
    print(f"\nDone. Output in {OUT_DIR}")
