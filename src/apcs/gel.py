"""Colony-PCR gel QC.

Two jobs, both of which happen before anything is pooled:

  1. Throw out wells that have no band, the wrong band, or more than one band.
     A well that fails here costs nothing; the same well sequenced costs a share
     of the submission and a share of every other well's read depth.
  2. Decide how much of each well to pool. Equal volumes give unequal read
     depth, which is the single commonest reason a well lands below min_depth.
     Band intensity is a good enough proxy for concentration to fix that.

Sizing uses a log(size) vs migration fit through the ladder lane.
"""
from __future__ import annotations
import numpy as np
from pathlib import Path
from dataclasses import dataclass, asdict

try:
    from PIL import Image
except ImportError:  # pragma: no cover
    Image = None


def load_image(path) -> np.ndarray:
    """Grayscale float array, 0..1, bands bright."""
    path = str(path)
    if path.lower().endswith((".tif", ".tiff")):
        try:
            import tifffile
            a = tifffile.imread(path)
        except ImportError:
            a = np.asarray(Image.open(path))
    else:
        if Image is None:
            raise RuntimeError("pillow is required: pip install pillow")
        a = np.asarray(Image.open(path))
    scale = np.iinfo(a.dtype).max if np.issubdtype(a.dtype, np.integer) else 1.0
    a = a.astype(float) / scale
    if a.ndim == 3:
        a = a[..., :3].mean(axis=2)
    if a.ndim != 2 or min(a.shape) < 30 or not np.isfinite(a).all() or a.max() > 1 or a.min() < 0:
        raise ValueError("gel must be a 2-D image (or RGB) with finite values in its native range")
    # Some scanners write dark bands on a light background.
    if np.median(a) > 0.5:
        a = 1.0 - a
    return a


def _smooth(x, w):
    if w < 3:
        return x
    k = np.ones(w) / w
    return np.convolve(x, k, mode="same")


def _peaks(x, min_dist, min_prom):
    out = []
    for i in range(1, len(x) - 1):
        if x[i] >= x[i - 1] and x[i] > x[i + 1]:
            lo = max(0, i - min_dist); hi = min(len(x), i + min_dist)
            prom = x[i] - min(x[lo:hi].min(), x[i])
            base = min(x[lo:i].min() if i > lo else x[i], x[i + 1:hi].min() if hi > i + 1 else x[i])
            if x[i] - base >= min_prom:
                out.append((i, x[i], x[i] - base))
    out.sort(key=lambda t: -t[2])
    keep = []
    for i, v, p in out:
        if all(abs(i - j) >= min_dist for j, _, _ in keep):
            keep.append((i, v, p))
    return sorted(keep)


def autocrop(img, pad=4):
    """Trim the dark surround so the gel fills the frame. Returns (img, (r0, c0))."""
    thr = np.percentile(img, 60) * 0.5 + img.max() * 0.02
    rows = np.where(img.max(axis=1) > thr)[0]
    cols = np.where(img.max(axis=0) > thr)[0]
    if len(rows) < 10 or len(cols) < 10:
        return img, (0, 0)
    c0, c1 = max(0, cols[0] - pad), min(img.shape[1], cols[-1] + pad)
    # Trimming to the last bright band would remove the bottom of the
    # migration region and invalidate detect_bands' vertical search range.
    return img[:, c0:c1], (0, c0)


def _grid_from_profile(prof, w):
    """Fit a regular comb to a 1-D column profile: returns (pitch, offset)."""
    p = prof - prof.mean()
    n = len(p)
    ac = np.correlate(p, p, mode="full")[n - 1:]
    lo, hi = max(8, int(w * 0.012)), int(w * 0.12)
    hi = min(hi, len(ac) - 1)
    if hi <= lo:
        return None, None
    pitch = int(lo + np.argmax(ac[lo:hi]))
    best, best_off = -1e18, 0
    for off in range(pitch):
        idx = np.arange(off, n, pitch)
        v = prof[idx].sum()
        if v > best:
            best, best_off = v, off
    return pitch, best_off


def find_lanes(img, n_lanes=None, min_width_frac=0.018, band_rows=(0.12, 0.65),
               well_rows=(0.0, 0.14)):
    """Lane windows as (left, right) column indices, one per comb position.

    The comb wells are visible whether or not a lane was loaded, so they give a
    regular grid that survives empty lanes, smeared lanes and uneven loading.
    Falls back to peak-picking on the migration region if no comb is visible.
    """
    h, w = img.shape
    wr = img[int(h * well_rows[0]):max(2, int(h * well_rows[1]))].mean(axis=0)
    wr = _smooth(wr - np.percentile(wr, 20), max(3, int(w * 0.003)))
    wr[wr < 0] = 0
    pitch, off = _grid_from_profile(wr, w)
    r0, r1 = int(h * band_rows[0]), int(h * band_rows[1])
    mig = _smooth(img[r0:r1].mean(axis=0), max(3, int(w * 0.006)))
    mig = mig - np.percentile(mig, 20)
    mig[mig < 0] = 0
    if pitch and pitch > 6:
        allc = list(range(off % pitch, w, pitch))
        thr_w, thr_m = wr.max() * 0.05, mig.max() * 0.03
        centres = [c for c in allc if wr[c] > thr_w or mig[c] > thr_m]
        if n_lanes and len(centres) > n_lanes:
            centres = sorted(sorted(centres, key=lambda c: -(wr[c] + mig[c]))[:n_lanes])
        if len(centres) >= 3:
            half = max(3, int(pitch * 0.38))
            return [(max(0, c - half), min(w, c + half)) for c in centres]
    pk = _peaks(mig, max(3, int(w * min_width_frac)), mig.max() * 0.15)
    if n_lanes:
        pk = sorted(sorted(pk, key=lambda t: -t[2])[:n_lanes])
    centres = [p[0] for p in pk]
    if len(centres) < 2:
        return []
    half = int(np.median(np.diff(centres)) * 0.35)
    return [(max(0, c - half), min(w, c + half)) for c in centres]


@dataclass
class Band:
    row: int
    intensity: float
    prominence: float
    size_bp: float | None = None
    saturated: bool = False


def lane_profile(img, lane, smooth_frac=0.004):
    l, r = lane
    prof = img[:, l:r].mean(axis=1)
    return _smooth(prof, max(3, int(img.shape[0] * smooth_frac)))


def detect_bands(img, lane, max_bands=4, prom_frac=0.10, band_rows=(0.08, 0.80)):
    """Peaks in one lane, restricted to the migration region so the loading well
    at the top of the lane and the plate edge are not called as bands.

    Rows returned are absolute image rows."""
    h = img.shape[0]
    r0, r1 = int(h * band_rows[0]), int(h * band_rows[1])
    prof = lane_profile(img, lane)[r0:r1]
    p = prof - np.percentile(prof, 20)
    p[p < 0] = 0
    noise = np.median(np.abs(np.diff(prof))) / (0.6745 * 2 ** 0.5)
    if p.max() <= max(img.max() * 0.01, noise * 6):
        return []
    pk = _peaks(p, max(4, int(h * 0.012)), p.max() * prom_frac)
    pk = sorted(pk, key=lambda t: -t[2])[:max_bands]
    l, r = lane
    bands = []
    for i, val, prom in sorted(pk):
        lo_i, hi_i = max(0, i - 10), min(len(p), i + 10)
        row = i + r0
        patch = img[max(0, row - 10):min(h, row + 10), l:r]
        bands.append(Band(row=int(row),
                          intensity=float(p[lo_i:hi_i].sum()),
                          prominence=float(prom),
                          saturated=bool((patch > 0.985).mean() > 0.02)))
    return bands


def calibrate(ladder_rows, ladder_sizes, deg=1):
    """log10(size) ~ poly(row). Returns a callable row -> bp."""
    rows = np.asarray(ladder_rows, float)
    bp = np.asarray(ladder_sizes, float)
    if (len(rows) != len(bp) or len(rows) < 3 or not np.isfinite(rows).all()
            or not np.isfinite(bp).all() or np.any(bp <= 0)
            or np.any(np.diff(rows) <= 0) or np.any(np.diff(bp) >= 0)):
        raise ValueError("ladder requires >=3 ordered band rows and matching descending positive sizes")
    sizes = np.log10(bp)
    deg = min(deg, len(rows) - 1)
    c = np.polyfit(rows, sizes, deg)
    return lambda row: float(10 ** np.polyval(c, row))


def analyse(path, expected_bp=None, n_lanes=None, ladder_lane=None, ladder_sizes=None,
            tolerance=0.20, sample_labels=None, secondary_frac=0.50, prom_frac=0.15,
            crop=True):
    """Full gel pass. Returns per-lane dicts with size, intensity, verdict and pool volume."""
    img = load_image(path)
    if (expected_bp is not None and expected_bp <= 0) or not 0 < tolerance < 1:
        raise ValueError("expected size must be positive; tolerance must be between 0 and 1")
    if not 0 < secondary_frac <= 1 or not 0 < prom_frac <= 1 or (n_lanes is not None and n_lanes < 2):
        raise ValueError("invalid band thresholds or lane count")
    if bool(ladder_sizes) != (ladder_lane is not None):
        raise ValueError("provide both --ladder-lane and --ladder-sizes")
    if crop:
        img, _ = autocrop(img)
    lanes = find_lanes(img, n_lanes)
    if not lanes:
        raise RuntimeError("no lanes detected - pass --lanes N or crop the image to the gel")
    if n_lanes is not None and len(lanes) != n_lanes:
        raise ValueError(f"detected {len(lanes)} lanes, expected {n_lanes}; crop/re-image and verify lane positions")
    if ladder_lane is not None and not 0 <= ladder_lane < len(lanes):
        raise ValueError("ladder lane is outside the detected gel")
    if sample_labels is not None and (len(sample_labels) != len(lanes) - (ladder_lane is not None)
                                      or len(set(sample_labels)) != len(sample_labels)):
        raise ValueError("provide one unique label for every sample lane, in image order")
    sizer = None
    if ladder_lane is not None and ladder_sizes:
        lb = detect_bands(img, lanes[ladder_lane], max_bands=len(ladder_sizes) + 2,
                          prom_frac=0.04)
        lb.sort(key=lambda b: b.row)
        sizes = sorted(ladder_sizes, reverse=True)      # big fragments migrate least
        if len(lb) != len(sizes):
            raise ValueError(f"detected {len(lb)} ladder bands for {len(sizes)} sizes; cannot pair them reliably")
        sizer = calibrate([b.row for b in lb], sizes)

    out = []
    sample_i = 0
    for i, lane in enumerate(lanes):
        if ladder_lane is not None and i == ladder_lane:
            out.append({"lane": i, "label": "ladder", "role": "ladder", "bands": [], "verdict": "-"})
            continue
        bands = detect_bands(img, lane, prom_frac=prom_frac)
        for b in bands:
            if sizer:
                b.size_bp = round(sizer(b.row))
        label = sample_labels[sample_i] if sample_labels is not None else f"lane{i}"
        sample_i += 1
        main = max(bands, key=lambda b: b.intensity) if bands else None
        verdict, note = "NO BAND", "nothing above background"
        if main:
            extra = [b for b in bands if b is not main
                     and b.intensity > secondary_frac * main.intensity]
            if expected_bp and main.size_bp:
                off = abs(main.size_bp - expected_bp) / expected_bp
                if off > tolerance:
                    verdict, note = "WRONG SIZE", f"{main.size_bp:.0f} bp vs {expected_bp} bp expected"
                elif extra:
                    verdict, note = "MULTIPLE BANDS", f"{len(extra)+1} bands - do not pool"
                else:
                    verdict, note = "OK", f"{main.size_bp:.0f} bp"
            elif extra:
                verdict, note = "MULTIPLE BANDS", f"{len(extra)+1} bands - do not pool"
            else:
                verdict = "REVIEW SIZE"
                note = f"{main.size_bp:.0f} bp" if main.size_bp else "single band, size not calibrated"
            if main.saturated:
                verdict = "REVIEW CLIPPING"
                note += "; possible clipping at image range limit; re-image before pooling"
            if sizer and not min(b.row for b in lb) <= main.row <= max(b.row for b in lb):
                verdict, note = "REVIEW SIZE", "band outside ladder calibration range"
        out.append({"lane": i, "label": label, "role": "sample",
                    "lane_bounds_px": [int(x) for x in lane],
                    "bands": [asdict(b) for b in bands],
                    "main_bp": main.size_bp if main else None,
                    "intensity": main.intensity if main else 0.0,
                    "saturated": bool(main.saturated) if main else False,
                    "verdict": verdict, "note": note})

    for r in out:
        r.update(pool_uL=0.0, pool_capped=False)

    ok = [r for r in out if r["verdict"] == "OK" and r["intensity"] > 0]
    if ok:
        ref = max(r["intensity"] for r in ok)
        for r in out:
            if r["verdict"] == "OK" and r["intensity"] > 0:
                v = 5.0 * ref / r["intensity"]
                r["pool_uL"] = round(min(v, 20.0), 1)
                r["pool_capped"] = v > 20.0
                if v > 20.0:
                    r.update(verdict="REVIEW VOLUME", pool_uL=0.0,
                             note="normalisation exceeds 20 uL cap; quantify or re-amplify")
            else:
                r["pool_uL"] = 0.0
                r["pool_capped"] = False
    return {"image": Path(path).name, "n_lanes": len(lanes), "calibrated": sizer is not None,
            "review_required": True,
            "note": "Heuristic image QC; verify all detected lanes and ladder assignments before pooling.",
            "lanes": out}
