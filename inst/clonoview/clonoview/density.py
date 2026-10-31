"""Density contours for a set of cells, returned as polygons in data coordinates.

WHY CONTOURS ARE COMPUTED HERE AND DRAWN THERE. The contour of one clone at one visit is
a few hundred numbers. Sending those to the browser, rather than a picture, is what lets
the viewer step or animate through a clone's timepoints without a server round trip per
frame: every frame of a clone arrives in a single request and scrubbing is a redraw.

WHAT THE LEVELS MEAN. Each contour encloses a stated share of THAT clone's cells -- the
0.5 ring contains half of them. This is the highest-density-region construction: evaluate
the kernel density at the cells themselves, and take the level at the corresponding
quantile of those values. It is chosen over evenly spaced density levels because the
number that comes out is one a reader can check ("half its cells are in here") rather
than an unlabelled gradient, and because it stays comparable between a clone with 40
cells and one with 800.

WHEN THERE ARE TOO FEW CELLS, THERE IS NO CONTOUR. The median clone unit in these
datasets has 18 cells and the minimum is 10. A kernel density estimate from 18 points
describes the bandwidth at least as much as the data: the ring would be an artefact with
a number attached, which is worse than no ring at all. Below `MIN_FOR_CONTOUR` this
returns no contours and says why, and the viewer draws the cells alone and labels the
panel as showing points only. The cells themselves are always drawn either way, at every
size -- they are the observation; the contour is a model of it.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt      # noqa: E402  (geometry only; nothing is rendered)
from scipy.ndimage import gaussian_filter   # noqa: E402

# Below this many cells, a contour is a statement about the smoothing kernel rather than
# about the clone. Set against the data: Wang's clone units run 10 to 825 cells with a
# median of 18, so a large share of single clones fall under it and are shown as points.
# That is the honest outcome, not a gap to be filled in.
MIN_FOR_CONTOUR = 30

# Fractions of the clone's own cells that each ring encloses, outermost first.
LEVELS = (0.9, 0.5, 0.25)

GRID = 140


def contours(xy, extent, levels=LEVELS, grid=GRID, min_cells=MIN_FOR_CONTOUR):
    """Density rings for `xy` (n, 2), as polygons in data coordinates.

    Returns {"levels": [...], "rings": [[polygon, ...], ...], "n": n, "reason": str|None}
    where each polygon is a flat [x0, y0, x1, y1, ...] list, and `rings[i]` are the
    polygons for `levels[i]`.
    """
    n = len(xy)
    out = dict(n=int(n), levels=[], rings=[], reason=None)
    if n < min_cells:
        out["reason"] = (f"{n} cells — below the {min_cells}-cell minimum for a density "
                         f"estimate; cells shown individually")
        return out

    # BINNED, NOT scipy.stats.gaussian_kde. The exact estimator costs O(cells x grid
    # points) inside a Python loop -- measured at 5.8 SECONDS for one 825-cell clone on a
    # 140x140 grid, which would make stepping through a clone's timepoints unusable and
    # is hopeless for the pooled multi-clone view where the cell count runs to tens of
    # thousands. Binning the cells and smoothing the histogram with a Gaussian filter is
    # the standard approximation: identical in the limit of a fine grid, and O(grid)
    # regardless of how many cells there are. Same bandwidth rule (Scott's), so the rings
    # are the same rings.
    sd = xy.std(axis=0, ddof=1) if n > 1 else np.zeros(2)
    if not np.all(np.isfinite(sd)) or np.all(sd <= 0):
        out["reason"] = "cells occupy no area; shown individually"
        return out
    # Scott's rule in 2 dimensions: n**(-1/(d+4)) = n**(-1/6).
    bw = np.where(sd > 0, sd, np.max(sd)) * float(n) ** (-1.0 / 6.0)

    # The grid spans the cells, not the whole embedding: a clone occupies a small part of
    # the map, and gridding the entire extent would spend almost every cell of it on
    # empty space and leave the clone itself a few pixels wide. Padded by 3 bandwidths so
    # the outer ring closes instead of being clipped at the edge.
    x0 = max(extent[0], xy[:, 0].min() - 3 * bw[0])
    x1 = min(extent[1], xy[:, 0].max() + 3 * bw[0])
    y0 = max(extent[2], xy[:, 1].min() - 3 * bw[1])
    y1 = min(extent[3], xy[:, 1].max() + 3 * bw[1])
    if not (x1 > x0 and y1 > y0):
        out["reason"] = "cells occupy no area; shown individually"
        return out

    # THE GRID IS SIZED TO THE BANDWIDTH, not fixed. A Gaussian filter costs time in
    # proportion to sigma measured in grid cells, and sigma is bandwidth over grid
    # spacing -- so a compact clone, whose grid spans little more than a few bandwidths,
    # ends up with an enormous kernel. Measured on a fixed 140-cell grid: a tight
    # 825-cell clone took 1.3 s while a diffuse 2,061-cell one took 11 ms, entirely
    # because of this. Solving for a constant sigma instead makes the cost flat and the
    # amount of smoothing per bandwidth identical for every clone, which is also what
    # makes two clones' rings comparable.
    target_sigma = 8.0
    nx = int(np.clip(round(target_sigma * (x1 - x0) / bw[0]), 48, grid))
    ny = int(np.clip(round(target_sigma * (y1 - y0) / bw[1]), 48, grid))

    gx = np.linspace(x0, x1, nx)
    gy = np.linspace(y0, y1, ny)
    H, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=[ny, nx],
                             range=[[y0, y1], [x0, x1]])
    sigma = (max(0.5, bw[1] / (gy[1] - gy[0])), max(0.5, bw[0] / (gx[1] - gx[0])))
    Z = gaussian_filter(H, sigma=sigma, mode="constant")
    if Z.max() <= 0:
        out["reason"] = "density is flat at this resolution; cells shown individually"
        return out
    X, Y = np.meshgrid(gx, gy)

    # Highest-density regions: the level enclosing fraction f of the cells is the
    # (1-f) quantile of the density AT the cells, read off the same grid the contours
    # come from, so the level and the ring cannot disagree.
    ix = np.clip(np.searchsorted(gx, xy[:, 0]) - 1, 0, nx - 1)
    iy = np.clip(np.searchsorted(gy, xy[:, 1]) - 1, 0, ny - 1)
    at_points = Z[iy, ix]
    want = []
    for f in levels:
        lv = float(np.quantile(at_points, 1.0 - f))
        if lv > Z.min():
            want.append((f, lv))
    if not want:
        out["reason"] = "density is flat at this resolution; cells shown individually"
        return out

    # matplotlib's contour generator, used for its geometry only -- nothing is rendered.
    # `allsegs` is the stable accessor; ContourSet.collections is deprecated and gone in
    # newer matplotlib, so it is not used here. Imported at module scope, not per call:
    # the first import cost ~970 ms and showed up as a slow first contour.
    fig = plt.figure(figsize=(1, 1))
    ax = fig.add_subplot()
    try:
        order = np.argsort([lv for _, lv in want])
        lvls = [want[i][1] for i in order]
        fracs = [want[i][0] for i in order]
        # Contour requires strictly increasing levels; identical ones collapse.
        keep = [0] + [i for i in range(1, len(lvls)) if lvls[i] > lvls[i - 1]]
        lvls = [lvls[i] for i in keep]
        fracs = [fracs[i] for i in keep]
        cs = ax.contour(X, Y, Z, levels=lvls)
        segs = cs.allsegs
    finally:
        plt.close(fig)

    for f, polys in zip(fracs, segs):
        rings = []
        for p in polys:
            if len(p) >= 3:
                rings.append([round(float(v), 4) for v in np.asarray(p).ravel()])
        out["levels"].append(round(float(f), 3))
        out["rings"].append(rings)
    return out
