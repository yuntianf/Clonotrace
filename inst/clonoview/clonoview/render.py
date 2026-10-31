"""Server-side rasters for the cell panel.

WHY THE RIGHT PANEL IS A PICTURE AND THE LEFT ONE IS NOT. The clone embedding has about
1,500 points, so it is sent to the browser as coordinates and drawn there -- every point
stays individually hoverable and clickable. The cell embedding has 150,150 (Mathew) to
273,508 (Wang); shipping that as JSON would be tens of megabytes and would not redraw
smoothly in Safari over an SSH tunnel. So the cell panel's *background* is rendered here
once per (colour mode, selection), cached on disk, and reused.

Everything that changes as the user interacts -- a clone's cells, its density contours,
its arrows, the direct labels -- is drawn in the browser on top of that backdrop. Which
is what makes scrubbing a clone through time a client-side redraw with no round trip.

THE EXTENT IS PINNED. Each raster is drawn to exactly `meta.extent.cell`, with the axes
filling the whole figure and no margins, so browser pixel (0,0) is the top-left of that
extent and nothing has to be negotiated between the two sides. If this file ever stops
using the stored extent, every overlay silently shifts.

Labels are NOT baked in. They are drawn client-side from centroids the API returns, so
they can be toggled, stay crisp at any zoom, and never end up half-rendered into a
cached PNG.
"""

import hashlib
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt      # noqa: E402

from . import colour as C            # noqa: E402

# Point size falls with cell count so a dense map does not become a solid slab and a
# sparse one does not vanish. The floor matters more than it looks: markers below about
# 2 pt^2 at this dpi are sub-pixel, and antialiasing then spreads each one into a pale
# tint of itself, so a fully saturated palette renders as pastel and the structure
# washes out. Compared at 0.45 / 1.5 / 3.0 / 5.0 on Wang's 273,508 cells -- 3.0 is where
# the colours hold. Cost is flat across that range (0.26-0.37 s), so there is nothing to
# buy by going smaller.
def _point_size(n, width):
    scale = (width / 1400.0) ** 2
    if n > 200_000:
        return 3.0 * scale
    if n > 100_000:
        return 4.0 * scale
    return 6.0 * scale


def _fig(extent, width):
    x0, x1, y0, y1 = extent
    aspect = (y1 - y0) / (x1 - x0)
    height = max(1, int(round(width * aspect)))
    dpi = 100.0
    fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_axis_off()
    ax.set_facecolor(C.SURFACE)
    fig.patch.set_facecolor(C.SURFACE)
    return fig, ax, height


def cache_key(**kw):
    s = "|".join(f"{k}={kw[k]}" for k in sorted(kw))
    return hashlib.sha1(s.encode()).hexdigest()[:16]


def background(b, mode="celltype", select=None, profile=None, mana=False, width=1400,
               cache_dir=None):
    """Render the cell panel backdrop.

    mode "celltype"  colour by the dataset's own cell annotation
    mode "profile"   colour by the cell's Clonotrace profile (argmax of its profile
                     probability vector -- see the confidence caveat in the README)
    mode "plain"     no colour at all; the map as pure context

    COLOURING AND FILTERING ARE SEPARATE, deliberately. `select` isolates one class of
    whatever `mode` colours by. `profile` restricts which cells are drawn at all,
    whatever the colouring. Keeping them apart is what allows the question the two
    panels exist to answer: "show me profile 7's cells, coloured by cell type" -- which
    is how you find out what a profile is made of. Tying the profile choice to the
    colour mode would make that view unreachable.

    Isolating is also the reliable way to answer "where is EM.GZMK" on Mathew, whose 16
    cell types cannot be separated by colour alone at any assignment, so there it is the
    primary tool rather than a convenience.
    """
    key = cache_key(name=b.name, mode=mode, select=select, profile=profile,
                    mana=mana, width=width)
    path = os.path.join(cache_dir, f"bg_{key}.png") if cache_dir else None
    if path and os.path.exists(path):
        return path

    extent = b.meta["extent"]["cell"]
    fig, ax, height = _fig(extent, width)
    xy = b.cxy
    s = _point_size(len(xy), width)

    if mode == "profile":
        values, colours = b.c_profile, b.profile_colours
        classes = sorted(colours)
    elif mode == "celltype":
        values = b.c_type
        pal = b.cell_type_colours
        colours = {i: pal[i] for i in range(len(b.cell_types))}
        classes = sorted(colours)
    else:
        values, colours, classes = None, {}, []

    if values is None:
        ax.scatter(xy[:, 0], xy[:, 1], s=s, c=C.CONTEXT, linewidths=0, marker=".")
    else:
        # Context first, coloured on top, so an assigned cell is never buried under an
        # unassigned one.
        if select is None:
            ctx = values < 0 if mode == "celltype" else values <= 0
        else:
            ctx = values != select
        if profile is not None:
            # Cells outside the chosen profile stay on the map as context rather than
            # being deleted from it. Dropping them would rescale nothing (the extent is
            # pinned) but would remove the landscape the selection is meant to be read
            # against, which is the entire point of showing a profile on the cell map.
            ctx = ctx | (b.c_profile != profile)
        ax.scatter(xy[ctx, 0], xy[ctx, 1], s=s, c=C.CONTEXT, linewidths=0, marker=".")

        # The coloured cells go down in ONE scatter, in a shuffled order, rather than
        # one class at a time. Drawing class by class puts whichever class is drawn last
        # on top of every overlap, so in a mixed region the map shows the highest-
        # numbered class present rather than what is actually there -- a systematic bias
        # that reads as biology. A fixed permutation makes occlusion unbiased: each
        # class covers a share of the overlap proportional to how many cells it has
        # there. Seeded, so the same bundle always renders identically.
        # Exactly the complement of the context mask. Deriving it independently is how
        # the profile filter got silently dropped once: `ctx` had the filter ORed in
        # but the selection recomputed itself from `values` alone, so asking for one
        # profile returned the whole map coloured. There is only one mask.
        idx = np.flatnonzero(~ctx)
        if len(idx):
            rng = np.random.default_rng(0)
            rng.shuffle(idx)
            cols = np.array([colours[int(k)] for k in classes])
            lut = {k: i for i, k in enumerate(classes)}
            ci = np.array([lut[int(v)] for v in values[idx]])
            # Under a filter the selection is figure and the rest is ground, so the
            # selected cells are drawn larger. Without this a diffuse profile -- one
            # whose clones spread across the whole map, which is a real and common kind
            # rather than a defect -- lays a thin scatter over 273k context points and
            # is close to invisible at the same size.
            sel_s = s * 2.6 if (select is not None or profile is not None) else s
            ax.scatter(xy[idx, 0], xy[idx, 1], s=sel_s, c=cols[ci], linewidths=0,
                       marker=".")

    # MANAscore last, so it is never buried. Red, as asked, and drawn LARGER WITH A
    # WHITE HALO rather than relying on the colour alone: the categorical palettes
    # already contain reds (Liu's first cell type is #94180d), so on a dense map hue
    # alone would not reliably separate a called cell from an ordinary one. The halo
    # makes the mark legible over any underlying colour, which is what an overlay has to
    # do -- it sits on top of whatever the panel is already coloured by.
    if mana and b.has_manascore:
        m = b.c_mana == 1
        if m.any():
            ax.scatter(xy[m, 0], xy[m, 1], s=s * 3.4, c=C.MANA_RED, linewidths=s * 0.16,
                       edgecolors="#ffffff", marker="o", zorder=5)

    if path:
        os.makedirs(cache_dir, exist_ok=True)
        fig.savefig(path, dpi=100, facecolor=C.SURFACE)
        plt.close(fig)
        return path
    plt.close(fig)
    return None


def clone_background(b, width=900, cache_dir=None, coloured=True):
    """The clone panel is drawn in the browser, but a rasterised copy is kept for export
    and as a fallback. Small enough that it costs nothing."""
    key = cache_key(name=b.name, panel="clone", coloured=coloured, width=width)
    path = os.path.join(cache_dir, f"cl_{key}.png") if cache_dir else None
    if path and os.path.exists(path):
        return path
    fig, ax, _ = _fig(b.meta["extent"]["clone"], width)
    if coloured:
        cols = b.profile_colours
        for k in sorted(cols):
            m = b.u_profile == k
            ax.scatter(b.uxy[m, 0], b.uxy[m, 1], s=14, c=cols[k], linewidths=0)
    else:
        ax.scatter(b.uxy[:, 0], b.uxy[:, 1], s=14, c=C.CONTEXT, linewidths=0)
    os.makedirs(cache_dir, exist_ok=True)
    fig.savefig(path, dpi=100, facecolor=C.SURFACE)
    plt.close(fig)
    return path
