"""Build the offline export: a folder that runs from a double-click, and its zip.

WHAT IS AND IS NOT COPIED. app.js, table.js, app.css and link.js are copied from
`static/` at build time and never edited in `export/`. They are the same files the server
serves, and that is the point: the served viewer is the reference the offline one is
checked against, so they have to be one front end with one seam in it (`api()`), not two
that have to be kept in step by hand. Anything under `export/` is a build artefact and is
safe to delete. The three files that exist only offline -- api.js, density.js, raster.js
-- live in `export/app/` in the source tree because they have no served counterpart.

THE HTML IS REWRITTEN, NOT FORKED, for the same reason. `static/index.html` refers to
`/static/app.css` and friends by absolute path, which resolves to the filesystem root when
the page is opened from a folder; the rewrite makes those relative and inserts the extra
script tags the offline build needs. Two hand-maintained copies of a 116-line document
would drift on the first change to either.
"""

import json
import os
import re
import shutil
import textwrap
import zipfile

import numpy as np

from . import bundle as bd
from . import manifest as mf
from . import pack

HERE = os.path.dirname(os.path.abspath(__file__))
VIZ_ROOT = os.path.dirname(HERE)
STATIC = os.path.join(VIZ_ROOT, "static")
# The hand-written half of the offline front end. It used to live in `export/app/`, the
# same directory the build wrote into, so source and build artefact were interleaved and
# the build could only run in a writable source tree. Shipped inside the package, it is
# copied out to the destination like everything else, and `export/` becomes purely
# output -- which is what lets the destination be anywhere, including a read-only
# install serving a folder the user named.
OFFLINE = os.path.join(VIZ_ROOT, "offline")

# Copied from static/ on every build. Sole source of truth is static/.
SHARED = ["app.css", "app.js", "table.js", "link.js", "markers.js"]

# Copied from offline/ into the export's app/ directory. These have no served
# counterpart: they are what replaces the server when the page is opened from a folder.
OFFLINE_ONLY = ["data.js", "density.js", "api.js", "raster.js"]

# Copied from offline/ to the top of the export. The self check is shipped WITH the
# export on purpose: it is the only way the person receiving it can find out whether
# their browser will run the viewer, without having to describe a symptom to somebody who
# cannot see their screen.
EXPORT_PAGES = ["selftest.html", "selftest_peer.html"]

# Order matters: data.js defines the loader the router calls, and the catalogue has to be
# in place before app.js boots and asks for the dataset list.
SCRIPTS = ["app/data.js", "app/density.js", "app/api.js", "app/raster.js",
           "data/index.js", "app/link.js", "app/markers.js"]


def _rewrite_html(src, out, entry):
    with open(src) as fh:
        html = fh.read()
    html = html.replace('href="/static/app.css"', 'href="app/app.css"')
    tags = "\n".join(f'<script src="{s}"></script>' for s in SCRIPTS)
    # The page's own script goes last, after the offline machinery it depends on.
    # The shared scripts are re-emitted under app/ by the SCRIPTS list above, so their
    # served-path tags are stripped here. An absolute /static/... path in the export
    # would resolve to the filesystem root and silently load nothing.
    html = re.sub(r'<script src="/static/link\.js"></script>\s*\n?', "", html)
    html = re.sub(r'<script src="/static/markers\.js"></script>\s*\n?', "", html)
    html = html.replace(f'<script src="/static/{entry}"></script>',
                        f'{tags}\n<script src="app/{entry}"></script>')
    if "app/" + entry not in html:
        raise SystemExit(f"{src}: could not find the <script> tag for {entry} to rewrite")
    with open(out, "w") as fh:
        fh.write(html)


README = """ClonoView — an interactive viewer for the Clonotrace results
===========================================================

TO OPEN IT
----------
Double-click  ClonoView.html

That is the whole of it. There is nothing to install, no Python, no R, and it does not
use the internet -- everything it shows is in the "data" folder beside it. It works in
Chrome, Edge, Firefox and Safari.

If it opens as a page of text instead of a picture, your computer opened the file in a
text editor. Right-click it, choose "Open with", and pick a web browser.

IF ANYTHING LOOKS WRONG, OPEN  selftest.html  FIRST
---------------------------------------------------
It is in the same folder, and it checks -- in your browser, on your machine -- that the
data files load, that the cell map can be drawn, and that the two tabs can talk to each
other. Each line says what to do if it fails. It takes a few seconds and it will usually
tell you the answer faster than describing the symptom to someone else would.

FIRST THING TO KNOW
-------------------
Keep the three things together. ClonoView.html needs the "app" and "data" folders that
came with it, sitting in the same place. Moving the HTML file on its own to the desktop
will give you a page that loads and then says it cannot find its data.

WHAT YOU ARE LOOKING AT
-----------------------
Two panels, side by side, always in the same places:

  LEFT   the clone embedding. One dot is one CLONE UNIT -- a single clone, in a single
         patient, at a single timepoint, that had at least the run's minimum number of
         cells. Two dots for the same clone at two timepoints are two dots, and the
         arrows between them are that clone moving.

  RIGHT  the cell embedding. One dot is one cell. This is the ordinary expression map:
         it is where the cell types live.

Pick a dataset from the menu at the top left. The two "table" buttons open a second tab
-- put it beside this one, because clicking a row there changes the view here, and the
point is to watch it change.

Click a clone (in the table, or directly on the left panel) and the time controls appear
along the bottom: step through the timepoints, or press Play. The right panel then shows
that clone's own cells, and the rings around them enclose 90%, 50% and 25% of them.

Click a patient row instead and you get all of that patient's clones at once.

Clicking an entry in the legend shows only that class and greys out the rest. That is
the reliable way to find one cell type on a crowded map, and on a dataset with more than
about a dozen cell types it is the intended way, because no set of that many colours can
be told apart reliably. The viewer says so under the legend rather than pretending
otherwise.

THINGS THAT LOOK LIKE BUGS AND ARE NOT
--------------------------------------
* A clone can have cells at a timepoint and no dot on the left panel. That is a real
  distinction, not a gap: Clonotrace only places a clone on the clone map at timepoints
  where it passed the minimum cell count. The cells are still there and still drawn; the
  timebar marks the frame as having no clone-map position. How common this is differs a
  lot between studies -- it is most clones in some and a minority in others -- so it is
  given per dataset at the end of this file rather than described here.

* Small clones get no rings. Below thirty cells a density estimate describes the
  smoothing more than it describes the clone, so the viewer draws the cells alone and
  says why.

* A dataset can have one timepoint. Where a study is a single sample per subject there
  is no movement through time to show, and the time controls correctly have one step.

WHERE THIS CAME FROM
--------------------
Built from a Clonotrace run. Every number shown is read from that run's outputs -- the
viewer computes layout and density, and nothing else.
"""


def _clones_below_threshold(data_root, name):
    """How many clones have at least one timepoint with cells but no clone-map position.

    Counted here rather than read from meta.json because meta.json does not carry it:
    `counts` is entirely cell-level apart from the totals. It is derived rather than
    stored on purpose -- everything it needs is already in the bundle, so storing it
    would be a second copy of a fact that could go stale against the arrays it came from.

    A cell qualifies when it is on a tracked clone (`lineage >= 0`) whose frame at that
    cell's own timepoint is not a clone unit (`clone_unit < 0`), which is exactly the
    situation the README describes. The clones are then the distinct lineages over those
    cells.
    """
    with np.load(os.path.join(data_root, name, "cells.npz")) as c:
        lineage, unit = c["lineage"], c["clone_unit"]
    return int(np.unique(lineage[(lineage >= 0) & (unit < 0)]).size)


def _dataset_notes(data_root, names):
    """The 'what is in this export' section, read from the bundles themselves.

    This used to be three paragraphs written by hand about the four datasets that
    happened to be prepared at the time, which is wrong in a shipped tool: an export
    built from somebody else's run would have described datasets it did not contain.
    Everything below comes out of each bundle's own meta.json, so the README describes
    what was actually packed -- including the study's own `notes`, which is where a
    caveat like "these are numbered biopsies, not successive visits" belongs.
    """
    out = ["\n\nWHAT IS IN THIS EXPORT\n----------------------\n"]
    for name in names:
        with open(os.path.join(data_root, name, "meta.json")) as fh:
            d = json.load(fh)
        n_vis = len(d.get("visits") or [])
        when = ("one timepoint" if n_vis <= 1 else
                f"{n_vis} timepoints" +
                (", which are not a time course" if d.get("kind") == "pre_post" else ""))
        c = d.get("counts") or {}

        def n(k, word):
            v = c.get(k, 0)
            return f"{v:,} {word}{'' if v == 1 else 's'}"

        body = [f"{n('cells', 'cell')}, {n('clone_units', 'clone unit')} over "
                f"{n('lineages', 'clone')}, {n('patients', 'subject')}, {when}."]
        # The figure behind the "cells at a timepoint with no dot" note above, and it has
        # to be counted over CLONES, because that is what the note is about. The obvious
        # substitute -- counts.cells_below_threshold, which is already in meta.json --
        # answers a different question and cannot stand in for this one: across the four
        # datasets developed against, the cell count ranges over a factor of twelve while
        # the share of clones affected ranges over a factor of two and a half, so no cell
        # count lets a reader tell 31% of clones from 74%. That distinction is the whole
        # point of the note, and getting it wrong is not academic: the sentence this
        # replaces said "most clones", which was true of three of those datasets and
        # false of the fourth at 31%.
        below = c.get("cells_below_threshold", 0)
        n_aff, n_lin = _clones_below_threshold(data_root, name), c.get("lineages", 0)
        if n_aff and n_lin:
            body.append(f"{n_aff:,} of its {n_lin:,} clones ({100 * n_aff / n_lin:.0f}%) "
                        f"have at least one such timepoint, accounting for {below:,} of "
                        f"the cells drawn with no clone-map position. The run's "
                        f"threshold was {d.get('expanded_min')} cells.")
        note = (d.get("notes") or "").strip()
        if note:
            body.append(note)
        out.append(f"{d.get('label', d['name'])} ({d['name']})\n")
        # Every paragraph wrapped to the same width as the hand-written text above it,
        # so a long study name or a four-timepoint summary does not run off the page in
        # a file people open in Notepad.
        for para in body:
            out.append("  " + "\n  ".join(textwrap.wrap(para, 76)) + "\n")
        out.append("\n")
    return "".join(out)


def build(names=None, data_root=None, out=None, make_zip=True, zip_path=None):
    data_root = data_root or mf.data_root()
    out = out or os.path.join(os.path.dirname(data_root), "export")
    names = names or [d["name"] for d in bd.available(data_root)]
    if not names:
        raise SystemExit(f"no prepared datasets in {data_root} to export.")

    app_dir = os.path.join(out, "app")
    os.makedirs(app_dir, exist_ok=True)

    for f in OFFLINE_ONLY:
        src = os.path.join(OFFLINE, f)
        if not os.path.exists(src):
            raise SystemExit(f"{src} is missing -- it is source, not a build artefact, "
                             f"and cannot be regenerated by this script.")
        shutil.copy2(src, os.path.join(app_dir, f))
    for f in EXPORT_PAGES:
        src = os.path.join(OFFLINE, f)
        if not os.path.exists(src):
            raise SystemExit(f"{src} is missing -- it is source, not a build artefact, "
                             f"and cannot be regenerated by this script.")
        shutil.copy2(src, os.path.join(out, f))
    for f in SHARED:
        shutil.copy2(os.path.join(STATIC, f), os.path.join(app_dir, f))

    _rewrite_html(os.path.join(STATIC, "index.html"),
                  os.path.join(out, "ClonoView.html"), "app.js")
    _rewrite_html(os.path.join(STATIC, "table.html"),
                  os.path.join(out, "table.html"), "table.js")

    print("packing datasets:")
    data_out = os.path.join(out, "data")
    for n in names:
        r = pack.pack_dataset(data_root, n, data_out)
        print(f"  {r['name']:8s} {r['cells']:>8,} cells  {r['bytes'] / 1e6:6.2f} MB")
    pack.pack_index(data_root, data_out, set(names))

    with open(os.path.join(out, "README.txt"), "w") as fh:
        fh.write(README + _dataset_notes(data_root, names))

    if not make_zip:
        return None

    zpath = zip_path or os.path.join(os.path.dirname(out), "ClonoView-offline.zip")
    os.makedirs(os.path.dirname(os.path.abspath(zpath)), exist_ok=True)
    # Everything is nested under one folder inside the zip. Unzipping a flat archive
    # scatters ClonoView.html, app/ and data/ into whatever folder the person happened to
    # be in, and the three have to stay together.
    root = "ClonoView"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for base, _dirs, files in os.walk(out):
            for f in sorted(files):
                p = os.path.join(base, f)
                rel = os.path.relpath(p, out)
                z.write(p, os.path.join(root, rel))
    print(f"\n  {zpath}  {os.path.getsize(zpath) / 1e6:.1f} MB")
    return zpath


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Build the offline ClonoView export.")
    ap.add_argument("names", nargs="*", help="datasets to include; default all prepared")
    ap.add_argument("--data", default=None, help="directory of prepared bundles")
    ap.add_argument("--out", default=None, help="directory to build the export in")
    ap.add_argument("--zip", dest="zip_path", default=None, help="path for the zip")
    ap.add_argument("--no-zip", action="store_true")
    a = ap.parse_args(argv)
    build(names=a.names or None, data_root=a.data, out=a.out,
          make_zip=not a.no_zip, zip_path=a.zip_path)


if __name__ == "__main__":
    main()
