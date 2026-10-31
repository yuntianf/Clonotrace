"""The ClonoView HTTP server -- Python standard library only.

NO FRAMEWORK ON PURPOSE. The cluster's analysis environment has numpy, pandas, scipy and
matplotlib but no Flask, FastAPI, Dash, Bokeh or Plotly, and installing into it is both
slow and a good way to break the environment other work depends on. `http.server` plus a
dictionary of routes is a few hundred lines and has no such cost, so the viewer runs
anywhere Python does.

WHAT GOES OVER THE WIRE. Small and structural things go as JSON: the clone embedding
(about 1,500 points), the patient table, one clone's cells. The cell embedding's
backdrop goes as a cached PNG, because 150k-273k points is too much JSON to redraw
smoothly over an SSH tunnel. See render.py for why that split is drawn where it is.

Bind address defaults to 127.0.0.1. This node has no public route, so reaching the
viewer means an SSH tunnel, and binding to localhost keeps it off the cluster's internal
network rather than quietly serving to every host in the VPC.
"""

import argparse
import gzip
import hmac
import json
import mimetypes
import os
import posixpath
import secrets
import socket
import sys
import threading
import time
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import urllib.error
import urllib.request
from urllib.parse import urlparse, parse_qs

import numpy as np

from . import bundle as bd
from . import colour as C
from . import manifest as mf
from . import density as dens
from . import render as rd

HERE = os.path.dirname(os.path.abspath(__file__))
VIZ_ROOT = os.path.dirname(HERE)
# Read-only, and shipped with the code: `static/` sits beside this package inside the
# installed R package, so it is found relative to the module.
STATIC = os.path.join(VIZ_ROOT, "static")
# Written to, so it must NOT be. An installed package directory is read-only on a managed
# R installation and shared between users on a cluster, so bundles and the pidfile live in
# the user's cache -- see manifest.data_root(). Both are still overridable: `--data` for
# the bundles, CLONOVIEW_DATA for both.
DATA = mf.data_root()
RUN_DIR = os.path.dirname(DATA)

# matplotlib's pyplot state is global, so rasters are rendered one at a time even though
# the server is threaded. They are cached on disk, so this is contended only on the
# first request for each (mode, selection).
RENDER_LOCK = threading.Lock()


def _np_json(o):
    """json.dumps refuses numpy scalars, and every route builds its payload out of
    numpy arrays. Converting at the source is still the rule -- this is the backstop so
    one missed int() returns a number instead of a 500 from deep inside the encoder."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(f"not JSON serialisable: {type(o).__name__}")


class Api:
    """Route handlers. Each returns a JSON-serialisable object, or ("file", path)."""

    def __init__(self, data_root):
        self.data_root = data_root

    def cache_dir(self, name):
        return os.path.join(self.data_root, name, "cache")

    # ---- catalogue -----------------------------------------------------------
    def datasets(self, q):
        return dict(datasets=bd.available(self.data_root))

    def dataset(self, q, name):
        b = bd.load(self.data_root, name)
        pal = b.cell_type_colours
        nprof = len(b.profile_ids)
        return dict(
            meta=b.meta,
            patients=b.names["patients"],
            cell_types=b.cell_types,
            profile_colours={str(k): v for k, v in b.profile_colours.items()},
            cell_type_colours={str(i): pal[i] for i in range(len(pal))},
            # The viewer states what the palette actually measures rather than implying
            # a separation it does not have. At Mathew's 16 cell types no assignment
            # clears the gates, so the UI leans on labels and isolate instead.
            cell_type_palette_quality=C.cell_type_palette_quality(len(b.cell_types)),
            profile_palette_quality=dict(
                (lambda m: dict(worst_cvd=round(m["worst_cvd"], 2),
                                worst_normal=round(m["worst_normal"], 2),
                                low_contrast=m["low_contrast"]))(
                    C.measure(C.PROFILE_COLOURS[:nprof]))),
            colours=dict(surface=C.SURFACE, ink=C.INK, ink2=C.INK_2, muted=C.MUTED,
                         context=C.CONTEXT, grid=C.GRID, time=C.SEQ_TIME,
                         mana=C.MANA_RED),
            manascore=b.meta.get("manascore"),
            axes=b.axes(),
            centroids=dict(
                profile=b.centroids(b.c_profile, b.profile_ids),
                cell_type=b.centroids(b.c_type, range(len(b.cell_types))),
            ),
        )

    def centroids(self, q, name):
        """Direct-label positions, optionally for one profile's cells only.

        These have to be recomputed under a filter rather than reused. With the map
        restricted to profile 7, a cell-type label parked at that type's centroid over
        ALL cells can easily sit where profile 7 has no cells at all -- the label would
        then point at empty space and name something that is not there. Profiles are
        distributions over the whole expression map, so this is the common case, not an
        edge case.
        """
        b = bd.load(self.data_root, name)
        prof = (q.get("profile") or [None])[0]
        prof = int(prof) if prof not in (None, "", "null") else None
        keep = None if prof is None else (b.c_profile == prof)

        def cent(values, labels):
            out = {}
            for lab in labels:
                m = values == lab
                if keep is not None:
                    m = m & keep
                if int(m.sum()) == 0:
                    continue
                out[int(lab)] = [float(np.median(b.cxy[m, 0])),
                                 float(np.median(b.cxy[m, 1])), int(m.sum())]
            return out

        return dict(profile=cent(b.c_profile, b.profile_ids),
                    cell_type=cent(b.c_type, range(len(b.cell_types))),
                    filtered_to_profile=prof)

    # ---- the clone panel -----------------------------------------------------
    def clones(self, q, name):
        b = bd.load(self.data_root, name)
        # Which clone layout to draw. The clone identities, profiles and every count are
        # identical across layouts -- only the coordinates and the extent change -- so a
        # layout switch is a redraw of the same points, never a different selection.
        layout = (q.get("layout") or ["default"])[0]
        uxy = b.uxy_for(layout)
        return dict(
            layout=layout,
            x=[round(float(v), 4) for v in uxy[:, 0]],
            y=[round(float(v), 4) for v in uxy[:, 1]],
            profile=[int(v) for v in b.u_profile],
            patient=[int(v) for v in b.u_patient],
            visit=[int(v) for v in b.u_visit],
            lineage=[int(v) for v in b.u_lineage],
            n_cells=[int(v) for v in b.u_n],
            # -1 means this dataset has no MANAscore, or this clone was never called.
            ptrc=[int(v) for v in b.u_ptrc],
            mana_usable=[int(v) for v in b.u_mana_usable],
            extent=b.clone_extent(layout),
        )

    # ---- the clone table -----------------------------------------------------
    def lineages(self, q, name):
        """Rows for the clone table, optionally restricted to some patients.

        Each row is one lineage: its label, its total cells, and which visits it is a
        clone unit at. `visits_with_cells` includes visits where it fell below the run's
        threshold and therefore has no clone-embedding position -- those are shown
        differently rather than hidden, because the cells are real.
        """
        b = bd.load(self.data_root, name)
        want = q.get("patient")
        pats = set(int(v) for v in want) if want else None
        rows = []
        for i, l in enumerate(b.names["lineages"]):
            if pats is not None and l["patient"] not in pats:
                continue
            units = np.flatnonzero(b.u_lineage == i)
            a, z = b.cells_of_lineage(i)
            vis = b.c_visit[a:z]
            rows.append(dict(
                lineage=i, label=l["label"], patient=l["patient"],
                receptor=l["receptor"], n_cells=l["n_cells"],
                n_units=int(len(units)),
                visits=[int(v) for v in sorted(b.u_visit[units])],
                profiles=[int(p) for p in b.u_profile[units[np.argsort(b.u_visit[units])]]],
                visits_with_cells=[int(v) for v in sorted(set(vis[vis >= 0].tolist()))],
            ))
        rows.sort(key=lambda r: (-r["n_units"], -r["n_cells"]))
        return dict(rows=rows, n=len(rows))

    def lineage(self, q, name, idx):
        """One lineage: its frames through time, and the cells behind each frame."""
        b = bd.load(self.data_root, name)
        i = int(idx)
        l = b.names["lineages"][i]
        layout = (q.get("layout") or ["default"])[0]
        axis = (q.get("axis") or ["time"])[0]
        fixed = int((q.get("fixed") or ["-1"])[0])
        if axis == "tissue" and fixed < 0:
            # Default to the visit this lineage actually has most of its cells at, so
            # "show me this clone by site" works with one click and lands somewhere with
            # data rather than on an empty first visit.
            aa, zz = b.cells_of_lineage(i)
            vv = b.c_visit[aa:zz]
            fixed = int(np.bincount(vv[vv >= 0], minlength=len(b.visits)).argmax()) \
                if int((vv >= 0).sum()) else 0
        frames = b.lineage_frames(i, layout=layout, axis=axis, fixed=fixed)
        a, z = b.cells_of_lineage(i)
        for f in frames:
            if f["cell_slice"]:
                s, e = f["cell_slice"]
                f["cells"] = dict(
                    x=[round(float(v), 4) for v in b.cxy[s:e, 0]],
                    y=[round(float(v), 4) for v in b.cxy[s:e, 1]],
                    cell_type=[int(v) for v in b.c_type[s:e]],
                    profile=[int(v) for v in b.c_profile[s:e]],
                )
            else:
                f["cells"] = dict(x=[], y=[], cell_type=[], profile=[])
            # Contours for every frame in this one response. Each is 11-17 ms and a
            # clone has at most a handful of visits, so the whole time course costs well
            # under a tenth of a second -- which is what lets the viewer step and animate
            # through time as a redraw, with no request per frame.
            #
            # FROM THE STORED COORDINATES, NOT THE ROUNDED ONES. This used to read the
            # x and y it had just written into f["cells"], which are rounded to four
            # decimals -- a size optimisation for the JSON going over the tunnel. Feeding
            # them back into the density estimate let that optimisation reach the
            # statistics: on a small clone the rounding moves cells between histogram
            # bins, which moves the level read off the density at the cells, and the 50%
            # ring came out about 0.8% different in area on a 38-cell frame of Mathew's
            # lineage 226. It also disagreed with the pooled multi-clone view four
            # methods below, which always contoured the stored coordinates -- so the same
            # clone could be drawn with two different rings depending on which view it was
            # reached from. Found by comparing this against the JavaScript port, which had
            # no reason to round and therefore did not.
            if f["cell_slice"]:
                s, e = f["cell_slice"]
                f["density"] = dens.contours(b.cxy[s:e], b.meta["extent"]["cell"])
            else:
                f["density"] = dict(n=0, levels=[], rings=[], reason="no cells")
            f.pop("cell_slice", None)
        return dict(lineage=i, label=l["label"], patient=l["patient"],
                    patient_id=b.names["patients"][l["patient"]],
                    receptor=l["receptor"], n_cells=l["n_cells"], frames=frames,
                    axis=axis, fixed=fixed,
                    min_cells_for_contour=dens.MIN_FOR_CONTOUR)

    # ---- a set of clones, through time --------------------------------------
    # Cap on cells shipped per frame. A whole patient at one visit runs to ~15,000 cells
    # and a multi-patient selection to far more; past this many the browser is drawing
    # more marks than the panel has pixels, so the cost buys nothing. Above the cap the
    # frame is subsampled and says by how much -- the density contours are always
    # computed from ALL the cells, before any subsampling, so the shape being read is
    # never the shape of a sample.
    MAX_CELLS_PER_FRAME = 12000

    def cset(self, q, name):
        """Every clone of a chosen set of patients, frame by frame through the schedule.

        Returns one frame per visit: where those clones sit on the clone map, the arrows
        from their positions at the previous visit they held one, all of their cells,
        and the density those cells form. This is the multi-clone counterpart of
        /api/lineage -- and like it, the whole time course comes back in ONE response so
        stepping through time is a redraw rather than a request.
        """
        b = bd.load(self.data_root, name)
        layout = (q.get("layout") or ["default"])[0]
        uxy = b.uxy_for(layout)
        axis = (q.get("axis") or ["time"])[0]
        fixed = int((q.get("fixed") or ["-1"])[0])
        pats = [int(v) for v in (q.get("patient") or [])]
        lins = set(int(v) for v in (q.get("lineage") or []))
        # Captured BEFORE the patient expansion below widens `lins`, because these are
        # the hand-picked clones and that distinction is what the markers encode.
        sel_lins = sorted(lins)
        if pats:
            lins |= set(int(i) for i, l in enumerate(b.names["lineages"])
                        if l["patient"] in pats)
        if not lins:
            return dict(patients=pats, lineages=[], frames=[], n_lineages=0,
                        layout=layout, axis=axis, fixed=fixed,
                        note="no clones selected")
        lin_arr = np.fromiter(sorted(lins), dtype=np.int64)

        in_set_units = np.isin(b.u_lineage, lin_arr)
        in_set_cells = np.isin(b.c_lineage, lin_arr)

        # Where each lineage sits at each visit, so arrows can span a visit the lineage
        # missed rather than being dropped at the gap.
        pos = {}
        for k in np.flatnonzero(in_set_units):
            pos.setdefault(int(b.u_lineage[k]), {})[int(b.u_visit[k])] = k

        # Which axis the frames run along. See Bundle.lineage_frames for why the two are
        # not symmetric: clone units exist per visit, never per tissue, so stepping
        # through anatomic positions holds the clone map still (same units, same arrows,
        # all belonging to the one fixed visit) and varies only the cells.
        if axis == "tissue" and fixed < 0:
            vv = b.c_visit[in_set_cells]
            fixed = int(np.bincount(vv[vv >= 0], minlength=len(b.visits)).argmax()) \
                if int((vv >= 0).sum()) else 0
        if axis == "tissue":
            steps = [(fixed, ti) for ti in range(len(b.tissues))]
        else:
            ti = fixed if 0 <= fixed < len(b.tissues) else -1
            steps = [(vi, ti) for vi in range(len(b.visits))]

        rng = np.random.default_rng(0)
        frames = []
        for vi, ti in steps:
            um = in_set_units & (b.u_visit == vi)
            uidx = np.flatnonzero(um)
            cm = in_set_cells & (b.c_visit == vi)
            if ti >= 0:
                cm = cm & (b.c_tissue == ti)
            cidx = np.flatnonzero(cm)
            # Never sampled at this site: skip rather than draw an empty panel that
            # reads as "no cells were found here".
            if axis == "tissue" and len(cidx) == 0:
                continue

            # Contours from every cell, before any subsampling.
            dmap = (dens.contours(b.cxy[cidx], b.meta["extent"]["cell"])
                    if len(cidx) else dict(n=0, levels=[], rings=[], reason="no cells"))

            shown, sub = cidx, False
            if len(cidx) > self.MAX_CELLS_PER_FRAME:
                shown = rng.choice(cidx, self.MAX_CELLS_PER_FRAME, replace=False)
                shown.sort()
                sub = True

            arrows = []
            for li, byv in (pos.items() if axis != "tissue" else []):
                if vi not in byv:
                    continue
                prev = [v for v in byv if v < vi]
                if not prev:
                    continue
                pv = max(prev)
                a, z = pos[li][pv], byv[vi]
                arrows.append(dict(
                    lineage=li,
                    x1=round(float(uxy[a, 0]), 4), y1=round(float(uxy[a, 1]), 4),
                    x2=round(float(uxy[z, 0]), 4), y2=round(float(uxy[z, 1]), 4),
                    from_visit=pv, p1=int(b.u_profile[a]), p2=int(b.u_profile[z]),
                    # True when the lineage held no position at the visit(s) between --
                    # the arrow spans a gap rather than a single step, and the viewer
                    # draws it dashed so the two are never confused.
                    skipped=vi - pv > 1,
                ))

            frames.append(dict(
                visit=vi,
                code=(b.tissues[ti] if axis == "tissue" else b.visits[vi]),
                label=(b.tissues[ti] if axis == "tissue"
                       else b.meta["visits"][vi]["label"]),
                **(dict(tissue=ti, tissue_label=b.tissues[ti]) if ti >= 0 else {}),
                n_units=int(len(uidx)), n_cells=int(len(cidx)),
                n_lineages=int(len(np.unique(b.u_lineage[uidx]))),
                units=dict(
                    x=[round(float(v), 4) for v in uxy[uidx, 0]],
                    y=[round(float(v), 4) for v in uxy[uidx, 1]],
                    profile=[int(v) for v in b.u_profile[uidx]],
                    lineage=[int(v) for v in b.u_lineage[uidx]],
                    n_cells=[int(v) for v in b.u_n[uidx]],
                    ptrc=[int(v) for v in b.u_ptrc[uidx]],
                ),
                cells=dict(
                    x=[round(float(v), 2) for v in b.cxy[shown, 0]],
                    y=[round(float(v), 2) for v in b.cxy[shown, 1]],
                ),
                cells_shown=int(len(shown)), subsampled=bool(sub),
                density=dmap,
                arrows=arrows,
            ))

        return dict(
            patients=pats,
            patient_ids=[b.names["patients"][p] for p in pats],
            # The clones the user picked BY HAND, in a stable order, with their names.
            # These are the ones the viewer gives a distinct marker to, so the order
            # must not depend on request order -- a reload would otherwise reshuffle
            # which clone is a square and which a triangle. Clones pulled in by a
            # patient selection are deliberately NOT here: there can be hundreds, and
            # a per-clone marker for hundreds of clones says nothing.
            lineages_selected=sel_lins,
            lineages_selected_labels=[b.names["lineages"][i]["label"]
                                      for i in sel_lins],
            axis=axis, fixed=fixed,
            n_lineages=len(lins),
            n_units=int(in_set_units.sum()),
            n_cells=int(in_set_cells.sum()),
            frames=frames,
            min_cells_for_contour=dens.MIN_FOR_CONTOUR,
            max_cells_per_frame=self.MAX_CELLS_PER_FRAME,
        )

    # ---- rasters -------------------------------------------------------------
    def bg(self, q, name):
        b = bd.load(self.data_root, name)
        mode = (q.get("mode") or ["celltype"])[0]
        if mode not in ("celltype", "profile", "plain"):
            raise ValueError(f"unknown mode {mode!r}")
        def opt_int(k):
            v = (q.get(k) or [None])[0]
            return int(v) if v not in (None, "", "null") else None

        sel, prof = opt_int("select"), opt_int("profile")
        mana = (q.get("mana") or ["0"])[0] in ("1", "true", "True")
        width = int((q.get("width") or [1400])[0])
        width = max(200, min(2400, width))
        with RENDER_LOCK:
            path = rd.background(b, mode=mode, select=sel, profile=prof, mana=mana,
                                 width=width, cache_dir=self.cache_dir(name))
        return ("file", path)


ROUTES = [
    ("/api/datasets", "datasets", 0),
    ("/api/dataset/", "dataset", 1),
    ("/api/centroids/", "centroids", 1),
    ("/api/clones/", "clones", 1),
    ("/api/set/", "cset", 1),
    ("/api/lineages/", "lineages", 1),
    ("/api/lineage/", "lineage", 2),
    ("/api/bg/", "bg", 1),
]


class Handler(BaseHTTPRequestHandler):
    server_version = "ClonoView"
    api = None
    token = None            # None means no auth (loopback only)

    def _authorised(self, q):
        """Token check, for when the server is not bound to loopback.

        These nodes have public IPv4 addresses, so binding to 0.0.0.0 puts the viewer
        wherever the security group allows -- possibly the open internet. The data is
        unpublished and patient-derived, so an open port is not acceptable even if the
        security group probably blocks it; "probably" is not a control. Same scheme
        Jupyter uses: a random token in the URL, swapped for a cookie on first load so
        the query string does not have to be repeated on every asset and API call.

        compare_digest, not ==, so the check does not leak the token's prefix through
        its own timing.
        """
        if self.token is None:
            return True
        got = (q.get("key") or [None])[0]
        if got is None:
            cookie = self.headers.get("Cookie") or ""
            for part in cookie.split(";"):
                k, _, v = part.strip().partition("=")
                if k == "cv_token":
                    got = v
                    break
        return got is not None and hmac.compare_digest(got, self.token)

    def log_message(self, fmt, *args):
        # Time and client address, both of which have already been needed once: with
        # neither, a request made by a curl check on this box is indistinguishable in the
        # log from one made by a browser through the tunnel, which is exactly the
        # question to ask when a click appears to do nothing.
        if self.server.verbose:
            sys.stderr.write("  %s %s %s\n" % (time.strftime("%H:%M:%S"),
                                               self.address_string(), fmt % args))

    def _send(self, code, body, ctype, extra=None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, obj, code=200):
        body = json.dumps(obj, default=_np_json).encode()
        # A whole patient's clone set runs to several hundred kilobytes of coordinates,
        # and this is usually reached down an SSH tunnel. JSON of rounded numbers
        # compresses about fourfold, and gzip costs a few milliseconds, so it is the
        # cheapest latency available. Only above a size where the header would dominate.
        if len(body) > 4096 and "gzip" in (self.headers.get("Accept-Encoding") or ""):
            body = gzip.compress(body, 5)
            return self._send(code, body, "application/json",
                              {"Content-Encoding": "gzip", "Vary": "Accept-Encoding"})
        self._send(code, body, "application/json")

    def do_HEAD(self):
        self.do_GET()

    def do_GET(self):
        u = urlparse(self.path)
        path, q = u.path, parse_qs(u.query)
        if not self._authorised(q):
            return self._send(403, b"forbidden: this ClonoView needs its access link. "
                                   b"Use the full URL that was printed when it started, "
                                   b"including the ?key=... part.", "text/plain")
        try:
            for prefix, fn, nargs in ROUTES:
                if nargs == 0 and path == prefix:
                    return self._dispatch(fn, q, [])
                if nargs and path.startswith(prefix):
                    rest = [p for p in path[len(prefix):].split("/") if p]
                    if len(rest) == nargs:
                        return self._dispatch(fn, q, rest)
            if path in ("/", "/index.html"):
                return self._static("index.html")
            if path.startswith("/static/"):
                return self._static(path[len("/static/"):])
            self._json(dict(error="not found", path=path), 404)
        except KeyError as e:
            self._json(dict(error=f"no such dataset: {e}"), 404)
        except (ValueError, IndexError) as e:
            self._json(dict(error=str(e)), 400)
        except Exception as e:                                  # noqa: BLE001
            traceback.print_exc()
            self._json(dict(error=f"{type(e).__name__}: {e}"), 500)

    def _dispatch(self, fn, q, args):
        out = getattr(self.api, fn)(q, *args)
        if isinstance(out, tuple) and out and out[0] == "file":
            with open(out[1], "rb") as fh:
                body = fh.read()
            # Rasters are content-addressed by their cache key, so a long cache is safe
            # and keeps the tunnel quiet when toggling back and forth.
            return self._send(200, body, "image/png",
                              {"Cache-Control": "public, max-age=86400"})
        self._json(out)

    def _static(self, rel):
        rel = posixpath.normpath("/" + rel).lstrip("/")
        p = os.path.join(STATIC, rel)
        if not os.path.abspath(p).startswith(STATIC) or not os.path.exists(p):
            return self._json(dict(error="not found", path=rel), 404)
        ctype = mimetypes.guess_type(p)[0] or "application/octet-stream"
        with open(p, "rb") as fh:
            body = fh.read()
        extra = {"Cache-Control": "no-cache"}
        if self.token and rel == "index.html":
            # Trade the token in the URL for a cookie on the first page load, so the
            # page's own asset and API requests carry it without every one of them
            # having to append the query string. SameSite=Strict: this is a viewer on a
            # private link, nothing should ever be able to make requests to it from
            # another site.
            extra["Set-Cookie"] = (f"cv_token={self.token}; Path=/; SameSite=Strict; "
                                   f"Max-Age=86400")
        self._send(200, body, ctype, extra)


def port_in_use(host, port):
    """Is something already listening, and does it look like a ClonoView?

    This exists because of a real hour lost: a ClonoView left running from an earlier
    session kept the port, so a browser reached THAT process. It served the current
    HTML and JavaScript -- static files are read from disk per request -- while running
    the old Python in memory, so a route added since simply 404ed and one feature
    silently did nothing. A stale server is worse than a dead one, and the only cheap
    defence is to say so loudly at startup.
    """
    with socket.socket() as s:
        # SO_REUSEADDR, because ThreadingHTTPServer sets allow_reuse_address and this
        # check must not be stricter than the thing it is checking for. Without it, a
        # port left in TIME_WAIT by a server that has just exited fails to bind here
        # while binding fine for the real server -- so restarting immediately after a
        # stop reported "in use by something that is not ClonoView" and refused to
        # start, on a port that was genuinely free.
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        s.settimeout(0.4)
        try:
            s.bind((host if host != "0.0.0.0" else "127.0.0.1", port))
            return None
        except OSError:
            pass
    probe = "127.0.0.1" if host in ("0.0.0.0", "127.0.0.1", "localhost") else host
    try:
        with urllib.request.urlopen(f"http://{probe}:{port}/api/datasets",
                                    timeout=1.5) as r:
            return "clonoview" if r.status in (200, 403) else "other"
    except urllib.error.HTTPError as e:
        return "clonoview" if e.code in (401, 403) else "other"
    except Exception:                                       # noqa: BLE001
        return "other"


def free_port(host, start):
    for p in range(start, start + 200):
        with socket.socket() as s:
            try:
                s.bind((host, p))
                return p
            except OSError:
                continue
    raise RuntimeError("no free port in range")


def public_ipv4():
    """This instance's public address, or None. Used only to print a usable URL."""
    try:
        req = urllib.request.Request(
            "http://169.254.169.254/latest/api/token", method="PUT",
            headers={"X-aws-ec2-metadata-token-ttl-seconds": "60"})
        tok = urllib.request.urlopen(req, timeout=2).read().decode()
        req = urllib.request.Request(
            "http://169.254.169.254/latest/meta-data/public-ipv4",
            headers={"X-aws-ec2-metadata-token": tok})
        return urllib.request.urlopen(req, timeout=2).read().decode().strip() or None
    except Exception:
        return None


def serve(data_root=DATA, host="127.0.0.1", port=0, verbose=True, open_msg=True,
          token=None):
    port = port or free_port(host, 8770)
    loopback = host in ("127.0.0.1", "localhost", "::1")
    # A token is mandatory off loopback and pointless on it. These nodes carry public
    # IPv4 addresses, so "bound to 0.0.0.0" can mean "on the internet" depending on a
    # security group this process cannot see; the data is unpublished and patient-
    # derived, so it does not get served without one.
    if token is None and not loopback:
        token = secrets.token_urlsafe(18)
    Handler.api = Api(data_root)
    Handler.token = token
    srv = ThreadingHTTPServer((host, port), Handler)
    srv.verbose = verbose
    srv.token = token
    ds = bd.available(data_root)
    if open_msg:
        node = socket.gethostname()
        q = f"?key={token}" if token else ""
        print(f"\nClonoView serving {len(ds)} dataset(s): "
              f"{', '.join(d['name'] for d in ds) or '(none prepared)'}")
        print(f"  bound to {host}:{port} on {node}")
        if loopback:
            print("\n  Loopback only. From your laptop, forward the port, then open"
                  f"\n  http://localhost:{port}/ :\n"
                  f"\n    ssh -L {port}:localhost:{port} <user>@<cluster-address>\n")
        else:
            ip = public_ipv4()
            print("\n  Open this link (it carries the access token -- treat it as a "
                  "password):\n")
            if ip:
                print(f"    http://{ip}:{port}/{q}\n")
            print(f"    http://{node}:{port}/{q}       (from inside the cluster network)\n")
            print("  If the public link does not load, the port is closed at the "
                  "firewall\n  and an SSH tunnel is the way in.\n")
        print("Ctrl-C here to stop.\n", flush=True)
    return srv


def main(argv=None):
    ap = argparse.ArgumentParser(description="Serve the ClonoView viewer.")
    ap.add_argument("dataset", nargs="?",
                    help="dataset to open on load; omit to choose in the app")
    ap.add_argument("--data", default=DATA, help="directory of prepared bundles")
    ap.add_argument("--port", type=int, default=0, help="0 picks a free port from 8770")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--token", default=None,
                    help="access token; generated automatically when the "
                         "bind address is not loopback")
    a = ap.parse_args(argv)
    if a.dataset and not os.path.isdir(os.path.join(a.data, a.dataset)):
        sys.exit(f"dataset {a.dataset!r} is not prepared in {a.data}. Run:\n"
                 f"    python -m clonoview.prepare {a.dataset}")

    if a.port:
        busy = port_in_use(a.host, a.port)
        if busy == "clonoview":
            sys.exit(
                f"\nA ClonoView is ALREADY running on port {a.port} of "
                f"{socket.gethostname()}.\n\n"
                f"  It is probably left over from an earlier session. It will keep\n"
                f"  serving the current web pages while running whatever code it was\n"
                f"  started with, so newer features can appear to be broken.\n\n"
                f"  Stop it and start again:\n"
                f"      ./clonoview.sh --stop\n"
                f"      ./clonoview.sh --port {a.port}\n")
        if busy:
            sys.exit(f"\nPort {a.port} on {socket.gethostname()} is in use by something "
                     f"that is not ClonoView.\nPick another: ./clonoview.sh --port "
                     f"{a.port + 1}\n")
    srv = serve(a.data, a.host, a.port, verbose=not a.quiet, token=a.token)

    # The SERVER writes the pidfile, not the launcher, and only now -- after the socket
    # is bound. The launcher used to write it from $! before knowing whether the start
    # succeeded, so a run that aborted (say, because the port was already taken)
    # overwrote the pid of the healthy server it had just refused to replace, and
    # --stop then reported that pid was "not a clonoview server" and left the real one
    # running. Writing it here means the file only ever names a process that got as far
    # as listening.
    # PER HOST, and in the user's cache rather than beside the code. The hostname keeps
    # two nodes sharing one networked home directory apart: without it, starting the
    # viewer on a login node overwrites the compute node's entry, and --stop then signals
    # that pid on the WRONG machine, where it belongs to somebody else's process.
    os.makedirs(RUN_DIR, exist_ok=True)
    pidfile = os.path.join(RUN_DIR, f"clonoview.{socket.gethostname().split('.')[0]}.pid")
    try:
        with open(pidfile, "w") as fh:
            fh.write(str(os.getpid()))
    except OSError:
        pidfile = None
    if a.dataset:
        print(f"opening with dataset: {a.dataset}  "
              f"(http://localhost:{srv.server_port}/?dataset={a.dataset})\n", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    finally:
        if pidfile and os.path.exists(pidfile):
            try:
                with open(pidfile) as fh:
                    if fh.read().strip() == str(os.getpid()):
                        os.remove(pidfile)     # only if it is still ours
            except OSError:
                pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
