"""Load a prepared bundle and answer the questions the viewer asks of it.

Everything here is a slice or a bincount over the arrays `prepare.py` wrote. The one
structural thing worth knowing is that cells are stored sorted by (lineage, visit), so
`cells_of_lineage` is a slice rather than a scan -- which is what makes scrubbing a
clone through time fast enough to animate.
"""

import json
import os
import threading

import numpy as np

from . import colour as C


class Bundle:
    def __init__(self, path):
        self.path = path
        self.name = os.path.basename(path)
        with open(os.path.join(path, "meta.json")) as fh:
            self.meta = json.load(fh)
        with open(os.path.join(path, "names.json")) as fh:
            self.names = json.load(fh)
        c = np.load(os.path.join(path, "cells.npz"))
        u = np.load(os.path.join(path, "clones.npz"))
        self.cxy = c["xy"]
        self.c_lineage, self.c_unit = c["lineage"], c["clone_unit"]
        self.c_visit, self.c_patient = c["visit"], c["patient"]
        self.c_type, self.c_profile = c["cell_type"], c["profile"]
        # Anatomic position, optional. -1 where a dataset has none, so every code path
        # below works unchanged on the four cohorts that were sampled at one site.
        self.c_tissue = c["tissue"] if "tissue" in c.files \
            else np.full(len(c["xy"]), -1, dtype=np.int8)
        self.c_prob, self.c_row = c["prof_prob"], c["source_row"]
        # MANAscore is optional per dataset. Three states per cell: 1 high, 0 scored and
        # not high, -1 never scored (it runs on the CD8 compartment only, so "never
        # scored" is the majority on both datasets that have it).
        self.c_mana = c["manascore_hi"] if "manascore_hi" in c.files \
            else np.full(len(self.cxy), -1, dtype=np.int8)
        self.c_mscore = c["manascore"] if "manascore" in c.files \
            else np.full(len(self.cxy), np.nan, dtype=np.float32)
        self.lin_off = c["lineage_offset"]
        self.uxy = u["xy"]
        # Alternative clone layouts, in the same clone order as `uxy` (prepare.py
        # reindexes them so a switch swaps coordinates under a fixed set of points).
        # `default` always exists and is `uxy`; anything else came from the manifest.
        self._layouts = {"default": self.uxy}
        for k in u.files:
            if k.startswith("xy_alt"):
                self._layouts[k[3:]] = u[k]
        self.u_profile, self.u_lineage = u["profile"], u["lineage"]
        self.u_visit, self.u_patient, self.u_n = u["visit"], u["patient"], u["n_cells"]
        self.u_ptrc = u["ptrc"] if "ptrc" in u.files \
            else np.full(len(self.uxy), -1, dtype=np.int8)
        self.u_mana_frac = u["mana_frac"] if "mana_frac" in u.files \
            else np.full(len(self.uxy), np.nan, dtype=np.float32)
        self.u_mana_usable = u["mana_usable"] if "mana_usable" in u.files \
            else np.ones(len(self.uxy), dtype=np.int8)
        self.visits = [v["code"] for v in self.meta["visits"]]
        self.profile_ids = [p["id"] for p in self.meta["profiles"]]
        self.cell_types = [t["name"] for t in self.meta["cell_types"]]
        self.tissues = [t["code"] for t in self.meta.get("tissues", [])]
        self._lock = threading.Lock()

    def layouts(self):
        """[{id, label, extent}] for the clone panel, default first."""
        return self.meta.get("clone_layouts") or [
            dict(id="default", label="UMAP", extent=self.meta["extent"]["clone"])]

    def uxy_for(self, layout):
        """Clone-unit coordinates under one layout. An unknown id falls back to the
        default rather than raising: a stale bookmark should show the map, not a 500."""
        return self._layouts.get(layout or "default", self.uxy)

    def clone_extent(self, layout):
        for L in self.layouts():
            if L["id"] == (layout or "default"):
                return L["extent"]
        return self.meta["extent"]["clone"]

    @property
    def has_manascore(self):
        return self.meta.get("manascore") is not None

    # ---- colours -------------------------------------------------------------
    @property
    def profile_colours(self):
        """Profile k is drawn in slot k-1 of the validated palette, the same colour it
        carries in the paper figures."""
        return {int(p): C.PROFILE_COLOURS[(int(p) - 1) % len(C.PROFILE_COLOURS)]
                for p in self.profile_ids}

    @property
    def cell_type_colours(self):
        return C.cell_type_palette(len(self.cell_types))

    # ---- slices --------------------------------------------------------------
    def cells_of_lineage(self, li):
        """(start, stop) into the cell arrays -- one contiguous slice."""
        return int(self.lin_off[li]), int(self.lin_off[li + 1])

    def lineage_frames(self, li, layout=None, axis="time", fixed=-1):
        """One frame per step along `axis`, for one lineage.

        TWO AXES, AND THEY ARE NOT SYMMETRIC. Time and anatomic position are both
        per-cell covariates, so the CELL panel varies along either one. The CLONE map
        does not: Clonotrace placed clone units at (patient, visit, receptor), so a
        lineage has a position per VISIT and never one per tissue. That asymmetry is real
        and is surfaced rather than papered over --

          axis="time"    frames are visits; `fixed` optionally restricts cells to one
                         anatomic position (-1 = all). Each frame carries its own clone
                         unit and position, so arrows between frames mean something.
          axis="tissue"  frames are anatomic positions at ONE visit, named by `fixed`.
                         Every frame carries the SAME clone unit -- the one at that visit
                         -- because there is no per-tissue position to carry. The left
                         panel therefore holds still while the right varies, which is the
                         honest rendering of "where did this clone's cells sit, by site,
                         at this timepoint".

        `axis="tissue"` requires a real visit in `fixed`. That is not a limitation of the
        implementation but of the storage it protects: cells are sorted by
        (lineage, visit, tissue), so one (visit, tissue) pair is a contiguous slice while
        "one tissue across all visits" is not, and a frame is served as a slice.
        """
        a, b = self.cells_of_lineage(li)
        uxy = self.uxy_for(layout)
        vis, tis = self.c_visit[a:b], self.c_tissue[a:b]
        units = {int(self.u_visit[k]): int(k)
                 for k in np.flatnonzero(self.u_lineage == li)}

        def frame(vi, ti, k):
            m = np.flatnonzero((vis == vi) & (tis == ti)) if ti >= 0 \
                else np.flatnonzero(vis == vi)
            f = dict(
                visit=vi, code=self.visits[vi],
                label=self.meta["visits"][vi]["label"],
                n_cells=int(len(m)),
                cell_slice=[a + int(m[0]), a + int(m[-1]) + 1] if len(m) else None,
                clone_unit=k,
                xy=[float(uxy[k, 0]), float(uxy[k, 1])] if k is not None else None,
                profile=int(self.u_profile[k]) if k is not None else None,
                below_threshold=k is None,
            )
            if ti >= 0:
                f["tissue"] = ti
                f["tissue_label"] = self.tissues[ti]
            return f

        if axis == "tissue":
            vi = int(fixed)
            if vi < 0 or vi >= len(self.visits):
                raise ValueError("axis='tissue' needs a visit index in `fixed`")
            k = units.get(vi)
            out = []
            for ti in range(len(self.tissues)):
                f = frame(vi, ti, k)
                if f["n_cells"] == 0:
                    continue          # never sampled here; do not draw an empty panel
                f["label"] = self.tissues[ti]
                f["code"] = self.tissues[ti]
                out.append(f)
            return out

        ti = int(fixed)
        ti = ti if 0 <= ti < len(self.tissues) else -1
        out = []
        for vi in range(len(self.visits)):
            k = units.get(vi)
            f = frame(vi, ti, k)
            if f["n_cells"] == 0 and k is None:
                continue
            out.append(f)
        return out

    def axes(self):
        """What the viewer may vary. Time always; anatomic position only where the
        cohort was sampled at more than one site."""
        out = [dict(id="time", label="Time",
                    steps=[dict(code=v["code"], label=v["label"])
                           for v in self.meta["visits"]])]
        if len(self.tissues) > 1:
            out.append(dict(id="tissue", label="Anatomic position",
                            steps=[dict(code=t["code"], label=t["label"],
                                        n_cells=t["n_cells"])
                                   for t in self.meta["tissues"]]))
        return out

    def centroids(self, values, labels):
        """Median position per class, for the direct labels. Median not mean: a diffuse
        profile's mean can land in a hole where none of its cells are."""
        out = {}
        for lab in labels:
            m = values == lab
            if int(m.sum()) == 0:
                continue
            out[int(lab)] = [float(np.median(self.cxy[m, 0])),
                             float(np.median(self.cxy[m, 1]))]
        return out


_CACHE = {}
_CACHE_LOCK = threading.Lock()


def load(root, name):
    with _CACHE_LOCK:
        if name not in _CACHE:
            p = os.path.join(root, name)
            if not os.path.isdir(p):
                raise KeyError(name)
            _CACHE[name] = Bundle(p)
        return _CACHE[name]


def available(root):
    out = []
    for n in sorted(os.listdir(root)):
        mp = os.path.join(root, n, "meta.json")
        if os.path.exists(mp):
            with open(mp) as fh:
                m = json.load(fh)
            out.append(dict(name=m["name"], label=m["label"], kind=m["kind"],
                            built=m.get("built"), counts=m["counts"],
                            notes=m.get("notes", "")))
    return out
