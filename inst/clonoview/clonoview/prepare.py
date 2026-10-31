"""Compile a Clonotrace run into a compact bundle the viewer can load instantly.

Why a bundle at all. The raw inputs are large and spread across trees -- Wang's
annotation file alone is 175 MB, its cell embedding 26 MB -- and the login shell is
capped at about 2 GB. Parsing them on every server start would be slow and would put
the viewer's memory use at the mercy of whatever else is in those files. So `prepare`
reads them once, streaming, keeps only the columns the viewer draws, and writes fixed
-width numpy arrays plus one small JSON file. The server then memory-maps a few
megabytes and starts in well under a second.

The layout, in `<out>/<name>/`:

  meta.json     everything small and stringy: visits, profiles, cell types, the patient
                table with its per-visit clone and cell counts, axis extents, and a
                `checks` block re-deriving every count from the arrays that were written
  names.json    patient / cell-type / lineage / clone-unit name tables
  cells.npz     one row per cell in the cell embedding
  clones.npz    one row per clone unit in the run

CELLS ARE SORTED, AND THAT IS LOad-BEARING. Cells belonging to a tracked lineage come
first, ordered by (lineage, visit); everything else follows. `lineage_offset` then gives
each lineage a contiguous slice, so drawing one clone's cells across all its timepoints
is a single array slice rather than a scan of 273,508 rows.

WHAT A LINEAGE IS, AND WHY IT IS NOT A CLONE UNIT. Clonotrace only models clones that
cleared the run's `expanded_min` cell threshold at a given visit -- those are the *clone
units*, and they are the points in the clone embedding. The same T-cell lineage,
(patient, receptor), typically also has cells at visits where it fell below that
threshold: 57,257 such cells in Mathew, 82,894 in Wang. Those cells are real and are
kept here, attached to the lineage with `clone_unit = -1`. The viewer draws them and
marks the frame as having no clone-embedding position, which is the honest rendering --
Clonotrace never placed that lineage on the map at that visit, so neither do we.
"""

import csv
import itertools
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict

import numpy as np

from . import manifest as mf

csv.field_size_limit(1 << 24)   # receptor strings can be long

# How R spells a missing number on the way out to a TSV.
NUM = {"NA": "nan", "NaN": "nan", "": "nan", "NULL": "nan"}


# ---- streaming TSV helpers ---------------------------------------------------
def stream(path, cols, header_offset=False, id_col=None, delim="\t"):
    """Yield dicts of the requested columns, one row at a time.

    `header_offset` handles a file whose header has no name for its row-id column, so
    every named field sits one place to the right of where the header says it does --
    Mathew's published annotation is written that way. `id_col=None` with an offset
    header means the row id is field 0 and is returned under the key "__id__".
    """
    with open(path, newline="") as fh:
        rdr = csv.reader(fh, delimiter=delim)
        hdr = next(rdr)
        shift = 1 if header_offset else 0
        idx = {}
        for c in cols:
            if c is None:
                continue
            if c not in hdr:
                raise KeyError(f"{path}: no column {c!r}; header has {hdr[:12]}"
                               f"{' ...' if len(hdr) > 12 else ''}")
            idx[c] = hdr.index(c) + shift
        want_id = header_offset and id_col is None
        for row in rdr:
            if not row:
                continue
            out = {c: row[i] for c, i in idx.items()}
            if want_id:
                out["__id__"] = row[0]
            yield out


def read_xy(path):
    """`<id> <x> <y>` -> ids, x, y.

    Deliberately tolerant, because the four datasets spell this file three different
    ways: `umap_1`/`umap_2` lowercase with a named id column (Wang, Mathew, Liu),
    `UMAP_1`/`UMAP_2` capitalised (GBM), and a header naming only the two coordinate
    columns so the row id sits in an unnamed leading field (GBM again). The offset is
    DETECTED, by comparing the header against the first row, rather than declared -- a
    file that changes shape is then caught here instead of silently loading coordinates
    one column to the left.
    """
    ids, xs, ys = [], [], []
    with open(path, newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        first = next(rdr)
        shift = len(first) - len(hdr)
        if shift not in (0, 1):
            raise ValueError(f"{path}: header has {len(hdr)} fields but rows have "
                             f"{len(first)}; cannot tell where the id column is")
        low = [h.strip().lower() for h in hdr]

        def find(*names):
            for n in names:
                if n in low:
                    return low.index(n) + shift
            raise KeyError(f"{path}: no coordinate column among {names}; header is {hdr}")

        ix, iy = find("umap_1", "umap1", "x"), find("umap_2", "umap2", "y")
        for row in itertools.chain([first], rdr):
            if row:
                ids.append(row[0])
                xs.append(row[ix])
                ys.append(row[iy])
    return ids, np.array(xs, dtype=np.float64), np.array(ys, dtype=np.float64)


def extent(x, y, pad=0.03):
    """Fixed drawing extent, with a small margin.

    The server renders the cell background once into a raster at exactly this extent and
    the browser draws contours over it, so the two only line up if the extent is pinned
    here and never recomputed per view.
    """
    x0, x1 = float(np.min(x)), float(np.max(x))
    y0, y1 = float(np.min(y)), float(np.max(y))
    mx, my = (x1 - x0) * pad, (y1 - y0) * pad
    return [x0 - mx, x1 + mx, y0 - my, y1 + my]


# ---- the per-cell profile ----------------------------------------------------
def read_profiles(run_dir, cell_index):
    """argmax profile and its probability, per cell.

    The columns are named X1..XK but are NOT written in numeric order -- Mathew's file
    heads X1, X3, X2, X4 -- so they are sorted by the integer in the name. Reading them
    positionally would relabel profiles 2 and 3 across the entire dataset.

    A cell absent from this file was dropped by Clonotrace's deviance filter or carries
    no clone; it gets profile 0, meaning none, and is drawn in the context grey.
    """
    path = os.path.join(run_dir, "cell_profile_prob.tsv")
    prof = np.zeros(len(cell_index), dtype=np.int8)
    prob = np.zeros(len(cell_index), dtype=np.float32)
    n_seen = n_unknown = n_allnan = 0
    with open(path, newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        cols = [(int(re.fullmatch(r"X(\d+)", c).group(1)), i)
                for i, c in enumerate(hdr) if re.fullmatch(r"X(\d+)", c)]
        cols.sort()
        order = [i for _, i in cols]
        labels = np.array([k for k, _ in cols], dtype=np.int8)
        for row in rdr:
            if not row:
                continue
            n_seen += 1
            j = cell_index.get(row[0])
            if j is None:
                n_unknown += 1
                continue
            # R writes a missing probability as the literal string "NA", so this cannot
            # go straight to float. Wang has five rows that are NA the whole way across
            # -- the ones S_shift_wang_melanoma/check_nan_propagates.R predicted. They
            # are confined to this file; the clone distance matrix has no NaN at all.
            v = np.array([NUM.get(row[i], row[i]) for i in order], dtype=np.float64)
            if not np.isfinite(v).any():
                # Known and audited: five Wang cells have an all-NaN probability row.
                # They stay profile 0 and are drawn grey rather than given a spurious
                # argmax of column one.
                n_allnan += 1
                continue
            v = np.nan_to_num(v, nan=-1.0)
            k = int(np.argmax(v))
            prof[j] = labels[k]
            prob[j] = max(0.0, float(v[k]))
    return prof, prob, dict(rows=n_seen, not_in_embedding=n_unknown, all_nan=n_allnan,
                            n_profiles=len(cols))


TRUEISH = {"TRUE", "True", "true", "T", "1", "1.0", "YES", "Y"}


def read_manascore(ms, cell_index, clone_ids):
    """Per-cell MANAscore calls and per-clone pTRC calls.

    THREE STATES, NOT TWO. A cell is high (1), scored-and-not-high (0), or **never
    scored** (-1). The distinction is not cosmetic: MANAscore runs on the CD8 compartment,
    so it covers 58,810 of Liu's 116,598 cells and 62,883 of GBM's 164,857. Collapsing
    "not scored" into "not high" would paint two thirds of the map as evidence of absence
    when it is absence of evidence. The viewer draws only the high cells and says how many
    were scored.

    The clone table is already keyed by the Clonotrace clone identifier, so nothing has to
    be matched here. `usable` records whether that clone's patient had a usable cutoff at
    all -- three patients across the two datasets did not, and a pTRC of 0 for those means
    "could not be called", not "called negative".
    """
    c = ms["cells"]
    hi = np.full(len(cell_index), -1, dtype=np.int8)
    score = np.full(len(cell_index), np.nan, dtype=np.float32)
    n_rows = n_unknown = 0
    cols = [c["cell"], c["hi"]] + ([c["score"]] if c.get("score") else [])
    for r in stream(c["path"], cols):
        n_rows += 1
        j = cell_index.get(r[c["cell"]])
        if j is None:
            n_unknown += 1
            continue
        hi[j] = 1 if r[c["hi"]] in TRUEISH else 0
        if c.get("score"):
            try:
                score[j] = float(NUM.get(r[c["score"]], r[c["score"]]))
            except ValueError:
                pass

    k = ms["clones"]
    idx = {cid: i for i, cid in enumerate(clone_ids)}
    ptrc = np.full(len(clone_ids), -1, dtype=np.int8)
    frac = np.full(len(clone_ids), np.nan, dtype=np.float32)
    usable = np.ones(len(clone_ids), dtype=np.int8)
    cols = [k["clone"], k["ptrc"]] + [k[x] for x in ("frac", "usable") if k.get(x)]
    n_clone_rows = n_clone_unmatched = 0
    for r in stream(k["path"], cols):
        n_clone_rows += 1
        i = idx.get(r[k["clone"]])
        if i is None:
            n_clone_unmatched += 1
            continue
        ptrc[i] = 1 if r[k["ptrc"]] in TRUEISH else 0
        if k.get("frac"):
            try:
                frac[i] = float(NUM.get(r[k["frac"]], r[k["frac"]]))
            except ValueError:
                pass
        if k.get("usable"):
            usable[i] = 1 if r[k["usable"]] in TRUEISH else 0

    stats = dict(
        cells_scored=int((hi >= 0).sum()), cells_high=int((hi == 1).sum()),
        cell_rows=n_rows, cell_rows_not_in_embedding=n_unknown,
        clone_rows=n_clone_rows, clone_rows_unmatched=n_clone_unmatched,
        clone_units_called=int((ptrc >= 0).sum()),
        clone_units_ptrc=int((ptrc == 1).sum()),
        clone_units_no_cutoff=int((usable == 0).sum()),
        note=ms.get("note", ""),
        sources=dict(cells=c["path"], clones=k["path"]),
    )
    return hi, score, ptrc, frac, usable, stats


# ---- main --------------------------------------------------------------------
def build(m, out_root, verbose=True):
    t0 = time.time()

    def say(*a):
        if verbose:
            print(*a, flush=True)

    rx = mf.clone_id_regex(m)
    one_patient = mf.patient_default(m)
    visits = list(m["visits"])
    visit_rank = {v: i for i, v in enumerate(visits)}
    aliases = dict(m.get("visit_aliases") or {})

    def canon_visit(raw):
        return aliases.get(raw, raw)

    # -- 1. clone units: profile, position, and the (patient, visit, receptor) triple --
    clusters = {r["clone"]: int(r["cluster"])
                for r in stream(os.path.join(m["run_dir"], "clone_clusters.tsv"),
                                ["clone", "cluster"])}
    emb_path = mf.clone_embedding_path(m)
    cids, cx, cy = read_xy(emb_path)
    say(f"  clone embedding {os.path.basename(emb_path)}: {len(cids)} clone units")
    if set(cids) != set(clusters):
        raise ValueError(f"clone embedding and clone_clusters disagree on membership: "
                         f"{len(set(cids) - set(clusters))} only in embedding, "
                         f"{len(set(clusters) - set(cids))} only in clusters")

    # Alternative clone layouts (e.g. a diffusion map beside the UMAP). Each is read in
    # the DEFAULT layout's clone order, not its own file order: the viewer switches
    # layouts by swapping a coordinate array under a fixed set of clone units, so a file
    # that happened to be written in a different order would silently move every point
    # onto the wrong clone. The reindex below is what makes the toggle safe, and the
    # membership check makes a mismatched file an error rather than a scramble.
    extra_layouts = {}
    for lbl, lpath in (m.get("clone_layouts") or {}).items():
        lids, lx, ly = read_xy(lpath)
        if set(lids) != set(cids):
            raise ValueError(
                f"clone layout {lbl!r} covers different clones than the default layout: "
                f"{len(set(lids) - set(cids))} only in {os.path.basename(lpath)}, "
                f"{len(set(cids) - set(lids))} only in the default")
        pos = {c: i for i, c in enumerate(lids)}
        take = [pos[c] for c in cids]
        extra_layouts[lbl] = (np.asarray(lx)[take], np.asarray(ly)[take])
        say(f"  clone layout {lbl!r} from {os.path.basename(lpath)}: "
            f"{len(lids)} clone units, reindexed to the default order")

    cu_patient, cu_visit, cu_receptor = [], [], []
    for c in cids:
        mt = rx.match(c)
        if not mt:
            raise ValueError(f"clone_id_regex did not match clone id {c!r}")
        # All three groups are optional, and each absence means something specific.
        # No `visit`: the dataset has no time dimension -- one tumour resection per
        # patient, or one barcoding library -- so every clone unit belongs to the single
        # declared visit. No `patient`: there is one sample, named by the manifest. No
        # `receptor`: the clone identifier is opaque and is the lineage key whole, which
        # is the lineage-barcoding case. Requiring the groups would have meant inventing
        # a fake timepoint or a fake subject inside the regex.
        gd = mt.groupdict()
        v = gd.get("visit") or visits[0]
        if v not in visit_rank:
            raise ValueError(f"clone {c!r} has visit {v!r}, which is not in the "
                             f"manifest's `visits` list {visits}")
        cu_patient.append(gd.get("patient") or one_patient)
        cu_visit.append(v)
        cu_receptor.append(gd.get("receptor") or c)

    # A lineage is one T-cell lineage in one patient, tracked across visits. This is the
    # object the viewer scrolls through time; a clone unit is that lineage at one visit.
    lin_key = list(dict.fromkeys(zip(cu_patient, cu_receptor)))
    lin_id = {k: i for i, k in enumerate(lin_key)}
    say(f"  {len(cids)} clone units over {len(lin_id)} lineages")

    # -- 2. cell embedding: defines the universe of drawable cells --
    cell_ids, ex, ey = read_xy(m["cell_umap"])
    cell_index = {c: i for i, c in enumerate(cell_ids)}
    n_cells = len(cell_ids)
    say(f"  cell embedding: {n_cells} cells")

    # -- 3. cell -> clone, then cell -> lineage / clone unit --
    cc = m["cell_clone"]
    cu_index = {c: i for i, c in enumerate(cids)}
    c_lineage = np.full(n_cells, -1, dtype=np.int32)
    c_unit = np.full(n_cells, -1, dtype=np.int32)
    c_visit_from_clone = np.full(n_cells, -1, dtype=np.int8)
    n_map = n_map_unknown = 0
    for r in stream(cc["path"], [cc["cell"], cc["clone"]]):
        j = cell_index.get(r[cc["cell"]])
        n_map += 1
        if j is None:
            n_map_unknown += 1
            continue
        clone = r[cc["clone"]]
        mt = rx.match(clone)
        if not mt:
            continue
        # Same optional groups as the clone-unit pass above, and resolved the same way,
        # so a lineage key built there is the key looked up here.
        gd = mt.groupdict()
        li = lin_id.get((gd.get("patient") or one_patient, gd.get("receptor") or clone))
        if li is None:
            continue            # lineage never cleared the threshold at any visit
        c_lineage[j] = li
        c_unit[j] = cu_index.get(clone, -1)
        c_visit_from_clone[j] = visit_rank.get(gd.get("visit") or visits[0], -1)
    say(f"  cell->clone: {n_map} rows, {int((c_lineage >= 0).sum())} cells on a tracked "
        f"lineage, of which {int((c_unit >= 0).sum())} are in a clone unit")

    # -- 4. annotation: cell type, patient, visit, patient-level fields --
    ann = m.get("cell_annotation")
    cell_type_name, c_type = ["(unannotated)"], np.full(n_cells, -1, dtype=np.int16)
    patient_name, c_patient = [], np.full(n_cells, -1, dtype=np.int16)
    c_visit = np.full(n_cells, -1, dtype=np.int8)
    # Anatomic position, a SECOND per-cell axis beside time. Optional: only a cohort
    # sampled at more than one site has one.
    tissue_name, c_tissue = [], np.full(n_cells, -1, dtype=np.int8)
    patient_fields = defaultdict(dict)
    field_conflicts = defaultdict(set)
    unmapped_visits = Counter()
    if ann:
        type_id, pat_id, tis_id = {}, {}, {}
        pf = dict(ann.get("patient_fields") or {})
        cols = [ann.get("cell"), ann.get("cell_type"), ann.get("patient"),
                ann.get("visit"), ann.get("tissue")] + list(pf.values())
        cols = [c for c in cols if c is not None or ann.get("cell") is None]
        for r in stream(ann["path"], cols, header_offset=ann.get("header_offset", False),
                        id_col=ann.get("cell"), delim=ann.get("delim", "\t")):
            cid = r["__id__"] if ann.get("cell") is None else r[ann["cell"]]
            j = cell_index.get(cid)
            if j is None:
                continue
            if ann.get("cell_type"):
                t = r[ann["cell_type"]]
                if t and t != "NA":
                    if t not in type_id:
                        type_id[t] = len(type_id)
                    c_type[j] = type_id[t]
            p = r[ann["patient"]] if ann.get("patient") else None
            if p and p != "NA":
                if p not in pat_id:
                    pat_id[p] = len(pat_id)
                c_patient[j] = pat_id[p]
                for disp, col in pf.items():
                    val = r[col]
                    if val and val != "NA":
                        # A patient-level field must be constant over that patient's
                        # cells -- it is a property of the person, and the viewer shows
                        # one value per patient. Where it is not, the last cell read used
                        # to win silently, so a response label could come out reversed
                        # with nothing to show for it. Collected here and raised below,
                        # so one run reports every field that disagrees rather than
                        # stopping at the first.
                        was = patient_fields[p].setdefault(disp, val)
                        if was != val:
                            field_conflicts[(p, disp)].update([was, val])
            if ann.get("tissue"):
                tv = r[ann["tissue"]]
                if tv and tv != "NA":
                    if tv not in tis_id:
                        tis_id[tv] = len(tis_id)
                    c_tissue[j] = tis_id[tv]
            if ann.get("visit"):
                v = canon_visit(r[ann["visit"]])
                if v in visit_rank:
                    c_visit[j] = visit_rank[v]
                elif r[ann["visit"]] and r[ann["visit"]] != "NA":
                    unmapped_visits[r[ann["visit"]]] += 1
        cell_type_name = [t for t, _ in sorted(type_id.items(), key=lambda kv: kv[1])]
        patient_name = [p for p, _ in sorted(pat_id.items(), key=lambda kv: kv[1])]
        # Ordered by cell count, most-sampled site first, so the axis control reads
        # sensibly rather than in file order.
        if tis_id:
            cnt = np.bincount(c_tissue[c_tissue >= 0], minlength=len(tis_id))
            raw = [t for t, _ in sorted(tis_id.items(), key=lambda kv: kv[1])]
            rank = np.argsort(-cnt, kind="stable")
            remap = np.full(len(raw), -1, dtype=np.int8)
            for new, old in enumerate(rank):
                remap[old] = new
            c_tissue = np.where(c_tissue >= 0, remap[c_tissue], -1).astype(np.int8)
            tissue_name = [raw[i] for i in rank]
            say(f"  anatomic positions: {len(tissue_name)} "
                f"({', '.join(f'{tissue_name[i]} {int(cnt[rank[i]]):,}' for i in range(len(tissue_name)))})")
        say(f"  annotation: {len(cell_type_name)} cell types, {len(patient_name)} patients")

    if field_conflicts:
        shown = ", ".join(f"{p}.{d} is both {' and '.join(sorted(v))}"
                          for (p, d), v in sorted(field_conflicts.items())[:6])
        raise ValueError(
            f"{len(field_conflicts)} patient-level field(s) take more than one value "
            f"within one patient: {shown}. A `patient_fields` column must be constant "
            f"over each patient's cells; if it varies per cell it is not a patient-level "
            f"field.")

    if unmapped_visits:
        # Not a warning. An unmapped visit means the annotation and the clone identifier
        # disagree on how the same visit is spelled, which silently empties the patient
        # table; `visit_aliases` exists to fix it.
        raise ValueError(
            f"{sum(unmapped_visits.values())} cells carry an annotation visit that is "
            f"neither in `visits` nor in `visit_aliases`: "
            f"{dict(unmapped_visits.most_common(8))}. Add them to `visit_aliases`.")

    # Fall back to the clone identifier where the annotation has no opinion. This is what
    # lets a dataset with no annotation file still populate a patient table.
    fill = (c_visit < 0) & (c_visit_from_clone >= 0)
    c_visit[fill] = c_visit_from_clone[fill]
    if not patient_name:
        pat_id = {}
        for i, p in enumerate(cu_patient):
            pat_id.setdefault(p, len(pat_id))
        patient_name = [p for p, _ in sorted(pat_id.items(), key=lambda kv: kv[1])]
        for j in np.flatnonzero(c_lineage >= 0):
            c_patient[j] = pat_id[lin_key[c_lineage[j]][0]]
    pat_index = {p: i for i, p in enumerate(patient_name)}

    for p in {lp for lp, _ in lin_key}:
        if p not in pat_index:
            raise ValueError(f"patient {p!r} appears in a clone identifier but not in "
                             f"the annotation's patient column -- the two sources name "
                             f"patients differently and the patient table would split")

    # -- 5. per-cell profile --
    c_prof, c_prob, prof_stats = read_profiles(m["run_dir"], cell_index)
    say(f"  profiles: {prof_stats['n_profiles']} classes, "
        f"{int((c_prof > 0).sum())} cells assigned, "
        f"{prof_stats['all_nan']} all-NaN rows left unassigned")

    # -- 5b. MANAscore, if this dataset has it ---------------------------------
    ms = m.get("manascore")
    if ms:
        c_mana, c_mscore, u_ptrc, u_frac, u_usable, mana_stats = read_manascore(
            ms, cell_index, cids)
        say(f"  manascore: {mana_stats['cells_high']} of "
            f"{mana_stats['cells_scored']} scored cells are high; "
            f"{mana_stats['clone_units_ptrc']} of {mana_stats['clone_units_called']} "
            f"clone units are pTRC"
            + (f" ({mana_stats['clone_units_no_cutoff']} had no usable patient cutoff)"
               if mana_stats["clone_units_no_cutoff"] else ""))
        if mana_stats["clone_rows_unmatched"]:
            say(f"    note: {mana_stats['clone_rows_unmatched']} rows of the MANAscore "
                f"clone table match no clone unit in this run")
    else:
        c_mana = np.full(n_cells, -1, dtype=np.int8)
        c_mscore = np.full(n_cells, np.nan, dtype=np.float32)
        u_ptrc = np.full(len(cids), -1, dtype=np.int8)
        u_frac = np.full(len(cids), np.nan, dtype=np.float32)
        u_usable = np.ones(len(cids), dtype=np.int8)
        mana_stats = None

    # -- 6. lineage display names: `<patient> #<k>`, k by size within the patient --
    lin_cells = np.zeros(len(lin_key), dtype=np.int64)
    if len(lin_key):
        np.add.at(lin_cells, c_lineage[c_lineage >= 0], 1)
    by_patient = defaultdict(list)
    for i, (p, rec) in enumerate(lin_key):
        by_patient[p].append(i)
    lin_label = [None] * len(lin_key)
    for p, idxs in by_patient.items():
        # Biggest lineage is #1. Ties break on the receptor string so the numbering is
        # identical on every rebuild.
        idxs.sort(key=lambda i: (-lin_cells[i], lin_key[i][1]))
        for k, i in enumerate(idxs, start=1):
            lin_label[i] = f"{p} #{k}"

    # -- 7. sort cells so each lineage is one contiguous slice --
    # TISSUE IS THE THIRD SORT KEY, and it is load-bearing exactly as the visit key is.
    # Frames are served as a contiguous [first, last+1] slice of the lineage's block, so
    # a frame can only be a slice if the cells of one (visit, tissue) pair sit together.
    # Sorting by (lineage, visit, tissue) keeps the existing (lineage, visit) contiguity
    # AND makes the finer pair contiguous, so both axes slice.
    tracked = c_lineage >= 0
    nt = max(len(tissue_name), 1) + 1
    key = ((np.where(tracked, c_lineage.astype(np.int64), len(lin_key)) * (len(visits) + 1)
            + np.where(c_visit >= 0, c_visit, len(visits)).astype(np.int64)) * nt
           + np.where(c_tissue >= 0, c_tissue, nt - 1).astype(np.int64))
    order = np.argsort(key, kind="stable")
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))

    xy = np.stack([ex[order], ey[order]], axis=1).astype(np.float32)
    c_lineage, c_unit = c_lineage[order], c_unit[order]
    c_visit, c_patient, c_type = c_visit[order], c_patient[order], c_type[order]
    c_tissue = c_tissue[order]
    c_prof, c_prob = c_prof[order], c_prob[order]
    c_mana, c_mscore = c_mana[order], c_mscore[order]

    n_tracked = int((c_lineage >= 0).sum())
    lineage_offset = np.zeros(len(lin_key) + 1, dtype=np.int64)
    if n_tracked:
        counts = np.bincount(c_lineage[:n_tracked], minlength=len(lin_key))
        lineage_offset[1:] = np.cumsum(counts)
    assert int(lineage_offset[-1]) == n_tracked, "lineage slices do not tile the cells"

    # -- 8. clone-unit arrays --
    cu_lin = np.array([lin_id[(p, r)] for p, r in zip(cu_patient, cu_receptor)],
                      dtype=np.int32)
    cu_arr = dict(
        xy=np.stack([cx, cy], axis=1).astype(np.float32),
        profile=np.array([clusters[c] for c in cids], dtype=np.int16),
        lineage=cu_lin,
        visit=np.array([visit_rank[v] for v in cu_visit], dtype=np.int8),
        patient=np.array([pat_index[p] for p in cu_patient], dtype=np.int16),
        n_cells=np.bincount(c_unit[c_unit >= 0], minlength=len(cids)).astype(np.int32),
        ptrc=u_ptrc, mana_frac=u_frac, mana_usable=u_usable,
    )
    layout_meta = [dict(id="default", label="UMAP",
                        extent=extent(cx, cy), source=emb_path)]
    for i, (lbl, (lx, ly)) in enumerate(extra_layouts.items()):
        cu_arr[f"xy_alt{i}"] = np.stack([lx, ly], axis=1).astype(np.float32)
        layout_meta.append(dict(id=f"alt{i}", label=lbl, extent=extent(lx, ly),
                                source=(m.get("clone_layouts") or {})[lbl]))

    # -- 9. patient table: clones and cells per visit --
    patients = []
    for pi, p in enumerate(patient_name):
        pm = c_patient == pi
        per_visit = {}
        for vi, v in enumerate(visits):
            cells_here = int((pm & (c_visit == vi)).sum())
            clones_here = int(((cu_arr["patient"] == pi) & (cu_arr["visit"] == vi)).sum())
            lin_here = int(len(np.unique(
                c_lineage[pm & (c_visit == vi) & (c_lineage >= 0)])))
            per_visit[v] = dict(cells=cells_here, clone_units=clones_here,
                                lineages=lin_here)
        patients.append(dict(
            id=p, fields=dict(patient_fields.get(p, {})),
            n_cells=int(pm.sum()),
            n_clone_units=int((cu_arr["patient"] == pi).sum()),
            n_lineages=len(by_patient.get(p, [])),
            per_visit=per_visit,
        ))

    # -- 10. write --
    out = os.path.join(out_root, m["name"])
    os.makedirs(out, exist_ok=True)
    # `source_row` is cell i's row in the cell embedding file. Cell barcodes are not
    # stored -- 150k-273k strings would dominate the bundle and the viewer never draws
    # them -- but this makes the identity recoverable, which is what lets the bundle be
    # checked against its sources and lets a selection be exported later.
    np.savez(os.path.join(out, "cells.npz"), xy=xy, lineage=c_lineage, clone_unit=c_unit,
             visit=c_visit, patient=c_patient, cell_type=c_type, profile=c_prof,
             tissue=c_tissue,
             prof_prob=c_prob, lineage_offset=lineage_offset,
             source_row=order.astype(np.int32),
             manascore_hi=c_mana, manascore=c_mscore)
    np.savez(os.path.join(out, "clones.npz"), **cu_arr)
    with open(os.path.join(out, "names.json"), "w") as fh:
        json.dump(dict(
            patients=patient_name, cell_types=cell_type_name,
            lineages=[dict(label=lin_label[i], patient=pat_index[lin_key[i][0]],
                           receptor=lin_key[i][1], n_cells=int(lin_cells[i]))
                      for i in range(len(lin_key))],
            clone_units=cids,
        ), fh)

    run_meta = {}
    rm = os.path.join(m["run_dir"], "run_meta.json")
    if os.path.exists(rm):
        with open(rm) as fh:
            run_meta = json.load(fh)

    prof_counts = Counter(int(v) for v in cu_arr["profile"])
    meta = dict(
        name=m["name"], label=m["label"], kind=m.get("kind", "sequence"),
        notes=m.get("notes", ""),
        built=time.strftime("%Y-%m-%d %H:%M:%S"),
        expanded_min=run_meta.get("expanded_min"),
        run_meta=run_meta,
        visits=[dict(code=v, label=(m.get("visit_labels") or {}).get(v, v),
                     n_cells=int((c_visit == i).sum()),
                     n_clone_units=int((cu_arr["visit"] == i).sum()))
                for i, v in enumerate(visits)],
        profiles=[dict(id=k, n_clone_units=prof_counts[k],
                       n_cells=int((c_prof == k).sum()))
                  for k in sorted(prof_counts)],
        cell_types=[dict(name=t, n_cells=int((c_type == i).sum()))
                    for i, t in enumerate(cell_type_name)],
        patients=patients,
        extent=dict(cell=extent(xy[:, 0], xy[:, 1]), clone=extent(cx, cy)),
        clone_layouts=layout_meta,
        tissues=[dict(code=t, label=t, n_cells=int((c_tissue == i).sum()))
                 for i, t in enumerate(tissue_name)],
        # WHICH (visit, tissue) PAIRS EXIST. In this cohort most do not -- four of the
        # seven sites were sampled at one visit only -- so a UI that offered the full
        # grid would offer mostly empty panels. Shipping the occupancy means the controls
        # can grey out what was never sampled instead of drawing a blank and letting the
        # reader infer "no cells here" when the truth is "never looked".
        visit_tissue=[[int(((c_visit == vi) & (c_tissue == ti)).sum())
                       for ti in range(len(tissue_name))]
                      for vi in range(len(visits))],
        counts=dict(cells=n_cells, cells_on_lineage=n_tracked,
                    cells_in_clone_unit=int((c_unit >= 0).sum()),
                    cells_below_threshold=int(((c_lineage >= 0) & (c_unit < 0)).sum()),
                    clone_units=len(cids), lineages=len(lin_key),
                    patients=len(patient_name)),
        sources=dict(run_dir=m["run_dir"], clone_embedding=emb_path,
                     cell_umap=m["cell_umap"], cell_clone=cc["path"],
                     cell_annotation=(ann or {}).get("path")),
        profile_read=prof_stats,
        manascore=mana_stats,
    )
    with open(os.path.join(out, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    say(f"  wrote {out} in {time.time() - t0:.1f}s")
    return out, meta


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    known = sorted(mf.available())
    ap.add_argument("dataset",
                    help="path to a manifest JSON, or a name on CLONOVIEW_MANIFESTS"
                         + (" (%s)" % ", ".join(known) if known else ""))
    ap.add_argument("--out", default=mf.data_root(),
                    help="where bundles are written (default: %(default)s)")
    a = ap.parse_args(argv)
    m = mf.load(a.dataset)
    print(f"preparing {m['label']} ({m['name']})", flush=True)
    out, meta = build(m, a.out)
    c = meta["counts"]
    print(f"\n{meta['label']}: {c['clone_units']} clone units over {c['lineages']} "
          f"lineages, {c['patients']} patients, {c['cells']} cells "
          f"({c['cells_in_clone_unit']} in a clone unit, {c['cells_below_threshold']} on "
          f"a tracked lineage but below the {meta['expanded_min']}-cell threshold)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
