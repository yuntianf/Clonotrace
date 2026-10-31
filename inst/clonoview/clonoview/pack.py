"""Pack a prepared bundle into JavaScript files a browser can load from `file://`.

WHY JAVASCRIPT AND NOT THE .npz. The export has to open by double-clicking an HTML file,
with no server and nothing installed. A page opened from `file://` is given an opaque
origin, and `fetch`, `XMLHttpRequest` and ES modules are all refused across it -- so the
viewer cannot read a data file the way the served version does. A classic
`<script src="...">` tag IS allowed from `file://` in Chrome, Edge and Firefox, and it is
the only loading mechanism that is. So every array is written into a .js file that assigns
itself into `window.CV_DATA`, and the page loads a dataset by appending a script tag.

WHAT IS DROPPED. `prepare.py` writes four per-cell arrays the viewer never reads:
`clone_unit`, `prof_prob`, `source_row` and the float `manascore`. They are provenance and
were kept in the served bundle because they cost nothing there. Here they are 12 bytes a
cell of pure download, so they are left out. The served bundle keeps them; this is an
export format, not a replacement.

WHAT IS NOT DROPPED: PRECISION. Coordinates are written as float32, exactly as they are
stored, and not quantised to int16. Quantising would have saved about 2 MB across the four
datasets, and at 0.0003 data units against 0.01 units per screen pixel it would have been
invisible on screen. It is still not done, because the JavaScript port is verified by
comparing it against this Python (see `selftest.html`), and a lossy step in the data would
put a floor under every comparison -- a real disagreement of a third of a pixel would then
be indistinguishable from the packing. Exactness here is worth more than the megabytes: it
means any deviation the self-check reports is a bug, with no threshold to argue about.
Integer columns ARE narrowed, but only where the narrowing is lossless and asserted so.
"""

import base64
import json
import os

import numpy as np

from . import bundle as bd
from . import colour as C
from . import manifest as mf

HERE = os.path.dirname(os.path.abspath(__file__))
VIZ_ROOT = os.path.dirname(HERE)

# numpy's native byte order is little-endian on x86 and on ARM as configured by every
# platform a browser runs on, and JavaScript's typed arrays use the platform's order with
# no way to ask for one. Rather than rely on both being little-endian by luck, the arrays
# are written little-endian explicitly and the loader asserts the platform matches.
DTYPES = {
    "f4": "<f4", "i4": "<i4", "i2": "<i2", "i1": "|i1", "u1": "|u1",
}

JS_VIEW = {
    "f4": "Float32Array", "i4": "Int32Array", "i2": "Int16Array",
    "i1": "Int8Array", "u1": "Uint8Array",
}


def _narrow(a, want):
    """Downcast `a` to `want`, refusing to do it if a single value would change.

    The narrowing is worth about 4 bytes a cell, which is 3 MB across the datasets, but
    it is exactly the kind of change that silently corrupts one dataset and not the
    others -- Wang has 703 lineages and would fit in an int16 whether or not the code
    were correct, while a dataset with 40,000 would wrap to negative and quietly relabel
    every cell. The assertion is what makes it safe to apply the same rule everywhere.
    """
    out = a.astype(DTYPES[want])
    if not np.array_equal(out.astype(np.int64), a.astype(np.int64)):
        raise ValueError(f"narrowing to {want} would change values "
                         f"(range {a.min()}..{a.max()})")
    return out


def _arr(a, kind):
    a = np.ascontiguousarray(a)
    if kind == "f4":
        a = a.astype(DTYPES["f4"])
    else:
        a = _narrow(a, kind)
    return dict(
        t=kind,
        n=int(a.shape[0]),
        cols=int(a.shape[1]) if a.ndim > 1 else 1,
        b64=base64.b64encode(a.tobytes()).decode("ascii"),
    )


def pack_dataset(root, name, out_dir):
    """Write `<out_dir>/<name>.js`. Returns a small record for the build log."""
    b = bd.load(root, name)
    n_cells = len(b.cxy)

    # The colours are resolved HERE, not in the browser. colour.py searches for a
    # cell-type palette by measuring perceptual distance under normal and colour-vision-
    # deficient simulation, over a cached set of candidates -- hundreds of lines that
    # exist to produce, for each dataset, one list of hex strings. Porting that search to
    # JavaScript would be pure risk for no gain, since the answer is fixed once the class
    # count is known. The measured quality is carried across with it, because the viewer
    # states those numbers rather than implying a separation it does not have.
    pal = b.cell_type_colours
    nprof = len(b.profile_ids)
    q = C.measure(C.PROFILE_COLOURS[:nprof])
    payload = dict(
        name=name,
        meta=b.meta,
        names=b.names,
        colours=dict(
            profile={str(k): v for k, v in b.profile_colours.items()},
            cell_type={str(i): pal[i] for i in range(len(pal))},
            cell_type_palette_quality=C.cell_type_palette_quality(len(b.cell_types)),
            profile_palette_quality=dict(worst_cvd=round(q["worst_cvd"], 2),
                                         worst_normal=round(q["worst_normal"], 2),
                                         low_contrast=q["low_contrast"]),
            ui=dict(surface=C.SURFACE, ink=C.INK, ink2=C.INK_2, muted=C.MUTED,
                    context=C.CONTEXT, grid=C.GRID, time=C.SEQ_TIME, mana=C.MANA_RED),
        ),
        arrays=dict(
            cxy=_arr(b.cxy, "f4"),
            c_lineage=_arr(b.c_lineage, "i2"),
            c_visit=_arr(b.c_visit, "i1"),
            c_patient=_arr(b.c_patient, "i1"),
            c_type=_arr(b.c_type, "i1"),
            c_profile=_arr(b.c_profile, "i1"),
            c_tissue=_arr(b.c_tissue, "i1"),
            c_mana=_arr(b.c_mana, "i1"),
            lin_off=_arr(b.lin_off, "i4"),
            uxy=_arr(b.uxy, "f4"),
            # Alternative clone layouts, same clone order as uxy. Packed under the ids
            # the meta advertises (`uxy_alt0`, ...) so the offline api.js can select one
            # exactly as the server does.
            **{f"uxy_{L['id']}": _arr(b.uxy_for(L["id"]), "f4")
               for L in b.layouts() if L["id"] != "default"},
            u_profile=_arr(b.u_profile, "i2"),
            u_lineage=_arr(b.u_lineage, "i2"),
            u_visit=_arr(b.u_visit, "i1"),
            u_patient=_arr(b.u_patient, "i1"),
            u_n=_arr(b.u_n, "i4"),
            u_ptrc=_arr(b.u_ptrc, "i1"),
            u_mana_usable=_arr(b.u_mana_usable, "i1"),
        ),
        has_manascore=bool(b.has_manascore),
        n_cells=int(n_cells),
        n_units=int(len(b.uxy)),
    )

    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{name}.js")
    # separators: json.dumps' default puts a space after every ':' and ',', which on a
    # file with 700k-element arrays of its own is megabytes of spaces.
    body = json.dumps(payload, separators=(",", ":"))
    with open(path, "w") as fh:
        fh.write("window.CV_DATA=window.CV_DATA||{};\n")
        fh.write(f"window.CV_DATA[{json.dumps(name)}]=")
        fh.write(body)
        fh.write(";\n")
    return dict(name=name, path=path, bytes=os.path.getsize(path),
                cells=n_cells, units=len(b.uxy))


def pack_index(root, out_dir, names):
    """The catalogue, which is what the page loads at startup.

    Kept separate from the datasets and tiny on purpose: the picker has to be populated
    before anything is chosen, and pulling a 4 MB dataset in to find out its label would
    make the first paint wait on the largest file in the export.
    """
    cat = [d for d in bd.available(root) if d["name"] in names]
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "index.js")
    with open(path, "w") as fh:
        fh.write("window.CV_INDEX=")
        fh.write(json.dumps(dict(datasets=cat, views=JS_VIEW),
                            separators=(",", ":")))
        fh.write(";\n")
    return path


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Pack bundles for the offline export.")
    ap.add_argument("names", nargs="*", help="datasets to pack (default: all prepared)")
    ap.add_argument("--data", default=mf.data_root(),
                    help="directory of prepared bundles (default: %(default)s)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(mf.data_root()),
                                                  "export", "data"),
                    help="where to write the packed data (default: %(default)s)")
    a = ap.parse_args(argv)

    names = a.names or [d["name"] for d in bd.available(a.data)]
    total = 0
    for n in names:
        r = pack_dataset(a.data, n, a.out)
        total += r["bytes"]
        print(f"  {r['name']:8s} {r['cells']:>8,} cells  {r['units']:>6,} units  "
              f"{r['bytes'] / 1e6:6.2f} MB")
    pack_index(a.data, a.out, set(names))
    print(f"  {'total':8s} {'':>8s}        {'':>6s}        {total / 1e6:6.2f} MB")


if __name__ == "__main__":
    main()
