"""Dataset manifests -- how to find the five things ClonoView needs in a Clonotrace run.

A Clonotrace run always produces the same objects, but every dataset stages its inputs
differently: one study's per-cell annotation sits in the run's own input directory,
another's lives in a separate published file whose header has no name for the row-id
column, and each encodes whatever sample structure it has inside the clone identifier
with a different syntax. A manifest is the small declarative file that reconciles those
differences, so `prepare.py` has exactly one code path.

The five things:

  1. clone embedding      run_dir/clone_embedding.tsv, or an override (see `clone_embedding`)
  2. clone -> profile     run_dir/clone_clusters.tsv
  3. cell embedding       `cell_umap`
  4. cell -> clone        `cell_clone`
  5. per-cell annotation  `cell_annotation` -- cell type, patient, visit, patient-level fields

and, alongside them, run_dir/cell_profile_prob.tsv for the per-cell profile.

MOST USERS SHOULD NEVER WRITE ONE OF THESE BY HAND. `Clonotrace::clonoview_bundle()`
takes the R objects a Clonotrace run already leaves in the session -- the cell embedding,
`cell_meta`, the `leiden_dis` result, `cell_profile_prob` -- writes them out in the layout
above, and generates the manifest. This file is the contract underneath that, and the way
in for a run staged on disk by something other than the R wrapper.

## Sample structure is optional

`clone_id_regex` pulls patient, visit and receptor out of the clone identifier, and all
three of its named groups are optional:

  * A longitudinal immune-repertoire study uses all three. `receptor` is what makes a
    clone trackable through time: one T-cell lineage is (patient, receptor), and it
    appears as a separate *clone unit* at each visit where it cleared the run's
    `expanded_min` cell threshold.
  * A study with one sample per patient omits `visit`; `visits` must then name the single
    timepoint, and every clone unit belongs to it.
  * A lineage-barcoding experiment -- one mouse, one library, barcodes rather than
    receptors -- omits the field entirely. The whole clone identifier is then the
    receptor, every clone unit belongs to the one sample named by `patient_default`, and
    the viewer's time control has a single step, which is the honest rendering: nothing
    in such a run moves through time.

The clone identifier's visit codes are canonical. Where the annotation file spells the
same visits differently -- clone ids saying `B` and `F1` where the annotation says
`Baseline` and `Follow Up 1` -- `visit_aliases` maps the annotation's spelling onto the
canonical code. Without it the two sources silently fail to join and every cell looks
like it belongs to no visit, so `prepare.py` treats an unmapped annotation visit as a
hard error rather than dropping the cell.
"""

import json
import os
import re

# Where prepared bundles are written and served from.
#
# This must not be inside the package. ClonoView ships inside the Clonotrace R package,
# whose installed directory is read-only on any managed R installation and shared between
# users on a cluster; the earlier version of this file put `data/` next to the source and
# would fail on the first `prepare` after a normal install. The R wrapper passes
# `tools::R_user_dir("Clonotrace", "cache")`; this env var and the XDG default are for
# running the Python directly.
def data_root():
    """The directory holding prepared bundles, one sub-directory per dataset."""
    env = os.environ.get("CLONOVIEW_DATA")
    if env:
        return os.path.abspath(os.path.expanduser(env))
    cache = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"),
                                                             ".cache")
    return os.path.join(cache, "clonoview", "data")


# Every field a manifest may carry. Anything else is a typo and is rejected loudly --
# a silently ignored key here would mean a dataset quietly loses its cell types.
FIELDS = {
    "name",              # short identifier, used as the bundle directory name
    "label",             # human-readable dataset name
    "run_dir",           # Clonotrace output directory
    "clone_embedding",   # optional override for run_dir/clone_embedding.tsv
    "clone_layouts",     # optional {label: path} of ALTERNATIVE clone layouts, offered
                         # in the viewer beside the default one. Same clone ids, same
                         # profiles -- only the coordinates differ.
    "cell_umap",         # cell embedding tsv: cell, umap_1, umap_2
    "cell_clone",        # {path, cell, clone}
    "clone_id_regex",    # optional; may name groups patient, visit, receptor
    "patient_default",   # sample name used when clone_id_regex has no `patient` group
    "cell_annotation",   # {path, cell, cell_type, patient, visit, tissue,
                         #  patient_fields, header_offset, delim}
                         # `tissue` is an optional SECOND per-cell axis, anatomic
                         # position, treated in the viewer exactly as time is: vary one
                         # while holding the other fixed.
    "manascore",         # optional {cells:{...}, clones:{...}, note} -- an independent
                         # per-cell score overlaid on both panels
    "visits",            # sampling schedule, in order (canonical codes)
    "visit_aliases",     # raw annotation value -> canonical code, where they differ
    "visit_labels",      # visit code -> display label
    "kind",              # "sequence" (a time course), "pre_post", or "single"
    "notes",             # free text shown in the app's dataset header
}
REQUIRED = {"name", "label", "run_dir", "cell_umap", "cell_clone", "visits"}

# Used when a manifest gives no `clone_id_regex`: the clone identifier is opaque and is
# taken whole as the receptor. This is the lineage-barcoding case.
OPAQUE_ID = r"^(?P<receptor>.+)$"
ID_GROUPS = {"patient", "visit", "receptor"}


def clone_id_regex(m):
    """The compiled clone-identifier parser for a manifest."""
    return re.compile(m.get("clone_id_regex") or OPAQUE_ID)


def patient_default(m):
    """The sample name for a manifest whose clone ids carry no patient."""
    return m.get("patient_default") or "sample"


def check(m):
    """Validate a manifest dict, returning it. Raises on anything that would produce a
    silently wrong bundle rather than a crash."""
    unknown = set(m) - FIELDS
    if unknown:
        raise ValueError(f"manifest {m.get('name')!r}: unknown field(s) "
                         f"{sorted(unknown)}; known fields are {sorted(FIELDS)}")
    missing = REQUIRED - set(m)
    if missing:
        raise ValueError(f"manifest {m.get('name')!r}: missing required field(s) "
                         f"{sorted(missing)}")

    groups = set(clone_id_regex(m).groupindex)
    stray = groups - ID_GROUPS
    if stray:
        raise ValueError(f"manifest {m['name']!r}: clone_id_regex names group(s) "
                         f"{sorted(stray)}, which nothing reads; the groups it may "
                         f"define are {sorted(ID_GROUPS)}")
    if not groups:
        raise ValueError(f"manifest {m['name']!r}: clone_id_regex defines no named "
                         f"groups, so it identifies nothing; name at least `receptor`, "
                         f"or omit the field to take the whole clone id as the receptor")
    # `visit` is optional. A dataset can be a single snapshot -- one tumour resection per
    # patient, or one lineage-barcoding library -- in which case the clone identifier has
    # no timepoint in it and every clone unit belongs to the one visit named in `visits`.
    # Requiring the group would have meant inventing a fake timepoint inside the regex.
    if "visit" not in groups and len(m["visits"]) != 1:
        raise ValueError(f"manifest {m['name']!r}: clone_id_regex has no `visit` group, "
                         f"so the dataset is a single snapshot and `visits` must name "
                         f"exactly one visit; it names {m['visits']}")

    if not m["visits"]:
        raise ValueError(f"manifest {m['name']!r}: `visits` is empty")
    if len(set(m["visits"])) != len(m["visits"]):
        raise ValueError(f"manifest {m['name']!r}: `visits` has duplicates")

    kind = m.get("kind", "sequence")
    if kind not in ("sequence", "pre_post", "single"):
        raise ValueError(f"manifest {m['name']!r}: kind must be 'sequence', 'pre_post' "
                         f"or 'single', got {kind!r}")

    for key in ("run_dir", "cell_umap"):
        if not os.path.exists(m[key]):
            raise FileNotFoundError(f"manifest {m['name']!r}: {key} does not exist: "
                                    f"{m[key]}")
    for key in ("cell_clone", "cell_annotation"):
        sub = m.get(key)
        if sub and not os.path.exists(sub["path"]):
            raise FileNotFoundError(f"manifest {m['name']!r}: {key}.path does not "
                                    f"exist: {sub['path']}")
    lay = m.get("clone_layouts")
    if lay:
        if not isinstance(lay, dict):
            raise ValueError(f"manifest {m['name']!r}: clone_layouts must be a "
                             f"{{label: path}} mapping")
        for lbl, path in lay.items():
            if not os.path.exists(path):
                raise FileNotFoundError(f"manifest {m['name']!r}: clone_layouts["
                                        f"{lbl!r}] does not exist: {path}")

    ms = m.get("manascore")
    if ms:
        for part in ("cells", "clones"):
            if part not in ms:
                raise ValueError(f"manifest {m['name']!r}: manascore needs both `cells` "
                                 f"and `clones`; missing {part!r}")
            if not os.path.exists(ms[part]["path"]):
                raise FileNotFoundError(f"manifest {m['name']!r}: manascore.{part}.path "
                                        f"does not exist: {ms[part]['path']}")
    return m


def clone_embedding_path(m):
    """The clone UMAP to use.

    Defaults to the run's own `clone_embedding.tsv`. Worth knowing what that is: it
    comes from `Clonotrace::leiden_dis`, which calls the R umap package with its
    defaults, n_neighbors 15 / min_dist 0.1. On the datasets checked here the k-NN graph
    over the same distances stays a single connected component down to k=10, so the
    islands that setting produces are a layout artefact rather than structure in the
    data, and re-laying the same distance matrix out at n_neighbors 30 / min_dist 0.4
    gives a map that matches the graph. A manifest can therefore point `clone_embedding`
    at a re-laid-out file. Profiles are unaffected either way, because they come from
    Leiden on a separate SNN graph over the same distances -- this moves points only.
    """
    return m.get("clone_embedding") or os.path.join(m["run_dir"], "clone_embedding.tsv")


def search_path():
    """Directories searched for manifests named by bare name, in order.

    `CLONOVIEW_MANIFESTS` is a colon-separated list, so a lab with a set of staged runs
    can keep their manifests together and refer to them as `prepare mystudy` rather than
    by path. There are no manifests built into this module: it used to carry four with
    absolute paths into one cluster's filesystem, which meant every install of the
    package advertised datasets nobody else could read.
    """
    raw = os.environ.get("CLONOVIEW_MANIFESTS", "")
    return [os.path.abspath(os.path.expanduser(p)) for p in raw.split(os.pathsep) if p]


def available():
    """Dataset names findable on the search path, with the manifest file for each."""
    found = {}
    for d in search_path():
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if fn.endswith(".json"):
                found.setdefault(fn[:-len(".json")], os.path.join(d, fn))
    return found


def load(name_or_path):
    """A path to a manifest JSON file, or a bare name found on the search path."""
    if os.path.exists(name_or_path) and not os.path.isdir(name_or_path):
        with open(name_or_path) as fh:
            return check(json.load(fh))
    found = available()
    if name_or_path in found:
        with open(found[name_or_path]) as fh:
            return check(json.load(fh))
    known = sorted(found) or ["none -- set CLONOVIEW_MANIFESTS to a directory of them"]
    raise KeyError(f"no manifest file at {name_or_path!r} and no dataset of that name "
                   f"on the search path (found: {', '.join(known)})")
