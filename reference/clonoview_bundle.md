# clonoview_bundle

Prepare a Clonotrace run for viewing in ClonoView.

## Usage

``` r
clonoview_bundle(
  cell_embedding,
  cell_meta,
  clone_cluster,
  cell_profile_prob,
  name,
  label = name,
  clone_column = "clone",
  cell_type_column = NULL,
  patient_column = NULL,
  visit_column = NULL,
  patient_fields = NULL,
  clone_id_regex = NULL,
  patient_default = "sample",
  visits = NULL,
  visit_labels = NULL,
  visit_aliases = NULL,
  kind = NULL,
  notes = "",
  expanded_min = NULL,
  python = NULL,
  verbose = TRUE
)
```

## Arguments

- cell_embedding:

  Cell embedding: a data frame or matrix with cells in rows, identified
  by row names or a \`cell\` column, and two coordinate columns.

- cell_meta:

  Per-cell table, cells in rows, identified by row names or a \`cell\`
  column. Must contain \`clone_column\`.

- clone_cluster:

  Clone embedding and profile assignment, as returned by
  \[leiden_dis()\] with \`if_umap = TRUE\`: clones in rows, with
  coordinates and a \`cluster\` column.

- cell_profile_prob:

  Cell-by-profile probability matrix, with cells as row names and
  profile numbers as column names. Sparse matrices are accepted.

- name:

  Short identifier for the dataset; becomes the bundle's directory name.

- label:

  Human-readable name shown in the viewer. Defaults to \`name\`.

- clone_column:

  Column of \`cell_meta\` holding the clone assignment.

- cell_type_column:

  Column of \`cell_meta\` holding the cell-type annotation.

- patient_column, visit_column:

  Columns of \`cell_meta\` holding patient and visit.

- patient_fields:

  Named character vector of further per-patient columns to show, as
  \`c("Display name" = "column_name")\`.

- clone_id_regex:

  Python regular expression parsing the clone identifier, naming any of
  the groups \`patient\`, \`visit\`, \`receptor\`. See details.

- patient_default:

  Sample name used when clone identifiers carry no patient.

- visits:

  Visit codes in sampling order. Defaults to the sorted unique values of
  \`visit_column\`, or a single visit if there is none.

- visit_labels:

  Named character vector, visit code to display label.

- visit_aliases:

  Named character vector mapping the spelling used in \`cell_meta\` onto
  the code used in the clone identifier, as \`c("Baseline" = "B")\`.

- kind:

  \`"sequence"\` for a time course, \`"pre_post"\` for numbered biopsies
  that are not successive visits, \`"single"\` for one timepoint.
  Inferred if not given.

- notes:

  Free text shown in the viewer's dataset header.

- expanded_min:

  The cell threshold the run used to call a clone expanded, recorded so
  the viewer can say how many cells fell below it.

- python:

  Python interpreter to use. Defaults to \[clonoview_python()\].

- verbose:

  Whether to print progress.

## Value

The dataset name, invisibly.

## Details

Takes the objects a Clonotrace run leaves in the session, writes them
out in the layout ClonoView reads, and compiles them into a bundle. Run
it once per dataset; \[clonoview()\] then serves whatever has been
bundled.

The four required arguments are the four things the viewer draws:

- \`cell_embedding\`:

  the cell UMAP, as passed to \[dimplot()\]

- \`cell_meta\`:

  the per-cell table, which must carry the clone assignment

- \`clone_cluster\`:

  the \[leiden_dis()\] result – clone coordinates and the profile each
  clone belongs to

- \`cell_profile_prob\`:

  the cell-by-profile probability matrix

Everything else describes study design, and all of it is optional. With
none of it supplied the run is treated as a single sample with no time
dimension, which is the right reading of a lineage-barcoding experiment:
there is one library, and the viewer's time control correctly has one
step.

**Sample structure lives in the clone identifier.** A clone in a
longitudinal study is a clone *at one timepoint*, so which patient and
which visit a clone belongs to has to be recoverable from its name.
\`clone_id_regex\` is how, and it must name the groups it extracts: with
clone identifiers spelled \`p12\|Baseline\|CASSLGQAYEQYF\`,
\`clone_id_regex =
"^(?P\<patient\>\[^\|\]+)\\\|(?P\<visit\>\[^\|\]+)\\\|(?P\<receptor\>.+)\$"\`
splits them. Any of the three groups may be left out; omit the argument
entirely and the clone identifier is taken whole, as one lineage in one
sample.

The \`patient_column\` and \`visit_column\` arguments say where the
*cell* table records the same facts. They are what puts a cell in the
patient table, and they are independent of the regex: cells that carry
no clone still have a patient.

## See also

\[clonoview()\], \[clonoview_datasets()\]

## Examples

``` r
if (FALSE) { # \dontrun{
# after the clone embedding and profile steps of the vignette
clonoview_bundle(cell_embedding = umap,
                 cell_meta = cell_meta,
                 clone_cluster = clone_cluster,
                 cell_profile_prob = cell_profile_prob,
                 name = "hematopoiesis",
                 cell_type_column = "Cell.type.annotation")
clonoview("hematopoiesis")
} # }
```
