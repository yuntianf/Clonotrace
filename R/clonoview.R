## ClonoView -- the interactive viewer bundled with Clonotrace.
##
## The viewer itself is a small Python program plus a browser front end, shipped inside
## this package under inst/clonoview and found with system.file(). The functions here are
## the whole of the R interface to it:
##
##   clonoview_bundle()   turn the objects a Clonotrace run leaves in the session into a
##                        prepared bundle
##   clonoview()          serve a prepared bundle and open it
##   clonoview_stop()     stop the server
##   clonoview_datasets() list what has been prepared
##
## WHY A SUBPROCESS AND NOT reticulate. The rest of the package calls Python in-process
## with source_python(), which is right for umap_from_knn(): it returns a value. The
## viewer's server is a loop that does not return until it is stopped, so running it
## in-process would freeze the R session for as long as the viewer is open. It is
## launched as a background process instead. The interpreter is still the one reticulate
## is configured with, so a user who has already run use_python() for the co-embedding
## step has nothing further to set up.


# ---- locations ---------------------------------------------------------------

#' @title clonoview_home
#' @description Directory holding the bundled ClonoView source.
#' @details
#' `options(clonoview.home = ...)` overrides it, which is how the viewer is run from a
#' working copy rather than from the installed package.
#' @return A path, or "" if the package was installed without it.
#' @keywords internal
clonoview_home = function() {
  getOption("clonoview.home", system.file("clonoview", package = "Clonotrace"))
}

#' @title clonoview_cache
#' @description Where ClonoView writes prepared bundles, staged inputs and its pidfile.
#' @details
#'
#' Not inside the installed package: that directory is read-only on a managed R
#' installation and shared between users on a cluster. The default is the standard
#' per-user cache for this package, overridable with `options(clonoview.cache = ...)`.
#'
#' @param sub Optional sub-directory, created if it does not exist.
#' @return A path.
#' @export
clonoview_cache = function(sub = NULL) {
  root = getOption("clonoview.cache", tools::R_user_dir("Clonotrace", "cache"))
  p = if (is.null(sub)) root else file.path(root, sub)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  normalizePath(p, winslash = "/", mustWork = FALSE)
}


# ---- the interpreter ---------------------------------------------------------

#' @title clonoview_python
#' @description The Python interpreter used to run ClonoView.
#' @details
#'
#' In order: `options(clonoview.python)`, then the interpreter reticulate is configured
#' with -- so `use_python()` from the co-embedding step covers this too -- then whatever
#' `python3` is on the PATH.
#'
#' ClonoView needs numpy, scipy and matplotlib. numpy and scipy are already in the
#' package's `requirements.txt`; matplotlib is the one addition the viewer makes, and it
#' is only used to rasterise the cell panel's backdrop.
#'
#' @param check Whether to verify the interpreter can import what the viewer needs.
#' @return Path to a Python interpreter.
#' @export
clonoview_python = function(check = TRUE) {
  py = getOption("clonoview.python")
  if (is.null(py) && requireNamespace("reticulate", quietly = TRUE)) {
    py = tryCatch(reticulate::py_exe(), error = function(e) NULL)
  }
  if (is.null(py) || !nzchar(py)) {
    py = Sys.which("python3")
    if (!nzchar(py)) py = Sys.which("python")
  }
  if (!nzchar(py)) {
    stop("no Python interpreter found. Set one with reticulate::use_python(\"...\"), ",
         "or options(clonoview.python = \"...\").", call. = FALSE)
  }
  if (check) {
    miss = system2(py, c("-c", shQuote("import numpy, scipy, matplotlib")),
                   stdout = FALSE, stderr = FALSE)
    if (miss != 0) {
      stop("the Python at ", py, " cannot import numpy, scipy and matplotlib, which ",
           "ClonoView needs.\n  Install them there:  ", py,
           " -m pip install numpy scipy matplotlib\n  Or point at another ",
           "interpreter with options(clonoview.python = \"...\").", call. = FALSE)
    }
  }
  unname(py)
}

# Run the bundled Python. PYTHONPATH rather than a working-directory change, because the
# package directory is not writable and the caller's working directory is not ours to
# move. Set through Sys.setenv rather than system2(env=), which is a no-op on Windows.
cv_python_run = function(args, wait = TRUE, stdout = "", stderr = "", py = NULL) {
  home = clonoview_home()
  if (!nzchar(home)) {
    stop("this Clonotrace installation has no inst/clonoview. Reinstall with ",
         "devtools::install_github(\"yuntianf/Clonotrace\").", call. = FALSE)
  }
  old = Sys.getenv(c("PYTHONPATH", "CLONOVIEW_DATA"), unset = NA)
  sep = .Platform$path.sep
  pp = Sys.getenv("PYTHONPATH")
  Sys.setenv(PYTHONPATH = if (nzchar(pp)) paste(home, pp, sep = sep) else home,
             CLONOVIEW_DATA = clonoview_cache("data"))
  on.exit({
    for (k in names(old)) {
      if (is.na(old[[k]])) Sys.unsetenv(k) else do.call(Sys.setenv, as.list(old[k]))
    }
  }, add = TRUE)
  system2(py %||% clonoview_python(check = FALSE), args,
          wait = wait, stdout = stdout, stderr = stderr)
}

`%||%` = function(a, b) if (is.null(a)) b else a


# ---- writing the run out -----------------------------------------------------

# Identifiers go into tab-separated files that the Python side reads with a plain CSV
# reader, so a tab or a newline inside one would shift every following field silently.
# A leading double quote would be read as an opening quote for the same reason. None of
# these occurs in a receptor, a barcode or a cell name, so hitting one means something
# upstream is wrong and it is worth stopping for.
cv_check_ids = function(x, what) {
  x = as.character(x)
  bad = grepl("[\t\r\n]", x) | startsWith(x, "\"")
  if (any(bad)) {
    stop(sum(bad), " ", what, " contain a tab, a newline or a leading quote, which ",
         "would corrupt the tab-separated file ClonoView reads. First: ",
         paste(sQuote(utils::head(x[bad], 3)), collapse = ", "), call. = FALSE)
  }
  x
}

cv_write_tsv = function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = TRUE, na = "NA")
  path
}

# Pull (id, x, y) out of whatever an embedding is held in: a data.frame or matrix, with
# the identifier either in a named column or in the row names, and the coordinates named
# any of the several ways the package's own functions produce them -- leiden_dis()
# returns umap1/umap2, a Seurat reduction gives UMAP_1/UMAP_2.
cv_xy = function(x, id_col, what) {
  df = as.data.frame(x, stringsAsFactors = FALSE)
  ids = if (id_col %in% names(df)) df[[id_col]] else rownames(df)
  if (is.null(ids)) {
    stop(what, " has neither a ", sQuote(id_col), " column nor row names, so there is ",
         "no way to tell which row is which.", call. = FALSE)
  }
  norm = tolower(gsub("[^A-Za-z0-9]", "", names(df)))
  ix = match(c("umap1", "x"), norm); ix = ix[!is.na(ix)]
  iy = match(c("umap2", "y"), norm); iy = iy[!is.na(iy)]
  if (!length(ix) || !length(iy)) {
    num = which(vapply(df, is.numeric, logical(1)))
    if (length(num) < 2) {
      stop(what, " has no pair of coordinate columns; looked for umap1/umap2 or x/y ",
           "and found columns ", paste(sQuote(names(df)), collapse = ", "), call. = FALSE)
    }
    ix = num[1]; iy = num[2]
  }
  out = data.frame(id = cv_check_ids(ids, paste(what, "identifiers")),
                   umap_1 = as.numeric(df[[ix[1]]]),
                   umap_2 = as.numeric(df[[iy[1]]]),
                   stringsAsFactors = FALSE)
  names(out)[1] = id_col
  keep = is.finite(out$umap_1) & is.finite(out$umap_2)
  if (!all(keep)) {
    warning(sum(!keep), " of ", nrow(out), " ", what, " rows have a missing coordinate ",
            "and are dropped.", call. = FALSE)
  }
  out[keep, , drop = FALSE]
}

# The per-cell profile probabilities, as the X1..XK table the viewer reads.
#
# The column names carry the profile numbers and are read by name, not by position: the
# argmax is taken over columns sorted by the integer in the name. That is deliberate on
# the Python side, because these matrices do not always come out in numeric order, and
# reading them positionally would relabel profiles across the whole dataset.
cv_profile_table = function(m, cells) {
  cn = colnames(m)
  if (is.null(cn)) {
    stop("cell_profile_prob has no column names, so there is no way to know which ",
         "profile each column is.", call. = FALSE)
  }
  k = sub("^X", "", cn)
  if (!all(grepl("^[0-9]+$", k))) {
    stop("cell_profile_prob column names must be profile numbers (\"1\", \"2\", ...) ",
         "or the same with an X prefix; got ",
         paste(sQuote(utils::head(cn, 6)), collapse = ", "), call. = FALSE)
  }
  rn = rownames(m)
  if (is.null(rn)) {
    stop("cell_profile_prob has no row names, so its rows cannot be matched to cells.",
         call. = FALSE)
  }
  unknown = setdiff(rn, cells)
  if (length(unknown)) {
    warning(length(unknown), " of ", nrow(m), " cell_profile_prob rows name cells that ",
            "are not in the cell embedding; they are ignored.", call. = FALSE)
  }
  out = as.data.frame(as.matrix(m), stringsAsFactors = FALSE)
  colnames(out) = paste0("X", k)
  cbind(cell = cv_check_ids(rn, "cell names"), out, stringsAsFactors = FALSE)
}


# ---- bundling ----------------------------------------------------------------

#' @title clonoview_bundle
#' @description Prepare a Clonotrace run for viewing in ClonoView.
#' @details
#'
#' Takes the objects a Clonotrace run leaves in the session, writes them out in the
#' layout ClonoView reads, and compiles them into a bundle. Run it once per dataset;
#' [clonoview()] then serves whatever has been bundled.
#'
#' The four required arguments are the four things the viewer draws:
#'
#' \describe{
#'   \item{`cell_embedding`}{the cell UMAP, as passed to [dimplot()]}
#'   \item{`cell_meta`}{the per-cell table, which must carry the clone assignment}
#'   \item{`clone_cluster`}{the [leiden_dis()] result -- clone coordinates and the
#'     profile each clone belongs to}
#'   \item{`cell_profile_prob`}{the cell-by-profile probability matrix}
#' }
#'
#' Everything else describes study design, and all of it is optional. With none of it
#' supplied the run is treated as a single sample with no time dimension, which is the
#' right reading of a lineage-barcoding experiment: there is one library, and the
#' viewer's time control correctly has one step.
#'
#' \strong{Sample structure lives in the clone identifier.} A clone in a longitudinal
#' study is a clone \emph{at one timepoint}, so which patient and which visit a clone
#' belongs to has to be recoverable from its name. `clone_id_regex` is how, and it must
#' name the groups it extracts: with clone identifiers spelled `p12|Baseline|CASSLGQAYEQYF`,
#' `clone_id_regex = "^(?P<patient>[^|]+)\\\\|(?P<visit>[^|]+)\\\\|(?P<receptor>.+)$"`
#' splits them. Any of the three groups may be left out; omit the argument entirely and
#' the clone identifier is taken whole, as one lineage in one sample.
#'
#' The `patient_column` and `visit_column` arguments say where the \emph{cell} table
#' records the same facts. They are what puts a cell in the patient table, and they are
#' independent of the regex: cells that carry no clone still have a patient.
#'
#' @param cell_embedding Cell embedding: a data frame or matrix with cells in rows,
#'   identified by row names or a `cell` column, and two coordinate columns.
#' @param cell_meta Per-cell table, cells in rows, identified by row names or a `cell`
#'   column. Must contain `clone_column`.
#' @param clone_cluster Clone embedding and profile assignment, as returned by
#'   [leiden_dis()] with `if_umap = TRUE`: clones in rows, with coordinates and a
#'   `cluster` column.
#' @param cell_profile_prob Cell-by-profile probability matrix, with cells as row names
#'   and profile numbers as column names. Sparse matrices are accepted.
#' @param name Short identifier for the dataset; becomes the bundle's directory name.
#' @param label Human-readable name shown in the viewer. Defaults to `name`.
#' @param clone_column Column of `cell_meta` holding the clone assignment.
#' @param cell_type_column Column of `cell_meta` holding the cell-type annotation.
#' @param patient_column,visit_column Columns of `cell_meta` holding patient and visit.
#' @param patient_fields Named character vector of further per-patient columns to show,
#'   as `c("Display name" = "column_name")`.
#' @param clone_id_regex Python regular expression parsing the clone identifier, naming
#'   any of the groups `patient`, `visit`, `receptor`. See details.
#' @param patient_default Sample name used when clone identifiers carry no patient.
#' @param visits Visit codes in sampling order. Defaults to the sorted unique values of
#'   `visit_column`, or a single visit if there is none.
#' @param visit_labels Named character vector, visit code to display label.
#' @param visit_aliases Named character vector mapping the spelling used in `cell_meta`
#'   onto the code used in the clone identifier, as `c("Baseline" = "B")`.
#' @param kind `"sequence"` for a time course, `"pre_post"` for numbered biopsies that
#'   are not successive visits, `"single"` for one timepoint. Inferred if not given.
#' @param notes Free text shown in the viewer's dataset header.
#' @param expanded_min The cell threshold the run used to call a clone expanded, recorded
#'   so the viewer can say how many cells fell below it.
#' @param python Python interpreter to use. Defaults to [clonoview_python()].
#' @param verbose Whether to print progress.
#'
#' @return The dataset name, invisibly.
#'
#' @examples
#' \dontrun{
#' # after the clone embedding and profile steps of the vignette
#' clonoview_bundle(cell_embedding = umap,
#'                  cell_meta = cell_meta,
#'                  clone_cluster = clone_cluster,
#'                  cell_profile_prob = cell_profile_prob,
#'                  name = "hematopoiesis",
#'                  cell_type_column = "Cell.type.annotation")
#' clonoview("hematopoiesis")
#' }
#'
#' @seealso [clonoview()], [clonoview_datasets()]
#' @export
clonoview_bundle = function(cell_embedding, cell_meta, clone_cluster, cell_profile_prob,
                            name, label = name,
                            clone_column = "clone",
                            cell_type_column = NULL,
                            patient_column = NULL, visit_column = NULL,
                            patient_fields = NULL,
                            clone_id_regex = NULL, patient_default = "sample",
                            visits = NULL, visit_labels = NULL, visit_aliases = NULL,
                            kind = NULL, notes = "", expanded_min = NULL,
                            python = NULL, verbose = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("clonoview_bundle() needs the jsonlite package.", call. = FALSE)
  }
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", name)) {
    stop("`name` becomes a directory name, so it must start with a letter or digit and ",
         "contain only letters, digits, dot, dash and underscore; got ", sQuote(name),
         call. = FALSE)
  }
  say = function(...) if (verbose) message(...)

  # -- the two embeddings --
  cells = cv_xy(cell_embedding, "cell", "cell_embedding")
  clones = cv_xy(clone_cluster, "clone", "clone_cluster")
  say("staging ", nrow(cells), " cells and ", nrow(clones), " clone units")

  # -- clone -> profile --
  cc_df = as.data.frame(clone_cluster, stringsAsFactors = FALSE)
  if (!"cluster" %in% names(cc_df)) {
    stop("clone_cluster has no `cluster` column. It should be the result of ",
         "leiden_dis(dismat, if_umap = TRUE), which returns coordinates plus the ",
         "profile each clone was assigned to.", call. = FALSE)
  }
  cl_ids = if ("clone" %in% names(cc_df)) as.character(cc_df$clone) else rownames(cc_df)
  clusters = data.frame(clone = cl_ids,
                        cluster = as.integer(as.character(cc_df$cluster)),
                        stringsAsFactors = FALSE)
  if (anyNA(clusters$cluster)) {
    stop("clone_cluster$cluster does not convert to whole numbers; ClonoView needs the ",
         "profile number, and a factor whose levels are not the numbers themselves ",
         "would be renumbered here.", call. = FALSE)
  }
  clusters = clusters[clusters$clone %in% clones$clone, , drop = FALSE]
  if (!setequal(clusters$clone, clones$clone)) {
    stop("the clone coordinates and the profile assignment do not cover the same ",
         "clones; ", length(setdiff(clones$clone, clusters$clone)),
         " have a position but no profile.", call. = FALSE)
  }

  # -- per-cell table --
  meta = as.data.frame(cell_meta, stringsAsFactors = FALSE)
  cell_ids = if ("cell" %in% names(meta)) as.character(meta$cell) else rownames(meta)
  if (is.null(cell_ids)) {
    stop("cell_meta has neither a `cell` column nor row names.", call. = FALSE)
  }
  if (!clone_column %in% names(meta)) {
    stop("cell_meta has no column ", sQuote(clone_column), "; it holds the clone each ",
         "cell was assigned to. Columns present: ",
         paste(sQuote(utils::head(names(meta), 12)), collapse = ", "), call. = FALSE)
  }
  if (length(patient_fields) &&
      (is.null(names(patient_fields)) || !all(nzchar(names(patient_fields))))) {
    stop("`patient_fields` must be named, as c(\"Response\" = \"response.phenotype\"): ",
         "the name is the label shown in the viewer, the value is the cell_meta column.",
         call. = FALSE)
  }
  need = c(cell_type = cell_type_column, patient = patient_column, visit = visit_column,
           patient_fields)
  for (col in need) {
    if (!col %in% names(meta)) {
      stop("cell_meta has no column ", sQuote(col), ".", call. = FALSE)
    }
  }

  ann = data.frame(cell = cv_check_ids(cell_ids, "cell names"), stringsAsFactors = FALSE)
  for (nm in names(need)) ann[[nm]] = as.character(meta[[need[[nm]]]])
  ann = ann[ann$cell %in% cells$cell, , drop = FALSE]

  cell_clone = data.frame(cell = as.character(cell_ids),
                          clone = as.character(meta[[clone_column]]),
                          stringsAsFactors = FALSE)
  cell_clone = cell_clone[!is.na(cell_clone$clone) & nzchar(cell_clone$clone) &
                            cell_clone$clone != "NA", , drop = FALSE]
  cv_check_ids(cell_clone$clone, "clone names")
  say("  ", nrow(cell_clone), " of ", nrow(meta), " cells carry a clone")
  orphan = setdiff(clones$clone, cell_clone$clone)
  if (length(orphan)) {
    warning(length(orphan), " clone units are on the clone map but named by no cell in ",
            sQuote(clone_column), ", so they will have no cells in the viewer. First: ",
            paste(sQuote(utils::head(orphan, 3)), collapse = ", "), call. = FALSE)
  }

  # -- study design --
  # A patient column on the cells with no patient group in the clone identifier is the
  # one combination that cannot work: the viewer keys lineages by (patient, receptor)
  # taken from the clone identifier, so every lineage would land on the synthetic
  # single sample while the cells sat under their real patients, and the patient table
  # would come out empty. Better to say so here than to let prepare fail on the join.
  if (!is.null(patient_column) && is.null(clone_id_regex)) {
    stop("`patient_column` says which patient each cell belongs to, but without a ",
         "`clone_id_regex` naming a `patient` group there is no way to tell which ",
         "patient each clone belongs to, and the two would not join.\n  Either pass a ",
         "regex that extracts the patient from the clone identifier, or drop ",
         "`patient_column` and let the run be a single sample.", call. = FALSE)
  }
  if (is.null(visits)) {
    visits = if (!is.null(visit_column)) {
      v = unique(stats::na.omit(ann$visit)); sort(v[nzchar(v)])
    } else "all"
  }
  visits = as.character(visits)
  if (is.null(kind)) kind = if (length(visits) == 1) "single" else "sequence"
  if (length(visits) > 1 && !grepl("visit", clone_id_regex %||% "", fixed = TRUE)) {
    stop("`visits` names ", length(visits), " timepoints, but the clone identifier is ",
         "not parsed for a `visit` group, so every clone would be placed at the first ",
         "one.\n  Add a `visit` group to `clone_id_regex`, or leave `visits` unset.",
         call. = FALSE)
  }

  # -- write it out --
  run = clonoview_cache(file.path("runs", name))
  cv_write_tsv(clones, file.path(run, "clone_embedding.tsv"))
  cv_write_tsv(clusters, file.path(run, "clone_clusters.tsv"))
  cv_write_tsv(cv_profile_table(cell_profile_prob, cells$cell),
               file.path(run, "cell_profile_prob.tsv"))
  cv_write_tsv(cells, file.path(run, "cell_umap.tsv"))
  cv_write_tsv(cell_clone, file.path(run, "cell_clone.tsv"))
  cv_write_tsv(ann, file.path(run, "cell_meta.tsv"))
  writeLines(jsonlite::toJSON(list(expanded_min = expanded_min), auto_unbox = TRUE,
                              null = "null"),
             file.path(run, "run_meta.json"))

  m = list(name = name, label = label, kind = kind, notes = notes,
           run_dir = run,
           cell_umap = file.path(run, "cell_umap.tsv"),
           cell_clone = list(path = file.path(run, "cell_clone.tsv"),
                             cell = "cell", clone = "clone"),
           visits = as.list(visits))
  if (!is.null(clone_id_regex)) m$clone_id_regex = clone_id_regex
  if (!identical(patient_default, "sample")) m$patient_default = patient_default
  if (ncol(ann) > 1) {
    a = list(path = file.path(run, "cell_meta.tsv"), cell = "cell")
    if (!is.null(cell_type_column)) a$cell_type = "cell_type"
    if (!is.null(patient_column)) a$patient = "patient"
    if (!is.null(visit_column)) a$visit = "visit"
    # The staged file names these columns by their display names, so the mapping the
    # manifest carries is display -> the column just written under that name.
    if (length(patient_fields)) {
      a$patient_fields = stats::setNames(as.list(names(patient_fields)),
                                         names(patient_fields))
    }
    m$cell_annotation = a
  }
  if (length(visit_labels)) m$visit_labels = as.list(visit_labels)
  if (length(visit_aliases)) m$visit_aliases = as.list(visit_aliases)
  mf = file.path(run, "manifest.json")
  writeLines(jsonlite::toJSON(m, auto_unbox = TRUE, null = "null"), mf)

  # -- compile --
  say("compiling the bundle")
  st = cv_python_run(c("-m", "clonoview.prepare", shQuote(mf),
                       "--out", shQuote(clonoview_cache("data"))),
                     stdout = if (verbose) "" else FALSE,
                     stderr = "", py = python)
  if (st != 0) {
    stop("clonoview.prepare failed (exit ", st, "). The staged inputs it read are in ",
         run, call. = FALSE)
  }
  say("bundled. View it with:  clonoview(\"", name, "\")")
  invisible(name)
}


# ---- serving -----------------------------------------------------------------

#' @title clonoview_example
#' @description The worked example that ships with the package: a real Clonotrace run.
#' @details
#'
#' Returns a list of the four objects [clonoview_bundle()] takes, so it doubles as the
#' clearest statement of what those arguments have to look like:
#'
#' \describe{
#'   \item{`cell_embedding`}{150,150 x 2, cells in row names}
#'   \item{`cell_meta`}{150,150 rows with `clone`, `cell_type`, `patient`, `visit` and
#'     `Response`; 94,390 of the cells carry a clone and the rest are context}
#'   \item{`clone_cluster`}{916 clone units with coordinates and a `cluster` column}
#'   \item{`cell_profile_prob`}{150,150 x 10}
#' }
#'
#' The data is the NSCLC anti-PD-1 plus JAK1-inhibitor trial of Mathew et al., eight
#' patients sampled over five treatment cycles, 16 cell types. Clone identifiers are
#' spelled `pt12cy1 TRB:CALRQGAGELFF;TRA:CAVSAERGNMLTF` -- patient and cycle run together,
#' then a space, then the receptor -- which is what the `clone_id_regex` in the example
#' below reads.
#'
#' @return A named list of four objects.
#' @seealso [clonoview_demo()], [clonoview_bundle()]
#' @export
clonoview_example = function() {
  f = system.file("extdata", "mathew_nsclc.rds", package = "Clonotrace")
  if (!nzchar(f)) {
    stop("this Clonotrace installation has no example data.", call. = FALSE)
  }
  readRDS(f)
}

#' @title clonoview_demo
#' @description Bundle the shipped example run and open it in ClonoView.
#' @details
#'
#' One call, from a fresh install, to a viewer in the browser. It is the fastest way to
#' find out whether ClonoView works on a given machine, and it exercises the whole path
#' rather than a stored result: the example ships as ordinary R objects and is put through
#' [clonoview_bundle()] here exactly as your own run would be. Bundling takes a few
#' seconds; it is skipped on later calls unless `rebuild = TRUE`.
#'
#' @param open Whether to open a browser. See [clonoview()].
#' @param rebuild Rebuild the bundle even if it is already prepared.
#' @param python Python interpreter to use. Defaults to [clonoview_python()].
#' @param verbose Whether to print progress.
#' @return The URL, invisibly.
#'
#' @examples
#' \dontrun{
#' clonoview_demo()
#' }
#'
#' @seealso [clonoview_example()], [clonoview_bundle()], [clonoview_stop()]
#' @export
clonoview_demo = function(open = interactive(), rebuild = FALSE, python = NULL,
                          verbose = TRUE) {
  name = "mathew_nsclc_example"
  if (rebuild || !name %in% clonoview_datasets()) {
    d = clonoview_example()
    clonoview_bundle(
      cell_embedding = d$cell_embedding, cell_meta = d$cell_meta,
      clone_cluster = d$clone_cluster, cell_profile_prob = d$cell_profile_prob,
      name = name, label = "Mathew NSCLC (example)",
      cell_type_column = "cell_type",
      patient_column = "patient", visit_column = "visit",
      patient_fields = c(Response = "Response"),
      clone_id_regex = "^(?P<patient>pt\\d+)(?P<visit>cy\\d+) (?P<receptor>.+)$",
      visits = c("cy1", "cy2", "cy4", "cy5", "cy6"),
      visit_labels = c(cy1 = "Cycle 1", cy2 = "Cycle 2", cy4 = "Cycle 4",
                       cy5 = "Cycle 5", cy6 = "Cycle 6"),
      expanded_min = 10L,
      notes = paste("Example data shipped with Clonotrace. Patients received two cycles",
                    "of pembrolizumab, then itacitinib plus pembrolizumab from the start",
                    "of cycle 3, then pembrolizumab alone from cycle 5, so the cycle 2 to",
                    "cycle 4 step is the one that spans the JAK1 inhibitor. Sampling is",
                    "irregular: only pt25 was ever sampled at cycle 5."),
      python = python, verbose = verbose)
  } else if (verbose) {
    message("using the example bundle already prepared; rebuild = TRUE to redo it.")
  }
  clonoview(name, open = open, python = python)
}

#' @title clonoview_export
#' @description Build a self-contained copy of the viewer that opens without R or Python.
#' @details
#'
#' Writes a folder holding the viewer and the data it draws, and zips it. The person
#' receiving it double-clicks `ClonoView.html`; there is nothing to install, and it does
#' not use the internet. This is the form to send to a collaborator or a reviewer, or to
#' attach as supplementary material.
#'
#' \strong{The export contains the data, not a picture of it.} Per-cell coordinates, cell
#' types, subject identifiers and every clone identifier -- which for an immune-repertoire
#' run means the receptor sequences -- are all in the folder, because the viewer computes
#' from them in the browser. Treat sending one as publishing that data.
#'
#' @param file Path for the zip file to write. The folder is built beside it.
#' @param datasets Datasets to include. Defaults to everything prepared.
#' @param zip Whether to zip the folder. FALSE leaves just the folder.
#' @param python Python interpreter to use. Defaults to [clonoview_python()].
#' @param verbose Whether to print progress.
#' @return The path written, invisibly.
#' @seealso [clonoview_bundle()], [clonoview()]
#' @export
clonoview_export = function(file = "ClonoView-offline.zip", datasets = NULL, zip = TRUE,
                            python = NULL, verbose = TRUE) {
  have = clonoview_datasets()
  if (!length(have)) {
    stop("nothing has been prepared yet. Build a bundle first with clonoview_bundle().",
         call. = FALSE)
  }
  bad = setdiff(datasets, have)
  if (length(bad)) {
    stop("not prepared: ", paste(sQuote(bad), collapse = ", "), ". Prepared: ",
         paste(sQuote(have), collapse = ", "), call. = FALSE)
  }
  file = normalizePath(file, winslash = "/", mustWork = FALSE)
  out = file.path(dirname(file), sub("\\.zip$", "", basename(file)))
  args = c("-m", "clonoview.build_export",
           "--data", shQuote(clonoview_cache("data")),
           "--out", shQuote(out))
  args = c(args, if (isTRUE(zip)) c("--zip", shQuote(file)) else "--no-zip")
  if (length(datasets)) args = c(args, shQuote(datasets))
  st = cv_python_run(args, stdout = if (verbose) "" else FALSE, stderr = "", py = python)
  if (st != 0) stop("clonoview.build_export failed (exit ", st, ").", call. = FALSE)
  path = if (isTRUE(zip)) file else out
  if (verbose) message("wrote ", path)
  invisible(path)
}

#' @title clonoview_datasets
#' @description Datasets that have been prepared with [clonoview_bundle()].
#' @return A character vector of dataset names.
#' @export
clonoview_datasets = function() {
  d = clonoview_cache("data")
  n = list.dirs(d, full.names = FALSE, recursive = FALSE)
  n[vapply(n, function(x) file.exists(file.path(d, x, "meta.json")), logical(1))]
}

# The port is chosen here rather than by the server, so that this function knows the URL
# without having to parse it back out of the server's output.
cv_free_port = function(from = 8770, tries = 40) {
  for (p in seq.int(from, from + tries - 1L)) {
    s = try(serverSocket(p), silent = TRUE)
    if (!inherits(s, "try-error")) {
      close(s)
      return(p)
    }
  }
  stop("no free port in ", from, "..", from + tries - 1L, call. = FALSE)
}

cv_state = new.env(parent = emptyenv())

#' @title clonoview
#' @description Serve a prepared dataset in the ClonoView viewer.
#' @details
#'
#' Starts the viewer in the background and returns, so the R session stays usable while
#' the browser is open. Stop it with [clonoview_stop()].
#'
#' The server binds to the loopback address, which means it is reachable from the machine
#' it runs on and from nowhere else. On a laptop that is all that is needed and the
#' browser opens straight onto it. On a cluster the way in is an SSH tunnel, and the
#' command to open one is printed; binding to the node's own address instead would put
#' unpublished data on every host that can route to it, for no gain over the tunnel.
#'
#' @param dataset Dataset to open on load. With more than one prepared and none named,
#'   the viewer opens on a chooser.
#' @param port Port to serve on. Defaults to the first free port from 8770.
#' @param open Whether to open a browser. Defaults to TRUE in an interactive session on
#'   a machine that has one.
#' @param python Python interpreter to use. Defaults to [clonoview_python()].
#' @return The URL, invisibly.
#' @seealso [clonoview_bundle()], [clonoview_stop()]
#' @export
clonoview = function(dataset = NULL, port = NULL, open = interactive(), python = NULL) {
  have = clonoview_datasets()
  if (!length(have)) {
    stop("nothing has been prepared yet. Build a bundle first with ",
         "clonoview_bundle().", call. = FALSE)
  }
  if (!is.null(dataset) && !dataset %in% have) {
    stop("dataset ", sQuote(dataset), " is not prepared. Prepared: ",
         paste(sQuote(have), collapse = ", "), call. = FALSE)
  }
  if (!is.null(cv_state$pid) && cv_state$pid %in% cv_running_pids()) {
    stop("a ClonoView started from this session is already running on port ",
         cv_state$port, ". Stop it first with clonoview_stop().", call. = FALSE)
  }
  port = port %||% cv_free_port()
  py = python %||% clonoview_python()
  log = file.path(clonoview_cache(), "clonoview.log")

  args = c("-m", "clonoview.server", "--port", port,
           "--data", shQuote(clonoview_cache("data")))
  if (!is.null(dataset)) args = c(args, shQuote(dataset))
  cv_python_run(args, wait = FALSE, stdout = log, stderr = log, py = py)

  url = paste0("http://localhost:", port, "/")
  # The server writes its pidfile only once it has bound, so waiting for the file is
  # waiting for a server that is actually listening rather than one that is about to
  # fail on a busy port.
  pid = NULL
  for (i in 1:50) {
    Sys.sleep(0.1)
    pid = cv_pid()
    if (!is.null(pid)) break
  }
  if (is.null(pid)) {
    stop("the viewer did not start. Its output:\n",
         paste(utils::tail(readLines(log, warn = FALSE), 20), collapse = "\n"),
         call. = FALSE)
  }
  cv_state$pid = pid
  cv_state$port = port

  message("ClonoView is serving at ", url)
  message("  If R is running on a remote machine, forward the port from your own:\n",
          "    ssh -L ", port, ":localhost:", port, " <user>@<host>")
  message("  Stop it with clonoview_stop().")
  if (isTRUE(open)) try(utils::browseURL(url), silent = TRUE)
  invisible(url)
}

# The pidfile is per host: this cache directory is often on a networked home directory
# shared by every node of a cluster, and without the hostname in the name, starting the
# viewer on one node would overwrite another node's entry -- so clonoview_stop() would
# then signal that pid on the wrong machine, where it belongs to somebody else.
cv_pidfile = function() {
  # The SHORT host name, because that is what the server writes: Python takes
  # socket.gethostname().split(".")[0] while R's nodename is the fully qualified name,
  # so leaving the domain on builds a path the server never writes to. The symptom is
  # not a missing file, it is clonoview() waiting for a server that started correctly
  # and is already serving, then reporting that it did not start -- and
  # clonoview_stop() never finding it afterwards.
  host = sub("\\..*$", "", Sys.info()[["nodename"]])
  file.path(clonoview_cache(), paste0("clonoview.", host, ".pid"))
}

cv_pid = function() {
  f = cv_pidfile()
  if (!file.exists(f)) return(NULL)
  p = suppressWarnings(as.integer(readLines(f, warn = FALSE)[1]))
  if (is.na(p)) NULL else p
}

cv_running_pids = function() {
  if (.Platform$OS.type == "windows") return(integer(0))
  suppressWarnings(as.integer(system2("ps", c("-o", "pid=", "-A"), stdout = TRUE)))
}

#' @title clonoview_stop
#' @description Stop the running ClonoView server.
#' @details
#'
#' Signals one process, named by the pidfile the server wrote once it had bound. It does
#' not match on the command line: a pattern like `pkill -f clonoview` also matches the
#' process doing the matching, and process ids are reused, so on a shared machine that
#' can signal somebody else's job.
#'
#' @return TRUE if a server was stopped, FALSE if none was running, invisibly.
#' @export
clonoview_stop = function() {
  pid = cv_pid()
  if (is.null(pid)) {
    message("no ClonoView server is recorded as running on this machine.")
    return(invisible(FALSE))
  }
  # Confirm it is still ours before signalling. On Unix the command line is checkable;
  # elsewhere fall back to the pid this session recorded, and refuse a pid we did not
  # start ourselves.
  ours = if (.Platform$OS.type == "windows") {
    identical(pid, cv_state$pid)
  } else {
    cmd = suppressWarnings(system2("ps", c("-o", "command=", "-p", pid), stdout = TRUE,
                                   stderr = FALSE))
    length(cmd) > 0 && any(grepl("clonoview.server", cmd, fixed = TRUE))
  }
  if (!ours) {
    message("process ", pid, " is not a ClonoView server; leaving it alone.")
    unlink(cv_pidfile())
    return(invisible(FALSE))
  }
  tools::pskill(pid)
  unlink(cv_pidfile())
  cv_state$pid = NULL
  message("stopped ClonoView (pid ", pid, ").")
  invisible(TRUE)
}
