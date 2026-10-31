# clonoview_export

Build a self-contained copy of the viewer that opens without R or
Python.

## Usage

``` r
clonoview_export(
  file = "ClonoView-offline.zip",
  datasets = NULL,
  zip = TRUE,
  python = NULL,
  verbose = TRUE
)
```

## Arguments

- file:

  Path for the zip file to write. The folder is built beside it.

- datasets:

  Datasets to include. Defaults to everything prepared.

- zip:

  Whether to zip the folder. FALSE leaves just the folder.

- python:

  Python interpreter to use. Defaults to \[clonoview_python()\].

- verbose:

  Whether to print progress.

## Value

The path written, invisibly.

## Details

Writes a folder holding the viewer and the data it draws, and zips it.
The person receiving it double-clicks \`ClonoView.html\`; there is
nothing to install, and it does not use the internet. This is the form
to send to a collaborator or a reviewer, or to attach as supplementary
material.

**The export contains the data, not a picture of it.** Per-cell
coordinates, cell types, subject identifiers and every clone identifier
– which for an immune-repertoire run means the receptor sequences – are
all in the folder, because the viewer computes from them in the browser.
Treat sending one as publishing that data.

## See also

\[clonoview_bundle()\], \[clonoview()\]
