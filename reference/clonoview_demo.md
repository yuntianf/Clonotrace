# clonoview_demo

Bundle the shipped example run and open it in ClonoView.

## Usage

``` r
clonoview_demo(
  open = interactive(),
  rebuild = FALSE,
  python = NULL,
  verbose = TRUE
)
```

## Arguments

- open:

  Whether to open a browser. See \[clonoview()\].

- rebuild:

  Rebuild the bundle even if it is already prepared.

- python:

  Python interpreter to use. Defaults to \[clonoview_python()\].

- verbose:

  Whether to print progress.

## Value

The URL, invisibly.

## Details

One call, from a fresh install, to a viewer in the browser. It is the
fastest way to find out whether ClonoView works on a given machine, and
it exercises the whole path rather than a stored result: the example
ships as ordinary R objects and is put through \[clonoview_bundle()\]
here exactly as your own run would be. Bundling takes a few seconds; it
is skipped on later calls unless \`rebuild = TRUE\`.

## See also

\[clonoview_example()\], \[clonoview_bundle()\], \[clonoview_stop()\]

## Examples

``` r
if (FALSE) { # \dontrun{
clonoview_demo()
} # }
```
