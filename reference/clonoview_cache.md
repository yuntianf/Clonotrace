# clonoview_cache

Where ClonoView writes prepared bundles, staged inputs and its pidfile.

## Usage

``` r
clonoview_cache(sub = NULL)
```

## Arguments

- sub:

  Optional sub-directory, created if it does not exist.

## Value

A path.

## Details

Not inside the installed package: that directory is read-only on a
managed R installation and shared between users on a cluster. The
default is the standard per-user cache for this package, overridable
with \`options(clonoview.cache = ...)\`.
