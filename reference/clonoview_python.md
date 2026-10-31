# clonoview_python

The Python interpreter used to run ClonoView.

## Usage

``` r
clonoview_python(check = TRUE)
```

## Arguments

- check:

  Whether to verify the interpreter can import what the viewer needs.

## Value

Path to a Python interpreter.

## Details

In order: \`options(clonoview.python)\`, then the interpreter reticulate
is configured with – so \`use_python()\` from the co-embedding step
covers this too – then whatever \`python3\` is on the PATH.

ClonoView needs numpy, scipy and matplotlib. numpy and scipy are already
in the package's \`requirements.txt\`; matplotlib is the one addition
the viewer makes, and it is only used to rasterise the cell panel's
backdrop.
