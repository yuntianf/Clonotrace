# clonoview_example

The worked example that ships with the package: a real Clonotrace run.

## Usage

``` r
clonoview_example()
```

## Value

A named list of four objects.

## Details

Returns a list of the four objects \[clonoview_bundle()\] takes, so it
doubles as the clearest statement of what those arguments have to look
like:

- \`cell_embedding\`:

  150,150 x 2, cells in row names

- \`cell_meta\`:

  150,150 rows with \`clone\`, \`cell_type\`, \`patient\`, \`visit\` and
  \`Response\`; 94,390 of the cells carry a clone and the rest are
  context

- \`clone_cluster\`:

  916 clone units with coordinates and a \`cluster\` column

- \`cell_profile_prob\`:

  150,150 x 10

The data is the NSCLC anti-PD-1 plus JAK1-inhibitor trial of Mathew et
al., eight patients sampled over five treatment cycles, 16 cell types.
Clone identifiers are spelled \`pt12cy1
TRB:CALRQGAGELFF;TRA:CAVSAERGNMLTF\` – patient and cycle run together,
then a space, then the receptor – which is what the \`clone_id_regex\`
in the example below reads.

## See also

\[clonoview_demo()\], \[clonoview_bundle()\]
