# cell_clone_coembed

Construct Cell-Cell Distance Matrix via Clone-Aware Coembedding

## Usage

``` r
cell_clone_coembed(
  cell_embedding,
  clone_embedding,
  cell_clone_prob,
  cell_k = 30,
  clone_k = 15
)
```

## Arguments

- cell_embedding:

  A numeric matrix of cell embeddings (rows = cells).

- clone_embedding:

  A numeric matrix of clone embeddings (rows = clones).

- cell_clone_prob:

  A numeric or sparse matrix of cell-by-clone membership probabilities.
  Row names must identify cells and column names must identify clones.

- cell_k:

  Integer. Number of neighbors to compute for cell-level kNN (default:
  30).

- clone_k:

  Integer. Number of neighbors to compute for clone-level kNN (default:
  15).

## Value

A symmetric sparse distance matrix (class \`dgCMatrix\`) encoding
clone-informed cell-cell distances.

## Details

Integrates cell and clone embeddings by using clone-level proximity to
reweight cell-cell distances. Generates a sparse distance matrix
suitable for downstream embedding or trajectory inference.
