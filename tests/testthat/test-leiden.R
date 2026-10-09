test_that("leiden_dis validates distance input and preserves names", {
  set.seed(1)
  coordinates <- matrix(rnorm(40), nrow = 10)
  distances <- as.matrix(stats::dist(coordinates))
  dimnames(distances) <- list(paste0("clone", 1:10), paste0("clone", 1:10))

  result <- leiden_dis(distances, k = 3, if_umap = FALSE)

  expect_s3_class(result, "factor")
  expect_identical(names(result), rownames(distances))
  expect_error(leiden_dis(distances[1:9, ], k = 3, if_umap = FALSE), "square")
  asymmetric <- distances
  asymmetric[1, 2] <- asymmetric[1, 2] + 1
  expect_error(leiden_dis(asymmetric, k = 3, if_umap = FALSE), "symmetric")
})
