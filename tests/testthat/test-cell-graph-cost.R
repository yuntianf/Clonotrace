test_that("SNN similarities are returned as shortest-path costs", {
  graph <- igraph::make_ring(4)
  graph <- igraph::set_edge_attr(graph, "weight", value = c(1, 2, 4, 8))

  corrected <- Clonotrace:::.clonotrace_snn_similarity_to_cost(graph)

  expect_equal(igraph::edge_attr(corrected, "weight"), exp(-c(1, 2, 4, 8)))
  expect_equal(igraph::edge_attr(graph, "weight"), c(1, 2, 4, 8))
  expect_true(all(diff(igraph::edge_attr(corrected, "weight")) < 0))
})

test_that("invalid SNN similarities are rejected", {
  graph <- igraph::make_ring(3)
  graph <- igraph::set_edge_attr(graph, "weight", value = c(1, -1, 2))

  expect_error(
    Clonotrace:::.clonotrace_snn_similarity_to_cost(graph),
    "finite, non-negative"
  )
})
