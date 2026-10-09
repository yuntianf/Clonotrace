test_that("long2square converts one-direction pair data", {
  long <- data.frame(
    group1 = c("clone1", "clone1", "clone2"),
    group2 = c("clone2", "clone3", "clone3"),
    dis = c(1, 2, 3)
  )

  result <- long2square(long, "group1", "group2", "dis", symmetric = TRUE)

  expect_type(result, "double")
  expect_identical(rownames(result), colnames(result))
  expect_equal(result["clone1", "clone2"], 1)
  expect_equal(result["clone2", "clone1"], 1)
  expect_equal(result["clone2", "clone3"], 3)
})

test_that("long2square preserves explicitly requested isolated nodes", {
  long <- data.frame(from = "a", to = "b", value = 1)
  result <- long2square(long, "from", "to", "value", nodes = c("a", "b", "c"))
  expect_identical(dim(result), c(3L, 3L))
  expect_true(all(is.na(result["c", ])))
  empty <- long2square(data.frame(), "from", "to", "value", nodes = c("b", "a"))
  expect_identical(dimnames(empty), list(c("a", "b"), c("a", "b")))
})

test_that("long2square rejects malformed and ambiguous pair tables", {
  duplicate <- data.frame(from = c("a", "a"), to = c("b", "b"), value = c(1, 1))
  conflict <- data.frame(from = c("a", "b"), to = c("b", "a"), value = c(1, 2))

  expect_error(long2square(duplicate, "from", "to", "value"), "Duplicate")
  expect_error(long2square(conflict, "from", "to", "value", symmetric = TRUE), "Conflicting")
  expect_error(long2square(data.frame(), "from", "to", "value"), "Missing required")
})
