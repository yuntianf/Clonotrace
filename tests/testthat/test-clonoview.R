test_that("ClonoView finds the current process on supported platforms", {
  expect_true(Sys.getpid() %in% Clonotrace:::cv_running_pids())
})

test_that("ClonoView preserves the requested dataset and ignores a stale pidfile", {
  skip_if(utils::packageVersion("testthat") < "3.2.0", "Requires testthat mocking support")
  state <- new.env(parent = emptyenv())
  pid_reads <- 0L
  testthat::local_mocked_bindings(
    cv_state = state,
    clonoview_datasets = function() c("first", "requested"),
    cv_free_port = function() 8770L,
    cv_python_run = function(...) 0L,
    cv_running_pids = function() 222L,
    cv_pid = function() {
      pid_reads <<- pid_reads + 1L
      if (pid_reads <= 2L) 111L else 222L
    },
    .package = "Clonotrace"
  )
  expect_message(url <- clonoview("requested", open = FALSE, python = "unused"),
                 "ClonoView is serving")
  expect_identical(url, "http://localhost:8770/?dataset=requested")
  expect_identical(state$pid, 222L)
  expect_gte(pid_reads, 3L)
})

test_that("ClonoView refuses a second server from the same session", {
  skip_if(utils::packageVersion("testthat") < "3.2.0", "Requires testthat mocking support")
  state <- new.env(parent = emptyenv())
  state$pid <- 222L
  state$port <- 8770L
  testthat::local_mocked_bindings(
    cv_state = state,
    clonoview_datasets = function() "example",
    cv_running_pids = function() 222L,
    cv_python_run = function(...) stop("A second process must not start"),
    .package = "Clonotrace"
  )
  expect_error(clonoview("example", open = FALSE), "already running")
})

test_that("viewer subprocesses use user caches and restore environment settings", {
  cache <- tempfile("clonoview-environment-")
  old_options <- options(clonoview.cache = cache)
  keys <- c("PYTHONDONTWRITEBYTECODE", "MPLCONFIGDIR", "PYTHONPATH", "CLONOVIEW_DATA")
  old_env <- Sys.getenv(keys, unset = NA_character_)
  on.exit({
    options(old_options)
    for (key in keys) {
      if (is.na(old_env[[key]])) Sys.unsetenv(key) else {
        do.call(Sys.setenv, stats::setNames(list(old_env[[key]]), key))
      }
    }
    unlink(cache, recursive = TRUE)
  }, add = TRUE)
  Sys.setenv(PYTHONDONTWRITEBYTECODE = "0")
  Sys.unsetenv("MPLCONFIGDIR")
  # A short R subprocess inspects the environment inherited by the Python runner;
  # this test does not require Python or optional viewer dependencies.
  expression <- 'cat(Sys.getenv("PYTHONDONTWRITEBYTECODE"), Sys.getenv("MPLCONFIGDIR"), Sys.getenv("CLONOVIEW_DATA"), sep = "\\n")'
  output <- Clonotrace:::cv_python_run(
    c("-e", shQuote(expression)), stdout = TRUE,
    py = file.path(R.home("bin"), "Rscript")
  )
  expect_identical(output[[1L]], "1")
  expect_identical(output[[2L]], normalizePath(file.path(cache, "matplotlib"), winslash = "/"))
  expect_identical(output[[3L]], normalizePath(file.path(cache, "data"), winslash = "/"))
  expect_identical(Sys.getenv("PYTHONDONTWRITEBYTECODE"), "0")
  expect_identical(Sys.getenv("MPLCONFIGDIR", unset = NA_character_), NA_character_)
})

test_that("profile columns retain their identities rather than their positions", {
  probability <- matrix(c(0.2, 0.8, 0.7, 0.3), 2, 2,
                        dimnames = list(c("cell1", "cell2"), c("2", "1")))
  table <- Clonotrace:::cv_profile_table(probability, rownames(probability))
  expect_identical(names(table), c("cell", "X2", "X1"))
  expect_equal(table$X2, probability[, "2"], ignore_attr = TRUE)
  expect_error(
    Clonotrace:::cv_profile_table(
      matrix(1, 2, 2, dimnames = list(c("a", "b"), c("response", "other"))),
      c("a", "b")
    ),
    "profile numbers"
  )
})

test_that("the shipped example bundles, serves, exports, and stops", {
  skip_if(Sys.getenv("CLONOVIEW_RUN_INTEGRATION") != "true",
          "Opt-in Python integration test; exercised by the ClonoView workflow")
  skip_if_not_installed("jsonlite")
  cache <- tempfile("clonoview-integration-")
  python <- Sys.getenv("CLONOVIEW_TEST_PYTHON", unset = Sys.which("python"))
  previous <- options(clonoview.cache = cache, clonoview.python = python)
  on.exit({
    clonoview_stop()
    options(previous)
    unlink(cache, recursive = TRUE)
  }, add = TRUE)
  home <- Clonotrace:::clonoview_home()
  before <- list.files(home, all.files = TRUE, recursive = TRUE)
  url <- clonoview_demo(open = FALSE, rebuild = TRUE, verbose = FALSE)
  expect_match(url, "[?]dataset=mathew_nsclc_example$")
  catalogue <- jsonlite::fromJSON(paste0(sub("[?].*$", "", url), "api/datasets"))
  expect_true("mathew_nsclc_example" %in% catalogue$datasets$name)
  meta <- jsonlite::fromJSON(paste0(sub("[?].*$", "", url), "api/dataset/mathew_nsclc_example"))
  expect_equal(meta$meta$counts$cells, 150150)
  expect_equal(meta$meta$counts$clone_units, 916)
  zip <- clonoview_export(file.path(cache, "offline.zip"), verbose = FALSE)
  expect_true(file.exists(zip))
  expect_true(file.exists(file.path(cache, "offline", "ClonoView.html")))
  expect_setequal(list.files(home, all.files = TRUE, recursive = TRUE), before)
  expect_true(clonoview_stop())
})
