# Static checks: every publication script parses, the manifest is
# consistent, and the reproduction tooling itself contains no machine-specific
# paths.

test_that("all R scripts in 02_Scripts/ and tools/ parse", {
  files <- c(list.files(file.path(ROOT, "02_Scripts"), "\\.R$", full.names = TRUE),
             list.files(file.path(ROOT, "tools"), "\\.R$", full.names = TRUE, recursive = TRUE))
  for (f in files) {
    res <- tryCatch({ parse(f, keep.source = FALSE); TRUE }, error = function(e) conditionMessage(e))
    expect_true(isTRUE(res), info = paste(basename(f), if (!isTRUE(res)) res))
  }
})

test_that("manifest: unique ids, scripts exist, well-formed checks", {
  ids <- vapply(PUBLICATION_MANIFEST, `[[`, "", "id")
  expect_false(anyDuplicated(ids) > 0)
  for (it in PUBLICATION_MANIFEST) {
    expect_true(file.exists(file.path(ROOT, it$script)), info = it$script)
    expect_true(length(it$checks) > 0, info = it$id)
    for (ck in it$checks) {
      expect_true(ck$type %in% c("table", "grid", "png"), info = it$id)
      expect_false(grepl("^([A-Za-z]:)?/", ck$reproduced), info = "reproduced paths must be relative")
    }
    for (inp in it$inputs) expect_true(inp$area %in% c("data", "orv"), info = inp$file)
  }
})

test_that("manifest: reference files are present in the working copy", {
  refs <- unlist(lapply(PUBLICATION_MANIFEST, function(it) vapply(it$checks, `[[`, "", "reference")))
  missing <- refs[!file.exists(file.path(ROOT, refs))]
  # On CI the checkout only contains committed files; references that the
  # author has not yet committed are reported rather than failing the build.
  if (nzchar(Sys.getenv("CI")) && length(missing))
    skip(paste("reference file(s) not committed:", paste(basename(missing), collapse = ", ")))
  expect_length(missing, 0)
})

test_that("reproduction tooling and tests contain no machine-specific paths", {
  files <- c(list.files(file.path(ROOT, "tools"), "\\.R$", full.names = TRUE, recursive = TRUE),
             list.files(file.path(ROOT, "tests"), "\\.R$", full.names = TRUE, recursive = TRUE))
  for (f in files) {
    txt <- readLines(f, warn = FALSE)
    # LEGACY_ROOTS in config.R is the single, documented place that names the
    # author's original locations (used to translate them, never to read them).
    # Lines marked "legacy-path-fixture" feed a legacy path INTO the sandbox
    # test on purpose.
    txt <- txt[!grepl("LEGACY_ROOTS|project = \"c:/|orv     = \"c:/|legacy-path-fixture|path_pattern <-", txt)]
    path_pattern <- paste0("[\"'][A-Za-z]:[/\\\\]|/Users/[a-z]|", "One", "Drive")
    hits <- grep(path_pattern, txt, value = TRUE)
    hits <- hits[!grepl("^\\s*#", hits)]
    expect_length(hits, 0)
    if (length(hits)) message(basename(f), ": ", paste(hits, collapse = " | "))
  }
})

test_that("upstream scripts set the documented random seeds", {
  seeds <- function(f) regmatches(readLines(file.path(ROOT, f), warn = FALSE),
                                  regexpr("set\\.seed\\([^)]*\\)", readLines(file.path(ROOT, f), warn = FALSE)))
  expect_true("set.seed(123)" %in% seeds("02_Scripts/01_setup.R"))
  expect_true("set.seed(20260420)" %in% seeds("02_Scripts/02_setup_age_props.R"))
  expect_true("set.seed(123)" %in% seeds("02_Scripts/06_brazil_travel_final_finite.R"))
  expect_true("set.seed(1)" %in% seeds("02_Scripts/03_brazil_all_draws_ori_v3.R"))
})
