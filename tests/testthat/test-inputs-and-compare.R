# Input validation and output comparison, using small synthetic files
# (clearly test-only; never used to reproduce results).

test_that("validate_rdata reports missing objects and columns", {
  f <- tempfile(fileext = ".RData")
  psa_df <- data.frame(draw = 1:3, state = "x")
  save(psa_df, file = f)
  spec <- list(psa_df = c("draw", "state", "daly_averted"), other = character(0))
  msg <- validate_rdata(f, spec)
  expect_match(msg, "lacks column\\(s\\): daly_averted")
  expect_match(msg, "object 'other' missing")
  expect_identical(validate_rdata(f, list(psa_df = c("draw", "state"))), "")
  expect_identical(validate_rdata(tempfile()), "file not found")
})

test_that("validate_rdata flags empty data frames and unreadable files", {
  f <- tempfile(fileext = ".RData")
  ar_summary <- data.frame(region = character(0), ar_med = numeric(0))
  save(ar_summary, file = f)
  expect_match(validate_rdata(f, list(ar_summary = c("region", "ar_med"))), "0 rows")
  bad <- tempfile(fileext = ".RData"); writeLines("not an RData file", bad)
  expect_match(suppressWarnings(validate_rdata(bad, list())), "cannot load")
})

test_that("check_inputs marks missing inputs as MISSING, never OK", {
  cfg <- list(root = tempdir(), data_dir = file.path(tempdir(), "no_such_dir"),
              orv_root = file.path(tempdir(), "no_orv"), out_dir = tempdir())
  res <- check_inputs(PUBLICATION_MANIFEST[1:2], cfg)
  expect_true(all(res$status == "MISSING"))
  expect_true(all(grepl("Needed for", res$detail)))
})

test_that("parse_numbers extracts estimate and uncertainty interval", {
  expect_equal(parse_numbers("5.19 (1.42–16.65)")[[1]], c(5.19, 1.42, 16.65))
  expect_equal(parse_numbers("0.00 (0.00-0.01)")[[1]], c(0, 0, 0.01))
})

test_that("compare_table: identical tables PASS, a changed cell FAILs with its difference", {
  skip_if_not_installed("writexl"); skip_if_not_installed("readxl")
  ref <- tempfile(fileext = ".xlsx"); rep <- tempfile(fileext = ".xlsx")
  tb <- data.frame(Age = c("18-64", "65+"), BRR = c("9.53 (2.51–34.22)", "0.08 (0.02–0.21)"))
  writexl::write_xlsx(tb, ref); writexl::write_xlsx(tb, rep)
  expect_identical(compare_table("t", rep, ref)$status, "PASS")
  tb2 <- tb; tb2$BRR[2] <- "0.09 (0.02–0.21)"
  writexl::write_xlsx(tb2, rep)
  r <- compare_table("t", rep, ref)
  expect_identical(r$status, "FAIL")
  expect_equal(r$n_mismatch, 1)
  expect_equal(r$max_abs_diff, 0.01, tolerance = 1e-12)
})

test_that("compare_grid: numeric tolerance, row order ignored", {
  skip_if_not_installed("writexl"); skip_if_not_installed("readxl")
  g <- data.frame(Setting = c("Low", "High"), `Age group` = "18-64", `Duration (days)` = 7,
                  `Entry day (day of year)` = c(1, 8), x = c(0.1, 0.2), check.names = FALSE)
  ref <- tempfile(fileext = ".xlsx"); rep <- tempfile(fileext = ".xlsx")
  writexl::write_xlsx(g, ref); writexl::write_xlsx(g[2:1, ], rep)
  expect_identical(compare_grid("g", rep, ref)$status, "PASS")
  g$x[1] <- 0.1 * (1 + 1e-6)
  writexl::write_xlsx(g, rep)
  expect_identical(compare_grid("g", rep, ref)$status, "FAIL")
  expect_identical(compare_grid("g", rep, ref, rel_tol = 1e-5)$status, "PASS")
})

test_that("compare_png: identical images PASS, altered pixels FAIL", {
  skip_if_not_installed("png")
  img <- array(0.5, c(10, 10, 3))
  a <- tempfile(fileext = ".png"); b <- tempfile(fileext = ".png")
  png::writePNG(img, a); png::writePNG(img, b)
  expect_identical(compare_png("p", a, b)$status, "PASS")
  img[1, 1, 1] <- 1; png::writePNG(img, b)
  r <- compare_png("p", a, b)
  expect_identical(r$status, "FAIL")
  expect_equal(r$n_mismatch, 1)
})

test_that("compare_one: missing reference is BLOCKED, missing output is FAIL", {
  expect_identical(compare_one("x", tempfile(), tempfile(), "png")$status, "BLOCKED")
  ref <- tempfile(); writeLines("x", ref)
  expect_identical(compare_one("x", tempfile(), ref, "png")$status, "FAIL")
})
