# The I/O sandbox must (1) translate the author's hard-coded paths,
# (2) read inputs from CHIK_DATA_DIR, and (3) never write outside
# CHIK_OUTPUT_DIR. Exercised end-to-end on a throwaway project in tempdir(),
# in a separate R process, exactly as tools/run_reproduction.R does.

test_that("sandbox redirects reads and writes and leaves the project untouched", {
  proj <- file.path(tempdir(), paste0("chik_sandbox_test_", Sys.getpid()))
  unlink(proj, recursive = TRUE)
  dir.create(file.path(proj, "02_Scripts"), recursive = TRUE)
  dir.create(file.path(proj, "06_Results"))
  file.create(file.path(proj, "CHIK_benefit_risk_modelling.Rproj"))
  data_dir <- file.path(proj, "elsewhere_data"); dir.create(data_dir)
  out_dir <- file.path(proj, "reproduced_outputs")

  setting_key <- c(A = "Low"); save(setting_key, file = file.path(data_dir, "setting_key.RData"))
  writeLines(c(
    'setwd("c:/Users/user/OneDrive/CHIK_benefit_risk")',     # legacy-path-fixture: author's original root
    'load("01_Data/setting_key.RData")',                       # relative data path
    'stopifnot(identical(setting_key, c(A = "Low")))',
    'write.csv(data.frame(x = 1), "06_Results/out.csv", row.names = FALSE)',
    'save(setting_key, file = "01_Data/cache.RData")',
    'load("01_Data/cache.RData")'                              # must read the reproduced copy
  ), file.path(proj, "02_Scripts", "toy.R"))

  old <- Sys.getenv(c("CHIK_PROJECT_ROOT", "CHIK_DATA_DIR", "CHIK_OUTPUT_DIR", "CHIK_ORV_ROOT"), unset = NA)
  on.exit({ for (n in names(old)) if (is.na(old[[n]])) Sys.unsetenv(n) else do.call(Sys.setenv, as.list(old[n])) },
          add = TRUE)
  Sys.setenv(CHIK_PROJECT_ROOT = proj, CHIK_DATA_DIR = data_dir, CHIK_OUTPUT_DIR = out_dir,
             CHIK_ORV_ROOT = proj)
  log <- tempfile()
  owd <- setwd(proj); on.exit(setwd(owd), add = TRUE)
  status <- system2(rscript_bin(), c(shQuote(file.path(ROOT, "tools/repro/run_script.R"), type = if (.Platform$OS.type == "windows") "cmd" else "sh"),
                                     "02_Scripts/toy.R"), stdout = log, stderr = log)
  expect_equal(status, 0, info = paste(readLines(log), collapse = "\n"))
  expect_true(file.exists(file.path(out_dir, "06_Results", "out.csv")))
  expect_true(file.exists(file.path(out_dir, "01_Data", "cache.RData")))
  expect_false(file.exists(file.path(proj, "06_Results", "out.csv")))
  expect_false(dir.exists(file.path(proj, "01_Data")))
  io <- read.csv(file.path(out_dir, "logs", "io_log.csv"))
  expect_true(any(io$op == "read" & grepl("elsewhere_data/setting_key.RData$", io$resolved)))
  expect_true(any(io$op == "read" & grepl("reproduced_outputs/01_Data/cache.RData$", io$resolved)))
  unlink(proj, recursive = TRUE)
})
