# =============================================================================
# tools/repro/config.R
#
# Path configuration for the reproduction tooling. Nothing here changes any
# scientific script; it only decides where inputs are read from and where
# reproduced outputs are written.
#
# Environment variables (all optional):
#   CHIK_PROJECT_ROOT   this repository (default: detected from the working dir)
#   CHIK_DATA_DIR       directory holding the 01_Data/ inputs
#                       (default: <root>/01_Data)
#   CHIK_ORV_ROOT       root of the companion CHIK_ORV_impact project, which
#                       holds the transmission-model fit, observed case data
#                       and outbreak-response week-sweep simulations
#                       (default: <root>/../CHIK_vaccine_impact/CHIK_ORV_impact)
#   CHIK_OUTPUT_DIR     where reproduced outputs are written
#                       (default: <root>/reproduced_outputs)
#
# The scientific scripts hard-code the author's original locations; these are
# listed in LEGACY_ROOTS and translated by tools/repro/sandbox.R.
# =============================================================================

chik_find_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (dir.exists(file.path(d, "02_Scripts")) && file.exists(file.path(d, "CHIK_benefit_risk_modelling.Rproj")))
      return(d)
    parent <- dirname(d)
    if (identical(parent, d)) return(NA_character_)
    d <- parent
  }
}

chik_norm <- function(p) {
  if (is.null(p) || !nzchar(p)) return(p)
  gsub("\\\\", "/", normalizePath(p, winslash = "/", mustWork = FALSE))
}

chik_config <- function() {
  root <- Sys.getenv("CHIK_PROJECT_ROOT", "")
  if (!nzchar(root)) root <- chik_find_root()
  if (is.na(root) || !dir.exists(root))
    stop("Cannot locate the project root. Run from inside the repository or set CHIK_PROJECT_ROOT.", call. = FALSE)
  root <- chik_norm(root)
  env_or <- function(var, default) {
    v <- Sys.getenv(var, "")
    chik_norm(if (nzchar(v)) v else default)
  }
  list(
    root     = root,
    data_dir = env_or("CHIK_DATA_DIR", file.path(root, "01_Data")),
    orv_root = env_or("CHIK_ORV_ROOT", file.path(dirname(root), "CHIK_vaccine_impact", "CHIK_ORV_impact")),
    out_dir  = env_or("CHIK_OUTPUT_DIR", file.path(root, "reproduced_outputs"))
  )
}

# Hard-coded locations used in the scientific scripts (author's machine).
# Matching is case-insensitive and slash-agnostic.
LEGACY_ROOTS <- list(
  project = "c:/Users/user/OneDrive/CHIK_benefit_risk",
  orv     = "c:/Users/user/OneDrive/CHIK_vaccine_impact/CHIK_ORV_impact"
)
