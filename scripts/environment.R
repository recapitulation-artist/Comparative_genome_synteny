#!/usr/bin/env Rscript
# Restore or verify the complete R dependency lock.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[1]
project_root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg))))
package_lock <- read.delim(file.path(project_root, "environment", "r_packages.lock.tsv"),
                          stringsAsFactors = FALSE, na.strings = "", check.names = FALSE)
project_library <- file.path(project_root, ".r_library")
Sys.setenv(TZ = "UTC")
if (dir.exists(project_library)) .libPaths(c(project_library, .libPaths()))

installed_version <- function(package) {
  path <- find.package(package, quiet = TRUE)
  if (!length(path)) return(NA_character_)
  read.dcf(file.path(path, "DESCRIPTION"), fields = "Version")[1, 1]
}

check_environment <- function() {
  if (getRversion() != "4.6.0") stop("R 4.6.0 is required; use the Conda lock.")
  actual <- vapply(package_lock$package, installed_version, character(1))
  bad <- is.na(actual) | actual != package_lock$version
  if (any(bad)) stop("Missing or mismatched packages: ", paste(package_lock$package[bad], collapse = ", "),
                      ". Run make setup.")
  if (!capabilities("cairo")) stop("R must support Cairo graphics.")
  invisible(TRUE)
}

restore_packages <- function() {
  if (getRversion() != "4.6.0") stop("Create the locked Conda environment first.")
  dir.create(project_library, showWarnings = FALSE)
  .libPaths(c(project_library, .Library), include.site = FALSE)
  options(timeout = 600)
  cran <- "https://cloud.r-project.org/src/contrib"
  for (i in seq_len(nrow(package_lock))) {
    p <- package_lock$package[i]
    version <- package_lock$version[i]
    if (identical(installed_version(p), version)) next
    archive <- sprintf("%s_%s.tar.gz", p, version)
    urls <- c(paste(cran, archive, sep = "/"), paste(cran, "Archive", p, archive, sep = "/"))
    dest <- tempfile(fileext = ".tar.gz")
    downloaded <- FALSE
    for (url in urls) {
      downloaded <- tryCatch(suppressWarnings(download.file(url, dest, mode = "wb", quiet = TRUE)) == 0,
                             error = function(e) FALSE)
      if (downloaded) break
    }
    if (!downloaded) stop("Cannot download locked package: ", archive)
    install.packages(dest, repos = NULL, type = "source", lib = project_library, quiet = TRUE)
    unlink(dest)
    if (!identical(installed_version(p), version)) stop("Installation failed: ", p, " ", version)
  }
  check_environment()
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (identical(args, "--install")) restore_packages() else {
    if (length(args)) stop("Usage: Rscript scripts/environment.R [--install]")
    check_environment()
  }
  cat(sprintf("Environment verified: R 4.6.0 and %d locked packages.\n", nrow(package_lock)))
}
