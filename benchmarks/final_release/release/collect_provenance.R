# Collect local release-candidate environment provenance.

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1L) args[[1L]] else "environment_provenance.txt"

connection <- file(output, open = "wt", encoding = "UTF-8")
sink(connection)
sink(connection, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(connection)
}, add = TRUE)

cat("R_VERSION=", R.version.string, "\n", sep = "")
cat("R_PLATFORM=", R.version$platform, "\n", sep = "")
cat("R_ARCH=", R.version$arch, "\n", sep = "")
cat("R_HOME=", R.home(), "\n", sep = "")
cat("PACKAGE_VERSION=", as.character(desc::desc_get_version(".")), "\n", sep = "")
cat("LOCALE=", paste(Sys.getlocale(), collapse = ";"), "\n", sep = "")
cat("NCORES_LOGICAL=", parallel::detectCores(logical = TRUE), "\n", sep = "")
cat("NCORES_PHYSICAL=", parallel::detectCores(logical = FALSE), "\n", sep = "")
cat("\nEXT_SOFT_VERSION\n")
print(extSoftVersion())
cat("\nCAPABILITIES\n")
print(capabilities())
cat("\nSESSION_INFO\n")
print(sessionInfo())
