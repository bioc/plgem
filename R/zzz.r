# ==========================================================================
# 'plgem' package initialization
# ==========================================================================

.onAttach <- function(lib, pkg) {
  pkgVersion <- packageDescription(pkg)$Version
  msg <- paste("\nWelcome to", pkg, "version", pkgVersion, "\n")
  packageStartupMessage(msg)
  Biobase::addVigs2WinMenu(pkg)
}
