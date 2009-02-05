# ==========================================================================
# 'plgem' package initialization
# ==========================================================================

.onLoad <- function(lib, pkg) {
  pkgVersion <- packageDescription(pkg)$Version
  msg <- paste("\nWelcome to", pkg, "version", pkgVersion, "\n")
  message(msg)
  Biobase::addVigs2WinMenu(pkg)
}
