# ==========================================================================
# 'plgem' package initialization
# ==========================================================================

.First.lib <- function(lib, pkg) {
  message(paste("\nWelcome to", pkg, "version", packageDescription(pkg)$Version,
    "\n"))
  addVigs2WinMenu(pkg)
}
