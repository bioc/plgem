# ==========================================================================
# 'plgem' package initialization
# ==========================================================================

.onLoad <- function(lib, pkg) {
  message(paste("\nWelcome to", pkg, "version", packageDescription(pkg)$Version,
    "\n"))
  addVigs2WinMenu(pkg)
}
