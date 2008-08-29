######################################################################
#
# zzz.r
#
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/mixtools package
#
# .First.lib is run when the package is loaded with library(mixtools)
#
######################################################################

.First.lib <- function(lib, pkg){
    library.dynam("mixtools", pkg, lib)
    DESCpath <- file.path(system.file(package="mixtools"), "DESCRIPTION")
    info <- read.dcf(DESCpath)
#   info <- read.dcf(file.path(lib, pkg, "DESCRIPTION")) does the same thing but
#   the other way works even outside the .First.lib function
    cat(pkg,'package version', info[,"Version"], '  Released', info[,"Date"], '\n')
    cat('Type help(package="mixtools") to get started.\n')
}

.Last.lib <- function(libpath){
  library.dynam.unload("mixtools",libpath)
}

