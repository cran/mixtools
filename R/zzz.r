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
    cat('mixtools package version 0.3.0,  released February 2008\n')
    cat('Type help(package="mixtools") to get started.\n')
}
