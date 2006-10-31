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
#    if(R.version$major=="1"){
#     ehelp <- help(package="mixtools")$info[[2]][[2]]
#     cat(paste("'",ehelp[4],"'\n",
#               "Version ",ehelp[2],
#               " created on ",ehelp[3],".\n", sep=""))
#    }else{
#     ehelp <- help(package="mixtools")$info[[2]]
#     cat(paste(substring(ehelp[4],first=16),"\n",
#               "Version ",substring(ehelp[2],first=16),
#               " created on ",
#                substring(ehelp[3],first=16),".\n", sep=""))
#    }
    cat('mixtools package, created April 2005\n')
    cat('Type help(package="mixtools") to get started.\n')
}
