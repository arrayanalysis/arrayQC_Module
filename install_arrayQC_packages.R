## Packages used for the following functions:
# - affy (?)
# - limma (read.maimages, normalizeWithinArrays, normalizeBetweenArrays
# - RColorBrewer (colorpallette functions)
# - gplots (heatmap.2)
# - geneplotter (MvA plots (Density))
# - shape

.checkPackages <- function(x) {
 #v0.01
 # This function will check if a required package is not installed on the system. If this is the case,
 # it will try to download the package from bioconductor.
 installedPackages <- rownames(installed.packages())
 a <- x %in% installedPackages
 if(sum(a) != length(x)) {
   cat(" -- The following package(s) are missing and will be installed: ")
   cat(x[!a])
   cat("\n\n")
   source("http://www.bioconductor.org/biocLite.R")
   biocLite(x[!a])
   cat(" -- All required arrayQC packages were successfully installed.\n")
 }
}

#.checkPackages(c("affy","limma","RColorBrewer","gplots","geneplotter", "shape", "plotrix"))
.checkPackages(c("affy","limma","RColorBrewer","gplots","geneplotter", "shape"))

