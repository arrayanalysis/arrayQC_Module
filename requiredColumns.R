## x.R
#
#
# Every function will first check and see which columns are required. If these column names / variable names
# are present, the basic functionality will be executed. If the optional columns are present, then the full
# function will be executed.
#
# For some functions there is no need to look for required column names, since they're rather specific themselves.
#
# Note that this function takes over one global variable name, and that is scriptpath!

.checkColumns <- function(functionName=NULL, requiredFile="requiredColumns.arrayQC", silent=FALSE, local.path = datapath, scriptpath = arrayQC.scriptpath) {
  error <- NULL
  if(is.null(functionName)) { error <- c(error, "- functionName\n") }
  if(!is.null(error)) { cat("[[WARNING]] .checkColumns stopped due to undefined variable(s):\n"); stop(cat(paste(error, sep=""))) }

  temp <- paste(local.path,"/",requiredFile, sep="")

#  if(!file.exists( temp )) { 
#    cat(paste("* Downloading", requiredFile, "...")) 
#    download.file(paste(arrayQC.scriptpath, requiredFile, sep=""), temp, quiet=TRUE)
#    cat(" done.\n")
#  }
  a <- as.matrix(read.table(file=paste(scriptpath, requiredFile, sep=""), sep="\t", row.names=1, colClasses="character", skip=9, fill=TRUE, stringsAsFactors=FALSE))
  if(silent == 0) { cat(paste("* Checking required / optional columns of ", functionName, "...\n", sep="")) }
  if(length(grep(functionName, rownames(a)))!=1) {
    stop(paste("[[ ERROR ]]", functionName, "can not be found in", requiredFile, "!\n"))
  } else {
    x <- list()
    b <- a[functionName,]
    y <- intersect( grep("[[", b, fixed=TRUE), grep("]]", b, fixed=TRUE) )
    if(length(y) == 0) {
      x[["required"]] <- as.vector(b[b[]!=""])
      x[["optional"]] <- NULL
    } else {
  ## Look for overlap between all values (without the optional value) and the values that are filled in
      x[["required"]] <- as.vector( intersect( b[setdiff(1:length(b), y)] ,b[b[]!=""] ) )
  ## Remove the "[[" and "]]" from the optional names
      b1 <- as.vector(gsub("[[","", b[y], fixed=TRUE))
      b2 <- gsub("]]","", b1, fixed=TRUE)
      x[["optional"]] <- as.vector( b2 )
    }
    return(x)
  }
}
