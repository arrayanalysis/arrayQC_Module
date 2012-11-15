
#CreateSummaryPlots <- function(RG) {
#
#  if(is.null(RG$summary)) stop("For summary plots to be created, a summary field must be present in the data object, call CalculateSummary first")
#  RGsummary <- RG$summary[1:(nrow(RG$summary)-4),]
#  totalArrays <- nrow(RGsummary)
#  npages <- ceiling(totalArrays/16)
#  x <- totalArrays/npages
#  ArraysPerPage <- ceiling(x)
#  hasRealBG <- !is.null(RG$other$rBGMeanSignal)|!is.null(RG$other$gBGMeanSignal)
#  hasManualFlags <- !is.null(RG$other$manualFlags)
#  
#  main.titles <- c("Foreground","Background",rep("Real Background",hasRealBG),"Number of Good Probes","Percentage of Good Probes",RG$QC.vector,rep("Manually Flagged",hasManualFlags))
#  y.axis.labels <- c(rep("Mean Signal",2+hasRealBG),"Number of reporters","% of reporters",rep("Number of reporters",length(RG$QC.vector)+hasManualFlags))
#  longest.filename <- max(nchar(main.titles))
#  filenames <- sub(" ","_", paste("Summary_", main.titles, sep=""))
#  
#  cat("Plotting Summary images\n")
#  pbar <- txtProgressBar(min=0, max=npages, char="*", width=20, style=3)
#  for (i in 1:npages) {
#    Sys.sleep(0.5)
#    if(i == 1) {k <- 1; k.max <- k + (ArraysPerPage-1)}
#    if(i > 1) {k <- k + ArraysPerPage; k.max <- k + (ArraysPerPage-1)}
#    if(k + (ArraysPerPage-1) >= totalArrays) {k.max <- totalArrays}
    
    if(npages > 1) {filename.suffix <- paste("_part_", i, ".png", sep="")} else {filename.suffix=".png"}
    
    for (j in 1:length(main.titles)) {
      png(paste(filenames[j], filename.suffix, sep=""), width=1600, height=1200, pointsize=20)
      par(las=2, mar=c(longest.filename/1.75,5,2,2))
      
      if (j == 1 | j == 2)
        if (RG$datatype == "both") {
          tmpObj <- rbind(RGsummary[,j],RGsummary[,j+2+hasRealBG])
          legend <- colnames(RG$summary)[c(j,j+2+hasRealBG)]
        } else {
          tmpObj <- RGsummary[,j]
          legend <- colnames(RG$summary)[j]
        }
      if (j == 3)
        if(hasRealBG)
          if (RG$datatype == "both") {
            tmpObj <- rbind(RGsummary[,j],RGsummary[,j+3])
            legend <- colnames(RG$summary)[c(j,j+3)]
          } else {
            tmpObj <- RGsummary[,j]
            legend <- colnames(RG$summary)[j]
          }
      if (j > 3 | (j == 3 & !hasRealBG)) {
        ext1 <- 2*(j>(length(main.titles)-length(RG$QC.vector)))
        ext2 <- (length(RG$QC.vector)-2)*(j>(length(main.titles)-length(RG$QC.vector)))
        tmpObj <- rbind(RGsummary[,j+2*(RG$datatype=="both")+(hasRealBG & RG$datatype=="both")+ext1],RGsummary[,j+2*(RG$datatype=="both")+(hasRealBG & RG$datatype=="both")+2+ext1+ext2])
        legend <- colnames(RG$summary)[c(j+2*(RG$datatype=="both")+(hasRealBG & RG$datatype=="both")+ext1,j+2*(RG$datatype=="both")+(hasRealBG & RG$datatype=="both")+2+ext1+ext2)]
      }
    if (j == length(main.titles) & hasManualFlags) {
      tmpObj <- tmpObj[2,]
    legend <- legend[2]
    }
    
      if (j == 1 | j == 2 | (j == 3 & hasRealBG == 1))
        col <- c(rep("salmon",RG$datatype!="green"),rep("springgreen3",RG$datatype!="red"))
      else
        col <- c("steelblue3","lightblue")
    
      if(!is.null(dim(tmpObj))) {
        plotData <- tmpObj[,k:k.max]
        labels <- colnames(tmpObj)[k:k.max]
      } else { #just one row
        plotData <- tmpObj[k:k.max]
        labels <- names(tmpObj)[k:k.max]
      }
      
      if(length(grep("percentage",tolower(main.titles[j])))>0) {
        y.lim <- c(0,112)
      } else {
        if(length(grep("number",tolower(main.titles[j])))>0) {
          if(RG$source!="genepix") {
            y.lim <- c(0,1.12*sum(RG$genes$ControlType==0))
          } else {
            y.lim <- c(0,1.12*sum(RG$genes$ControlType=="false"))
          }
        } else {
        y.lim <- c(0.92*min(tmpObj), 1.02*max(tmpObj)+0.1*(max(tmpObj)-min(tmpObj)))
        }
      }
    
      #barplot(height=plotData, beside=TRUE, xaxt="n", lend=1, lwd=700/length(plotData), ylim=y.lim, xlab="", main=main.titles[j], ylab=y.axis.labels[j], col=c("steelblue3","mediumblue"))
      barplot(height=plotData, beside=TRUE, main=main.titles[j], ylab=y.axis.labels[j], ylim=y.lim, xpd=FALSE, col=col,legend.text=legend)

      if(!is.null(dim(tmpObj))) {
        abline(h=mean(tmpObj[1,]),lty=2, col=col[1])
        abline(h=mean(tmpObj[2,]),lty=2, col=col[2])
      } else {
        abline(h=mean(tmpObj),lty=2, col=col[1])
      }
      #axis(1, labels=labels, at=1:length(labels))
      dev.off()
      setTxtProgressBar(pbar, i)
    }
    cat("\n")
  }
  cat("\n")
}

