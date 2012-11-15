# version 1.2.2

require(RColorBrewer)  ## Color Palette
require(affy)          ## plotDensity function
require(gplots)        ## heatmap.2 function
require(geneplotter)   ## densCols function
require(shape)         ## imageplot3by2Adp function
require(KernSmooth)    ## Required for barplots / MvA plots

##############################
#CreateSummaryPlots function##
##############################

CreateSummaryPlots <- function(RG) {
# This function is a bit more flexible in terms of available summary fields.
# - based on the datatype make a distinction between red/green images.
# - Updated the error detection text.
# - Improved removal of the 4 fixed columns (mean, median, min, max)
  error <- NULL
  if(is.null(RG$datatype)) { error <- c(error, "\n- a datatype field must be present in your data object with value \"green\", \"red\" or \"both\")" ) }
  if(is.null(RG$summary)) { error <- c(error, "\n- a summary field must be present in the data object, call CalculateSummary first") }

  if(!is.null(error)) {
    cat("[ERROR] CreateSummaryPlots could not be started due to the following errors:")
    stop(paste(error, "\n"))
  }

  datatype <- RG$datatype

  ## Remove the Minimum, Mean, MEdian and Maximum rows.
  temp <- c("Minimum", "Mean", "Median", "Maximum")
  temp2 <- rownames(RG$summary) %in% temp

  RGsummary <- RG$summary[!temp2,]
  rm(temp, temp2)

  totalArrays <- nrow(RGsummary)
  npages <- ceiling(totalArrays/16)
  x <- totalArrays/npages
  ArraysPerPage <- ceiling(x)
  hasRealBG <- !is.null(RG$other$rBGMeanSignal)|!is.null(RG$other$gBGMeanSignal)
  hasManualFlags <- !is.null(RG$other$manualFlags)

  main.titles <- c("Foreground","Background",rep("Real Background",hasRealBG),"Number of Good Probes","Percentage of Good Probes",RG$QC.vector,rep("Manually Flagged",hasManualFlags))
  y.axis.labels <- c(rep("Mean Signal",2+hasRealBG),"Number of reporters","% of reporters",rep("Number of reporters",length(RG$QC.vector)+hasManualFlags))
  longest.filename <- max(nchar(main.titles))
  filenames <- sub(" ","_", paste("Summary_", main.titles, sep=""))
  
  cat("Plotting Summary images\n")
  pbar <- txtProgressBar(min=0, max=npages, char="*", width=20, style=3)
  for (i in 1:npages) {
#   Sys.sleep(0.5)
    if(i == 1) {k <- 1; k.max <- k + (ArraysPerPage-1)}
    if(i > 1) {k <- k + ArraysPerPage; k.max <- k + (ArraysPerPage-1)}
    if(k + (ArraysPerPage-1) >= totalArrays) {k.max <- totalArrays}
    
    if(npages > 1) {filename.suffix <- paste("_part_", i, ".png", sep="")} else {filename.suffix=".png"}
    
    for (j in 1:length(main.titles)) {
      ## IF the main.title name is not present in the object, the data does not exist and should be skipped.
      if(length(grep(main.titles[j], colnames(RGsummary))) == 0 ) { next }
      png(paste(filenames[j], filename.suffix, sep=""), width=1600, height=1200, pointsize=20)
      par(las=2, mar=c(longest.filename/1.75,5,2,2), cex.axis=0.75)
      
      ## Foreground / Background plots
      if (j == 1 | j == 2) {
        if (RG$datatype == "both") {
          tmpObj <- rbind(RGsummary[,j],RGsummary[,j+2+hasRealBG])
          legend <- colnames(RG$summary)[c(j,j+2+hasRealBG)]
        } else {
          tmpObj <- RGsummary[,j]
          legend <- colnames(RG$summary)[j]
        }
      }
      ## Plot the realBG column if present
      if (j == 3 & hasRealBG == 1) {
        if (RG$datatype == "both") {
          tmpObj <- rbind(RGsummary[,j],RGsummary[,j+3])
          legend <- colnames(RG$summary)[c(j,j+3)]
        } else {
          tmpObj <- RGsummary[,j]
          legend <- colnames(RG$summary)[j]
        }
      }
      if ( j > 3 | (j == 3 & !hasRealBG == 1)) {
        ## check first if the name exists in your summary object. If not, skip image.
        ext1 <- 2*(j>(length(main.titles)-length(RG$QC.vector)))
        ext2 <- (length(RG$QC.vector)-2)*(j>(length(main.titles)-length(RG$QC.vector)))
        tmpObj <- rbind(RGsummary[,j+2*(RG$datatype=="both")+(hasRealBG & RG$datatype=="both")+ext1],RGsummary[,j+2*(RG$datatype=="both")+(hasRealBG & RG$datatype=="both")+2+ext1+ext2])
        legend <- colnames(RG$summary)[c(j+2*(RG$datatype=="both")+(hasRealBG & RG$datatype=="both")+ext1,j+2*(RG$datatype=="both")+(hasRealBG & RG$datatype=="both")+2+ext1+ext2)]
      }
      if (j == length(main.titles) & hasManualFlags) {
        tmpObj <- tmpObj[2,]
        legend <- legend[2]
      }
    
      if (j == 1 | j == 2 | (j == 3 & hasRealBG == 1)) {
        col <- c(rep("salmon",RG$datatype!="green"),rep("springgreen3",RG$datatype!="red"))
      } else {
        col <- c("steelblue3","lightblue")
      }

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
      barplot(height=plotData, beside=TRUE, main=main.titles[j], ylab=y.axis.labels[j], ylim=y.lim, xpd=FALSE, col=col,legend.text=legend, cex.lab=0.75)
#      axis(1, labels=labels, at=1:length(labels))

      if(!is.null(dim(tmpObj))) {
        abline(h=mean(tmpObj[1,]),lty=2, col=col[1])
        abline(h=mean(tmpObj[2,]),lty=2, col=col[2])
      } else {
        abline(h=mean(tmpObj),lty=2, col=col[1])
      }
      dev.off()
      setTxtProgressBar(pbar, i)
    }
    cat("\n")
  }
  cat("\n")
}




##############################
##imageplot3by2Adp function ##
##############################

## Adapted imageplot3by2 function

imageplot3by2Adp <- function (data.object, which.field, name, high, low, log.transform=FALSE, symm=FALSE, orientation=NULL, debug=0) {
  if(!class(data.object) %in% c("RGList","MAList","EListRaw")) stop("Object is not of RGList, ElistRaw, or MAList class")
    
  if(is.null(data.object$printer))
    stop("the printer field in the full object is empty, whereas this information is required for image plots")

  if(debug==1) {
    cat("Paste the following below to debug:\n")
    cat(paste("data.object <-", deparse(substitute(data.object)), "\nwhich.field <-", deparse(substitute(which.field)), "\n"))
    cat(paste("name <-", deparse(substitute(name)), "\nhigh <-", deparse(substitute(high)), "\nlow <-", deparse(substitute(low)), "\n"))
    cat(paste("log.transform <-", deparse(substitute(log.transform)), "\nsymm <-", deparse(substitute(symm)), "\norientation <-", deparse(substitute(orientation)), "\n"))
    stop("Happy debugging!")
  }
  cat(paste("    --", paste(rep("-", nchar(deparse(substitute(name)))), collapse=""), "--\n", sep="", collapse=""))
  cat(paste("    - ", deparse(substitute(name)), " -\n", sep="", collapse=""))
  cat(paste("    --", paste(rep("-", nchar(deparse(substitute(name)))), collapse=""), "--\n", sep="", collapse=""))

  cat("   * Preprocessing... ")
  prefix <- paste("VirtualArray", name, sep = "-")
  narrays <- ncol(which.field)
  npages <- ceiling(narrays/6)
  cnames <- colnames(which.field)

  #Calculating z-ranges and their mean and stdev
  minima <- NULL
  maxima <- NULL

  nr.of.well.above <- NULL
  nr.of.saturated <- NULL

  if(log.transform) {plot.field <- log(which.field,2)} else {plot.field <- which.field}
  ## Check if object contains binary values, predefined values or regular data:
  if((sum(!(names(table(plot.field)) %in% c(0,1))) > 0)) {
    binary <- FALSE
    #check whether only values from (0,-3,-5,-8) are present
    if(sum(names(table(plot.field)) %in% c(0,-3,-5,-8)) == length(names(table(plot.field)))) {
      QCField <- TRUE
      values <- as.numeric(names(table(plot.field)))
      number <- apply(plot.field,2,function(f) sum(f!=0,na.rm=TRUE))
      mean.number <- mean(number)
      stdev.number <- sd(number)
    } else {
      QCField <- FALSE      
      minima <- apply(plot.field,2,min,na.rm=TRUE)
      maxima <- apply(plot.field,2,max,na.rm=TRUE)
      mean.min <- mean(minima)
      stdev.min <- sd(minima)
      mean.max <- mean(maxima)
      stdev.max <- sd(maxima)
    }
  } else {
    binary <- TRUE
    QCField <- FALSE
    number <- apply(plot.field,2,sum,na.rm=TRUE)
    mean.number <- mean(number)
    stdev.number <- sd(number)
  }

  #Adjusting object for optimal visualisation (landscape orientation)
  printer.info <- data.object$printer
  printer.c <- printer.info$nspot.c
  printer.r <- printer.info$nspot.r

  if((printer.c > printer.r) & (printer.info$ngrid.c == 1) & (printer.info$ngrid.r == 1)) {
     printer.info$nspot.c <- printer.r
     printer.info$nspot.r <- printer.c
     plot.field <- plot.field[order(-data.object$genes[,"Col"], data.object$genes[,"Row"], decreasing=TRUE),]
   }

  #Plotting the array images
  cat("  ok.\n   * Plotting images:\n")
  pbar <- txtProgressBar(min=0, max=npages, char="*", width=20, style=3)
  dimensions <- c((14 * printer.info$nspot.c * printer.info$ngrid.c), (14 * printer.info$nspot.r * printer.info$ngrid.r))

  for (ipage in 1:npages) {
    i1 <- ipage * 6 - 5
    i2 <- min(ipage * 6, narrays)
    fn <- paste(prefix, "-", i1, "-", i2, ".png", sep = "")
    ## Adjusted size of image, dependant on the data dimensions
#    png(filename=fn, width=8.5 * 140, height=10 * 140, pointsize=20)
    png(filename=fn, width=dimensions[2], height=dimensions[1], pointsize=20)
    layout(rbind(
       matrix(data=c(rep(1,4),rep(2,4)), nrow=3, ncol=8, byrow=T),
       matrix(data=c(rep(3,4),rep(4,4)), nrow=3, ncol=8, byrow=T),
       matrix(data=c(rep(5,4),rep(6,4)), nrow=3, ncol=8, byrow=T)))

    for (i in i1:i2) {
      col.main <- "black"
      if (binary | QCField) {
        if (stdev.number > 0) {
          if(sum(plot.field[,i]!=0, na.rm=TRUE) <= mean.number-2*stdev.number) {col.main <- "red"}
          if(sum(plot.field[,i]!=0, na.rm=TRUE) >= mean.number+2*stdev.number) {col.main <- "red"}
        }
      } else { 
        if (stdev.max > 0 | stdev.min > 0) {
          if(min(plot.field[,i], na.rm=TRUE) <= mean.min-2*stdev.min) {col.main <- "red"}
          if(min(plot.field[,i], na.rm=TRUE) >= mean.min+2*stdev.min) {col.main <- "red"}
          if(max(plot.field[,i], na.rm=TRUE) <= mean.max-2*stdev.max) {col.main <- "red"}
          if(max(plot.field[,i], na.rm=TRUE) >= mean.max+2*stdev.max) {col.main <- "red"}                                                  
        }
      }
  
      limit <- 0.05
      if(symm) {
        zlim=c(-max(abs(quantile(plot.field[,i], c(limit,1-limit), na.rm=TRUE))),max(abs(quantile(plot.field[,i], c(limit,1-limit), na.rm=TRUE))))
      } else {
        zlim=quantile(plot.field[,i], c(limit,1-limit), na.rm=TRUE)
      }
  
      if(binary) zlim <-c(0, 1)
      if(QCField) zlim <- c((data.object$datatype=="green") * -3 + (data.object$datatype=="red") * -5 + (data.object$datatype=="both") * -8, 0)
      par(mar=(c(2,3,3,5) + 0.1))
      imageplot(plot.field[,i], printer.info, zlim=zlim, mar=c(2, 2, 4, 4), main=cnames[i], high=high, low=low, legend=TRUE, col.main=col.main, zerocenter=symm)
      # emptyplot()
      ColorFunction <- colorRampPalette(c(low,high))

      cex.size <- 0.6
      if( dimensions[1] > (3 * dimensions[2])) {
        pos.x <- c(0.87, 0.90)
        pos.y <- c(0.01, 0.99)
      } else {
        pos.x <- c(0.95, 0.96)  #0.91, 0.93
        pos.y <- c(0.01, 0.94)
      }

      if(binary) {
        colors2use <- ColorFunction(2)
        colorlegend(col=colors2use, zlim, left=FALSE, zval=c(zlim[1], zlim[2]), posy=pos.y, posx=pos.x, cex=cex.size, digit=1)
      } else {
        if(QCField) {
          colors2use <- ColorFunction(length(values))
          colorlegend(col=colors2use, zlim, left=FALSE, zval=values, posy=pos.y, posx=pos.x, cex=cex.size, digit=1)
        } else {
          colors2use <- ColorFunction(100)
          colorlegend(col=colors2use, zlim, left=FALSE, zval=c(zlim[1], mean(zlim), zlim[2]), posy=pos.y, posx=pos.x, cex=cex.size, digit=1)
        }
      }
    }
## Have to check if this is the proper position.
  setTxtProgressBar(pbar, ipage)
  dev.off()
  }
  cat("\n\n")
}

##############################
## Hierarchcluster function ##
##############################

HierarchCluster <- function(x, dist.method = "euclidean", clust.method = "ward", main, data.return=FALSE, ...) {
  a <- dist(t(x), method=dist.method)
  b <- hclust(a, method=clust.method)
  png(paste("Cluster_",main,".png", sep=""), width=1600, height=1200, pointsize=20)
  if(dim(x)[2] >= 100) par(cex=0.8, cex.axis=1.25, cex.lab=1.25, cex.main=1.5, cex.sub=1.25)
  plot(b, main=main, ...)
  dev.off()
  if(data.return==TRUE) {return(b)}
}

##############################
##      Heatmap functie     ##
##############################

CreateHeatMap <- function(data, main=NULL) {
  no.title <- ifelse(is.null(main), 1, 0)
  
  my.dist <- function(x) dist(x, method="euclidean")
  my.hclust <- function(d) hclust(d, method="ward")

  if(class(data)!= "matrix") {
    if(data$datatype=="both") {
      if(no.title == 1) { main1=paste("Heatmap Logratio ", deparse(substitute(data)), sep="") }
      if(no.title == 0) { 
        if(length(main)>1) { stop("main variable has more than 1 argument!\n") }; 
        main1 <- paste("Heatmap Logratio - ", main, sep="") 
      }

      png(paste("Heatmap Logratio ", deparse(substitute(data)),".png", sep=""), width=2000, height=1400, pointsize=18)
      crp <- cor(data$M, use="complete.obs")
      heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", margins=c(20,20), symm=TRUE, density.info="density", main=main1, dendrogram="row")

    #  heatmap(crp, sym=TRUE, hclustfun=function(w1) hclust(w1, method="ward"), main=main)
    #  legend("topleft",c(paste("min =",round(min(crp),2)),paste("max =",max(crp))),col=c("red","beige"),pch=15, cex=1.2)
      dev.off()

      if(no.title == 1) { main2 <- paste("Heatmap Average Intensity ", deparse(substitute(data)), sep="") }
      if(no.title == 0) { main2 <- paste("Heatmap Average Intensity - ", main, sep="") }

      png(paste("Heatmap Average Intensity ", deparse(substitute(data)),".png", sep=""), width=2000, height=1400, pointsize=18)
      crp <- cor(data$A, use="complete.obs")
      heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", margins=c(20,20), symm=TRUE, density.info="density", main=main2, dendrogram="row")
    #  heatmap(crp, sym=TRUE, hclustfun=function(w1) hclust(w1, method="ward"))
    #  legend("topleft",c(paste("min =",round(min(crp),2)),paste("max =",max(crp))),col=c("red","beige"),pch=15, cex=1.2)
      dev.off()
    } else {
      if(no.title == 1) { 
        main2 <- paste("Heatmap Estimated Intensity ", deparse(substitute(data)), sep="") }
      if(no.title == 0) { main2 <- paste("Heatmap Estimated Intensity - ", main, sep="") }

      png(paste("Heatmap Estimated Intensity ", deparse(substitute(data)),".png", sep=""), width=2000, height=1400, pointsize=18)
      crp <- cor(data$other$EST, use="complete.obs")
      heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", margins=c(20,20), symm=TRUE, density.info="density", main=main2, dendrogram="row")
      dev.off()
    }
  } else { # class is matrix
      if(no.title == 1) { main2 <- paste("Heatmap Estimated Intensity ", deparse(substitute(data)), sep="") }
      if(no.title == 0) { main2 <- paste("Heatmap Estimated Intensity - ", main, sep="") }

      png(paste("Heatmap Estimated Intensity ", deparse(substitute(data)),".png", sep=""), width=2000, height=1400, pointsize=18)
      crp <- cor(data, use="complete.obs")
      heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", margins=c(20,20), symm=TRUE, density.info="density", main=main2, dendrogram="row")
      dev.off()
  }
}

##############################
## CreatePCAplot function   ##
##############################

CreatePCAplot <- function(data, scaled_pca=TRUE, namesInPlot=FALSE){
  # PCA performed on reporters NOT containing NAs only
  # Scaled PCA by default
  
  if(class(data)!= "matrix") {
    if(data$datatype == "both") {
      pcaData <- data$M
    } else {
      pcaData <- data$other$EST
    }
  } else {
    pcaData <- data
  }

  pca1 <- NULL  
  try(pca1 <- prcomp(t(pcaData[apply(pcaData,1,function(r) {sum(is.na(r)) == 0}),]), retx=T, center=T, scale=scaled_pca),TRUE)
  if(!is.null(pca1)) {
    perc_expl1 <- round(((pca1$sdev[1:3]^2)/sum(pca1$sdev^2))*100,2)

    tmain <- paste("PCA analysis of", deparse(substitute(data)))
    plotColors <- rainbow(dim(pcaData)[2])
  arrayNames <- as.character(colnames(pcaData))
    cex.circle <- 1.3
    cex.text <- 0.7
    tcol <- "#444444"

    png(file = paste("PCAplot",deparse(substitute(data)),".png",sep=""), width=1200+(!namesInPlot)*400, height=1200, pointsiz=25)
    if(!namesInPlot) {
      layout(rbind(c(1,1,2,2,5),c(3,3,4,4,5)))
      par(cex.main=1.8, cex.lab=1.4, cex.axis=1.4)
    } else {
      layout(rbind(c(1,2),c(3,4)))
    }
  plot(pca1$x[,1],pca1$x[,2],cex=cex.circle,pch=17:0,
        col=plotColors,xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
        ylab=paste("PC2 (",perc_expl1[2],"%)",sep=""))
    if(namesInPlot) text(pca1$x[,1],pca1$x[,2], arrayNames,pos=4,cex=cex.text, 
        col=tcol)
    title(main=tmain, font=2)  
    plot(pca1$x[,1],pca1$x[,3],cex=cex.circle,pch=17:0,
        col=plotColors,xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
        ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
    if(namesInPlot) text(pca1$x[,1],pca1$x[,3], arrayNames,pos=4,cex=cex.text, 
        col=tcol)
    plot(pca1$x[,2],pca1$x[,3],cex=cex.circle,pch=17:0,
        col=plotColors,xlab=paste("PC2 (",perc_expl1[2],"%)",sep=""),
        ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
    if(namesInPlot) text(pca1$x[,2],pca1$x[,3], arrayNames,pos=4,cex=cex.text, 
        col=tcol)
    barplot((pca1$sdev^2)/sum(pca1$sdev^2),xlab="components",
        ylab="% of total variance explained")
    if(!namesInPlot) {
      plot(1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
      legend("topright",arrayNames,pch=17:0,col=plotColors,cex=cex.circle)
    }
    dev.off()
  } else {
    cat("Warning: pca on the",deparse(substitute(data)),"data set unsuccessful, possibly due to lack of memory\n\n")
  }
}

##############################
## CreateMAplots function   ##
##############################

CreateMAplots <- function(first=NULL, second=NULL, third=NULL, fourth=NULL, labels=NULL, postfix=NULL) {
  cat("* Preprocessing MA Plots ...")
  # Check if the user passed in the correct arguments
  which.norm <- !c(is.null(first),is.null(second),is.null(third),is.null(fourth))
  if(sum(which.norm) == 0) stop("At least one data object - MAList or matrix - has to be provided")
  which.class <- c(class(first),class(second),class(third),class(fourth))
  if(sum(which.class == "MAList") < sum(which.norm)) stop("only data objects of class MAList can be provided")
  if(length(labels) != sum(which.norm)) stop("The number of data object and labels to describe the type of normalization has to be equal")
  
  # Create subtexts for the graphs 
  subtexts1 <- paste(labels, " - Non-zero-weighted points (filtered)", sep="")
  subtexts2 <- paste(labels, " - All Points (not filtered)", sep="")

  #Create identifier for the normalization objects used
  objects <- c("first","second","third","fourth")[which.norm]

  #determine x and y limits and check whether objects have a weight field
  xlim <- c(0,0)
  ylim <- c(0,0)
  plot.weighted <- rep(FALSE, length(objects))
  for (j in 1:length(objects)) {
    rng <- range(as.matrix(get(objects[j])$A), na.rm=TRUE)
    if (rng[1] < xlim[1]) xlim[1] <- rng[1]
    if (rng[2] > xlim[2]) xlim[2] <- rng[2]
    rng <- range(as.matrix(get(objects[j])$M), na.rm=TRUE)
    if (rng[1] < ylim[1]) ylim[1] <- rng[1]
    if (rng[2] > ylim[2]) ylim[2] <- rng[2]
    if(!is.null(get(objects[j])$weights)) plot.weighted[j] <- TRUE
  }
  
  #plotting
  maxArrays <- ncol(get(objects[which.norm][1]))
  cat("  ok.\n* MA plotting progress:\n")
  pbar <- txtProgressBar(min=0, max=maxArrays, char="*", width=20, style=3)

  for (i in 1:maxArrays) {
    titleName <- colnames(get(objects[which.norm][1]))[i]
    
    if(sum(plot.weighted) > 0) {
      fileName = paste("MA_Plots_-_", titleName, "_-_Zero_Weighted", postfix, ".png", sep="", collapse=NULL)
      png(filename = fileName, width=1600, height=1200, pointsize=20)
      par(mar=c(8,4,4,3)+0.1, mfrow=c(1+(sum(which.norm) > 1),1+(sum(which.norm) > 2)))
      
      # Plot a MA.Plot (Zero-weighted) only when the normalization method is present in the which.norm vector and weights are present
      for (j in 1:length(objects)) {
        if(plot.weighted[j]) {
          # Putting A and M values in a matrix. Setting zero-weighted M-values to NA before plotting
          a <- as.matrix(get(objects[j])$A)[,i]
          m <- as.matrix(get(objects[j])$M)[,i]
          w <- as.matrix(get(objects[j])$weights)[,i]
          k <- is.na(w) | (w <= 0)
          m[k] <- NA
          
          #Determining color densities and plotting the MA-plot
          color <- densCols((get(objects[j])$A[,i]), (get(objects[j])$M[,i]))
          plot(a, m, xlim=xlim, ylim=ylim, col=color, cex=0.3, pch=16, xlab="A", ylab="M")
          title(main=titleName, sub=subtexts1[j], outer=FALSE)
          abline(0,0,col="blue")
        } else {
          cat(paste("Not plotting MA-Plot excluding zero-weighted spots for", deparse(substitute(get(objects[j]))), "- Reason: no weights available for this data object \n"))
          plot(0,type='n',xaxt='n',yaxt='n',xlab="",ylab="", bty='n')
        }
      }
      dev.off()
    }
    
    fileName = paste("MA_Plots_-_", titleName,"_-_All_Points", postfix, ".png", sep="", collapse=NULL)
    png(filename = fileName, width=1600, height=1200, pointsize=20)
    par(mar=c(8,4,4,3)+0.1, mfrow=c(1+(sum(which.norm) > 1),1+(sum(which.norm) > 2)))
    
    # Plot a MA.Plot (all points) only when the normalization method is present in the which.norm vector    
    for (j in 1:length(objects)) {
      # Putting A and M values in a matrix, to be consistent with the method in the zero-weighted plots 
      a <- as.matrix(get(objects[j])$A)[,i]
      m <- as.matrix(get(objects[j])$M)[,i]
      
      #Determining color densities and plotting the MA-plot
      color <- densCols((get(objects[j])$A[,i]), (get(objects[j])$M[,i]))
      plot(a, m, xlim=xlim, ylim=ylim, col=color, cex=0.3, pch=16, xlab="A", ylab="M")
      title(main=titleName, sub=subtexts2[j], outer=FALSE)
      abline(0,0,col="blue")
    }  
    dev.off()
    setTxtProgressBar(pbar, i)
    
    if(i == maxArrays) {cat("\n")}
  }
}


#################################
## CreateDensityPlots function ##
#################################

CreateDensityPlots <- function(x, name=NULL) {

  if((!class(x) %in% c("RGList","EListRaw","MAList")) & (!is.matrix(x))) stop("Only object of class RGList, EListRaw, MAList, or matrix can be handled")
  if(is.null(name)) stop("a name label must be provided")

  # Preparation of the number of pages needed
  maxArrays <- ncol(x)
  npages <- ceiling(maxArrays/9)
  y <- maxArrays/npages
  ArrayPerPage <- ceiling(y)
  
  # Loop that calls the plotting
  cat("Plotting images ")
  for (i in 1:npages) {
    ## k and k.max are the counters, with k being the start position and
    ## k.max being the end position
    if(i == 1) {k <- 1; k.max <- k + (ArrayPerPage-1)}
    if(i > 1) {k <- k + ArrayPerPage; k.max <- k + (ArrayPerPage-1)}
    if(k + (ArrayPerPage-1) >= maxArrays) { k.max <- maxArrays}
    
    if(npages > 1)
      outputname <- paste("DensityPlots_", name, "_Part_", i, ".png", sep="")
    else
      outputname <- paste("Density_Plots_", name, ".png")

    if(k == k.max) {
      PlotDensities(x, position=c(k,k.max), outputName=outputname, name=name)
    } else {
      PlotDensities(x, position=c(k,k.max), outputName=outputname, name=name)
    }

    cat(".")
    
  }
  cat(" ok.\n")
} 

##############################
##  PlotDensities function  ##
##############################

PlotDensities <- function (x, position, outputName=NULL, name=NULL) {
  maxArrays <- position[2]-position[1]+1
  if(is.null(outputName)) {stop(cat("PlotDensities() error. Please define the outputName variable in the function\n\n"))}
  if(is.null(name)) {stop(cat("PlotDensities() error. Please define the name variable in the function\n\n"))}
  
  labels.green <- NULL
  labels.red <- NULL
  labels.matrix <- NULL
  
  for (i in position[1]:position[2]) {
    arrayName <- colnames(x)[i]
    labels.green <- c(labels.green, paste(arrayName, "Green Signal", sep=" "))
    labels.red <- c(labels.red, paste(arrayName, "Red Signal", sep=" "))
    labels.matrix <- c(labels.matrix, arrayName)
  }
    
  greens <- brewer.pal(9, "Greens")[-1]
  reds <- brewer.pal(9, "Reds")[-1]
  blues <- brewer.pal(9, "Blues")[-1]
  extra.color <- brewer.pal(9, "RdBu")[-1]
    
  colors.green <- NULL
  colors.red <- NULL
  colors.blue <- NULL
  
  for (i in 1:(maxArrays)) {
    if(i > 8 ) {
      k <- 5
      l <- i-8
      colors.green <- c(colors.green, extra.color[k+l])
      colors.red   <- c(colors.red,   extra.color[k-l])
      colors.blue   <- c(colors.blue,   extra.color[k+l])
    } else {
      colors.green <- c(colors.green, greens[i])
      colors.red   <- c(colors.red,   reds[i])
      colors.blue   <- c(colors.blue,   blues[i])
    }
  }
  
  ## Preparing to Write to .png file
  
  background <- NULL
  background2 <- NULL
  corr.foreground <- NULL
  foreground <- NULL
  
  if(is.matrix(x)){
    #one channel data, but we don't know which channel, thus use blue for the plots
    color <- colors.blue
    label <- c(labels.matrix)
    corr.foreground <- x
    selection <- position[1]:position[2]
  } else {
    if(x$datatype == "both") {
      color <- c(colors.green,colors.red)
      label <- c(labels.green,labels.red)
      if(class(x) == "RGList") {
        if(!is.null(x$other$gBGMeanSignal) & !is.null(x$other$rBGMeanSignal))
          background2 <- cbind(log(x$other$gBGMeanSignal,2),log(x$other$rBGMeanSignal,2))
        foreground <- cbind(log(x$G,2),log(x$R,2))
        background <- cbind(log(x$Gb,2),log(x$Rb,2))
      } else {
        transform <- RG.MA(x)
        corr.foreground <- cbind(log(transform$G,2),log(transform$R,2))
      }
      selection <- c(position[1]:position[2],position[1]:position[2]+ncol(x))
    } else {
      if (x$datatype=="red") {
        color <- colors.red
        label <- labels.matrix
        if(class(x) == "EListRaw") {
          if(!is.null(x$other$rBGMeanSignal)) background2 <- log(x$other$rBGMeanSignal,2)
          foreground <- log(x$E,2)
          background <- log(x$Eb,2)
        } else {
          corr.foreground <- x$other$EST
        }
        selection <- position[1]:position[2]
      }
      if (x$datatype=="green") {
        color <- colors.green
        label <- labels.matrix
        if(class(x) == "EListRaw") {
          if(!is.null(x$other$gBGMeanSignal)) background2 <- log(x$other$gBGMeanSignal,2)
          foreground <- log(x$E,2)
          background <- log(x$Eb,2)
        } else {
          corr.foreground <- x$other$EST
        }
        selection <- position[1]:position[2]
      }
    }
  }

  png(file=outputName, width=1600, height=600+600*(is.null(corr.foreground)), pointsize=20)
  
  if(class(x) == "RGList" || class(x) == "EListRaw") {
    par(mfrow=c(2,2))
    if(x$source == "agilent") { 
      # Limits for the Spatial Detrend Signal Distribution
      if(!is.null(background2)) {
        x.min.sd <- 0
        x.max.sd <- ceiling(max(sapply(apply(background2, 2, density, na.rm=TRUE), function(z) max(z$x))))
        y.min.sd <- 0
        y.max.sd <- ceiling(max(sapply(apply(background2, 2, density, na.rm=TRUE), function(z) max(z$y))))
        # Plotting R/G Background distribution and setting name of main.title to Spatial Detrend Distibution, because this is made when plotting Gb and Rb
        plotDensity(as.matrix(background2[,selection]), col=color, xlim=c(x.min.sd,x.max.sd), lty=1:length(selection), lwd=2, ylim=c(y.min.sd,y.max.sd), main="Background Distribution", xlab=paste("Log2(Intensity) -",name), ylab="Density")
      } else {
        plot(0,type='n',xaxt='n',yaxt='n',xlab="",ylab="",bty='n')
      }
      main.title <- "Spatial Detrend Distribution"
    } else { 
      main.title <- "Background Distribution"
    }
    # Limits for the Foreground Signal Distribution
    if(!is.null(foreground)) {
      x.min.fg <- 0
      x.max.fg <- ceiling(max(sapply(apply(foreground, 2, density, na.rm=TRUE), function(z) max(z$x))))
      y.min.fg <- 0
      y.max.fg <- ceiling(max(sapply(apply(foreground, 2, density, na.rm=TRUE), function(z) max(z$y))))
      plotDensity(as.matrix(foreground[,selection]), col=color, xlim=c(x.min.fg,x.max.fg), lty=1:length(selection), lwd=2, ylim=c(y.min.fg,y.max.fg), main="Foreground Signal Distribution", xlab=paste("Log2(Intensity) -",name), ylab="Density")
    }
    
    # Limits for the Background Signal Distribution
    if(!is.null(background)) {
      x.min.bg <- 0
      x.max.bg <- ceiling(max(sapply(apply(background, 2, density, na.rm=TRUE), function(z) max(z$x))))
      y.min.bg <- 0
      y.max.bg <- ceiling(max(sapply(apply(background, 2, density, na.rm=TRUE), function(z) max(z$y))))
      plotDensity(as.matrix(background[,selection]), col=color, xlim=c(x.min.bg,x.max.bg), lty=1:length(selection), lwd=2, ylim=c(y.min.bg,y.max.bg),main=main.title, xlab=paste("Log2(Intensity) -",name), ylab="Density")
    }
  } else {
    par(mfrow=c(1,2))
    if(!is.null(corr.foreground)) {
      x.min.fg <- 0
      x.max.fg <- ceiling(max(sapply(apply(corr.foreground, 2, density, na.rm=TRUE), function(z) max(z$x))))
      y.min.fg <- 0
      y.max.fg <- ceiling(max(sapply(apply(corr.foreground, 2, density, na.rm=TRUE), function(z) max(z$y))))
      plotDensity(as.matrix(corr.foreground[,selection]), col=color, lty=1:length(selection), lwd=2, xlim=c(x.min.fg,x.max.fg), ylim=c(y.min.fg,y.max.fg), main="BG Corrected Foreground Signal Distribution", xlab=paste("Log2(Intensity) -",name), ylab="Density")
    }
  }
  
  plot.new()
  #if(is.null(corr.foreground)) plot.new()
  if((class(x) == "RGList" || class(x) == "EListRaw") & is.null(background2)) plot.new()

  legend("topright",label,col=color, lty=1:length(selection),text.col="black",cex=1-0.2*(!is.null(corr.foreground)), lwd=3)
  dev.off()
}

###################################
##    CreateBoxplot function     ##
###################################

## class1 contains the actual data you want to visualise (e.g. Agilent.RG, Agilent.MA.LOESS, etc)

CreateBoxplot <- function(class1, fileName, figTitle, y.axis="M", use.col=NULL, y.lim=NULL, non.zero.weight=FALSE, weights=NULL) {
  
  if((!class(class1) %in% c("RGList","EListRaw","MAList")) & (!is.matrix(class1))) stop("Only object of class RGList, EListRaw, MAList, or matrix can be handled")
  
  if(non.zero.weight == TRUE) { 
    a <- naZeroWeights(class1, weights=weights)
    suffix <- c("_High_Quality_Only")
  } else { 
    a <- class1 
    suffix <- c("_All_Reporters")
  }
  columnName <- colnames(class1)
  longest.filename <- max(nchar(colnames(class1)))
  
  if(class(class1)=="RGList" || class(class1)=="EListRaw") {
    if(class(class1)=="RGList") {
      a$R <- log(a$R,2)
      a$G <- log(a$G,2)
      ylim.min <- min(a$R, a$G, na.rm=TRUE)
      ylim.max <- max(a$R, a$G, na.rm=TRUE)
    } else {
      if(class1$datatype=="red") {
        a$R <- log(a$E,2)
        ylim.min <- min(a$R, na.rm=TRUE)
        ylim.max <- max(a$R, na.rm=TRUE)
      } else {
        a$G <- log(a$E,2)
        ylim.min <- min(a$G, na.rm=TRUE)
        ylim.max <- max(a$G, na.rm=TRUE)
      }
    }
    if(class1$datatype!="green") {
      fileName2 <- paste(fileName,"_Raw_Red_Signal", suffix, ".png", sep="")
      FigTitle2 <- paste(figTitle," Raw Red Signal", suffix, ".png", sep="")
      png(filename=fileName2, width=1600, height=1200, pointsize=20)
      par(mar = c(longest.filename/1.75,4,4,3)+0.1, las=2)
      if(dim(class1)[2] >= 100)
        par(cex.axis=0.7)
      else
        par(cex.axis=0.8)
      boxplot(as.data.frame(a$R), names=rep("",length(columnName)), ylim=c(ylim.min,ylim.max), ylab = "M", col=use.col)
      title(main=FigTitle2, cex.main=0.9, col.main="blue", font.main=4)
      axis(1,at=1:length(columnName),labels=columnName)
      abline(0,0,col="blue")
      dev.off()
    }
    if(class1$datatype!="red") {
      fileName2 <- paste(fileName,"_Raw_Green_Signal", suffix, ".png", sep="")
      FigTitle2 <- paste(figTitle," Raw Green Signal", suffix, ".png", sep="")
      png(filename=fileName2, width=1600, height=1200, pointsize=20)
      par(mar = c(longest.filename/1.75,4,4,3)+0.1, las=2)
      if(dim(class1)[2] >= 100)
        par(cex.axis=0.7)
      else
        par(cex.axis=0.8)
      boxplot(as.data.frame(a$G), names=rep("",length(columnName)), ylim=c(ylim.min,ylim.max), ylab = "M", col=use.col)
      title(main=FigTitle2, cex.main=0.9, col.main="blue", font.main=4)
      axis(1,at=1:length(columnName),labels=columnName)
      abline(0,0,col="blue")
      dev.off()
      }
  } else {
    if(class(class1)=="MAList") {
      if(class1$datatype=="both") {
        if(y.axis == "M") { 
          fileName <- paste(fileName, "_M", suffix, ".png", sep="") 
  #         Plot.This <- a$M~col(as.matrix(a$M))
          Plot.This <- as.data.frame(a$M)
          y.lab <- "M"
          lheight <- 0
          if(is.null(y.lim)) {  y.lim <- c(-15,15) }
        }
        if(y.axis == "A") { 
          fileName <- paste(fileName,"_A", suffix, ".png", sep="") 
  #         Plot.This <- a$A~col(a$A)
          Plot.This <- as.data.frame(a$A)
          y.lab <- "A"
          lheight <- mean(a$A,na.rm=TRUE)
          if(is.null(y.lim)) { y.lim <- c(-10,25) }
        }
      } else {
          fileName <- paste(fileName, "_EST", suffix, ".png", sep="") 
  #         Plot.This <- a$other$EST~col(as.matrix(a$other$EST))
          Plot.This <- as.data.frame(a$other$EST)
          lheight <- mean(a$other$EST,na.rm=TRUE)
          y.lab <- "estimated expression"
          if(is.null(y.lim)) { y.lim <- c(-10,25) }
      }
      png(filename=fileName, width=1600, height=1200, pointsize=20)
      par(mar=c(longest.filename/2.1,4,4,3)+0.1, las=2)
      if(dim(class1)[2] >= 100)
        par(cex.axis=0.7)
      else
        par(cex.axis=0.8)
      boxplot(Plot.This, names=rep("",length(columnName)), ylab = y.lab, cex.axis=0.8, ylim=y.lim, col=use.col)
      title(main=figTitle, cex.main=0.9, col.main="blue", font.main=4)
      axis(1,at=1:length(columnName),labels=columnName)
      abline(lheight, 0, col="blue", lty=2)
      dev.off()
    } else {
      if(class(class1)=="matrix") {
        fileName <- paste(fileName, "_EST", suffix, ".png", sep="") 
  #     Plot.This <- a~col(as.matrix(a))
        Plot.This <- as.data.frame(a)
        lheight <- mean(a,na.rm=TRUE)
        y.lab <- "estimated expression"
        if(is.null(y.lim)) { y.lim <- c(-10,25) }
        png(filename=fileName, width=1600, height=1200, pointsize=20)
        par(mar=c(longest.filename/2.1,4,4,3)+0.1, las=2)
        if(dim(class1)[2] >= 100)
          par(cex.axis=0.7)
        else
          par(cex.axis=0.8)
        boxplot(Plot.This, names=rep("",length(columnName)), ylab = y.lab, cex.axis=0.8, ylim=y.lim, col=use.col)
        title(main=figTitle, cex.main=0.9, col.main="blue", font.main=4)
        axis(1,at=1:length(columnName),labels=columnName)
        abline(lheight, 0, col="blue", lty=2)
        dev.off()
      }
    }
  } 
}

################################
##  naZeroWeights function    ##
################################

# In order to fit a boxplot all values that were weighted zero should be removed.
# This is not done in a standard boxplot output window, so we need to remove
# the zero-weighted values.

naZeroWeights <- function(x, weights=NULL) {
  if(class(x) != "matrix")
    if(is.null(x$weights))
      stop("for RGList or MAList objects a weights field must be present")
  for(i in 1:ncol(x)) {
    switch(class(x),
      RGList = {
        zeros <- which(x[,i]$weights == 0)
        x$R[zeros,i] <- NA
        x$G[zeros,i] <- NA
        x$Rb[zeros,i] <- NA
        x$Gb[zeros,i] <- NA
      },
      EListRaw = {
        zeros <- which(x[,i]$weights == 0)
        x$E[zeros,i] <- NA
        x$Eb[zeros,i] <- NA
      },
      MAList = {
        zeros <- which(x[,i]$weights == 0)
        x$M[zeros,i] <- NA
        x$A[zeros,i] <- NA
        if(!is.null(x$other$EST)) x$other$EST[zeros,i] <- NA
      },
      matrix = {
        if(is.null(weights)) stop("for matrix arguments - one-channel between normalizations - the weights parameter must be provided")
        zeros <- which(weights[,i] == 0)
        x[zeros,i] <- NA
      }, stop("Given data object is not of correct class")
    )
  }
  return(x)
}

######################################
##  boxplotOverview2color function  ##
######################################

boxplotOverview2color <- function (class1=NULL, class2=NULL, class3=NULL, class4=NULL, fileName=NULL, figTitles=NULL, y.axis="M", use.col=NULL, non.zero.weight=FALSE) {
  select <- c(!is.null(class1),!is.null(class2),!is.null(class3),!is.null(class4))
  if(sum(select) < 2) stop(cat("Minimum 2 classes need to be used up until four classes maximum"))

  if(sum(select) != length(figTitles)) stop("The same number of data objects and figure titles must be specified")
  
  objects <- c("class1", "class2", "class3", "class4")[select]
  
  classes <- NULL
  for (i in  objects) classes <- c(classes, class(get(i)))
  if(sum(classes!="MAList") > 0) stop("All arguments must be of type MAList")
  
  datatypes <- NULL
  for (i in objects) datatypes <- c(datatypes, get(i)$datatype)
  if(sum(datatypes!="both") > 0) stop("All arguments must be two-channel MALists")
  
  if(non.zero.weight)
    for(i in objects) assign(i, naZeroWeights(get(i)))
  
  png(file=paste(fileName, ".png", sep=""), width=1600, height=1200, pointsize=20)
  par(mfrow=c(2,1+(sum(select) > 2)))
  
  for(i in 1:length(objects)) {
    object <- get(objects[i])
    if(y.axis == "M") {
#     boxplot(object$M~col(as.matrix(object$M)), main=figTitles[i], ylim=c(-10,10), ylab="M", col=use.col)
      data2plot <- as.data.frame(object$M)
      #names will not fit on these combined plots, so remove them
      names(data2plot) <- NULL
      boxplot(data2plot, main=figTitles[i], ylim=c(-10,10), ylab="M", col=use.col)
      abline(0, 0, col="blue", lty=2)
    }
    if(y.axis == "A") {
#     boxplot(object$A~col(as.matrix(object$A)), main=figTitles[i], ylim=c(-10,25), ylab="A", col=use.col)
      data2plot <- as.data.frame(object$A)
      #names will not fit on these combined plots, so remove them
      names(data2plot) <- NULL
      boxplot(data2plot, main=figTitles[i], ylim=c(-10,25), ylab="A", col=use.col)
      abline(mean(object$A, na.rm=TRUE), 0, col="blue", lty=2)
    }
  }
  dev.off()
}

######################################
##  boxplotOverview1color function  ##
######################################

boxplotOverview1color <- function (class1=NULL, class2=NULL, class3=NULL, class4=NULL, fileName=NULL, figTitles=NULL, use.col=NULL, non.zero.weight=FALSE, weights=NULL) {
  select <- c(!is.null(class1),!is.null(class2),!is.null(class3),!is.null(class4))
  if(sum(select) < 2) stop(cat("Minimum 2 classes need to be used up until four classes maximum"))

  if(sum(select) != length(figTitles)) stop("The same number of data objects and figure titles must be specified")

  objects <- c("class1", "class2", "class3", "class4")[select]
  
  classes <- NULL
  for (i in  objects) classes <- c(classes, class(get(i)))
  if(sum(classes!="MAList" & classes!="matrix")>0) stop("All arguments must be of type MAList or matrix")
  
  datatypes <- NULL
  for (i in objects[classes == "MAList"]) datatypes <- c(datatypes, get(i)$datatype)
  if(sum(datatypes=="both") > 0) stop("All arguments must be one-channel MALists or matrices")

  if(non.zero.weight) { 
    if((sum(classes=="MAList") == 0) & is.null(weights)) stop("when only matrix objects are provided (e.g. between normalization objects) weights must be provided")
    if((sum(classes=="MAList") > 0) & !is.null(weights)) cat("WARNING: weights provided will only be used for matrix objects, not for MALists \n")
    if((sum(classes=="MAList") > 0) & is.null(weights)) {
      cat("WARNING: no weights provided for matrix object, these will be taken from the first MAList object with weights \n")
      for (i in objects[classes == "MAList"])
        if(is.null(weights) & !is.null(get(i)$weights))
          weights <- get(i)$weights
    }
    for(i in objects) assign(i, naZeroWeights(get(i), weights))
    # note that in this call the weights parameter will be ignored for MAList objects
  }
  
  png(file=paste(fileName, ".png", sep=""), width=1600, height=1200, pointsize=20)
  par(mfrow=c(2,1+(sum(select) > 2)))
  
  for (i in 1:length(objects)) {
    object <- get(objects[i])
    if (class(object) == "MAList") {
#     boxplot(object$other$EST~col(as.matrix(object$other$EST)), main=figTitles[i], ylim=c(-10,25), ylab="estimated intensity", col=use.col)
      data2plot <- as.data.frame(object$other$EST)
      #names will not fit on these combined plots, so remove them
      names(data2plot) <- NULL
      boxplot(data2plot, main=figTitles[i], ylim=c(-10,25), ylab="estimated intensity", col=use.col)
      abline(mean(object$other$EST, na.rm=TRUE), 0, col="blue", lty=2)
    } else {
      if (class(object) == "matrix") {
#       boxplot(object~col(as.matrix(object)), main=figTitles[i], ylim=c(-10,25), ylab="estimated intensity", col=use.col)
        data2plot <- as.data.frame(object)
        #names will not fit on these combined plots, so remove them
        names(data2plot) <- NULL
        boxplot(data2plot, main=figTitles[i], ylim=c(-10,25), ylab="estimated intensity", col=use.col)
        abline(mean(object, na.rm=TRUE), 0, col="blue", lty=2)
      } else {
        stop("all objects must be of data type 'MAList' or 'matrix'")
      }
    }
  }
  dev.off()
}

###########################################
##     END OF CREATEQCPLOTS SCRIPT       ##
###########################################

