# version 1.4.0
# - Fully compatible with first dev release

require(RColorBrewer)  ## Color Palette
require(affy)          ## plotDensity function
require(gplots)        ## heatmap.2 function
require(geneplotter)   ## densCols function
require(shape)         ## imageplot3by2Adp function
require(KernSmooth)    ## Required for barplots / MvA plots
require(Ringo)         ## Required for corPlot

maxSamples <- 65
cex.scale <-  0.125


#  legendTemp <- floor( dim(x[[2]])[2] / 65 )
#  cex.axis <- 1 - (legendTemp * 0.125)
#  png(file=paste(fileName, ".png", sep=""), width=1600+(800 * legendTemp), height=1200, pointsize=20)

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
  error <- NULL
  if(!class(data.object) %in% c("RGList","MAList","EListRaw")) { error <- c(error, "- Object is not of RGList, ElistRaw, or MAList class\n") }
  if(is.null(name)) { error <- c(error, "- name parameter not filled in\n") }
  if(!is.character(name)) { error <- c(error, "- name parameter is not a character string\n") }
  if(is.null(data.object$printer)) { error <- c(error, "- The printer field in the full object is empty, whereas this information is required for image plots") }
 
  if(!is.null(error)) {
    cat(paste("[ERROR] Please resolve the following issue(s):\n", paste(error, collapse="\n"), "\n", sep=""))
  }

  if(debug==1) {
    cat("Paste the following below to debug:\n")
    cat(paste("data.object <-", deparse(substitute(data.object)), "\nwhich.field <-", deparse(substitute(which.field)), "\n"))
    cat(paste("name <-", deparse(substitute(name)), "\nhigh <-", deparse(substitute(high)), "\nlow <-", deparse(substitute(low)), "\n"))
    cat(paste("log.transform <-", deparse(substitute(log.transform)), "\nsymm <-", deparse(substitute(symm)), "\norientation <-", deparse(substitute(orientation)), "\n"))
    stop("Happy debugging!")
  }
  cat(paste("  --", paste(rep("-", nchar(name)), collapse=""), "--\n", sep="", collapse=""))
  cat(paste("  - ", name, " -\n", sep="", collapse=""))
  cat(paste("  --", paste(rep("-", nchar(name)), collapse=""), "--\n", sep="", collapse=""))

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

HierarchCluster <- function(x, dist.method = "euclidean", clust.method = "ward", main, data.return=FALSE, image.width=NULL, image.height=1400, pointsize=25, ...) {
  a <- dist(t(x), method=dist.method)
  b <- hclust(a, method=clust.method)
  if(is.null(image.width)) {
     temp <- floor( dim(x)[2] / maxSamples )
     image.width= 1600 + (800 * temp)
  }
  png(paste("Cluster_",main,".png", sep=""), width=image.width, height=image.height, pointsize=pointsize)
  if(dim(x)[2] >= 100) par(cex=0.8, cex.axis=1.25, cex.lab=1.25, cex.main=1.5, cex.sub=1.25)
  plot(b, main=main, ...)
  dev.off()
  if(data.return==TRUE) {return(b)}
}

##############################
##      Heatmap function     ##
##############################

CreateHeatMap <- function(data, main=NULL, image.width=NULL, image.height=1400, pointsize=25) {
  no.title <- ifelse(is.null(main), 1, 0)
  list.used <- 0
  if(is.list(data)) { temp.main <- names(data); data <- data[[1]]; list.used <- 1 }
  
  if(is.null(image.width) | is.null(image.height) ) {
     temp <- floor( dim(data)[2] / maxSamples )
     image.width <- 2000 + (800 * temp)
     image.height <- 1400 + (560 * temp)
  }

  my.dist <- function(x) dist(x, method="euclidean")
  my.hclust <- function(d) hclust(d, method="ward")
  legend.keysize <- 1
  if(class(data) != "matrix") {
    if(data$datatype=="both") {
      if(no.title == 1) {
        if(list.used == 1) {
          main1 <- paste("Heatmap Logratio - ", temp.main, sep="") 
        } else {
          main1 <- paste("Heatmap Logratio - ", deparse(substitute(data)), sep="")
        }
      }
      if(no.title == 0) { 
        if(length(main)>1) { stop("main variable has more than 1 argument!\n") }; 
        main1 <- paste("Heatmap Logratio - ", main, sep="") 
      }

      png(paste(gsub(" ", "_", main1), ".png", sep=""), width=image.width, height=image.height, pointsize=pointsize)
        crp <- cor(data$M, use="complete.obs")
        heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", margins=c(20,20), symm=TRUE, density.info="density", main=main1, dendrogram="row", keysize=legend.keysize)
      dev.off()
      if(no.title == 1) { 
        if(list.used == 1) {
          main2 <- paste("Heatmap Average Intensity - ", temp.main, sep="") 
        } else { 
          main2 <- paste("Heatmap Average Intensity - ", deparse(substitute(data)), sep="")
        }
      }
      if(no.title == 0) { main2 <- paste("Heatmap Average Intensity - ", main, sep="") }

      png(paste(gsub(" ", "_", main2),".png", sep=""), width=image.width, height=image.height, pointsize=pointsize)
        crp <- cor(data$A, use="complete.obs")
        heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", margins=c(20,20), symm=TRUE, density.info="density", main=main2, dendrogram="row", keysize=legend.keysize)
      dev.off()
    } else {
      if(no.title == 1) { 
        if(list.used == 1) {  
          main2 <- paste("Heatmap Estimated Intensity - ", temp.main, sep="") 
        } else { 
          main2 <- paste("Heatmap Estimated Intensity - ", deparse(substitute(data)), sep="")
        }
      }
      if(no.title == 0) { main2 <- paste("Heatmap Estimated Intensity - ", main, sep="") }

      png(paste(gsub(" ","_", main2),".png", sep=""), width=image.width, height=image.height, pointsize=pointsize)
      crp <- cor(data$other$EST, use="complete.obs")
      heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", margins=c(20,20), symm=TRUE, density.info="density", main=main2, dendrogram="row", keysize=legend.keysize)
      dev.off()
    }
  } else { # object has not been generated within arrayQC
    if(no.title == 1) { main2 <- paste("Heatmap Estimated Intensity - ", deparse(substitute(data)), sep="")  }
    if(no.title == 0) { main2 <- paste("Heatmap Estimated Intensity - ", main, sep="") }
    png(paste(gsub(" ","_", main2), ".png", sep=""), width=2000, height=1400, pointsize=18)
    crp <- cor(data, use="complete.obs")
    heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", margins=c(20,20), symm=TRUE, density.info="density", main=main2, dendrogram="row", keysize=legend.keysize)
    dev.off()
  }
}

## CreateCorplot
CreateCorplot <- function(x, which.channel=NULL, data.type=NULL, fileName = NULL) {
  y <- NULL
  if(class(x) == "list") {
    if(is.null(data.type)) {
      data.type <- names(x)
    }
    x <- x[[1]]
  }

  checks <- c("MAList", "RGList", "EListRaw", "matrix")
  if(sum(class(x) %in% checks) == 1) {
    if(sum(class(x) %in% checks[1:3]) == 1) {
      ## Now we expect that which.channel is filled in!!
      if(is.null(which.channel)) { stop("For an object of the MAList/RGList/EListRaw class the which.channel variable needs to be filled in!") }
      if(is.null(fileName)) {
        if(!is.null(data.type)) { 
          fileName <- paste("Correlation_Plot_", data.type, "_", which.channel, ".png", sep="")
        } else {
          cat("Please provide the type.data parameter if you would like to add a specific title to your filename!\n")
          fileName <- paste("Correlation_Plot_", which.channel,".png")
        }
      }
      y <- x[[which.channel]]
    }
    if( class(x) == "matrix" ) {
      y <- x
      if(is.null(fileName)) {
        if(is.null(data.type)) { 
          fileName <- "Correlation_Plot.png"
        } else { 
          fileName <- paste("Correlation_Plot_", data.type, ".png", sep="")
        }
      }
    }
  }
  if(is.null(y)) { stop(paste("object x is not of the following class:\n- ", paste(checks, collapse=" / "), sep="")) }
  png(file=fileName, width= (1400 * ceiling( dim(y)[2] / 15) ), height=(1400 * ceiling( dim(y)[2] / 15) ), pointsize=25)
  corPlot(y, useSmoothScatter=FALSE)
  dev.off()
}


##############################
## CreatePCAplot function   ##
##############################

CreatePCAplot <- function(data, main=NULL, scaled_pca=TRUE, namesInPlot=FALSE){
  # PCA performed on reporters NOT containing NAs only
  # Scaled PCA by default
  if(class(data) == "list") { 
    if(length(data) > 1) { stop("List with more than 1 element found!\n") }
    main.temp <- names(data)
    data <- data[[1]]
  } else {
    main.temp <- main
  }
  
  if(class(data)!= "matrix") {
    if(data$datatype == "both") {
      pcaData <- data$M
    } else {
      pcaData <- data$other$EST
    }
  } else {
    pcaData <- data
  }

  if(is.null(main.temp)) { baseName <- deparse(substitute(data)) } else { baseName <- main.temp }

  pca1 <- NULL  
  try(pca1 <- prcomp(t(pcaData[apply(pcaData,1,function(r) {sum(is.na(r)) == 0}),]), retx=T, center=T, scale=scaled_pca),TRUE)
  if(!is.null(pca1)) {
    perc_expl1 <- round(((pca1$sdev[1:3]^2)/sum(pca1$sdev^2))*100,2)
    
    tmain <- paste("PCA analysis of", baseName)

    plotColors <- rainbow(dim(pcaData)[2])
    arrayNames <- as.character(colnames(pcaData))
    cex.circle <- 1.3
    cex.text <- 0.6
    tcol <- "#444444"

    png(file = paste("PCAplot_",baseName,".png",sep=""), width=1200+(!namesInPlot)*400, height=1200, pointsize=25)
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
    cat("Warning: pca on the",baseName,"data set unsuccessful, possibly due to lack of memory\n\n")
  }
}
#################################
## CreateDensityPlots function ##
#################################

CreateDensityPlots <- function(x, name=NULL) {
  if(is.list(x) && !class(x) %in% c("RGList", "EListRaw", "MAList")) {
    name <- names(x)
    x <- x[[1]]
    # If maximum value of x > 1.000 then the data needs to be log transformed.
  }
  if((!class(x) %in% c("RGList","EListRaw","MAList")) & (!is.matrix(x))) stop("Only object of class RGList, EListRaw, MAList, or matrix can be handled")
  if(is.null(name)) stop("The \"name\" parameter must be provided\n")
  
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
    
    if(npages > 1) {
      outputname <- paste("DensityPlots_", name, "_Part_", i, ".png", sep="")
    } else {
      outputname <- paste("Density_Plots_", name, ".png")
    }

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

  png(file=outputName, width=1200, height=600+600*(is.null(corr.foreground)), pointsize=20)
  
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
    layout( t(as.matrix(c(1,1,2)) ) )
    if(!is.null(corr.foreground)) {
      x.min.fg <- 0
      x.max.fg <- ceiling(max(sapply(apply(corr.foreground, 2, density, na.rm=TRUE), function(z) max(z$x))))
      y.min.fg <- 0
      y.max.fg <- ceiling(max(sapply(apply(corr.foreground, 2, density, na.rm=TRUE), function(z) max(z$y))))
     if(class(x) == "list") { lab.title <- paste(names(x), " Foreground Signal Distribution") } else { lab.title <- "BG Corrected Foreground Signal Distribution" }
      plotDensity(as.matrix(corr.foreground[,selection]), col=color, lty=1:length(selection), lwd=2, xlim=c(x.min.fg,x.max.fg), ylim=c(y.min.fg,y.max.fg), main=lab.title, xlab=paste("Log2(Intensity) -",name), ylab="Density")
    }
  }
  
  plot.new()
  #if(is.null(corr.foreground)) plot.new()
  if((class(x) == "RGList" || class(x) == "EListRaw") & is.null(background2)) plot.new()

  legend("topright",label,col=color, lty=1:length(selection),text.col="black",cex=1-0.2*(!is.null(corr.foreground)), lwd=3)
  dev.off()
}

################################
##  naZeroWeights function    ##
################################

# In order to fit a boxplot all values that were weighted zero should be removed.
# This is not done in a standard boxplot output window, so we need to remove
# the zero-weighted values.

naZeroWeights <- function(x, weights=NULL) {
  if(class(x) != "matrix")
    if(is.null(weights))
      stop("Please provide a weights matrix")
    if(dim(x)[2] != dim(weights)[2] | dim(x)[1] != dim(weights)[1] )
      stop("Object dimension and weight matrix dimensions do not match!")

  for(i in 1:ncol(x)) {
    switch(class(x),
      RGList = {
        zeros <- which(weights[,i] == 0)
        x$R[zeros,i] <- NA
        x$G[zeros,i] <- NA
        x$Rb[zeros,i] <- NA
        x$Gb[zeros,i] <- NA
      },
      EListRaw = {
        zeros <- which(weights[,i] == 0)
        x$E[zeros,i] <- NA
        x$Eb[zeros,i] <- NA
      },
      MAList = {
        zeros <- which(weights[,i] == 0)
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
##  boxplotOverview function  ##
######################################

boxplotOverview <- function(x, fileName=NULL, figTitles=NULL, groupcols=NULL, use.weights=FALSE, weights = NULL, max.characters=20, y.axis=NULL, y.lim=NULL) {
  ## This function assumes that x is a list coming from arrayQC consisting of four or five data matrices / MALists.
  error <- NULL
  if(use.weights) { tempTitle <- "[HIGH QUALITY SPOTS ONLY]" } else { tempTitle <- "[ALL]" }
  if(is.null(fileName)) { fileName <- "BoxplotOverview" }  ## Check fileName
  
  ## Check if y-axis is filled in properly (in combination with y.lim)
  check <- 0
  if(is.null(y.axis)) { 
    check <- 1 
  } else {
    switch(y.axis, "M" = { 
     check <- 1
      fileName <- paste(fileName, "_M", sep="")
      if(is.null(y.lim)) { y.lim=c(-6, 6) }
    }, "A" = { 
      check <- 1
      fileName <- paste(fileName, "_AverageSignal", sep="") 
    })
  }
  if(check == 0) { error <- c(error, "- y.axis: parameter is not set to NULL, \"M\" or \"A\"") }
  rm(check)

  if(use.weights == TRUE) { fileName <- paste(fileName, "_weighted", sep="") }
  if(is.null(y.lim)) { y.lim <- c(0,20) }

  if(!is.null(figTitles)) { 
    if( length(figTitles) != length(x) ) { 
      error <- c(error, "- figTitles: data object size does not correspond with figure title length") 
    } else {
      figTitles <- paste( figTitles, rep(tempTitle, length(x)) )
      figTitles <- gsub("."," + ", figTitles, fixed=TRUE)
    }
  } else {
    ## Creating Figure Titles if none are supplied
    figTitles <- paste( names(x), rep(tempTitle, length(x)) )
    figTitles <- gsub("."," + ", figTitles, fixed=TRUE)
  }
  if(is.list(x)) {  
    if(length(x) == 5) { ## Happens for two-color arrays (5 normalizations)
      normalizations <- c("BGCORRECTED", "LOESS", "LOESS.QUANTILE", "LOESS.AQUANTILE") ## Normalization values
      ## Check if 4 normalizations are present:
      if( sum( normalizations %in% names(x) ) != 4 ) { 
        error <- c(error, "- x: object is a list with 5 elements, but does not contain all recommended normalization results") 
      } else { ## if everything is correct, make the proper data object:
        selected <- which( names(x) %in% normalizations == FALSE )
        x[[selected]] <- NULL
      }
    } else { ## Single-channel arrays
      if(length(x) != 4) { error <- c(error, "- x: object is a list with more or less than 4 data fields") }
      temp <- as.vector(lapply(x, class))
      if(sum(temp!="matrix") > 0) { error <- c(error, "- All elements of the list must be of the matrix class.") }
      rm(temp)
    }
  } else {
    error <- c(error, "- x: object is not a list") 
  }

  ## Check the maximum characters of the values assigned to the legend:
  a <- nchar(colnames(x[[1]]))
  b <- a[] >= max.characters
  if(sum(b, na.rm=TRUE) > 0) {
    error <- c(error, paste("- The following sample descriptions are too long ( >= ", max.characters, " characters) :\n", sep=""))
    error <- c(error, paste("-", a[b], "\n"))
  }


  if(!is.null(error)) { stop(paste("[WARNING] The following error(s) occurred:\n", paste(error, collapse="\n"), "\n\nPlease adress the above issues!\n")) }

  ## Image legend will support up to 65 sample names per row. If more samples are present, the image needs to widen up.
  legendTemp <- floor( dim(x[[2]])[2] / 65 )
  cex.axis <- 1 - (legendTemp * 0.125)
  png(file=paste(fileName, ".png", sep=""), width=1600+(800 * legendTemp), height=1200, pointsize=20)
    image.layout <- rbind( c(1,1,2,2,5), c(3,3,4,4,5) )
    layout(image.layout)
    for(i in 1:length(x)) {
      ## Which value to plot
      if( use.weights == 1 ) { 
        if(is.null(weights)) { weights <- x[[i]]$weights }
        data <- naZeroWeights( x[[i]], weights=weights ) 
      } else { data <- x[[i]] }

      if( is.null(y.axis) ) { data <- data   } else { 
        if( y.axis == "M" ) { data <- data$M }
        if( y.axis == "A" ) { data <- data$A }
      }
      boxplot(data, names=c(1:dim(x[[i]])[2]), ylim=y.lim, main=figTitles[i], cex.axis=cex.axis, las=2)
    }
    plot(0,type='n',xaxt='n',yaxt='n',xlab="",ylab="", bty='n')
    legend("topright", paste( c(1:dim(x[[1]])[2]), ": ", colnames(x[[1]]), sep=""), ncol=(legendTemp + 1), box.lwd = 0,box.col = "white",bg = "white")
  dev.off()

}


##############################
## CreateMAplots function   ##
##############################
createMAplots <- function(x, lab=NULL, weight = NULL, postfix=NULL, image.width=1600, image.height=1200, pointsize=25, y.lim=NULL, x.lim=NULL) {
  cat("* Preprocessing MA Plots ...")
  error <- NULL
  if( class(x) != "list" ) { error <- c(error, "- object x is not a list class!") }
  ## Check if each component of x is an MAlist
  for(i in 1:length(x)) {
    if(class(x[[i]]) != "MAList") { error <- c(error, paste("- Element ", names(x)[i], " is not of MAList class!", sep="")) }
  }
  ## Only allow lists with 2 or 5 elements (single channel vs dual channel only).
  if(! (length(x) == 2 | length(x) == 5) ) { error <- c(error, "- lists must consist of 2 or 5 elements!") }
  ## If 5 elements are supplied, extract only the more usefull normalizations:
  if(length(x) == 5) {
    normalizations <- c("BGCORRECTED", "LOESS", "LOESS.QUANTILE", "LOESS.AQUANTILE") ## Normalization values
    ## Check if 4 normalizations are present:
    if( sum( normalizations %in% names(x) ) != 4 ) { 
      error <- c(error, "- x: object is a list with 5 elements, but does not contain all recommended normalization results") 
    } else { ## if everything is correct, make the proper data object:
      selected <- which( names(x) %in% normalizations == FALSE )
      x[[selected]] <- NULL
    }
  }

  ## Check if labels exist
  if(is.null(lab)) { lab <- names(x) }
  if(sum( lab %in% names(x) ) != length(lab) ) { error <- c(error, "- label names do not match names of object x!") }

  
  subtext1 <- paste( lab, " - High quality spots (filtered)", sep="")
  subtext2 <- paste( lab, " - All spots (not filtered)", sep="")
  
  ## Determine x and y limits (based on full experiment) and check whether objects have a weight field

  tempM <- tempA <- NULL
  for(i in 1:length(x)) {
    tempM <- c( tempM, range(x[[i]]$M, na.rm=TRUE) )
    tempA <- c( tempA, range(x[[i]]$A, na.rm=TRUE) )
  }
  if(is.null(x.lim)) { x.lim <- c( min( floor(tempA) ) - 1, max( ceiling(tempA) ) + 1 ) }
  if(is.null(y.lim)) { y.lim <- c( min( floor(tempM) ) - 1, max( ceiling(tempM) ) + 1 ) }

  cat("  ok.\n* MA plotting progress:\n")
  #### PLOTTING
  #############
  maxArrays <- ncol(x[[1]])
  pbar <- txtProgressBar(min=0, max=maxArrays, char="*", width=20, style=3)
  for(i in 1:maxArrays) {
    ### WEIGHTED PLOT
    if(!is.null(weights)) {
      titleName <- colnames(x[[1]])[i]
      fileName <- paste("MA_Plots_-_", titleName, "_-_Zero_Weighted", postfix, ".png", sep="", collapse=NULL)
      png(filename = fileName, width=image.width, height=image.height, pointsize=pointsize)
      par( mar=c(8,4,4,3)+0.1, mfrow=c(1+ (length(x) == 2 | length(x) == 4), 1+ ( length(x) == 4)) )
      for(j in 1:length(x)) {
        m <- x[[j]]$M[,i]
        a <- x[[j]]$A[,i]
        w <- weight[,i]
        k <- is.na(w) | (w <= 0) 
        m[k] <- NA
        a[k] <- NA
        #Determining color densities and plotting the MA-plot
        color <- densCols(a, m)
        plot(a, m, xlim=x.lim, ylim=y.lim, col=color, cex=0.3, pch=16, xlab="A", ylab="M")
        title(main=titleName, sub=subtext1[j], outer=FALSE)
        abline(0,0,col="darkgrey", lty=2)
      }
      dev.off()
    }
    ### NORMAL PLOT
    titleName <- colnames(x[[1]])[i]
    fileName <- paste("MA_Plots_-_", titleName,"_-_All_Points",postfix,".png", sep="", collapse=NULL)
    png(filename = fileName, width=image.width, height=image.height, pointsize=pointsize)
    par( mar=c(8,4,4,3)+0.1, mfrow=c(1+ (length(x) == 2 | length(x) == 4), 1+ ( length(x) == 4)) )
    for(j in 1:length(x)) {
      m <- x[[j]]$M[,i]
      a <- x[[j]]$A[,i]
      #Determining color densities and plotting the MA-plot
      color <- densCols(a, m)
      plot(a, m, xlim=x.lim, ylim=y.lim, col=color, cex=0.3, pch=16, xlab="A", ylab="M")
      title(main=titleName, sub=subtext2[j], outer=FALSE)
      abline(0,0,col="darkgrey", lty=2)
    }
    dev.off()
    setTxtProgressBar(pbar, i)
    if(i == maxArrays) {cat("\n")}
  }      
}



###########################################
##     END OF CREATEQCPLOTS SCRIPT       ##
###########################################
