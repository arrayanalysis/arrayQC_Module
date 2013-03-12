## One and two channel QC - Dept. of Bioinformatics - BiGCaT - Maastricht University - the Netherlands
# Version: 1.3.1
# Last adjustment: 25-01-2013
## Note: Changed xBGMeanSignal to xBGMediansignal, to reflect the xMedianSignal parameter (it basically means this).
## need to check if any other functions are affected.

#####################
## START OF SCRIPT ##
#####################
## VALIDATING USER ARGUMENTS
workpath  <- paste(datapath, "qualityControl/", sep="")
if(!is.null(number.of.blocks)) {
  if(is.null(nblock.row) || is.null(nblock.row)) { stop("Please define your nrow.block and ncol.block in the main script!") }
  blocks <-list(total=number.of.blocks, nblock.row=nblock.row, nblock.col=nblock.col)
} else {
  blocks <- NULL
}

dataformat <- tolower(dataformat)
datatype <- tolower(datatype)

dataformat <- match.arg(dataformat, c("fes", "genepix", "generic"))
  if(dataformat == "fes") { package <- "agilent" }
  if(dataformat == "genepix") { package <- "genepix" }
  if(dataformat == "generic") { package <- "generic" } 
datatype <- match.arg(datatype, c("two-channel","green","red"))

## Check if all necessary packages are installed
source(paste(arrayQC.scriptpath, "install_arrayQC_packages.R", sep=""))

cat("* Entering data directory...")
setwd(datapath)
cat(" done.\n")

cat("* Creating output directory...")
dir.create(workpath, showWarnings = FALSE)
cat(" done.\n")


library(RColorBrewer)
library(limma)
library(affy)

# Reading specific R-scripts (use the reload()-function if you want to reload the scripts later on, e.g. after modifications)
reload <- function() {
  source(paste(arrayQC.scriptpath,"requiredColumns.R", sep=""))
  source(paste(arrayQC.scriptpath,"ReadFiles.R", sep=""))
  source(paste(arrayQC.scriptpath,"QC.R", sep=""))
  source(paste(arrayQC.scriptpath,"CreateQCPlots.R", sep=""))
  source(paste(arrayQC.scriptpath,"functions.EListRaw.R", sep=""))
}

reload()

###############
## IMPORTANT ##
###############
# Depending on the usage of the Feature Extraction Softare or GenePix Pro software the appropriate columns need to be checked
# manually to see if the columns are really present. If not, adjust the variable accordingly. A short summary of the needed
# columns follows below.

## Necessary columns to perform an analysis:
## -----------------------------------------
# - r/gMeanSignal (Mean Foreground signal)
# - r/gBGUsed (Measured or Estimated Backgroundsignal (usually the latter - See "Spatial Detrending" in Agilent manual)

## Necessary for Extensive Quality Control:
## ----------------------------------------
# - r/gNumPix (Number of pixels detected in one spot - These pixels are used to extract other data (intensities))
# - r/gIsWellAboveBG (Is my signal above 2.6x SD of the BG signal?)
# - r/gIsSaturated (Do I have more than 75% saturated pixels in my spot) 
# - r/gIsFound (Predecessor check to r/gIsWellAboveBackground, can actually be omitted)
# - r/gMedianSignal (Need this part for the mean / median ratio calculation)

## Additional Information
## ----------------------
# - LogRatio (Calculated intensity from rProcessed / gProcessed intensities - Not a necessity here)

## Annotation information
## ----------------------
# - FeatureNum (Number of the reporter (i.e. NOT the name!)
# - Row (row position of the reporter)
# - Col (col position of the reporter)
# - ProbeName (Agilent specific reporter ID)
# - Controltype (1 - Positive control; 0 - No control; -1 - Negative control)
# - GeneName (target gene of the reporter)
# - Description (Full name of the GeneName - Also includes possible crossreference
# - SystematicName (No Idea - Have to look this up)

## For "generic", a minimum number of columns is needed for arrayQC to run.
## All information is not always available to perform spot-specific QC, and as such
## the scripts should automatically detect this.

## It will be difficult to guess which columns are needed, but for basic arrayQC functionality the needed
## columns should be:
## - R/G Foreground + R/G Background [Intensities]
## - Reporter ID + Gene Annotation   [Annotation]
## - Row / Column                    [Array Layout + Block Detection]

#####
## CREATING / READING COLUMN HEADER FILE
#####
default.columnHeaderFile <- "columnHeaders.arrayQC"
if(exists("columnHeaderFile")) { 
  if(is.null(columnHeaderFile)) { columnHeaderFile <- default.columnHeaderFile }
} else {
  columnHeaderFile <- default.columnHeaderFile
}

if(.FileExists(columnHeaderFile, xpath=datapath)==0) {
  GenerateHeaderFile(xpath=datapath, fileName = columnHeaderFile, datatype = datatype, package=package)
  switch(package,
    "generic" = {
      cat(paste("\n-> A new default column header file", columnHeaderFile, " has been generated in", datapath, ".\n\n"))
      cat("WARNING!!\n   Please fill in the correct column names under the ColumnHeader column and restart arrayQC!\n")
      cat("     --> Note: if you use a different filename, then please set the columnHeaderFile parameter in the original arrayQC.R script!\n")
      cat("         Example:  columnHeaderFile <- \"new.filename.txt\"\n")     
      stop("arrayQC halted!") 
    },
    "agilent" = {
      cat("\n -> Creating default FES profile for use within arrayQC.\n")
    },
    "genepix" = {
      cat("\n -> Creating default GenePix Pro Format for use within arrayQC.\n")
    })
}

cat("\n -> Reading profile file for current arrayQC run.\n")
arrayQC.description <- parseHeaderFile(xpath = datapath, fileName = columnHeaderFile, datatype = datatype, package = package)

columns <- arrayQC.description$columns

other.columns <- arrayQC.description$other.columns
annotation <- arrayQC.description$annotation

## Pre-check to see if basic files (spottypes and description) are present.
# - If spottypes file is not present --> continue
# - If description file is not present --> stop

if(!file.exists(paste(datapath, "spottypes.txt", sep=""))) { 
  spotfile <- NULL
} else {
  spotfile <- "spottypes.txt"
}

if(!file.exists(paste(datapath, "description.txt", sep=""))) { 
  experimentalDescr <- NULL
  cat(paste("\n>>WARNING<<\n\nThe \"description.txt\" file was not located in", datapath, "and is necessary for arrayQC to proceed.\n"))
  cat("If this file has a different name, please enter the correct name below:\n  ")
  experimentalDescr <- scan(n=1, what="character")
  check <- .FileExists(experimentalDescr, xpath = datapath)
  if(check == 0 || nchar(experimentalDescr) < 2) { stop(paste("\nFile \"", datapath, experimentalDescr, "\" could not be found. Terminating arrayQC script.", sep="")) }
  if(check == 1) { cat(paste(experimentalDescr, "successfully located. Proceeding with arrayQC...\n")) }
} else {
  experimentalDescr <- "description.txt"
}

## Reading in scanner output files and save this in the RG-object (makes use of read.maimages()):

setwd(workpath)
RG <- ReadFiles(description.file=experimentalDescr, spottypes.file=spotfile, data.path=datapath, columns=columns, other=other.columns, annotation=annotation, blocks=blocks, source=package, use.description=TRUE, save.backup=TRUE, manual.flags=useAsManualFlags, controlType.value=controlType.value, arrayQC.path=arrayQC.scriptpath)

## Sometimes the software doesn't provide information about the Blocknumber, whereas this information
#  is present in the GAL file. .addBlockInfo will recreate these numbers and add them to the annotation
#  field.
if(!is.null(number.of.blocks)) {
  blocks$nspot.row <- max(RG$genes$Row) / blocks$nblock.row
  blocks$nspot.col <- max(RG$genes$Col) / blocks$nblock.col
}

## If provided, use manual flags (note that the manualFlags columns has been made binary by the ReadFiles function)
if(!is.null(RG$other$manualFlags))
  RG <- removeManualFlaggedSpots(RG)

## Perform this step if you want to update annotation information. Depending on the annotation file, you might need to change 
## fields that need to be read in.
## Note: this code only applies when the order of reporters in the file is the same as in the data object!
# upd.annotation <- read.delim("annotation.tab", as.is=TRUE)
# RG$genes$Description <- upd.annotation$Description
# RG$genes$Ensembl_ID <- upd.annotation$Ensembl_ID 

## Performing makeLimmaCompatible step
RG <- makeLimmaCompatible(RG)

## Performing QualityControl

switch(RG$source, "agilent" = {
  switch(RG$datatype,"red" =   {  RG <- QualityControl(RG, num.pix=NULL, low.pix=TRUE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-5,-5,-9,-5), criteria=c(-5,-5,-5,-5))    },
                     "green" = {  RG <- QualityControl(RG, num.pix=NULL, low.pix=TRUE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-3,-3,-9,-3), criteria=c(-3,-3,-3,-3))    }, 
                               {  RG <- QualityControl(RG, num.pix=NULL, low.pix=TRUE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-3,-3,-9,-3), criteria=c(-3,-3,-8,-8)) })
  }, "genepix" = {
  switch(RG$datatype, "red" =   { RG <- QualityControl(RG, num.pix=NULL, low.pix=FALSE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-5,-9,-5), criteria=c(-5,-5,-5))    },
                      "green" = { RG <- QualityControl(RG, num.pix=NULL, low.pix=FALSE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-3,-9,-3), criteria=c(-3,-3,-3))    },
                                { RG <- QualityControl(RG, num.pix=NULL, low.pix=FALSE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-3,-9,-3), criteria=c(-3,-8,-8))    })
  }, cat("SpotQC skipped due to non-agilent / non-genepix file format. Please run the QualityControl() function manually to include this information!\n")
)

cat("*------------\n| Generating Summary files\n*------------\n")
## AddSummary
RG <- AddSummary(RG)

## CreateSummaryPlots creates barcharts from the information in the summary object.
CreateSummaryPlots(RG)


cat("*------------\n| Data Normalization\n*------------\n")
## MA are MA-objects, MA2 are matrices (single-channel only)
MA <- MA2 <- list()

if(RG$datatype == "both") {
  MA[["BGCORRECTED"]] <- normalizeWithinArrays(RG, RG$printer, method="none", bc.method="subtract", offset=0, weights=NULL)
  MA[["LOESS"]] <- normalizeWithinArrays(RG, RG$printer, method="loess", iterations=4, bc.method="subtract", offset=0, weights=RG$other$weights.norm)
  MA[["LOESS.SCALED"]] <- normalizeBetweenArrays(MA[["LOESS"]], method="scale")
  MA[["LOESS.QUANTILE"]] <- normalizeBetweenArrays(MA[["LOESS"]], method="quantile")
  MA[["LOESS.AQUANTILE"]] <- normalizeBetweenArrays(MA[["LOESS"]], method="Aquantile")
} else {
  ## normalizeWithinArrays does NOT support EListRaw objects, as such we have to make use of RG2
  RG2 <- RG
  RG2$R <- RG2$G <- RG2$E
  RG2$Rb <- RG2$Gb <- RG2$Eb
  switch(datatype, "red" = {
    RG2$G[] <- apply(RG2$R,1,mean,na.rm=TRUE)
    RG2$Gb[] <- apply(RG2$Rb,1,mean,na.rm=TRUE)
  }, "green" = {
    RG2$R[] <- apply(RG2$G,1,mean,na.rm=TRUE)
    RG2$Rb[] <- apply(RG2$Gb,1,mean,na.rm=TRUE)
  })
  MA[["BGCORRECTED"]] <- normalizeWithinArrays(RG2, RG2$printer, method="none", bc.method="subtract", offset=0, weights=NULL)
  switch(datatype, "red" = {   MA[["BGCORRECTED"]]$other$EST <- MA[["BGCORRECTED"]]$A + MA[["BGCORRECTED"]]$M/2 },
                   "green" = { MA[["BGCORRECTED"]]$other$EST <- MA[["BGCORRECTED"]]$A - MA[["BGCORRECTED"]]$M/2 })
  MA[["LOESS"]] <- normalizeWithinArrays(RG2, RG2$printer, method="loess", iterations=4, bc.method="subtract", offset=0, weights=RG2$other$weights.norm) 
  switch(datatype, "red" = {   MA[["LOESS"]]$other$EST <- MA[["LOESS"]]$A + MA[["LOESS"]]$M/2 },
                   "green" = { MA[["LOESS"]]$other$EST <- MA[["LOESS"]]$A - MA[["LOESS"]]$M/2 })

  MA2[["BGCORRECTED"]] <- MA[["BGCORRECTED"]]$other$EST
  MA2[["LOESS"]] <- MA[["BGCORRECTED"]]$other$EST
  MA2[["SCALED"]] <- normalizeBetweenArrays(MA[["BGCORRECTED"]]$other$EST, method="scale")
  MA2[["QUANTILE"]] <- normalizeBetweenArrays(MA[["BGCORRECTED"]]$other$EST, method="quantile")
  rm(RG2)
}
cat("status: OK\n\n")

if(plotVirtualImages == 1) {
  ## Plotting of virtual array images
  cat("*------------\n| Virtual Images\n*------------\n")
  switch(RG$datatype, "both" = {
    cat("* Based on red channel ...\n")
    imageplot3by2Adp(RG, RG$R, "RedSignal", high="red", low="yellow")
    imageplot3by2Adp(RG, RG$Rb, "RedBckgrndSignal", high="red", low="yellow")
    cat("* Based on green channel ...\n")
    imageplot3by2Adp(RG, RG$G, "GreenSignal", high="#00aa00", low="yellow")
    imageplot3by2Adp(RG, RG$Gb, "GreenBckgrndSignal", high="#00aa00", low="yellow")
  }, "red" = {
    cat("* Based on red channel ...\n")
    imageplot3by2Adp(RG, RG$E, "RedSignal", high="red", low="yellow")
    imageplot3by2Adp(RG, RG$Eb, "RedBckgrndSignal", high="red", low="yellow")
  }, "green" = {
    cat("* Based on green channel ...\n")
    imageplot3by2Adp(RG, RG$E, "GreenSignal", high="#00aa00", low="yellow")
    imageplot3by2Adp(RG, RG$Eb, "GreenBckgrndSignal", high="#00aa00", low="yellow")
  })

  cat("* Based on normalized data ...\n") 
  if(RG$datatype == "both") {
    for(i in 1:length(MA)) {
      imageplot3by2Adp(MA[[i]], MA[[i]]$M, paste(names(MA)[i],"_LogRatio", sep=""), high="blue", low="yellow", symm=TRUE)
      imageplot3by2Adp(MA[[i]], MA[[i]]$A, paste(names(MA)[i],"_AverageIntensity", sep=""), high="blue", low="yellow", symm=TRUE)
    }
  } else {
    temp <- names(MA2)
    for(i in 1:length(MA2)) {
      if(temp[i] == "LOESS") {
        imageplot3by2Adp(MA[["LOESS"]], MA2[[i]], paste(temp[i],"_EstimatedSignal", sep=""), high="blue", low="yellow", symm=TRUE)
      } else {
        imageplot3by2Adp(MA[["BGCORRECTED"]], MA2[[i]], paste(temp[i],"_EstimatedSignal", sep=""), high="blue", low="yellow", symm=TRUE)
      }
    }
    rm(temp)
  }
  cat("* Based on various (quality) parameters ...\n") 
  if(!is.null(RG$other$rNumPix)) { imageplot3by2Adp(RG, RG$other$rNumPix, "RedNumberPixels", high="red", low="yellow") }
  if(!is.null(RG$other$rIsWellAboveBG)) { imageplot3by2Adp(RG, RG$other$rIsWellAboveBG, "RedWellAboveBckgrnd", high="red", low="white") }
  if(!is.null(RG$other$rIsSaturated)) { imageplot3by2Adp(RG, RG$other$rIsSaturated, "RedSaturated", high="red", low="white") }
  if(!is.null(RG$other$gNumPix)) { imageplot3by2Adp(RG, RG$other$gNumPix, "GreenNumberPixels", high="#00aa00", low="yellow") }
  if(!is.null(RG$other$gIsWellAboveBG)) { imageplot3by2Adp(RG, RG$other$gIsWellAboveBG, "GreenWellAboveBckgrnd", high="green", low="white") }
  if(!is.null(RG$other$gIsSaturated)) { imageplot3by2Adp(RG, RG$other$gIsSaturated, "GreenSaturated", high="green", low="white") }
  if(!is.null(RG$other$manualFlags)) imageplot3by2Adp(RG, RG$other$manualFlags, "ManualFlags", high="sienna", low="orange")
  if(!is.null(RG$other$LowNumPix)) imageplot3by2Adp(RG, RG$other$LowNumPix, "LowNumPix", high="orange", low="sienna")
  if(!is.null(RG$other$LowMeanMedianRatio)) imageplot3by2Adp(RG, RG$other$LowMeanMedianRatio, "LowMeanMedianRatio", high="orange", low="sienna")
  if(!is.null(RG$other$NOTWellAboveBG)) imageplot3by2Adp(RG, RG$other$NOTWellAboveBG, "NOTWellAboveBG", high="orange", low="sienna")
  if(!is.null(RG$other$Saturated)) imageplot3by2Adp(RG, RG$other$Saturated, "Saturated", high="orange", low="sienna")
  if(!is.null(RG$other$weights.norm)) imageplot3by2Adp(RG, RG$other$weights.norm, "WeightsForNormalization", high="purple2", low="white")
  if(!is.null(RG$weights)) imageplot3by2Adp(RG, RG$weights, "arrayQC_Weights", high="purple2", low="white")
}
cat(" done.\n\nstatus: OK\n\n")

if( dim(RG)[2] < 3 ) {
  cat("[[ CLUSTER ]] No clustering was performed due to too few samples!\n")
  plotClust <- 0
}

if(plotClust == 1) {
  if(!exists(as.character(substitute(cluster.distance)))) { cluster.distance <- "euclidean" }
  if(!exists(as.character(substitute(cluster.method)))) { cluster.method <- "ward" }
  cat("*------------\n| Clustering\n*------------\n")
  ## Cluster plots - by default: cluster.distance="euclidean", cluster.method="ward". Can be changed by other methods available in ?hdist and ?hclust
  if(RG$datatype=="both") {
    for(i in 1:length(MA)) {
      switch(names(MA)[i], "RAW" = { ######### NEED TO CHECK IF THIS IS STILL VALID!!!
         addTitle <- "RAW Data"
      }, addTitle <- paste( names(MA)[i], "Normalized Data"))
      HierarchCluster(MA[[i]]$M, main=addTitle, dist.method = cluster.distance, clust.method = cluster.method)
      rm(addTitle)
    }
  } else {
    for(i in 1:length(MA2)) {
      switch(names(MA2)[i], "BGCORRECTED" = {
        addTitle <- "BGCorrected Data)"
      }, addTitle <- paste(names(MA2)[i], "Normalized Data") })
      HierarchCluster(MA2[[i]], main=addTitle, dist.method = cluster.distance, clust.method = cluster.method)
      rm(addTitle)
    }
  }
  cat("status: OK\n\n")
}

if(plotCor == 1) {
  cat("*------------\n| Correlation plots \n*------------\n")
  switch(datatype, "both" = {
    CreateCorplot(RG, which.channel="R", data.type="RAW_Red_data")
    CreateCorplot(RG, which.channel="G", data.type="RAW_Green_data")
  }, "green" = {
    CreateCorplot(RG, which.channel="E", data.type="RAW_Green_data") 
  }, "red" = {
    CreateCorplot(RG, which.channel="E", data.type="RAW_Red_data")
  })
  if(datatype == "both") {
    for(i in 1:length(MA)) {
      CreateCorplot(MA[i], which.channel="M")
      CreateCorplot(MA[i], which.channel="A")
      cat(".")
    }
  } else {
    for(i in 1:length(MA2)) {
      CreateCorplot(MA2[i])
      cat(".")
    }
  }
  cat("\nstatus: OK\n\n")
}

## Heatmaps
cat("*------------\n| Heatmaps \n*------------\n")
if(plotHeatmap == 1) {
  for(i in 1:length(MA)) {
    CreateHeatMap(MA[i])
  }
}


## PCA plots
if(plotPCA == 1) {
  cat("*------------\n| PCA plots\n*------------\n")
  switch(datatype, "two-channel" = {
    for(i in 1:length(MA)) {
      CreatePCAplot(MA[i])
    }
  }, {
    for(i in 1:length(MA2)) {
      CreatePCAplot(MA2[i])
    }
  })
  cat("status: OK\n\n")
}


## Density Plots
if(plotDensity == 1) {
  cat("*------------\n| Density Plots\n*------------\n")
#  CreateDensityPlots(RG, name="Raw Data") ## We seem to have some issues with EListRaw objects. RAW data doesn't provide extra information that BGCORRECTED doesn't show (I think), so adding a comment line (for now)
  switch(datatype, "two-channel" = {
    for(i in 1:length(MA)) { CreateDensityPlots(MA[i]) }
  }, {
    for(i in 1:length(MA2)) { CreateDensityPlots(MA2[i]) }
  })
}

if(plotBoxplot == 1) {
  cat("*------------\n| Boxplots\n*------------\n")
  boxplotOverview
}




## Boxplots
cat("*------------\n| Boxplots\n*------------\n")
if(RG$datatype=="both") {
  CreateBoxplot(MA.RAW.W, "Boxplots_Before_Normalization", "Unweighted BoxPlots Before Normalization (ALL)", non.zero.weight=FALSE)
  CreateBoxplot(MA.LOESS.W, "Boxplots_LOESS_Normalization", "Unweighted BoxPlots After LOESS Normalization (ALL)", non.zero.weight=FALSE)
  CreateBoxplot(MA.LOESS.SCALED.W, "Boxplots_LOESS_Normalization_And_Scaled", "Unweighted BoxPlots After LOESS + SCALED (ALL)", non.zero.weight=FALSE)
  CreateBoxplot(MA.LOESS.QUANTILE.W, "Boxplots_LOESS_Normalization_And_Quantile", "Unweighted BoxPlots After LOESS + QUANTILE (ALL)", non.zero.weight=FALSE)
  CreateBoxplot(MA.LOESS.AQUANTILE.W, "Boxplots_LOESS_Normalization_And_AQuantile", "Unweighted BoxPlots After LOESS + AQUANTILE (ALL)", non.zero.weight=FALSE)
  if(!is.null(MA.RAW.W$weights))
    CreateBoxplot(MA.RAW.W, "Boxplots_Before_Normalization", "Weighted BoxPlots Before Normalization (FILTERED)", non.zero.weight=TRUE)
  if(!is.null(MA.LOESS.W$weights))
    CreateBoxplot(MA.LOESS.W, "Boxplots_LOESS_Normalization", "Weighted BoxPlots After LOESS Normalization (FILTERED)", non.zero.weight=TRUE)
  if(!is.null(MA.LOESS.SCALED.W$weights))
    CreateBoxplot(MA.LOESS.SCALED.W, "Boxplots_LOESS_Normalization_And_Scaled", "Weighted BoxPlots After LOESS + SCALED (FILTERED)", non.zero.weight=TRUE)
  if(!is.null(MA.LOESS.QUANTILE.W$weights))
    CreateBoxplot(MA.LOESS.QUANTILE.W, "Boxplots_LOESS_Normalization_And_Quantile", "Weighted BoxPlots After LOESS + QUANTILE (FILTERED)", non.zero.weight=TRUE)
  if(!is.null(MA.LOESS.AQUANTILE.W$weights))
    CreateBoxplot(MA.LOESS.AQUANTILE.W, "Boxplots_LOESS_Normalization_And_AQuantile", "Weighted BoxPlots After LOESS + AQUANTILE (FILTERED)", non.zero.weight=TRUE)
} else {
  CreateBoxplot(MA.RAW.W, "Boxplots_Before_Normalization", "Unweighted BoxPlots Before Normalization (ALL)", non.zero.weight=FALSE)
  CreateBoxplot(MA.LOESS.W, "Boxplots_LOESS_Normalization", "Unweighted BoxPlots After LOESS Normalization (ALL)", non.zero.weight=FALSE)
  CreateBoxplot(MA.SCALED, "Boxplots_Scaled_Normalization", "Unweighted BoxPlots After SCALED (ALL)", non.zero.weight=FALSE)
  CreateBoxplot(MA.QUANTILE, "Boxplots_Quantile_Normalization", "Unweighted BoxPlots After QUANTILE (ALL)", non.zero.weight=FALSE)
  if(!is.null(MA.RAW.W$weights))
    CreateBoxplot(MA.RAW.W, "Boxplots_Before_Normalization", "Weighted BoxPlots Before Normalization (FILTERED)", non.zero.weight=TRUE)
  if(!is.null(MA.LOESS.W$weights))
    CreateBoxplot(MA.LOESS.W, "Boxplots_LOESS_Normalization", "Weighted BoxPlots After LOESS Normalization (FILTERED)", non.zero.weight=TRUE)
  if(!is.null(MA.RAW.W$weights))
    CreateBoxplot(MA.SCALED, "Boxplots_Scaled_Normalization", "Weighted BoxPlots After SCALED (FILTERED)", non.zero.weight=TRUE, weights=MA.RAW.W$weights)
  if(!is.null(MA.RAW.W$weights))
    CreateBoxplot(MA.QUANTILE, "Boxplots_Quantile_Normalization", "Weighted BoxPlots After QUANTILE (FILTERED)", non.zero.weight=TRUE, weights=MA.RAW.W$weights)
}

## BoxplotOverview.
if(RG$datatype=="both") {
  boxplotOverview2color(class1=MA.RAW.W, class2=MA.LOESS.W, class3=MA.LOESS.SCALED.W, class4=MA.LOESS.QUANTILE.W, figTitles=c("Unnormalized data (ALL)","LOESS Normalization (ALL)","LOESS + SCALED (ALL)","LOESS + QUANTILE (ALL)"), non.zero.weight=FALSE, fileName="Boxplots_Overview_MValue_ALL_Reporters", y.axis="M")
  boxplotOverview2color(class1=MA.RAW.W, class2=MA.LOESS.W, class3=MA.LOESS.SCALED.W, class4=MA.LOESS.QUANTILE.W, figTitles=c("Unnormalized data (FILTERED)","LOESS Normalization (FILTERED)","LOESS + SCALED (FILTERED)","LOESS + AQUANTILE (FILTERED)"), non.zero.weight=TRUE, fileName="Boxplots_Overview_MValue_FILTERED_Reporters", y.axis="M")
  boxplotOverview2color(class1=MA.RAW.W, class2=MA.LOESS.W, class3=MA.LOESS.SCALED.W, class4=MA.LOESS.QUANTILE.W, figTitles=c("Unnormalized data (ALL)","LOESS Normalization (ALL)","LOESS + SCALED (ALL)","LOESS + QUANTILE (ALL)"), non.zero.weight=FALSE, fileName="Boxplots_Overview_AValue_ALL_Reporters", y.axis="A")
  boxplotOverview2color(class1=MA.RAW.W, class2=MA.LOESS.W, class3=MA.LOESS.SCALED.W, class4=MA.LOESS.QUANTILE.W, figTitles=c("Unnormalized data (FILTERED)","LOESS Normalization (FILTERED)","LOESS + SCALED (FILTERED)","LOESS + QUANTILE (FILTERED)"), non.zero.weight=TRUE, fileName="Boxplots_Overview_AValue_FILTERED_Reporters", y.axis="A")
} else {
  boxplotOverview1color(class1=MA.RAW.W, class2=MA.LOESS.W, class3=MA.SCALED, class4=MA.QUANTILE, figTitles=c("Unnormalized data (ALL)","LOESS Normalization (ALL)","SCALED Normalization (ALL)","QUANTILE Normalization (ALL)"), non.zero.weight=FALSE, fileName="Boxplots_Overview_ESTValue_ALL_Reporters")
  boxplotOverview1color(class1=MA.RAW.W, class2=MA.LOESS.W, class3=MA.SCALED, class4=MA.QUANTILE, figTitles=c("Unnormalized data (FILTERED)","LOESS Normalization (FILTERED)","SCALED Normalization (FILTERED)","QUANTILE Normalization (FILTERED)"), non.zero.weight=TRUE, fileName="Boxplots_Overview_ESTValue_FILTERED_Reporters")
}

## MA plots
cat("*------------\n| MvA plots\n*------------\n")
if(RG$datatype=="both") {
  CreateMAplots(first=MA.RAW.W, second=MA.LOESS.W, third=MA.LOESS.SCALED.W, fourth=MA.LOESS.QUANTILE.W, labels=c("Raw intensities", "LOESS normalization","LOESS + SCALED normalization", "LOESS + QUANTILE normalization")) 
} else {
  CreateMAplots(first=MA.RAW.W, second=MA.LOESS.W, third=NULL, fourth=NULL, labels=c("Raw intensities", "LOESS normalization")) 
}
cat("status: OK\n\n")


## saving normalized data, LOESS as default for dual channel data, Quantile for single channel data
if(RG$datatype=="both") {
  SaveData(MA.LOESS.W, save.RG=FALSE, save.MA=TRUE, save.QC=TRUE, detailed.QC=TRUE, raw=RG)
} else {
  SaveData(MA.QUANTILE, save.RG=TRUE, save.MA=FALSE, save.QC=TRUE, detailed.QC=TRUE, raw=RG)
} 

cat(paste("\n\n#-- ARRAYQC COMPLETED --# \nAll QC files generated by this script can be found in ", workpath, " !\n\n", sep=""))

###################
## END OF SCRIPT ##
###################
