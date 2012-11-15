## One and two channel QC - Dept. of Bioinformatics - BiGCaT - Maastricht University - the Netherlands
# Version: 1.1
# Last adjustment: 12-10-2010
## Function to look into: txtProgressBar

## Note: Changed xBGMeanSignal to xBGMediansignal, to reflect the xMedianSignal parameter (it basically means this).
## need to check if any other functions are affected.

## Suggestion: change softwarepackage to data.format? Fit better for generic.

#####################
## START OF SCRIPT ##
#####################

## VALIDATING USER ARGUMENTS
softwarepackage <- tolower(softwarepackage)
datatype <- tolower(datatype)

softwarepackage <- match.arg(softwarepackage, c("fes", "genepix", "generic"))
  if(softwarepackage == "fes") { package <- "agilent" }
  if(softwarepackage == "genepix") { package <- "genepix" }
  if(softwarepackage == "generic") { package <- "generic" } 
datatype <- match.arg(datatype, c("two-channel","green","red"))

cat("* Entering data directory...")
setwd(datapath)
cat(" done.\n")

cat("* Creating output directory...")
dir.create(workpath, showWarnings = FALSE)
cat(" done.\n")


## READING IN LIBRARIES
library(RColorBrewer)
library(limma)
library(affy)
#library(plotrix)

# Reading specific R-scripts (use the reload()-function if you want to reload the scripts later on, e.g. after modifications)
reload <- function() {
  source(paste(scriptpath,"ReadFiles.R", sep=""))
  source(paste(scriptpath,"QC.R", sep=""))
  source(paste(scriptpath,"CreateQCPlots.R", sep=""))
  source(paste(scriptpath,"functions.EListRaw.R", sep=""))
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

if (softwarepackage == "fes") {

  switch(datatype,
    "red" = {
      columns <- list( 
        R = "rMeanSignal", 
        Rb = "rBGUsed")
    },
    "green" = {
      columns <- list( 
        G = "gMeanSignal", 
        Gb = "gBGUsed")
    },
    {
      columns <- list( 
        R = "rMeanSignal", 
        G = "gMeanSignal", 
        Rb = "rBGUsed", 
        Gb = "gBGUsed")
    }
  )

  switch(datatype,
    "red" = {
      other.columns <- list(
        rBGMedianSignal = "rBGMeanSignal",
        rMedianSignal = "rMedianSignal",
        rNumPix = "rNumPix",
        rIsWellAboveBG = "rIsWellAboveBG",
        rIsSaturated = "rIsSaturated",
        manualFlags = "IsManualFlag")
    },
    "green" = {
      other.columns <- list(
        gBGMeanSignal = "gBGMeanSignal",
        gMedianSignal = "gMedianSignal",
        gNumPix = "gNumPix",
        gIsWellAboveBG = "gIsWellAboveBG",
        gIsSaturated = "gIsSaturated",
        manualFlags = "IsManualFlag")
    },
    {
      other.columns <- list(
        rBGMeanSignal = "rBGMeanSignal",
        gBGMeanSignal = "gBGMeanSignal",
        rMedianSignal = "rMedianSignal",
        gMedianSignal = "gMedianSignal",
        rNumPix = "rNumPix",
        gNumPix = "gNumPix",
        rIsWellAboveBG = "rIsWellAboveBG",
        gIsWellAboveBG = "gIsWellAboveBG",
        rIsSaturated = "rIsSaturated",
        gIsSaturated = "gIsSaturated",
        manualFlags = "IsManualFlag")
    }
  )

  annotation <- list(
    FeatureNum = "FeatureNum",
    Row = "Row",
    Col = "Col",
    ProbeName = "ProbeName",
    ControlType = "ControlType",
    GeneName = "GeneName",
    Description = "GeneName",
    SystematicName = "SystematicName")
}

if (softwarepackage == "genepix") {

  switch(datatype,
    "red" = {
      columns <- list( 
        R = "F635 Mean", 
        Rb = "B635 Mean") 
    },
    "green" = {
      columns <- list( 
        G = "F532 Mean", 
        Gb = "B532 Mean")
    },
    {
      columns <- list( 
        R = "F635 Mean", 
        G = "F532 Mean", 
        Rb = "B635 Mean", 
        Gb = "B532 Mean")
    }
  )

  switch(datatype,
    "red" = {
      other.columns <- list(
        rMedianSignal = "F635 Median",
        rBGSDev = "B635 SD",
        rPercSatPix = "F635 % Sat.",
        if(!is.null(useAsManualFlags))
          manualFlags = "Flags")
    },
    "green" = {
      other.columns <- list(
        gMedianSignal = "F532 Median", 
        gBGSDev = "B532 SD",
        gPercSatPix = "F532 % Sat.",
        manualFlags = "Flags")
    },
    {
      other.columns <- list(
        rMedianSignal = "F635 Median",
        gMedianSignal = "F532 Median", 
        rBGSDev = "B635 SD",
        gBGSDev = "B532 SD",
        rPercSatPix = "F635 % Sat.",
        gPercSatPix = "F532 % Sat.",
		    if(!is.null(useAsManualFlags))
          manualFlags = "Flags")
    }
  )

  annotation <- list(
    Block = "Block",
    Row = "Row",
    Col = "Column",
    ProbeName = "ID",
    GeneName = "Name")
}

## For "generic", a minimum number of columns is needed for arrayQC to run.
## All information is not always available to perform spot-specific QC, and as such
## the scripts should automatically detect this.

## It will be difficult to guess which columns are needed, but for basic arrayQC functionality the needed
## columns should be:
## - R/G Foreground + R/G Background [Intensities]
## - Reporter ID + Gene Annotation   [Annotation]
## - Row / Column                    [Array Layout + Block Detection]


####################
## For now, the "generic" format is required to trigger the next part.
## The idea would be that arrayQC uses a file where the proper column names are defined.
## For "fes" and "genepix" the data formats are already available for download at the arrayQC server
## (http://svn.bigcat.unimaas.nl/r-packages/arrayQC/dataformats/)

## For generic file formats, a local file will be created, where the user has to fill in specific fields
## OFcourse, it is not possible to fill in all fields. Only thing we DO need to do for every function, is to
## detect if the required columns are present (only for Generic)

## By removing the first if-statement, the idea would be that if "fes" or "genepix" would be chosen and the file
## does not exist on your local hard drive, then it will be downloaded to your hard drive. Once that is done, there
## are two possibilities: 1) continue with the default (recommended) settings; 2) stop arrayQC, fill in the correct
## fields and restart arrayQC.
## If the file DOES exist, it will just be parsed and tested just as the generic file.

## Ofcourse the aim of arrayQC is to make automated workflows possible. As such, users are able to use
## the same file for other purposes. We'll implement another variable in the main script called generic.settings.file
	
## Remove the next if-statement if files are ready
if(softwarepackage == "generic") {
  default.description.file <- "description_generic_format.arrayQC"
  if(is.null(description.file)) { description.file <- default.description.file }

  if(CheckFile(description.file, xpath=datapath)==0) {
    switch(softwarepackage,
      "generic" = {
        cat(paste("\n\n-> A new default description file", default.description.file, " has been generated in", datapath, ".\n"))
        cat("Feel free to edit and/or rename this file. \n\n\tNote: if a different filename is used, please refer to this file \n\tby adjusting the  >> description.file << parameter in the original\n\tarrayQC.R script!\n\n\tExample: description.file <- \"new.filename.txt\"")
        GenerateDescriptionFile(xpath=datapath, fileName = default.description.file, datatype = datatype)
      },

      "fes" = {
        cat("\n -> Retrieving default Feature Extraction Format from arrayQC server.\n")
        download.file(paste(scriptpath, "dataformats/description_fes_format.arrayQC", sep=""), 
                      paste(datapath, "description_fes_format.arrayQC", sep=""), 
                      method="wget")
        cat("\n -> File grabbed successfully. Do you want to continue with the default format settings? (recommended)\n")
        switch(
## Bug: Menu doesn't seem to do anything when "No" has been chosen...
          menu(c("Yes", "No")+1, 
            temp <- parseDescriptionFile(xpath = datapath, fileName = "description_fes_format.arrayQC", datatype = datatype, softwarepackage = softwarepackage),
            stop("\n\n>> arrayQC interrupted by user's request <<\n"))
          )
      },

      "genepix" = {
        cat("\n -> Retrieving default GenePix Pro Format from arrayQC server.\n")
        download.file(paste(scriptpath, "dataformats/description_genepix_format.arrayQC", sep=""), 
                      paste(datapath, "description_genepix_format.arrayQC", sep=""), 
                      method="wget")
        switch(
          menu(c("Yes", "No")+1, 
            temp <- parseDescriptionFile(xpath = datapath, fileName = "description_fes_format.arrayQC", datatype = datatype, softwarepackage = softwarepackage),
            stop("\n\n>> arrayQC interrupted by user's request <<\n"))
          )
      }

    )

  } else {
    temp <- parseDescriptionFile(xpath = datapath, fileName = description.file, datatype = datatype, softwarepackage = softwarepackage)
  }
  columns <- temp$columns
  other.columns <- temp$other.columns
  annotation <- temp$annotation
  rm(temp)
}

## Pre-check to see if basic files (spottypes and description) are present.
# - If spottypes file is not present --> continue
# - If description file is not present --> stop

if(!file.exists(paste(datapath, "spottypes.txt", sep=""))) { 
  spotfile <- NULL
} else {
  spotfile <- "spottypes.txt"
}

if(!file.exists(paste(datapath, "description.txt", sep=""))) { 
  descfile <- NULL
  cat(paste("\n>>WARNING<<\n\nThe \"description.txt\" file was not located in", datapath, "and is necessary for arrayQC to proceed.\n"))
  cat("If this file has a different name, please enter the correct name below:\n  ")
  descfile <- scan(n=1, what="character")
  check <- CheckFile(descfile, xpath = datapath)
  if(check == 0 || nchar(descfile) < 2) { stop(paste("\nFile \"", datapath, descfile, "\" could not be found. Terminating arrayQC script.", sep="")) }
  if(check == 1) { cat(paste(descfile, "successfully located. Proceeding with arrayQC...\n")) }
} else {
  descfile <- "description.txt"
}

## Reading in scanner output files and save this in the RG-object (makes use of read.maimages()):

setwd(workpath)
RG <- ReadFiles(description.file=descfile, spottypes.file=spotfile, data.path=datapath, columns=columns, other=other.columns, annotation=annotation, blocks=NULL, source=package, use.description=TRUE, save.backup=TRUE, manual.flags=useAsManualFlags)

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
if(RG$source != "genepix") {
  switch(RG$datatype,
    "red" = {
      RG <- QualityControl(RG, num.pix=NULL, low.pix=TRUE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-5,-5,-9,-5), criteria=c(-5,-5,-5,-5))
    },
    "green" = {
      RG <- QualityControl(RG, num.pix=NULL, low.pix=TRUE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-3,-3,-9,-3), criteria=c(-3,-3,-3,-3))
    },
    {
      RG <- QualityControl(RG, num.pix=NULL, low.pix=TRUE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-3,-3,-9,-3), criteria=c(-3,-3,-8,-8))
    }
  )
} else {
  switch(RG$datatype,
    "red" = {
      RG <- QualityControl(RG, num.pix=NULL, low.pix=FALSE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-5,-9,-5), criteria=c(-5,-5,-5))
    },
    "green" = {
      RG <- QualityControl(RG, num.pix=NULL, low.pix=FALSE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-3,-9,-3), criteria=c(-3,-3,-3))
    },
    {
      RG <- QualityControl(RG, num.pix=NULL, low.pix=FALSE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=c(-3,-9,-3), criteria=c(-3,-8,-8))
    }
  )
}

cat("*------------\n| Generating Summary files\n*------------\n")
## AddSummary
RG <- AddSummary(RG)

## CreateSummaryPlots creates barcharts from the information in the summary object.
CreateSummaryPlots(RG)

## Normalization -- Current set up is by using normalization weights only. Set weights=NULL to remove normalization weights.
cat("*------------\n| Data Normalization\n*------------\n")
if(RG$datatype=="both") {
  MA.RAW.W <- normalizeWithinArrays(RG, RG$printer, method="none", bc.method="subtract", offset=0, weights=NULL)
  MA.LOESS.W <- normalizeWithinArrays(RG, RG$printer, method="loess", iterations=4, bc.method="subtract", offset=0, weights=RG$other$weights.norm)
  MA.LOESS.SCALED.W <- normalizeBetweenArrays(MA.LOESS.W, method="scale")
  MA.LOESS.QUANTILE.W <- normalizeBetweenArrays(MA.LOESS.W, method="quantile")
  MA.LOESS.AQUANTILE.W <- normalizeBetweenArrays(MA.LOESS.W, method="Aquantile")
} else {
  RG2<-RG
  class(RG2) <- "RGList"
  if(RG$datatype=="red") {
    RG2$R <- RG2$G <- RG2$E
    RG2$Rb <- RG2Gb <- RG2$Eb
    RG2$G[] <- apply(RG2$R,1,mean,na.rm=TRUE)
    RG2$Gb[] <- apply(RG2$Rb,1,mean,na.rm=TRUE)
    MA.RAW.W <- normalizeWithinArrays(RG2, RG2$printer, method="none", bc.method="subtract", offset=0, weights=NULL)
    MA.RAW.W$other$EST <- MA.RAW.W$A+MA.RAW.W$M/2
    MA.LOESS.W <- normalizeWithinArrays(RG2, RG2$printer, method="loess", iterations=4, bc.method="subtract", offset=0, weights=RG2$other$weights.norm)
    MA.LOESS.W$other$EST <- MA.LOESS.W$A+MA.LOESS.W$M/2
  }
  if(RG$datatype=="green") {
    RG2$G <- RG2$R <- RG2$E
    RG2$Gb <- RG2$Rb <- RG2$Eb
    RG2$R[] <- apply(RG2$G,1,mean,na.rm=TRUE)
    RG2$Rb[] <- apply(RG2$Gb,1,mean,na.rm=TRUE)
    MA.RAW.W <- normalizeWithinArrays(RG2, RG2$printer, method="none", bc.method="subtract", offset=0, weights=NULL)
    MA.RAW.W$other$EST <- MA.RAW.W$A-MA.RAW.W$M/2
    MA.LOESS.W <- normalizeWithinArrays(RG2, RG2$printer, method="loess", iterations=4, bc.method="subtract", offset=0, weights=RG2$other$weights.norm)
    MA.LOESS.W$other$EST <- MA.LOESS.W$A-MA.LOESS.W$M/2
  }
  MA.SCALED <- normalizeBetweenArrays(MA.RAW.W$other$EST, method="scale")
  MA.QUANTILE <- normalizeBetweenArrays(MA.RAW.W$other$EST, method="quantile")
  #following codes are very similar to doing normalizeWithinArrays with method="none" and bc.method="none"
  #followed by normlizeBetweenArrays on the $EST field. However, it does not allow for bc correction
  #so we still keep the workaround
  #MA.SCALED2 <- normalizeBetweenArrays(RG, method="scale")
  #MA.QUANTILE2 <- normalizeBetweenArrays(RG, method="quantile")
  #other combinations of methods, not used by default
  #MA.LOESS.SCALED.W <- normalizeBetweenArrays(MA.LOESS.W$other$EST, method="scale")
  #MA.LOESS.QUANTILE.W <- normalizeBetweenArrays(MA.LOESS.W$other$EST, method="quantile")
  rm(RG2)
}

## Plotting of virtual array images
cat("status: OK\n\n")
cat("*------------\n| Virtual Images\n*------------\n")
switch(RG$datatype,
  both = {
    cat("* Based on red channel ...\n")
    imageplot3by2Adp(RG, RG$R, "RedSignal", high="red", low="yellow")
    imageplot3by2Adp(RG, RG$Rb, "RedBckgrndSignal", high="red", low="yellow")
    cat("* Based on green channel ...\n")
    imageplot3by2Adp(RG, RG$G, "GreenSignal", high="#00aa00", low="yellow")
    imageplot3by2Adp(RG, RG$Gb, "GreenBckgrndSignal", high="#00aa00", low="yellow")
  },
  red = {
    cat("* Based on red channel ...\n")
    imageplot3by2Adp(RG, RG$E, "RedSignal", high="red", low="yellow")
    imageplot3by2Adp(RG, RG$Eb, "RedBckgrndSignal", high="red", low="yellow")
  },
  green = {
    cat("* Based on green channel ...\n")
    imageplot3by2Adp(RG, RG$E, "GreenSignal", high="#00aa00", low="yellow")
    imageplot3by2Adp(RG, RG$Eb, "GreenBckgrndSignal", high="#00aa00", low="yellow")
  }
)

if(RG$datatype!="green") {
  if(!is.null(RG$other$rNumPix)) imageplot3by2Adp(RG, RG$other$rNumPix, "RedNumberPixels", high="red", low="yellow")
  if(!is.null(RG$other$rIsWellAboveBG)) imageplot3by2Adp(RG, RG$other$rIsWellAboveBG, "RedWellAboveBckgrnd", high="red", low="white")
  if(!is.null(RG$other$rIsSaturated)) imageplot3by2Adp(RG, RG$other$rIsSaturated, "RedSaturated", high="red", low="white")
}

if(RG$datatype!="red") {
  if(!is.null(RG$other$gNumPix)) imageplot3by2Adp(RG, RG$other$gNumPix, "GreenNumberPixels", high="#00aa00", low="yellow")
  if(!is.null(RG$other$gIsWellAboveBG)) imageplot3by2Adp(RG, RG$other$gIsWellAboveBG, "GreenWellAboveBckgrnd", high="green", low="white")
  if(!is.null(RG$other$gIsSaturated)) imageplot3by2Adp(RG, RG$other$gIsSaturated, "GreenSaturated", high="green", low="white")
}

cat(" done.\n")
cat("* Based on normalized data ...\n") 
if(RG$datatype=="both") {
  imageplot3by2Adp(MA.RAW.W, MA.RAW.W$M, "RAW_LogRatio", high="blue", low="yellow", symm=TRUE)
  imageplot3by2Adp(MA.RAW.W, MA.RAW.W$A, "RAW_AverageIntensity", high="blue", low="yellow", symm=FALSE)
  imageplot3by2Adp(MA.LOESS.W, MA.LOESS.W$M, "LOESS_LogRatio", high="blue", low="yellow", symm=TRUE)
  imageplot3by2Adp(MA.LOESS.W, MA.LOESS.W$A, "LOESS_AverageIntensity", high="blue", low="yellow", symm=FALSE)
} else {
  imageplot3by2Adp(MA.RAW.W, MA.RAW.W$other$EST, "EstimatedSignal", high="blue", low="yellow", symm=FALSE)
  imageplot3by2Adp(MA.LOESS.W, MA.LOESS.W$other$EST, "LOESS_EstimatedSignal", high="blue", low="yellow", symm=FALSE)
  imageplot3by2Adp(MA.RAW.W, MA.SCALED, "Scaled_EstimatedSignal", high="blue", low="yellow", symm=FALSE)
  imageplot3by2Adp(MA.RAW.W, MA.QUANTILE, "Quantile_EstimatedSignal", high="blue", low="yellow", symm=FALSE)
}

if(!is.null(RG$other$manualFlags)) imageplot3by2Adp(RG, RG$other$manualFlags, "ManualFlags", high="sienna", low="orange")
#if(!is.null(RG$other$LowNumPix)) imageplot3by2Adp(RG, RG$other$LowNumPix, "LowNumPix", high="orange", low="sienna")
if(!is.null(RG$other$LowMeanMedianRatio)) imageplot3by2Adp(RG, RG$other$LowMeanMedianRatio, "LowMeanMedianRatio", high="orange", low="sienna")
if(!is.null(RG$other$NOTWellAboveBG)) imageplot3by2Adp(RG, RG$other$NOTWellAboveBG, "NOTWellAboveBG", high="orange", low="sienna")
#if(!is.null(RG$other$Saturated)) imageplot3by2Adp(RG, RG$other$Saturated, "Saturated", high="orange", low="sienna")
if(!is.null(RG$other$weights.norm)) imageplot3by2Adp(RG, RG$other$weights.norm, "WeightsForNormalization", high="purple2", low="white")
if(!is.null(RG$weights)) imageplot3by2Adp(RG, RG$weights, "Weights", high="purple2", low="white")

cat(" done.\n\nstatus: OK\n\n")
cat("*------------\n| Clustering\n*------------\n")
## Cluster plots - by default: dist.method="euclidean", clust.method="ward". Can be changed by other methods available in ?hdist and ?hclust
if(RG$datatype=="both") {
  HierarchCluster(MA.RAW.W$M, main="Raw Data (RAW)")
  HierarchCluster(MA.LOESS.W$M, main="Normalized Data (LOESS)")
  HierarchCluster(MA.LOESS.SCALED.W$M, main="Normalized Data (LOESS+SCALED)")
  HierarchCluster(MA.LOESS.QUANTILE.W$M, main="Normalized Data (LOESS+QUANTILE)")
  HierarchCluster(MA.LOESS.AQUANTILE.W$M, main="Normalized Data (LOESS+AQUANTILE)")
} else {
  HierarchCluster(MA.RAW.W$other$EST, main="Raw Data (RAW)")
  HierarchCluster(MA.LOESS.W$other$EST, main="Normalized Data (LOESS)")
  HierarchCluster(MA.SCALED, main="Normalized Data (SCALED)")
  HierarchCluster(MA.QUANTILE, main="Normalized Data (QUANTILE)")
}

cat("status: OK\n\n")
cat("*------------\n| Correlation plots / heatmaps\n*------------\n")
## Heatmaps
CreateHeatMap(MA.RAW.W)
CreateHeatMap(MA.LOESS.W)
if(RG$datatype=="both") {
  CreateHeatMap(MA.LOESS.SCALED.W)
  CreateHeatMap(MA.LOESS.QUANTILE.W)
  CreateHeatMap(MA.LOESS.AQUANTILE.W)
} else {
  CreateHeatMap(MA.SCALED)
  CreateHeatMap(MA.QUANTILE)
}
cat("status: OK\n\n")

## PCA plots
cat("*------------\n| PCA plots\n*------------\n")
CreatePCAplot(MA.RAW.W)
CreatePCAplot(MA.LOESS.W)
if(RG$datatype=="both") {
  CreatePCAplot(MA.LOESS.SCALED.W)
  CreatePCAplot(MA.LOESS.QUANTILE.W)
  CreatePCAplot(MA.LOESS.AQUANTILE.W)
} else {
  CreatePCAplot(MA.SCALED)
  CreatePCAplot(MA.QUANTILE)
}
cat("status: OK\n\n")

## MA plots
cat("*------------\n| MvA plots\n*------------\n")
if(RG$datatype=="both") {
  CreateMAplots(first=MA.RAW.W, second=MA.LOESS.W, third=MA.LOESS.SCALED.W, fourth=MA.LOESS.QUANTILE.W, labels=c("Raw intensities", "LOESS normalization","LOESS + SCALED normalization", "LOESS + QUANTILE normalization")) 
} else {
  CreateMAplots(first=MA.RAW.W, second=MA.LOESS.W, third=NULL, fourth=NULL, labels=c("Raw intensities", "LOESS normalization")) 
}
cat("status: OK\n\n")

## Density Plots
cat("*------------\n| Density Plots\n*------------\n")
CreateDensityPlots(RG, name = "Raw Data")
CreateDensityPlots(MA.RAW.W, name = "Raw BG Corrected Data")
CreateDensityPlots(MA.LOESS.W, name = "LOESS Normalized Data")
if(RG$datatype=="both") {
    CreateDensityPlots(MA.LOESS.SCALED.W, name = "LOESS and SCALED Data")
    CreateDensityPlots(MA.LOESS.QUANTILE.W, name = "LOESS and QUANTILE Data")
    CreateDensityPlots(MA.LOESS.AQUANTILE.W, name = "LOESS and AQUANTILE Data")
} else {
  CreateDensityPlots(MA.SCALED, name = "SCALED Normalized Data")
  CreateDensityPlots(MA.QUANTILE, name = "QUANTILE Normalized Data")
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

## saving normalized data, LOESS as default for dual channel data, Quantile for single channel data
if(RG$datatype=="both") {
  SaveData(MA.LOESS.W, save.RG=FALSE, save.MA=TRUE, save.QC=TRUE, detailed.QC=TRUE, raw=RG)
} else {
  SaveData(MA.QUANTILE, save.RG=TRUE, save.MA=FALSE, save.QC=TRUE, detailed.QC=TRUE, raw=RG)
} 

cat(paste("\n\n[[COMPLETED]]\nAll QC files generated by this scriptcan be found in ", workpath, " !\n\n", sep=""))

###################
## END OF SCRIPT ##
###################
