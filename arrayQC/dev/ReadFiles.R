## NEED ONE MAJOR HAUL OF SCRIPTPATH! Should be renamed to arrayQC.scriptpath to avoid confusion with other scripts!

################################
##     ReadFiles function     ##
################################
## 08-08-2012
#  - fixed wrong reference (RG$source) in makeLimmaCompatible()
## 25/07/2012
#  - Added .NAcheck for generic formats

## 16/05/2012
#  - ReadFiles() adjustments:
#  -- If an annotation column is missing, ReadFiles() will now check if it's a required annotation field.
#     If it's not required, the script will continue.

## 01/05/2012
#  - Implemented first part of major 'general' applicability of data.
#  - Added required fields parameter that will be initiated during the initial arrayQC check.

## 25/01/2012
#  - ReadFiles() adjustments:
#    - "Generic" files may possible contain additional rows of NA information. ReadFiles will now immediately correct for this
#      if this is occurring.
#    - Optimized detection of block columns
#      - max(row) * max(col) * total.blocks == total datapoints (5% margin) --> overwrite row and col
#      - max(row) * max(col) == total datapoints --> keep row, check if block column is numeric

## 14/12/2011 (Stan - v1.2.2
## - Added .checkReqColumns() for each function name.
## - ReadFiles now supports ControlType for "generic" formats
##   users must provide a controlType.value parameter in the main script mentioning
##   all possible names of the controls. These will then be flagged and removed
##   from subsequent analysis.
## 30/11/2011 (Stan) - v1.2.1
## - Rename and bugfixed .FileExists() and GenerateDescriptionFile()

## 28/07/2011 (Stan)
## - Added switch-parameter for calculating block size / orientation

## 13/12/2010 (Stan)
## -- Worked out working with 'files' that contain the columns that arrayQC wants to use.
## -- Works now for "generic" format, semi-support for "fes" and "genepix" output as well.
## -- Idea would be that if user selects "fes", it will check if a file is on the local harddrive with settings,
##    if not, the file will be downloaded and stored on the local hard drive. Downside of this approach is that
##    arrayQC needs an internetconnection for it's first real run!

## Hehe, LOL...
## read.maimages maakt nu geen gebruik meer van channels, maar van green.only <- TRUE/FALSE parameter... (R 2.11.1)
## Dwz dat arrayQC op die manier niet meer backwards compatibel is met een oudere R versie (of we moeten dit opvangen met een if-statement)

## Sommige stukken aangepast, werken niet optimaal:
## - nspots deel, enkel van toepassing op niet-FES arrays
## - checken van kolomnamen in bestand (gaat voor FES mis)

## To do:
## if datatype <- "dual-channel" and only one channel is present, give a specific warning!

require(limma)

.NAcheck <- function(x) {
  countNA <- apply(x, 1, function(x) { sum(is.na(x)) })
  selectNA <- which(countNA[] == dim(x)[2])
  if(length(selectNA) != 0) { x[-c(selectNA),] }
  return(x)
}

.FileExists <- function(x, xpath=NULL) {
  if(is.null(x)) { 
    check <- FALSE 
  } else {
    if(is.null(xpath)) { xpath <- getwd() }
    if(file.exists(paste(xpath, x, sep=""))) { 
      check <- TRUE 
    } else { 
      cat(paste(" -->> Unable to locate file ", x, " in ", xpath, "\n", sep="")); check <- FALSE 
    }
  }
return(check)
}

GenerateHeaderFile <- function(xpath = NULL, fileName = "columnHeaders.arrayQC", datatype = NULL, package = NULL) {
## (25-07-2012)
# - set default value of fileName to "columnHeaders.arrayQC"  
  if(is.null(datatype)) { stop("Define datatype!") }
  if(is.null(xpath)) { xpath <- getwd() }
  temp.names <- c("MedianSignal", "BGMedianSignal", "BGSDev", "PercSatPix", "NumPix", "IsSaturated", "IsWellAboveBG")
  temp.reqA <- c("Y", "Y", "N", "N", "N", "N", "N", "N", "N")
  temp.descr1 <- c("foreground intensity (i.e. mean)", 
                  "background intensity (i.e. mean)", 
                  "foreground intensity (i.e. median)",
                  "background intensity (i.e. median)",
                  "background SD value", "percentage of pixels saturated",
                  "number of used spot pixels to estimate the signal",
                  "flag column for saturated (1) or unsaturated (0) reporters",
                  "flag column for well measured (1) or badly measured (0) reporters")
  temp.descr2 <- c("fes/genepix", "fes/genepix", "fes/genepix", "fes/genepix", "genepix", "genepix", "fes", "fes", "fes")

 aQC <- NULL
  ## Green channel
 aQC$Names       <- c(aQC$Names, "G", "Gb", paste("g", temp.names, sep=""))
 aQC$Required    <- c(aQC$Required, temp.reqA)
 aQC$Description1 <- c(aQC$Description1, paste("Green", temp.descr1))
 aQC$Description2 <- c(aQC$Description2, temp.descr2) 
  ## Red channel
 aQC$Names       <- c(aQC$Names, "R", "Rb", paste("r", temp.names, sep=""))
 aQC$Required    <- c(aQC$Required, temp.reqA)
 aQC$Description1 <- c(aQC$Description1, paste("Red", temp.descr1))
 aQC$Description2 <- c(aQC$Description2, temp.descr2)
  ## Annotation
 aQC$Names        <- c(aQC$Names, "---", "FeatureNum", "Block", "Row", "Col", "ProbeName", "ControlType", "GeneName", "Description", "SystematicName", "manualFlags")
 aQC$Required     <- c(aQC$Required, "---", "N","N", "Y","Y","Y","N","N","N","N","N")
 aQC$Description1 <- c(aQC$Description1, "---", "Number of the reporter (often a combination of ([row] * [col]) * blocknr", 
                               "Block number", "Row number of the reporter", 
                               "Column number of the reporter",
                               "Platform-specific reporter ID",
                               "Numerical column - Is reporter a control reporter (i.e. 0 = no control)",
                               "Target gene name",
                               "Gene description",
                               "Additional gene name / annotation",
                               "Column used for manual flagging")
  aQC$Description2 <- c(aQC$Description2, "", "fes/genepix", "", rep("fes/genepix",7), "")


  x <- matrix(ncol=4, data=cbind(c(""), aQC$Required, aQC$Description1, aQC$Description2), nrow=length(aQC$Names), 
                dimnames=list(Fields=aQC$Names, Title=c("columnHeader","Required","fieldDescription", "packageCompatibility")))

  switch(package, "genepix" = { 
     x[,1] <- c("F532 Mean", "B532 Mean", "F532 Median", "", "B532 SD", 
                "F532 % Sat.", "", "", "", "F635 Mean", "B635 Mean", 
                "F635 Median", "", "B635 SD", "F635 % Sat.", "", "" ,"", "", 
                "RefNumber", "", "Row", "Column", "ID" , "ControlType",
                "Name", "Description", "GeneName", "Flags")                      
  }, "agilent" = {
    x[,1] <- c("gMeanSignal","gBGUsed","gMedianSignal","gBGMedianSignal","","","gNumPix","gIsSaturated",
               "gIsWellAboveBG","rMeanSignal","rBGUsed","rMedianSignal","rBGMedianSignal","","","rNumPix",
               "rIsSaturated","rIsWellAboveBG","","FeatureNum","", "Row","Col","ProbeName","ControlType","GeneName","Description",
               "SystematicName","IsManualFlag")
  })
  write.table(x, file=paste(xpath, fileName, sep=""), col.names=NA, quote=FALSE, sep="\t")
}

parseHeaderFile <- function(xpath=NULL, fileName="columnHeaders.arrayQC", datatype=NULL, package=NULL) {
  ## xpath= 		location of the description file
  ## fileName=		annotation header file name
  ## datatype=		"green","red" or "dual-channel"
  ## package=	        "agilent", "genepix" or "generic"
  match.arg(package, c("agilent", "genepix", "generic"))
  match.arg(datatype, c("green", "red", "two-channel"))

  funct.error <- NULL ## Captures faulty function paramter
  if(is.null(xpath)) { funct.error <- c(funct.error, "xpath") }
  if(is.null(datatype)) { funct.error <- c(funct.error, "datatype") }
  if(is.null(package)) { funct.error <- c(funct.error, "package") }
  if(!is.null(funct.error)) { stop(paste(" -->> Error! The following function parameter is invalid / missing:", funct.error)) }

  x <- as.matrix(read.table(file=paste(xpath, fileName, sep=""), row.names=1, sep="\t", header=TRUE))
  required <- list()
## Which columns are required? That depends on two things: a) the type of dataset (here it's generic, but
## let's assume that in a near future we'll extend this functionality to the 'agilent' and 'genepix' output
## as well.; b) the datatype.
  switch(package,
    "generic" = {
      required$annotation <- c("Row","Col","ProbeName")

      if(datatype == "green" || datatype == "two-channel") {
        required$columns <- c(required$columns, "G", "Gb")
        required$other.columns <- c(required$other.columns, "gBGSDev")
      }
      if(datatype == "red" || datatype == "two-channel") {
        required$columns <- c(required$columns, "R", "Rb")
        required$other.columns <- c(required$other.columns, "rBGSDev")
      }
    },
    "agilent" = {
      required$annotation <- c("FeatureNum", "Row","Col","ProbeName", "ControlType", "GeneName", "Description", "SystematicName", "manualFlags")

      if(datatype == "green" || datatype == "two-channel") {
        required$columns <- c(required$columns, "G", "Gb")
        required$other.columns <- c(required$other.columns, "gBGMedianSignal", "gMedianSignal", "gNumPix", "gIsWellAboveBG", "gIsSaturated")
      }
      if(datatype == "red" || datatype == "two-channel") {
        required$columns <- c(required$columns, "R", "Rb")
        required$other.columns <- c(required$other.columns, "rBGMedianSignal", "rMedianSignal", "rNumPix", "rIsWellAboveBG", "rIsSaturated")
      }
    },
    "genepix" = {
      required$annotation <- c("FeatureNum", "Row","Col","ProbeName", "ControlType", "GeneName", "Description", "SystematicName", "manualFlags")

      if(datatype == "green" || datatype == "two-channel") {
        required$columns <- c(required$columns, "G", "Gb")
        required$other.columns <- c(required$other.columns, "gMedianSignal", "gBGSDev", "gPercSatPix")
      }
      if(datatype == "red" || datatype == "two-channel") {
        required$columns <- c(required$columns, "R", "Rb")
        required$other.columns <- c(required$other.columns, "rMedianSignal", "rBGSDev", "rPercSatPix")
      }
    })

  y <- NULL
  y$descriptionFile <- x
  switch(datatype,
## Is "red" and "green" obsolete with the new $E object? If it is, Lars will need to adapt this part.
## i.e. R, G, Rb, and Gb need to be changed to E and Eb.
    "red" = {
      y$columns <- list( 
        R  = x["R",1],
        Rb = x["Rb",1]
      )
      y$other.columns <- list(
        rMedianSignal = x["rMedianSignal",1],
        rBGMedianSignal = x["rBGMedianSignal",1],
        rBGSDev = x["rBGSDev",1],
        rPercSatPix = x["rPercSatPix",1],
        rNumPix = x["rNumPix",1],
        rIsSaturated = x["rIsSaturated",1],
        rIsWellAboveBG = x["rIsWellAboveBG",1]
      )
    },
    "green" = {
      y$columns <- list( 
        G  = x["G",1],
        Gb = x["Gb",1]
      )
      y$other.columns <- list(
        gMedianSignal = x["gMedianSignal",1],
        gBGMedianSignal = x["gBGMedianSignal",1],
        gBGSDev = x["gBGSDev",1],
        gPercSatPix = x["gPercSatPix",1],
        gNumPix = x["gNumPix",1],
        gIsSaturated = x["gIsSaturated",1],
        gIsWellAboveBG = x["gIsWellAboveBG",1]
      )
    },"two-channel" = {
      y$columns <- list( 
        R = x["R",1], 
        G = x["G",1], 
        Rb = x["Rb",1], 
        Gb = x["Gb",1]
      )
      y$other.columns <- list(
        rMedianSignal = x["rMedianSignal",1],
        rBGMedianSignal = x["rBGMedianSignal",1],
        rBGSDev = x["rBGSDev",1],
        rPercSatPix = x["rPercSatPix",1],
        rNumPix = x["rNumPix",1],
        rIsSaturated = x["rIsSaturated",1],
        rIsWellAboveBG = x["rIsWellAboveBG",1],
        gMedianSignal = x["gMedianSignal",1],
        gBGMedianSignal = x["gBGMedianSignal",1],
        gBGSDev = x["gBGSDev",1],
        gPercSatPix = x["gPercSatPix",1],
        gNumPix = x["gNumPix",1],
        gIsSaturated = x["gIsSaturated",1],
        gIsWellAboveBG = x["gIsWellAboveBG",1]
      )
    }
  )
  y$annotation <- list(
    FeatureNum = x["FeatureNum",1],
    Block = x["Block",1],
    Row = x["Row",1],
    Col = x["Col",1],
    ProbeName = x["ProbeName",1],
    ControlType = x["ControlType",1],
    GeneName = x["GeneName",1],
    Description = x["Description",1],
    SystematicName = x["SystematicName",1]
  )

  ## Removing all values that are not filled in.
  a <- unlist(y$other.columns)
  b <- which(a[]=="")
  if(length(b) > 0) { y$other.columns[b] <- NULL }
  if(length(a) == length(b)) { y$other.columns <- NULL }
  rm(a, b)

  a <- unlist(y$annotation)
  b <- which(a[]=="")
  if(length(b) > 0) { y$annotation[b] <- NULL } 
  cat(" [INFO] arrayQC will work with the intensity values derived from columns:\n")
  for(i in 1:length(y$columns)) { cat(paste(" -> ", y$columns[i], "\n", sep="")) }
  return(y)
}

ReadFiles <- function(description.file=NULL, spottypes.file=NULL, data.path=NULL, columns=NULL, other=NULL, annotation=NULL, blocks=NULL, source=NULL, use.description=FALSE, save.backup=TRUE, manual.flags=NULL, debug.parameter=FALSE, controlType.value=NULL, arrayQC.path="http://svn.bigcat.unimaas.nl/r-packages/arrayQC/") {
#-- To add for "generic" support:
#   number of lines to skip prior to reading in the files!
  error <- NULL
  if(is.null(description.file)) {error <- c(error,"description.file")}
  if(is.null(data.path)) {error <- c(error,"data.path")}
  if(is.null(columns)) {error <- c(error,"columns")}
  if(is.null(annotation)) {error <- c(error,"annotation")}
  if(is.null(annotation$Row)) {error <- c(error,"Row")}
  if(is.null(annotation$Col)) {error <- c(error,"Col")}
  if(is.null(source)) {error <- c(error,"source")}
  if(is.null(arrayQC.path)) {error <- c(error,"arrayQC.path")}

  if(!is.null(error)) {stop(paste(" The following parameters were not defined correctly: ", paste(error, collapse=", "), "\n", sep=""))}
  
#-- Debugging purposes

  if(debug.parameter==1) {
    cat("Paste the following in your R session window to debug the current function:\n")
    cat(" description.file <- experimentalDescr\n spottypes.file <- spotfile\n data.path <- datapath\n")
    cat(" columns <- columns\n other <- other.columns\n annotation <- annotation\n blocks <- NULL\n")
    cat(" source <- softwarepackage\n use.description <- TRUE\n save.backup <- TRUE\n manual.flags <- useAsManualFlags\n")
    cat(" controlType.value <- controlType.value\narrayQC.path <- scriptpath")
    stop("[/StartDebug]\n----")
  }

  # Check for blocks parameter:
  # - is it present?
  # - is it a list?
  # - does it contain at least 3 values?
  if(!is.null(blocks) & class(blocks) != "list" & length(blocks)<=2) {
    error <- NULL
    if(class(blocks) != "list") {  error <- c(error, "- blocks is no list!\n") }
    if(length(blocks)<= 2) { error <- c(error, "- blocks does not contain at least 3 values!\n") }
    stop(paste("The following errors were encountered while checking the blocks parameter:\n", error, sep=""))
  }
  #check which channels are present
  present <- 0
  if(!is.null(columns$R))  { present <- present + 1 }
  if(!is.null(columns$Rb)) { present <- present + 1 }
  if(!is.null(columns$G))  { present <- present + 10 }
  if(!is.null(columns$Gb)) { present <- present + 10 }
  #check whether any channel is present
  if(present==0) { stop("red, green or both channels have to be specified\n") }
  #check whether background and channel are both present, or both not
  if(present %in% c(1,10,11,12,21)) { stop("both foreground and background values have to be specified for each channel present\n") }
  datatype <- "unknown"
  if(present == 2) { cat("\nDataset type: single color red channel\n\n");datatype <- "red"}
  if(present == 20) {cat("\nDataset type: single color green channel\n\n");datatype <- "green"}
  if(present == 22) {cat("\nDataset type: dual color\n\n");datatype <- "both"}


  #just in case check
  if((!(present %in% c(0,1,2,10,11,12,20,21,22))) | datatype=="unknown") { stop("incorrect channel check value\n") }

  # Check if default files are present (description.txt and spottypes.txt)
  cat(" * Checking for target and spottype files:\n")
  cat("   - Spottype file  ")
  if(!is.null(spottypes.file) & file.exists(file=paste(data.path, spottypes.file, sep=""))) {
    cat(paste(spottypes.file, "located\n"))
  } else {
    if(is.null(spottypes.file)) {
      cat(paste(spottypes.file, "undefined.\n"))
    } else {
      cat(paste(spottypes.file, "not found!\n"))
    }
  }

  cat("   - Description file: ")
  if(!is.null(description.file) & file.exists(file=paste(data.path, description.file, sep=""))) {
    cat(paste(description.file, "located\n"))
  } else {
    if(is.null(description.file)) {
      cat(paste(description.file, "undefined.\n"))
    } else {
      cat(paste(description.file, "not found!\n"))
    }    
    stop(" -- Script halted -- \n")
  }
 
  cat(" * Reading target file ...");

  targets <- readTargets(file=paste(data.path, description.file, sep=""), sep="\t", quote="\"");
  cat(" ok\n");
  
  #check whether target file is ok
  #first convert to correct case if needed
  names(targets)[grep("filename",tolower(names(targets)))] <- "FileName"
  names(targets)[grep("description",tolower(names(targets)))] <- "Description"
  
  if(is.null(targets$FileName)) { stop("Targets file must contain a FileName column") }
  #we will for now not check whether there are Cy3 and/or Cy5 columns, as the user may wish to name these otherwise
  if(is.null(targets$Description) & use.description) { 
    stop("Targets file must contain a Description column if you want to use this as header")
  }

  #check for manualFlags field, when manual.flags values are provided
  if(is.null(other$manualFlags)  & !is.null(manual.flags)) {
    stop("manual.flags are provided, without providing a manualFlags column, no flagging can be performed")
  }
  na.files <- NULL
  wd <- getwd()
  setwd(data.path)
  for (f in targets$FileName) {
    if(!file.exists(f))
      na.files <- c(na.files, f)
  }
  if(length(na.files)>0)
    stop("The following files cannot be found: ", paste(na.files,collapse=", "))

  #select the number of channels
  channels <- 2
  if((present == 2) | (present == 20)) { channels <- 1 }


  #check whether all columns to be included, are present in the (first) data file
  #make use of the fact that, if not found, read.maimages will just not include it in the object returned
  #change to errors stopping the code later (after debugging that no incorrect calls are made)
  cat(" * Checking column names (reading first file for verification)... ")

  ## 13/12/2010
  ## Extraction Version number of R to use the proper function parameters (i.e. channels or green.only)
  ## Actually, this needs to be done on the limma version, but I have no idea which version this would be...
  ver    <- as.character(R.version["minor"])
  subver <- strsplit(ver, ".", fixed=TRUE)[[1]][1]

## Redoing the read-part:
## - Read in the first xxx lines
## - If one match is found in a line, continue in that same line searching for other possible matches.
## - If at least 4 matches are found, then the correct header row has been located.

## We'll do this with the read.columns function and a small loop.

#  cat("\n * Additional column check ...\n")
  .addColumns <- function(x, y) {
     x <- c(x, as.vector(unlist(y)))
     return(x)
  }
  required.columns <- NULL
  required.columns <- .addColumns(required.columns, columns)
  required.columns <- .addColumns(required.columns, other)
  required.columns <- .addColumns(required.columns, annotation)

  required.colnames <- paste("columns$", names(columns), sep="")
  if(!is.null(other)) { required.colnames <- c(required.colnames, paste("other$", names(other), sep="")) }
  if(!is.null(annotation)) {  required.colnames <- c(required.colnames, paste("annotation$", names(annotation), sep="")) }

  maxLines <- 100
  file1 <- readLines(targets$FileName[1], n=maxLines)

  temp2 <- matrix(ncol=1+length(required.columns), data=NA, nrow=maxLines)
  for(i in 1:maxLines) {
    temp <- gsub("\"","",file1[i])
    temp <- strsplit(temp,"\t")[[1]]
    temp2[i,1] <- sum(required.columns %in% temp)
    temp2[i,2:dim(temp2)[2]] <- required.columns %in% temp
    rm(temp)
  }

  if(max(temp2[,1], na.rm=TRUE) == length(required.columns)) { 
    rm(temp2)
    cat(" ok\n") 
  } else {
    xxx <- which(temp2[,1] == max(temp2[,1], na.rm=TRUE))
    if(length(xxx)!=1) { stop("\n   [[ERROR]]\n No specific / multiple rows could be associated with the given column names!\n") }
    cat(paste("\n   [[ERROR]]\n  ReadFiles() was able to match ", max(temp2[,1], na.rm=TRUE), " given column names. The following column names were NOT found:\n"))
    check.required <- 1
## Check which row contains this value
    header.row <- which(temp2[,2] == max(temp2[,2]))
    mispos <- temp2[header.row,2:dim(temp2)[2]]
##### FUTURE UPDATE
## Need to implment a check here that looks for all 'optional' values. (grab them with .checkColumns). Then check if all missing values are optional, then proceed.
## ALso need to include that the objects themselves should be altered. If optional and not present, it should be removed from the object in general.
    reqCols <- .checkColumns(functionName="annotation", silent=TRUE, location=arrayQC.path)[["required"]]
    missing.col <- required.columns[mispos[]==0]
    missing.colname <- required.colnames[mispos[]==0]
    reqValue <- missing.col %in% reqCols
    ## If all values in reqValue are FALSE, then just proceed with the QC procedure
    if(sum(reqValue) == 0) { 
      cat(paste("*", missing.colname, ":", missing.col, "\n"))
      check.required <- 0
      cat("  --> These columns were all identified as NON-ESSENTIAL and as such arrayQC will continue to read through the files...\n")
    }
    if(check.required == 1) {
      if(xxx == 0) {
        cat(paste(" *", missing.colname, ":", required.columns, "\n"))
      } else {
        cat(paste(" *", missing.colname, ":", required.columns[!temp2[xxx,2:dim(temp2)[2]]], "\n"))
        Sys.sleep(5)
        temp <- strsplit(file1[xxx],"\t")[[1]]
  #      temp <- paste("*", temp)
    ## Have to fix this later, but the problem is that while creating a matrix and there are a surplus of columns, the data starts
    ## to repeat itself. Creating an empty matrix and filling it up row by row.
  #      temp <- matrix(data=temp, ncol=4, nrow=ceiling(length(temp)/4), byrow=TRUE)
  #     mat1 <- matrix(data="", ncol=4, nrow=ceiling(length(temp)/4))
  #     if(ceiling(length(temp)/4) == length(temp)/4) {
  #       temp <- matrix(data=temp, ncol=4, nrow=ceiling(length(temp)/4), byrow=TRUE)
  #     } else {
  #       for(i in 1:ceiling(length(temp)/4)) {
  #         start.pos <- (4*(i-1))+1
  #         end.pos <- 4*i
  #         if(i == ceiling(length(temp)/4)) {
  #           remaining <- length(temp) - (4 * (i-1))
  #           mat1[1:remaining,i] <- temp[ start.pos : (start.pos + remaining)]
  #         } else {
  #           mat1[1:4,i] <- temp[ start.pos : end.pos ]        
  #         }
  #       }
  #    }
  
        cat("\n\n   [[INFO]]\nThe following column names were found (and can be chosen):\n")
        cat(paste(" *", temp, "\n", sep=""))
        rm(temp, temp2)
        stop("ReadFiles stopped due to unmatching column names!")
      }
    }
  }
  cat(" * Reading microarray text files ...\n");
  if(subver > 10) {
    x <- read.maimages(targets$FileName, source=source, path=data.path, green.only=ifelse(channels==1, TRUE, FALSE), columns=columns, other.columns=other, annotation=annotation)
  } else {
    x <- read.maimages(targets$FileName, source=source, path=data.path, channels=channels, columns=columns, other.columns=other, annotation=annotation)
  }

  # Changing column names of x$other  and x$genes into names of other.columns and annotation respectively. 
  # This is needed if the user specified column names are different from the default ones. 
  # Since the script assumes these 'default' column names in upcoming functions, it is needed that the names are 
  # changed back to the defaults in this case
  if(!is.null(other)) { names(x$other) <- names(other[!sapply(other,is.null)]) }
  names(x$genes) <- names(annotation[!sapply(annotation,is.null)])

  if(source=="generic") {
    ## It is possible that in a generic format there might be additional lines in the end
    #  of the datafile that contain no expression data. These lines need to be removed prior
    #  to continue. This to ensure that the data dimensions are kept constant!
    if(channels == 1) {
      x$E  <- .NAcheck(x$E)
      x$Eb <- .NAcheck(x$Eb)
    } else {
      x$G  <- .NAcheck(x$G)
      x$R  <- .NAcheck(x$R)
      x$Gb <- .NAcheck(x$Gb)
      x$Rb <- .NAcheck(x$Rb)
    }
    x$genes <- .NAcheck(x$genes)
    if(!is.null(x$other)) {
      for(i in 1:length(names(x$other))) {
        x$other[i][[1]] <- .NAcheck(x$other[i][[1]])
      }
    }
  }

  ## For red channel data no E or Eb is being used... We're forcing R to do this.
  ## One thing that needs to be added is how to remove subparts of the data (i.e. rm(x$R))
  if(datatype == "red") {
    if(!is.null(x$R)) {
      x$E <- x$R
      x$Eb <- x$Rb
      #x$R <- x$Rb <- NULL
    }
  }
  if(datatype == "green") {
    if(!is.null(x$G)) {
      x$E <- x$G
      x$Eb <- x$Gb
      #x$G <- x$Gb <- NULL
    }
  }

  setwd(wd)

  warningMessages <- NULL

  #if no feature numbers are present, create artificial ones 1:number of features
  if(is.null(annotation$FeatureNum)) { c(warningMessages, "[[WARNING]] No feature number column defined, creating artificial numbers 1:number of features")}
  if(is.null(x$genes$FeatureNum)) {x$genes$FeatureNum<-1:length(x$genes$Row)}
  #if no ControlType column is given, consider all probesets to be real genes
  if(is.null(annotation$ControlType)) {
    warningMessages <- c(warningMessages, "[[WARNING]] No control type column defined, no probesets marked as controls")
    if(is.null(x$genes$ControlType)) {if(source=="genepix") {x$genes$ControlType <- "false"} else {x$genes$ControlType <- 0}}
    useControls <- 0  
  } else {
    useControls <- 1
  }
  
  
  #add information about datatype
  x$datatype <- datatype

  if(source == "genepix") { 
    warningMessages <- c(warningMessages, "[[WARNING]] With GenePix data, it is not possible to perform criterium 'low.pix' in upcoming Quality Control. \n")
  } else {
    if((is.null(x$other$rNumPix) & (x$datatype!="green")) | (is.null(other$gNumPix) & (x$datatype!="red"))) {
      warningMessages <- c(warningMessages, "[[WARNING]] No spot sizes given. Unable to perform criterium 'low.pix' in upcoming Quality Control.")
    }    
  }
  
  if((is.null(x$other$rMedianSignal) & (x$datatype!="green")) | (is.null(other$gMedianSignal) & (x$datatype!="red"))) {
    warningMessages <- c(warningMessages, "[[WARNING]] No Median Signal given. Unable to perform criterium 'mean.vs.median' in upcoming Quality Control. ")
  }    
  
  if((is.null(x$other$rIsWellAboveBG) & (x$datatype!="green")) | (is.null(x$other$gIsWellAboveBG) & (x$datatype!="red"))) {
    if((is.null(x$other$rBGSDev) & (x$datatype!="green")) | (is.null(x$other$gBGSDev) & (x$datatype!="red"))) {
      warningMessages <- c(warningMessages, "[[WARNING]] No pixel standard deviation given. Unable to perform criterium 'not.above.bg' in upcoming Quality Control.")
    }
  }
  if((is.null(x$other$rIsSaturated) & (x$datatype!="green")) | (is.null(x$other$gIsSaturated) & (x$datatype!="red"))) {
    if((is.null(x$other$rPercSatPix) & (x$datatype!="green")) | (is.null(x$other$gPercSatPix) & (x$datatype!="red"))) {
      warningMessages <- c(warningMessages, "[[WARNING]] No percentage of saturated pixels given. Unable to perform criterium 'saturated' in upcoming Quality Control.")
    }
  }    
  
  if(use.description) {
    if(x$datatype=="both") {
      colnames(x$R) <- colnames(x$G) <- colnames(x$Rb) <- colnames(x$Gb) <- targets$Description
    } else {
      colnames(x$E) <- colnames(x$Eb) <- targets$Description
    }
    if(!is.null(x$Rb.real)) colnames(x$Rb.real) <- targets$Description
    if(!is.null(x$Gb.real)) colnames(x$Gb.real) <- targets$Description
    for (n in names(x$other)) eval(parse("",-1,paste("colnames(x$other$",n,") <- targets$Description",sep="")))
  }

  ## Convert manualFlags (if present) to a standard format (1 means flagged, 0 means not flagged) 
  ## The standard Agilent FES IsManualFlag column does not need to be converted (this is in the standard boolean format)
  if(!is.null(x$other$manualFlags) & is.null(manual.flags)) {
    #agilent Flags will normally be boolean, just in case other non-zero values have been added, set to 1 as well
    if(source == "agilent") {
      x$other$manualFlags[x$other$manualFlags != 0] <- 1
    }
    #the standard genepix column will contain negative numbers for flagged spots
    #note that it is NOT recommended to use these as such, as spots will be flagged automatically as well using the same column
    #so we discard flagging in such case
    #this combination of conditions will never occur unless you modify the main script
    if(source == "genepix") {
      x$other$manualFlags <- NULL
      warningMessages <- c(warningMessages, "[[WARNING]] Manual flags column provided, but no manual flag values given (for Genepix data): no manual flags set")
    }
  } else {
    #when manual flags have been provided, consider these to be the values to be flagged
    #this means: when you want to add values for FES instead of replace, also provide the default 1 in the list of manual.flags
    if(!is.null(x$other$manualFlags)  & !is.null(manual.flags)) {
      bad.flagged <- matrix(data=FALSE, nrow=dim(x$other$manualFlags)[1], ncol=dim(x$other$manualFlags)[2])
      for (i in manual.flags) {
        if(i == "neg") {
          bad.flagged <- x$other$manualFlags<0 | bad.flagged
        } else if(i == "pos") {
          bad.flagged <- x$other$manualFlags>0 | bad.flagged
        } else {
          bad.flagged <- x$other$manualFlags==manual.flags[i] | bad.flagged
        }
      }
      x$other$manualFlags[] <- 0
      x$other$manualFlags[bad.flagged] <- 1
    }
  }

  ## New block-code
  cat(" * Checking block information")
  if(!is.null(blocks)) {
    temp.names <- colnames(x$genes)
    ## There are two situations possible: either block information is present, and the maxmium column and row number are low
    ## or the row and column numbering is increasing.
    a <- max(x$genes$Row)
    b <- max(x$genes$Col)

    if(sum(is.na(a)) > 0 | sum(is.na(b)) > 0) { stop("[[ ERROR ]] Row or column contains an NA value!") }

    if( dim(x$genes)[1] == a * b * blocks$total) {
      ## Some functionalities will not work with the current Row and Col setup. We will rewrite the Row and Col values
      #  in such a way that they follow each other up in a horizontal way.
      blocks$nspot.row <- a
      blocks$nspot.col <- b
      x$genes[,"Row_Original"] <- x$genes[,"Row"]
      x$genes[,"Col_Original"] <- x$genes[,"Col"]
      x$genes <- .adjustRowCols(x$genes, blocks=blocks)


      if(sum("Block" %in% temp.names) != 1) {
        x$genes <- .addBlockInfo(x$genes, blocks=blocks, plot.blocks=FALSE)
      }    


    } else {
      blocks$nspot.row <- max(x$genes$Row) / blocks$nblock.row
      blocks$nspot.col <- max(x$genes$Col) / blocks$nblock.col
      if(sum("Block" %in% temp.names) != 1) {
        x$genes <- .addBlockInfo(x$genes, blocks=blocks, plot.blocks=FALSE)
      }    
    }
    cat(paste("  ", max(x$genes$Block, na.rm=TRUE), " blocks found!\n", sep=""))

    x$printer <- list(ngrid.r=blocks$nblock.row, ngrid.c=blocks$nblock.col, 
                      nspot.r=blocks$nspot.row, nspot.c=blocks$nspot.col)
    cat(paste("  ", max(x$genes$Block, na.rm=TRUE), " blocks found!\n", sep=""))

    rm(temp.names)
  } else {
    x$printer <- list(ngrid.r=1, ngrid.c=1, nspot.r=max(x$genes$Row, na.rm=TRUE), nspot.c=max(x$genes$Col, na.rm=TRUE))
    cat("\n  None found\n")
  }
  ## Place targets into RG object
  x$targets <- targets
  
  ## Import spottypes (if provided) and add this information to the RGList
  if(!is.null(spottypes.file))
    x$genes$Status <- controlStatus(readSpotTypes(paste(data.path,spottypes.file,sep="")),x$genes)
 
  if(save.backup) {
    RG_Unprocessed <- x
    save(RG_Unprocessed, file="RG_Unprocessed.RData")
    rm(RG_Unprocessed)
  }
  # Remove all values on R/G/Rb/Gb that are associated with controls
  if(useControls == 1) {
    switch(source, "agilent" = { err.passed  <- (x$genes$ControlType[] != 0) },
                   "genepix" = { err.passed  <- (x$genes$ControlType[] != "false") },
                   "generic" = {
      if(is.null(controlType.value)) {
        stop("controlType.value not defined!")
      } else {
        temp <- matrix(data=NA, ncol=length(controlType.value), nrow=dim(x$genes$ControlType)[1])
        for(i in 1:length(controlType.value)) {
          temp[,i] <- x$genes$ControlType[] == controlType.value[i]
        }
        err.passed <- apply(temp, 1, function(x) { sum(x)[] >= 1 })
        rm(temp)
      }
                   })
    if(x$datatype=="both") {
      x$R[err.passed,]  <- NA
      x$G[err.passed,]  <- NA
      x$Rb[err.passed,] <- NA
      x$Gb[err.passed,] <- NA
    } else {
      x$E[err.passed,]  <- NA
      x$Eb[err.passed,] <- NA
    }
    if(!is.null(x$Rb.real))
      x$Rb.real[err.passed,] <- NA
    if(!is.null(x$Gb.real))
      x$Gb.real[err.passed,] <- NA
    if(!is.null(x$Eb.real))
      x$Eb.real[err.passed,] <- NA
    if(!is.null(x$other$rBGMedianSignal))
      x$other$rBGMedianSignal[err.passed,] <- NA
    if(!is.null(x$other$gBGMedianSignal))
      x$other$gBGMedianSignal[err.passed,] <- NA
    if(!is.null(x$other$rMedianSignal))
      x$other$rMedianSignal[err.passed,] <- NA
    if(!is.null(x$other$gMedianSignal))
      x$other$gMedianSignal[err.passed,] <- NA
  }
  if(!is.null(warningMessages)) {
    cat("--- arrayQC WARNINGS ---\n")
    cat(paste(warningMessages, collapse="\n"))
    cat("\n")
  }
  return(x)

}

############################
# removeManualFlaggedSpots #
############################

removeManualFlaggedSpots <- function(x) {
  if(is.null(x$other$manualFlags)) stop("the manualFlags field is empty, no flagging can be done")
  err.passed <- x$other$manualFlags == 1

  if(x$datatype=="both") {
    x$R[err.passed]  <- NA
    x$G[err.passed]  <- NA
    x$Rb[err.passed] <- NA
    x$Gb[err.passed] <- NA
  } else {
    x$E[err.passed]  <- NA
    x$Eb[err.passed] <- NA
  }

  if(!is.null(x$Rb.real))
    x$Rb.real[err.passed] <- NA
  if(!is.null(x$Gb.real))
    x$Gb.real[err.passed] <- NA
  if(!is.null(x$Eb.real))
    x$Eb.real[err.passed] <- NA
  if(!is.null(x$other$rBGMedianSignal))
    x$other$rBGMedianSignal[err.passed] <- NA
  if(!is.null(x$other$gBGMedianSignal))
    x$other$gBGMedianSignal[err.passed] <- NA
  if(!is.null(x$other$rMedianSignal))
    x$other$rMedianSignal[err.passed] <- NA
  if(!is.null(x$other$gMedianSignal))
    x$other$gMedianSignal[err.passed] <- NA

  return(x)
}

################################
##makeLimmaCompatible function##
################################

makeLimmaCompatible <- function(x) {
  
  if(x$source=="agilent") {
  
    cat("Executing makeLimmaCompatible script ...\n");
    
    #find out which feature numbers are missing from x
    feat <- x$genes$FeatureNum
    allfeat <- c(min(feat):max(feat))
    miss <- setdiff(allfeat, feat)
    
    # if miss has at least 1 value, correct for this, otherwise, don't do anything
    if(length(miss) >= 1) {
      #take part of the data set to use as mold (will be emptied)
      new <- x[1:length(miss),]
    
      #this is required to exist
      new$genes[] <- NA
    
      if(!is.null(new$genes$ProbeName ))
        new$genes$ProbeName <- "AgilentFE_deleted"
      if(!is.null(new$genes$GeneName ))
        new$genes$GeneName <- "AgilentFE_deleted"
      if(!is.null(new$genes$SystematicName ))
        new$genes$SystematicName <- "AgilentFE_deleted"
    
      #this is required to exist
      new$genes$ControlType <- 2
    
      if(x$datatype=="both") {
        new$R[] <- new$Rb[] <- NA
        new$G[] <- new$Gb[] <- NA
      } else {
        new$E[] <- new$Eb[] <- NA
      }
    
      if(!is.null(new$weights))
        new$weights[] <- 0
    
      for(field in names(new$other)) 
        eval(parse("",-1,paste("new$other$", field, "[] <- NA",sep="")))
    
      if(!is.null(new$Rb.real))
        new$Rb.real[] <- NA
      if(!is.null(new$Gb.real))
        new$Gb.real[] <- NA
      if(!is.null(new$Eb.real))
        new$Eb.real[] <- NA
    
      new$genes[,"FeatureNum"] <- miss
      
      # Add all the empty rows to the agilent object
      new <- rbind(x, new)
    
      # order all rows based on feature number
      new <- new[order(new$genes$FeatureNum),]
      rownames(new$genes) <- c(1:dim(new)[1])

      # Fill in the row and columns for the missing features, assuming that the featurenumbers are sorted by row first, then by column
      new$genes$Col <- rep(1:max(new$genes[,"Col"], na.rm=TRUE), max(new$genes[,"Row"], na.rm=TRUE))
      new$genes$Row <- sort(rep(1:max(new$genes[,"Row"], na.rm=TRUE), max(new$genes[,"Col"], na.rm=TRUE)))
      cat(paste("[INFO] Object dimensions altered. New dimensions:", dim(new$genes)[1], "(rows) -",dim(new$genes)[2], " (cols)\n"))
    } else {
      cat("[INFO] Object dimensions are OK, no alterations made!\n")
      new <- x
    }
    return(new)
  } else {
    
    cat("\n[[WARNING]] makeLimmaCompatible is only needed for Agilent Feature Extraction data sets.\nThe data object is returned as is.\n")
    
    return(x)
    
  }   
}

.addBlockInfo <- function(x, blocks=NULL, plot.blocks=FALSE) {
  error.block <- NULL
  temp.colnames <- colnames(x)
  if(sum( c("Row","Col") %in% temp.colnames) != 2) {
    if(sum("Row" %in% temp.colnames) == 0) { error.block <- c(error.block, "\n  - x$Row not present") }
    if(sum("Col" %in% temp.colnames) == 0) { error.block <- c(error.block, "\n  - x$Cow not present") }
  }
  if( sum( "Block" %in% temp.colnames ) == 1) { error.block <- c(error.block, "\n  - BLOCK information already present in main object!!") }
  if( max(x$Col, na.rm=TRUE) %% blocks$nblock.col != 0) { error.block <- c(error.block, "\n  - Incorrect nblock.col values!") }
  if( max(x$Row, na.rm=TRUE) %% blocks$nblock.row != 0) { error.block <- c(error.block, "\n  - Incorrect nblock.row values!") }

  if(!is.null(error.block)) { stop(" [ERROR] Incorrect values used for blocks. Please check the following variables:", error.block, "\n") }

  temp.block <- rep(NA, dim(x)[1])

  #-- Selecting blocks one by one.

  colFactor <- blocks$nspot.col
  rowFactor <- blocks$nspot.row
  blocknr <- 0

  for(i in 1:blocks$nblock.row) {
    for(j in 1:blocks$nblock.col) {
     blocknr <- blocknr + 1
     colStart <- ((j-1) * blocks$nspot.col) + 1
     colEnd   <- (j * blocks$nspot.col)
     rowStart <- ((i-1) * blocks$nspot.row) + 1
     rowEnd   <- (i * blocks$nspot.row)

     grab <- NULL
     for(k in 1:blocks$nspot.col) {
       grab$Cols <- c(grab$Cols, which(x$Col[] == c(colStart:colEnd)[k]))
     }
      
     for(l in 1:blocks$nspot.row) {
       grab$Rows <- c(grab$Rows, which(x$Row[] == c(rowStart:rowEnd)[l]))
     }
     selectedBlock <- sort(intersect(grab$Cols, grab$Rows))
     temp.block[selectedBlock] <- blocknr
     rm(grab, selectedBlock)
    }
  }
  if(plot.blocks == TRUE) {  
    matplot(x[,"Col"], -x[,"Row"], pch="*", ylab="Row #", xlab="Column #")
    for(i in 1:blocks$total) {
      matplot(x[which(temp.block[]==i),"Col"], -x[which(temp.block[]==i),"Row"], pch="*", add=T, col=rainbow(blocks$total)[i])  
    }
  }
  Block <- temp.block
  x <- cbind(x, Block)
  return(x)
}

.adjustRowCols <- function(y, blocks=NULL) {
  ##  For some dataformats the rows and column maximum values correspond to a block size.
  #   This may give further problems in the downstream workflow. To solve this we'll convert
  #   the rows and numbers to a more 'continuous' value, as if the blocks never existed.
  
  #   The direction is from left to right, so blocks$nblock.col will be used to set the maximum
  #   number of columns.
  #   In practice, for an 19x20 block (rowxblock format) the first block will retain it's values,
  #   however the values for the 2nd block will lie between 1-19 (row) and 21-40 (col) untill the maximum
  #   number of columns have been found.


  error <- NULL
  if(is.null(blocks)) { error <- c(error, "blocks is NULL") }
  if( sum(c("Col", "Row") %in% colnames(y)) != 2) { error <- c(error, "Col or Row not present in your object!") }

  block.size <- blocks$nspot.row * blocks$nspot.col

  row.counts <- 1
  col.counts <- 1

  ## Counter information:
  #  First counter goes through each block (length: blocksize)
  #  Second counter goes through all rows
  #  Third counter goes through all columns
  row.counter <- 1
  col.counter <- 1

  for( ii in 1:blocks$total ) {
    Block.end <- ((ii - 1) * block.size) + block.size
    Block.start <- ((ii - 1) * block.size) + 1
    selected <- c(Block.start: Block.end)
#    cat(paste("Block Row:", row.counter, " - Block Col:", col.counter, "\n"))

    if(row.counter != 1) {
      y[selected,"Row"] <- y[selected,"Row"] + ( blocks$nspot.row * (row.counter - 1) )
    }
    if( col.counter != 1) {
      y[selected,"Col"] <- y[selected,"Col"] + ( blocks$nspot.col * (col.counter - 1) )
    }

    if( (col.counter + 1) <= blocks$nblock.col ) {
      col.counter <- col.counter + 1
    } else {
      col.counter <- 1
      if( (row.counter + 1) <= blocks$nblock.row ) {
        row.counter <- row.counter + 1
      } else {
        row.counter <- 1
      }
    }


    rm(Block.end, Block.start, selected)
  }
  return(y)
}


#######################################
##     END OF READFILES SCRIPT       ##
#######################################
