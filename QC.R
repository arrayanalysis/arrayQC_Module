###############################
##  QualityControl Function  ##
###############################

QualityControl <- function (x=NULL, num.pix=NULL, low.pix=TRUE, mean.vs.median=TRUE, not.above.bg=TRUE, saturated=TRUE, criteria.norm=NULL, criteria=NULL) {

  if (is.null(x))
    stop("no data object specified")  
  if (sum(low.pix+mean.vs.median+not.above.bg+saturated) != length(criteria.norm))
    stop("the number of QC tests does not correspond to the number of weight criteria for normalization")
  if (sum(low.pix+mean.vs.median+not.above.bg+saturated) != length(criteria))
    stop("the number of QC tests does not correspond to the number of weight criteria for further processing")
    
  cat("-----------------------------------\n")
  cat("   Q U A L I T Y   C O N T R O L\n") 
  cat("-----------------------------------\n") 
  if(low.pix & x$source == "genepix") {stop("Unable to perform low.pix criterium for GenePix data. Please set low.pix to FALSE. \n")}
  if(low.pix & ((is.null(x$other$rNumPix) & (x$datatype!="green")) | (is.null(x$other$gNumPix) & (x$datatype!="red")))) {stop("No Number-of-pixels fields found in data-object. Please make sure that columns that contain red and green number of pixels are read from the datafiles. \n")}
  if(mean.vs.median & ((is.null(x$other$rMedianSignal) & (x$datatype!="green")) | (is.null(x$other$gMedianSignal) & (x$datatype!="red")))) {stop("No MedianSignal fields found in data-object. Please make sure that columns that contain red and green MedianSignal are read from the datafiles. \n")}
  if(not.above.bg & ((is.null(x$other$rIsWellAboveBG) & (x$datatype!="green")) | (is.null(x$other$gIsWellAboveBG) & (x$datatype!="red")))) {
    if((is.null(x$other$rBGSDev) & (x$datatype!="green")) | (is.null(x$other$gBGSDev) & (x$datatype!="red"))) {
      stop("Unable to perform quality control for Well Above Background. Please make sure that columns WellAboveBG (in case of Feature Extraction Software) or BGPixSDev (in case of GenePix) are read from the datafiles. \n")
    }
  }  
  if(saturated & ((is.null(x$other$rIsSaturated) & (x$datatype!="green")) | (is.null(x$other$gIsSaturated) & (x$datatype!="red")))) {
    if((is.null(x$other$rPercSatPix) & (x$datatype!="green")) | (is.null(x$other$gPercSatPix) & (x$datatype!="red"))) {
      stop("Unable to perform quality control for Saturation. Please make sure that columns IsSaturated (in case of Feature Extraction Software) or F[wavelength]%Sat. (in case of GenePix) are read from the datafiles. \n")
    }
  }  
  
  cat(" The following QC tests will be performed:\n")
  if(low.pix) {
    cat(" * Number of Pixels criteria")
    max.pix <- max(c(x$other$gNumPix, x$other$rNumPix), na.rm=TRUE)
    if(is.null(num.pix)) {
      num.pix <- ceiling(max.pix * 0.75)
    }
  cat(paste(" (Max pixels: ",  max.pix, " Cut-off: ", num.pix, ")\n", sep=""))
  }
  
  if(mean.vs.median) {cat(" * Mean vs Median ratio > 0.9\n")}
  if(not.above.bg) {cat(" * Signal Well Above Background\n")}
  if(saturated) {cat(" * Spot saturation\n")}
  cat("\n\n")
  
  #reset fields in case the function had been run before
  x$other$LowNumPix <- x$other$LowMeanMedianRatio <- x$other$NOTWellAboveBG <- x$other$Saturated <- NULL
  #also reset summary field, as it will not be adapted automatically in THIS function
  x$summary <- NULL
  
  QC.vector <- NULL
  tests <- rep(1,4)
  if(low.pix) {QC.vector <- c(QC.vector, "LowNumPix")} else {tests[1]  <- 0}
  if(mean.vs.median) {QC.vector <- c(QC.vector, "LowMeanMedianRatio")} else {tests[2] <- 0}
  if(not.above.bg) {QC.vector <- c(QC.vector, "NOTWellAboveBG")} else {tests[3] <- 0}
  if(saturated) {QC.vector <- c(QC.vector, "Saturated")} else {tests[4] <- 0}
  
  if(x$datatype == "green" & sum(criteria == -5) > 0) {
    cat("WARNING: -5 used in criteria while having Cy3 array, all -5 values set to -3\n")
  criteria[criteria == -5] <- -3
  }
  if(x$datatype == "red" & sum(criteria == -3) > 0) {
    cat("WARNING: -3 used in criteria while having Cy5 array, all -3 values set to -5\n")
  criteria[criteria == -3] <- -5
  }
  if(x$datatype == "green" & sum(criteria.norm == -5) > 0) {
    cat("WARNING: -5 used in criteria.norm while having Cy3 array, all -5 values set to -3\n")
  criteria.norm[criteria.norm == -5] <- -3
  }
  if(x$datatype == "red" & sum(criteria.norm == -3) > 0) {
    cat("WARNING: -3 used in criteria.norm while having Cy5 array, all -3 values set to -5\n")
  criteria.norm[criteria.norm == -3] <- -5
  }
  
  x$QC.vector <- QC.vector
  x$criteria.norm <- criteria.norm
  x$criteria <- criteria
  
  # Calculating IsWellAboveBG 0 or 1 if standard deviation columns are present 
  if(not.above.bg & ((is.null(x$other$rIsWellAboveBG) & (x$datatype!="green")) | (is.null(x$other$gIsWellAboveBG) & (x$datatype!="red")))) {
    if(!is.null(x$other$rBGSDev) & !is.null(x$other$gBGSDev) & (x$datatype=="both")) {
      rBGLimit <- x$Rb+x$other$rBGSDev*2.6
      gBGLimit <- x$Gb+x$other$gBGSDev*2.6
      x$other$rIsWellAboveBG <- (x$R>rBGLimit)+0
      x$other$gIsWellAboveBG <- (x$G>gBGLimit)+0  
    }
    if(!is.null(x$other$rBGSDev) & (x$datatype=="red")) {
      rBGLimit <- x$Eb+x$other$rBGSDev*2.6
      x$other$rIsWellAboveBG <- (x$E>rBGLimit)+0
    }
    if(!is.null(x$other$gBGSDev) & (x$datatype=="green")) {
      gBGLimit <- x$Eb+x$other$gBGSDev*2.6
      x$other$gIsWellAboveBG <- (x$E>gBGLimit)+0  
    }
  }  
  
  # Calculating Saturated
  if(saturated & ((is.null(x$other$rIsSaturated) & (x$datatype!="green")) | (is.null(x$other$gIsSaturated) & (x$datatype!="red")))) {
    if(!is.null(x$other$rPercSatPix) & (x$datatype!="green")) {
      x$other$rIsSaturated <- (x$other$rPercSatPix>50)+0
    }
    if(!is.null(x$other$gPercSatPix) & (x$datatype!="red")) {
      x$other$gIsSaturated <- (x$other$gPercSatPix>50)+0
    }
  }
  
  x <- QCFields(x, num.pix=num.pix, tests=tests)
   
  tests <- c("Number of Pixels","Mean vs Median","Well Above Background","Saturation")
  
  cat("Weights for normalization are created based on the following criteria:\n")
  cat(paste(tests[c(low.pix,mean.vs.median,not.above.bg,saturated)],  collapse=", "), ":", paste(criteria.norm, collapse=","), "\n")
  x$other$weights.norm <- CreateWeights(x, norm=TRUE)
  
  cat("Weights for further processing are created based on the following criteria:\n")
  cat(paste(tests[c(low.pix,mean.vs.median,not.above.bg,saturated)], collapse=", "),":", paste(criteria, collapse=","), "\n")
  x$weights <- CreateWeights(x, norm=FALSE)
  
  #create weights.analysis field for possible later use, copy structure from weights field and overwrite with NAs
  x$other$weights.analysis <- x$weights
  x$other$weights.analysis[] <- NA
  
  return(x)
  }

###########################
##   QCFields function   ##
########################### 

QCFields <- function(x, num.pix=NULL, tests=NULL) {
  
  if((class(x)!="RGList") && (class(x)!="EListRaw")) stop("Given argument is NOT of RGList or EListRaw class!")

  Control <- array(dim = c(length(x$QC.vector),dim(x)[2],dim(x)[1]));
  dimnames(Control) <- list(x$QC.vector, colnames(x), x$genes$FeatureNum);
  
  cat(paste("* Performing gene-specific QC on", dim(x)[1], "reporters (for", dim(x)[2], "microarrays). \n "))
  
  current_index <- 1
  
  if(x$datatype=="both") {
    signalIsNA <- is.na(x$R)&is.na(x$G)
  } else {
    signalIsNA <- is.na(x$E)
  }
  
  if(tests[1]==1){
  eval(parse("",-1,paste("x$other$",x$QC.vector[current_index],"<-",ifelse(x$datatype!="red","((x$other$gNumPix<num.pix)*-3)+",""),ifelse(x$datatype!="green","((x$other$rNumPix<num.pix)*-5)+",""),"ifelse(signalIsNA,NA,0)",sep="")))
  current_index<-current_index+1
  }
  if(tests[2]==1){
    mm.cutoff.r <- mm.cutoff.g <- 0.9
    if(RG$datatype!="fes") {
      if(RG$datatype=="both") {
        mmratioG <- ((RG$other$gMedianSignal/RG$G)*(RG$other$gMedianSignal<=RG$G))+((RG$G/RG$other$gMedianSignal)*(RG$other$gMedianSignal>RG$G))
        mmratioG <- sort(mmratioG[!is.na(mmratioG)])
        mm.cutoff.g <- min(0.9,max(0.75,mmratioG[0.05*(length(mmratioG))]))
        mmratioR <- ((RG$other$rMedianSignal/RG$R)*(RG$other$rMedianSignal<=RG$R))+((RG$R/RG$other$rMedianSignal)*(RG$other$rMedianSignal>RG$R))
        mmratioR <- sort(mmratioR[!is.na(mmratioR)])
        mm.cutoff.r <- min(0.9,max(0.75,mmratioR[0.05*(length(mmratioR))]))
      }
      if(RG$datatype=="green") {
        mmratioG <- ((RG$other$gMedianSignal/RG$E)*(RG$other$gMedianSignal<=RG$E))+((RG$E/RG$other$gMedianSignal)*(RG$other$gMedianSignal>RG$E))
        mmratioG <- sort(mmratioG[!is.na(mmratioG)])
        mm.cutoff.g <- min(0.9,max(0.75,mmratioG[0.05*(length(mmratioG))]))
      }
      if(RG$datatype=="red") {
        mmratioR <- ((RG$other$rMedianSignal/RG$E)*(RG$other$rMedianSignal<=RG$E))+((RG$E/RG$other$rMedianSignal)*(RG$other$rMedianSignal>RG$E))
        mmratioR <- sort(mmratioR[!is.na(mmratioR)])
        mm.cutoff.r <- min(0.9,max(0.75,mmratioR[0.05*(length(mmratioR))]))
      }
    }
    eval(parse("",-1,paste("x$other$",x$QC.vector[current_index],"<-",ifelse(x$datatype=="both","((((x$other$gMedianSignal/x$G)<mm.cutoff.g)|((x$G/x$other$gMedianSignal)<mm.cutoff.g))*-3)+",""),
                                                                      ifelse(x$datatype=="both","((((x$other$rMedianSignal/x$R)<mm.cutoff.r)|((x$R/x$other$rMedianSignal)<mm.cutoff.r))*-5)+",""),
                                                                      ifelse(x$datatype=="green","((((x$other$gMedianSignal/x$E)<mm.cutoff.g)|((x$E/x$other$gMedianSignal)<mm.cutoff.g))*-3)+",""),
                                                                      ifelse(x$datatype=="red","((((x$other$rMedianSignal/x$E)<mm.cutoff.r)|((x$E/x$other$rMedianSignal)<mm.cutoff.r))*-5)+",""),"0",sep="")))
    current_index<-current_index+1
  }
  if(tests[3]==1){
    eval(parse("",-1,paste("x$other$",x$QC.vector[current_index],"<-",ifelse(x$datatype!="red","((x$other$gIsWellAboveBG==0)*-3)+",""),ifelse(x$datatype!="green","((x$other$rIsWellAboveBG==0)*-5)+",""),"ifelse(signalIsNA,NA,0)",sep="")))
    current_index<-current_index+1
  }
  if(tests[4]==1){
    eval(parse("",-1,paste("x$other$",x$QC.vector[current_index],"<-",ifelse(x$datatype!="red","((x$other$gIsSaturated==1)*-3)+",""),ifelse(x$datatype!="green","((x$other$rIsSaturated==1)*-5)+",""),"ifelse(signalIsNA,NA,0)",sep="")))
    current_index<-current_index+1
  }
    
  cat(" ok.\n\n")
  return(x)
}

##############################
##  CreateWeights function  ##
##############################

CreateWeights <- function(x, norm=FALSE) {

  if(norm) {criteria <- x$criteria.norm} else {criteria <- x$criteria}
  
  if(length(criteria) != length(x$QC.vector))
    stop("Criteria and QC.vector are not of equal size")
  
  cat("Creating Weights") 
  
  weight <- matrix(data=NA, ncol=dim(x)[2], nrow=dim(x)[1])
  colnames(weight) <- dimnames(x)[[2]]
  total <- matrix(data=0, ncol=dim(x)[2], nrow=dim(x)[1])
  
  for(j in 1:length(x$QC.vector))
    total <- total+eval(parse("",-1,paste("x$other$",x$QC.vector[j],">criteria[j]",sep="")))
  weight <- (total==length(criteria))+0
  
  if(!is.null(x$other$manualFlags)) weight[x$other$manualFlags == 1] <- 0
  
  weight[is.na(weight)] <- 0
 
  cat("\n\n")
  return(weight)
}

###########################
##  AddSummary function  ##
###########################

# AddSummary summarizes all information available on each microarray, such as foreground signal, measured background signal,
# estimated background signal (Spatial Detrending), etc.
AddSummary <- function(x) {

  if((class(x)!="RGList") && (class(x)!="EListRaw")) stop("Given argument is NOT of RGList or EListRaw class!")

  cat("#- SUMMARY GENERATOR FOR DATA OBJECT -#\n\n")
    
  summary <- NULL
  colnames.summary <- NULL
  if(x$datatype=="both") {
    rownames.summary <- colnames(x$R)
  } else {
    rownames.summary <- colnames(x$E)
  }
  rownames.summary <- c(rownames.summary,"Minimum","Mean","Median","Maximum")

  if(x$source=="genepix") {
    isControl <- x$genes$ControlType!="false"
  } else {
    isControl <- x$genes$ControlType!=0
  }
  if(x$source=="agilent") {
    est <- " Est. "
  } else {
    est <- " "
  }
  
  if(x$datatype=="both") {
    cat("* Creating Summary for Red Foreground/Background signals")
    Rf <- apply(x$R, 2, function(x) { mean(x, na.rm=TRUE) })
    Rf <- c(Rf,min(Rf),mean(Rf),median(Rf),max(Rf))
    Rb <- apply(x$Rb, 2, function(x) { mean(x, na.rm=TRUE) })
    Rb <- c(Rb,min(Rb),mean(Rb),median(Rb),max(Rb))
    summary <- cbind(Rf,Rb)
    colnames.summary <- c("Red Foreground",paste("Red",est,"Background",sep=""))
    if(!is.null(x$other$rBGMeanSignal)) {
      Rb.real <- apply(x$other$rBGMeanSignal, 2, function(x) mean(x, na.rm=TRUE))
      Rb.real <- c(Rb.real,min(Rb.real),mean(Rb.real),median(Rb.real),max(Rb.real))
      summary <- cbind(summary,Rb.real)
      colnames.summary <- c(colnames.summary,"Red Real Background")
    }
    cat(" - ok\n")
    cat("* Creating Summary for Green Foreground/Background signals")
    Gf <- apply(x$G, 2, function(x) { mean(x, na.rm=TRUE) })
    Gf <- c(Gf,min(Gf),mean(Gf),median(Gf),max(Gf))
    Gb <- apply(x$Gb, 2, function(x) { mean(x, na.rm=TRUE) })
    Gb <- c(Gb,min(Gb),mean(Gb),median(Gb),max(Gb))
    summary<-cbind(summary,Gf,Gb)
    colnames.summary <- c(colnames.summary,"Green Foreground",paste("Green",est,"Background",sep=""))
    if(!is.null(x$other$gBGMeanSignal)) {
      Gb.real <- apply(x$other$gBGMeanSignal, 2, function(x) mean(x, na.rm=TRUE))
      Gb.real <- c(Gb.real,min(Gb.real),mean(Gb.real),median(Gb.real),max(Gb.real))
      summary <- cbind(summary,Gb.real)
      colnames.summary <- c(colnames.summary,"Green Real Background")
    }
    cat(" - ok\n")
  }

  if(x$datatype=="red") {
    cat("* Creating Summary for Red Foreground/Background signals")
    Rf <- apply(x$E, 2, function(x) { mean(x, na.rm=TRUE) })
    Rf <- c(Rf,min(Rf),mean(Rf),median(Rf),max(Rf))
    Rb <- apply(x$Eb, 2, function(x) { mean(x, na.rm=TRUE) })
    Rb <- c(Rb,min(Rb),mean(Rb),median(Rb),max(Rb))
    summary <- cbind(Rf,Rb)
    colnames.summary <- c("Red Foreground",paste("Red",est,"Background",sep=""))
    if(!is.null(x$other$rBGMeanSignal)) {
      Rb.real <- apply(x$other$rBGMeanSignal, 2, function(x) mean(x, na.rm=TRUE))
      Rb.real <- c(Rb.real,min(Rb.real),mean(Rb.real),median(Rb.real),max(Rb.real))
      summary <- cbind(summary,Rb.real)
      colnames.summary <- c(colnames.summary,"Red Real Background")
    }
    cat(" - ok\n")
  }

  if(x$datatype=="green") {
    cat("* Creating Summary for Green Foreground/Background signals")
    Gf <- apply(x$E, 2, function(x) { mean(x, na.rm=TRUE) })
    Gf <- c(Gf,min(Gf),mean(Gf),median(Gf),max(Gf))
    Gb <- apply(x$Eb, 2, function(x) { mean(x, na.rm=TRUE) })
    Gb <- c(Gb,min(Gb),mean(Gb),median(Gb),max(Gb))
    summary<-cbind(summary,Gf,Gb)
    colnames.summary <- c(colnames.summary,"Green Foreground",paste("Green",est,"Background",sep=""))
    if(!is.null(x$other$gBGMeanSignal)) {
      Gb.real <- apply(x$other$gBGMeanSignal, 2, function(x) mean(x, na.rm=TRUE))
      Gb.real <- c(Gb.real,min(Gb.real),mean(Gb.real),median(Gb.real),max(Gb.real))
      summary <- cbind(summary,Gb.real)
      colnames.summary <- c(colnames.summary,"Green Real Background")
    }
    cat(" - ok\n")
  }
  
  if(!is.null(x$other$weights.norm)) {
    cat("* Creating Summary for High Quality Spots (Normalization step)")
    weight.values.norm <- apply(x$other$weights.norm[!isControl,], 2, function(x) { sum(x,na.rm=TRUE)} )
    weight.values.norm <- c(weight.values.norm,min(weight.values.norm),round(mean(weight.values.norm),0),round(median(weight.values.norm),0),max(weight.values.norm))
    summary <- cbind(summary,weight.values.norm,round((weight.values.norm/sum(!isControl,na.rm=TRUE))*100,2))
    colnames.summary <- c(colnames.summary,"Number of Good Probes (Norm)","Percentage Good Probes (Norm)")
    cat(" - ok\n")
  }
        
  if(!is.null(x$weights)) {
    cat("* Creating Summary for High Quality Spots (Further processing)")
    weight.values <- apply(x$weights[!isControl,], 2, function(x) { sum(x,na.rm=TRUE)} )
    weight.values <- c(weight.values,min(weight.values),round(mean(weight.values),0),round(median(weight.values),0),max(weight.values))
    summary <- cbind(summary,weight.values,round((weight.values/sum(!isControl,na.rm=TRUE))*100,2))
    colnames.summary <- c(colnames.summary,"Number of Good Probes (Further processing)","Percentage Good Probes (Further processing)" )
    cat(" - ok\n")
  }
  
  if(!is.null(x$QC.vector)) {
    cat("* Creating Criteria-specific Summary for High Quality Spots (Normalization step)")
    spec.values.norm <- NULL
    for (i in 1:length(x$QC.vector)) {
      eval(parse("",-1,paste(x$QC.vector[i],".temp <- x$other$",x$QC.vector[i],"; ",x$QC.vector[i],".temp[isControl", ifelse(!is.null(x$other$manualFlags), " | x$other$manualFlags==1", ""),"] <- NA",sep="")))
      spec.values.norm.col <- apply(eval(parse("",-1,paste(x$QC.vector[i],".temp",sep=""))), 2, function(y) { sum(y<=x$criteria.norm[i],na.rm=TRUE)} )
      spec.values.norm.col <- c(spec.values.norm.col,min(spec.values.norm.col),round(mean(spec.values.norm.col),0),round(median(spec.values.norm.col),0),max(spec.values.norm.col))
      spec.values.norm <- cbind(spec.values.norm,spec.values.norm.col)
    }
    summary <- cbind(summary,spec.values.norm)
    for (n in x$QC.vector) { colnames.summary <- c(colnames.summary,paste(n," (Norm)", sep="")) }
    cat(" - ok\n")
  }
              
  if(!is.null(x$QC.vector)) {
    cat("* Creating Criteria-specific Summary for High Quality Spots (Further processing)")
    spec.values <- NULL
    for (i in 1:length(x$QC.vector)) {
      eval(parse("",-1,paste(x$QC.vector[i],".temp <- x$other$",x$QC.vector[i],"; ",x$QC.vector[i],".temp[isControl", ifelse(!is.null(x$other$manualFlags), " | x$other$manualFlags==1", ""),"] <- NA",sep="")))
      spec.values.col <- apply(eval(parse("",-1,paste(x$QC.vector[i],".temp",sep=""))), 2, function(y) { sum(y<=x$criteria[i],na.rm=TRUE)} )
      spec.values.col <- c(spec.values.col,min(spec.values.col),round(mean(spec.values.col),0),round(median(spec.values.col),0),max(spec.values.col))
      spec.values <- cbind(spec.values,spec.values.col)
    }
    summary <- cbind(summary,spec.values)
    for (n in x$QC.vector) { colnames.summary <- c(colnames.summary,paste(n," (Proc)", sep="")) }
    cat(" - ok\n")
  }
  
  if(!is.null(x$other$manualFlags)) {
    cat("* Creating Summary for Manually Flagged Spots")
    flagged.values <- apply(x$other$manualFlags[!isControl,], 2, function(x) { sum(x,na.rm=TRUE)} )
    flagged.values <- c(flagged.values,min(flagged.values),round(mean(flagged.values),0),round(median(flagged.values),0),max(flagged.values))
    summary <- cbind(summary,flagged.values)
    colnames.summary <- c(colnames.summary,"Number of Manually Flagged Probes")
    cat(" - ok\n")
  }
      
  colnames(summary) <- colnames.summary
  rownames(summary) <- rownames.summary
  
  x$summary <- summary

  cat("[- Ready -]\n\n")
  return(x)
}

########################
## SaveData function  ##
########################

SaveData <- function(x, save.RG=TRUE, save.MA=TRUE, save.QC=TRUE, detailed.QC=TRUE, raw=NULL, store.output=FALSE) {
  if(is.null(x)) stop("data object (MAList or matrix) must be provided")
  if(class(x)!="MAList" & class(x)!="matrix") stop("data object must be either of class MAList or of class matrix")
  if(class(x)=="matrix" & is.null(raw)) stop("when saving a matrix object, also provide raw data object (i.e. RG) for addition of gene information")
  if((save.QC | detailed.QC) & is.null(raw)) stop("when saving QC information, also provide raw data object for adding this information")
  if(!is.null(raw) & (class(raw)!="RGList") & class(raw)!="EListRaw") stop("the raw data object should be of class RGList or EListRaw")
  #if raw is not given (e.g. x is of class MAList and no QC needed to be saved), copy x into raw to get uniform code for both cases
  if(is.null(raw)) raw <- x
  
  outputData <- NULL
  a <- list()

  if(!is.null(raw$genes$FeatureNum)) outputData <- cbind(outputData, FeatureNum=raw$genes$FeatureNum)
  if(!is.null(raw$genes$ProbeName)) outputData <- cbind(outputData, ProbeName=raw$genes$ProbeName)
  if(!is.null(raw$genes$ControlType)) outputData <- cbind(outputData, ControlType=raw$genes$ControlType)
  if(!is.null(raw$genes$GeneName)) outputData <- cbind(outputData, GeneName=raw$genes$GeneName)
  if(!is.null(raw$genes$Description)) outputData <- cbind(outputData, Description=raw$genes$Description)
  if(!is.null(raw$genes$Status)) outputData <- cbind(outputData, Status=raw$genes$Status)
  
  if(save.RG) {
    if(class(x)=="matrix") {
      cat("---One-channel data, only estimates for the", raw$datatype,"channel will be saved---\n")
      ESTval <- x
      colnames(ESTval) <- paste(colnames(ESTval),ifelse(raw$datatype=="green","G","R"))
      outputData <- cbind(outputData, ESTval)
      rm(ESTval)
    } else { #MAList
      if (raw$datatype=="both") {
        RGval <- cbind(x$A+0.5*x$M, x$A-0.5*x$M)
        colnames(RGval) <- paste(colnames(RGval),rep(c("R","G"),rep(dim(x$A)[2],2)))
        RGval <- RGval[,rep(1:dim(x$A)[2],rep(2,dim(x$A)[2]))+c(0,dim(x$A)[2])]
        outputData <- cbind(outputData, RGval)
        rm(RGval)
      } else { #single channel
        cat("---One-channel data, only estimates for the", raw$datatype,"channel will be saved---\n")
        ESTval <- x$other$EST
        colnames(ESTval) <- paste(colnames(ESTval),ifelse(raw$datatype=="green","G","R"))
        outputData <- cbind(outputData, ESTval)
        a[["RG"]] <- outputData
        rm(ESTval)
      }
    }
  }

  if(save.MA) {
    if(class(x)=="matrix") {
      cat("WARNING: one-channel data, no MA values can be saved\n")
      outputData <- cbind(outputData, x)
    } else  { #MAList
      if (raw$datatype=="both") {
        MAval <- cbind(x$M, x$A)
        colnames(MAval) <- paste(colnames(MAval),rep(c("M","A"),rep(dim(x$A)[2],2)))
        outputData <- cbind(outputData, MAval)
        rm(MAval)
      } else { #single channel
        cat("WARNING: one-channel data, no MA values can be saved\n")
      }
    }
  }
  
  if(save.QC) {
    for(v in raw$QC.vector) {
      if(raw$datatype!="both") {
        eval(parse("",-1,paste("outputData <- cbind(outputData, \"num ", v, ifelse(raw$datatype=="green"," G"," R"), "\"=apply(raw$other$",v,",1,function(y) sum(y!=0)))",sep="")))
      } else { #save flag totals per channel
        eval(parse("",-1,paste("outputData <- cbind(outputData, \"num ", v, " R\"=apply(raw$other$",v,",1,function(y) sum(y<=-5)), \"num ", v, " G\"=apply(raw$other$",v,",1,function(y) sum((y==-3) | (y==-8))))",sep="")))
      }
    }
    if(!is.null(raw$other$manualFlags)) {
      if(sum(!is.na(raw$other$manualFlags) & raw$other$manualFlags!=0)>0) {
        outputData <- cbind(outputData, "num manualFlags"=apply(raw$other$manualFlags,1,sum))
      }
    }
    outputData <- cbind(outputData, "num flagged norm"=apply(1-raw$other$weights.norm,1,sum))
    outputData <- cbind(outputData, "num flagged"=apply(1-raw$weights,1,sum))
    if(sum(!is.na(raw$other$weights.analysis))>0) {
      outputData <- cbind(outputData, "num flagged analysis"=apply(1-raw$other$weights.analysis,1,sum))
    }
  }
  
  if(detailed.QC) {  #note that in this case, raw will neccesarily be an RGList
    outputQC <- NULL
    if(!is.null(raw$genes$FeatureNum)) outputQC <- cbind(outputQC,FeatureNum=raw$genes$FeatureNum)
    if(!is.null(raw$genes$ProbeName)) outputQC <- cbind(outputQC,ProbeName=raw$genes$ProbeName)
    if(!is.null(raw$genes$ControlType)) outputQC <- cbind(outputQC,ControlType=raw$genes$ControlType)
    if(!is.null(raw$genes$GeneName)) outputQC <- cbind(outputQC,GeneName=raw$genes$GeneName)
    if(!is.null(raw$genes$Description)) outputQC <- cbind(outputQC,Description=raw$genes$Description)
    if(!is.null(raw$genes$Status)) outputQC <- cbind(outputQC, Status=raw$genes$Status)
 
    for(v in raw$QC.vector) {
      if(raw$datatype!="both") {
        eval(parse("",-1,paste("QCval <- (raw$other$",v,"<0)+0", sep="")))
        colnames(QCval) <- paste(colnames(QCval),v,ifelse(raw$datatype=="green","G","R"))
        outputQC <- cbind(outputQC, QCval)
        rm(QCval)
      } else {#save flags per channel
        eval(parse("",-1,paste("QCval <- cbind((raw$other$",v,"<=-5)+0, ((raw$other$",v,"==-3) | (raw$other$",v,"==-8))+0)", sep="")))
        colnames(QCval) <- paste(colnames(QCval),v,rep(c("R","G"),rep(dim(raw)[2],2)))
        QCval <- QCval[,rep(1:dim(raw)[2],rep(2,dim(raw)[2]))+c(0,dim(raw)[2])]
        outputQC <- cbind(outputQC, QCval)
        rm(QCval)
      }
    }
    if(!is.null(raw$other$manualFlags)) {
      if(sum(!is.na(raw$other$manualFlags) & raw$other$manualFlags!=0)>0) {
        outputQC <- cbind(outputQC, raw$other$manualFlags)
        colnames(outputQC)[(length(colnames(outputQC))-dim(raw$other$manualFlags)[2]+1):(length(colnames(outputQC)))] <- paste(colnames(outputQC)[(length(colnames(outputQC))-dim(raw$other$manualFlags)[2]+1):(length(colnames(outputQC)))],"manualFlags")
      }
    }
    outputQC <- cbind(outputQC, raw$other$weights.norm)
    colnames(outputQC)[(length(colnames(outputQC))-dim(raw$other$weights.norm)[2]+1):(length(colnames(outputQC)))] <- paste(colnames(outputQC)[(length(colnames(outputQC))-dim(raw$other$weights.norm)[2]+1):(length(colnames(outputQC)))],"flagged norm")
    outputQC <- cbind(outputQC, raw$weights)
    colnames(outputQC)[(length(colnames(outputQC))-dim(raw$weights)[2]+1):(length(colnames(outputQC)))] <- paste(colnames(outputQC)[(length(colnames(outputQC))-dim(raw$weights)[2]+1):(length(colnames(outputQC)))],"flagged")
    if(sum(!is.na(raw$other$weights.analysis))>0) {
      outputQC <- cbind(outputQC, raw$other$weights.analysis)
      colnames(outputQC)[(length(colnames(outputQC))-dim(raw$other$weights.analysis)[2]+1):(length(colnames(outputQC)))] <- paste(colnames(outputQC)[(length(colnames(outputQC))-dim(raw$other$weights.analysis)[2]+1):(length(colnames(outputQC)))],"flagged analysis")
    }
  }

  write.table(outputData, file="arrayQC_normData.table", quote=FALSE, sep="\t", row.names=FALSE)
  if(detailed.QC) { write.table(outputQC, file="arayQC_flags.table", quote=FALSE, sep="\t", row.names=FALSE) }
  
  if(store.output == 1) {
    a <- list()
    a[["output"]] <-  outputData
    if(detailed.QC) { a[["QC"]] <- outputQC }
    return(a)
  }
}

########################
##  END OF QC SCRIPT  ##
########################
