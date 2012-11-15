## One and two channel QC - Dept. of Bioinformatics / Toxicogenomics - Maastricht University - the Netherlands
# Version: 1.3.0-DEV
# Last adjustment: 13-Nov-2012

###########################
## USER-DEFINED SETTINGS ##
###########################

## Enter the directories to be used, use forward slashes (/) or double backward slashed (\\) as separator
#  and end with a forward slash
datapath  <- ""
scriptpath <- "/home/stan/Desktop/SVN/r-packages/arrayQC/dev/"
# scriptpath <- "http://svn.bigcat.unimaas.nl/r-packages/arrayQC/dev/"


## Select the data format of the microarray data files
# "fes" (Feature Extraction Software), "genepix" or "generic" are valid here
dataformat <- ""

## Select the type of data set
# "two-channel", "red" or "green" are valid options
datatype <- ""

## If the microarray contains more than one block, please change the values below accordingly.
# 
number.of.blocks <- NULL
nblock.row <- NULL
nblock.col <- NULL


## columnHeader file. Change this value to the name (and location) of the file containing the 
#  column headers that will be read in. For standard Agilent and GenePix runs this variable
#  should be empty.
columnHeaderFile <- NULL

## Give one or more values - if any - that you used as manual flagging values within the data set
#  examples:
#  NULL       (default, no flagging for Genepix, all non-zero values are flags for Agilent FES)
#  -1         (-1 as a value)
#  c(-1,-10)  (-1 and -10 as values, the c() is used to group the values)
#  You can also use the words "neg" and "pos" for all negative or all positive values
#  example:   c("neg",10)
useAsManualFlags <- NULL

## For generic formats: please supply the controlType variable below with the value(s)
#  corresponding with the positive and negative controls used in the ControlType column.
controlType.value <- NULL

## Variables that are set to run the script. Only when arrayQC.mode is set to "local"
## Nog even over nadenken hoe dit wordt aangepast door de webservice...
if(arrayQC.mode == "local") {
  CheckVirtualImages <- TRUE
  CheckBoxplot <- TRUE
  CheckHeatmap <- TRUE
  CheckClustering <- TRUE
  CheckPCA <- TRUE
  CheckCorplot <- TRUE
  CheckMvA <- TRUE
}


## If you want to run arrayQC locally, download the files from the scriptpath mentioned above
#  and refer to these files in the scriptpath variable. Recommended for experienced R programmers.

source(paste(scriptpath, "run_arrayQC.R", sep=""))


###########################
##    TROUBLESHOOTING    ##
###########################
## 1) I get a 'could not allocate memory vector' error. What should I do?
# Answer: increase your available memory in R by setting memory.limit(2000) prior to running the script
# (here, 2000 corresponds with 2Gb of physical RAM)

