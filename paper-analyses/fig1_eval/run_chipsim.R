library("ChIPsim")

######################################################
# params
args = commandArgs(trailingOnly=TRUE)
nReads = as.integer(args[1])
genome = args[2]
bedfile = args[3]
outprefix = args[4]
readLen = as.integer(args[5])

######################################################
# Define a "binding" feature for TFBD
bindingFeature <- function(start, length=200, shape=1, scale=20, enrichment=5, r=1.5){
 stopifnot(r > 1)
 avgWeight <- shape * scale * enrichment
 lowerBound <- ((r - 1) * avgWeight)
 weight <- 1
 params <- list(start = start, length = length, weight = weight, overlap = 0)
 class(params) <- c("Binding", "SimulatedFeature")
 params
}

# TODO - background as non-binding sites. ugh I have to define myself??

# Update "features" function to read from the BED file
getpeakfeatures <- function(..., maxTail=0.01, compoundFeatures=list("Binding")) {
   features = list()

   fnum = 1
   con = file(bedfile, "r")
   while(TRUE) {
     line = readLines(con, n=1)
     if (length(line) == 0) {
       break
     }
     # TODO need to use length
     start = as.integer(strsplit(line, '\t')[[1]][2])
     end = as.integer(strsplit(line, '\t')[[1]][3])
     f = bindingFeature(start)
     features[[fnum]] = f
     fnum = fnum + 1
   }
   close(con)
   head(features, 100)
}

constRegion <- function(weight, length) rep(weight, length)
featureDensity.Binding <- function(feature, ...) constRegion(feature$weight, feature$length)	

######################################################
# Custom functions for writing reads
randomQuality <- function(read, ...) {
  paste(sample(unlist(strsplit(rawToChar(as.raw(64:104)),"")),
  nchar(read), replace = TRUE), collapse="")
}

myWriteReads <- function(readPos, readNames, sequence, quality, file, ...) {
  fileName <- file
  file <- file(file, open="w", blocking=FALSE)
  for(i in 1:length(sequence))
  pos2fastq(readPos=readPos[[i]], names=readNames[[i]], sequence[[i]], quality=quality, file=file, ...)
  close(file)
  fileName
}

######################################################
# Set up main functions item
myFunctions <- defaultFunctions()
myFunctions$features = getpeakfeatures
myFunctions$readSequence = myWriteReads

myControl = defaultControl(readSequence=list(qualityFun=randomQuality, readLen=readLen))
simulated <- simChIP(nReads, genome, file=outprefix, functions=myFunctions, verbose=TRUE, load=FALSE, control=myControl)
