#####################
#### GeneCycle ######
#####################

library(readr)
library(plyr)
# Change directory name depending on if it's local or on cluster
if (dir.exists("/scratch/cluster/monthly/dlaloum/")==TRUE){
  main.dir <- "/scratch/cluster/monthly/dlaloum/Documents/cost_theory_workingspace"
} else {
  main.dir <- "~/Documents/cost_theory_workingspace"
}

install.packages(paste(main.dir, "scripts/GeneCycle", sep = "/"), repos = NULL, type = "source")
library(GeneCycle)

file.dir <- paste(paste(main.dir, "DATA/arabidopsis/leaves/transcriptome", sep = "/"), "/", sep="")
tissue.file <- "transcript_level.txt"


### SCRIPT ###
system(paste("echo Running the file :", tissue.file, sep=" "))

file.name <- paste(file.dir, tissue.file, sep = "")
raw.dataset <- read.table(file.name, head=TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, fill=TRUE)


# time.points
time.points <- colnames(raw.dataset)
time.points <- time.points[grep("CT|ZT|LD|DD|LL", time.points)]
time.points <- parse_number(time.points)

# Number of replicates for each timepoint
if ( length(duplicated(time.points) == FALSE) != length(duplicated(time.points) == FALSE) ){
  print("ERROR: The number of replicates isn't identical between the timepoints")
} else { replicates.number <- as.numeric(unique(plyr::count(time.points)$freq)) }


# If there is several replicates per time-point, we are going to transform replicates as new days measurements 
# Indeed, robust.spectrum function of GeneCycle can not deal with several replicates per time-point

if ( any(duplicated(time.points)) ) {
  GeneCycle.time.points <- time.points[duplicated(time.points)==FALSE]
  nb.of.replicates <- unique(plyr::count(time.points)$freq)
  GeneCycle.seq <- seq(from=2, to=ncol(raw.dataset), by=nb.of.replicates)
  # interval.time
  interval.time <- GeneCycle.time.points[3]-GeneCycle.time.points[2]
  
  GeneCycle.raw.dataset <- raw.dataset[ , c(1, GeneCycle.seq)]
  for (i in 1:(nb.of.replicates-1)) {
    GeneCycle.raw.dataset <- cbind(GeneCycle.raw.dataset, raw.dataset[ , c(GeneCycle.seq+i)])
  }
  
  unique.time.points <- seq(from = min(GeneCycle.time.points), by = interval.time, length.out = length(time.points))
  length(unique.time.points)
  
  ### GeneCycle :
  raw.dataset.IDs <- as.data.frame(GeneCycle.raw.dataset[,1])
  raw.dataset.transposed <- t(GeneCycle.raw.dataset[,-1])
  
} else {
  # interval.time
  unique.time.points <- unique(time.points)
  
  ### GeneCycle :
  raw.dataframe <- raw.dataset
  #raw.dataframe[raw.dataframe=="-"] <- NA
  #raw.dataframe <- na.omit(raw.dataframe)
  raw.dataset.IDs <- raw.dataframe$ID
  
  raw.dataset.transposed <- t(raw.dataframe[,-1])
  #raw.dataset.transposed <- apply(raw.dataset.transposed, c(1,2), FUN = as.numeric)
}

#raw.dataset.transposed <- raw.dataset.transposed + abs(min(raw.dataset.transposed, na.rm = TRUE))
raw.dataset.transposed <- raw.dataset.transposed+0.0001 

## Lets use the robust regression based approach (Ahdesmaki et al. 2007) 
## with robust.spectrum function. Because we known the periodicity time.
## (computation will take a lot of time depending on how many permutations are used per time series and time series length).
GeneCycle.output <- robust.spectrum(x = raw.dataset.transposed,
                                    t = unique.time.points,
                                    periodicity.time = 24,
                                    #noOfPermutations = 150,
                                    algorithm = "regression")

GeneCycle <- data.frame(ID = raw.dataset.IDs, 
                        default.pvalue = GeneCycle.output)
hist(GeneCycle$default.pvalue, breaks = 150)
colnames(GeneCycle)[1] <- "ID"

#GeneCycle <- subset(GeneCycle, default.pvalue =! NA)

write.table(GeneCycle, paste(file.dir, "transcript_GeneCycle.txt", sep = ""), row.names = F, quote = F, sep = "\t")




