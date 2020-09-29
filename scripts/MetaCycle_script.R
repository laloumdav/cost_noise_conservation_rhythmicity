#################
### METACYCLE ###
#################
library(readr)
library(plyr)
library(MetaCycle)


if (dir.exists("/scratch/cluster/monthly/dlaloum/")==TRUE){
  main.dir <- "/scratch/cluster/monthly/dlaloum/Documents/cost_theory_workingspace"
} else {
  main.dir <- "~/Documents/cost_theory_workingspace"
}

file.dir <- paste(paste(main.dir, "DATA/cyanobacteria/transcriptome", sep = "/"), "/", sep="")
tissue.file <- "transcript_level.txt"

file.name <- paste(file.dir, tissue.file, sep = "")
raw.dataset <- read.table(file.name, head=TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, fill=TRUE)

time.points <- colnames(raw.dataset)
time.points <- time.points[grep("CT|ZT|LD|DD|LL", time.points)]
time.points <- parse_number(time.points)

### MetaCycle : 
metaCycle.output.default <- meta2d(infile=file.name, filestyle="txt",
                             outputFile=FALSE,
                             minper = 23, maxper = 25,
                             timepoints=time.points,
                             cycMethod=c("ARS", "LS", "JTK"), outRawData=FALSE) 

meta.default <- metaCycle.output.default$meta

ARS <- meta.default[ , c("CycID", "ARS_pvalue", "ARS_period", "ARS_adjphase", "ARS_amplitude")]
colnames(ARS) <- c("ID", "default.pvalue", "period", "phase", "amplitude")
#JTK <- meta.default[ , c("CycID", "JTK_pvalue", "JTK_period", "JTK_adjphase", "JTK_amplitude")]
#colnames(JTK) <- c("ID", "default.pvalue", "period", "phase", "amplitude")
LS <- meta.default[ , c("CycID", "LS_pvalue", "LS_period", "LS_adjphase", "LS_amplitude")]
colnames(LS) <- c("ID", "default.pvalue", "period", "phase", "amplitude")
meta2d <- meta.default[ , c("CycID", "meta2d_pvalue", "meta2d_period", "meta2d_phase", "meta2d_AMP")]
colnames(meta2d) <- c("ID", "default.pvalue", "period", "phase", "amplitude")




if (exists("ARS")){ write.table(ARS, paste(file.dir, "transcript_ARS.txt", sep = ""), row.names = F, quote = F, sep = "\t") }
if (exists("JTK")){ write.table(JTK, paste(file.dir, "transcript_JTK.txt", sep = ""), row.names = F, quote = F, sep = "\t") }
if (exists("LS")){ write.table(LS, paste(file.dir, "transcript_LS.txt", sep = ""), row.names = F, quote = F, sep = "\t") }
if (exists("meta2d")){ write.table(meta2d, paste(file.dir, "transcript_meta2d.txt", sep = ""), row.names = F, quote = F, sep = "\t") }


