###
### Pre-treatment Arabidopsis data ###
###
file <- read.csv("~/Documents/DATA/Arabidopsis/WT_normalized_protein_abundance.csv",skip = 1, h=T)
file.tmp <- file[, c(1,44:73)]

# Sometimes, several genes exist for 1 identification: 
# because related proteins are very similar and only one or a few peptides 
# were identified that it is impossible to tell which of the three proteins it is from. 
# Then the Progenesis software put all three of them in. 

# We are going to ignore them:
file.tmp <- file.tmp[-grep(";", file.tmp$Accession), ]

# Lets consider biological replicates as new cycles: 
#colnames(file.tmp) <- gsub("A|B|C|D|E", "", colnames(file.tmp))
#colnames(file.tmp) <- gsub("wt.", "CT", colnames(file.tmp))
vector.i <- c(1,6,11,16,21,26)
vector <- NULL
for(i in 1:5){
  vector <- c(vector, vector.i+i)
}
vector <- c(1, vector)

file.tmp <- file.tmp[, vector]
time.points <- seq(from = 12,to = 32*4, by = 4)
time.points <- paste("CT", time.points, sep="")
colnames(file.tmp) <- c("ID", time.points)

write.table(file.tmp, "~/Documents/DATA/Arabidopsis/prot_level.txt", row.names = F, quote = F, sep = "\t")


###
### Pre-treatment Ostreococcus tauri data ###
###
file <- read.csv("~/Documents/DATA/Ostreococcus_tauri/protein_level.csv",skip = 2, h=T)
# Retreive protein abundance (constant protein abundance)
file.tmp <- file[, c(5,2)]
colnames(file.tmp) <- c("Gene.ID", "abundance")
write.table(file.tmp, "~/Documents/DATA/Ostreococcus_tauri/prot_abundance.txt", row.names = F, quote = F, sep = "\t")


file.tmp <- file[, c(5,7:36)]
# Lets consider biological replicates as new cycles: 
vector.i <- c(1,6,11,16,21,26)
vector <- NULL
for(i in 1:5){
  vector <- c(vector, vector.i+i)
}
vector <- c(1, vector)

file.tmp <- file.tmp[, vector]
time.points <- seq(from = 0,to = 20*6-1, by = 4)
time.points <- paste("ZT", time.points, sep="")
colnames(file.tmp) <- c("ID", time.points)

write.table(file.tmp, "~/Documents/DATA/Ostreococcus_tauri/prot_level.txt", row.names = F, quote = F, sep = "\t")

# Transcripts level
transcript.data <- read.table("~/Documents/DATA/Ostreococcus_tauri/GSE16422_transcript_timeseries.txt",sep = "\t", head=T, fill=T, check.names = F)

# Modify correctly gene.ID in the UniProt file: 
file <- read.table("~/Documents/DATA/Ostreococcus_tauri/uniprot.txt",sep = "\t", head=T, fill=T, check.names = F)
file$Gene.Name <- gsub("\\ .*", "", file$Gene.Name)
write.table(file, "~/Documents/DATA/Ostreococcus_tauri/uniprot.txt", row.names = F, quote = F, sep = "\t")


file <- read.csv("~/Downloads/pgen.1004047.s006.csv",head=T, fill=T, check.names = F)
length(unique(file$`Gene name`))


file$Gene.Name <- gsub("\\ .*", "", file$Gene.Name)
write.table(file, "~/Documents/DATA/Ostreococcus_tauri/uniprot.txt", row.names = F, quote = F, sep = "\t")



