### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
# Mouse ######
### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
main.dir <- "~/Documents/DATA/Mus_Musculus/Liver_Protein_Mauvoisin/"
file.dir <- main.dir
species <- "mouse"
tissue <- "liver"

### SCRIPT ###
file.name <- paste(file.dir, "protein_level_raw.csv", sep="")
raw.dataset <- read.csv(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Gene.ID and Protein.ID : 
file.dir <- "~/Documents/DATA/Genes_ID/"
species.genesID.protID.file <- paste(file.dir, "GeneID_GeneName_TranscriptID_ProtID_MOUSE.txt", sep="")
geneID.proteinID <- read.table(species.genesID.protID.file, head=TRUE, fill = TRUE)
geneID.proteinID <- unique(subset(geneID.proteinID, Protein.ID != ""))
geneID.proteinID <- geneID.proteinID[, c("Gene.ID", "Gene.Name")]

raw.dataset.prot <- raw.dataset[, c(1,4, 6:21)]
raw.dataset.mRNA <- raw.dataset[, c(1,4, 28:51)]
colnames(raw.dataset.prot) <- gsub("nZT", "ZT", colnames(raw.dataset.prot))
colnames(raw.dataset.prot)[2] <- "Protein.Name"
colnames(raw.dataset.prot)[1] <- "Gene.Name"
colnames(raw.dataset.mRNA)[2] <- "Protein.Name"

raw.dataset.prot <- raw.dataset[, c(1,4, 6:21)]
raw.dataset.mRNA <- raw.dataset[, c(1,4, 28:51)]
colnames(raw.dataset.prot) <- gsub("nZT", "ZT", colnames(raw.dataset.prot))
colnames(raw.dataset.prot)[2] <- "Protein.Name"
colnames(raw.dataset.mRNA)[2] <- "Protein.Name"

raw.dataset.prot <- raw.dataset.prot[,-1]
raw.dataset.mRNA <- raw.dataset.mRNA[,-1]

main.dir <- "~/Documents/DATA/Mus_Musculus/Liver_Protein_Mauvoisin/"
file.dir <- main.dir
write.table(raw.dataset.prot, paste(file.dir, "protein_level.txt", sep="/"), row.names = F, quote = F, sep = "\t")
write.table(raw.dataset.mRNA, paste(file.dir, "mRNA_level.txt", sep="/"), row.names = F, quote = F, sep = "\t")


# raw protein data
# original timeseries data (from http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001211) 
protein.abundance.dataset <- read.table(paste(proteome.dir, "proteinGroups.txt", sep="/"), sep="\t", fill= TRUE, h=T,  check.names = FALSE, stringsAsFactors = FALSE)
normalized.count.values <- grep("Ratio H/L coun", colnames(protein.abundance.dataset))
protein.abundance.dataset <- protein.abundance.dataset[, c(2, 6, 7, normalized.count.values)]
protein.abundance.dataset <- protein.abundance.dataset[, -4]
colnames(protein.abundance.dataset)[1:3] <- c("Uniprot.ID", "Protein.Name", "Gene.Name")
colnames(protein.abundance.dataset) <- gsub("Ratio H/L count ", "", colnames(protein.abundance.dataset))

protein.level <- read.table(paste(proteome.dir, "protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#protein.level <- protein.level[, 1:2]

merged.data <- merge(protein.level, protein.abundance.dataset, by=c("Protein.Name", "Gene.Name"), sort = FALSE)
counts <- plyr::ddply(merged.data, .(merged.data$Protein.Name, merged.data$Gene.Name), nrow)
names(counts) <- c("Protein.Name", "Gene.Name", "Freq")
merged.data <- merge(merged.data, counts, by=c("Protein.Name", "Gene.Name"), sort = FALSE)

merged.data <- merged.data[, c("Protein.Name", "Gene.Name", "Freq")]

protein.level.merged <- merge(protein.level, merged.data, by=c("Protein.Name", "Gene.Name"), sort = FALSE)
protein.level.merged <- subset(protein.level.merged, Freq ==1)[-ncol(protein.level.merged)]

write.table(protein.level.merged, paste(proteome.dir, "protein_level.txt", sep="/"), quote = F, row.names = F, sep = "\t")


## Transcripts ##
# original timeseries data (from http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD001211) 
protein.abundance.dataset <- read.table(paste(proteome.dir, "proteinGroups.txt", sep="/"), sep="\t", fill= TRUE, h=T,  check.names = FALSE, stringsAsFactors = FALSE)
normalized.count.values <- grep("Ratio H/L coun", colnames(protein.abundance.dataset))
protein.abundance.dataset <- protein.abundance.dataset[, c(2, 6, 7, normalized.count.values)]
protein.abundance.dataset <- protein.abundance.dataset[, -4]
colnames(protein.abundance.dataset)[1:3] <- c("Uniprot.ID", "Protein.Name", "Gene.Name")
colnames(protein.abundance.dataset) <- gsub("Ratio H/L count ", "", colnames(protein.abundance.dataset))

raw.data <- read.csv(paste(proteome.dir, "complete_data_from_original_paper.csv", sep="/"), head=TRUE)[, -2]
colnames(raw.data)[1:3] <- c("Gene.Name", "Uniprot.ID", "Protein.Name")

merged.data <- merge(raw.data, protein.abundance.dataset, by=c("Protein.Name", "Gene.Name"), sort = FALSE)
#library("plyr")
counts <- plyr::ddply(merged.data, .(merged.data$Protein.Name, merged.data$Gene.Name), nrow)
names(counts) <- c("Protein.Name", "Gene.Name", "Freq")

merged.data <- merge(merged.data, counts, by=c("Protein.Name", "Gene.Name"), sort = FALSE)

merged.data <- merged.data[, c("Protein.Name", "Gene.Name", "Freq")]

protein.level.merged <- merge(raw.data, merged.data, by=c("Protein.Name", "Gene.Name"), sort = FALSE)
protein.level.merged <- subset(protein.level.merged, Freq ==1)[-ncol(protein.level.merged)]
protein.level.merged <- protein.level.merged[, c("Protein.Name", "Gene.Name", "Uniprot.ID", "mRNA.ID", "ZT00", "ZT02", "ZT04", "ZT06", "ZT08", paste("ZT", seq(10, 46, by=2), sep=""))]

write.table(protein.level.merged, paste(proteome.dir, "protein_harmonicRegression.txt", sep="/"), quote = F, row.names = F, sep = "\t")

### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 


### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### Krill LD data (transcripts):
file <- read.table("~/Documents/cost_theory_workingspace/DATA/krill/GSE94756_series_matrix.txt", h=T, fill=T,check.names = FALSE, stringsAsFactors = FALSE)
# remove duplicated time-point LD06
file <- file[, c(-ncol(file), -ncol(file)+1, -ncol(file)+2)]
#change.order <- c(1, 
#                  seq(from = 2, to = ncol(file), by = 3), 
#                  seq(from = 3, to = ncol(file), by = 3), 
#                  seq(from = 4, to = ncol(file), by = 3))
#file <- file[, change.order]
file <- file[, c(1, ncol(file), ncol(file)-1, ncol(file)-2, 2:(ncol(file)-3))]
column.names <- rep(c("LD03", "LD06", "LD09", "LD12", "LD15", "LD18", "LD21", "LD24"), each=3)
colnames(file) <- c("ID", column.names)

file[file=="null"] <- NA
# Keep only transcripts which have at least 10 data:
file <- file[rowSums(is.na(file)) < ncol(file)-11, ]

file <- file[-nrow(file), ]

write.table(file, "~/Documents/cost_theory_workingspace/DATA/krill/transcript_level.txt", row.names = F, quote = F, sep = "\t")
### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 



### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### Cyanobacteria data (transcriptome):
cyano.transcript <- read.table("~/Documents/cost_theory_workingspace/DATA/cyanobacteria/transcriptome/GSE14225_series_matrix.txt", h=T, fill=T)
wt.timeseries <- c(paste("GSM35640", 1:9, sep=""), paste("GSM3564", 10:22, sep=""))
colnames(cyano.transcript)[1] <- "ID"
cyano.transcript <- cyano.transcript[, c("ID", wt.timeseries)]
# We remove the 2 last time-points of rep2 since LL44 is missing 
colunm.names <- paste("LL", c(seq(4, 48, by=4), seq(4, 40, by=4)), sep="")
colnames(cyano.transcript) <- c("ID", colunm.names)

write.table(cyano.transcript, "~/Documents/cost_theory_workingspace/DATA/cyanobacteria/transcriptome/transcript_level.txt", sep = "\t", quote = F, row.names = F)
### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 


### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### Cyanobacteria data (proteome):
cyano.prot <- read.csv("~/Documents/cost_theory_workingspace/DATA/cyanobacteria/proteome/initial_data.csv", head=TRUE, skip=4, fill=TRUE)
cyano.prot <- cyano.prot[, c(1, 3, 4, grep("X", colnames(cyano.prot)))]
cyano.prot[cyano.prot == "-"] <- NA
cyano.prot$counting.na <- apply(cyano.prot[,c(-1,-2,-3)], 1, FUN = function(x){sum(is.na(x))})
# Remove data with NA at some time-points
cyano.prot <- subset(cyano.prot, counting.na<1)[, -ncol(cyano.prot)]
colnames(cyano.prot)[4:ncol(cyano.prot)] <- gsub("X", "LD", colnames(cyano.prot)[4:ncol(cyano.prot)])
colnames(cyano.prot)[4:ncol(cyano.prot)] <- gsub(".hours.LD", "", colnames(cyano.prot)[4:ncol(cyano.prot)])
cyano.prot <- cyano.prot[-nrow(cyano.prot), ]
write.table(cyano.prot[, 1:3], "~/Documents/cost_theory_workingspace/DATA/cyanobacteria/proteome/proteinNames_initial_data.txt", sep = "\t", quote = F, row.names = F)
colnames(cyano.prot)[1] <- c("ID")
write.table(cyano.prot[, c(-2,-3)], "~/Documents/cost_theory_workingspace/DATA/cyanobacteria/proteome/protein_level.txt", sep = "\t", quote = F, row.names = F)
#...
cyano.prot <- read.table("~/Documents/cost_theory_workingspace/DATA/cyanobacteria/proteome/protein_level.txt", head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
cyano.prot[cyano.prot == "-"] <- NA
cyano.prot$counting.na <- apply(cyano.prot[,-1], 1, FUN = function(x){sum(is.na(x))})
cyano.prot <- subset(cyano.prot, counting.na<1)[, -ncol(cyano.prot)]
write.table(cyano.prot, "~/Documents/cost_theory_workingspace/DATA/cyanobacteria/proteome/protein_level.txt", sep = "\t", quote = F, row.names = F)
### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 


### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ###
### arabidopsis "whole" transcriptome data :
file <- read.csv("~/Documents/cost_theory_workingspace/DATA/arabidopsis/whole/transcriptome/processed_initial_data.txt", 
                 sep="\t", head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
file <- file[-1,]
file <- file[, c(1, grep("LD", colnames(file)))]
file <- file[, 1:13]
colnames(file) <- gsub("col_", "LD", colnames(file))
colnames(file) <- gsub("_LDHC", "", colnames(file))
colnames(file)[1:4] <- c("ID", "LD00", "LD04", "LD08")

arabidopsis.IDs <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/whole/transcriptome/GeneIDs.txt", head=TRUE, fill = TRUE)
file <- merge(arabidopsis.IDs[, c("Gene.ID", "AFFY.ATH1")], file, by.x= "AFFY.ATH1", by.y = "ID")

# THEN: We remove prob.ID if it is affected to several Gene.ID (long process) 
raw.dataset.merged <- file
column.name <- "AFFY.ATH1"
rows.to.remove <- c()
for (i in 1:nrow(raw.dataset.merged)) {
  if ((i %in% rows.to.remove)==FALSE) {
    prob.ID <- raw.dataset.merged[i, column.name]
    if (length(unique(raw.dataset.merged$Gene.ID[grep(prob.ID, raw.dataset.merged[, column.name])])) != 1){
      rows.to.remove <- c(rows.to.remove, grep(prob.ID, raw.dataset.merged[, column.name]))
    }
  }
}

rows.to.remove <- unique(rows.to.remove)
# So, we remove prob.IDs refered to several Gene.IDs:
if (is.null(rows.to.remove)==FALSE){
  raw.dataset.merged <- raw.dataset.merged[-rows.to.remove, ]
}
#length(unique(raw.dataset.merged$AFFY.ATH1))

write.table(raw.dataset.merged[,2:1], "~/Documents/cost_theory_workingspace/DATA/arabidopsis/whole/transcriptome/probID_GeneID.txt", row.names = F, quote = F, sep = "\t")

file <- raw.dataset.merged[,2:ncol(raw.dataset.merged)]
colnames(file)[1] <- "ID"
write.table(file, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/whole/transcriptome/transcript_level.txt", row.names = F, quote = F, sep = "\t")
### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 


### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### Arabidopsis (Leaves: transcriptome):
file <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/transcriptome/GSE3416_series_matrix.txt", 
                   head=TRUE, fill=T)
time.points <- rep(c("LD00", "LD04", "LD08", "LD12", "LD16", "LD20"), each=3)
colnames(file) <- c("ID", time.points)
# consider biological replicates as new cycles
ordered.columns <- c(1, seq(2, 19, by=3), seq(3, 19, by=3), seq(4, 19, by=3))
file <- file[, ordered.columns]
time.points <- paste("LD", c("00", "04", "08", seq(12, 68, by=4)), sep="")
colnames(file) <- c("ID", time.points)

write.table(file, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/transcriptome/transcript_level.txt", row.names = F, quote = F, sep = "\t")
### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 


### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ###
library(sjmisc)
### arabidopsis HISTONES MODIFICATIONS Data, Leaves
file.1 <- read.csv("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/epigenetics/H3K4me3_initial_data.csv", 
                   head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, skip=1)
file.2 <- read.csv("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/epigenetics/H3K9ac_initial_data.csv", 
                   head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
genes.positions <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/epigenetics/arabidopsis_genes_positions.txt", 
                              head=TRUE, fill = TRUE, sep="\t")

file.1.pos <- file.1[, 1:4]
file.2.pos <- file.2[, 1:4]

kept.data.1 <- NULL
for (p in 1:nrow(file.1.pos)) {
  peak.start <- file.1.pos[p, "Start"]
  peak.stop <- file.1.pos[p, "Stop"]
  chr <- as.numeric(gsub("Chr", "", file.1.pos[p, "Chromosome"]))
  
  overlap.genes <- genes.positions[as.numeric(genes.positions$Chromosome) == chr &
                                     as.numeric(genes.positions$Transcription.start.site) <= as.numeric(peak.start) &
                                     as.numeric(genes.positions$Transcript.end) >= as.numeric(peak.stop), ]
  # Consider the sequence as part of the promoter if included into -1500 to +700 base from the TSS
  overlap.promoter <- genes.positions[as.numeric(genes.positions$Chromosome) == chr &
                                        as.numeric(genes.positions$Transcription.start.site)-1500 <= as.numeric(peak.start) &
                                        as.numeric(genes.positions$Transcription.start.site)+700 >= as.numeric(peak.stop), ]
  if (length(unique(overlap.genes$Gene.ID)) == 1) {
    if (is_empty(overlap.promoter)) {
      overlap.genes <- cbind(overlap.genes[, 1:3], file.1[p, ])
      overlap.genes$promoter.region <- rep(FALSE, nrow(overlap.genes))
      kept.data.1 <- rbind(kept.data.1, overlap.genes)
    } else {
      overlap.genes <- cbind(overlap.genes[, 1:3], file.1[p, ])
      overlap.genes$promoter.region <- rep(TRUE, nrow(overlap.genes))
      kept.data.1 <- rbind(kept.data.1, overlap.genes)
    }
  } else if (length(unique(overlap.promoter$Gene.ID)) == 1) {
    overlap.promoter <- cbind(overlap.promoter[, 1:3], file.1[p, ])
    overlap.promoter$promoter.region <- rep(TRUE, nrow(overlap.promoter))
    kept.data.1 <- rbind(kept.data.1, overlap.promoter)
  }
}

kept.data.2 <- NULL
for (p in 1:nrow(file.2.pos)) {
  peak.start <- file.2.pos[p, "Start"]
  peak.stop <- file.2.pos[p, "Stop"]
  chr <- as.numeric(gsub("Chr", "", file.2.pos[p, "Chromosome"]))
  
  overlap.genes <- genes.positions[as.numeric(genes.positions$Chromosome) == chr &
                                     as.numeric(genes.positions$Transcription.start.site) <= as.numeric(peak.start) &
                                     as.numeric(genes.positions$Transcript.end) >= as.numeric(peak.stop), ]
  # Consider the sequence as part of the promoter if included into -1500 to +700 base from the TSS
  overlap.promoter <- genes.positions[as.numeric(genes.positions$Chromosome) == chr &
                                        as.numeric(genes.positions$Transcription.start.site)-1500 <= as.numeric(peak.start) &
                                        as.numeric(genes.positions$Transcription.start.site)+700 >= as.numeric(peak.stop), ]
  if (length(unique(overlap.genes$Gene.ID)) == 1) {
    if (is_empty(overlap.promoter)) {
      overlap.genes <- cbind(overlap.genes[, 1:3], file.2[p, ])
      overlap.genes$promoter.region <- rep(FALSE, nrow(overlap.genes))
      kept.data.2 <- rbind(kept.data.2, overlap.genes)
    } else {
      overlap.genes <- cbind(overlap.genes[, 1:3], file.2[p, ])
      overlap.genes$promoter.region <- rep(TRUE, nrow(overlap.genes))
      kept.data.2 <- rbind(kept.data.2, overlap.genes)
    }
  } else if (length(unique(overlap.promoter$Gene.ID)) == 1) {
    overlap.promoter <- cbind(overlap.promoter[, 1:3], file.2[p, ])
    overlap.promoter$promoter.region <- rep(TRUE, nrow(overlap.promoter))
    kept.data.2 <- rbind(kept.data.2, overlap.promoter)
  }
}

#### #### #### 
#### H3K4me3 
#### #### #### 
# Lets keep data only per unique gene, with TRUE or FALSE for the promoter information 
# In case of multiples peaks for one gene, we calculate the mean of the intensity at each time-point
kept.data.1_unique.genes <- unique(kept.data.1[, c(1, 8:(ncol(kept.data.1)))])
kept.data.1_unique.genes <- aggregate(.~Gene.ID+promoter.region,
                  data = kept.data.1_unique.genes,
                  FUN = mean)
# Lets consider biological replicates as new cycles
# We remove ZT24 to get regular intervals
kept.data.1_unique.genes <- kept.data.1_unique.genes[, c(1, 2, 
                                                         seq(3, ncol(kept.data.1_unique.genes)-4, by=3),
                                                         seq(4, ncol(kept.data.1_unique.genes)-4, by=3),
                                                         seq(5, ncol(kept.data.1_unique.genes), by=3))]
colnames(kept.data.1_unique.genes) <- c("ID", "promoter.region", paste("ZT", seq(0, 72, by=3), sep=""))

#### #### #### 
#### H3K9ac 
#### #### #### 
kept.data.2_unique.genes <- unique(kept.data.2[, c(1, 8:(ncol(kept.data.2)))])
kept.data.2_unique.genes <- aggregate(.~Gene.ID+promoter.region,
                                      data = kept.data.2_unique.genes,
                                      FUN = mean)
kept.data.2_unique.genes <- kept.data.2_unique.genes[, c(1, 2,
                                                         seq(3, ncol(kept.data.2_unique.genes)-4, by=3),
                                                         seq(4, ncol(kept.data.2_unique.genes)-4, by=3),
                                                         seq(5, ncol(kept.data.2_unique.genes), by=3))]
colnames(kept.data.2_unique.genes) <- c("ID", "promoter.region", paste("ZT", seq(0, 72, by=3), sep=""))

write.table(kept.data.1_unique.genes, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/epigenetics/H3K4me3_data.txt", row.names = F, quote = F, sep = "\t")
write.table(kept.data.2_unique.genes, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/epigenetics/H3K9ac_data.txt", row.names = F, quote = F, sep = "\t")





### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### Ostreococcus data (transcriptome):
file <- read.table("~/Documents/cost_theory_workingspace/DATA/ostreococcus/transcriptome/GSE16422_transcript_timeseries.txt", 
                   head=TRUE, sep="\t", check.names = FALSE, stringsAsFactors = FALSE, fill=T)
file <- subset(file, CDS_ID_II != "****")
file <- file[, c(3, grep("GSM", colnames(file)))]
# remove the duplicated LD09 (GSM412697, GSM412706 and GSM412715):
file <- file[, c(1, 2:09, 11:18, 20:27)]
time.points <- rep(c("LD09", "LD12", "LD15", "LD18", "LD21", "LD24", "LD03", "LD06"), 3)
colnames(file) <- c("ID", time.points)
file <- file[, c(1, 8, 9, 2:7, 16, 17, 10:15, 24, 25, 18:23)]
file <- file[, order(colnames(file))]

# Remove rows with to much NA (let say more than 7 NA) : 
file$nb.NA <- apply(file[,-1], 1, FUN = function(x){sum(is.na(x))}) # compte le nombre de NA pour chaque ligne
#count(dataset$nb.NA)
file <- subset(file, nb.NA<=7)[, -ncol(file)]

time.points <- rep(c("LD03", "LD06", "LD09", "LD12", "LD15", "LD18", "LD21", "LD24"), each=3)
colnames(file) <- c("ID", time.points)

#length(unique(file$ID))

write.table(file, "~/Documents/cost_theory_workingspace/DATA/ostreococcus/transcriptome/transcript_level.txt", row.names = F, quote = F, sep = "\t")

# Lets do the same for the ProbIDs-GeneIDs correspondance file
# Gene.ID and Protein.ID :
geneIDs <- read.table("~/Documents/cost_theory_workingspace/DATA/ostreococcus/transcriptome/GSE16422_transcript_timeseries.txt", 
                      head=TRUE, sep="\t", check.names = FALSE, stringsAsFactors = FALSE, fill=T)
geneIDs <- geneIDs[!(geneIDs$CDS_ID_II == "****"),]
geneIDs <- geneIDs[, c(2, 3, grep("GSM", colnames(geneIDs)))]
# remove the duplicated LD09 (GSM412697, GSM412706 and GSM412715):
geneIDs <- geneIDs[, c(1, c(1, 2:09, 11:18, 20:27)+1)]
geneIDs$counting.na <- apply(geneIDs[, c(-1, -2)], 1, FUN = function(x){sum(is.na(x))})
geneIDs <- subset(geneIDs, counting.na<=7)
geneIDs <- geneIDs[, 1:2]

write.table(geneIDs, "~/Documents/cost_theory_workingspace/DATA/ostreococcus/transcriptome/transcriptIDs.txt", row.names = F, quote = F, sep = "\t")
### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 

### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### ostreococcus data (proteome):
file <- read.csv("~/Documents/cost_theory_workingspace/DATA/ostreococcus/proteome/protein_level.csv", head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, skip=2)
file <- file[,c(5, grep("LD.R", colnames(file)))]
column.names <- c("ID", rep(c("ZT00", "ZT04", "ZT08", "ZT12", "ZT16", "ZT20"), each=5))
file <- file[, 1:length(column.names)]
colnames(file) <- column.names

write.table(file, "~/Documents/cost_theory_workingspace/DATA/ostreococcus/proteome/protein_level_bis.txt", row.names = F, quote = F, sep = "\t")
### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 






### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### Mouse (Cartilage: proteome):
inital.data <- read.csv("~/Documents/cost_theory_workingspace/DATA/mouse/cartilage/proteome/initial_data.csv", head=TRUE, fill=TRUE)
#ARS.data <- inital.data[, c("UNIPROTKB", "ARS_pvalue")]
#colnames(ARS.data) <- c("ID", "default.pvalue")
#write.table(ARS.data, "~/Documents/cost_theory_workingspace/DATA/mouse/cartilage/proteome/protein_ARS.txt", row.names = F, quote = F, sep = "\t")
prot.data <- inital.data[, c(ncol(inital.data)-2, ncol(inital.data)-1, grep("Time_point", colnames(inital.data)))]
colnames(prot.data) <- gsub("Time_point_", "ZT", colnames(prot.data))
colnames(prot.data)[1:4] <- c("Uniprot.ID", "Gene.Name", "ZT02", "ZT06")

#proteome.data <- merge(prot.data[, c("Uniprot.ID", "Gene.Name", "Protein.Description")], cartilage.proteome.data, by="Uniprot.ID")
#write.table(proteome.data, "~/Documents/cost_theory_workingspace/DATA/mouse/cartilage/proteome/proteome_data.txt", row.names = F, quote = F, sep = "\t")
write.table(prot.data, "~/Documents/cost_theory_workingspace/DATA/mouse/cartilage/proteome/protein_level.txt", row.names = F, quote = F, sep = "\t")


### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### Mouse (Tendon: proteome):
inital.data <- read.csv("~/Documents/cost_theory_workingspace/DATA/mouse/tendon/proteome/initial_data.csv", head=TRUE, fill=TRUE)
prot.data <- inital.data[, c(2, grep("CT", colnames(inital.data)))]
colnames(prot.data)[1:3] <- c("ID", "CT03", "CT07")
write.table(prot.data, "~/Documents/cost_theory_workingspace/DATA/mouse/tendon/proteome/protein_level.txt", row.names = F, quote = F, sep = "\t")



## ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### 
### Mouse (Forebrain: proteome):
inital.data <- read.csv("~/Documents/cost_theory_workingspace/DATA/mouse/forebrain/proteome/initial_data.csv", head=TRUE, fill=TRUE, skip=8)
prot.data <- inital.data[, c(1, grep("ZT", colnames(inital.data)))]
prot.data <- prot.data[, c(1, 20, 2, 6, 8, 11, 14, 18,
                           3, 7, 9, 12, 15, 19, 4)]
colnames(prot.data) <- gsub("_1", "", colnames(prot.data))
colnames(prot.data) <- gsub("_2", "", colnames(prot.data))
colnames(prot.data) <- gsub("_3", "", colnames(prot.data))
colnames(prot.data) <- c("ID", paste("ZT", seq(20, 72, by=4), sep=""))
write.table(prot.data, "~/Documents/cost_theory_workingspace/DATA/mouse/forebrain/proteome/protein_level.txt", row.names = F, quote = F, sep = "\t")








