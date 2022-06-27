#####################
### ARABIDOPSIS ###
#####################

species <- "arabidopsis"
tissue <- "leaves"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, sep="/")
proteome.dir <- paste(main.dir, tissue, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, tissue, "transcriptome", sep="/")
epigenetics.dir <- paste(main.dir, tissue, "epigenetics", sep="/")
noise.dir <- paste(main.dir, "root/single-cell", sep="/")
  
### Proteomic data ### 
# original timeseries data
protein.level <- read.table(paste(proteome.dir, "protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level) <- c("Protein.ID", paste("protein", colnames(protein.level)[-1], sep="_"))
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
protein.genecycle.data <- data.frame(Protein.ID = protein.genecycle.data$Protein.ID,
                                     protein_rhythm.pvalue = protein.genecycle.data$default.pvalue)
# phase detection using Lomb Scargle (LS)
protein.LS.data <- read.table(paste(proteome.dir, "protein_LS.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
protein.LS.data <- data.frame(Protein.ID = protein.LS.data$ID,
                              protein_phase = protein.LS.data$phase)
# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(main.dir, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

proteome.data <- merge(protein.level, protein.genecycle.data, by = "Protein.ID")
proteome.data <- merge(protein.LS.data, proteome.data, by = "Protein.ID")
#proteome.data <- merge(proteome.data, paxdb.leaf, by.x="ID", by.y="Protein.ID", all.x = TRUE)
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by = "Protein.ID", all.x = TRUE)

if (length(proteome.data$Protein.ID) == length(unique(proteome.data$Protein.ID))) {
  print("OK")
} else { print("Warning: protein.IDs are not unique ...") }

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

### Transcriptomic data (Leaves data...) ###
# Gene.ID and Protein.ID :
geneIDs.file <- paste(main.dir, "GeneIDs.txt", sep="/")
geneIDs <- read.table(geneIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
geneIDs <- unique(subset(geneIDs, Protein.ID != ""))
geneIDs <- unique(geneIDs[, c("AFFY.ATH1", "Protein.ID")])

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.level) <- c("ID", paste("RNA", colnames(transcript.level)[-1], sep="_"))
# calculate mean of transcripts levels between all timepoints:
transcript.level$RNA_mean.level <- apply(transcript.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$RNA_max.level <- apply(transcript.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.genecycle.data)[2] <- "RNA_rhythm.pvalue"
# phase detection using Lomb Scargle (LS)
transcript.LS.data <- read.table(paste(transcriptome.dir, "transcript_LS.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
transcript.LS.data <- data.frame(ID = transcript.LS.data$ID,
                                 RNA_phase = transcript.LS.data$phase)

transcriptome.data <- merge(transcript.LS.data, transcript.genecycle.data, by="ID", check.names = FALSE)
transcriptome.data <- merge(transcript.level, transcriptome.data, by="ID", check.names = FALSE)
transcriptome.column.names <- colnames(transcriptome.data)[-1]
transcriptome.data <- merge(geneIDs, transcriptome.data, by.x="AFFY.ATH1", by.y="ID", all.y = TRUE)
#write.table(unique(transcriptome.data[,-2]), paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

TOT.data <- unique(merge(proteome.data, transcriptome.data, by="Protein.ID", all.x = TRUE, all.y = TRUE, check.names = FALSE))

# THEN: We remove AFFY.ATH1 if it is affected to several Protein.ID
raw.dataset.merged <- TOT.data
column.name <- "AFFY.ATH1"
rows.to.remove <- c()
for (i in 1:nrow(raw.dataset.merged)) {
  if ((i %in% rows.to.remove)==FALSE) {
    prob.ID <- raw.dataset.merged[i, column.name]
    if (length(unique(raw.dataset.merged$Protein.ID[grep(prob.ID, raw.dataset.merged[, column.name])])) != 1){
      rows.to.remove <- c(rows.to.remove, grep(prob.ID, raw.dataset.merged[, column.name]))
    }
  }
}

rows.to.remove <- unique(rows.to.remove)
# So, we remove prob.IDs refered to several Gene.IDs:
if (is.null(rows.to.remove)==FALSE){
  raw.dataset.merged <- raw.dataset.merged[-rows.to.remove, ]
}

TOT.data <- raw.dataset.merged
# Remove cases where there are several proteins.IDs for a given gene (to simplify the analysis)
#TOT.data <- TOT.data[!(duplicated(TOT.data$Protein.ID) | duplicated(TOT.data$Protein.ID, fromLast=TRUE)), ]
# Affect NA (to keep all proteins data) into transcripts data to cases where there is 1 prob.ID (Affymetrix probe set from the transcriptomic study) for several proteins.IDs
#TOT.data[duplicated(TOT.data$AFFY.ATH1) | duplicated(TOT.data$AFFY.ATH1, fromLast=TRUE), c(grep("AFFY.ATH1", colnames(TOT.data)) : ncol(TOT.data))] <- NA

#write.csv(TOT.data, paste(main.dir, tissue, "tot_data.csv", sep="/"), quote = F, row.names = F, col.names = colnames(TOT.data))

colnames(TOT.data)[(ncol(TOT.data)-length(transcriptome.column.names)+1) : ncol(TOT.data)] <- transcriptome.column.names


# Epigenetics data
H3K4me3.data <- read.table(paste(epigenetics.dir, "H3K4me3_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
H3K9ac.data <- read.table(paste(epigenetics.dir, "H3K9ac_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(H3K4me3.data)[-1] <- paste("H3K4me3", colnames(H3K4me3.data)[-1], sep="_")
colnames(H3K9ac.data)[-1] <- paste("H3K9ac", colnames(H3K9ac.data)[-1], sep="_")
# rhythm detection using GeneCycle
H3K4me3.genecycle.data <- read.table(paste(epigenetics.dir, "H3K4me3_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(H3K4me3.genecycle.data)[2] <- "H3K4me3_rhythm.pvalue"
H3K9ac.genecycle.data <- read.table(paste(epigenetics.dir, "H3K9ac_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(H3K9ac.genecycle.data)[2] <- "H3K9ac_rhythm.pvalue"
# phase detection using Lomb Scargle (LS)
H3K4me3.LS.data <- read.table(paste(epigenetics.dir, "H3K4me3_LS.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
H3K9ac.LS.data <- read.table(paste(epigenetics.dir, "H3K9ac_LS.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
# tot
H3K4me3.data.tmp <- H3K4me3.genecycle.data
H3K4me3.data.tmp$H3K4me3_phase <- H3K4me3.LS.data$phase
H3K4me3.data <- merge(H3K4me3.data, H3K4me3.data.tmp, by = "ID")

H3K9ac.data.tmp <- H3K9ac.genecycle.data
H3K9ac.data.tmp$H3K9ac_phase <- H3K9ac.LS.data$phase
H3K9ac.data <- merge(H3K9ac.data, H3K9ac.data.tmp, by = "ID")

epigenetic.data <- merge(H3K4me3.data, H3K9ac.data, by = "ID")

## TOT
TOT.data$ID <- gsub("\\..*", "", TOT.data$Protein.ID)
TOT.data_complete <- merge(TOT.data, epigenetic.data, by = "ID", all = TRUE)


# Noise data (F* from Baroso et al.)
noise.data <- read.table(paste(noise.dir, "F*_data.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
noise.data <- noise.data[, c("ID", "Variance", "sd", "noise.adj.sd", "Fano.Factor", "Fstar1", "Fstar4", "Coeff.Variation.squared", "distance.to.median")]

## TOT
TOT.data_complete <- merge(TOT.data_complete, noise.data, by = "ID", all.x = TRUE)

#colnames(TOT.data_complete)
write.csv(TOT.data_complete, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/tot_data.csv", quote = F, row.names = F)
write.table(TOT.data_complete, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/tot_data.txt", sep = "\t", quote = F, row.names = F)





#####################
### CYANOBACTERIA ###
#####################

species <- "cyanobacteria"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

### Proteomic data ### 
# original timeseries data
protein.level <- read.table(paste(proteome.dir, "protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level) <- c("ID", paste("protein", colnames(protein.level)[-1], sep="_"))
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.genecycle.data)[2] <- "protein_rhythm.pvalue"

# Gene.IDs mapping :
geneIDs.file <- paste(main.dir, "UniProt_Ensembl_IDs.txt", sep="/")
geneIDs <- read.table(geneIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
# protein lengh
proteome.data <- merge(protein.level, protein.genecycle.data, by="ID")
proteome.data <- merge(proteome.data, geneIDs, by.x="ID", by.y="UniProt.ID", all.x = TRUE)
colnames(proteome.data)[1] <- "UniProt.ID"

# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(main.dir, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by.x = "UniProt.ID", by.y = "Protein.ID", all.x = TRUE)


if (length(proteome.data$UniProt.ID) == length(unique(proteome.data$UniProt.ID))) {
  print("OK")
} else { print("Warning: protein.IDs are not unique ...") }

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

### Transcriptomic data ###

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.level) <- c("ID", paste("RNA", colnames(transcript.level)[-1], sep="_"))
# calculate mean of transcripts levels between all timepoints:
transcript.level$RNA_mean.level <- apply(transcript.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$RNA_max.level <- apply(transcript.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#transcript.genecycle.data <- transcript.genecycle.data[, c("ID", "default.pvalue")]
colnames(transcript.genecycle.data)[2] <- "RNA_rhythm.pvalue"

transcriptome.data <- merge(transcript.level, transcript.genecycle.data, by="ID", check.names = FALSE)
#write.table(transcriptome.data, paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

TOT.data <- merge(proteome.data, transcriptome.data, by.x="Transcript.ID", by.y="ID", all = TRUE)

write.csv(TOT.data, paste(main.dir, "tot_data.csv", sep="/"), quote = F, row.names = F)
write.table(TOT.data, paste(main.dir, "tot_data.txt", sep="/"), sep = "\t", quote = F, row.names = F)






#####################
### OSTREOCOCCUS ###
#####################

species <- "ostreococcus"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

### Proteomic data ### 
# original timeseries data
protein.level <- read.table(paste(proteome.dir, "protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level) <- c("ID", paste("protein", colnames(protein.level)[-1], sep="_"))
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.genecycle.data)[2] <- "protein_rhythm.pvalue"

# Uniprot.IDs mapping :
uniprotIDs <- read.table(paste(proteome.dir, "uniprot.txt", sep="/"), head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
uniprotIDs <- unique(uniprotIDs[, c("Entry", "Gene.Name")])
# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(main.dir, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa.synth.average.cost.per.protein <- merge(uniprotIDs, aa.synth.average.cost.per.protein, by.x="Entry", by.y="Protein.ID")

proteome.data <- merge(protein.level, protein.genecycle.data, by.x="ID", by.y="Protein.ID")
proteome.data <- proteome.data[proteome.data$ID != "", ]
colnames(proteome.data)[1] <- "Protein.ID"
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by.x="Protein.ID", by.y= "Gene.Name", all.x = TRUE)
proteome.data <- proteome.data[!(duplicated(proteome.data$Protein.ID) | duplicated(proteome.data$Protein.ID, fromLast=TRUE)), ]

if (length(proteome.data$Protein.ID) == length(unique(proteome.data$Protein.ID))) {
  print("OK")
} else { print("Warning: protein.IDs are not unique ...") }

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

### Transcriptomic data ###

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.level) <- c("ID", paste("RNA", colnames(transcript.level)[-1], sep="_"))
column.names <- colnames(transcript.level)
# Transform all normalized values in positives values : 
transcript.level[, -1] <- transcript.level[, -1] + abs(min(transcript.level[, -1], na.rm = TRUE))
transcript.level.raw <- transcript.level
# calculate mean of transcripts levels between all timepoints:
transcript.level$RNA_mean.level <- apply(transcript.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$RNA_max.level <- apply(transcript.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#transcript.genecycle.data <- transcript.genecycle.data[, c("ID", "default.pvalue")]
colnames(transcript.genecycle.data)[2] <- "RNA_rhythm.pvalue"
# Since there are several data per "geneID", we normalize using Brown normalization
library(EmpiricalBrownsMethod)
genes.counts <- plyr::count(transcript.genecycle.data$ID)
genes.in.several.copies <- subset(genes.counts, freq >1)[,"x"]

close <- function(x, value, tol=NULL){
  if(!is.null(tol)){
    x[abs(x-value) <= tol]
  } else {
    x[order(abs(x-value))]
  }
}

original.data <- NULL
pvalues <- NULL
transcript.genecycle.data$pvalue_brown <- NA
for (k in 1:nrow(transcript.genecycle.data)) {
  gene.id <- transcript.genecycle.data$ID[k]
  
  if (gene.id %in% genes.in.several.copies) {
    original.data[[k]] <- transcript.level.raw[grep(gene.id, transcript.level.raw$ID),-1]
    # empiricalBrownsMethod() does not deal with missing values => We replace missing values by the mean of time-points around
    for (l in 1:nrow(original.data[[k]])){
      if ( any(is.na(original.data[[k]][l, ])) ) {
        values.available <- setdiff(1:ncol(original.data[[k]]), which(is.na(original.data[[k]][l, ])))
        for (n in which(is.na(original.data[[k]][l, ])) ) {
          closest.values <- close(values.available, value=n)[1:2]
          original.data[[k]][l, n] <- mean(original.data[[k]][l, closest.values[1]], original.data[[k]][l, closest.values[2]])
        }
      }
    }
    
    pvalues[[k]] <- transcript.genecycle.data[grep(gene.id, transcript.genecycle.data$ID), "RNA_rhythm.pvalue"]
    empirical.browns <- empiricalBrownsMethod(data_matrix=original.data[[k]], p_values=pvalues[[k]], extra_info=TRUE)
    pvalue.brown <- empirical.browns[[1]]
    transcript.genecycle.data$pvalue_brown[k] <- pvalue.brown
  } else {
    transcript.genecycle.data$pvalue_brown[k] <- transcript.genecycle.data$RNA_rhythm.pvalue[k]
  }
}

#write.table(transcript.genecycle.data, paste(transcriptome.dir, "transcript_GeneCycle_brownAggregation.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

transcriptome.data <- unique(transcript.level)
transcriptome.data$RNA_rhythm.pvalue_original <- transcript.genecycle.data$RNA_rhythm.pvalue
transcriptome.data$RNA_rhythm.pvalue_Brown <- transcript.genecycle.data$pvalue_brown
transcriptome.data$RNA_rhythm.pvalue <- transcriptome.data$RNA_rhythm.pvalue_Brown
transcriptome.column.names <- colnames(transcriptome.data)[-1]
#write.table(transcriptome.data, paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

# Gene.IDs mapping :
geneIDs.file <- paste(transcriptome.dir, "transcriptIDs.txt", sep="/")
geneIDs <- read.table(geneIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)

transcriptome.data <- cbind(geneIDs, transcriptome.data[, -1])

# Uniprot.IDs mapping :
uniprotIDs.file <- paste(main.dir, "GeneIDs.tab", sep="/")
uniprotIDs <- read.table(uniprotIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
uniprotIDs <- unique(uniprotIDs[, c("Entry", "Gene.Names")])

proteome.data$Protein.ID <- gsub("OT_ostta", "Ot", proteome.data$Protein.ID)

proteome.data$Uniprot.ID <- NA
for (i in 1:nrow(proteome.data)) {
  protID <- proteome.data[i, "Protein.ID"]
  uniprot.ID <- uniprotIDs[grep(protID, uniprotIDs$Gene.Names), "Entry"]
  if (length(uniprot.ID == 1)) {
    proteome.data$Uniprot.ID[i] <- uniprot.ID
  }
}

TOT.data.1 <- merge(proteome.data, transcriptome.data, by.x="Protein.ID", by.y="PROBES_ID", check.names = FALSE)
TOT.data.2 <- merge(proteome.data, transcriptome.data, by.x="Protein.ID", by.y="CDS_ID_II", check.names = FALSE)
TOT.data <- unique(rbind(TOT.data.1[, -grep("CDS_ID_II", colnames(TOT.data.1))], 
                  TOT.data.2[, -grep("PROBES_ID", colnames(TOT.data.2))]))

# Remove cases where there are several proteins.IDs for a given gene (to simplify the analysis)
TOT.data <- TOT.data[!(duplicated(TOT.data$Protein.ID) | duplicated(TOT.data$Protein.ID, fromLast=TRUE)), ]

colnames(TOT.data)[(ncol(TOT.data)-length(transcriptome.column.names)+1) : ncol(TOT.data)] <- transcriptome.column.names

write.csv(TOT.data, paste(main.dir, "tot_data.csv", sep="/"), quote = F, row.names = F)
write.table(TOT.data, paste(main.dir, "tot_data.txt", sep="/"), sep = "\t", quote = F, row.names = F)







#####################
###### MOUSE ########
#####################
### liver
species <- "mouse"
tissue <- "liver"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

### Proteomic data ### 
protein.level <- read.table(paste(proteome.dir, "raw_protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level)[grep("ZT|LD|CT|LL", colnames(protein.level))] <- paste("protein", grep("ZT|LD|CT|LL", colnames(protein.level), value = TRUE), sep="_")
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, c(-1, -2, -3)], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[, c(-1, -2, -3)], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using Harmonic Regression (Since GeneCycle seemed to not really detect rhythmicity on this dataset..)
protein.harmonic.regression.data <- read.table(paste(proteome.dir, "protein_harmonicRegression.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
protein.harmonic.regression.data <- protein.harmonic.regression.data[, -ncol(protein.harmonic.regression.data)]
colnames(protein.harmonic.regression.data)[ncol(protein.harmonic.regression.data)] <- "protein_rhythm.pvalue"

proteome.data <- protein.level
proteome.data$protein_rhythm.pvalue <- protein.harmonic.regression.data$protein_rhythm.pvalue

### Transcriptomic data ###

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.level)[grep("ZT|LD|CT|LL", colnames(transcript.level))] <- paste("RNA", grep("ZT|LD|CT|LL", colnames(transcript.level), value = TRUE), sep="_")
# calculate mean of transcripts levels between all timepoints:
transcript.level$RNA_mean.level <- apply(transcript.level[, c(1:4)*(-1)], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$RNA_max.level <- apply(transcript.level[, c(1:4)*(-1)], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#transcript.genecycle.data <- transcript.genecycle.data[, c("ID", "default.pvalue")]
colnames(transcript.genecycle.data)[ncol(transcript.genecycle.data)] <- "RNA_rhythm.pvalue"

transcriptome.data <- transcript.level
transcriptome.data$RNA_rhythm.pvalue <- transcript.genecycle.data$RNA_rhythm.pvalue

TOT.data <- cbind(proteome.data, transcriptome.data[, c(-1,-2,-3)])

# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste("~/Documents/cost_theory_workingspace/DATA/", species, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
colnames(aa.synth.average.cost.per.protein)[1] <- "Uniprot.ID"

# If there are several Protein.IDs (Uniprot.IDs) for the same data, we keep only the one "rewieved" by Uniprot
is.integer0 <- function(x) {
  is.integer(x) && length(x) == 0L 
}

TOT.data$protein.length <- NA
TOT.data$aa.synthesis.average.cost_Akashi <- NA
TOT.data$aa.synthesis.average.cost_Wagner <- NA
for (i in 1:nrow(TOT.data)) {
    proteins <- strsplit(as.character(TOT.data[i, "Uniprot.ID"]), split = ";")
    # Search if one or more of them are within Uniprot "REVIEWED" proteins
    match.proteins <- grep(paste(unlist(proteins),collapse="|"), aa.synth.average.cost.per.protein[, "Uniprot.ID"], value=FALSE)
    
    if ( length(match.proteins) == 1 ) {
      TOT.data[i, c("protein.length", "aa.synthesis.average.cost_Akashi", "aa.synthesis.average.cost_Wagner")] <- aa.synth.average.cost.per.protein[match.proteins, c("protein.length", "aa.synthesis.average.cost_Akashi", "aa.synthesis.average.cost_Wagner")]
    } #else if ( is.integer0(match.proteins) ) {} 
  }

genes.unreviewed.by.uniprot <- TOT.data[is.na(TOT.data$protein.length), ]
print(paste("Warning:", nrow(genes.unreviewed.by.uniprot), "genes are being lost (out of", nrow(TOT.data), "genes) [affect NA values in the .csv file]", sep=" "))


#write.table(TOT.data[, c(1:ncol(proteome.data), ncol(TOT.data)-2, ncol(TOT.data)-1, ncol(TOT.data))], 
#            paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)
#write.table(unique(TOT.data[, c((ncol(proteome.data)+1) : (ncol(TOT.data)-3))]), 
#            paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

#write.csv(TOT.data, paste(main.dir, "tot_data.csv", sep="/"), quote = F, row.names = F)
#write.table(TOT.data, paste(main.dir, "tot_data.txt", sep="/"), sep = "\t", quote = F, row.names = F)






#####################
### cartilage
species <- "mouse"
tissue <- "cartilage"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
### Proteomic data ### 
# original timeseries data
protein.level <- read.table(paste(proteome.dir, "protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, c(-1, -2)], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[, c(-1, -2)], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using ARS
protein.ARS.data <- read.table(paste(proteome.dir, "protein_ARS.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.genecycle.data)[2] <- "protein_rhythm.pvalue"

# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste("~/Documents/cost_theory_workingspace/DATA/", species, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
colnames(aa.synth.average.cost.per.protein)[1] <- "Uniprot.ID"

proteome.data <- protein.level
proteome.data$protein_rhythm.pvalue <- protein.ARS.data$default.pvalue
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by="Uniprot.ID", all.x = TRUE)

if (length(proteome.data$Uniprot.ID) == length(unique(proteome.data$Uniprot.ID))) {
  print("OK")
} else { print("Warning: protein.IDs are not unique ...") }

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)


#####################
### tendon
species <- "mouse"
tissue <- "tendon"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
### Proteomic data ### 
# original timeseries data
protein.level <- read.table(paste(proteome.dir, "protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level)[1] <- "Gene.Name"
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.genecycle.data)[2] <- "protein_rhythm.pvalue"

# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste("~/Documents/cost_theory_workingspace/DATA/", species, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
colnames(aa.synth.average.cost.per.protein)[1] <- "Uniprot.ID"

uniprot <- read.table(paste(main.dir, "uniprot_id.txt", sep = "/"), sep = "\t", fill=TRUE, h=T)
uniprot <- uniprot[, c("Entry", "Gene.Name")]
colnames(uniprot) <- c("Uniprot.ID", "Gene.Name")
#inital.data <- read.csv(paste(proteome.dir, "initial_data.csv", sep="/"), head=TRUE, fill=TRUE)

proteome.data <- protein.level
proteome.data$protein_rhythm.pvalue <- protein.genecycle.data$protein_rhythm.pvalue
proteome.data <- proteome.data[proteome.data$Gene.Name != "", ]
proteome.data <- merge(uniprot, proteome.data, by="Gene.Name", all.y = TRUE)
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by="Uniprot.ID", all.x = TRUE)

#proteome.data <- proteome.data[!(duplicated(proteome.data$Protein.ID) | duplicated(proteome.data$Protein.ID, fromLast=TRUE)), ]

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)





#####################
### forebrain
species <- "mouse"
tissue <- "forebrain" 
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
# original timeseries data
protein.level <- read.table(paste(proteome.dir, "protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level)[1] <- "Gene.Name"
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.genecycle.data)[2] <- "protein_rhythm.pvalue"

# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste("~/Documents/cost_theory_workingspace/DATA/", species, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
colnames(aa.synth.average.cost.per.protein)[1] <- "Uniprot.ID"

uniprot <- read.table(paste(main.dir, "uniprot_id.txt", sep = "/"), sep = "\t", fill=TRUE, h=T)
uniprot <- uniprot[, c("Entry", "Gene.Name")]
colnames(uniprot) <- c("Uniprot.ID", "Gene.Name")
#inital.data <- read.csv(paste(proteome.dir, "initial_data.csv", sep="/"), head=TRUE, fill=TRUE)

proteome.data <- protein.level
proteome.data$protein_rhythm.pvalue <- protein.genecycle.data$protein_rhythm.pvalue
proteome.data <- proteome.data[proteome.data$Gene.Name != "", ]
proteome.data <- merge(uniprot, proteome.data, by="Gene.Name", all.y = TRUE)
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by="Uniprot.ID", all.x = TRUE)

#proteome.data <- proteome.data[!(duplicated(proteome.data$Protein.ID) | duplicated(proteome.data$Protein.ID, fromLast=TRUE)), ]

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)










##################################################################
##################################################################
##################################################################
##################################################################
#         idem with z-scores values for expressions
##################################################################
##################################################################
##################################################################
##################################################################

#####################
### ARABIDOPSIS ###
#####################

species <- "arabidopsis"
tissue <- "leaves"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, sep="/")
proteome.dir <- paste(main.dir, tissue, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, tissue, "transcriptome", sep="/")
epigenetics.dir <- paste(main.dir, tissue, "epigenetics", sep="/")

### Proteomic data ### 
# original timeseries data
protein.level <- read.table(paste(proteome.dir, "protein_zscore_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level) <- c("Protein.ID", paste("protein", colnames(protein.level)[-1], sep="_"))
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
protein.genecycle.data <- data.frame(Protein.ID = protein.genecycle.data$Protein.ID,
                                     protein_rhythm.pvalue = protein.genecycle.data$default.pvalue)
# phase detection using Lomb Scargle (LS)
protein.LS.data <- read.table(paste(proteome.dir, "protein_LS.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
protein.LS.data <- data.frame(Protein.ID = protein.LS.data$ID,
                              protein_phase = protein.LS.data$phase)
# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(proteome.dir, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

proteome.data <- merge(protein.level, protein.genecycle.data, by = "Protein.ID")
proteome.data <- merge(protein.LS.data, proteome.data, by = "Protein.ID")
#proteome.data <- merge(proteome.data, paxdb.leaf, by.x="ID", by.y="Protein.ID", all.x = TRUE)
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by = "Protein.ID", all.x = TRUE)

if (length(proteome.data$Protein.ID) == length(unique(proteome.data$Protein.ID))) {
  print("OK")
} else { print("Warning: protein.IDs are not unique ...") }

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

### Transcriptomic data (Leaves data...) ###
# Gene.ID and Protein.ID :
geneIDs.file <- paste(main.dir, "GeneIDs.txt", sep="/")
geneIDs <- read.table(geneIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
geneIDs <- unique(subset(geneIDs, Protein.ID != ""))
geneIDs <- unique(geneIDs[, c("AFFY.ATH1", "Protein.ID")])

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_zscore_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.level) <- c("ID", paste("RNA", colnames(transcript.level)[-1], sep="_"))
# calculate mean of transcripts levels between all timepoints:
transcript.level$RNA_mean.level <- apply(transcript.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$RNA_max.level <- apply(transcript.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.genecycle.data)[2] <- "RNA_rhythm.pvalue"
# phase detection using Lomb Scargle (LS)
transcript.LS.data <- read.table(paste(transcriptome.dir, "transcript_LS.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
transcript.LS.data <- data.frame(ID = transcript.LS.data$ID,
                                 RNA_phase = transcript.LS.data$phase)

transcriptome.data <- merge(transcript.LS.data, transcript.genecycle.data, by="ID", check.names = FALSE)
transcriptome.data <- merge(transcript.level, transcriptome.data, by="ID", check.names = FALSE)
transcriptome.column.names <- colnames(transcriptome.data)[-1]
transcriptome.data <- merge(geneIDs, transcriptome.data, by.x="AFFY.ATH1", by.y="ID", all.y = TRUE)
#write.table(unique(transcriptome.data[,-2]), paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

TOT.data <- unique(merge(proteome.data, transcriptome.data, by="Protein.ID", all.x = TRUE, all.y = TRUE, check.names = FALSE))

# THEN: We remove AFFY.ATH1 if it is affected to several Protein.ID
raw.dataset.merged <- TOT.data
column.name <- "AFFY.ATH1"
rows.to.remove <- c()
for (i in 1:nrow(raw.dataset.merged)) {
  if ((i %in% rows.to.remove)==FALSE) {
    prob.ID <- raw.dataset.merged[i, column.name]
    if (length(unique(raw.dataset.merged$Protein.ID[grep(prob.ID, raw.dataset.merged[, column.name])])) != 1){
      rows.to.remove <- c(rows.to.remove, grep(prob.ID, raw.dataset.merged[, column.name]))
    }
  }
}

rows.to.remove <- unique(rows.to.remove)
# So, we remove prob.IDs refered to several Gene.IDs:
if (is.null(rows.to.remove)==FALSE){
  raw.dataset.merged <- raw.dataset.merged[-rows.to.remove, ]
}

TOT.data <- raw.dataset.merged
# Remove cases where there are several proteins.IDs for a given gene (to simplify the analysis)
#TOT.data <- TOT.data[!(duplicated(TOT.data$Protein.ID) | duplicated(TOT.data$Protein.ID, fromLast=TRUE)), ]
# Affect NA (to keep all proteins data) into transcripts data to cases where there is 1 prob.ID (Affymetrix probe set from the transcriptomic study) for several proteins.IDs
#TOT.data[duplicated(TOT.data$AFFY.ATH1) | duplicated(TOT.data$AFFY.ATH1, fromLast=TRUE), c(grep("AFFY.ATH1", colnames(TOT.data)) : ncol(TOT.data))] <- NA

#write.csv(TOT.data, paste(main.dir, tissue, "tot_data.csv", sep="/"), quote = F, row.names = F, col.names = colnames(TOT.data))

colnames(TOT.data)[(ncol(TOT.data)-length(transcriptome.column.names)+1) : ncol(TOT.data)] <- transcriptome.column.names


# Epigenetics data
H3K4me3.data <- read.table(paste(epigenetics.dir, "H3K4me3_zscore_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
H3K9ac.data <- read.table(paste(epigenetics.dir, "H3K9ac_zscore_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(H3K4me3.data)[-1] <- paste("H3K4me3", colnames(H3K4me3.data)[-1], sep="_")
colnames(H3K9ac.data)[-1] <- paste("H3K9ac", colnames(H3K9ac.data)[-1], sep="_")
# rhythm detection using GeneCycle
H3K4me3.genecycle.data <- read.table(paste(epigenetics.dir, "H3K4me3_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(H3K4me3.genecycle.data)[2] <- "H3K4me3_rhythm.pvalue"
H3K9ac.genecycle.data <- read.table(paste(epigenetics.dir, "H3K9ac_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(H3K9ac.genecycle.data)[2] <- "H3K9ac_rhythm.pvalue"
# phase detection using Lomb Scargle (LS)
H3K4me3.LS.data <- read.table(paste(epigenetics.dir, "H3K4me3_LS.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
H3K9ac.LS.data <- read.table(paste(epigenetics.dir, "H3K9ac_LS.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
# tot
H3K4me3.data.tmp <- H3K4me3.genecycle.data
H3K4me3.data.tmp$H3K4me3_phase <- H3K4me3.LS.data$phase
H3K4me3.data <- merge(H3K4me3.data, H3K4me3.data.tmp, by = "ID")

H3K9ac.data.tmp <- H3K9ac.genecycle.data
H3K9ac.data.tmp$H3K9ac_phase <- H3K9ac.LS.data$phase
H3K9ac.data <- merge(H3K9ac.data, H3K9ac.data.tmp, by = "ID")

epigenetic.data <- merge(H3K4me3.data, H3K9ac.data, by = "ID")

## TOT
TOT.data$ID <- gsub("\\..*", "", TOT.data$Protein.ID)
TOT.data_complete <- merge(TOT.data, epigenetic.data, by = "ID", all = TRUE)


# Noise data (F* from Baroso et al.)
noise.data <- read.table(paste(main.dir, "F*_data.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
noise.data <- noise.data[, c("ID", "Fstar", "Variance", "Noise")]

## TOT
TOT.data_complete <- merge(TOT.data_complete, noise.data, by = "ID", all.x = TRUE)

#colnames(TOT.data_complete)
write.csv(TOT.data_complete, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/tot_data_zscore.csv", quote = F, row.names = F)
write.table(TOT.data_complete, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/tot_data_zscore.txt", sep = "\t", quote = F, row.names = F)





#####################
### CYANOBACTERIA ###
#####################

species <- "cyanobacteria"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

### Proteomic data ### 
# original timeseries data
protein.level <- read.table(paste(proteome.dir, "raw_protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level) <- c("ID", paste("protein", colnames(protein.level)[-1], sep="_"))
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.genecycle.data)[2] <- "protein_rhythm.pvalue"

# Gene.IDs mapping :
geneIDs.file <- paste(main.dir, "UniProt_Ensembl_IDs.txt", sep="/")
geneIDs <- read.table(geneIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(proteome.dir, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

proteome.data <- merge(protein.level, protein.genecycle.data, by="ID")
proteome.data <- merge(proteome.data, geneIDs, by.x="ID", by.y="UniProt.ID", all.x = TRUE)
colnames(proteome.data)[1] <- "UniProt.ID"
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by="Ensembl.Genome.Protein.ID", all.x = TRUE)


if (length(proteome.data$UniProt.ID) == length(unique(proteome.data$UniProt.ID))) {
  print("OK")
} else { print("Warning: protein.IDs are not unique ...") }

write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

### Transcriptomic data ###

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_zscore_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.level) <- c("ID", paste("RNA", colnames(transcript.level)[-1], sep="_"))
# calculate mean of transcripts levels between all timepoints:
transcript.level$RNA_mean.level <- apply(transcript.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$RNA_max.level <- apply(transcript.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#transcript.genecycle.data <- transcript.genecycle.data[, c("ID", "default.pvalue")]
colnames(transcript.genecycle.data)[2] <- "RNA_rhythm.pvalue"

transcriptome.data <- merge(transcript.level, transcript.genecycle.data, by="ID", check.names = FALSE)
#write.table(transcriptome.data, paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

TOT.data <- merge(proteome.data, transcriptome.data, by.x="Transcript.ID", by.y="ID", all = TRUE)

write.csv(TOT.data, paste(main.dir, "tot_data_zscore.csv", sep="/"), quote = F, row.names = F)
write.table(TOT.data, paste(main.dir, "tot_data_zscore.txt", sep="/"), sep = "\t", quote = F, row.names = F)






#####################
### OSTREOCOCCUS ###
#####################

species <- "ostreococcus"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

### Proteomic data ### 
# original timeseries data
protein.level <- read.table(paste(proteome.dir, "protein_zscore_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level) <- c("ID", paste("protein", colnames(protein.level)[-1], sep="_"))
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.genecycle.data)[2] <- "protein_rhythm.pvalue"

# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(proteome.dir, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa.synth.average.cost.per.protein <- unique(aa.synth.average.cost.per.protein)
# protein abundance data from PaxDB
#paxdb <- read.table(paste(proteome.dir, "/", species, "_paxdb.txt", sep=""), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#colnames(paxdb)[2] <- "paxdb.protein.abundance"
# since there is several abundances for a same protein (gene for ostreococcus), we aggregate them by the mean
#paxdb <- aggregate(paxdb, by = list(paxdb$Protein.ID), FUN = mean)[-1, c(1,3)]
#colnames(paxdb)[1] <- "Protein.ID"

proteome.data <- merge(protein.level, protein.genecycle.data, by.x="ID", by.y="Protein.ID")
proteome.data <- proteome.data[proteome.data$ID != "", ]
#proteome.data <- merge(proteome.data, paxdb, by.x="ID", by.y="Protein.ID", all.x = TRUE)
colnames(proteome.data)[1] <- "Protein.ID"
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by="Protein.ID", all.x = TRUE)
proteome.data <- proteome.data[!(duplicated(proteome.data$Protein.ID) | duplicated(proteome.data$Protein.ID, fromLast=TRUE)), ]

if (length(proteome.data$Protein.ID) == length(unique(proteome.data$Protein.ID))) {
  print("OK")
} else { print("Warning: protein.IDs are not unique ...") }

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

### Transcriptomic data ###

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_zscore_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.level) <- c("ID", paste("RNA", colnames(transcript.level)[-1], sep="_"))
column.names <- colnames(transcript.level)
# Transform all normalized values in positives values : 
transcript.level[, -1] <- transcript.level[, -1] + abs(min(transcript.level[, -1], na.rm = TRUE))
transcript.level.raw <- transcript.level
# calculate mean of transcripts levels between all timepoints:
transcript.level$RNA_mean.level <- apply(transcript.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$RNA_max.level <- apply(transcript.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#transcript.genecycle.data <- transcript.genecycle.data[, c("ID", "default.pvalue")]
colnames(transcript.genecycle.data)[2] <- "RNA_rhythm.pvalue"
# Since there are several data per "geneID", we normalize using Brown normalization
library(EmpiricalBrownsMethod)
genes.counts <- plyr::count(transcript.genecycle.data$ID)
genes.in.several.copies <- subset(genes.counts, freq >1)[,"x"]

close <- function(x, value, tol=NULL){
  if(!is.null(tol)){
    x[abs(x-value) <= tol]
  } else {
    x[order(abs(x-value))]
  }
}

original.data <- NULL
pvalues <- NULL
transcript.genecycle.data$pvalue_brown <- NA
for (k in 1:nrow(transcript.genecycle.data)) {
  gene.id <- transcript.genecycle.data$ID[k]
  
  if (gene.id %in% genes.in.several.copies) {
    original.data[[k]] <- transcript.level.raw[grep(gene.id, transcript.level.raw$ID),-1]
    # empiricalBrownsMethod() does not deal with missing values => We replace missing values by the mean of time-points around
    for (l in 1:nrow(original.data[[k]])){
      if ( any(is.na(original.data[[k]][l, ])) ) {
        values.available <- setdiff(1:ncol(original.data[[k]]), which(is.na(original.data[[k]][l, ])))
        for (n in which(is.na(original.data[[k]][l, ])) ) {
          closest.values <- close(values.available, value=n)[1:2]
          original.data[[k]][l, n] <- mean(original.data[[k]][l, closest.values[1]], original.data[[k]][l, closest.values[2]])
        }
      }
    }
    
    pvalues[[k]] <- transcript.genecycle.data[grep(gene.id, transcript.genecycle.data$ID), "RNA_rhythm.pvalue"]
    empirical.browns <- empiricalBrownsMethod(data_matrix=original.data[[k]], p_values=pvalues[[k]], extra_info=TRUE)
    pvalue.brown <- empirical.browns[[1]]
    transcript.genecycle.data$pvalue_brown[k] <- pvalue.brown
  } else {
    transcript.genecycle.data$pvalue_brown[k] <- transcript.genecycle.data$RNA_rhythm.pvalue[k]
  }
}

#write.table(transcript.genecycle.data, paste(transcriptome.dir, "transcript_GeneCycle_brownAggregation.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

transcriptome.data <- unique(transcript.level)
transcriptome.data$RNA_rhythm.pvalue_original <- transcript.genecycle.data$RNA_rhythm.pvalue
transcriptome.data$RNA_rhythm.pvalue_Brown <- transcript.genecycle.data$pvalue_brown
transcriptome.column.names <- colnames(transcriptome.data)[-1]
#write.table(transcriptome.data, paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

# Gene.IDs mapping :
geneIDs.file <- paste(transcriptome.dir, "transcriptIDs.txt", sep="/")
geneIDs <- read.table(geneIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)

transcriptome.data <- cbind(geneIDs, transcriptome.data[, -1])

# Uniprot.IDs mapping :
uniprotIDs.file <- paste(main.dir, "GeneIDs.tab", sep="/")
uniprotIDs <- read.table(uniprotIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
uniprotIDs <- unique(uniprotIDs[, c("Entry", "Gene.Names")])

proteome.data$Protein.ID <- gsub("OT_ostta", "Ot", proteome.data$Protein.ID)

proteome.data$Uniprot.ID <- NA
for (i in 1:nrow(proteome.data)) {
  protID <- proteome.data[i, "Protein.ID"]
  uniprot.ID <- uniprotIDs[grep(protID, uniprotIDs$Gene.Names), "Entry"]
  if (length(uniprot.ID == 1)) {
    proteome.data$Uniprot.ID[i] <- uniprot.ID
  }
}

TOT.data.1 <- merge(proteome.data, transcriptome.data, by.x="Protein.ID", by.y="PROBES_ID", check.names = FALSE)
TOT.data.2 <- merge(proteome.data, transcriptome.data, by.x="Protein.ID", by.y="CDS_ID_II", check.names = FALSE)
TOT.data <- unique(rbind(TOT.data.1[, -grep("CDS_ID_II", colnames(TOT.data.1))], 
                         TOT.data.2[, -grep("PROBES_ID", colnames(TOT.data.2))]))

# Remove cases where there are several proteins.IDs for a given gene (to simplify the analysis)
TOT.data <- TOT.data[!(duplicated(TOT.data$Protein.ID) | duplicated(TOT.data$Protein.ID, fromLast=TRUE)), ]

colnames(TOT.data)[(ncol(TOT.data)-length(transcriptome.column.names)+1) : ncol(TOT.data)] <- transcriptome.column.names

write.csv(TOT.data, paste(main.dir, "tot_data_zscore.csv", sep="/"), quote = F, row.names = F)
write.table(TOT.data, paste(main.dir, "tot_data_zscore.txt", sep="/"), sep = "\t", quote = F, row.names = F)







#####################
###### MOUSE ########
#####################
### liver
species <- "mouse"
tissue <- "liver"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

### Proteomic data ###
protein.level <- read.table(paste(proteome.dir, "protein_zscore_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.level)[grep("ZT|LD|CT|LL", colnames(protein.level))] <- paste("protein", grep("ZT|LD|CT|LL", colnames(protein.level), value = TRUE), sep="_")
# calculate mean of protein levels between all timepoints:
protein.level$protein_mean.level <- apply(protein.level[, c(-1, -2, -3)], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$protein_max.level <- apply(protein.level[, c(-1, -2, -3)], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using Harmonic Regression (Since GeneCycle seemed to not really detect rhythmicity on this dataset..)
protein.harmonic.regression.data <- read.table(paste(proteome.dir, "protein_harmonicRegression.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
protein.harmonic.regression.data <- protein.harmonic.regression.data[, -ncol(protein.harmonic.regression.data)]
colnames(protein.harmonic.regression.data)[ncol(protein.harmonic.regression.data)] <- "protein_rhythm.pvalue"

proteome.data <- protein.level
proteome.data$protein_rhythm.pvalue <- protein.harmonic.regression.data$protein_rhythm.pvalue

### Transcriptomic data ###

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_zscore_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(transcript.level)[grep("ZT|LD|CT|LL", colnames(transcript.level))] <- paste("RNA", grep("ZT|LD|CT|LL", colnames(transcript.level), value = TRUE), sep="_")
# calculate mean of transcripts levels between all timepoints:
transcript.level$RNA_mean.level <- apply(transcript.level[, c(1:4)*(-1)], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$RNA_max.level <- apply(transcript.level[, c(1:4)*(-1)], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#transcript.genecycle.data <- transcript.genecycle.data[, c("ID", "default.pvalue")]
colnames(transcript.genecycle.data)[ncol(transcript.genecycle.data)] <- "RNA_rhythm.pvalue"

transcriptome.data <- transcript.level
transcriptome.data$RNA_rhythm.pvalue <- transcript.genecycle.data$RNA_rhythm.pvalue

TOT.data <- cbind(proteome.data, transcriptome.data[, c(-1,-2,-3)])

# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(proteome.dir, "/", species, "_UniprotReviewed_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa.synth.average.cost.per.protein <- unique(aa.synth.average.cost.per.protein)

# If there are several Protein.IDs (Uniprot.IDs) for the same data, we keep only the one "rewieved" by Uniprot
is.integer0 <- function(x) {
  is.integer(x) && length(x) == 0L 
}

TOT.data$protein.length <- NA
TOT.data$aa.synthesis.average.cost_Akashi <- NA
TOT.data$aa.synthesis.average.cost_Wagner <- NA
TOT.data$Mass <- NA
for (i in 1:nrow(TOT.data)) {
  proteins <- strsplit(as.character(TOT.data[i, "Uniprot.ID"]), split = ";")
  # Search if one or more of them are within Uniprot "REVIEWED" proteins
  match.proteins <- grep(paste(unlist(proteins),collapse="|"), aa.synth.average.cost.per.protein[, "Uniprot.ID"], value=FALSE)
  
  if ( length(match.proteins) == 1 ) {
    TOT.data[i, c("protein.length", "aa.synthesis.average.cost_Akashi", "aa.synthesis.average.cost_Wagner", "Mass")] <- aa.synth.average.cost.per.protein[match.proteins, c("protein.length", "aa.synthesis.average.cost_Akashi", "aa.synthesis.average.cost_Wagner", "Mass")]
  } #else if ( is.integer0(match.proteins) ) {} 
}

genes.unreviewed.by.uniprot <- TOT.data[is.na(TOT.data$protein.length), ]
print(paste("Warning:", nrow(genes.unreviewed.by.uniprot), "genes are being lost (out of", nrow(TOT.data), "genes) [affect NA values in the .csv file]", sep=" "))

# protein abundance data from PaxDB: lets ignore these data otherwize we would need to retreive mapping between Uniprot.IDs and Ensembl.IDs for proteins leading to loose information 
#paxdb <- read.table(paste(proteome.dir, "/", species, "_", tissue, "_paxdb.txt", sep=""), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#colnames(paxdb)[2] <- "paxdb.protein.abundance"


#write.table(TOT.data[, c(2:ncol(proteome.data), ncol(TOT.data)-2, ncol(TOT.data)-1, ncol(TOT.data))], 
#            paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)
#write.table(unique(TOT.data[, c((ncol(proteome.data)+1) : (ncol(TOT.data)-3))]), 
#            paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

write.csv(TOT.data, paste(main.dir, "tot_data_zscore.csv", sep="/"), quote = F, row.names = F)
write.table(TOT.data, paste(main.dir, "tot_data_zscore.txt", sep="/"), sep = "\t", quote = F, row.names = F)







