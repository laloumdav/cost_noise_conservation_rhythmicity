#####################
### ARABIDOPSIS ###
#####################

species <- "arabidopsis"
tissue <- "leaves"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/")

tot.data <- read.csv(paste(main.dir,"tot_data.csv", sep="/"), head=TRUE)
data.nb <- nrow(tot.data)
tot.data <- na.omit(tot.data)
print(paste(data.nb-nrow(tot.data), "over", nrow(tot.data), "genes data are being lost", sep=" "))

### transcripts ###
transcriptome.data <- tot.data[, (1+grep("aa.synthesis.average.cost", colnames(tot.data))):ncol(tot.data)]
# retreive replicates for each time-point
timepoint.names <- paste("LD", rep(c("00", "04", "08", "12", "16", "20"), each=3), sep="")
colnames(transcriptome.data)[grep("LD", colnames(transcriptome.data))] <- timepoint.names

col.names.IDs.nb = 1

t.transcriptome.data <- as.data.frame(t(transcriptome.data[, col.names.IDs.nb*-1]))
t.transcriptome.data$col.names <- colnames(transcriptome.data)[col.names.IDs.nb*-1]
t.transcriptome.data <- aggregate(.~col.names, data = t.transcriptome.data, FUN="mean")

transcriptome.data.aggregated <- as.data.frame(t(t.transcriptome.data[, col.names.IDs.nb*-1]))
colnames(transcriptome.data.aggregated) <- t.transcriptome.data$col.names
transcriptome.data.aggregated[, colnames(transcriptome.data)[col.names.IDs.nb]] <- transcriptome.data[, col.names.IDs.nb]
# re-order columns 
transcriptome.data.aggregated <- transcriptome.data.aggregated[, c(ncol(transcriptome.data.aggregated), 1:(ncol(transcriptome.data.aggregated)-1))]

### proteins ###
proteome.data <- tot.data[, 1:grep("aa.synthesis.average.cost", colnames(tot.data))]
# retreive replicates for each time-point
timepoint.names <- paste("LL", rep(c("12", "16", "20", "00", "04", "08"), 5), sep="")
colnames(proteome.data)[grep("LL", colnames(proteome.data))] <- timepoint.names

col.names.IDs.nb = 1

t.proteome.data <- as.data.frame(t(proteome.data[, col.names.IDs.nb*-1]))
t.proteome.data$col.names <- colnames(proteome.data)[col.names.IDs.nb*-1]
t.proteome.data <- aggregate(.~col.names, data = t.proteome.data, FUN="mean")

proteome.data.aggregated <- as.data.frame(t(t.proteome.data[, col.names.IDs.nb*-1]))
colnames(proteome.data.aggregated) <- t.proteome.data$col.names
proteome.data.aggregated[, colnames(proteome.data)[col.names.IDs.nb]] <- proteome.data[, col.names.IDs.nb]
# re-order columns 
proteome.data.aggregated <- proteome.data.aggregated[, c(ncol(proteome.data.aggregated), 1:(ncol(proteome.data.aggregated)-1))]


### at which time-point the protein is max expressed ?
timePointMaxExpr <- function(x) {
  return( min(grep(max(x), x)) )
}

timepoint.names <- colnames(transcriptome.data.aggregated)[grep("LD", colnames(transcriptome.data.aggregated))]
time.point.nb <- apply(transcriptome.data.aggregated[, timepoint.names], MARGIN = 1, FUN = timePointMaxExpr)
transcriptome.data.aggregated$RNA.timepoint.max.expr <- timepoint.names[time.point.nb]

timepoint.names <- colnames(proteome.data.aggregated)[grep("LL", colnames(proteome.data.aggregated))]
time.point.nb <- apply(proteome.data.aggregated[, timepoint.names], MARGIN = 1, FUN = timePointMaxExpr)
proteome.data.aggregated$protein.timepoint.max.expr <- timepoint.names[time.point.nb]

tot.data <- cbind(proteome.data.aggregated, transcriptome.data.aggregated)

subset.rhythmic.RNA <- subset(tot.data$max.protein.level, RNA.rhythm.pvalue <=0.01 
                              & protein.rhythm.pvalue>0.1 
                              & )

tot.data$approx.Q = (tot.data$aa.synthesis.average.cost * tot.data$max.protein.level) / tot.data$mean.RNA.level
  
rhythmic <- subset(tot.data, protein.rhythm.pvalue <= 0.01)
non.rhythmic <- subset(tot.data, protein.rhythm.pvalue > 0.2)
non.rhythmic <- non.rhythmic[sample(1:nrow(non.rhythmic), replace = F, size = nrow(rhythmic)), ]

wilcox.test(rhythmic$approx.Q, non.rhythmic$approx.Q)

timepoint.RNA <- "LD08"
timepoint.protein <- "LL08"
subset.rhythmic.prot <- subset(tot.data, protein.rhythm.pvalue <=0.01
                               & RNA.rhythm.pvalue >0.5
                               & protein.timepoint.max.expr == timepoint.protein)
plot(log10(tot.data[[timepoint.RNA]]), log10(tot.data[[timepoint.protein]])/log10(tot.data[[timepoint.RNA]]), xlim = c(0, 5), ylim = c(1,9))
par(new=TRUE)
plot(log10(subset.rhythmic.prot[[timepoint.RNA]]), log10(subset.rhythmic.prot[[timepoint.protein]])/log10(subset.rhythmic.prot[[timepoint.RNA]]), 
     col="blue", xlim = c(0, 5), ylim = c(1,9))


plot(log10(tot.data[[timepoint.RNA]]), log10(tot.data[[timepoint.protein]]), xlim = c(0, 5), ylim = c(1,9))
par(new=TRUE)
plot(log10(subset.rhythmic.prot[[timepoint.RNA]]), log10(subset.rhythmic.prot[[timepoint.protein]]), col="blue", xlim = c(0, 5), ylim = c(1,9))


par(new=TRUE)
plot(log(subset.rhythmic.RNA$LL04), log(subset.rhythmic.RNA$LD04), col="red")












### Transcriptomic data (Leaves data...) ###
# Gene.ID and Protein.ID :
geneIDs.file <- paste(main.dir, "GeneIDs.txt", sep="/")
geneIDs <- read.table(geneIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
geneIDs <- unique(subset(geneIDs, Protein.ID != ""))
geneIDs <- unique(geneIDs[, c("AFFY.ATH1", "Protein.ID")])


groups <- as.factor(rbinom(32, n = 5, prob = 0.4))
tapply(groups, groups, length) #- is almost the same as
table(groups)
tapply()


install.packages("remotes")
remotes::install_github("brooksandrew/Rsenal")
library(Rse)
require('data.table')
## establishing variables to aggregate on
lengthVs <- c('Sepal.Length', 'Petal.Length')
widthVs <- c('Sepal.Width', 'Petal.Width')
## aggregating using 2 different functions and identifying columns to aggregate by variable names
irisAgg1 <- smartAgg(df=iris, by='Species', 'mean', lengthVs, 'sum', widthVs)
## aggregating using 2 dimensions ("Specied" and "randthing")
iris$randthing <- as.character(sample(1:5, nrow(iris), replace=T))
irisAgg2 <- smartAgg(df=iris, by=c('Species', 'randthing'), 'mean', lengthVs, 'sum', widthVs, catN=T, printAgg=T)
## aggregating variables by column number
irisAgg3 <- smartAgg(df=iris, by=c('Species', 'randthing'), 'mean', 1:2, 'sum', 3:4, catN=T, printAgg=T)
## use anonymous functions
data(mtcars)
smartAgg(mtcars, by='cyl', function(x) sum(x*100), c('drat', 'mpg', 'disp'))
## use anonymous functions with more than 1 argument.  Uses the provided variables for all unassigned arguments in anonymous function
smartAgg(mtcars, by='cyl', function(x,y='carb') sum(x*y), c('drat', 'mpg', 'disp'))
with(mtcars[mtcars$cyl==6,], c(sum(drat*carb), sum(mpg*carb), sum(disp*carb)))
## with anonymous functions with more than 1 argument.
## Example of possible unintended behavior - the user-provided variable is used for both and x and y in this example.
smartAgg(mtcars, by='cyl', function(x,y) sum(x*y), c('drat', 'mpg', 'disp'))
with(mtcars[mtcars$cyl==6,], c(sum(drat*drat), sum(mpg*mpg), sum(carb*carb)))







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
# calculate mean of protein levels between all timepoints:
protein.level$mean.protein.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$max.protein.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.genecycle.data)[2] <- "protein.rhythm.pvalue"

# Gene.IDs mapping :
geneIDs.file <- paste(main.dir, "UniProt_Ensembl_IDs.txt", sep="/")
geneIDs <- read.table(geneIDs.file, head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(proteome.dir, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa.synth.average.cost.per.protein <- merge(geneIDs, aa.synth.average.cost.per.protein, by.x = "Ensembl.Genome.Protein.ID", by.y = "Protein.ID")
aa.synth.average.cost.per.protein <- aa.synth.average.cost.per.protein[, -1]

proteome.data <- merge(protein.level, protein.genecycle.data, by="ID")
colnames(proteome.data)[1] <- "UniProt.ID"
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by="UniProt.ID", all.y = TRUE)


if (length(proteome.data$UniProt.ID) == length(unique(proteome.data$UniProt.ID))) {
  print("OK")
} else { print("Warning: protein.IDs are not unique ...") }

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

### Transcriptomic data ###

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
# calculate mean of transcripts levels between all timepoints:
transcript.level$mean.RNA.level <- apply(transcript.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$max.RNA.level <- apply(transcript.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#transcript.genecycle.data <- transcript.genecycle.data[, c("ID", "default.pvalue")]
colnames(transcript.genecycle.data)[2] <- "RNA.rhythm.pvalue"

transcriptome.data <- merge(transcript.level, transcript.genecycle.data, by="ID", check.names = FALSE)
#write.table(transcriptome.data, paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

TOT.data <- merge(proteome.data, transcriptome.data, by.x="Transcript.ID", by.y="ID", all.x = TRUE)
length(unique(TOT.data$UniProt.ID))

#write.csv(TOT.data, paste(main.dir, "tot_data.csv", sep="/"), quote = F, row.names = F)





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
# calculate mean of protein levels between all timepoints:
protein.level$mean.protein.level <- apply(protein.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$max.protein.level <- apply(protein.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
protein.genecycle.data <- read.table(paste(proteome.dir, "protein_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(protein.genecycle.data)[2] <- "protein.rhythm.pvalue"

# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(proteome.dir, "/", species, "_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa.synth.average.cost.per.protein <- unique(aa.synth.average.cost.per.protein)
# protein abundance data from PaxDB
paxdb <- read.table(paste(proteome.dir, "/", species, "_paxdb.txt", sep=""), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
colnames(paxdb)[2] <- "paxdb.protein.abundance"
# since there is several abundances for a same protein (gene for ostreococcus), we aggregate them by the mean
paxdb <- aggregate(paxdb, by = list(paxdb$Protein.ID), FUN = mean)[-1, c(1,3)]
colnames(paxdb)[1] <- "Protein.ID"

proteome.data <- merge(protein.level, protein.genecycle.data, by.x="ID", by.y="Protein.ID")
proteome.data <- proteome.data[proteome.data$ID != "", ]
proteome.data <- merge(proteome.data, paxdb, by.x="ID", by.y="Protein.ID", all.x = TRUE)
colnames(proteome.data)[1] <- "Protein.ID"
proteome.data <- merge(proteome.data, aa.synth.average.cost.per.protein, by="Protein.ID", all.x = TRUE)
proteome.data <- proteome.data[!(duplicated(proteome.data$Protein.ID) | duplicated(proteome.data$Protein.ID, fromLast=TRUE)), ]

if (length(proteome.data$Protein.ID) == length(unique(proteome.data$Protein.ID))) {
  print("OK")
} else { print("Warning: protein.IDs are not unique ...") }

#write.table(proteome.data, paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

### Transcriptomic data ###

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
column.names <- colnames(transcript.level)
# Transform all normalized values in positives values : 
transcript.level[, -1] <- transcript.level[, -1] + abs(min(transcript.level[, -1], na.rm = TRUE))
transcript.level.raw <- transcript.level
# calculate mean of transcripts levels between all timepoints:
transcript.level$mean.RNA.level <- apply(transcript.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$max.RNA.level <- apply(transcript.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#transcript.genecycle.data <- transcript.genecycle.data[, c("ID", "default.pvalue")]
colnames(transcript.genecycle.data)[2] <- "RNA.rhythm.pvalue"
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
    
    pvalues[[k]] <- transcript.genecycle.data[grep(gene.id, transcript.genecycle.data$ID), "RNA.rhythm.pvalue"]
    empirical.browns <- empiricalBrownsMethod(data_matrix=original.data[[k]], p_values=pvalues[[k]], extra_info=TRUE)
    pvalue.brown <- empirical.browns[[1]]
    transcript.genecycle.data$pvalue_brown[k] <- pvalue.brown
  } else {
    transcript.genecycle.data$pvalue_brown[k] <- transcript.genecycle.data$RNA.rhythm.pvalue[k]
  }
}

#write.table(transcript.genecycle.data, paste(transcriptome.dir, "transcript_GeneCycle_brownAggregation.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

transcriptome.data <- unique(transcript.level)
transcriptome.data$RNA.rhythm.pvalue_original <- transcript.genecycle.data$RNA.rhythm.pvalue
transcriptome.data$RNA.rhythm.pvalue_Brown <- transcript.genecycle.data$pvalue_brown
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

#write.csv(TOT.data, paste(main.dir, "tot_data.csv", sep="/"), quote = F, row.names = F)








#####################
###### MOUSE ########
#####################

species <- "mouse"
tissue <- "liver"
main.dir <- paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/")
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

### Proteomic data ### 
# original timeseries data (raw proteins level values)
protein.level <- read.table(paste(proteome.dir, "raw_protein_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
# calculate mean of protein levels between all timepoints:
protein.level$mean.protein.level <- apply(protein.level[, c(-1, -2, -3)], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of protein levels in all timepoints:
protein.level$max.protein.level <- apply(protein.level[, c(-1, -2, -3)], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using Harmonic Regression (Since GeneCycle seemed to not really detect rhythmicity on this dataset..)
protein.harmonic.regression.data <- read.table(paste(proteome.dir, "protein_harmonicRegression.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
protein.harmonic.regression.data <- protein.harmonic.regression.data[, -ncol(protein.harmonic.regression.data)]
colnames(protein.harmonic.regression.data)[ncol(protein.harmonic.regression.data)] <- "protein.rhythm.pvalue"

proteome.data <- protein.level
proteome.data$protein.rhythm.pvalue <- protein.harmonic.regression.data$protein.rhythm.pvalue

### Transcriptomic data ###

# original timeseries data
transcript.level <- read.table(paste(transcriptome.dir, "transcript_level.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
# calculate mean of transcripts levels between all timepoints:
transcript.level$mean.RNA.level <- apply(transcript.level[, c(1:4)*(-1)], 1, FUN = function(x){mean(x, na.rm = TRUE)})
# calculate mean of the 2 maximum values of transcripts levels in all timepoints:
transcript.level$max.RNA.level <- apply(transcript.level[, c(1:4)*(-1)], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
# rhythm detection using GeneCycle
transcript.genecycle.data <- read.table(paste(transcriptome.dir, "transcript_GeneCycle.txt", sep="/"), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#transcript.genecycle.data <- transcript.genecycle.data[, c("ID", "default.pvalue")]
colnames(transcript.genecycle.data)[ncol(transcript.genecycle.data)] <- "RNA.rhythm.pvalue"

transcriptome.data <- transcript.level
transcriptome.data$RNA.rhythm.pvalue <- transcript.genecycle.data$RNA.rhythm.pvalue

TOT.data <- cbind(proteome.data, transcriptome.data[, c(-1,-2,-3)])

# protein lenght and AA synthesis average cost calculated for each protein :
aa.synth.average.cost.per.protein <- read.table(paste(proteome.dir, "/", species, "_UniprotReviewed_AA_synthesis_average_cost_per_protein.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
aa.synth.average.cost.per.protein <- unique(aa.synth.average.cost.per.protein)

# If there are several Protein.IDs (Uniprot.IDs) for the same data, we keep only the one "rewieved" by Uniprot
is.integer0 <- function(x) {
  is.integer(x) && length(x) == 0L 
}

TOT.data$protein.length <- NA
TOT.data$aa.synthesis.average.cost <- NA
for (i in 1:nrow(TOT.data)) {
  proteins <- strsplit(as.character(TOT.data[i, "Uniprot.ID"]), split = ";")
  # Search if one or more of them are within Uniprot "REVIEWED" proteins
  match.proteins <- grep(paste(unlist(proteins),collapse="|"), aa.synth.average.cost.per.protein[, "Uniprot.ID"], value=FALSE)
  
  if ( length(match.proteins) == 1 ) {
    TOT.data[i, c("protein.length", "aa.synthesis.average.cost")] <- aa.synth.average.cost.per.protein[match.proteins, c("protein.length", "aa.synthesis.average.cost")]
  } #else if ( is.integer0(match.proteins) ) {} 
}

genes.unreviewed.by.uniprot <- TOT.data[is.na(TOT.data$protein.length), ]
print(paste("Warning:", nrow(genes.unreviewed.by.uniprot), "genes are being lost (out of", nrow(TOT.data), "genes) [affect NA values in the .csv file]", sep=" "))

# protein abundance data from PaxDB: lets ignore these data otherwize we would need to retreive mapping between Uniprot.IDs and Ensembl.IDs for proteins leading to loose information 
#paxdb <- read.table(paste(proteome.dir, "/", species, "_", tissue, "_paxdb.txt", sep=""), head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
#colnames(paxdb)[2] <- "paxdb.protein.abundance"


#write.table(TOT.data[, c(2:ncol(proteome.data), ncol(TOT.data)-1, ncol(TOT.data))], 
#            paste(proteome.dir, "proteome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)
#write.table(unique(TOT.data[, c((ncol(proteome.data)+1) : (ncol(TOT.data)-2))]), 
#            paste(transcriptome.dir, "transcriptome_data.txt", sep = "/"), sep = "\t", quote = F, row.names = F)

#write.csv(TOT.data, paste(main.dir, "tot_data.csv", sep="/"), quote = F, row.names = F)

colnames(proteome.data)
colnames(transcriptome.data)







# Total cost for a steady state protein level
#protein.abundance.and.cost$tot.cost.ss.protein.level <- protein.abundance.and.cost$protein.abundance * protein.abundance.and.cost$protein.length * (protein.abundance.and.cost$aa.synthesis.average.cost - 1)



