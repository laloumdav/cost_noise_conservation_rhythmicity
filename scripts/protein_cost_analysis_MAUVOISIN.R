library(ggplot2)

main.dir <- "~/Documents/cost_theory_workingspace/DATA/mauvoisin/"
file.dir <- main.dir
species <- "mouse"
tissue <- "liver"

### SCRIPT ###
file.name <- paste(file.dir, "protein_level.txt", sep="")
raw.dataset <- read.table(file.name, head=TRUE, fill=TRUE, sep = "\t")

# Gene.ID and Protein.ID : 
species.genesID.protID.file <- paste(file.dir, "GeneID_GeneName_TranscriptID_ProtID_MOUSE.txt", sep="")
geneID.proteinID <- read.table(species.genesID.protID.file, head=TRUE, fill = TRUE)
geneID.proteinID <- unique(subset(geneID.proteinID, Protein.ID != ""))
geneID.proteinID$Protein.ID <- gsub("\\..*", "", geneID.proteinID$Protein.ID)

# Circad vs Not-Circad proteins:
circad.tissue.dir <- paste(file.dir, "protein_GeneCycle.txt", sep="")
circad.proteins <- read.table(circad.tissue.dir, head=TRUE, fill=TRUE, check.names = FALSE, sep = "\t", stringsAsFactors = FALSE)

random.circad.proteins <- circad.proteins[sample(1:nrow(circad.proteins), size = 400, replace = F), ]
nrow(random.circad.proteins[random.circad.proteins$default.pvalue<=0.005,])

# Expression level per transcript in the dataset of the circadian experiment where circadian data come from :
#express.level.per.transcript.file <- paste(file.dir, "mRNA_level.txt", sep="")
#express.level.per.transcript <- read.table(express.level.per.transcript.file, head=TRUE, fill=TRUE, sep="\t")
#express.level.per.transcript$mediane.transcript.expr <- apply(express.level.per.transcript[,-1], 1, median)
#hist(log(circad.proteins$mediane.transcript.expr), breaks = 150)

# cost per protein :
file.name <- paste(file.dir, "mouse_AA_synthesis_average_cost_per_protein.txt", sep="")
aa.synthesis.average.cost.per.protein <- read.table(file.name, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Protein abundance data available :
protein.abundance.file.dir <- paste(file.dir, "liver.txt", sep="")
protein.abundance.dataset <- read.table(protein.abundance.file.dir, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
colnames(protein.abundance.dataset)[2] <- "protein.abundance"

protein.abundance.dataset <- read.table("~/Documents/DATA/Mus_Musculus/Liver_Protein_Mauvoisin/Combined_WT/proteinGroups.txt", sep="\t", fill= TRUE, h=T,  check.names = FALSE, stringsAsFactors = FALSE)
normalized.count.values <- grep("Ratio H/L coun", colnames(protein.abundance.dataset))
protein.abundance.dataset <- protein.abundance.dataset[, c(2, 6, 7, normalized.count.values)]
protein.abundance.dataset <- protein.abundance.dataset[, -4]
colnames(protein.abundance.dataset)[1:3] <- c("Uniprot.ID", "Protein.Name", "Gene.Name")
colnames(protein.abundance.dataset) <- gsub("Ratio H/L count ", "", colnames(protein.abundance.dataset))

reviewed.status <- read.table("~/Downloads/reviewed_mouse.tab", h=T, fill = T, check.names = F, stringsAsFactors = F)
length(unique(reviewed.status$Entry))

is.integer0 <- function(x) {
  is.integer(x) && length(x) == 0L 
  }
# If there are several Protein.IDs (Uniprot.IDs) for the same data, we keep only those "rewieved" by Uniprot
for (i in 1:nrow(protein.abundance.dataset)) {
  if (grepl(";", protein.abundance.dataset[i, "Uniprot.ID"])) {
    
    protein.group <- strsplit(as.character(protein.abundance.dataset[i, "Uniprot.ID"]), split = ";")
    # Search if one or more of them are within Uniprot "REVIEWED" proteins
    match.proteins.ids <- grep(paste(unlist(protein.group),collapse="|"), reviewed.status[, "Uniprot.ID"], value=FALSE)
    
    if (is.integer0(match.proteins.ids)) {
      
    } else if (length(match.proteins.ids) == 1) {
      
    } else if (length(match.proteins.ids) > 1) {
      
    } else {
      print("... Error ... ??")
    }
    
    
    length(grep(unlist(protein.group), reviewed.status[, "Entry"]))
    grep("O70423", reviewed.status[, "Entry"])
    reviewed.status[8801,]
    class(reviewed.status[, "Entry"])
    
    any(unlist(protein.group) %in% reviewed.status[, "Entry"])
    
    class(raw.data$Protein.IDs)
    raw.data.bis <- data.frame(V1 = rep(raw.data$Protein.names, sapply(raw.data.bis, length)), V2 = unlist(raw.data.bis))
    test3 <- merge(reviewed.status, raw.data.bis, by="Entry", by.y="V2")
    
    
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

reviewed.status <- read.table("~/Downloads/reviewed_mouse.txt", header = T, sep = "\t", fill = TRUE)
reviewed.status[grep("P49443", reviewed.status$Entry),]
length(unique(reviewed.status$Protein.names))
reviewed.status <- unique(reviewed.status[, c(4,7)])




protein.abundance.dataset_bis <- read.table("~/Documents/DATA/Mus_Musculus/Liver_Protein_Mauvoisin/Combined_WT/peptides.txt", sep="\t", fill= TRUE, h=T,  check.names = FALSE, stringsAsFactors = FALSE)


nb.timepoints <- ncol(protein.abundance.dataset)-1
protein.abundance.dataset$protein.abundance <- apply(protein.abundance.dataset[,-1], 1, FUN = function(x){mean(max(x), max(x[x!=max(x)]))})
#protein.abundance.dataset$protein.abundance <- apply(protein.abundance.dataset[,-1], 1, max)

protein.abundance.dataset <- protein.abundance.dataset[, c("Gene.Name", "protein.abundance")]
protein.abundance.dataset$protein.abundance <- protein.abundance.dataset$protein.abundance + abs(min(protein.abundance.dataset$protein.abundance))
protein.abundance.dataset <- merge(protein.abundance.dataset, geneID.proteinID, by="Gene.Name")
protein.abundance.dataset <- unique(protein.abundance.dataset[, c("Protein.ID", "protein.abundance")])

# Merging
protein.abundance.and.cost <- merge(protein.abundance.dataset, aa.synthesis.average.cost.per.protein, by="Protein.ID")
protein.abundance.and.cost <- unique(protein.abundance.and.cost)

# Total cost for steady state protein level 
protein.abundance.and.cost$total.protein.cost <- protein.abundance.and.cost$protein.abundance * protein.abundance.and.cost$protein.length * (protein.abundance.and.cost$aa.synthesis.average.cost - 1)


# Merging
################
#### IMPORTANT STEP ####
all.data <- merge(protein.abundance.and.cost, geneID.proteinID, by="Protein.ID")
all.data <- merge(circad.proteins, all.data, by="Gene.Name", all.x = TRUE)
# Remove proteins with no cost of expression :
all.data$total.protein.cost <- all.data$total.protein.cost +1
all.data <- subset(all.data, total.protein.cost > 1)
#all.data <- merge(all.data, express.level.per.transcript, by="Transcript.ID")

hist(log(all.data$total.protein.cost), breaks = 150)

all.data$default.pvalue <- p.adjust(all.data$default.pvalue, method = "fdr")
hist(all.data$default.pvalue, breaks = 150)
#############
### ANOVA ###
#############
# effect of total protein cost into the p-value representing the rhythmicity nature
anova.1factor <- aov(default.pvalue ~ log(total.protein.cost), data = all.data)
summary(anova.1factor)
# Calculate the Adjusted R-squared
anova.data <- all.data[, c("default.pvalue", "total.protein.cost")]
anova.data <- na.exclude(anova.data)
lm.1factor <- lm(default.pvalue ~ log(total.protein.cost), data = anova.data)
summary(lm.1factor)

# effect of protein level into the p-value representing the rhythmicity nature
anova.1factor <- aov(default.pvalue ~ log(protein.abundance), data = all.data)
summary(anova.1factor)
# Calculate the Adjusted R-squared
anova.data <- all.data[, c("default.pvalue", "protein.abundance")]
anova.data <- na.exclude(anova.data)
lm.1factor <- lm(default.pvalue ~ log(protein.abundance), data = anova.data)
summary(lm.1factor)

# effect of each factors needed to calculate the tot protein cost into the p-value representing the rhythmicity nature
anova.3factors <- aov(default.pvalue ~ log(aa.synthesis.average.cost)*log(protein.abundance)*log(protein.length), data = all.data)
summary(anova.3factors)

anova.data <- all.data[, c("default.pvalue", "aa.synthesis.average.cost", "protein.abundance", "protein.length")]
anova.data <- na.exclude(anova.data)
is.na(anova.data) <- sapply(anova.data, is.infinite)
lm.3factors <- lm(default.pvalue ~ log(aa.synthesis.average.cost)*log(protein.abundance)*log(protein.length), data = anova.data)
summary(lm.3factors)

# => non pas le coût par unité de protéine mais le coût lié à l'abondance nécessaire de ces protéines en steady-state

# Verification des hypotheses
# Contradictions avec les hypothèses ? Points aberrants ?
par(mfrow=c(2,2))
plot(anova.2factors)
#shapiro.test(resid(anova.2factors))
bartlett.test(log(all.data$total.protein.cost), all.data$circadian.gene)
#bartlett.test(log(all.data$total.protein.cost), all.data$transcript.expression.level)

# La variable circadian.gene influence significativement la variable total.protein.cost en raison de l’effet croisé significatif. 
# On vérifie ceci par un test de Fisher de 
# (H0) log(total.protein.cost)~transcript.expression.level , contre
# (H1) log(total.protein.cost)~transcript.expression.level*circadian.gene
anova(aov(log(total.protein.cost) ~ mediane.transcript.expr, data = all.data), anova.2factors)
# Ce test confirme qu’on peut rejeter (H0) 
# et ainsi conclure que le facteur circadian.gene est pertinent dans le modèle

ks.test(all.data[all.data$circadian.gene==TRUE, "total.protein.cost"], all.data[all.data$circadian.gene==FALSE, "total.protein.cost"])
ggplot(all.data, aes(log(total.protein.cost), fill = all.data[[ "circadian.gene" ]])) + 
  geom_density(alpha = 0.2) +
  geom_vline(xintercept = log(median(all.data[all.data$circadian.gene==TRUE, "weighted.total.protein.cost"])), linetype="dotted", 
             color = "blue", size=1.5) +
  geom_vline(xintercept = log(median(all.data[all.data$circadian.gene==F, "weighted.total.protein.cost"])), linetype="dotted", 
             color = "red", size=1.5)

subset.first.rhythmic.proteins <- all.data[order(all.data$default.pvalue),]
subset.first.rhythmic.proteins <- subset.first.rhythmic.proteins[1:500, ]
plot(subset.first.rhythmic.proteins$default.pvalue, log(subset.first.rhythmic.proteins$weighted.total.protein.cost))
test.lm <- lm(default.pvalue ~ weighted.total.protein.cost, data = subset.first.rhythmic.proteins)
test.lm$df.residual
ggplot(all.data, aes(default.pvalue, log(weighted.total.protein.cost))) +
  geom_point() +
  stat_smooth(method = lm)

ks.test(all.data[all.data$circadian.gene==TRUE, "weighted.total.protein.cost"], all.data[all.data$circadian.gene==FALSE, "weighted.total.protein.cost"])
ggplot(all.data, aes(log(weighted.total.protein.cost), fill = all.data[[ "circadian.gene" ]])) + 
  geom_density(alpha = 0.2)

# MEME CHOSE  mais en fonction du cout moyen par AA
anova.1factor <- aov(log(aa.synthesis.average.cost) ~ circadian.gene, data = all.data)
summary(anova.1factor)
anova.2factors <- aov(log(aa.synthesis.average.cost) ~ mediane.transcript.expr*circadian.gene, data = all.data)
summary(anova.2factors)

plot <- ggplot(all.data, aes(log(aa.synthesis.average.cost), fill = all.data[[ "circadian.gene" ]])) + 
  geom_density(alpha = 0.2)









all.data.timepoints <- cbind(raw.dataset, circad.proteins[, c(-1,-2)])

all.data <- merge(protein.abundance.and.cost, geneID.proteinID, by="Protein.ID")
all.data.timepoints <- merge(all.data.timepoints, all.data, by="Gene.Name", all.x = TRUE)


time.points <- colnames(raw.dataset)[c(-1, -2)]

for (i in 1:length(time.points)) {
  all.data.timepoints[[time.points[i]]] <- all.data.timepoints[[time.points[i]]] * all.data.timepoints$protein.length * all.data.timepoints$aa.synthesis.average.cost
}

anova.1factor <- aov(default.pvalue ~ log(ZT00), data = subset(all.data.timepoints,  (!is.na(all.data.timepoints[,"ZT00"])) & (!is.na(all.data.timepoints[,"default.pvalue"]))))
summary(anova.1factor)







#### 
#### 
#### 
# Gene.ID and Protein.ID : 
geneID.EntrezID <- read.table("~/Documents/DATA/Genes_ID/GeneID_EntrezGeneID_Mouse.txt", head=TRUE, fill = TRUE)
####
all.data <- subset(all.data, protein.abundance > 5)

all.data <- merge(all.data, geneID.EntrezID, by="Gene.ID")

circad.proteins <- subset(all.data, default.pvalue <= 0.01)
non.circad.proteins <- subset(all.data, default.pvalue > 0.8)
non.circad.proteins <- non.circad.proteins[sample(1:nrow(non.circad.proteins), size = nrow(circad.proteins)), ]

hist(log(circad.proteins$total.protein.cost), breaks = 150)
hist(log(non.circad.proteins$total.protein.cost), breaks = 150)
t.test(log(circad.proteins$total.protein.cost), log(non.circad.proteins$total.protein.cost), alternative = "greater")
t.test(circad.proteins$total.protein.cost, non.circad.proteins$total.protein.cost, alternative = "greater")

t.test(circad.proteins$total.protein.cost*24, non.circad.proteins$total.protein.cost*24, alternative = "greater")

t.test(log(circad.proteins$protein.abundance), log(non.circad.proteins$protein.abundance), alternative = "greater")
ks.test(log(circad.proteins$protein.abundance), log(non.circad.proteins$protein.abundance))
ks.test(circad.proteins$protein.abundance, non.circad.proteins$protein.abundance)

ks.test(log(circad.proteins$total.protein.cost), log(non.circad.proteins$total.protein.cost), alternative = "greater")
ks.test(circad.proteins$total.protein.cost, non.circad.proteins$total.protein.cost, alternative = "greater")
ks.test(circad.proteins$total.protein.cost, non.circad.proteins$total.protein.cost)

ks.test(log(circad.proteins$total.protein.cost*24), log(non.circad.proteins$total.protein.cost*24))




write.table(as.data.frame(unique(circad.proteins$Gene.ID)), "~/Documents/circad_proteins_genes.txt", row.names = F, quote = F, col.names = F)
write.table(as.data.frame(unique(non.circad.proteins$Gene.ID)), "~/Documents/non_circad_proteins_genes.txt", row.names = F, quote = F, col.names = F)

write.table(as.data.frame(unique(circad.proteins$Transcript.ID)), "~/Documents/circad_proteins_genes.txt", row.names = F, quote = F, col.names = F)
write.table(as.data.frame(unique(non.circad.proteins$Transcript.ID)), "~/Documents/non_circad_proteins_genes.txt", row.names = F, quote = F, col.names = F)


#source("https://bioconductor.org/biocLite.R")
#biocLite("goProfiles")
library(goProfiles)

require(goProfiles)
data(prostateIds)
welsh.MF <- basicProfile(welsh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db")
singh.MF <- basicProfile(singh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db")
welsh.singh.MF <- mergeProfilesLists(welsh.MF, singh.MF, profNames=c("nonRhythmic_nb", "Rhythmic_nb"))
printProfiles(welsh.singh.MF, percentage=TRUE)

compared.welsh.singh.01.MF <- compareGeneLists(welsh01EntrezIDs[1:100], singh01EntrezIDs[1:100], onto="MF", level=2, orgPackage="org.Hs.eg.db")
plotProfiles(welsh.singh.MF, percentage=T, aTitle="Welsh vs Singh", legend=T)

class(welsh.singh.MF$MF[[3]])

panther.data <- read.table("~/Downloads/panther.txt", sep="\t", head=TRUE, fill = TRUE)
colnames(panther.data) <- gsub("non_circad_proteins_genes.txt", "nonRhythmic", colnames(panther.data))
colnames(panther.data) <- gsub("circad_proteins_genes.txt", "Rhythmic", colnames(panther.data))
panther.data <- panther.data[, c("GO.biological.process.complete", "nonRhythmic_220", "nonRhythmic_expected", "nonRhythmic_fold.Enrichment", 
                                 "Rhythmic_86", "Rhythmic_expected", "Rhythmic_fold.Enrichment")]
colnames(panther.data) <- c("GO.biological.process.complete", "nonRhythmic_nb", "nonRhythmic_expected", "nonRhythmic_fold.Enrichment", 
                            "Rhythmic_nb", "Rhythmic_expected", "Rhythmic_fold.Enrichment")
library(gtools)
Description <- unlist(strsplit.description)[odd(1:length(strsplit.description))]
Description <- trimws(Description, "l")
Description <- trimws(Description, "r")
GOID <- unlist(strsplit.description)[!odd(1:length(strsplit.description))]
panther.data$Description <- Description
panther.data$GOID <- GOID

welsh.MF$MF <- panther.data[, c("Description", "GOID", "nonRhythmic_nb")]
colnames(welsh.MF$MF)[3] <- "Frequency"
singh.MF$MF <- panther.data[, c("Description", "GOID", "Rhythmic_nb")]
colnames(singh.MF$MF)[3] <- "Frequency"
welsh.singh.MF <- mergeProfilesLists(welsh.MF, singh.MF, profNames=c("nonRhythmic_nb", "Rhythmic_nb"))
plotProfiles(welsh.singh.MF, percentage=T, aTitle="Welsh vs Singh", legend=F)

class(welsh.singh.MF$MF)
class(welsh.MF$MF) <- c("BasicGOProfile", "data.frame")
class(singh.MF$MF$Description) <- c("BasicGOProfile", "data.frame")

welsh.MF$MF$Description <- as.factor(panther.data$Description[1:19])
welsh.MF$MF$GOID <- as.factor(panther.data$GOID[1:19])
welsh.MF$MF$Frequency <- panther.data$nonRhythmic_nb[1:19]
singh.MF$MF$Description <- as.factor(panther.data$Description[1:19])
singh.MF$MF$GOID <- as.factor(panther.data$GOID[1:19])
singh.MF$MF$Frequency <- panther.data$Rhythmic_nb[1:19]







egIDs <- c("839235", "838362", "838961", "837091", "837455", "837543")

## use select to quickly translate these into TAIR IDs, and then grab that column of IDs back out.
## (You may find it more convenient to just start with the TAIR IDs that you said were in your file, but I don't have those here)
tairIDs <- as.character(select(test.1, keys=egIDs, cols="TAIR", keytype="ENTREZID")[[2]])

## THEN call basicProfile function and pass in tair IDs instead...
## Now when it calls mget on the GO mapping, it will actually get some matches.
basicProfile(tairIDs, idType="Entrez", onto ="ANY", level=2, orgPackage="org.At.tair.db", ordúLSE)



