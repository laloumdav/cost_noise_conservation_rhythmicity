# *************
############## ARABIDOPSIS Leaves ###################
# *************
dnds <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/dnds_with_lyrata.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- unique(dnds[!is.na(dnds$DNDS), ])
dnds[is.infinite(dnds$DNDS), ] = 1
#length(unique(dnds$Gene))
hist(log2(dnds$DNDS), breaks = 700)
hist(dnds$DNDS, breaks = 700)

#*** PROTEIN ****#
prot.data <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/proteome/proteome_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
prot.data$Protein.ID <- gsub("\\..*", "", prot.data$Protein.ID)
prot.data <- merge(prot.data, dnds, by.x="Protein.ID", by.y="Gene.ID")
head(sort(prot.data$DNDS, decreasing = T))
hist(prot.data$DNDS, breaks = 150)
hist(-log(prot.data$DNDS), breaks = 150)

rProt <- subset(prot.data, protein_rhythm.pvalue <= .01)
nrProt <- subset(prot.data, protein_rhythm.pvalue > .8)

# Nb of genes for which we have protein data: 
nrow(prot.data)
# Nb of rhythmic prot : 
nrow(rProt)
# Nb of non-rhythmic prot : 
nrow(nrProt)

t.test(rProt$DNDS, nrProt$DNDS)

# residuals
plot(log(prot.data$protein_mean.level), log(prot.data$DNDS))
cor.test(log(prot.data$protein_mean.level), log(prot.data$DNDS))

lmTest <- lm(log(protein_mean.level) ~ log(DNDS), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .01]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.8]
t.test(rProt, nrProt)

lmTest <- lm(log(DNDS) ~ log(protein_mean.level), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .01]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.8]
t.test(rProt, nrProt)


#*** RNA ****#
geneIDs <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/GeneIDs.txt", head=TRUE, fill = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
geneIDs <- unique(subset(geneIDs, Protein.ID != ""))
geneIDs <- unique(geneIDs[, c("AFFY.ATH1", "Protein.ID")])
# RNA
RNA.data <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/transcriptome/transcriptome_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.data <- merge(RNA.data, geneIDs, by="AFFY.ATH1")
RNA.data$Protein.ID <- gsub("\\..*", "", RNA.data$Protein.ID)
RNA.level <- unique(RNA.data[, c("Protein.ID", grep("RNA_LD", colnames(RNA.data), value = TRUE))])

# Aggregate multiple data for one gene 
RNA.rhythm <- unique(RNA.data[, c("Protein.ID", "RNA_rhythm.pvalue")])

# 2) Brown normalization:
genes.counts <- plyr::count(RNA.rhythm$Protein.ID)
genes.in.several.copies <- subset(genes.counts, freq >1)[,"x"]
library(EmpiricalBrownsMethod)
RNA.rhythm.unique <- data.frame(Protein.ID = unique(RNA.rhythm$Protein.ID),
                                RNA_rhythm.pvalue = NA)
for (k in 1:nrow(RNA.rhythm)) {
  gene.id <- RNA.rhythm$Protein.ID[k]
  
  if (gene.id %in% genes.in.several.copies) {
    original.data <- subset(RNA.level, Protein.ID == gene.id)[, -1]
    # modify the original.data[[k]] adding 1e-16 at the first column if all columns have 0 values
    # because empiricalBrownsMethod() does not work otherwise
    # => this step is not needed if these data have been removed initialy
    for (l in 1:nrow(original.data)){
      if (sum(original.data[l, ])==0){
        original.data[l, 1] <- 1e-16
      }
    }
    
    pvalues <- subset(RNA.rhythm, Protein.ID == gene.id)$RNA_rhythm.pvalue
    empirical.browns <- empiricalBrownsMethod(data_matrix=original.data, p_values=pvalues, extra_info=TRUE)
    
    RNA.rhythm.unique[grep(gene.id, RNA.rhythm.unique$Protein.ID), "RNA_rhythm.pvalue"] <- empirical.browns[[1]] 
  } else {
    RNA.rhythm.unique[grep(gene.id, RNA.rhythm.unique$Protein.ID), "RNA_rhythm.pvalue"] <- RNA.rhythm$RNA_rhythm.pvalue[k]
  }
}

# Aggregate multiple data for one gene 
RNA.level <- aggregate(RNA.level[, -1], by=list(RNA.level$Protein.ID), mean)
colnames(RNA.level)[1] <- "ID"
RNA.level$RNA_mean.level <- apply(RNA.level[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.level$RNA_max.level <- apply(RNA.level[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.data <- merge(RNA.level, RNA.rhythm.unique, by.x="ID", by.y="Protein.ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")
# remove extreme value 
RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .8)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > 0.8]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > 0.8]
t.test(rRNA, nrRNA)

# Both rhythmic
tot.data <- unique(merge(prot.data[, c("Protein.ID", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level")], 
                  RNA.data[, c("ID", "RNA_rhythm.pvalue", "RNA_mean.level", "RNA_max.level", "DNDS")], by.x="Protein.ID", by.y="ID"))

rProt <- subset(tot.data, protein_rhythm.pvalue <= .05)
nrow(rProt)
rProt.rRNA <- subset(rProt, RNA_rhythm.pvalue <= .01)
nrow(rProt.rRNA)
rProt.nrRNA <- subset(rProt, RNA_rhythm.pvalue > .3)
nrow(rProt.nrRNA)

t.test(rProt.rRNA$DNDS, rProt.nrRNA$DNDS)


rRNA <- subset(tot.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
rRNA.rProt <- subset(rRNA, protein_rhythm.pvalue <= .05)
nrow(rRNA.rProt)
rRNA.nrProt <- subset(rRNA, protein_rhythm.pvalue > .3)
nrow(rRNA.nrProt)

t.test(rRNA.rProt$DNDS, rRNA.nrProt$DNDS)



# *************
############## MOUSE Liver ###################
# *************
dnds <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.Name", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1

#*** PROTEIN ****#
prot.data <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/proteome/proteome_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
prot.level <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/proteome/proteome_level.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
prot.level$protein_mean.level <- apply(prot.level[, c(-1, -2, -3, -4)], 1, FUN = function(x){mean(x, na.rm = TRUE)})
prot.level$protein_max.level <- apply(prot.level[, c(-1, -2, -3, -4)], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
prot.data <- cbind(prot.data[, c("Gene.Name", "protein_rhythm.pvalue")], prot.level[, c("protein_mean.level", "protein_max.level")])
#gene.ids <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/Ensembl_geneName_Affy.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
#gene.ids <- unique(gene.ids[, c("Gene.ID", "Gene.Name")])
#prot.data <- merge(gene.ids, prot.data, by="Gene.Name")
prot.data <- merge(prot.data, dnds, by="Gene.Name")
head(sort(prot.data$DNDS, decreasing = T))
hist(prot.data$DNDS, breaks = 150)
hist(-log(prot.data$DNDS), breaks = 150)

rProt <- subset(prot.data, protein_rhythm.pvalue <= .01)
nrow(rProt)
nrProt <- subset(prot.data, protein_rhythm.pvalue > .7)
nrow(nrProt)
t.test(rProt$DNDS, nrProt$DNDS)

# residuals
plot(log(prot.data$protein_mean.level), log(prot.data$DNDS))
cor.test(log(prot.data$protein_mean.level), log(prot.data$DNDS))

lmTest <- lm(log(protein_mean.level) ~ log(DNDS), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .01]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.7]
t.test(rProt, nrProt, na.action=na.omit)

lmTest <- lm(log(DNDS) ~ log(protein_mean.level), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .01]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.7]
t.test(rProt, nrProt, na.action=na.omit)

#*** RNA ****#
dnds <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1

RNA.data <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/transcriptome/transcriptome_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
gene.ids <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/Ensembl_geneName_Affy.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
gene.ids <- unique(gene.ids[, c("Gene.ID", "AFFY.ATH1")])
RNA.data <- unique(merge(RNA.data, gene.ids, by.x="mRNA.ID", by.y="AFFY.ATH1"))
RNA.data <- merge(RNA.data, dnds, by="Gene.ID")
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .001)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .2)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .001]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .2]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .001]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .2]
t.test(rRNA, nrRNA)

# Both rhythmic
tot.data <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/tot_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
gene.ids <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/Ensembl_geneName_Affy.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
gene.ids <- unique(gene.ids[, c("Gene.ID", "AFFY.ATH1")])
tot.data <- unique(merge(tot.data, gene.ids, by.x="mRNA.ID", by.y="AFFY.ATH1"))
tot.data <- merge(tot.data, dnds, by="Gene.ID")
head(sort(tot.data$DNDS, decreasing = T))

rProt <- subset(tot.data, protein_rhythm.pvalue <= .05)
nrow(rProt)
rProt.rRNA <- subset(rProt, RNA_rhythm.pvalue <= .01)
nrow(rProt.rRNA)
rProt.nrRNA <- subset(rProt, RNA_rhythm.pvalue > .3)
nrow(rProt.nrRNA)

t.test(rProt.rRNA$DNDS, rProt.nrRNA$DNDS)


rRNA <- subset(tot.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
rRNA.rProt <- subset(rRNA, protein_rhythm.pvalue <= .05)
nrow(rRNA.rProt)
rRNA.nrProt <- subset(rRNA, protein_rhythm.pvalue > .3)
nrow(rRNA.nrProt)

t.test(rRNA.rProt$DNDS, rRNA.nrProt$DNDS)

# *************
############## MOUSE Tendon ###################
# *************
dnds <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.Name", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1

#*** PROTEIN ****#
prot.data <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/tendon/proteome/proteome_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
prot.data <- prot.data[, c("Gene.Name", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
#gene.ids <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/Ensembl_geneName_Affy.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
#gene.ids <- unique(gene.ids[, c("Gene.ID", "Gene.Name")])
#prot.data <- merge(gene.ids, prot.data, by="Gene.Name")
prot.data <- unique(merge(prot.data, dnds, by="Gene.Name"))
head(sort(prot.data$DNDS, decreasing = T))
hist(prot.data$DNDS, breaks = 150)
hist(-log(prot.data$DNDS), breaks = 150)

rProt <- subset(prot.data, protein_rhythm.pvalue <= .05)
nrow(rProt)
nrProt <- subset(prot.data, protein_rhythm.pvalue > .4)
nrow(nrProt)
t.test(rProt$DNDS, nrProt$DNDS)

# residuals
plot(log(prot.data$protein_mean.level), log(prot.data$DNDS))
cor.test(log(prot.data$protein_mean.level), log(prot.data$DNDS))

lmTest <- lm(log(protein_mean.level) ~ log(DNDS), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .05]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.4]
t.test(rProt, nrProt, na.action=na.omit)

lmTest <- lm(log(DNDS) ~ log(protein_mean.level), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .05]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.4]
t.test(rProt, nrProt, na.action=na.omit)


# *************
############## MOUSE Tendon ###################
# *************
dnds <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.Name", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1

#*** PROTEIN ****#
prot.data <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/forebrain/proteome/proteome_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
prot.data <- prot.data[, c("Gene.Name", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
#gene.ids <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/Ensembl_geneName_Affy.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
#gene.ids <- unique(gene.ids[, c("Gene.ID", "Gene.Name")])
#prot.data <- merge(gene.ids, prot.data, by="Gene.Name")
prot.data <- unique(merge(prot.data, dnds, by="Gene.Name"))
head(sort(prot.data$DNDS, decreasing = T))
hist(prot.data$DNDS, breaks = 150)
hist(-log(prot.data$DNDS), breaks = 150)

rProt <- subset(prot.data, protein_rhythm.pvalue <= .05)
nrow(rProt)
nrProt <- subset(prot.data, protein_rhythm.pvalue > .4)
nrow(nrProt)
t.test(rProt$DNDS, nrProt$DNDS)

# residuals
plot(log(prot.data$protein_mean.level), log(prot.data$DNDS))
cor.test(log(prot.data$protein_mean.level), log(prot.data$DNDS))

lmTest <- lm(log(protein_mean.level) ~ log(DNDS), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .05]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.4]
t.test(rProt, nrProt, na.action=na.omit)

lmTest <- lm(log(DNDS) ~ log(protein_mean.level), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .05]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.4]
t.test(rProt, nrProt, na.action=na.omit)


# *************
############## MOUSE Cartilage ###################
# *************
dnds <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.Name", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1

#*** PROTEIN ****#
prot.data <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/cartilage/proteome/proteome_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
prot.data <- prot.data[, c("Gene.Name", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
#gene.ids <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/Ensembl_geneName_Affy.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
#gene.ids <- unique(gene.ids[, c("Gene.ID", "Gene.Name")])
#prot.data <- merge(gene.ids, prot.data, by="Gene.Name")
prot.data <- unique(merge(prot.data, dnds, by="Gene.Name"))
head(sort(prot.data$DNDS, decreasing = T))
hist(prot.data$DNDS, breaks = 150)
hist(-log(prot.data$DNDS), breaks = 150)

rProt <- subset(prot.data, protein_rhythm.pvalue <= .05)
nrow(rProt)
nrProt <- subset(prot.data, protein_rhythm.pvalue > .4)
nrow(nrProt)
t.test(rProt$DNDS, nrProt$DNDS)

# residuals
plot(log(prot.data$protein_mean.level), log(prot.data$DNDS))
cor.test(log(prot.data$protein_mean.level), log(prot.data$DNDS))

lmTest <- lm(log(protein_mean.level) ~ log(DNDS), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .05]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.4]
t.test(rProt, nrProt, na.action=na.omit)

lmTest <- lm(log(DNDS) ~ log(protein_mean.level), prot.data)
rProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue <= .05]
nrProt <- lmTest$residuals[prot.data$protein_rhythm.pvalue > 0.4]
t.test(rProt, nrProt, na.action=na.omit)



# *************
############## MOUSE Lung ###################
# *************
dnds <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1
hist(log2(dnds$DNDS), breaks = 50)
hist(dnds$DNDS, breaks = 50)

RNA.data <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/lung/lung.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
# Aggregate multiple data for one gene 
RNA.data <- aggregate(RNA.data[, -1], by=list(RNA.data$ID), mean)
colnames(RNA.data)[1] <- "ID"
RNA.data$RNA_mean.level <- apply(RNA.data[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.data$RNA_max.level <- apply(RNA.data[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.rhythm <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/lung/normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.rhythm <- data.frame(ID = subset(RNA.rhythm, algorithm == "GeneCycle")$ID,
                         RNA_rhythm.pvalue = subset(RNA.rhythm, algorithm == "GeneCycle")$pvalue_brown)
RNA.data <- merge(RNA.data, RNA.rhythm, by="ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")

hist(RNA.data$RNA_rhythm.pvalue, breaks = 80)
head(sort(RNA.data$DNDS, decreasing = T))
# remove extreme value 
RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .001)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .8)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .001]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .001]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)


# *************
############## MOUSE Muscle ###################
# *************
dnds <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1
hist(log2(dnds$DNDS), breaks = 50)
hist(dnds$DNDS, breaks = 50)

RNA.data <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/muscle/muscle.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
# Aggregate multiple data for one gene 
RNA.data <- aggregate(RNA.data[, -1], by=list(RNA.data$ID), mean)
colnames(RNA.data)[1] <- "ID"
RNA.data$RNA_mean.level <- apply(RNA.data[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.data$RNA_max.level <- apply(RNA.data[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.rhythm <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/muscle/normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.rhythm <- data.frame(ID = subset(RNA.rhythm, algorithm == "GeneCycle")$ID,
                         RNA_rhythm.pvalue = subset(RNA.rhythm, algorithm == "GeneCycle")$pvalue_brown)
RNA.data <- merge(RNA.data, RNA.rhythm, by="ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")

hist(RNA.data$RNA_rhythm.pvalue, breaks = 80)
head(sort(RNA.data$DNDS, decreasing = T))
# remove extreme value 
RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .4)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)


# *************
############## MOUSE Heart ###################
# *************
dnds <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1
hist(log2(dnds$DNDS), breaks = 50)
hist(dnds$DNDS, breaks = 50)

RNA.data <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/heart/heart.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
# Aggregate multiple data for one gene 
RNA.data <- aggregate(RNA.data[, -1], by=list(RNA.data$ID), mean)
colnames(RNA.data)[1] <- "ID"
RNA.data$RNA_mean.level <- apply(RNA.data[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.data$RNA_max.level <- apply(RNA.data[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.rhythm <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/heart/normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.rhythm <- data.frame(ID = subset(RNA.rhythm, algorithm == "GeneCycle")$ID,
                         RNA_rhythm.pvalue = subset(RNA.rhythm, algorithm == "GeneCycle")$pvalue_brown)
RNA.data <- merge(RNA.data, RNA.rhythm, by="ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")

hist(RNA.data$RNA_rhythm.pvalue, breaks = 80)
head(sort(RNA.data$DNDS, decreasing = T))
# remove extreme value 
RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .4)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)


# *************
############## MOUSE kidney ###################
# *************
dnds <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1
hist(log2(dnds$DNDS), breaks = 50)
hist(dnds$DNDS, breaks = 50)

RNA.data <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/kidney/kidney.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
# Aggregate multiple data for one gene 
RNA.data <- aggregate(RNA.data[, -1], by=list(RNA.data$ID), mean)
colnames(RNA.data)[1] <- "ID"
RNA.data$RNA_mean.level <- apply(RNA.data[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.data$RNA_max.level <- apply(RNA.data[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.rhythm <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/kidney/normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.rhythm <- data.frame(ID = subset(RNA.rhythm, algorithm == "GeneCycle")$ID,
                         RNA_rhythm.pvalue = subset(RNA.rhythm, algorithm == "GeneCycle")$pvalue_brown)
RNA.data <- merge(RNA.data, RNA.rhythm, by="ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")

hist(RNA.data$RNA_rhythm.pvalue, breaks = 80)
head(sort(RNA.data$DNDS, decreasing = T))
# remove extreme value 
RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .4)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)



# *************
############## MOUSE Aorta ###################
# *************
dnds <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/dnds_with_rat.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1
hist(log2(dnds$DNDS), breaks = 50)
hist(dnds$DNDS, breaks = 50)

RNA.data <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/aorta/aorta.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
# Aggregate multiple data for one gene 
RNA.data <- aggregate(RNA.data[, -1], by=list(RNA.data$ID), mean)
colnames(RNA.data)[1] <- "ID"
RNA.data$RNA_mean.level <- apply(RNA.data[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.data$RNA_max.level <- apply(RNA.data[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.rhythm <- read.table("~/Documents/rhythm_detection_benchmark/DATA/mouse_microarray/aorta/normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.rhythm <- data.frame(ID = subset(RNA.rhythm, algorithm == "GeneCycle")$ID,
                         RNA_rhythm.pvalue = subset(RNA.rhythm, algorithm == "GeneCycle")$pvalue_brown)
RNA.data <- merge(RNA.data, RNA.rhythm, by="ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")

hist(RNA.data$RNA_rhythm.pvalue, breaks = 80)
head(sort(RNA.data$DNDS, decreasing = T))
# remove extreme value 
RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .4)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)




# *************
############## RAT Lung ###################
# *************
dnds <- read.table("~/Documents/rhythm_detection_benchmark/DATA/rat/dnds_with_mouse.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1
hist(log2(dnds$DNDS), breaks = 50)
hist(dnds$DNDS, breaks = 50)

RNA.data <- read.table("~/Documents/rhythm_detection_benchmark/DATA/rat/lung/lung.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
# Aggregate multiple data for one gene 
RNA.data <- aggregate(RNA.data[, -1], by=list(RNA.data$ID), mean)
colnames(RNA.data)[1] <- "ID"
RNA.data$RNA_mean.level <- apply(RNA.data[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.data$RNA_max.level <- apply(RNA.data[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.rhythm <- read.table("~/Documents/rhythm_detection_benchmark/DATA/rat/lung/normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.rhythm <- data.frame(ID = subset(RNA.rhythm, algorithm == "GeneCycle")$ID,
                         RNA_rhythm.pvalue = subset(RNA.rhythm, algorithm == "GeneCycle")$pvalue_brown)
RNA.data <- merge(RNA.data, RNA.rhythm, by="ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")

hist(RNA.data$RNA_rhythm.pvalue, breaks = 80)
head(sort(RNA.data$DNDS, decreasing = T))
# remove extreme value 
#RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
#head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .001)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .4)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .001]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .001]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)


# *************
############## Anopheles head ###################
# *************
dnds <- read.table("~/Documents/rhythm_detection_benchmark/DATA/anopheles/dnds_with_aedes.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1
hist(log2(dnds$DNDS), breaks = 50)
hist(dnds$DNDS, breaks = 50)

RNA.data <- read.table("~/Documents/rhythm_detection_benchmark/DATA/anopheles/head/head.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
# Aggregate multiple data for one gene 
RNA.data <- aggregate(RNA.data[, -1], by=list(RNA.data$ID), mean)
colnames(RNA.data)[1] <- "ID"
RNA.data$RNA_mean.level <- apply(RNA.data[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.data$RNA_max.level <- apply(RNA.data[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.rhythm <- read.table("~/Documents/rhythm_detection_benchmark/DATA/anopheles/head/normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.rhythm <- data.frame(ID = subset(RNA.rhythm, algorithm == "GeneCycle")$ID,
                         RNA_rhythm.pvalue = subset(RNA.rhythm, algorithm == "GeneCycle")$pvalue_brown)
RNA.data <- merge(RNA.data, RNA.rhythm, by="ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")

hist(RNA.data$RNA_rhythm.pvalue, breaks = 80)
head(sort(RNA.data$DNDS, decreasing = T))
# remove extreme value 
RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .4)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)


# *************
############## Anopheles Body ###################
# *************
dnds <- read.table("~/Documents/rhythm_detection_benchmark/DATA/anopheles/dnds_with_aedes.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1
hist(log2(dnds$DNDS), breaks = 50)
hist(dnds$DNDS, breaks = 50)

RNA.data <- read.table("~/Documents/rhythm_detection_benchmark/DATA/anopheles/body/body.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
# Aggregate multiple data for one gene 
RNA.data <- aggregate(RNA.data[, -1], by=list(RNA.data$ID), mean)
colnames(RNA.data)[1] <- "ID"
RNA.data$RNA_mean.level <- apply(RNA.data[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.data$RNA_max.level <- apply(RNA.data[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.rhythm <- read.table("~/Documents/rhythm_detection_benchmark/DATA/anopheles/body/normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.rhythm <- data.frame(ID = subset(RNA.rhythm, algorithm == "GeneCycle")$ID,
                         RNA_rhythm.pvalue = subset(RNA.rhythm, algorithm == "GeneCycle")$pvalue_brown)
RNA.data <- merge(RNA.data, RNA.rhythm, by="ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")

hist(RNA.data$RNA_rhythm.pvalue, breaks = 80)
head(sort(RNA.data$DNDS, decreasing = T))
# remove extreme value 
RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .4)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)


# *************
############## Aedes Head ###################
# *************
dnds <- read.table("~/Documents/rhythm_detection_benchmark/DATA/aedes/dnds_with_anopheles.txt", h=T, fill=T, sep="\t")
dnds <- unique(dnds[, c("Gene.ID", "DN", "DS")])
dnds$DNDS = dnds$DN / dnds$DS 
dnds <- dnds[!is.na(dnds$DNDS), ]
dnds[is.infinite(dnds$DNDS), ] = 1
hist(log2(dnds$DNDS), breaks = 50)
hist(dnds$DNDS, breaks = 50)

RNA.data <- read.table("~/Documents/rhythm_detection_benchmark/DATA/aedes/head/head.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
# Aggregate multiple data for one gene 
RNA.data <- aggregate(RNA.data[, -1], by=list(RNA.data$ID), mean)
colnames(RNA.data)[1] <- "ID"
RNA.data$RNA_mean.level <- apply(RNA.data[, -1], 1, FUN = function(x){mean(x, na.rm = TRUE)})
RNA.data$RNA_max.level <- apply(RNA.data[,-1], 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})

RNA.rhythm <- read.table("~/Documents/rhythm_detection_benchmark/DATA/aedes/head/normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
RNA.rhythm <- data.frame(ID = subset(RNA.rhythm, algorithm == "GeneCycle")$ID,
                         RNA_rhythm.pvalue = subset(RNA.rhythm, algorithm == "GeneCycle")$pvalue_brown)
RNA.data <- merge(RNA.data, RNA.rhythm, by="ID")
RNA.data <- merge(RNA.data, dnds, by.x="ID", by.y="Gene.ID")

hist(RNA.data$RNA_rhythm.pvalue, breaks = 80)
head(sort(RNA.data$DNDS, decreasing = T))
# remove extreme value 
RNA.data <- RNA.data[which(RNA.data$DNDS != sort(RNA.data$DNDS, decreasing = T)[1]), ]
head(sort(RNA.data$DNDS, decreasing = T))

rRNA <- subset(RNA.data, RNA_rhythm.pvalue <= .01)
nrow(rRNA)
nrRNA <- subset(RNA.data, RNA_rhythm.pvalue > .4)
nrow(nrRNA)
t.test(rRNA$DNDS, nrRNA$DNDS)

# residuals
plot(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))
cor.test(log(RNA.data$RNA_mean.level), log(RNA.data$DNDS))

lmTest <- lm(log(RNA_mean.level) ~ log(DNDS), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)

lmTest <- lm(log(DNDS)~log(RNA_mean.level), RNA.data)
rRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue <= .01]
nrRNA <- lmTest$residuals[RNA.data$RNA_rhythm.pvalue > .4]
t.test(rRNA, nrRNA)


























##########################
#######################
###################
################
#############
#########
####
#
rProt$group <- "rProt"
nrProt$group <- "nrProt"
kept.data <- rbind(rProt, nrProt)

ggplot(kept.data, aes(group, -log2(DNDS), fill = group)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = .5) +
  #facet_wrap(~ species, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=6)+
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=16), 
    axis.title.y = element_text(size=19), 
    axis.text.x = element_text(size=19), 
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    #panel.spacing = unit(.01, "lines"),
    #axis.line.x = element_line(colour = "grey"),
    #axis.line.y = element_line(colour = "grey"),
    legend.key = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "grey32"),
    strip.text = element_blank(),
    legend.position = "none")
####




rRNA$group <- "rRNA"
nrRNA$group <- "nrRNA"
kept.data <- rbind(rRNA, nrRNA)

ggplot(kept.data, aes(group, -log2(DNDS), fill = group)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = .5) +
  #facet_wrap(~ species, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=6)+
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=16), 
    axis.title.y = element_text(size=19), 
    axis.text.x = element_text(size=19), 
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    #panel.spacing = unit(.01, "lines"),
    #axis.line.x = element_line(colour = "grey"),
    #axis.line.y = element_line(colour = "grey"),
    legend.key = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "grey32"),
    strip.text = element_blank(),
    legend.position = "none")


####




kept.data <- subset(tot.data, (!is.na(tot.data[, "protein_rhythm.pvalue"]) & !is.na(tot.data[, "RNA_rhythm.pvalue"]) & !is.na(tot.data[, "DNDS"])))
kept.data <- kept.data[, c("ID", "protein_rhythm.pvalue", "RNA_rhythm.pvalue", "noise.adj.sd", "Fstar1", "Fstar4", "distance.to.median", "DNDS")]
# order genes by their rhythm p-value
kept.data <- kept.data[order(kept.data$protein_rhythm.pvalue), ]
kept.data$genes.order <- 1:nrow(kept.data)
kept.data$group <- "first 50%"
kept.data$group[round(nrow(kept.data)/2):nrow(kept.data)] <- "last 50%"

hist(kept.data$DNDS, breaks = 150)

test <- kept.data[kept.data$DNDS < quantile(kept.data$DNDS, 0.98) &
                    kept.data$DNDS > quantile(kept.data$DNDS, 0.02), ]
test <- test[(test$protein_rhythm.pvalue <= .05) | (test$RNA_rhythm.pvalue <= .02), ]
hist(-log(test$DNDS), breaks = 150)

plot1 <- ggplot(kept.data, aes(x=-log2(protein_rhythm.pvalue), y=-log2(RNA_rhythm.pvalue), color=DNDS)) +
  geom_point(size=1) +
  #geom_smooth(aes(group = group), size = .3, se = FALSE, color="black", span = 0.8) +
  scale_fill_manual(values=c("seagreen4", "#E7B800")) +
  #scale_color_gradientn(colours = rainbow(5)) +
  scale_color_gradient2(midpoint = .18, 
                        low = "darkorchid4", mid = "grey", high = "springgreen4", space = "Lab" ) +
  #scale_size_continuous(range=c(1,10)) +
  labs(y = expression(paste("rhythm ", italic("p"), "-value (RNA) (log)")),
       x = expression(paste("rhythm ", italic("p"), "-value (protein) (log)")),
       color="dn/ds") +
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=9), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_text(size=11),
    panel.background = element_blank(),
    legend.position = "top",
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    legend.key = element_blank(),
    legend.text = element_text(size=10), 
    axis.line = element_line(colour = "grey32"))
plot1 <- ggMarginal(plot1, size = 2, type = "histogram", col = "grey70", fill = "grey")
plot1



rRNA <- subset(kept.data, RNA_rhythm.pvalue <= .07)
rRNA$group <- rep("rRNA", nrow(rRNA))
nrRNA <- subset(kept.data, RNA_rhythm.pvalue > .8)
nrRNA$group <- rep("nrRNA", nrow(nrRNA))
t.test(rRNA$DNDS, nrRNA$DNDS)

rProt <- subset(prot.data, RNA_rhythm.pvalue <= .01)
rProt$group <- rep("rRNA", nrow(rProt))
nrProt <- subset(prot.data, RNA_rhythm.pvalue > .3)
nrProt$group <- rep("nrRNA", nrow(nrProt))
t.test(-log2(rProt$DNDS), -log2(nrProt$DNDS))



kept.data <- rbind(rProt, nrProt)

ggplot(kept.data, aes(group, -log2(DNDS), fill = group)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = 1) +
  #facet_wrap(~ species, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=2.4)+
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=6), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_text(size=11), 
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(.07, "lines"),
    #axis.line.x = element_line(colour = "grey"),
    #axis.line.y = element_line(colour = "grey"),
    legend.key = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "grey32"),
    strip.text = element_blank(),
    legend.position = "top",
    plot.caption = element_text(size=7, face = "italic", hjust=0))


