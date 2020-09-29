library(ggplot2)

paralogs <- read.table("Downloads/dgd_Mmu_all_v71.tsv", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)

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

prot.data.paralogs <- merge(prot.data, paralogs, by.x="Gene.Name", by.y="Name")
# remove paralogs that do not have data for at least 2 of them
nb.paralogs <- plyr::count(prot.data.paralogs$group_id)
prot.data.paralogs <- prot.data.paralogs[prot.data.paralogs$group_id %in% nb.paralogs[nb.paralogs$freq > 1, "x"], ]
# Aggregate by the mean in case of several data (dNdS)
prot.data.paralogs <- aggregate(prot.data.paralogs[, c( "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level", "DNDS", "group_id")], by = list(prot.data.paralogs$Gene.Name), FUN = mean)
# Again remove paralogs that do not have data for at least 2 of them
nb.paralogs <- plyr::count(prot.data.paralogs$group_id)
prot.data.paralogs <- prot.data.paralogs[prot.data.paralogs$group_id %in% nb.paralogs[nb.paralogs$freq > 1, "x"], ]
# order genes by their rhythm p-value
row.nb <- order(prot.data.paralogs$protein_rhythm.pvalue)
prot.data.paralogs[row.nb, "rhythmic.order"] <- 1:nrow(prot.data.paralogs)
plyr::count(prot.data.paralogs$group_id)

TOT.prot.data.paralogs <- NULL
for (i in 1:unique(prot.data.paralogs$group_id)) {
  paralog.group.id <- unique(prot.data.paralogs$group_id)[i]
  paralog.group <- subset(prot.data.paralogs, group_id==paralog.group.id)
  paralog.group.added <- data.frame(prot.data.paralogs[i, ], 
                     rhythmic.order.paralogs = prot.data.paralogs[prot.data.paralogs$group_id == paralog.group, "rhythmic.order"],
                     DNDS.paralogs = prot.data.paralogs[prot.data.paralogs$group_id == paralog.group, "DNDS"])[-1, ]
  
  TOT.prot.data.paralogs <- rbind(TOT.prot.data.paralogs, paralog.group.added)
}




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

RNA.data.paralogs <- merge(RNA.data, paralogs, by.x="Gene.ID", by.y="ENS_ID")
# remove paralogs that do not have data for at least 2 of them
nb.paralogs <- plyr::count(RNA.data.paralogs$group_id)
RNA.data.paralogs <- RNA.data.paralogs[RNA.data.paralogs$group_id %in% nb.paralogs[nb.paralogs$freq > 1, "x"], ]
# Aggregate by the mean in case of several data (dNdS)
RNA.data.paralogs <- aggregate(RNA.data.paralogs[, c( "RNA_rhythm.pvalue", "RNA_mean.level", "RNA_max.level", "DNDS", "group_id")], by = list(RNA.data.paralogs$Name), FUN = mean)
# Again remove paralogs that do not have data for at least 2 of them
nb.paralogs <- plyr::count(RNA.data.paralogs$group_id)
RNA.data.paralogs <- RNA.data.paralogs[RNA.data.paralogs$group_id %in% nb.paralogs[nb.paralogs$freq > 1, "x"], ]
# order genes by their rhythm p-value
row.nb <- order(RNA.data.paralogs$RNA_rhythm.pvalue)
RNA.data.paralogs[row.nb, "rhythmic.order"] <- 1:nrow(RNA.data.paralogs)
plyr::count(RNA.data.paralogs$group_id)

TOT.rhythmic.RNA.data.paralogs <- NULL
nb.rhythmic.paralog.gp = 0
nb.at.least.2.rhythmic.paralog = 0
for (i in 1:length(unique(RNA.data.paralogs$group_id))) {
  paralog.group.id <- unique(RNA.data.paralogs$group_id)[i]
  paralog.group <- subset(RNA.data.paralogs, group_id==paralog.group.id)
  
  if(any(paralog.group$RNA_rhythm.pvalue <= 0.05)) {
    if(length(which(paralog.group$RNA_rhythm.pvalue <= 0.05)) >1) {
      nb.at.least.2.rhythmic.paralog = nb.at.least.2.rhythmic.paralog +1
    }
    paralog.group <- data.frame(group_id=paralog.group.id,
                              rhythmic.order.1=paralog.group[1, "rhythmic.order"],
                              rhythmic.order.paralog=paralog.group[-1, "rhythmic.order"],
                              DNDS.1=paralog.group[1, "DNDS"],
                              DNDS.paralog=paralog.group[-1, "DNDS"])
  TOT.rhythmic.RNA.data.paralogs <- rbind(TOT.rhythmic.RNA.data.paralogs, paralog.group)
  nb.rhythmic.paralog.gp = nb.rhythmic.paralog.gp+1
  }
}
# Percentage of paralogs groups which have at least one of them rhythmic (p<0.05)
nb.rhythmic.paralog.gp/length(unique(RNA.data.paralogs$group_id))
# Percentage of paralogs groups which have at least 2 of them rhythmic (p<0.05) among all paralogs groups
nb.at.least.2.rhythmic.paralog/length(unique(RNA.data.paralogs$group_id))
# Percentage of paralogs groups which have at least 2 of them rhythmic (p<0.05) among paralogs genes with at least 1 of them rhythmic
nb.at.least.2.rhythmic.paralog/nb.rhythmic.paralog.gp

RNA.data.paralogs$group_id <- as.character(RNA.data.paralogs$group_id)
test <- lm(DNDS ~ group_id, RNA.data.paralogs)
summary(test)
test <- lm(RNA_rhythm.pvalue*DNDS ~ group_id, RNA.data.paralogs)
summary(test)


ggplot(TOT.rhythmic.RNA.data.paralogs, aes(x=rhythmic.order.1, y=rhythmic.order.paralog, group=rhythmic.order.1)) +
  geom_boxplot()


