# *** functions *** #

# Function to calculate Tau :
# for tissue i: xi=xi/max(xi)
# For a total of n tissues: tau=sum(1-xi)/(n-1)
ftau <- function(x) {
  if(!all(is.na(x)))
  {
    res <- sum(x, na.rm=TRUE)
    res <- res/(length(x)-1)
  } else {
    res <- NA
  }
  return(res)
}

trendFunction <- function(obs.value, theoretic.value, signif) {
  if ((obs.value < theoretic.value) & (signif <= 0.05)) {
    return("0>") 
  } else if ((obs.value > theoretic.value) & (signif <= 0.05)) { 
    return("0<") 
  } else if ((obs.value > theoretic.value) & (signif > 0.05)) { 
    return("0=<") 
  } else { 
    return("0>=") }
}

signiFunction <- function(pval) {
  p.formated <- NULL
  for (i in 1:length(pval)) {
    if (pval[i] < 2.2e-16) {
      p.formated[i] <- "< 2.2e-16" 
    } else if (pval[i] >= 0.01) {
      p.formated[i] <- round(pval[i], 4)
    } else { 
      p.formated[i] <- formatC(pval[i], format = "e", digits = 1)
    }
  }
  return(p.formated)
}
# *** functions *** #


library(reshape2)
library(ggplot2)
library(ggpubr)


table.tot.1 <- NULL
table.tot.2 <- NULL

main.dir <- "~/Documents/rhythm_detection_benchmark"
####################################
####################################
#######  MOUSE 
####################################
####################################
species <- "mouse_microarray"
tissue.list <- c("liver", "lung", "kidney", "muscle", "aorta", "adrenal_gland", "brain_stem", "brown_adipose", "cerebellum", "heart", "white_adipose")
p.val <- "default.pvalue"
######
#if (grepl("mouse", species)){ species.name <- "mouse" } else { species.name <- species }

full.dataframe <- data.frame()
for (t in 1:length(tissue.list)) {
  tissue <- tissue.list[t]
  
  file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")
  metacycle.file.dir <- paste(file.dir, "MetaCycle_on_aggregated_data", sep="")
  
  normalized.pval <- read.table(paste(file.dir, "normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if ("pvalue_brown" %in% colnames(normalized.pval)) {
    normalized.pval <- subset(normalized.pval, algorithm == "GeneCycle")[, c("ID", "pvalue_brown")]
    colnames(normalized.pval) <- c("ID", "GeneCycle.pvalue")
  } else {
    normalized.pval <- subset(normalized.pval, algorithm == "GeneCycle")[, -grep("algorithm", colnames(normalized.pval))]
    colnames(normalized.pval) <- c("ID", "GeneCycle.pvalue")
  }
  
  raw.file <- paste(file.dir, tissue, ".txt", sep = "")
  raw.data <- read.table(raw.file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  raw.data <- aggregate(raw.data[-1], by = raw.data[1], FUN = mean)
  raw.data$mean.RNA.level <- apply(raw.data[, grep("CT|ZT|LD|LL", colnames(raw.data))], 1, FUN = function(x){round(mean(x, na.rm = TRUE), digits = 1)})
  raw.data$max.RNA.level <- apply(raw.data[, grep("CT|ZT|LD|LL", colnames(raw.data))], 1, FUN = function(x){round(mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T)), digits = 1)})
  
  ### Z-score normalization of each dataset
  m.mean = mean(raw.data$mean.RNA.level, na.rm = TRUE)
  m.max = mean(raw.data$max.RNA.level, na.rm = TRUE)
  sd.mean = sd(raw.data$mean.RNA.level, na.rm = TRUE)
  sd.max = sd(raw.data$max.RNA.level, na.rm = TRUE)
  
  raw.data$zscore.mean.expr.value <- (raw.data$mean.RNA.level - m.mean) / sd.mean
  raw.data$zscore.max.expr.value <- (raw.data$max.RNA.level - m.max) / sd.max
  # to get positive values
  raw.data$zscore.mean.expr.value <- raw.data$zscore.mean.expr.value + abs(min(raw.data$zscore.mean.expr.value, na.rm = TRUE))
  raw.data$zscore.max.expr.value <- raw.data$zscore.max.expr.value + abs(min(raw.data$zscore.max.expr.value, na.rm = TRUE))
  # And get a score from 0 to 1 :
  raw.data$zscore.mean.expr.value <- raw.data$zscore.mean.expr.value / max(raw.data$zscore.mean.expr.value, na.rm = TRUE)
  raw.data$zscore.max.expr.value <- raw.data$zscore.max.expr.value / max(raw.data$zscore.max.expr.value, na.rm = TRUE)
  
  full.dataframe.tmp <- merge(normalized.pval, raw.data[, c("ID", "zscore.mean.expr.value", "zscore.max.expr.value", "mean.RNA.level", "max.RNA.level")], by="ID")
  full.dataframe.tmp$tissue <- rep(tissue, nrow(full.dataframe.tmp))
  
  full.dataframe <- rbind(full.dataframe, full.dataframe.tmp)
  
}

# Rhythmic transcripts p<0.01  /  Not Rhythmic transcripts p>0.5
RNA_rhythm.pvalue_threshold.inf = 0.01
RNA_rhythm.pvalue_threshold.sup = 0.5

full.dataframe$rhythmic <- full.dataframe$GeneCycle.pvalue <= RNA_rhythm.pvalue_threshold.inf
full.dataframe$non.rhythmic <- full.dataframe$GeneCycle.pvalue > RNA_rhythm.pvalue_threshold.sup

tot.data <- data.frame(ID=NA,
                       mean.rhythmic.tissues=NA, 
                       mean.non.rhythmic.tissues=NA,
                       max.rhythmic.tissues=NA, 
                       max.non.rhythmic.tissues=NA,
                       tau.mean=NA,
                       tau.max=NA, 
                       nb.rhythmic.tissues=NA)
for (i in 1:length(unique(full.dataframe$ID))) {
  gene.id <- unique(full.dataframe$ID)[i]
  sub.full.dataframe <- subset(full.dataframe, ID==gene.id)
  
  if (nrow(sub.full.dataframe) < 2) {
    next
  }
  
  tau.mean = ftau(sub.full.dataframe$zscore.mean.expr.value)
  tau.max = ftau(sub.full.dataframe$zscore.max.expr.value)
  
  non.rhythmic.tissues.list <- sub.full.dataframe[which(sub.full.dataframe$non.rhythmic==TRUE), "tissue"]
  rhythmic.tissues.list <- sub.full.dataframe[which(sub.full.dataframe$rhythmic==TRUE), "tissue"]
  if (length(non.rhythmic.tissues.list) != 0 && length(rhythmic.tissues.list) != 0) {
    mean.rhythmic.tissues.tmp <- mean(sub.full.dataframe[which(sub.full.dataframe$rhythmic==TRUE), "zscore.mean.expr.value"], na.rm=TRUE)
    mean.non.rhythmic.tissues.tmp <- mean(sub.full.dataframe[which(sub.full.dataframe$non.rhythmic==TRUE), "zscore.mean.expr.value"], na.rm=TRUE)
    
    max.rhythmic.tissues.tmp <- mean(sub.full.dataframe[which(sub.full.dataframe$rhythmic==TRUE), "zscore.max.expr.value"], na.rm=TRUE)
    max.non.rhythmic.tissues.tmp <- mean(sub.full.dataframe[which(sub.full.dataframe$non.rhythmic==TRUE), "zscore.max.expr.value"], na.rm=TRUE)
    
    tot.data <- rbind(tot.data, data.frame(ID=gene.id,
                                           mean.rhythmic.tissues=mean.rhythmic.tissues.tmp, 
                                           mean.non.rhythmic.tissues=mean.non.rhythmic.tissues.tmp,
                                           max.rhythmic.tissues=max.rhythmic.tissues.tmp, 
                                           max.non.rhythmic.tissues=max.non.rhythmic.tissues.tmp,
                                           tau.mean=tau.mean,
                                           tau.max=tau.max,
                                           nb.rhythmic.tissues=length(rhythmic.tissues.list)))
  } else {
    tot.data <- rbind(tot.data, data.frame(ID=gene.id,
                                           mean.rhythmic.tissues=NA, 
                                           mean.non.rhythmic.tissues=NA,
                                           max.rhythmic.tissues=NA, 
                                           max.non.rhythmic.tissues=NA,
                                           tau.mean=tau.mean,
                                           tau.max=tau.max,
                                           nb.rhythmic.tissues=length(rhythmic.tissues.list)))
  }
}
tot.data <- tot.data[-1, ]

full.dataframe.tot <- merge(full.dataframe, tot.data[, c("ID", "tau.mean", "tau.max", "nb.rhythmic.tissues")], by="ID")

# check
if (nrow(full.dataframe.tot) != nrow(full.dataframe)) {
  print("Warning ! there is an issue ...")
}



####################################
####################################
#######  MOUSE (Proteomics) 
####################################
####################################
species <- "mouse"
main.dir <- "~/Documents/cost_theory_workingspace/DATA/mouse"
######

####
#### With 3 tissues 
####
liver.proteome.data <- read.table(paste(main.dir, "liver/proteome/proteome_data.txt", sep="/"), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep="\t")
tendon.proteome.data <- read.table(paste(main.dir, "tendon/proteome/proteome_data.txt", sep="/"), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep="\t")
forebrain.proteome.data <- read.table(paste(main.dir, "forebrain/proteome/proteome_data.txt", sep="/"), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep="\t")

liver.proteome.data <- liver.proteome.data[, c("Protein.Name", "Gene.Name", "Uniprot.ID", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
liver.proteome.data$tissue <- rep("liver", nrow(liver.proteome.data))
tendon.proteome.data <- tendon.proteome.data[, c("Uniprot.ID", "Gene.Name", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
tendon.proteome.data$tissue <- rep("tendon", nrow(tendon.proteome.data))
forebrain.proteome.data <- forebrain.proteome.data[, c("Uniprot.ID", "Gene.Name", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
forebrain.proteome.data$tissue <- rep("forebrain", nrow(forebrain.proteome.data))

data.list <- list(liver.proteome.data, 
                  tendon.proteome.data,
                  forebrain.proteome.data)
for (i in 1:length(data.list)) {
  raw.dataset.values <- data.list[[i]] 
  proteome.data <- data.list[[i]]  
  ### Z-score normalization of each dataset
  m.mean = mean(raw.dataset.values$protein_mean.level, na.rm = TRUE)
  m.max = mean(raw.dataset.values$protein_max.level, na.rm = TRUE)
  sd.mean = sd(raw.dataset.values$protein_mean.level, na.rm = TRUE)
  sd.max = sd(raw.dataset.values$protein_max.level, na.rm = TRUE)
  
  proteome.data$zscore.mean.expr.value <- (raw.dataset.values$protein_mean.level - m.mean) / sd.mean
  proteome.data$zscore.max.expr.value <- (raw.dataset.values$protein_max.level - m.max) / sd.max
  # to get positive values to then calculate a coherent total expression cost
  proteome.data$zscore.mean.expr.value <- proteome.data$zscore.mean.expr.value + abs(min(proteome.data$zscore.mean.expr.value, na.rm = TRUE))
  proteome.data$zscore.max.expr.value <- proteome.data$zscore.max.expr.value + abs(min(proteome.data$zscore.max.expr.value, na.rm = TRUE))
  # And get a score from 0 to 1 :
  data.list[[i]]$zscore.mean.expr.value <- proteome.data$zscore.mean.expr.value / max(proteome.data$zscore.mean.expr.value, na.rm = TRUE)
  data.list[[i]]$zscore.max.expr.value <- proteome.data$zscore.max.expr.value / max(proteome.data$zscore.max.expr.value, na.rm = TRUE)
}
liver.proteome.data <- data.list[[1]]
tendon.proteome.data <- data.list[[2]]
forebrain.proteome.data <- data.list[[3]]

common.protein.data <- intersect(intersect(liver.proteome.data$Gene.Name, 
                                           tendon.proteome.data$Gene.Name),
                                 forebrain.proteome.data$Gene.Name)

liver.proteome.data <- unique(liver.proteome.data[liver.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level")]
liver.proteome.data$tissue <- rep("liver", nrow(liver.proteome.data))
tendon.proteome.data <- unique(tendon.proteome.data[tendon.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level")]
tendon.proteome.data$tissue <- rep("tendon", nrow(tendon.proteome.data))
forebrain.proteome.data <- unique(forebrain.proteome.data[forebrain.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level")]
forebrain.proteome.data$tissue <- rep("forebrain", nrow(forebrain.proteome.data))

full.protein.dataframe <- rbind(liver.proteome.data, 
                        tendon.proteome.data, 
                        forebrain.proteome.data)
# 
protein_rhythm.pvalue_threshold.inf = 0.02
protein_rhythm.pvalue_threshold.sup = 0.2

full.protein.dataframe$rhythmic <- full.protein.dataframe$protein_rhythm.pvalue <= protein_rhythm.pvalue_threshold.inf
full.protein.dataframe$non.rhythmic <- full.protein.dataframe$protein_rhythm.pvalue > protein_rhythm.pvalue_threshold.sup

gene.id <- read.table("~/Documents/DATA/Genes_ID/GeneID_GeneName_MOUSE.txt", h=T, fill=T)

full.protein.dataframe <- merge(gene.id, full.protein.dataframe, by="Gene.Name")
full.protein.dataframe <- merge(full.protein.dataframe, tot.data, by.x="Gene.ID", by.y="ID")
full.protein.dataframe <- unique(full.protein.dataframe[, c("Gene.Name", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level", "tissue", "rhythmic", "non.rhythmic", "tau.mean", "tau.max")])

# check
if (nrow(full.protein.dataframe) != nrow(full.protein.dataframe)) {
  print("Warning ! there is an issue ...")
}

for (i in 1:length(unique(full.protein.dataframe$tissue))) {
  sub.dataframe.tissue <- subset(full.protein.dataframe, tissue == unique(full.protein.dataframe$tissue)[i])
  
  # Replace p-values = 0 by the min of p-values divided by 10
  min.value = min(sub.dataframe.tissue[sub.dataframe.tissue$protein_rhythm.pvalue != 0, "protein_rhythm.pvalue"])
  sub.dataframe.tissue[sub.dataframe.tissue$protein_rhythm.pvalue == 0, "protein_rhythm.pvalue"] <- min.value/10
  # get only positives values of expression levels
  min.value = min(sub.dataframe.tissue$protein_mean.level)
  sub.dataframe.tissue$protein_mean.level <- sub.dataframe.tissue$protein_mean.level + abs(min.value)
  # Replace protein expr = 0 by the min of protein expr divided by 1000
  min.value = min(sub.dataframe.tissue[sub.dataframe.tissue$protein_mean.level != 0, "protein_mean.level"])
  sub.dataframe.tissue[sub.dataframe.tissue$protein_mean.level == 0, "protein_mean.level"] <- min.value/1000
  
  # See respective influence of gene expression and tissue-specificity into the rhythmic expression
  reg.lin <- lm(formula = -log10(protein_rhythm.pvalue) ~ log(protein_mean.level) + tau.mean, data = sub.dataframe.tissue)                      
  reg.lin
  summary.reg.lin <- summary(reg.lin)
  
  table.tot.tmp <- NULL
  table.tot.tmp <- data.frame(species = "mouse",
                              omics = "proteins", 
                              technique = "SILAC",
                              tissue = unique(full.protein.dataframe$tissue)[i],
                              `nb tissues` = length(unique(full.protein.dataframe$tissue)),
                              `nb genes` = length(unique(full.dataframe$Gene.Name)),
                              `rhythmic` = paste("p<", RNA_rhythm.pvalue_threshold.inf, 
                                                 " (", nrow(subset(sub.dataframe.tissue, rhythmic == TRUE)), " genes)", sep=""),
                              `non-rhythmic` = paste("p>", RNA_rhythm.pvalue_threshold.sup, 
                                                     " (", nrow(subset(sub.dataframe.tissue, non.rhythmic == TRUE)), " genes)", sep=""),
                              parameters = "-log10(rhythm.p-value) ~ log(mean.Protein.level) + tau")
  table.tot.tmp <- rbind(table.tot.tmp, NA, NA)
  table.tot.tmp$terms <- row.names(summary.reg.lin$coefficients)
  table.tot.tmp[, c("Estimate", "Std. Error", "t value", "Signif")] <- cbind(round(summary.reg.lin$coefficients[, 1:3], 3),
                                                                             signiFunction(summary.reg.lin$coefficients[, 4]))
  table.tot.tmp$`Adj R-squared` = round(summary.reg.lin$adj.r.squared, 3)
  table.tot.tmp$`F-statistic` = round(summary.reg.lin$fstatistic[1], 3)
  
  table.tot.1 <- rbind(table.tot.1, table.tot.tmp)
  
  
  # See influence of the interaction between gene expression level and tissue-specificity into the rhythmic expression
  reg.lin <- lm(formula = -log10(protein_rhythm.pvalue) ~ log(protein_mean.level) * tau.mean, data = sub.dataframe.tissue)                      
  reg.lin
  summary.reg.lin <- summary(reg.lin)
  
  table.tot.tmp <- NULL
  table.tot.tmp <- data.frame(species = "mouse",
                              omics = "proteins", 
                              technique = "SILAC",
                              tissue = unique(full.protein.dataframe$tissue)[i],
                              `nb tissues` = length(unique(full.protein.dataframe$tissue)),
                              `nb genes` = length(unique(full.dataframe$Gene.Name)),
                              `rhythmic` = paste("p<", RNA_rhythm.pvalue_threshold.inf, 
                                                 " (", nrow(subset(sub.dataframe.tissue, rhythmic == TRUE)), " genes)", sep=""),
                              `non-rhythmic` = paste("p>", RNA_rhythm.pvalue_threshold.sup, 
                                                     " (", nrow(subset(sub.dataframe.tissue, non.rhythmic == TRUE)), " genes)", sep=""),
                              parameters = "-log10(rhythm.p-value) ~ log(mean.Protein.level) * tau")
  table.tot.tmp <- rbind(table.tot.tmp, NA, NA, NA)
  table.tot.tmp$terms <- row.names(summary.reg.lin$coefficients)
  table.tot.tmp[, c("Estimate", "Std. Error", "t value", "Signif")] <- cbind(round(summary.reg.lin$coefficients[, 1:3], 3),
                                                                             signiFunction(summary.reg.lin$coefficients[, 4]))
  table.tot.tmp$`Adj R-squared` = round(summary.reg.lin$adj.r.squared, 3)
  table.tot.tmp$`F-statistic` = round(summary.reg.lin$fstatistic[1], 3)
  
  table.tot.2 <- rbind(table.tot.2, table.tot.tmp)
}



####
#### With 4 tissues 
####
liver.proteome.data <- read.table(paste(main.dir, "liver/proteome/proteome_data.txt", sep="/"), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep="\t")
cartilage.proteome.data <- read.table(paste(main.dir, "cartilage/proteome/proteome_data.txt", sep="/"), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep="\t")
tendon.proteome.data <- read.table(paste(main.dir, "tendon/proteome/proteome_data.txt", sep="/"), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep="\t")
forebrain.proteome.data <- read.table(paste(main.dir, "forebrain/proteome/proteome_data.txt", sep="/"), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE, sep="\t")

liver.proteome.data <- liver.proteome.data[, c("Protein.Name", "Gene.Name", "Uniprot.ID", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
liver.proteome.data$tissue <- rep("liver", nrow(liver.proteome.data))
tendon.proteome.data <- tendon.proteome.data[, c("Uniprot.ID", "Gene.Name", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
tendon.proteome.data$tissue <- rep("tendon", nrow(tendon.proteome.data))
forebrain.proteome.data <- forebrain.proteome.data[, c("Uniprot.ID", "Gene.Name", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
forebrain.proteome.data$tissue <- rep("forebrain", nrow(forebrain.proteome.data))
cartilage.proteome.data <- cartilage.proteome.data[, c("Uniprot.ID", "Gene.Name", "protein_mean.level", "protein_max.level", "protein_rhythm.pvalue")]
cartilage.proteome.data$tissue <- rep("tendon", nrow(cartilage.proteome.data))

data.list <- list(liver.proteome.data,
                  tendon.proteome.data,
                  forebrain.proteome.data,
                  cartilage.proteome.data)
for (i in 1:4) {
  raw.dataset.values <- data.list[[i]] 
  proteome.data <- data.list[[i]]  
  ### Z-score normalization of each dataset
  m.mean = mean(raw.dataset.values$protein_mean.level, na.rm = TRUE)
  m.max = mean(raw.dataset.values$protein_max.level, na.rm = TRUE)
  sd.mean = sd(raw.dataset.values$protein_mean.level, na.rm = TRUE)
  sd.max = sd(raw.dataset.values$protein_max.level, na.rm = TRUE)
  
  proteome.data$zscore.mean.expr.value <- (raw.dataset.values$protein_mean.level - m.mean) / sd.mean
  proteome.data$zscore.max.expr.value <- (raw.dataset.values$protein_max.level - m.max) / sd.max
  # to get positive values to then calculate a coherent total expression cost
  proteome.data$zscore.mean.expr.value <- proteome.data$zscore.mean.expr.value + abs(min(proteome.data$zscore.mean.expr.value, na.rm = TRUE))
  proteome.data$zscore.max.expr.value <- proteome.data$zscore.max.expr.value + abs(min(proteome.data$zscore.max.expr.value, na.rm = TRUE))
  # And get a score from 0 to 1 :
  data.list[[i]]$zscore.mean.expr.value <- proteome.data$zscore.mean.expr.value / max(proteome.data$zscore.mean.expr.value, na.rm = TRUE)
  data.list[[i]]$zscore.max.expr.value <- proteome.data$zscore.max.expr.value / max(proteome.data$zscore.max.expr.value, na.rm = TRUE)
}
liver.proteome.data <- data.list[[1]]
tendon.proteome.data <- data.list[[2]]
forebrain.proteome.data <- data.list[[3]]
cartilage.proteome.data <- data.list[[4]]

common.protein.data <- unique(intersect(intersect(intersect(liver.proteome.data$Gene.Name, 
                                                            tendon.proteome.data$Gene.Name),
                                                  forebrain.proteome.data$Gene.Name),
                                        cartilage.proteome.data$Gene.Name))

liver.proteome.data <- unique(liver.proteome.data[liver.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level")]
liver.proteome.data$tissue <- rep("liver", nrow(liver.proteome.data))
tendon.proteome.data <- unique(tendon.proteome.data[tendon.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level")]
tendon.proteome.data$tissue <- rep("tendon", nrow(tendon.proteome.data))
forebrain.proteome.data <- unique(forebrain.proteome.data[forebrain.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level")]
forebrain.proteome.data$tissue <- rep("forebrain", nrow(forebrain.proteome.data))
cartilage.proteome.data <- unique(cartilage.proteome.data[cartilage.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level")]
cartilage.proteome.data$tissue <- rep("cartilage", nrow(cartilage.proteome.data))

full.protein.dataframe <- rbind(liver.proteome.data, 
                        tendon.proteome.data, 
                        forebrain.proteome.data,
                        cartilage.proteome.data)

# 
protein_rhythm.pvalue_threshold.inf = 0.02
protein_rhythm.pvalue_threshold.sup = 0.2

full.protein.dataframe$rhythmic <- full.protein.dataframe$protein_rhythm.pvalue <= protein_rhythm.pvalue_threshold.inf
full.protein.dataframe$non.rhythmic <- full.protein.dataframe$protein_rhythm.pvalue > protein_rhythm.pvalue_threshold.sup

gene.id <- read.table("~/Documents/DATA/Genes_ID/GeneID_GeneName_MOUSE.txt", h=T, fill=T)

full.protein.dataframe <- merge(gene.id, full.protein.dataframe, by="Gene.Name")
full.protein.dataframe <- merge(full.protein.dataframe, tot.data, by.x="Gene.ID", by.y="ID")
full.protein.dataframe <- unique(full.protein.dataframe[, c("Gene.Name", "protein_rhythm.pvalue", "protein_mean.level", "protein_max.level", "tissue", "rhythmic", "non.rhythmic", "tau.mean", "tau.max")])

# check
if (nrow(full.protein.dataframe) != nrow(full.protein.dataframe)) {
  print("Warning ! there is an issue ...")
}

for (i in 1:length(unique(full.protein.dataframe$tissue))) {
  sub.dataframe.tissue <- subset(full.protein.dataframe, tissue == unique(full.protein.dataframe$tissue)[i])
  
  # Replace p-values = 0 by the min of p-values divided by 10
  min.value = min(sub.dataframe.tissue[sub.dataframe.tissue$protein_rhythm.pvalue != 0, "protein_rhythm.pvalue"])
  sub.dataframe.tissue[sub.dataframe.tissue$protein_rhythm.pvalue == 0, "protein_rhythm.pvalue"] <- min.value/10
  # get only positives values of expression levels
  min.value = min(sub.dataframe.tissue$protein_mean.level)
  sub.dataframe.tissue$protein_mean.level <- sub.dataframe.tissue$protein_mean.level + abs(min.value)
  # Replace protein expr = 0 by the min of protein expr divided by 1000
  min.value = min(sub.dataframe.tissue[sub.dataframe.tissue$protein_mean.level != 0, "protein_mean.level"])
  sub.dataframe.tissue[sub.dataframe.tissue$protein_mean.level == 0, "protein_mean.level"] <- min.value/1000
  
  # See respective influence of gene expression and tissue-specificity into the rhythmic expression
  reg.lin <- lm(formula = -log10(protein_rhythm.pvalue) ~ log(protein_mean.level) + tau.mean, data = sub.dataframe.tissue)                      
  reg.lin
  summary.reg.lin <- summary(reg.lin)
  
  table.tot.tmp <- NULL
  table.tot.tmp <- data.frame(species = "mouse",
                              omics = "proteins", 
                              technique = "SILAC",
                              tissue = unique(full.protein.dataframe$tissue)[i],
                              `nb tissues` = length(unique(full.protein.dataframe$tissue)),
                              `nb genes` = length(unique(full.dataframe$Gene.Name)),
                              `rhythmic` = paste("p<", RNA_rhythm.pvalue_threshold.inf, 
                                                 " (", nrow(subset(sub.dataframe.tissue, rhythmic == TRUE)), " genes)", sep=""),
                              `non-rhythmic` = paste("p>", RNA_rhythm.pvalue_threshold.sup, 
                                                     " (", nrow(subset(sub.dataframe.tissue, non.rhythmic == TRUE)), " genes)", sep=""),
                              parameters = "-log10(rhythm.p-value) ~ log(mean.Protein.level) + tau")
  table.tot.tmp <- rbind(table.tot.tmp, NA, NA)
  table.tot.tmp$terms <- row.names(summary.reg.lin$coefficients)
  table.tot.tmp[, c("Estimate", "Std. Error", "t value", "Signif")] <- cbind(round(summary.reg.lin$coefficients[, 1:3], 3),
                                                                             signiFunction(summary.reg.lin$coefficients[, 4]))
  table.tot.tmp$`Adj R-squared` = round(summary.reg.lin$adj.r.squared, 3)
  table.tot.tmp$`F-statistic` = round(summary.reg.lin$fstatistic[1], 3)
  
  table.tot.1 <- rbind(table.tot.1, table.tot.tmp)
  
  
  # See influence of the interaction between gene expression level and tissue-specificity into the rhythmic expression
  reg.lin <- lm(formula = -log10(protein_rhythm.pvalue) ~ log(protein_mean.level) * tau.mean, data = sub.dataframe.tissue)                      
  reg.lin
  summary.reg.lin <- summary(reg.lin)
  
  table.tot.tmp <- NULL
  table.tot.tmp <- data.frame(species = "mouse",
                              omics = "proteins", 
                              technique = "SILAC",
                              tissue = unique(full.protein.dataframe$tissue)[i],
                              `nb tissues` = length(unique(full.protein.dataframe$tissue)),
                              `nb genes` = length(unique(full.dataframe$Gene.Name)),
                              `rhythmic` = paste("p<", RNA_rhythm.pvalue_threshold.inf, 
                                                 " (", nrow(subset(sub.dataframe.tissue, rhythmic == TRUE)), " genes)", sep=""),
                              `non-rhythmic` = paste("p>", RNA_rhythm.pvalue_threshold.sup, 
                                                     " (", nrow(subset(sub.dataframe.tissue, non.rhythmic == TRUE)), " genes)", sep=""),
                              parameters = "-log10(rhythm.p-value) ~ log(mean.Protein.level) * tau")
  table.tot.tmp <- rbind(table.tot.tmp, NA, NA, NA)
  table.tot.tmp$terms <- row.names(summary.reg.lin$coefficients)
  table.tot.tmp[, c("Estimate", "Std. Error", "t value", "Signif")] <- cbind(round(summary.reg.lin$coefficients[, 1:3], 3),
                                                                             signiFunction(summary.reg.lin$coefficients[, 4]))
  table.tot.tmp$`Adj R-squared` = round(summary.reg.lin$adj.r.squared, 3)
  table.tot.tmp$`F-statistic` = round(summary.reg.lin$fstatistic[1], 3)
  
  table.tot.2 <- rbind(table.tot.2, table.tot.tmp)
}



write.csv(table.tot.1, "~/Documents/cost_theory_workingspace/TableS3a_bis.csv", quote = F, row.names = FALSE)
write.csv(table.tot.2, "~/Documents/cost_theory_workingspace/TableS3b_bis.csv", quote = F, row.names = FALSE)
library(xtable)
print(xtable(table.tot.1, type = "latex"), file = "~/Documents/cost_theory_workingspace/TableS3a_bis.tex", include.rownames=FALSE)
print(xtable(table.tot.2, type = "latex"), file = "~/Documents/cost_theory_workingspace/TableS3b_bis.tex", include.rownames=FALSE)


