# *** functions *** #
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


table.tot <- NULL
delta.values <- NULL

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
  
  full.dataframe.tmp <- merge(normalized.pval, raw.data[, c("ID", "zscore.mean.expr.value", "zscore.max.expr.value")], by="ID")
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
                       max.non.rhythmic.tissues=NA)
for (i in 1:length(unique(full.dataframe$ID))) {
  gene.id <- unique(full.dataframe$ID)[i]
  sub.full.dataframe <- subset(full.dataframe, ID==gene.id)
  
  if (nrow(sub.full.dataframe) < 2) {
    next
  }
  
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
                                           max.non.rhythmic.tissues=max.non.rhythmic.tissues.tmp))
  }
}
tot.data <- tot.data[-1, ]

#tot.data$delta.mean <- tot.data$mean.rhythmic.tissues - tot.data$mean.non.rhythmic.tissues
tot.data$delta.mean <- tot.data$max.rhythmic.tissues - tot.data$max.non.rhythmic.tissues
tot.data <- tot.data[,-1]
tot.data$tissues.nb <- rep(paste(length(tissue.list), " tissues", "\n",
                                 nrow(tot.data), " transcripts", sep=""), nrow(tot.data))

#hist(tot.data$delta.mean, breaks = 150)
ttest <- t.test(tot.data$delta.mean, mu=0)
ttest

table.tot <- rbind(table.tot, data.frame(species = "mouse",
                                         omics = "transcript", 
                                         technique = "microarray",
                                         `nb tissues` = length(tissue.list),
                                         `nb genes` = length(unique(full.dataframe$ID)),
                                         parameter = "max expression level",
                                         `rhythmic` = paste("p<", RNA_rhythm.pvalue_threshold.inf, sep=""),
                                         `non-rhythmic` = paste("p>", RNA_rhythm.pvalue_threshold.sup, sep=""),
                                         `Mean of Delta tissues` = round(ttest$estimate[1], 4),
                                         t = round(ttest$statistic, 2),
                                         df = round(ttest$parameter, 2),
                                         lower = round(ttest$conf.int[1], 4),
                                         upper = round(ttest$conf.int[2], 4), 
                                         Signif = signiFunction(ttest$p.value)
                                         #,trend = trendFunction(obs.value = ttest$estimate[1], theoretic.value = ttest$null.value, signif = ttest$p.value)
                                         , check.names = FALSE))
delta.values <- rbind(delta.values, 
                      data.frame(species = "mouse, RNA, 11 tissues",
                                 delta = tot.data$delta.mean))



####################################
####################################
#######  DROSOPHILA 
####################################
####################################
species <- "drosophila"
tissue.list <- c("body", "head", "heart")
p.val <- "default.pvalue"
######

full.dataframe <- data.frame()
for (t in 1:length(tissue.list)) {
  tissue <- tissue.list[t]
  
  file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")
  
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
  
  full.dataframe.tmp <- merge(normalized.pval, raw.data[, c("ID", "zscore.mean.expr.value", "zscore.max.expr.value")], by="ID")
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
                       max.non.rhythmic.tissues=NA)
for (i in 1:length(unique(full.dataframe$ID))) {
  gene.id <- unique(full.dataframe$ID)[i]
  sub.full.dataframe <- subset(full.dataframe, ID==gene.id)
  
  if (nrow(sub.full.dataframe) < 2) {
    next
  }
  
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
                                           max.non.rhythmic.tissues=max.non.rhythmic.tissues.tmp))
  }
}
tot.data <- tot.data[-1, ]

#tot.data$delta.mean <- tot.data$mean.rhythmic.tissues - tot.data$mean.non.rhythmic.tissues
tot.data$delta.mean <- tot.data$max.rhythmic.tissues - tot.data$max.non.rhythmic.tissues

tot.data <- tot.data[,-1]
tot.data$tissues.nb <- rep(paste(length(tissue.list), " tissues", "\n",
                                 nrow(tot.data), " transcripts", sep=""), nrow(tot.data))
#hist(tot.data$delta.mean, breaks = 150)
ttest <- t.test(tot.data$delta.mean, mu=0)
ttest

table.tot <- rbind(table.tot, data.frame(species = "drosophila",
                                         omics = "transcript", 
                                         technique = "RNAseq",
                                         `nb tissues` = length(tissue.list),
                                         `nb genes` = length(unique(full.dataframe$ID)),
                                         parameter = "max expression level",
                                         `rhythmic` = paste("p<", RNA_rhythm.pvalue_threshold.inf, sep=""),
                                         `non-rhythmic` = paste("p>", RNA_rhythm.pvalue_threshold.sup, sep=""),
                                         `Mean of Delta tissues` = round(ttest$estimate[1], 4),
                                         t = round(ttest$statistic, 2),
                                         df = round(ttest$parameter, 2),
                                         lower = round(ttest$conf.int[1], 4),
                                         upper = round(ttest$conf.int[2], 4), 
                                         Signif = signiFunction(ttest$p.value)
                                         #,trend = trendFunction(obs.value = ttest$estimate[1], theoretic.value = ttest$null.value, signif = ttest$p.value)
                                         , check.names = FALSE))

delta.values <- rbind(delta.values, 
                      data.frame(species = "drosophila, RNA, 3 tissues",
                                 delta = tot.data$delta.mean))







####################################
####################################
#######  ANOPHELES 
####################################
####################################
species <- "anopheles"
tissue.list <- c("body", "head")
p.val <- "default.pvalue"
######

full.dataframe <- data.frame()
for (t in 1:length(tissue.list)) {
  tissue <- tissue.list[t]
  
  file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")
  
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
  
  full.dataframe.tmp <- merge(normalized.pval, raw.data[, c("ID", "zscore.mean.expr.value", "zscore.max.expr.value")], by="ID")
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
                       max.non.rhythmic.tissues=NA)
for (i in 1:length(unique(full.dataframe$ID))) {
  gene.id <- unique(full.dataframe$ID)[i]
  sub.full.dataframe <- subset(full.dataframe, ID==gene.id)
  
  if (nrow(sub.full.dataframe) < 2) {
    next
  }
  
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
                                           max.non.rhythmic.tissues=max.non.rhythmic.tissues.tmp))
  }
}
tot.data <- tot.data[-1, ]

#tot.data$delta.mean <- tot.data$mean.rhythmic.tissues - tot.data$mean.non.rhythmic.tissues
tot.data$delta.mean <- tot.data$max.rhythmic.tissues - tot.data$max.non.rhythmic.tissues

tot.data <- tot.data[,-1]
tot.data$tissues.nb <- rep(paste(length(tissue.list), " tissues", "\n",
                                 nrow(tot.data), " transcripts", sep=""), nrow(tot.data))
#hist(tot.data$delta.mean, breaks = 150)
ttest <- t.test(tot.data$delta.mean, mu=0)
ttest

table.tot <- rbind(table.tot, data.frame(species = "anopheles",
                                         omics = "transcript", 
                                         technique = "microarray",
                                         `nb tissues` = length(tissue.list),
                                         `nb genes` = length(unique(full.dataframe$ID)),
                                         parameter = "max expression level",
                                         `rhythmic` = paste("p<", RNA_rhythm.pvalue_threshold.inf, sep=""),
                                         `non-rhythmic` = paste("p>", RNA_rhythm.pvalue_threshold.sup, sep=""),
                                         `Mean of Delta tissues` = round(ttest$estimate[1], 4),
                                         t = round(ttest$statistic, 2),
                                         df = round(ttest$parameter, 2),
                                         lower = round(ttest$conf.int[1], 4),
                                         upper = round(ttest$conf.int[2], 4), 
                                         Signif = signiFunction(ttest$p.value)
                                         #, trend = trendFunction(obs.value = ttest$estimate[1], theoretic.value = ttest$null.value, signif = ttest$p.value)
                                         , check.names = FALSE))

delta.values <- rbind(delta.values, 
                      data.frame(species = "anopheles, RNA, 2 tissues",
                                 delta = tot.data$delta.mean))





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

liver.proteome.data <- unique(liver.proteome.data[liver.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue")]
liver.proteome.data$tissue <- rep("liver", nrow(liver.proteome.data))
tendon.proteome.data <- unique(tendon.proteome.data[tendon.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue")]
tendon.proteome.data$tissue <- rep("tendon", nrow(tendon.proteome.data))
forebrain.proteome.data <- unique(forebrain.proteome.data[forebrain.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue")]
forebrain.proteome.data$tissue <- rep("forebrain", nrow(forebrain.proteome.data))

full.dataframe <- rbind(liver.proteome.data, 
                        tendon.proteome.data, 
                        forebrain.proteome.data)
  
# 
protein_rhythm.pvalue_threshold.inf = 0.05
protein_rhythm.pvalue_threshold.sup = 0.2

full.dataframe$rhythmic <- full.dataframe$protein_rhythm.pvalue <= protein_rhythm.pvalue_threshold.inf
full.dataframe$non.rhythmic <- full.dataframe$protein_rhythm.pvalue > protein_rhythm.pvalue_threshold.sup

tot.data <- data.frame(ID=NA,
                       mean.rhythmic.tissues=NA, 
                       mean.non.rhythmic.tissues=NA,
                       max.rhythmic.tissues=NA, 
                       max.non.rhythmic.tissues=NA)
for (i in 1:length(unique(full.dataframe$Gene.Name))) {
  gene.id <- unique(full.dataframe$Gene.Name)[i]
  sub.full.dataframe <- subset(full.dataframe, Gene.Name == gene.id)
  
  if (nrow(sub.full.dataframe) < 2) {
    next
  }
  
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
                                           max.non.rhythmic.tissues=max.non.rhythmic.tissues.tmp))
  }
}
tot.data <- tot.data[-1, ]

#tot.data$delta.mean <- tot.data$mean.rhythmic.tissues - tot.data$mean.non.rhythmic.tissues
tot.data$delta.mean <- tot.data$max.rhythmic.tissues - tot.data$max.non.rhythmic.tissues
tot.data <- tot.data[,-1]
tot.data$tissues.nb <- rep(paste("3 tissues", "\n",
                                 length(unique(full.dataframe$Gene.Name)), " genes", sep=""), nrow(tot.data))

#hist(tot.data$delta.mean, breaks = 150)
ttest <- t.test(tot.data$delta.mean, mu=0)
ttest

table.tot <- rbind(table.tot, data.frame(species = "mouse",
                                         omics = "protein", 
                                         technique = "SILAC",
                                         `nb tissues` = "3",
                                         `nb genes` = length(unique(full.dataframe$Gene.Name)),
                                         parameter = "max expression level",
                                         `rhythmic` = paste("p<", protein_rhythm.pvalue_threshold.inf, sep=""),
                                         `non-rhythmic` = paste("p>", protein_rhythm.pvalue_threshold.sup, sep=""),
                                         `Mean of Delta tissues` = round(ttest$estimate[1], 4),
                                         t = round(ttest$statistic, 2),
                                         df = round(ttest$parameter, 2),
                                         lower = round(ttest$conf.int[1], 4),
                                         upper = round(ttest$conf.int[2], 4), 
                                         Signif = signiFunction(ttest$p.value)
                                         #,trend = trendFunction(obs.value = ttest$estimate[1], theoretic.value = ttest$null.value, signif = ttest$p.value)
                                         , check.names = FALSE))

delta.values <- rbind(delta.values, 
                      data.frame(species = "mouse, protein, 3 tissues",
                                 delta = tot.data$delta.mean))


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

liver.proteome.data <- unique(liver.proteome.data[liver.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue")]
liver.proteome.data$tissue <- rep("liver", nrow(liver.proteome.data))
tendon.proteome.data <- unique(tendon.proteome.data[tendon.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue")]
tendon.proteome.data$tissue <- rep("tendon", nrow(tendon.proteome.data))
forebrain.proteome.data <- unique(forebrain.proteome.data[forebrain.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue")]
forebrain.proteome.data$tissue <- rep("forebrain", nrow(forebrain.proteome.data))
cartilage.proteome.data <- unique(cartilage.proteome.data[cartilage.proteome.data$Gene.Name %in% common.protein.data, ])[, c("Gene.Name", "zscore.mean.expr.value", "zscore.max.expr.value", "protein_rhythm.pvalue")]
cartilage.proteome.data$tissue <- rep("cartilage", nrow(cartilage.proteome.data))

full.dataframe <- rbind(liver.proteome.data, 
                        tendon.proteome.data, 
                        forebrain.proteome.data,
                        cartilage.proteome.data)

protein_rhythm.pvalue_threshold.inf = 0.05
protein_rhythm.pvalue_threshold.sup = 0.2

full.dataframe$rhythmic <- full.dataframe$protein_rhythm.pvalue <= protein_rhythm.pvalue_threshold.inf
full.dataframe$non.rhythmic <- full.dataframe$protein_rhythm.pvalue > protein_rhythm.pvalue_threshold.sup

tot.data <- data.frame(ID=NA,
                       mean.rhythmic.tissues=NA, 
                       mean.non.rhythmic.tissues=NA,
                       max.rhythmic.tissues=NA, 
                       max.non.rhythmic.tissues=NA)
for (i in 1:length(unique(full.dataframe$Gene.Name))) {
  gene.id <- unique(full.dataframe$Gene.Name)[i]
  sub.full.dataframe <- subset(full.dataframe, Gene.Name==gene.id)
  
  if (nrow(sub.full.dataframe) < 2) {
    next
  }
  
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
                                           max.non.rhythmic.tissues=max.non.rhythmic.tissues.tmp))
  }
}
tot.data <- tot.data[-1, ]

#tot.data$delta.mean <- tot.data$mean.rhythmic.tissues - tot.data$mean.non.rhythmic.tissues
tot.data$delta.mean <- tot.data$max.rhythmic.tissues - tot.data$max.non.rhythmic.tissues
tot.data <- tot.data[,-1]
tot.data$tissues.nb <- rep(paste("4 tissues", "\n",
                                 length(unique(full.dataframe$Gene.Name)), " genes", sep=""), nrow(tot.data))

#hist(tot.data$delta.mean, breaks = 150)
ttest <- t.test(tot.data$delta.mean, mu=0)
ttest

table.tot <- rbind(table.tot, data.frame(species = "mouse",
                                         omics = "protein", 
                                         technique = "SILAC",
                                         `nb tissues` = "4",
                                         `nb genes` = length(unique(full.dataframe$Gene.Name)),
                                         parameter = "max expression level",
                                         `rhythmic` = paste("p<", protein_rhythm.pvalue_threshold.inf, sep=""),
                                         `non-rhythmic` = paste("p>", protein_rhythm.pvalue_threshold.sup, sep=""),
                                         `Mean of Delta tissues` = round(ttest$estimate[1], 4),
                                         t = round(ttest$statistic, 2),
                                         df = round(ttest$parameter, 2),
                                         lower = round(ttest$conf.int[1], 4),
                                         upper = round(ttest$conf.int[2], 4), 
                                         Signif = signiFunction(ttest$p.value)
                                         #,trend = trendFunction(obs.value = ttest$estimate[1], theoretic.value = ttest$null.value, signif = ttest$p.value)
                                         , check.names = FALSE))

delta.values <- rbind(delta.values, 
                      data.frame(species = "mouse, protein, 4 tissues",
                                 delta = tot.data$delta.mean))

write.csv(table.tot, "~/Documents/cost_theory_workingspace/Table1.csv", quote = F, row.names = FALSE)
library(xtable)
print(xtable(table.tot, type = "latex", digits=4), file = "~/Documents/cost_theory_workingspace/Table1.tex", include.rownames=FALSE)


##### Histograms
# Ordering plots : 
ordered.group <- c("mouse, RNA, 11 tissues",
                   "mouse, protein, 3 tissues",
                   "anopheles, RNA, 2 tissues",
                   "mouse, protein, 4 tissues",
                   "drosophila, RNA, 3 tissues")
delta.values$species <- factor(delta.values$species, levels=ordered.group)

HistogramsPlots <- ggplot(delta.values, aes(x=delta, fill=species)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', 
                 alpha=0.6, na.rm = TRUE, bins = 60) + 
  facet_wrap( ~ species, scales = "free", nrow = 3) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  labs(x = TeX("$\\textbf{\\delta}$"), fill="genes group") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.position = "none",
        plot.caption = element_text(size=11, face = "italic", hjust=0))

