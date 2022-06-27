##################################
# Z-score normalization of data #
##################################
# Normalize mean and max expression values (over time-points) to obtain a null distribution
# With z-score
# ........
# ........
# ........
########################

species.list <- c("arabidopsis", "cyanobacteria", "ostreococcus")
tissue.list <- c("leaves", NA, NA)

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species.list, tissue.list, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

for (y in 1:length(species.list)) {
  proteome.data <- read.table(paste(proteome.dir[y], "protein_level.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  transcriptome.data <- read.table(paste(transcriptome.dir[y], "transcript_level.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  
  raw.dataset.values <- proteome.data[, grep("LL|LD|ZT|CT", colnames(proteome.data))]
   
  nb.timepoints <- ncol(raw.dataset.values)
  m.mean = apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  M.mean = mean(m.mean, na.rm = TRUE)
  sd.mean = sd(m.mean, na.rm = TRUE)
  zscore.raw.dataset.values <- raw.dataset.values
  for (i in 1:nrow(raw.dataset.values)) {
    for (j in 1:ncol(raw.dataset.values)) {
      zscore.raw.dataset.values[i, j] = (nb.timepoints*raw.dataset.values[i, j] - M.mean)/sd.mean
    }
  }
  zscore.raw.dataset.values <- zscore.raw.dataset.values + abs(min(zscore.raw.dataset.values, na.rm = T))
  #zscore.raw.dataset.values$mean <- apply(zscore.raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  #raw.dataset.values$mean <- apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  #plot(zscore.raw.dataset.values$mean, raw.dataset.values$mean)
  
  #zscore.raw.dataset.values <- zscore.raw.dataset.values / max(zscore.raw.dataset.values, na.rm = TRUE)
  zscore.raw.dataset.values <- cbind(as.data.frame(proteome.data[, -grep("LL|LD|ZT|CT", colnames(proteome.data))]),
                                     zscore.raw.dataset.values)
  colnames(zscore.raw.dataset.values)[-grep("LL|LD|ZT|CT", colnames(proteome.data))] <- colnames(proteome.data)[-grep("LL|LD|ZT|CT", colnames(proteome.data))]
  
  write.table(zscore.raw.dataset.values, paste(proteome.dir[y], "protein_zscore_level.txt", sep="/"), sep = "\t", quote = F, row.names = F)
  
  raw.dataset.values <- transcriptome.data[, grep("LL|LD|ZT|CT", colnames(transcriptome.data))]
  
  nb.timepoints <- ncol(raw.dataset.values)
  m.mean = apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  M.mean = mean(m.mean, na.rm = TRUE)
  sd.mean = sd(m.mean, na.rm = TRUE)
  zscore.raw.dataset.values <- raw.dataset.values
  for (i in 1:nrow(raw.dataset.values)) {
    for (j in 1:ncol(raw.dataset.values)) {
      zscore.raw.dataset.values[i, j] = (nb.timepoints*raw.dataset.values[i, j] - M.mean)/sd.mean
    }
  }
  zscore.raw.dataset.values <- zscore.raw.dataset.values + abs(min(zscore.raw.dataset.values, na.rm = T))
  #zscore.raw.dataset.values$mean <- apply(zscore.raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  #raw.dataset.values$mean <- apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  #plot(zscore.raw.dataset.values$mean, raw.dataset.values$mean)
  
  #zscore.raw.dataset.values <- zscore.raw.dataset.values / max(zscore.raw.dataset.values, na.rm = TRUE)
  zscore.raw.dataset.values <- cbind(as.data.frame(transcriptome.data[, -grep("LL|LD|ZT|CT", colnames(transcriptome.data))]),
                                     zscore.raw.dataset.values)
  colnames(zscore.raw.dataset.values)[-grep("LL|LD|ZT|CT", colnames(transcriptome.data))] <- colnames(transcriptome.data)[-grep("LL|LD|ZT|CT", colnames(transcriptome.data))]
  
  write.table(zscore.raw.dataset.values, paste(transcriptome.dir[y], "transcript_zscore_level.txt", sep="/"), sep = "\t", quote = F, row.names = F)
}


#### special case for mouse liver 
species.list <- "mouse"
tissue.list <- "liver"

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species.list, tissue.list, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

for (y in 1:length(species.list)) {
  proteome.data <- read.table(paste(proteome.dir[y], "raw_protein_level.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  transcriptome.data <- read.table(paste(transcriptome.dir[y], "transcript_level.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  
  raw.dataset.values <- proteome.data[, grep("LL|LD|ZT|CT", colnames(proteome.data))]
  
  nb.timepoints <- ncol(raw.dataset.values)
  m.mean = apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  M.mean = mean(m.mean, na.rm = TRUE)
  sd.mean = sd(m.mean, na.rm = TRUE)
  zscore.raw.dataset.values <- raw.dataset.values
  for (i in 1:nrow(raw.dataset.values)) {
    for (j in 1:ncol(raw.dataset.values)) {
      zscore.raw.dataset.values[i, j] = (nb.timepoints*raw.dataset.values[i, j] - M.mean)/sd.mean
    }
  }
  zscore.raw.dataset.values <- zscore.raw.dataset.values + abs(min(zscore.raw.dataset.values, na.rm = T))
  #zscore.raw.dataset.values$mean <- apply(zscore.raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  #raw.dataset.values$mean <- apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  #plot(zscore.raw.dataset.values$mean, raw.dataset.values$mean)
  
  zscore.raw.dataset.values <- zscore.raw.dataset.values / max(zscore.raw.dataset.values, na.rm = TRUE)
  zscore.raw.dataset.values <- cbind(as.data.frame(proteome.data[, -grep("LL|LD|ZT|CT", colnames(proteome.data))]),
                                     zscore.raw.dataset.values)
  colnames(zscore.raw.dataset.values)[-grep("LL|LD|ZT|CT", colnames(proteome.data))] <- colnames(proteome.data)[-grep("LL|LD|ZT|CT", colnames(proteome.data))]
  
  write.table(zscore.raw.dataset.values, paste(proteome.dir[y], "protein_zscore_level.txt", sep="/"), sep = "\t", quote = F, row.names = F)
  
  raw.dataset.values <- transcriptome.data[, grep("LL|LD|ZT|CT", colnames(transcriptome.data))]
  
  nb.timepoints <- ncol(raw.dataset.values)
  m.mean = apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  M.mean = mean(m.mean, na.rm = TRUE)
  sd.mean = sd(m.mean, na.rm = TRUE)
  zscore.raw.dataset.values <- raw.dataset.values
  for (i in 1:nrow(raw.dataset.values)) {
    for (j in 1:ncol(raw.dataset.values)) {
      zscore.raw.dataset.values[i, j] = (nb.timepoints*raw.dataset.values[i, j] - M.mean)/sd.mean
    }
  }
  zscore.raw.dataset.values <- zscore.raw.dataset.values + abs(min(zscore.raw.dataset.values, na.rm = T))
  #zscore.raw.dataset.values$mean <- apply(zscore.raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  #raw.dataset.values$mean <- apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  #plot(zscore.raw.dataset.values$mean, raw.dataset.values$mean)
  
  zscore.raw.dataset.values <- zscore.raw.dataset.values / max(zscore.raw.dataset.values, na.rm = TRUE)
  zscore.raw.dataset.values <- cbind(as.data.frame(transcriptome.data[, -grep("LL|LD|ZT|CT", colnames(transcriptome.data))]),
                                     zscore.raw.dataset.values)
  colnames(zscore.raw.dataset.values)[-grep("LL|LD|ZT|CT", colnames(transcriptome.data))] <- colnames(transcriptome.data)[-grep("LL|LD|ZT|CT", colnames(transcriptome.data))]
  
  write.table(zscore.raw.dataset.values, paste(transcriptome.dir[y], "transcript_zscore_level.txt", sep="/"), sep = "\t", quote = F, row.names = F)
}




### Idem other tissues proteomic data of Mouse:
species.list <- c("mouse", "mouse", "mouse")
tissue.list <- c("cartilage", "tendon", "forebrain")

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species.list, tissue.list, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")

for (y in 1:length(species.list)) {
  proteome.data <- read.table(paste(proteome.dir[y], "protein_level.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  
  raw.dataset.values <- proteome.data[, grep("LL|LD|ZT|CT", colnames(proteome.data))]
  
  nb.timepoints <- ncol(raw.dataset.values)
  m.mean = apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
  M.mean = mean(m.mean, na.rm = TRUE)
  sd.mean = sd(m.mean, na.rm = TRUE)
  zscore.raw.dataset.values <- raw.dataset.values
  for (i in 1:nrow(raw.dataset.values)) {
    for (j in 1:ncol(raw.dataset.values)) {
      zscore.raw.dataset.values[i, j] = (nb.timepoints*raw.dataset.values[i, j] - M.mean)/sd.mean
    }
  }
  zscore.raw.dataset.values <- zscore.raw.dataset.values + abs(min(zscore.raw.dataset.values, na.rm = T))
  zscore.raw.dataset.values <- zscore.raw.dataset.values / max(zscore.raw.dataset.values, na.rm = TRUE)
  zscore.raw.dataset.values <- cbind(as.data.frame(proteome.data[, -grep("LL|LD|ZT|CT", colnames(proteome.data))]),
                                     zscore.raw.dataset.values)
  colnames(zscore.raw.dataset.values)[-grep("LL|LD|ZT|CT", colnames(proteome.data))] <- colnames(proteome.data)[-grep("LL|LD|ZT|CT", colnames(proteome.data))]
  
  write.table(zscore.raw.dataset.values, paste(proteome.dir[y], "protein_zscore_level.txt", sep="/"), sep = "\t", quote = F, row.names = F)
}


### Idem for epigenetics data of arabidopsis
# H3K4me3 #
H3K4me3.data <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/epigenetics/H3K4me3_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
raw.dataset.values <- H3K4me3.data[, grep("LL|LD|ZT|CT", colnames(H3K4me3.data))]

nb.timepoints <- ncol(raw.dataset.values)
m.mean = apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
M.mean = mean(m.mean, na.rm = TRUE)
sd.mean = sd(m.mean, na.rm = TRUE)
zscore.raw.dataset.values <- raw.dataset.values
for (i in 1:nrow(raw.dataset.values)) {
  for (j in 1:ncol(raw.dataset.values)) {
    zscore.raw.dataset.values[i, j] = (nb.timepoints*raw.dataset.values[i, j] - M.mean)/sd.mean
  }
}
zscore.raw.dataset.values <- zscore.raw.dataset.values + abs(min(zscore.raw.dataset.values, na.rm = T))
#zscore.raw.dataset.values$mean <- apply(zscore.raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
#raw.dataset.values$mean <- apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
#plot(zscore.raw.dataset.values$mean, raw.dataset.values$mean)

#zscore.raw.dataset.values <- zscore.raw.dataset.values / max(zscore.raw.dataset.values, na.rm = TRUE)
zscore.raw.dataset.values <- cbind(as.data.frame(H3K4me3.data[, -grep("LL|LD|ZT|CT", colnames(H3K4me3.data))]),
                                   zscore.raw.dataset.values)
colnames(zscore.raw.dataset.values)[-grep("LL|LD|ZT|CT", colnames(H3K4me3.data))] <- colnames(H3K4me3.data)[-grep("LL|LD|ZT|CT", colnames(H3K4me3.data))]

write.table(zscore.raw.dataset.values, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/epigenetics/H3K4me3_zscore_level.txt", sep = "\t", quote = F, row.names = F)

# H3K9ac #
H3K9ac.data <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/epigenetics/H3K9ac_data.txt", sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
raw.dataset.values <- H3K9ac.data[, grep("LL|LD|ZT|CT", colnames(H3K9ac.data))]

nb.timepoints <- ncol(raw.dataset.values)
m.mean = apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
M.mean = mean(m.mean, na.rm = TRUE)
sd.mean = sd(m.mean, na.rm = TRUE)
zscore.raw.dataset.values <- raw.dataset.values
for (i in 1:nrow(raw.dataset.values)) {
  for (j in 1:ncol(raw.dataset.values)) {
    zscore.raw.dataset.values[i, j] = (nb.timepoints*raw.dataset.values[i, j] - M.mean)/sd.mean
  }
}
zscore.raw.dataset.values <- zscore.raw.dataset.values + abs(min(zscore.raw.dataset.values, na.rm = T))
#zscore.raw.dataset.values$mean <- apply(zscore.raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
#raw.dataset.values$mean <- apply(raw.dataset.values, 1, FUN = function(x) {mean(x, na.rm = TRUE)})
#plot(zscore.raw.dataset.values$mean, raw.dataset.values$mean)
#zscore.raw.dataset.values <- zscore.raw.dataset.values / max(zscore.raw.dataset.values, na.rm = TRUE)
zscore.raw.dataset.values <- cbind(as.data.frame(H3K9ac.data[, -grep("LL|LD|ZT|CT", colnames(H3K9ac.data))]),
                                   zscore.raw.dataset.values)
colnames(zscore.raw.dataset.values)[-grep("LL|LD|ZT|CT", colnames(H3K9ac.data))] <- colnames(H3K9ac.data)[-grep("LL|LD|ZT|CT", colnames(H3K9ac.data))]

write.table(zscore.raw.dataset.values, "~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/epigenetics/H3K9ac_zscore_level.txt", sep = "\t", quote = F, row.names = F)



