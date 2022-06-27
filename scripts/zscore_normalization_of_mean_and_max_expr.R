##################################
# Z-score normalization of data #
##################################
# Normalize mean and max expression values (over time-points) to obtain a null distribution
# With z-score
# Such as: z-score of gene i = Zi = (mi-m)/sd
# with:
# mean expression of full dataset = m
# standard deviation of full dataset = sd
# mean or max expression values of time-points for the gene i = mi
########################

species.list <- c("arabidopsis", "cyanobacteria", "ostreococcus", "mouse")
tissue.list <- c("leaves", NA, NA, "liver")

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species.list, tissue.list, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

for (y in 1:length(species.list)) {
  proteome.data <- read.table(paste(proteome.dir[y], "proteome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE)
  
  raw.dataset.values <- proteome.data[, grep("LL|LD|ZT|CT", colnames(proteome.data))]
  raw.dataset.values$mean.expr.value = apply(raw.dataset.values, 1, FUN = function(x) {return(mean(x, na.rm=TRUE))})
  raw.dataset.values$max.expr.value = apply(raw.dataset.values, 1, FUN = function(x){mean(max(x, na.rm = T), max(x[x!=max(x, na.rm = T)], na.rm = T))})
  
  m.mean = mean(raw.dataset.values$mean.expr.value, na.rm = TRUE)
  m.max = mean(raw.dataset.values$max.expr.value, na.rm = TRUE)
  sd.mean = sd(raw.dataset.values$mean.expr.value, na.rm = TRUE)
  sd.max = sd(raw.dataset.values$max.expr.value, na.rm = TRUE)
  
  proteome.data$zscore.mean.expr.value <- (raw.dataset.values$mean.expr.value - m.mean) / sd.mean
  proteome.data$zscore.max.expr.value <- (raw.dataset.values$max.expr.value - m.max) / sd.max
  # to get positive values to then calculate a coherent total expression cost
  proteome.data$zscore.mean.expr.value <- proteome.data$zscore.mean.expr.value + abs(min(proteome.data$zscore.mean.expr.value, na.rm = TRUE))
  proteome.data$zscore.max.expr.value <- proteome.data$zscore.max.expr.value + abs(min(proteome.data$zscore.max.expr.value, na.rm = TRUE))
  
  #hist(log(raw.dataset.values$mean.expr.value), breaks = 150)
  #hist(log(raw.dataset.values$max.expr.value), breaks = 150)
  #hist(log(raw.dataset.values$zscore.mean.expr.value), breaks = 150)
  #hist(log(raw.dataset.values$zscore.max.expr.value), breaks = 150)
  
  write.table(proteome.data, paste(proteome.dir[y], "proteome_data.txt", sep="/"), sep = "\t", quote = F, row.names = F)
}

