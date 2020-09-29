library(ggplot2)
library(ggridges)
library(officer)
library(magrittr)
library(ggpubr)

species.list <- c("arabidopsis", "cyanobacteria", "ostreococcus", "mouse")
tissue.list <- c("leaves", NA, NA, "liver")

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species.list, tissue.list, sep="/"))

percent.rhythmic = 0.15


##############
transcriptome.data <- list()
proteome.data <- list()
tot.data <- data.frame(max.expr.value=NA, 
                       mean.expr.value=NA, 
                       rhythm.pvalue=NA,
                       nb.genes=NA,
                       data=NA, 
                       species=NA,
                       rhythmic.order=NA,
                       rhythmic=NA,
                       random.group=NA)
for (i in 1:length(species.list)) {
  tot.data.tmp <- read.table(paste(main.dir[i], "tot_data.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  
  transcriptome.data[[i]] <- tot.data.tmp[, grep("RNA", colnames(tot.data.tmp))]
  transcriptome.data[[i]] <- unique(subset(transcriptome.data[[i]], !is.na(RNA_mean.level)))
  proteome.data[[i]] <- tot.data.tmp[, grep("protein", colnames(tot.data.tmp))]
  proteome.data[[i]] <- unique(subset(proteome.data[[i]], !is.na(protein_mean.level)))
  
  tot.data.transcripts <- data.frame(max.expr.value = transcriptome.data[[i]]$RNA_max.level, 
                                     mean.expr.value = transcriptome.data[[i]]$RNA_mean.level, 
                                     rhythm.pvalue = transcriptome.data[[i]][, grep("^RNA_rhythm.pvalue$|^RNA_rhythm.pvalue_Brown$", colnames(transcriptome.data[[i]]))],
                                     nb.genes = nrow(transcriptome.data[[i]]),
                                     data = rep("transcripts", nrow(transcriptome.data[[i]])), 
                                     species = rep(species.list[i], nrow(transcriptome.data[[i]])))
  tot.data.proteins <- data.frame(max.expr.value = proteome.data[[i]]$protein_max.level, 
                                  mean.expr.value = proteome.data[[i]]$protein_mean.level, 
                                  rhythm.pvalue = proteome.data[[i]]$protein_rhythm.pvalue,
                                  nb.genes = nrow(proteome.data[[i]]),
                                  data = rep("proteins", nrow(proteome.data[[i]])), 
                                  species = rep(species.list[i], nrow(proteome.data[[i]])))
  
  # order genes by their rhythm p-value
  row.nb <- order(tot.data.transcripts$rhythm.pvalue)
  tot.data.transcripts[row.nb, "rhythmic.order"] <- 1:nrow(tot.data.transcripts)
  row.nb <- order(tot.data.proteins$rhythm.pvalue)
  tot.data.proteins[row.nb, "rhythmic.order"] <- 1:nrow(tot.data.proteins)
  # Lets consider the first n percent (percent.rhythmic) as rhythmic genes
  tot.data.transcripts$rhythmic <- (tot.data.transcripts$rhythmic.order/nrow(tot.data.transcripts) <= percent.rhythmic)
  rhythmic.group.size.transcripts = length(which(tot.data.transcripts$rhythmic == TRUE))
  tot.data.proteins$rhythmic <- (tot.data.proteins$rhythmic.order/nrow(tot.data.proteins) <= percent.rhythmic)
  rhythmic.group.size.proteins = length(which(tot.data.proteins$rhythmic == TRUE))
  
  # Lets take a random group of genes of size equal to rhythmic genes group size and affect TRUE randomly (size of TRUE genes = rhythmic.group.size) 
  random.true <- sample(1:nrow(tot.data.transcripts), size=rhythmic.group.size.transcripts, replace = FALSE)
  tot.data.transcripts[random.true, "random.group"] <- TRUE
  random.true <- sample(1:nrow(tot.data.proteins), size=rhythmic.group.size.proteins, replace = FALSE)
  tot.data.proteins[random.true, "random.group"] <- TRUE
  
  
  tot.data.tmp <- rbind(tot.data.transcripts, tot.data.proteins)
  tot.data <- rbind(tot.data, tot.data.tmp)
}
tot.data <- tot.data[-1, ]




proteome.data[[i]] <- read.table(paste(main.dir[i], "tot_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE)
tot.data.proteins <- data.frame(max.expr.value = proteome.data[[i]]$protein_max.level, 
                                mean.expr.value = proteome.data[[i]]$protein_mean.level, 
                                aa.synthesis.average.cost = proteome.data[[i]]$aa.synthesis.average.cost,
                                protein_rhythm.pvalue = proteome.data[[i]]$protein_rhythm.pvalue,
                                RNA_rhythm.pvalue = proteome.data[[i]]$RNA_rhythm.pvalue,
                                nb.genes = nrow(proteome.data[[i]]),
                                data = rep("proteins", nrow(proteome.data[[i]])), 
                                species = rep(species.list[i], nrow(proteome.data[[i]])))

row.nb <- order(tot.data.proteins$rhythm.pvalue)
tot.data.proteins[row.nb, "rhythmic.order"] <- 1:nrow(tot.data.proteins)
# Lets consider the first n percent (percent.rhythmic) as rhythmic genes
tot.data.proteins$rhythmic <- (tot.data.proteins$rhythmic.order/nrow(tot.data.proteins) <= percent.rhythmic)
rhythmic.group.size.proteins = length(which(tot.data.proteins$rhythmic == TRUE))
# Lets take a random group of genes of size equal to rhythmic genes group size and affect TRUE randomly (size of TRUE genes = rhythmic.group.size) 
random.true <- sample(1:nrow(tot.data.proteins), size=rhythmic.group.size.proteins, replace = FALSE)
tot.data.proteins[random.true, "random.group"] <- TRUE

tot.data.proteins <- subset(tot.data.proteins, (!is.na(tot.data.proteins[, "aa.synthesis.average.cost"])))
t.test(tot.data.proteins[tot.data.proteins$protein_rhythm.pvalue<=0.02 &
                           tot.data.proteins$RNA_rhythm.pvalue>0.1, "mean.expr.value"], 
       tot.data.proteins[tot.data.proteins$protein_rhythm.pvalue>0.1 |
                           tot.data.proteins$RNA_rhythm.pvalue<=0.02, "mean.expr.value"])
wilcox.test(tot.data.proteins[tot.data.proteins$protein_rhythm.pvalue<=0.02 &
                                tot.data.proteins$RNA_rhythm.pvalue>0.1, "mean.expr.value"], 
            tot.data.proteins[tot.data.proteins$protein_rhythm.pvalue>0.1 |
                                tot.data.proteins$RNA_rhythm.pvalue<=0.02, "mean.expr.value"])

nrRNA-rProt
nrRNA-nrProt

### PLOTS ### 
subset.rhythmic.data <- subset(tot.data, rhythmic == TRUE)
subset.rhythmic.data$median.mean.expr.value <- NA
subset.rhythmic.data$median.max.expr.value <- NA
subset.random.data <- subset(tot.data, random.group == TRUE)
for (i in 1:length(species.list)) {
  subset.data.transcripts.tmp <- subset(subset.random.data, species == species.list[i] & data == "transcripts")
  subset.data.proteins.tmp <- subset(subset.random.data, species == species.list[i] & data == "proteins")
  
  subset.random.data[subset.random.data$species == species.list[i] & subset.random.data$data == "transcripts", "median.mean.expr.value"] <- median(subset.data.transcripts.tmp$mean.expr.value, na.rm = TRUE)
  subset.random.data[subset.random.data$species == species.list[i] & subset.random.data$data == "proteins", "median.mean.expr.value"] <- median(subset.data.proteins.tmp$mean.expr.value, na.rm = TRUE)
  
  subset.random.data[subset.random.data$species == species.list[i] & subset.random.data$data == "transcripts", "median.max.expr.value"] <- median(subset.data.transcripts.tmp$max.expr.value, na.rm = TRUE)
  subset.random.data[subset.random.data$species == species.list[i] & subset.random.data$data == "proteins", "median.max.expr.value"] <- median(subset.data.proteins.tmp$max.expr.value, na.rm = TRUE)
}
subset.tot.data <- rbind(subset.random.data, subset.rhythmic.data)
subset.tot.data$rhythmic.vs.random <- c( rep("random", nrow(subset.random.data)), rep("rhythmic", nrow(subset.rhythmic.data)))


### 1. max ###
##############
maxExprHistogramsPerGroup <- ggplot(subset.tot.data, aes(x=log10(max.expr.value), fill=rhythmic.vs.random)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.6, na.rm = TRUE, bins = 40) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "max expression (log)", fill="genes group", caption=paste("proportion of rhythmic or random genes group = ", percent.rhythmic*100, "%", sep="")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0))

### 2. mean ###
##############
meanExprHistogramsPerGroup <- ggplot(subset.tot.data, aes(x=log10(mean.expr.value), fill=rhythmic.vs.random)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.6, na.rm = TRUE, bins = 40) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "mean expression (log)", fill="genes group", caption=paste("proportion of rhythmic or random genes group = ", percent.rhythmic*100, "%", sep="")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0))

##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##max_expr_histograms_per_group##") %>%
  body_add_gg(value = maxExprHistogramsPerGroup, width = 5, height = 7) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##mean_expr_histograms_per_group##") %>%
  body_add_gg(value = meanExprHistogramsPerGroup, width = 5, height = 7) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")








####################################################################################
### Comparison of mean values of expression between both groups of genes ###########
####################################################################################
alternative.wilcox.option <- "two.sided"  # choice between "two.sided", "less", "greater"

caption.text <- paste("Box-plot log scaled",
                      #"\n",
                      paste("Wilcox test (alternative='", alternative.wilcox.option, "'):", sep=""),
                      "ns: p > 0.05",
                      "*: p <= 0.05",
                      "**: p <= 0.01",
                      "***: p <= 0.001",
                      "****: p <= 0.0001",
                      #"\n",
                      paste("proportion of rhythmic or random genes group = ", percent.rhythmic*100, "%", sep=""),
                      sep="\n")


### PLOTS ### 
subset.rhythmic.data <- subset(tot.data, rhythmic == TRUE)
subset.rhythmic.data$median.mean.expr.value <- NA
subset.rhythmic.data$median.max.expr.value <- NA
subset.random.data <- subset(tot.data, random.group == TRUE)
for (i in 1:length(species.list)) {
  subset.data.transcripts.tmp <- subset(subset.random.data, species == species.list[i] & data == "transcripts")
  subset.data.proteins.tmp <- subset(subset.random.data, species == species.list[i] & data == "proteins")
  
  subset.random.data[subset.random.data$species == species.list[i] & subset.random.data$data == "transcripts", "median.mean.expr.value"] <- median(subset.data.transcripts.tmp$mean.expr.value, na.rm = TRUE)
  subset.random.data[subset.random.data$species == species.list[i] & subset.random.data$data == "proteins", "median.mean.expr.value"] <- median(subset.data.proteins.tmp$mean.expr.value, na.rm = TRUE)
  
  subset.random.data[subset.random.data$species == species.list[i] & subset.random.data$data == "transcripts", "median.max.expr.value"] <- median(subset.data.transcripts.tmp$max.expr.value, na.rm = TRUE)
  subset.random.data[subset.random.data$species == species.list[i] & subset.random.data$data == "proteins", "median.max.expr.value"] <- median(subset.data.proteins.tmp$max.expr.value, na.rm = TRUE)
}
subset.tot.data <- rbind(subset.random.data, subset.rhythmic.data)
subset.tot.data$rhythmic.vs.random <- c( rep("random", nrow(subset.random.data)), rep("rhythmic", nrow(subset.rhythmic.data)))



### 1. max ###
##############
maxMeanExprComparedPerGroup <- ggboxplot(subset.tot.data, x = "rhythmic.vs.random", 
                                         y = "max.expr.value", yscale = "log10",
                                         color = "rhythmic.vs.random", palette = c("#00AFBB", "#E7B800"),
                                         width = 0.4,
                                         #add = "jitter", 
                                         short.panel.labs = T) + 
  geom_hline(aes(yintercept = median.max.expr.value), linetype = 2, color="grey39") +
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = alternative.wilcox.option), 
                     label="p.signif" , label.x.npc = 0.4, label.y.npc = 0.9, na.rm = TRUE) +
  labs(y = "max expression level", color="genes group:", caption=caption.text) +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=6), 
        axis.text.x = element_text(size=10), 
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0),
        legend.position = "right",
        strip.background = element_rect(fill="gray44", colour="black", size = 0.1), 
        strip.text = element_text(size=10, colour="white"))

### 2. mean ###
##############
meanMeanExprComparedPerGroup <- ggboxplot(subset.tot.data, x = "rhythmic.vs.random", 
                                          y = "mean.expr.value", yscale = "log10",
                                          color = "rhythmic.vs.random", palette = c("#00AFBB", "#E7B800"),
                                          width = 0.4,
                                          #add = "jitter", 
                                          short.panel.labs = T) + 
  geom_hline(aes(yintercept = median.mean.expr.value), linetype = 2, color="grey39") +
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = alternative.wilcox.option), 
                     label="p.signif" , label.x.npc = 0.4, label.y.npc = 0.9, na.rm = TRUE) +
  labs(y = "mean expression level", color="genes group:", caption=caption.text) +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=6), 
        axis.text.x = element_text(size=10), 
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0),
        legend.position = "right",
        strip.background = element_rect(fill="gray44", colour="black", size = 0.1), 
        strip.text = element_text(size=10, colour="white"))


##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##max_expr_compared_between_groups##") %>%
  body_add_gg(value = maxMeanExprComparedPerGroup, width = 5.5, height = 8.5) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##mean_expr_compared_between_groups##") %>%
  body_add_gg(value = meanMeanExprComparedPerGroup, width = 5.5, height = 8.5) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")


#rhythm.data <- subset(tot.data, rhythmic == TRUE)
#random.data <- subset(tot.data, random.group == TRUE)
#random.data <- random.data[sample(1:nrow(random.data), replace = F, size = nrow(rhythm.data)), ]
#shapiro.test(rhythm.data$mean.expr.value[sample(1:length(rhythm.data$mean.expr.value), replace = FALSE, size = 5000)])


##############
# Also save it as an R object:
caption.text <- paste("Box-plot log scaled",
                      paste("Wilcox test (alternative='", alternative.wilcox.option, "'):", sep=""),
                      sep="\n")

### 1. max ###
##############
maxMeanExprComparedPerGroup <- ggboxplot(subset.tot.data, x = "rhythmic.vs.random", 
                                         y = "max.expr.value", yscale = "log10", outlier.size=0.1,
                                         fill = "rhythmic.vs.random", palette = c("seagreen4", "#E7B800"),
                                         width = 0.4, alpha=0.8, size=0.2,
                                         #add = "jitter", 
                                         short.panel.labs = T) + 
  geom_hline(aes(yintercept = median.max.expr.value), linetype = 2, color="grey39", size=0.2) +
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = alternative.wilcox.option), 
                     label="p.signif" , label.x.npc = 0.4, label.y.npc = 0.9, na.rm = TRUE, size=2.2) +
  labs(y = "max expression level", fill="genes group:" #, caption=caption.text
  ) +
  theme(axis.text.y = element_text(size=5), 
        axis.title.y = element_text(size=7), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.caption = element_text(size=5, face = "italic", hjust=0),
        legend.position = "right",
        legend.text = element_text(size=5), 
        legend.title = element_text(size=5), 
        strip.background = element_blank(),
        panel.background = element_rect(fill = "grey93"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size=6, face = "bold", margin = margin(0)),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(size = 0.1))

### 2. mean ###
##############
meanMeanExprComparedPerGroup <- ggboxplot(subset.tot.data, x = "rhythmic.vs.random", 
                                          y = "mean.expr.value", yscale = "log10", outlier.size=0.1,
                                          fill = "rhythmic.vs.random", palette = c("seagreen4", "#E7B800"),
                                          width = 0.4, alpha=0.8, size=0.2,
                                          #add = "jitter", 
                                          short.panel.labs = T) + 
  geom_hline(aes(yintercept = median.mean.expr.value), linetype = 2, color="grey39") +
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = alternative.wilcox.option), 
                     label="p.signif" , label.x.npc = 0.4, label.y.npc = 0.9, na.rm = TRUE, size=2.2) +
  labs(y = "mean expression level", fill="genes group:" #, caption=caption.text
  ) +
  theme(axis.text.y = element_text(size=5), 
        axis.title.y = element_text(size=7), 
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.caption = element_text(size=5, face = "italic", hjust=0),
        legend.position = "right",
        legend.text = element_text(size=5), 
        legend.title = element_text(size=5), 
        strip.background = element_blank(),
        panel.background = element_rect(fill = "grey93"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size=6, face = "bold", margin = margin(0)),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(size = 0.1))


out.rds.name <- paste("~/Documents/cost_theory_workingspace/rds_objects", "maxMeanExprComparedPerGroup.rds", sep="/")
saveRDS(maxMeanExprComparedPerGroup, file = out.rds.name)
out.rds.name <- paste("~/Documents/cost_theory_workingspace/rds_objects", "meanMeanExprComparedPerGroup.rds", sep="/")
saveRDS(meanMeanExprComparedPerGroup, file = out.rds.name)



