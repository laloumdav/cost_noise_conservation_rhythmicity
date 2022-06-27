library(ggplot2)
library(ggridges)
library(officer)
library(magrittr)
library(ggpubr)
library(latex2exp)

species.list <- c("arabidopsis", "cyanobacteria", "ostreococcus", "mouse")
tissue.list <- c("leaves", NA, NA, "liver")

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species.list, tissue.list, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

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
                       rhythmic=NA)
for (i in 1:length(species.list)) {
  transcriptome.data[[i]] <- read.table(paste(transcriptome.dir[i], "transcriptome_data.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "proteome_data.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  
  tot.data.transcripts <- data.frame(max.expr.value = transcriptome.data[[i]]$RNA_max.level, 
                                     mean.expr.value = transcriptome.data[[i]]$RNA_mean.level, 
                                     rhythm.pvalue = transcriptome.data[[i]]$RNA_rhythm.pvalue,
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
  #random.true <- sample(1:nrow(tot.data.transcripts), size=rhythmic.group.size.transcripts, replace = FALSE)
  #tot.data.transcripts[random.true, "random.group"] <- TRUE
  #random.true <- sample(1:nrow(tot.data.proteins), size=rhythmic.group.size.proteins, replace = FALSE)
  #tot.data.proteins[random.true, "random.group"] <- TRUE
  
  tot.data.tmp <- rbind(tot.data.transcripts, tot.data.proteins)
  tot.data <- rbind(tot.data, tot.data.tmp)
}
tot.data <- tot.data[-1, ]


### PLOTS ### 
subset.rhythmic.data <- subset(tot.data, rhythmic == TRUE)
subset.rhythmic.data$median.mean.expr.value <- NA
subset.rhythmic.data$median.max.expr.value <- NA
subset.nonrhythmic.data <- subset(tot.data, rhythmic == FALSE)
for (i in 1:length(species.list)) {
  subset.data.transcripts.tmp <- subset(subset.nonrhythmic.data, species == species.list[i] & data == "transcripts")
  subset.data.proteins.tmp <- subset(subset.nonrhythmic.data, species == species.list[i] & data == "proteins")
  
  subset.nonrhythmic.data[subset.nonrhythmic.data$species == species.list[i] & subset.nonrhythmic.data$data == "transcripts", "median.mean.expr.value"] <- median(subset.data.transcripts.tmp$mean.expr.value, na.rm = TRUE)
  subset.nonrhythmic.data[subset.nonrhythmic.data$species == species.list[i] & subset.nonrhythmic.data$data == "proteins", "median.mean.expr.value"] <- median(subset.data.proteins.tmp$mean.expr.value, na.rm = TRUE)

  subset.nonrhythmic.data[subset.nonrhythmic.data$species == species.list[i] & subset.nonrhythmic.data$data == "transcripts", "median.max.expr.value"] <- median(subset.data.transcripts.tmp$max.expr.value, na.rm = TRUE)
  subset.nonrhythmic.data[subset.nonrhythmic.data$species == species.list[i] & subset.nonrhythmic.data$data == "proteins", "median.max.expr.value"] <- median(subset.data.proteins.tmp$max.expr.value, na.rm = TRUE)
}
subset.tot.data <- rbind(subset.nonrhythmic.data, subset.rhythmic.data)
subset.tot.data$genes.gp <- c( rep("non-rhythmic", nrow(subset.nonrhythmic.data)), rep("rhythmic", nrow(subset.rhythmic.data)))


### 1. max ###
##############
maxExprHistogramsPerGroup <- ggplot(subset.tot.data, aes(x=log10(max.expr.value), fill=genes.gp)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.6, na.rm = TRUE, bins = 40) + 
  scale_fill_manual(values = c("seagreen4", "#E7B800")) + 
  scale_color_manual(values = c("seagreen4", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "max expression (log)", fill="genes group", caption=paste("proportion of rhythmic genes group = ", percent.rhythmic*100, "%", sep="")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0))

### 2. mean ###
##############
meanExprHistogramsPerGroup <- ggplot(subset.tot.data, aes(x=log10(mean.expr.value), fill=genes.gp)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.6, na.rm = TRUE, bins = 40) + 
  scale_fill_manual(values = c("seagreen4", "#E7B800")) + 
  scale_color_manual(values = c("seagreen4", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "mean expression (log)", fill="genes group", caption=paste("proportion of rhythmic genes group = ", percent.rhythmic*100, "%", sep="")) +
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
caption.text <- paste("Box-plot log10 scaled",
                      paste("rhythmic genes: first ", percent.rhythmic*100, "%", sep=""), sep="\n")

### 1. max ###
##############
sub.table <- subset(subset.tot.data, data == "proteins")
maxMean_ProteinExpr_ComparedPerGroup <- ggplot(sub.table, aes(species, log10(max.expr.value), fill = genes.gp)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = 1) +
  facet_wrap(~ species, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=2.4) +
  geom_hline(aes(yintercept = log10(median.max.expr.value)), linetype = 2, color="grey39", size=.2) +
  scale_fill_manual(values=c("seagreen4", "#E7B800")) +
  labs(y = TeX("$\\textbf{_{max}N_{p}}$"),
       fill="genes group:",
       caption=caption.text) +
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=6), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_text(size=11), 
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(.07, "lines"),
    legend.key = element_blank(),
    #axis.line.x = element_line(colour = "grey"),
    #axis.line.y = element_line(colour = "grey"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "grey32"),
    strip.text = element_blank(),
    legend.position = "top",
    plot.caption = element_text(size=7, face = "italic", hjust=0))

sub.table <- subset(subset.tot.data, data == "transcripts")
maxMean_RNAexpr_ComparedPerGroup <- ggplot(sub.table, aes(species, log10(max.expr.value), fill = genes.gp)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = 1) +
  facet_wrap(~ species, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=2.4) +
  geom_hline(aes(yintercept = log10(median.max.expr.value)), linetype = 2, color="grey39", size=.2) +
  scale_fill_manual(values=c("seagreen4", "#E7B800")) +
  labs(y = TeX("$\\textbf{_{max}N_{RNA}}$"),
       fill="genes group:",
       caption=caption.text) +
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=6), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_text(size=11), 
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(.07, "lines"),
    legend.key = element_blank(),
    #axis.line.x = element_line(colour = "grey"),
    #axis.line.y = element_line(colour = "grey"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "grey32"),
    strip.text = element_blank(),
    legend.position = "top",
    plot.caption = element_text(size=7, face = "italic", hjust=0))

### 2. mean ###
##############
sub.table <- subset(subset.tot.data, data == "proteins")
meanMean_ProteinExpr_ComparedPerGroup <- ggplot(sub.table, aes(species, log10(mean.expr.value), fill = genes.gp)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = 1) +
  facet_wrap(~ species, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=2.4) +
  geom_hline(aes(yintercept = log10(median.mean.expr.value)), linetype = 2, color="grey39", size=.2) +
  scale_fill_manual(values=c("seagreen4", "#E7B800")) +
  labs(y = TeX("$\\textbf{_{mean}N_{p}}$"),
       fill="genes group:",
       caption=caption.text) +
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=6), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_text(size=11), 
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(.07, "lines"),
    legend.key = element_blank(),
    #axis.line.x = element_line(colour = "grey"),
    #axis.line.y = element_line(colour = "grey"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "grey32"),
    strip.text = element_blank(),
    legend.position = "top",
    plot.caption = element_text(size=7, face = "italic", hjust=0))

sub.table <- subset(subset.tot.data, data == "transcripts")
meanMean_RNAexpr_ComparedPerGroup <- ggplot(sub.table, aes(species, log10(mean.expr.value), fill = genes.gp)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = 1) +
  facet_wrap(~ species, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=2.4) +
  geom_hline(aes(yintercept = log10(median.mean.expr.value)), linetype = 2, color="grey39", size=.2) +
  scale_fill_manual(values=c("seagreen4", "#E7B800")) +
  labs(y = TeX("$\\textbf{_{mean}N_{RNA}}$"),
       fill="genes group:",
       caption=caption.text) +
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=6), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_text(size=11), 
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(.07, "lines"),
    legend.key = element_blank(),
    #axis.line.x = element_line(colour = "grey"),
    #axis.line.y = element_line(colour = "grey"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "grey32"),
    strip.text = element_blank(),
    legend.position = "top",
    plot.caption = element_text(size=7, face = "italic", hjust=0))

##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##max_expr_compared_between_groups##") %>%
  body_add_gg(value = maxMean_RNAexpr_ComparedPerGroup, width = 5.5, height = 8.5) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##mean_expr_compared_between_groups##") %>%
  body_add_gg(value = meanMean_RNAexpr_ComparedPerGroup, width = 5.5, height = 8.5) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")



out.rds.name <- paste("~/Documents/cost_theory_workingspace/rds_objects", "proteins_maxMeanExprComparedPerGroup.rds", sep="/")
saveRDS(maxMean_ProteinExpr_ComparedPerGroup, file = out.rds.name)
out.rds.name <- paste("~/Documents/cost_theory_workingspace/rds_objects", "proteins_meanMeanExprComparedPerGroup.rds", sep="/")
saveRDS(meanMean_ProteinExpr_ComparedPerGroup, file = out.rds.name)

pdf("~/Documents/cost_theory_workingspace/RNA_meanMeanExprComparedPerGroup.pdf", width = 7, height = 4) # or other device
print(meanMean_RNAexpr_ComparedPerGroup)
dev.off()

pdf("~/Documents/cost_theory_workingspace/RNA_maxMeanExprComparedPerGroup.pdf", width = 7, height = 4) # or other device
print(maxMean_RNAexpr_ComparedPerGroup)
dev.off()
