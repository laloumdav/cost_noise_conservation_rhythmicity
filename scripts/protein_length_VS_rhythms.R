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
tot.data <- data.frame(average.AA.synth.cost=NA, 
                       protein.length=NA,
                       rhythm.pvalue=NA,
                       nb.genes=NA,
                       data=NA, 
                       species=NA,
                       rhythmic.order=NA,
                       rhythmic=NA)
for (i in 1:length(species.list)) {
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "proteome_data.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  
  tot.data.proteins <- data.frame(average.AA.synth.cost = proteome.data[[i]]$aa.synthesis.average.cost_Akashi, 
                                  protein.length = proteome.data[[i]]$protein.length, 
                                  rhythm.pvalue = proteome.data[[i]]$protein_rhythm.pvalue,
                                  nb.genes = nrow(proteome.data[[i]]),
                                  data = rep("proteins", nrow(proteome.data[[i]])), 
                                  species = rep(species.list[i], nrow(proteome.data[[i]])))
  
  # order genes by their rhythm p-value
  row.nb <- order(tot.data.proteins$rhythm.pvalue)
  tot.data.proteins[row.nb, "rhythmic.order"] <- 1:nrow(tot.data.proteins)
  # Lets consider the first n percent (percent.rhythmic) as rhythmic genes
  tot.data.proteins$rhythmic <- (tot.data.proteins$rhythmic.order/nrow(tot.data.proteins) <= percent.rhythmic)
  rhythmic.group.size = length(which(tot.data.proteins$rhythmic == TRUE))
  # Lets take a random group of genes of size equal to rhythmic genes group size and affect TRUE randomly (size of TRUE genes = rhythmic.group.size) 
  #random.true <- sample(1:nrow(tot.data.proteins), size=rhythmic.group.size, replace = FALSE)
  #tot.data.proteins[random.true, "random.group"] <- TRUE
  
  # re-order columns 
  tot.data.proteins <- tot.data.proteins[, colnames(tot.data)]
  
  tot.data <- rbind(tot.data, tot.data.proteins)
}
tot.data <- tot.data[-1,]


### PLOTS ### 
subset.rhythmic.data <- subset(tot.data, rhythmic == TRUE)
subset.rhythmic.data$median.protein.length <- NA
subset.nonrhythmic.data <- subset(tot.data, rhythmic == FALSE)
for (i in 1:length(species.list)) {
  subset.data.tmp <- subset(subset.nonrhythmic.data, species == species.list[i])
  subset.nonrhythmic.data[subset.nonrhythmic.data$species == species.list[i], "median.protein.length"] <- median(subset.data.tmp$protein.length, na.rm = TRUE)
}
subset.tot.data <- rbind(subset.nonrhythmic.data, subset.rhythmic.data)
subset.tot.data$genes.gp <- c( rep("non-rhythmic", nrow(subset.nonrhythmic.data)), rep("rhythmic", nrow(subset.rhythmic.data)))


### 2. protein length ###
##############

proteinLengthHistPerGroup <- ggplot(subset.tot.data, aes(x=log2(protein.length), fill=genes.gp)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.6, na.rm = TRUE, bins = 40) + 
  scale_fill_manual(values = c("seagreen4", "#E7B800")) + 
  scale_color_manual(values = c("seagreen4", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  labs(x = "protein length (log)", fill="genes group", caption=paste("proportion of rhythmic or random genes group = ", percent.rhythmic*100, "%", sep="")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0),
        strip.text = element_text(size=14))

proteinLengthQQplot <- ggplot(subset.tot.data, aes(sample=log2(protein.length))) + 
  stat_qq(aes(color=genes.gp), size=0.1) + 
  stat_qq_line(aes(color=genes.gp)) +
  scale_color_manual("genes group", values = c("seagreen4", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  labs(caption=paste("proportion of rhythmic or random genes group = ", percent.rhythmic*100, "%", sep="")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0),
        strip.text = element_text(size=14))

# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/supplementary_information/proteinLength_hist_perGroup.pdf", height = 3, width = 9)
print(proteinLengthHistPerGroup)
dev.off()
pdf("~/Documents/cost_theory_workingspace/supplementary_information/proteinLength_QQplot.pdf", height = 3, width = 9)
print(proteinLengthQQplot)
dev.off()

##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##protein_length_histograms_per_group##") %>%
  body_add_gg(value = proteinLengthHistPerGroup, width = 7, height = 2) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##protein_length_QQplots_per_group##") %>%
  body_add_gg(value = proteinLengthQQplot, width = 7, height = 2) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")


####################################################################################
### Comparison of protein lengths between both groups of genes ###########
####################################################################################
caption.text <- paste("Box-plot log2 scaled",
                      paste("rhythmic genes: first ", percent.rhythmic*100, "%", sep=""), sep="\n")

proteinLengthComparedPerGroup <- ggplot(subset.tot.data, aes(species, log2(protein.length), fill = genes.gp)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = .9) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=2.4) +
  scale_fill_manual(values=c("seagreen4", "#E7B800")) +
  labs(y = TeX("$\\textbf{L_{p}}$"), 
       fill="genes group:", 
       caption=caption.text) +
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=9), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_text(size=11), 
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    legend.key = element_blank(),
    axis.line = element_line(colour = "grey32"),
    plot.caption = element_text(size=7, face = "italic", hjust=0))

out.rds.name <- paste("~/Documents/cost_theory_workingspace/rds_objects", "proteinLengthComparedPerGroup.rds", sep="/")
saveRDS(proteinLengthComparedPerGroup, file = out.rds.name)

