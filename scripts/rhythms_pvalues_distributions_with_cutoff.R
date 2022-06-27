############################################
###### rhythms pvalues distributions #######
############################################
library(ggplot2)
library(ggridges)
library(officer)
library(magrittr)
library(gridExtra)
library(grid)

species <- c("arabidopsis", "cyanobacteria", "ostreococcus", "mouse")
tissue <- c("leaves", NA, NA, "liver")

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

transcriptome.data <- list()
proteome.data <- list()
tot.data <- data.frame(pvalue=NA, 
                       fdr=NA,
                       nb.rhythmic=NA,
                       nb.genes=NA,
                       data=NA, 
                       species=NA)
for (i in 1:length(species)) {
  transcriptome.data[[i]] <- read.table(paste(transcriptome.dir[i], "transcriptome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE, sep="\t")
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "proteome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE, sep="\t")
  
  tot.data.transcripts <- data.frame(pvalue=transcriptome.data[[i]]$RNA_rhythm.pvalue, 
                                     fdr = p.adjust(transcriptome.data[[i]]$RNA_rhythm.pvalue, method = "fdr"),
                                     nb.rhythmic = length(which(transcriptome.data[[i]]$RNA_rhythm.pvalue <= 0.01)),
                                     nb.genes = nrow(transcriptome.data[[i]]),
                                     data=rep("transcripts", nrow(transcriptome.data[[i]])), 
                                     species=rep(species[i], nrow(transcriptome.data[[i]])))
  tot.data.proteins <- data.frame(pvalue=proteome.data[[i]]$protein_rhythm.pvalue, 
                                  fdr = p.adjust(proteome.data[[i]]$protein_rhythm.pvalue, method = "fdr"),
                                  nb.rhythmic = length(which(proteome.data[[i]]$protein_rhythm.pvalue <= 0.01)),
                                  nb.genes = nrow(proteome.data[[i]]),
                                  data=rep("proteins", nrow(proteome.data[[i]])), 
                                  species=rep(species[i], nrow(proteome.data[[i]])))
  tot.data.tmp <- rbind(tot.data.transcripts, tot.data.proteins)
  tot.data <- rbind(tot.data, tot.data.tmp)
} 

tot.data <- tot.data[-1,]

tot.data$percentage.rhythmic <- tot.data$nb.rhythmic / tot.data$nb.genes
tot.data$text.percentage.rhythmic <- paste(tot.data$data, "\n* ",round(tot.data$percentage.rhythmic*100, digits = 1), "% (", tot.data$nb.genes, ")", sep="")

### PLOTS ###
pvalDistribPlot <- ggplot(tot.data, aes(x=pvalue)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.8, na.rm = TRUE, bins = 100) + 
  facet_wrap(species ~ text.percentage.rhythmic, scales = "free", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "rhythm detection p-value (original paper)") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))

# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/supplementary_information/pval_distrib.pdf", width = 4, height = 8.5)
print(pvalDistribPlot)
dev.off()

# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##pval_distrib##") %>%
  body_add_gg(value = pvalDistribPlot, width = 4, height = 8) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##pval_distrib_text##") %>%
  body_add_par("*percentage of rhythmic genes with p-value < 0.01 (among total of genes)") %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")


