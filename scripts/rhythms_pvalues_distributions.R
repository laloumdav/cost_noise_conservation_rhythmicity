############################################
###### rhythms pvalues distributions #######
############################################
library(ggplot2)
library(ggridges)
library(officer)
library(magrittr)

species <- c("arabidopsis", "cyanobacteria", "ostreococcus", "mouse")
tissue <- c("leaves", NA, NA, "liver")

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

transcriptome.data <- list()
proteome.data <- list()
tot.data <- data.frame(pvalue=NA, 
                       fdr=NA,
                       data=NA, 
                       species=NA)
for (i in 1:length(species)) {
  transcriptome.data[[i]] <- read.table(paste(transcriptome.dir[i], "transcriptome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE, sep="\t")
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "proteome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE, sep="\t")
  
  tot.data.transcripts <- data.frame(pvalue=transcriptome.data[[i]]$RNA_rhythm.pvalue, 
                                     fdr = p.adjust(transcriptome.data[[i]]$RNA_rhythm.pvalue, method = "fdr"),
                                     data=rep("transcripts", nrow(transcriptome.data[[i]])), 
                                     species=rep(species[i], nrow(transcriptome.data[[i]])))
  tot.data.proteins <- data.frame(pvalue=proteome.data[[i]]$protein_rhythm.pvalue, 
                                  fdr = p.adjust(proteome.data[[i]]$protein_rhythm.pvalue, method = "fdr"),
                                  data=rep("proteins", nrow(proteome.data[[i]])), 
                                  species=rep(species[i], nrow(proteome.data[[i]])))
  tot.data.tmp <- rbind(tot.data.transcripts, tot.data.proteins)
  tot.data <- rbind(tot.data, tot.data.tmp)
} 

tot.data <- tot.data[-1,]
  
# Special case for Ostreococcus transcriptome where there are several data for unique geneID => Brown normalization
transcriptome.ostreococcus <- read.table(paste(transcriptome.dir[grep("ostreococcus", species)], "transcriptome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE)
transcriptome.ostreococcus <- data.frame(pvalue=transcriptome.ostreococcus$RNA_rhythm.pvalue_Brown, 
                                data=rep("unique transcripts", nrow(transcriptome.ostreococcus)), 
                                species=rep("ostreococcus", nrow(transcriptome.ostreococcus)))

### PLOTS ###
pvalDistribPlot <- ggplot(tot.data, aes(x=pvalue)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.8, na.rm = TRUE, bins = 100) + 
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))

fdrDistribPlot <- ggplot(tot.data, aes(x=fdr)) + 
  geom_density() + 
  facet_wrap(species ~ data, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  labs(x = "rhythm detection p-value", y = "density") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))

pvalDistribBrownPlot <- ggplot(transcriptome.ostreococcus, aes(x=pvalue)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.8, na.rm = TRUE, bins = 100) + 
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  labs(x = "rhythm detection p-value", y = "density") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))

# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/supplementary_information/pval_distrib.pdf", width = 4)
print(pvalDistribPlot)
dev.off()
pdf("~/Documents/cost_theory_workingspace/supplementary_information/pval_distrib_brown_ostreococcus.pdf", width = 2.5, height = 2.5)
print(pvalDistribBrownPlot)
dev.off()

# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##pval_distrib##") %>%
  body_add_gg(value = pvalDistribPlot, width = 4) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##fdr_distrib##") %>%
  body_add_gg(value = fdrDistribPlot, width = 4) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##pval_distrib_brown##") %>%
  body_add_gg(value = pvalDistribBrownPlot, width = 2.4, height = 2.4) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")
 

