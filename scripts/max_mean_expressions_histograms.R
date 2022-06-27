library(ggplot2)
library(ggridges)
library(officer)
library(magrittr)

species <- c("arabidopsis", "cyanobacteria", "ostreococcus", "mouse")
tissue <- c("leaves", NA, NA, "liver")

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")


##############
transcriptome.data <- list()
proteome.data <- list()
tot.data <- data.frame(max.expr.value=NA, 
                       mean.expr.value=NA, 
                       data=NA, 
                       species=NA)
for (i in 1:length(species)) {
  transcriptome.data[[i]] <- read.table(paste(transcriptome.dir[i], "transcriptome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE, sep="\t")
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "proteome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE, sep="\t")
  
  tot.data.transcripts <- data.frame(max.expr.value=transcriptome.data[[i]]$RNA_max.level, 
                                     mean.expr.value=transcriptome.data[[i]]$RNA_mean.level, 
                                     data=rep("transcripts", nrow(transcriptome.data[[i]])), 
                                     species=rep(species[i], nrow(transcriptome.data[[i]])))
  tot.data.proteins <- data.frame(max.expr.value=proteome.data[[i]]$protein_max.level, 
                                  mean.expr.value=proteome.data[[i]]$protein_mean.level, 
                                  data=rep("proteins", nrow(proteome.data[[i]])), 
                                  species=rep(species[i], nrow(proteome.data[[i]])))
  tot.data.tmp <- rbind(tot.data.transcripts, tot.data.proteins)
  tot.data <- rbind(tot.data, tot.data.tmp)
} 

tot.data <- tot.data[-1,]


### PLOTS ### 

### 1. max ###
maxExprHistograms <- ggplot(tot.data, aes(x=log10(max.expr.value))) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.8, na.rm = TRUE, bins = 60) + 
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "max expression (log)") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))

### 2. mean ###
meanExprHistograms <- ggplot(tot.data, aes(x=log10(mean.expr.value))) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.8, na.rm = TRUE, bins = 60) + 
  facet_wrap(species ~ data, scales = "free", ncol = 2) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "mean expression (log)") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))

# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/supplementary_information/max_expr.pdf", width = 4)
print(maxExprHistograms)
dev.off()
pdf("~/Documents/cost_theory_workingspace/supplementary_information/mean_expr.pdf", width = 4)
print(meanExprHistograms)
dev.off()

# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##max_expr_histograms##") %>%
  body_add_gg(value = maxExprHistograms, width = 4) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##mean_expr_histograms##") %>%
  body_add_gg(value = meanExprHistograms, width = 4) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")




