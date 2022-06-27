library(ggplot2)
library(ggridges)
library(officer)
library(magrittr)

Dir <- "~/Documents/cost_theory_workingspace"

species <- c("arabidopsis", "cyanobacteria", "ostreococcus", "mouse")
tissue <- c("leaves", NA, NA, "liver")

main.dir <- gsub("/NA", "", paste(Dir, "DATA", species, tissue, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

transcriptome.data <- list()
proteome.data <- list()
tot.data <- data.frame(rhythm.pvalue=NA, 
                       data=NA,
                       species=NA)
for (i in 1:length(species)) {
  transcriptome.data[[i]] <- read.table(paste(transcriptome.dir[i], "transcript_GeneCycle.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "protein_GeneCycle.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  
  tot.data.transcripts <- data.frame(rhythm.pvalue = transcriptome.data[[i]]$default.pvalue,
                                     data=rep(paste(tissue[i], "transcripts", sep="_"), nrow(transcriptome.data[[i]])), 
                                     species=rep(species[i], nrow(transcriptome.data[[i]])))
  tot.data.proteins <- data.frame(rhythm.pvalue = proteome.data[[i]]$default.pvalue,
                                  data=rep(paste(tissue[i], "proteins", sep="_"), nrow(proteome.data[[i]])), 
                                  species=rep(species[i], nrow(proteome.data[[i]])))
  
  tot.data.tmp <- rbind(tot.data.transcripts, tot.data.proteins)
  tot.data <- rbind(tot.data, tot.data.tmp)
}
tot.data <- tot.data[-1,]

# Add other tissues proteomics data from Mice
species <- c("mouse", "mouse", "mouse")
tissue <- c("cartilage", "tendon", "forebrain")
main.dir <- gsub("/NA", "", paste(Dir, "DATA", species, tissue, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")

proteome.data <- list()
for (i in 1:length(species)) {
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "protein_GeneCycle.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  
  tot.data.proteins <- data.frame(rhythm.pvalue = proteome.data[[i]]$default.pvalue,
                                  data=rep(paste(tissue[i], "proteins", sep="_"), nrow(proteome.data[[i]])), 
                                  species=rep(species[i], nrow(proteome.data[[i]])))
  
  tot.data <- rbind(tot.data, tot.data.proteins)
}

# Add epigenetics data
arabidopsis.H3K9ac <- read.table(paste(Dir, "DATA/arabidopsis/leaves/epigenetics/H3K9ac_GeneCycle.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
arabidopsis.H3K9ac <- data.frame(rhythm.pvalue = arabidopsis.H3K9ac$default.pvalue,
                                data=rep("leaves_H3K9ac", nrow(arabidopsis.H3K9ac)), 
                                species=rep(species[1], nrow(arabidopsis.H3K9ac)))
arabidopsis.H3K4me3 <- read.table(paste(Dir, "DATA/arabidopsis/leaves/epigenetics/H3K4me3_GeneCycle.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
arabidopsis.H3K4me3 <- data.frame(rhythm.pvalue = arabidopsis.H3K4me3$default.pvalue,
                                 data=rep("leaves_H3K4me3", nrow(arabidopsis.H3K4me3)), 
                                 species=rep(species[1], nrow(arabidopsis.H3K4me3)))

tot.data <- rbind(tot.data, arabidopsis.H3K9ac, arabidopsis.H3K4me3)

tot.data$data <- gsub("NA_", "", tot.data$data)

### PLOTS ###
rhythmPvalDistrib <- ggplot(tot.data, aes(x=rhythm.pvalue)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.8, na.rm = TRUE, bins = 100) + 
  facet_wrap(species ~ data, scales = "free", ncol = 2, labeller = label_parsed) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "rhythm detection p-value (GeneCycle)")

########
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##rhythm_Pvalues_Distribution##") %>%
  body_add_gg(value = rhythmPvalDistrib, width = 4.5, height = 9) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")


########
# Special case for Mouse Proteomics (Liver, Cartilage, Tendon)
mouse.liver.prot.harmonicregr <- read.table(paste(Dir, "DATA/mouse/liver/proteome/protein_harmonicRegression.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
mouse.liver.prot.harmonicregr <- data.frame(rhythm.pvalue = mouse.liver.prot.harmonicregr$Pvalue.protein,
                                            data = rep("`Mouse Liver / proteins`", nrow(mouse.liver.prot.harmonicregr)),
                                            method = rep("method: harmonic_regression", nrow(mouse.liver.prot.harmonicregr)))
mouse.cartilage.prot <- read.csv(paste(Dir, "DATA/mouse/cartilage/proteome/initial_data.csv", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE)
mouse.cartilage.prot <- data.frame(rhythm.pvalue = mouse.cartilage.prot$ARS_pvalue,
                                            data = rep("`Mouse Cartilage / proteins`", nrow(mouse.cartilage.prot)),
                                            method = rep("method: ARS", nrow(mouse.cartilage.prot)))
mouse.tendon.prot <- read.csv(paste(Dir, "DATA/mouse/tendon/proteome/initial_data.csv", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE)
mouse.tendon.prot <- data.frame(rhythm.pvalue = mouse.tendon.prot$ARS_pvalue,
                                   data = rep("`Mouse Tendon / proteins`", nrow(mouse.tendon.prot)),
                                   method = rep("method: ARS", nrow(mouse.tendon.prot)))
tot.data <- rbind(mouse.liver.prot.harmonicregr, mouse.cartilage.prot, mouse.tendon.prot)

specialCasesPvalDistrib <- ggplot(tot.data, aes(x=rhythm.pvalue)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.8, na.rm = TRUE, bins = 100) + 
  facet_wrap(data ~ method, scales = "free", ncol = 2, labeller = label_parsed) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "rhythm detection p-value (original paper)") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))

# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/supplementary_information/specialCase_pval_distrib.pdf", width = 4, height = 4.3)
print(specialCasesPvalDistrib)
dev.off()

########
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##special_cases_Pval_Distrib##") %>%
  body_add_gg(value = specialCasesPvalDistrib, width = 4.5) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")
