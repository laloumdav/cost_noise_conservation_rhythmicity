####################################################################################
### Comparison of averaged AA synthesis cost between both groups of genes ###########
### From AA synth costs of Akashi and Gojobori VS those of Wagner  ###########
####################################################################################
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
                       average.AA.synth.cost_Akashi=NA, 
                       average.AA.synth.cost_Wagner=NA, 
                       protein.length=NA,
                       rhythm.pvalue=NA,
                       nb.genes=NA,
                       data=NA, 
                       species=NA,
                       rhythmic.order=NA,
                       rhythmic=NA)
for (i in 1:length(species.list)) {
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "proteome_data.txt", sep="/"), sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  
  tot.data.proteins <- data.frame(max.expr.value = proteome.data[[i]]$protein_max.level, 
                                  mean.expr.value = proteome.data[[i]]$protein_mean.level, 
                                  average.AA.synth.cost_Akashi = proteome.data[[i]]$aa.synthesis.average.cost_Akashi,
                                  average.AA.synth.cost_Wagner = proteome.data[[i]]$aa.synthesis.average.cost_Wagner, 
                                  protein.length = proteome.data[[i]]$protein.length, 
                                  rhythm.pvalue = proteome.data[[i]]$protein_rhythm.pvalue,
                                  nb.genes = nrow(proteome.data[[i]]),
                                  data = rep("proteins", nrow(proteome.data[[i]])), 
                                  species = rep(species.list[i], nrow(proteome.data[[i]])),
                                  random.group = FALSE)
  
  # order genes by their rhythm p-value
  row.nb <- order(tot.data.proteins$rhythm.pvalue)
  tot.data.proteins[row.nb, "rhythmic.order"] <- 1:nrow(tot.data.proteins)
  # Lets consider the first n percent (percent.rhythmic) as rhythmic genes
  tot.data.proteins$rhythmic <- (tot.data.proteins$rhythmic.order/nrow(tot.data.proteins) <= percent.rhythmic)
  rhythmic.group.size = length(which(tot.data.proteins$rhythmic == TRUE))
  
  # re-order columns 
  tot.data.proteins <- tot.data.proteins[, colnames(tot.data)]
  
  tot.data <- rbind(tot.data, tot.data.proteins)
}
tot.data <- tot.data[-1,]


###  ###  ###  ###  ###  ### 
### PLOTS ### 
subset.tot.data <- rbind(tot.data[, c("species", "rhythmic")], tot.data[, c("species", "rhythmic")])
subset.tot.data[subset.tot.data$rhythmic == TRUE, "rhythmic"] <- "rhythmic"
subset.tot.data[subset.tot.data$rhythmic == FALSE, "rhythmic"] <- "non-rhythmic"
subset.tot.data$average.AA.synth.cost <- c(tot.data$average.AA.synth.cost_Akashi, tot.data$average.AA.synth.cost_Wagner)
subset.tot.data$AA.cost.data <- c(rep("Akashi&Gojobori", nrow(tot.data)), rep("Wagner", nrow(tot.data)))


####################################################################################
### Comparison of averaged AA synthesis cost between both groups of genes ###########
####################################################################################
caption.text <- paste("rhythmic proteins: first ", percent.rhythmic*100, "%", sep="")

AAsynthCostComparedPerGroup <- ggplot(subset.tot.data, aes(species, average.AA.synth.cost, fill = rhythmic)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = .9) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=2.4) +
  facet_wrap(~ AA.cost.data, scales = "fixed", nrow = 1) +
  scale_fill_manual(values=c("seagreen4", "#E7B800")) +
  labs(y = TeX("$\\textbf{\\bar{c}_{AA}}$"), 
       fill="genes group:", 
       caption=caption.text) +
  theme(#aspect.ratio = .9,
    axis.text.y = element_text(size=9), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_text(size=9), 
    axis.title.x = element_blank(),
    #panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey97", colour = "grey90"), 
    axis.line = element_line(colour = "grey32"),
    panel.spacing = unit(.03, "lines"),
    legend.position = "top",
    plot.caption = element_text(size=7, face = "italic", hjust=0))

# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/additional_files/AAsynthCost_compared_perGroup.pdf", height = 6, width = 8)
print(AAsynthCostComparedPerGroup)
dev.off()

##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##AA_synth_cost_Compared_per_group_and_data##") %>%
  body_add_gg(value = AAsynthCostComparedPerGroup, width = 7, height = 4) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")
