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
                       average.AA.synth.cost=NA, 
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
                                  average.AA.synth.cost = proteome.data[[i]]$aa.synthesis.average.cost_Wagner, # Use Akashi and Gojobori AA synth costs data
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
# Total cost per gene calculated with the mean of expression values of all timepoints normalized by Z-score
tot.data$total.cost = tot.data$mean.expr.value * tot.data$protein.length * (tot.data$average.AA.synth.cost - 1)
# To get scaled graphs, lets have scores: 
#for (i in 1:length(species.list)) {
#  tot.data[tot.data$species == species.list[i], "total.cost"] = tot.data[tot.data$species == species.list[i], "total.cost"] / max(tot.data[tot.data$species == species.list[i], "total.cost"], na.rm=TRUE)
#}

###  ###  ###  ###  ###  ### 
### PLOTS ### 
subset.rhythmic.data <- subset(tot.data, rhythmic == TRUE)
subset.rhythmic.data$median.average.AA.synth.cost <- NA
subset.rhythmic.data$median.total.cost <- NA
subset.nonrhythmic.data <- subset(tot.data, rhythmic == FALSE)
for (i in 1:length(species.list)) {
  subset.data.tmp <- subset(subset.nonrhythmic.data, species == species.list[i])
  subset.nonrhythmic.data[subset.nonrhythmic.data$species == species.list[i], "median.average.AA.synth.cost"] <- median(subset.data.tmp$average.AA.synth.cost, na.rm = TRUE)
  subset.nonrhythmic.data[subset.nonrhythmic.data$species == species.list[i], "median.total.cost"] <- median(subset.data.tmp$total.cost, na.rm = TRUE)
}
subset.tot.data <- rbind(subset.nonrhythmic.data, subset.rhythmic.data)
subset.tot.data$genes.gp <- c( rep("non-rhythmic", nrow(subset.nonrhythmic.data)), rep("rhythmic", nrow(subset.rhythmic.data)))


### 1. averaged AA synthesis cost ###
##############

AAsynthCostHistPerGroup <- ggplot(subset.tot.data, aes(x=average.AA.synth.cost, fill=genes.gp)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.6, na.rm = TRUE, bins = 40) + 
  scale_fill_manual(values = c("seagreen4", "#E7B800")) + 
  scale_color_manual(values = c("seagreen4", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  labs(x = TeX("$\\textbf{\\bar{c}_{AA}}$ (Wagner, 2005)"), fill="genes group", caption=paste("proportion of rhythmic genes group = ", percent.rhythmic*100, "%", sep="")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0),
        strip.text = element_text(size=14))


#shapiro.test
#dat_text <- unique(subset.tot.data[, c("data", "species", "rhythmic.vs.random", "rhythmic.shapiro.test", "random.shapiro.test")])
#library(reshape2)
#dat_text <- melt(dat_text, id=c("data", "species", "rhythmic.vs.random"), value.name = "shapiro.test")
#dat_text$shapiro.test = round(dat_text$shapiro.test,2)

#subset.tot.data <- melt(subset.tot.data, id=colnames(subset.tot.data)[-grep("shapiro", colnames(subset.tot.data))], value.name = "shapiro.test")

AAsynthCostQQplot <- ggplot(subset.tot.data, aes(sample=average.AA.synth.cost)) + 
  stat_qq(aes(color=genes.gp), size=0.1) + 
  stat_qq_line(aes(color=genes.gp)) +
  scale_color_manual("genes group", values = c("seagreen4", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  labs(caption=paste("proportion of rhythmic genes group = ", percent.rhythmic*100, "%", sep="")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        plot.caption = element_text(size=11, face = "italic", hjust=0),
        strip.text = element_text(size=14))


# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/supplementary_information/AAsynthCost_hist_perGroup.pdf", height = 3, width = 9)
print(AAsynthCostHistPerGroup)
dev.off()
pdf("~/Documents/cost_theory_workingspace/supplementary_information/AAsynthCost_QQplot.pdf", height = 3, width = 9)
print(AAsynthCostQQplot)
dev.off()

##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##AA_synth_cost_histograms_per_group##") %>%
  body_add_gg(value = AAsynthCostHistPerGroup, width = 7, height = 2) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##AA_synth_cost_QQplots_per_group##") %>%
  body_add_gg(value = AAsynthCostQQplot, width = 7, height = 2) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")

####################################################################################
### Comparison of averaged AA synthesis cost between both groups of genes ###########
####################################################################################
caption.text <- paste("rhythmic proteins: first ", percent.rhythmic*100, "%", sep="")

AAsynthCostComparedPerGroup <- ggplot(subset.tot.data, aes(species, average.AA.synth.cost, fill = genes.gp)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = .9) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=2.4) +
  scale_fill_manual(values=c("seagreen4", "#E7B800")) +
  #labs(y = TeX("$\\textbf{\\bar{c}_{AA}}$ (Akashi and Gojobori, 2002)"), 
  labs(y = TeX("$\\textbf{\\bar{c}_{AA}}$ (Wagner, 2005)"), 
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

out.rds.name <- paste("~/Documents/cost_theory_workingspace/rds_objects", "AAsynthCostComparedPerGroup.rds", sep="/")
saveRDS(AAsynthCostComparedPerGroup, file = out.rds.name)





### 2. total cost ###
##############

####################################################################################
### Comparison of the total cost between both groups of genes ###########
####################################################################################
caption.text <- paste("Box-plot log2 scaled",
                      paste("rhythmic genes: first ", percent.rhythmic*100, "%", sep=""), sep="\n")

totalCostComparedPerGroup <- ggplot(subset.tot.data, aes(species, log2(total.cost), fill = genes.gp)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = 1) +
  facet_wrap(~ species, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", #method.args = list(alternative = alternative.wilcox.option), 
                     label="p.format" , #label.x.npc = 0.4, label.y.npc = 0.95, na.rm = TRUE
                     size=2.4) +
  geom_hline(aes(yintercept = log2(median.total.cost)), linetype = 2, color="grey39", size=.2) +
  scale_fill_manual(values=c("seagreen4", "#E7B800")) +
  labs(y = TeX("$\\textbf{C_{p}}$"),
       fill="genes group:",
       caption=caption.text) +
  theme(#aspect.ratio = .9,
        axis.text.y = element_text(size=6), 
        axis.title.y = element_text(size=12), 
        axis.text.x = element_text(size=11), 
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(.07, "lines"),
        #axis.line.x = element_line(colour = "grey"),
        #axis.line.y = element_line(colour = "grey"),
        legend.key = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "grey32"),
        strip.text = element_blank(),
        legend.position = "top",
        plot.caption = element_text(size=7, face = "italic", hjust=0))

out.rds.name <- paste("~/Documents/cost_theory_workingspace/rds_objects", "totalCostComparedPerGroup.rds", sep="/")
saveRDS(totalCostComparedPerGroup, file = out.rds.name)

