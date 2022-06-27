library(ggplot2)
library(ggridges)
library(officer)
library(magrittr)
library(ggpubr)

species.list <- c("arabidopsis", "cyanobacteria", "ostreococcus", "mouse")
tissue.list <- c("leaves", NA, NA, "liver")

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species.list, tissue.list, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

percent.rhythmic = 0.15


##############
transcriptome.data <- list()
proteome.data <- list()
tot.data <- data.frame(Mass=NA, 
                       protein.length=NA,
                       rhythm.pvalue=NA,
                       nb.genes=NA,
                       data=NA, 
                       species=NA,
                       rhythmic.order=NA,
                       rhythmic=NA,
                       random.group=NA,
                       rhythmic.shapiro.test=NA,
                       random.shapiro.test=NA)
for (i in 1:length(species.list)) {
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "proteome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE)
  
  tot.data.proteins <- data.frame(Mass = proteome.data[[i]]$Mass, 
                                  protein.length = proteome.data[[i]]$protein.length, 
                                  rhythm.pvalue = proteome.data[[i]]$protein.rhythm.pvalue,
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
  # Lets take a random group of genes of size equal to rhythmic genes group size and affect TRUE randomly (size of TRUE genes = rhythmic.group.size) 
  random.true <- sample(1:nrow(tot.data.proteins), size=rhythmic.group.size, replace = FALSE)
  tot.data.proteins[random.true, "random.group"] <- TRUE
  # Retriece p.val of shapiro.test
  tot.data.proteins$rhythmic.shapiro.test <- shapiro.test(subset(tot.data.proteins, rhythmic == TRUE)$Mass)$p.value
  tot.data.proteins$random.shapiro.test <- shapiro.test(subset(tot.data.proteins, random.group == TRUE)$Mass)$p.value
  
  # re-order columns 
  tot.data.proteins <- tot.data.proteins[, colnames(tot.data)]
  
  tot.data <- rbind(tot.data, tot.data.proteins)
}
tot.data <- tot.data[-1,]


### PLOTS ### 
subset.rhythmic.data <- subset(tot.data, rhythmic == TRUE)
subset.rhythmic.data$median.Mass <- NA
subset.random.data <- subset(tot.data, random.group == TRUE)
for (i in 1:length(species.list)) {
  subset.data.tmp <- subset(subset.random.data, species == species.list[i])
  subset.random.data[subset.random.data$species == species.list[i], "median.Mass"] <- median(subset.data.tmp$Mass, na.rm = TRUE)
}
subset.tot.data <- rbind(subset.random.data, subset.rhythmic.data)
subset.tot.data$rhythmic.vs.random <- c( rep("random", nrow(subset.random.data)), rep("rhythmic", nrow(subset.rhythmic.data)))


### 1. molecular weight ###
##############

MolWeightHistPerGroup <- ggplot(subset.tot.data, aes(x=log(Mass), fill=rhythmic.vs.random)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.6, na.rm = TRUE, bins = 40) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  labs(x = "molecular weight (log)", fill="genes group", caption=paste("proportion of rhythmic or random genes group = ", percent.rhythmic*100, "%", sep="")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0))

#shapiro.test
#dat_text <- unique(subset.tot.data[, c("data", "species", "rhythmic.vs.random", "rhythmic.shapiro.test", "random.shapiro.test")])
#library(reshape2)
#dat_text <- melt(dat_text, id=c("data", "species", "rhythmic.vs.random"), value.name = "shapiro.test")
#dat_text$shapiro.test = round(dat_text$shapiro.test,2)

#subset.tot.data <- melt(subset.tot.data, id=colnames(subset.tot.data)[-grep("shapiro", colnames(subset.tot.data))], value.name = "shapiro.test")

MolWeightQQplot <- ggplot(subset.tot.data, aes(sample=log(Mass))) + 
  stat_qq(aes(color=rhythmic.vs.random), size=0.1) + 
  stat_qq_line(aes(color=rhythmic.vs.random)) +
  scale_color_manual("genes group", values = c("#00AFBB", "#E7B800")) +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  labs(caption=paste("molecular weight (log)", percent.rhythmic*100, "%", sep="")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0))



####################################################################################
### Comparison of averaged AA synthesis cost between both groups of genes ###########
####################################################################################
alternative.wilcox.option <- "two.sided"  # choice between "two.sided", "less", "greater"

caption.text <- paste("Box-plot log scaled",
                      "\n",
                      paste("Wilcox test (alternative='", alternative.wilcox.option, "'):", sep=""),
                      "ns: p > 0.05",
                      "*: p <= 0.05",
                      "**: p <= 0.01",
                      "***: p <= 0.001",
                      "****: p <= 0.0001",
                      "\n",
                      paste("proportion of rhythmic or random genes group = ", percent.rhythmic*100, "%", sep=""),
                      sep="\n")

MolWeightComparedPerGroup <- ggboxplot(subset.tot.data, x = "rhythmic.vs.random", 
                                         y = "Mass", yscale = "log2",
                                         color = "rhythmic.vs.random", palette = c("#00AFBB", "#E7B800"),
                                         width = 0.4,
                                         #add = "jitter", 
                                         short.panel.labs = T) + 
  geom_hline(aes(yintercept = median.Mass), linetype = 2, color="grey39") +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = alternative.wilcox.option), 
                     label="p.signif" , label.x.npc = 0.6, label.y.npc = 0.9, na.rm = TRUE) +
  labs(y = "molecular weight", color="genes group:", caption=caption.text) +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size=6), 
        axis.text.x = element_text(size=10), 
        axis.title.y = element_text(size=12), 
        plot.caption = element_text(size=11, face = "italic", hjust=0))

##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##molecular_weight_histograms_per_group##") %>%
  body_add_gg(value = MolWeightHistPerGroup, width = 7, height = 2) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##molecular_weight_QQplots_per_group##") %>%
  body_add_gg(value = MolWeightQQplot, width = 7, height = 2) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_reach(keyword = "##molecular_weight_compared_between_groups##") %>%
  body_add_gg(value = MolWeightComparedPerGroup, width = 6.8, height = 5.9) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")



##############
# Also save it as an R object:
caption.text <- paste(paste("Wilcox test (alternative='", alternative.wilcox.option, "'):", sep=""),
                      sep="\n")

MolWeightComparedPerGroup <- ggboxplot(subset.tot.data, x = "rhythmic.vs.random", 
                                         y = "Mass", yscale = "log2",
                                         fill = "rhythmic.vs.random", palette = c("seagreen4", "#E7B800"),
                                         width = 0.4, alpha=0.8, size=0.2, outlier.size=0.1,
                                         #add = "jitter", 
                                         short.panel.labs = T) + 
  geom_hline(aes(yintercept = median.Mass), linetype = 2, color="grey39") +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  stat_compare_means(method = "wilcox.test", method.args = list(alternative = alternative.wilcox.option), 
                     label="p.signif" , label.x.npc = 0.6, label.y.npc = 0.9, na.rm = TRUE, size=2.2) +
  labs(y = "molecular weight", fill="genes group:", caption=caption.text) +
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

out.rds.name <- paste("~/Documents/cost_theory_workingspace/rds_objects", "MolWeightComparedPerGroup.rds", sep="/")
saveRDS(MolWeightComparedPerGroup, file = out.rds.name)


