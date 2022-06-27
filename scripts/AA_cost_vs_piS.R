library(ggplot2)
library(ggridges)
library(officer)
library(magrittr)
library(ggpubr)
library(latex2exp)
library(readxl)
library(ggrepel)
library(ggpmisc)

main.dir <- "~/Documents/cost_theory_workingspace"

#############
#
#############
average.AA.synth.cost.files <- list.files(paste(main.dir, "DATA/protein_AA_composition/with_pi_s/", sep="/"))
average.AA.synth.cost.files <- paste(main.dir, "DATA/protein_AA_composition/with_pi_s", grep("AA_synthesis_average_cost_per_protein.txt", average.AA.synth.cost.files, value = TRUE), sep="/")

# retrieve pi s data from Supp Table 2 of Romiguier et al. 2014
pi.s <- read_excel(paste(main.dir, "DATA/protein_AA_composition/with_pi_s/piS.xlsx", sep="/"), skip = 2)
pi.s <- as.data.frame(sapply(pi.s, function(x) gsub("\"", "", x)))
colnames(pi.s) <- gsub("\"", "", colnames(pi.s))
pi.s[, 7:(ncol(pi.s)-1)] <- sapply(pi.s[, 7:(ncol(pi.s)-1)], as.numeric)

##############
tot.data <- data.frame()
median.cost <- data.frame()
for (n in 1:length(average.AA.synth.cost.files)) {
  AA.cost.data <- read.table(average.AA.synth.cost.files[n], sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  Species.name <- gsub("~/Documents/cost_theory_workingspace/DATA/protein_AA_composition/with_pi_s/", "", average.AA.synth.cost.files[n])
  Species.name <- gsub("_AA_synthesis_average_cost_per_protein.txt", "", Species.name)
  AA.cost.data$Species <- Species.name
  
  merged.data <- merge(AA.cost.data, pi.s, by="Species")
  tot.data <- rbind(tot.data, merged.data)
  
  median.cost.value_Wagner <- median(merged.data$aa.synthesis.average.cost_Wagner, na.rm = TRUE)
  median.cost.value_Akashi <- median(merged.data$aa.synthesis.average.cost_Akashi, na.rm = TRUE)
  protein.coding.genome.size <- as.numeric(unique(AA.cost.data$protein.coding.genome.size))
  median.cost.values <- cbind(pi.s[pi.s$Species == Species.name, ], median.cost.value_Wagner, median.cost.value_Akashi, protein.coding.genome.size)
  median.cost <- rbind(median.cost, median.cost.values)
}
tot.data$Species <- gsub("_", " ", tot.data$Species)
median.cost$Species <- gsub("_", " ", median.cost$Species)

##############
tot.data$piS <- as.factor(format(tot.data$piS, digits = 2))

piS_AAcostPlot <- ggplot(median.cost, aes(x = piS, y = median.cost.value_Wagner, color = Phylum)) +
  geom_point(size=2) +
  geom_label_repel(aes(label = Species),
                   box.padding   = 0.5, point.padding = 0.1, label.size = .02, 
                   force = 1.7, size = 5,
                   segment.color = 'grey', show.legend=FALSE) +
  labs(x = TeX("synonymous nucleotide diversity ($\\pi_{S}$)"),
       y = TeX("median $\\textbf{\\bar{c}_{AA}}$ (Wagner, 2005)")) +
  #scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  theme(axis.text.y = element_text(size=13), 
        axis.title.y = element_text(size=14), 
        axis.text.x = element_text(size=13), 
        axis.title.x = element_text(size=14), 
        legend.position="top",
        #legend.title = element_blank(), 
        legend.text = element_text(size=14, face = "bold"), 
        axis.line = element_line(colour = "grey32")) 

pdf("~/Documents/cost_theory_workingspace/Fig3_bis.pdf", width = 10, height = 7) # or other device
print(piS_AAcostPlot)
dev.off()


ggplot(median.cost, aes(x=-log(Propagule_Size), y=median.cost.value_Wagner, color = Phylum)) +
  geom_point(size=2) +
  geom_label_repel(aes(label = Species),
                   box.padding   = 0.5, point.padding = 0.1, label.size = .02, 
                   force = 1.7, size = 2.3,
                   segment.color = 'grey', show.legend=FALSE) 
ggplot(median.cost, aes(x=log(Longevity), y=median.cost.value_Wagner, color = Phylum)) +
  geom_point(size=2) +
  geom_label_repel(aes(label = Species),
                   box.padding   = 0.5, point.padding = 0.1, label.size = .02, 
                   force = 1.7, size = 2.3,
                   segment.color = 'grey', show.legend=FALSE) 
ggplot(median.cost, aes(x=log(Adult_Size), y=median.cost.value_Wagner, color = Phylum)) +
  geom_point(size=2) +
  geom_label_repel(aes(label = Species),
                   box.padding   = 0.5, point.padding = 0.1, label.size = .02, 
                   force = 1.7, size = 2.3,
                   segment.color = 'grey', show.legend=FALSE) 
ggplot(median.cost, aes(x=log(Body_Mass), y=median.cost.value_Wagner, color = Phylum)) +
  geom_point(size=2) +
  geom_label_repel(aes(label = Species),
                   box.padding   = 0.5, point.padding = 0.1, label.size = .02, 
                   force = 1.7, size = 2.3,
                   segment.color = 'grey', show.legend=FALSE) 

ggplot(median.cost, aes(x=-log(Body_Mass/protein.coding.genome.size), y=median.cost.value_Wagner)) +
  geom_point(size=2) + 
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  geom_label_repel(aes(label = Species),
                   box.padding   = 0.5, point.padding = 0.1, label.size = .02, 
                   force = 1.7, size = 2.3,
                   segment.color = 'grey', show.legend=FALSE) 

test <- lm(median.cost$median.cost.value_Wagner ~ median.cost$Body_Mass/median.cost$protein.coding.genome.size)
summary(test)

plot1 <- ggplot(median.cost, aes(x = piS, y = median.cost.value_Wagner, color = Social)) +
  geom_point(size=2) +
  geom_label_repel(aes(label = Species),
                   box.padding   = 0.5, point.padding = 0.1, label.size = .02, 
                   force = 1.7, size = 2.3,
                   segment.color = 'grey', show.legend=FALSE) +
  labs(x = TeX("synonymous nucleotide diversity ($\\pi_{S}$)"),
       y = TeX("median $\\textbf{\\bar{c}_{AA}}$ (Wagner, 2005)")) +
  scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  theme(axis.text.y = element_text(size=9), 
        axis.title.y = element_text(size=12), 
        axis.text.x = element_text(size=11), 
        axis.title.x = element_text(size=12), 
        legend.position="top",
        legend.title = element_blank(), 
        legend.text = element_text(size=12, face = "bold"), 
        axis.line = element_line(colour = "grey32")) 

plot2 <- ggplot(tot.data, aes(x = piS, y = aa.synthesis.average.cost_Wagner, fill = Species)) +
  geom_boxplot(outlier.size = .1, size=0.1, width = .9) +
  facet_grid(cols = vars(Social), scales = "free_x") +
  labs(x = TeX("synonymous nucleotide diversity ($\\pi_{S}$)"),
       y = TeX("$\\textbf{\\bar{c}_{AA}}$ (Wagner, 2005)")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme(axis.text.y = element_text(size=9), 
        axis.title.y = element_text(size=12), 
        axis.text.x = element_text(size=9), 
        axis.title.x = element_text(size=12), 
        legend.title = element_text(size=11, face = "bold"),
        legend.key = element_blank(),
        axis.line = element_line(colour = "grey32"),
        plot.caption = element_text(size=7, face = "italic", hjust=0))

plot3 <- ggplot(median.cost, aes(x=piS, y=median.cost.value_Wagner, color=Social, shape=Social))+
  geom_point(size=2) + 
  geom_smooth(method=lm, se=FALSE, fullrange=FALSE)+
  scale_color_manual(values=c('#E69F00', '#56B4E9'))+
  annotate('text', 0.043, 25, label=social.lm.text, 
           hjust=1, size=3.5, color='#56B4E9') +
  annotate('text', 0.04, 26, label=nonsocial.lm.text, 
           hjust=1, size=3.5, color='#E69F00') +
  labs(x = TeX("synonymous nucleotide diversity ($\\pi_{S}$)"),
       y = TeX("median $\\textbf{\\bar{c}_{AA}}$ (Wagner, 2005)")) +
  theme(axis.text.y = element_text(size=9), 
        axis.title.y = element_text(size=12), 
        axis.text.x = element_text(size=11), 
        axis.title.x = element_text(size=12), 
        legend.title = element_blank(), 
        legend.key = element_blank(),
        legend.position="none",
        axis.line = element_line(colour = "grey32"),
        plot.caption = element_text(size=7, face = "italic", hjust=0))

############################
##############
##############
##########################################
# make the complete figure :
library(cowplot)
library(ggpubr)
library(magick)

fig <- plot_grid(plot_grid(plot1, plot3,
                           labels = c("a", "c"),
                           label_size = 17), 
                 plot_grid(plot2, NA,
                           labels = c("b", ""),
                            rel_widths = c(1.5, .5),
                            label_size = 17), 
                 nrow = 2)

pdf("~/Documents/cost_theory_workingspace/Fig3.pdf", width = 11, height = 8) # or other device
print(fig)
dev.off()






