######################################################################
######################################################################
### Gene expression variation over time, of some genes ###############
######################################################################
######################################################################
library(ggplot2)
require(RColorBrewer)
library(esquisse)

main.dir <- "~/Documents/cost_theory_workingspace/"
dataset <- read.csv(paste(main.dir, "DATA/theoretical_genes.csv", sep=""), head=TRUE)
#write.csv(data, paste(main.dir, "DATA/theoretical_genes.csv", sep=""), quote = F, row.names = F)

# Relative expression calcul
t.dataset <- as.data.frame(t(dataset[,-1]))
colnames(t.dataset) <- dataset$theoretical_gene
#t.dataset$tot_expr <- apply(t.dataset, 1, sum)
t.dataset <- t.dataset / apply(t.dataset, 1, sum)

# stack function manner : 
final.table <- stack(dataset, select = -theoretical_gene)
final.table.relativ <- stack(as.data.frame(t(t.dataset)))
final.table <- cbind(final.table.relativ[,1], final.table)
final.table$theoretical_gene <- rep(dataset$theoretical_gene, ncol(dataset)-1)
colnames(final.table) <- c("relativ_expr", "expr_value", "time_point", "theoretical_gene")
final.table$time_point <- gsub("CT|ZT", "", final.table$time_point)
final.table$time_point <- as.numeric(final.table$time_point)


########################
# Normalize all data to follow a null distribution
# With z-score
# Such as: for gene i (I mean one microarray data for instance, not normalized per unique gene), we have j time-points
# mean expression of gene.i = mi = sum(gene.ij)/j
# mean expression of full dataset = m
# standard deviation of full dataset = sd
# In our case, z-score of gene.i is: Zi=(mi-m)/sd
# Then, we normalize each data of gene.i
# For this, we will substract xi to each data of gene.i: gene.ij.normalized = gene.ij - xi
# such as: sum(gene.ij.normalized)=Zi
# => sum(gene.ij) - j*xi = Zi
# => xi = (sum(gene.ij)-Zi)/j = (j*mi - Zi)/j
# => xi = mi - (Zi/j)
########################
m = mean(as.numeric(unlist(dataset[,-1])), na.rm = TRUE)
sd = sd(as.numeric(unlist(dataset[,-1])), na.rm = TRUE)
normalized.raw.dataset.genes <- dataset
for (i in 1:nrow(dataset)) {
  gene.i <- as.numeric(dataset[i, -1])
  mi = mean(gene.i, na.rm = TRUE)
  Zi = (mi-m)/sd
  xi = mi - (Zi/length(gene.i))
  for (j in 1:length(gene.i)) {
    gene.ij <- gene.i[j]
    gene.ij.normalized = gene.ij - xi
    
    normalized.raw.dataset.genes[i, 1+j] <- gene.ij.normalized
  }
}
normalized.raw.dataset.genes[,-1] <- normalized.raw.dataset.genes[,-1] + abs(round(min(normalized.raw.dataset.genes[,-1])))
normalized.raw.dataset.genes$mean.expr <- apply(normalized.raw.dataset.genes[,-1], 1, mean)
  
# stack function manner : 
final.table.zscore <- stack(normalized.raw.dataset.genes, select = -theoretical_gene)
final.table <- cbind(final.table.zscore[,1], final.table)
colnames(final.table)[1] <- "zscore_expr"


ggplot(final.table, aes(x = time_point)) +
  geom_line(aes(y=expr_value), color = '#0c4c8a') +
  labs(x = 'time-point (hours)', y = 'gene expression') +
  coord_cartesian(ylim = c(0, 500)) +
  facet_wrap( ~ theoretical_gene, scales = "free") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y =element_text(size=6),
        axis.title.x = element_text(face="bold", hjust=0.5),
        axis.title.y = element_text(face="bold", hjust=0.5),
        legend.text = element_text(size = 11),
        legend.background = element_rect(color = "black", linetype = "solid", size = 0.2),
        legend.key = element_rect(colour = "white", fill = NA), 
        legend.title = element_text(face="bold", hjust=0.5, size=12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggplot(final.table, aes(x = time_point)) +
  geom_line(aes(y=relativ_expr), color = 'red') +
  labs(x = 'time-point (hours)', y = 'gene expression') +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap( ~ theoretical_gene, scales = "free") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y =element_text(size=6),
        axis.title.x = element_text(face="bold", hjust=0.5),
        axis.title.y = element_text(face="bold", hjust=0.5),
        legend.text = element_text(size = 11),
        legend.background = element_rect(color = "black", linetype = "solid", size = 0.2),
        legend.key = element_rect(colour = "white", fill = NA), 
        legend.title = element_text(face="bold", hjust=0.5, size=12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

ggplot(final.table, aes(x = time_point)) +
  geom_line(aes(y=zscore_expr), color = 'red') +
  labs(x = 'time-point (hours)', y = 'gene expression') +
  coord_cartesian(ylim = c(0, 500)) +
  facet_wrap( ~ theoretical_gene, scales = "free") +
  theme(axis.text.x=element_text(size=8),
        axis.text.y =element_text(size=6),
        axis.title.x = element_text(face="bold", hjust=0.5),
        axis.title.y = element_text(face="bold", hjust=0.5),
        legend.text = element_text(size = 11),
        legend.background = element_rect(color = "black", linetype = "solid", size = 0.2),
        legend.key = element_rect(colour = "white", fill = NA), 
        legend.title = element_text(face="bold", hjust=0.5, size=12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))



