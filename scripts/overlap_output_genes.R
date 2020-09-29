library(venn)
library(wesanderson)


if (dir.exists("/scratch/cluster/monthly/dlaloum/")==TRUE){
  main.dir <- "/scratch/cluster/monthly/dlaloum/Documents/cost_theory_workingspace"
} else {
  main.dir <- "~/Documents/cost_theory_workingspace"
}

#N = 100

file.dir <- paste(paste(main.dir, "DATA/mauvoisin", sep = "/"), "/", sep="")

# Original paper algorithm:
file.name <- paste(file.dir, "complete_data_from_original_paper.csv", sep = "")
original <- read.csv(file.name, check.names = FALSE, stringsAsFactors = FALSE, fill=TRUE)
original <- original[, c("Gene.names", "Protein.names", "Pvalue.protein")]
colnames(original) <- c("Gene.Name", "Protein.Name", "pvalue")
original$algorithm <- rep("paper", nrow(original))
original$protein_number <- 1:nrow(original)
# GeneCycle:
file.name <- paste(file.dir, "protein_GeneCycle.txt", sep = "")
GeneCycle <- read.table(file.name, head=TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, fill=TRUE)
GeneCycle <- GeneCycle[, c("Protein.Name", "default.pvalue")]
colnames(GeneCycle)[2] <- "pvalue"
GeneCycle$Gene.Name <- original$Gene.Name
GeneCycle <- GeneCycle[, c("Gene.Name", "Protein.Name", "pvalue")]
GeneCycle$algorithm <- rep("GeneCycle", nrow(GeneCycle))
GeneCycle$protein_number <- 1:nrow(GeneCycle)

# Mix data
mixed.data <- rbind(original, GeneCycle)

algorithm.list <- unique(mixed.data$algorithm)


# Function to unlist intersections of rhythmic proteins between algorithms which keep proteins.number 
FromList <- function(input.list, elements) {
  data <- unlist(lapply(input.list, function(x) { x <- as.vector(match(elements, x)) } ))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input.list), byrow = F))
  names(data) <- names(input.list)
  data$elements <- elements
  return(data)
}



N.max = 4000
intersect.ratio.list <- data.frame(`number of top rhythmic proteins`=1:N.max, `propotion in common`=NA)
for (n in 1:N.max) {
  N = n
  
  # keep first "N" rhythmic proteins ordered:
  list.rhythmic.proteins.ordered <- list()
  for (i in 1:length(algorithm.list)){
    data <- subset(mixed.data, algorithm == algorithm.list[i])
    # subset of first "N" rhythmic proteins
    data <- data[order(data[, "pvalue"]), ]
    data <- data[1:N, ]
    list.rhythmic.proteins.ordered[[algorithm.list[i]]] <- unique(data$protein_number)
  }
  
  # Take "N" proteins kept with the highest threshold
  proteins.number <- unique(unlist(list.rhythmic.proteins.ordered))
  
  intersect.table.N.rhythmicproteins <- FromList(list.rhythmic.proteins.ordered, proteins.number)
  intersection.length <- length(which(intersect.table.N.rhythmicproteins$paper == intersect.table.N.rhythmicproteins$GeneCycle))
  # Ratio of the intersection:
  intersect.ratio <- intersection.length/N
  intersect.ratio.list$propotion.in.common[N] <- intersect.ratio
}


ggplot(intersect.ratio.list, aes(x=number.of.top.rhythmic.proteins, y=propotion.in.common)) +
  geom_line() 
  geom_smooth(aes(x=genes_ordered, y=proportion, col=algorithm), method = "loess", span = 0.06, se = FALSE, size=0.7) +
  geom_point(data = cutoff.table, aes(col=algorithm, fill=algorithm), na.rm = TRUE, 
             size = 2.5, shape=23, color="black") +
  theme_light(base_size = 0.25) +
  scale_fill_manual(values = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue", Naive="black")) +
  scale_color_manual(values = c(JTK=colors.list[1], LS=colors.list[6], empJTK=colors.list[7], ARS=colors.list[4], RAIN=colors.list[5], meta2d=colors.list[2], GeneCycle="blue", Naive="black")) +
  guides(colour = guide_legend(override.aes = list(shape = NA, size = 2))) +
  labs(x = "x\nnumber of orthologs detected rhythmic", y = paste("proportion which are\nrhythmic in", species2.name, tissue2, sep=" ")) +
  theme(axis.text.x = element_text(size=12, hjust = 0.5, vjust = 0.5, margin = margin(l = 20)),
        axis.title.x = element_text(size=14, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size=12, hjust = 0.5, vjust = 0.5, margin = margin(l = 20)),
        axis.title.y = element_text(size=14, hjust = 0.5, vjust = 0.5),
        legend.position = "none")





################
# Venn diagram #
################
colors.list <- c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest2"))
#title <- paste(species, tissue, sep = " ")
#subtitle.1 <- paste(p.val, "<", threshold.1, sep = " ")
#subtitle.2 <- paste(p.val, "<", threshold.2, sep = " ")
#tot.genes <- paste("tot genes =", tot.genes.nb, sep = " ")
#p.adj <- paste("p.adj.method =", p.adj.method, sep = " ")
# order for colors for the venn diagram
col.algo <- colnames(intersect.table.N.rhythmicproteins)
col.algo <- col.algo[col.algo %in% c("GeneCycle", "paper")]
col.algo[col.algo == "paper"] <- colors.list[3]
col.algo[col.algo == "GeneCycle"] <- "blue"


###  Exit to powerpoint file
venn(intersect.table.N.rhythmicproteins[colnames(intersect.table.N.rhythmicproteins) != "elements"], 
                    ilabels = TRUE, 
                    zcolor = col.algo, 
                    opacity = 0.5,
                    cexsn = 0.9,
                    cexil=1)



