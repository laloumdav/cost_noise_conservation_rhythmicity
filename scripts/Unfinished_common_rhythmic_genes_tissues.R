# Installation du package "plyr" pour utiliser la fonction count() :
library(plyr, ggplot2)
library(ggplot2)
library(reshape) 
library(scales)
library("RColorBrewer")

### Etablissement de la fonction multiplot :
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


main.dir <- "~/Documents/rhythm_detection_benchmark"

######
species <- "mouse_microarray"
tissue.list <- c("liver", "lung", "kidney", "muscle", "aorta", "brain_stem", "adrenal_gland", "brown_adipose",
                 "cerebellum", "heart", "white_adipose")

p.val <- "default.pvalue"
algorithm <- "GeneCycle"

full.dataframe <- data.frame()
for (t in 1:length(tissue.list)) {
  tissue <- tissue.list[t]
  file.dir <- paste(paste(main.dir, "DATA", species, tissue, sep = "/"), "/", sep="")
  normalized.pval <- read.table(paste(file.dir, "normalized_default.pvalue/normalized_pvalue_per_gene.txt", sep=""), head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  if ("pvalue_brown" %in% colnames(normalized.pval)) {
    normalized.pval <- subset(normalized.pval, algorithm == "GeneCycle")[, c("ID", "pvalue_brown")]
    colnames(normalized.pval) <- c("ID", tissue)
  } else {
    normalized.pval <- subset(normalized.pval, algorithm == "GeneCycle")[, -grep("algorithm", colnames(normalized.pval))]
    colnames(normalized.pval) <- c("ID", tissue)
  }
  
  if (t==1) {
    full.dataframe <- normalized.pval
  } else {
    full.dataframe <- merge(full.dataframe, normalized.pval)
  }
}

full.dataframe$rhythmic_count <- apply(full.dataframe[, -1], 1, FUN = function(x){return(length(which(x <= 0.01)))})
# keep genes rhythmic in at least 1 tissue
full.dataframe <- subset(full.dataframe, rhythmic_count >=1)
full.dataframe[full.dataframe$rhythmic_count == 1, "rhythmic_count"] <- "tissue-specific"
full.dataframe[full.dataframe$rhythmic_count == 2, "rhythmic_count"] <- "2 tissues"
full.dataframe[full.dataframe$rhythmic_count == 3, "rhythmic_count"] <- "3 tissues"
full.dataframe[full.dataframe$rhythmic_count >3, "rhythmic_count"] <- "common"

# Kegg circadian genes
kegg.circadian.genes <- read.table(paste(main.dir, "DATA", species, "kegg_circadian_rhythmc_genes.txt", sep = "/"), 
                                   head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# Add clock genes info
full.dataframe[which(full.dataframe$ID %in% kegg.circadian.genes$Gene.ID), "clock gene"] <- TRUE
# remove gene.IDs
full.dataframe <- full.dataframe[, -1]
# number of rhythmic genes per tissue
full.dataframe["number of rhythmic genes", c(-ncol(full.dataframe)+1, -ncol(full.dataframe))] <- apply(full.dataframe[, c(-ncol(full.dataframe)+1, -ncol(full.dataframe))], 2, function(x){return(length(which(x<=0.01)))})
# Apply it to tissues
for (i in 1:nrow(full.dataframe)) {
  full.dataframe[i, which(full.dataframe[i, -ncol(full.dataframe)] <= 0.01)] <- full.dataframe[i, "rhythmic_count"]
  full.dataframe[i, which(full.dataframe[i, -ncol(full.dataframe)] > 0.01)] <- NA
}

full.dataframe

ggplot(Tableau_plot,aes(x = Tissue, y = `Number of rhythmic genes`, fill = specificity,width=.6)) +
  geom_bar(position = "fill",stat = "identity") +
  geom_bar(stat="identity")+
  theme(legend.text = element_text(size=13,face="bold")) +
  coord_flip() +
  theme(axis.text.y = element_text(size=13,face="bold")) +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5, size = 12)) +
  scale_fill_manual(values=c(` tissue specific`="plum3", `2 tissues`="palegreen2", `3 tissues`="indianred1", `Common`="steelblue2")) +
  ggtitle("Distribution of the tissue-specificity\nof circadian genes in Mouse") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  geom_point(data=Tableau_plot,aes(Tissue,valeur.Ordonnee, size=Number.of.CoreClock.Genes)) +
  scale_size_continuous(range = c(-0.8,7)) 






# Liste de tous les genes circadien dans l'organisme: 
List_All_circad <- NULL
#Circad_Tissue_List <- list()
for(i in 1:Nbtissue)
{ 
  filename <- paste(Listtissues[i],sep="")  
  Mouse_Tissue_1_stats <- read.table(filename, header = T)
  Mouse_Tissue_1_stats <- Mouse_Tissue_1_stats[,c("EnsemblID", "JTK_q.value")]
  Mouse_Tissue_1_signif <- subset(Mouse_Tissue_1_stats, JTK_q.value<0.05) 
  Mouse_Tissue_Circad <- as.data.frame(unique(Mouse_Tissue_1_signif$EnsemblID))
  colnames(Mouse_Tissue_Circad) <- "Circad.GeneID"
  #Circad_Tissue_List[[i]] <- Mouse_Tissue_Circad
  
  List_All_circad <- rbind(List_All_circad, Mouse_Tissue_Circad)
}

# TEST 
Mouse_Tissue_1_stats[grep("ENSMUSG00000028957", Mouse_Tissue_1_stats$EnsemblID),]
length(unique(Mouse_Tissue_1_stats$EnsemblID))
Gene_protcoding <- merge(Correspondance_ID, Mouse_Tissue_1_stats, by.x="Gene.ID", by.y="EnsemblID")
Gene_protcoding_bis <- Gene_protcoding[,c("Gene.ID", "Gene.Type")]
Gene_protcoding_bis <- subset(Gene_protcoding_bis, Gene.Type=="protein_coding")
length(unique(Gene_protcoding_bis$Gene.ID))
head(Gene_protcoding_bis)

# Nombre de fois que chaque gene circadien est retrouve dans 1 tissu different :
library(plyr)
Count_List <- count(List_All_circad, 'Circad.GeneID')
# Densite (divise par le nombre total de tissus) :
Count_List$Density <- Count_List$freq/14
# Nombre total de genes circadiens : 
a <- nrow(Count_List)
# Nombre de genes circadiens commun en fonction du nombre de tissus :
Freq_List <- count(Count_List, 'freq')
colnames(Freq_List)[1] <- "Nb of tissues"
colnames(Freq_List)[2] <- "Nb of CircadGenes"
# Genes rhythmic in all tissues: 
Genes_rhythmic_AllTissues <- Count_List[grep(14, Count_List$freq), "Circad.GeneID"]

# Pour chaque gene de la liste par tissu, on regarde le nombre de "count" qui leur sont associes :
List_count_totale <- matrix()
List_count_CoreClock <- matrix()
for(i in 1:Nbtissue) { 
  filename <- paste(Listtissues[i],sep="") 
  Mouse_Tissue_1_stats <- read.table(filename, header = T)
  Mouse_Tissue_1_stats <- Mouse_Tissue_1_stats[,c("EnsemblID", "JTK_q.value")]
  Mouse_Tissue_1_signif <- subset(Mouse_Tissue_1_stats, JTK_q.value<0.05) 
  Mouse_Tissue_Circad <- as.data.frame(unique(Mouse_Tissue_1_signif$EnsemblID))
  colnames(Mouse_Tissue_Circad) <- "Circad.GeneID"
  
  # On Merge avec la count-liste :
  count_List_par_tissue <- merge(Mouse_Tissue_Circad,Count_List, by="Circad.GeneID")
  # On change le nom de la colonne freq, car sinon il fait la somme apr??s : 
  colnames(count_List_par_tissue)[2] <- "frequency"
  # Et on fait la proportion de genes commun : 
  Freq_list_par_tissue <- count(count_List_par_tissue,'frequency' )
  colnames(Freq_list_par_tissue)[1] <- "Nb.of.tissues"
  colnames(Freq_list_par_tissue)[2] <- "Nb.of.CircadGenes" 
  # On considere si un gene au moins dans 4 tissus alors il est "commun": 
  count_tot_common_genes <-  sum(Freq_list_par_tissue$Nb.of.CircadGenes[match(4:14,Freq_list_par_tissue$Nb.of.tissues)],na.rm = TRUE)
  Freq_list_par_tissue <- as.data.frame(Freq_list_par_tissue[1:4,2])
  colnames(Freq_list_par_tissue) <- Listtissues_Bis_Title[i]
  Freq_list_par_tissue[4,] <- count_tot_common_genes
  row.names(Freq_list_par_tissue) <- c(" tissue specific","2 tissues","3 tissues","Common")
  # On cre un fichier avec juste les "core clock genes" present dans ce tissus et avec le nombre de tissus tot dans lesquels on les retrouve,
  # et on lui fait le meme procede :
  CoreClock_genes_par_tissue <- merge(count_List_par_tissue,Mouse_BioSystem_Circadian_Genes, by.x = "Circad.GeneID", by.y = "Gene.ID" )
  Freq_CoreClock_genes_par_tissue <- count(CoreClock_genes_par_tissue,'frequency' )
  colnames(Freq_CoreClock_genes_par_tissue)[1] <- "Nb.of.tissues"
  colnames(Freq_CoreClock_genes_par_tissue)[2] <- "Nb.of.CircadGenes" 
  
  value1 <- Freq_CoreClock_genes_par_tissue$Nb.of.CircadGenes[match(1,Freq_CoreClock_genes_par_tissue$Nb.of.tissues)]
  value2 <- Freq_CoreClock_genes_par_tissue$Nb.of.CircadGenes[match(2,Freq_CoreClock_genes_par_tissue$Nb.of.tissues)]
  value3 <- Freq_CoreClock_genes_par_tissue$Nb.of.CircadGenes[match(3,Freq_CoreClock_genes_par_tissue$Nb.of.tissues)]
  Tot_Commum <-  sum(Freq_CoreClock_genes_par_tissue$Nb.of.CircadGenes[match(4:14,Freq_CoreClock_genes_par_tissue$Nb.of.tissues)],na.rm = TRUE)
  # On fait le tableau final pour les core clock genes:
  a <- c(value1,value2,value3,Tot_Commum)
  row_names <- c(" tissue specific","2 tissues","3 tissues","Common")
  FreqFinal_CoreClock_tissue <- data.frame(a, row.names = row_names)
  FreqFinal_CoreClock_tissue[is.na(FreqFinal_CoreClock_tissue)] <- 0
  colnames(FreqFinal_CoreClock_tissue) <- Listtissues_Bis_Title[i]
  
  # On fait la liste totale : 
  List_count_totale <- cbind(List_count_totale,Freq_list_par_tissue)
  List_count_CoreClock <- cbind(List_count_CoreClock,FreqFinal_CoreClock_tissue)
}

List_count_totale <- List_count_totale[,-1]
Reverse_List_count_totale <- as.data.frame(t(List_count_totale))
# On ajoute une colonne avec le Nb total de genes circadiens par tissus :
Reverse_List_count_totale$Tot.Genes <- Reverse_List_count_totale$`2 tissues`+ Reverse_List_count_totale$`3 tissues`+ Reverse_List_count_totale$Common+ Reverse_List_count_totale$` tissue specific`
# On garde la liste des tissus avec leur nombre total de genes circadiens pour d'autres analyses:
Tot_CircadGenes_par_tissue <- as.data.frame(Reverse_List_count_totale[5])
colnames(Tot_CircadGenes_par_tissue) <- "Tot.CircadGenes"
# On range dans l'ordre en fonction du Nombre totale de genes circadiens : 
Reverse_List_count_totale_ordered <- Reverse_List_count_totale[order(Reverse_List_count_totale$Tot.Genes,decreasing = FALSE),]
Reverse_List_count_totale_ordered <- Reverse_List_count_totale_ordered[,-5]
# On cr?? un tableau avec les futurs abscisses ou ordonnees des points :
Dots_Ordonnees <- matrix()


Dots_Ordonnees$tissueSpeOrdonnee <- Reverse_List_count_totale_ordered$`2 tissues` + Reverse_List_count_totale_ordered$`3 tissues` + Reverse_List_count_totale_ordered$Common + Reverse_List_count_totale_ordered$` tissue specific`/2
Dots_Ordonnees$tissues2.Ordonnee <- Reverse_List_count_totale_ordered$`2 tissues`/2 + Reverse_List_count_totale_ordered$`3 tissues` + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$tissues3.Ordonnee <- Reverse_List_count_totale_ordered$`3 tissues`/2 + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$Common.Ordonnee <- Reverse_List_count_totale_ordered$Common/2
Dots_Ordonnees <- as.data.frame(Dots_Ordonnees)
Dots_Ordonnees <- Dots_Ordonnees[,-1]
rownames(Dots_Ordonnees) <- rownames(Reverse_List_count_totale_ordered)
Dots_Ordonnees <- as.data.frame(t(Dots_Ordonnees))
Dots_Ordonnees <- melt(cbind(Dots_Ordonnees, ind = rownames(Dots_Ordonnees)), id.vars = c('ind'))
colnames(Dots_Ordonnees) <- c("Ordonnes", "Tissue", "valeur.Ordonnee")
# Puis seulement maintenant on r??verse puis on transform par melt() :
Reverse_List_count_totale_ordered <- as.data.frame(t(Reverse_List_count_totale_ordered))
Tableau_plot <- melt(cbind(Reverse_List_count_totale_ordered, ind = rownames(Reverse_List_count_totale_ordered)), id.vars = c('ind'))
# On ajoute les valeur des ordonnees pour les futurs points noirs : 
Tableau_plot$valeur.Ordonnee <- Dots_Ordonnees[,3]
colnames(Tableau_plot) <- c("specificity","Tissue","Number of rhythmic genes","valeur.Ordonnee")


# On fait pareil avec la partie des core clock genes : 
List_count_CoreClock <- List_count_CoreClock[,-1]
Reverse_List_count_CoreClock <- as.data.frame(t(List_count_CoreClock))
# On ajoute une colonne avec le Nb total de genes circadiens par tissus :
Reverse_List_count_CoreClock$Tot.Genes <- Reverse_List_count_CoreClock$`2 tissues`+ Reverse_List_count_CoreClock$`3 tissues`+ Reverse_List_count_CoreClock$Common+ Reverse_List_count_CoreClock$` tissue specific`
Reverse_CoreClock_totale_ordered <- Reverse_List_count_CoreClock[order(Reverse_List_count_CoreClock$Tot.Genes,decreasing = FALSE),]
Reverse_CoreClock_totale_ordered <- Reverse_CoreClock_totale_ordered[,-5]
# Puis seulement maintenant on r??verse puis on transform par melt() :
Reverse_CoreClock_totale_ordered <- as.data.frame(t(Reverse_CoreClock_totale_ordered))
Tableau_CoreClock <- melt(cbind(Reverse_CoreClock_totale_ordered, ind = rownames(Reverse_CoreClock_totale_ordered)), id.vars = c('ind'))

# On ajoute alors les valeurs des core clock genes pour les futurs points noirs : 
Tableau_plot$Number.of.CoreClock.Genes <- Tableau_CoreClock[,3]
# On enleve le spleen : 
Tableau_plot <- Tableau_plot[c(-4:-1),]

TissueSpe_Plot <- 
  ggplot(Tableau_plot,aes(x = Tissue, y = `Number of rhythmic genes`,fill = specificity,width=.6)) +
  geom_bar(position = "fill",stat = "identity") +
  geom_bar(stat="identity")+
  theme(legend.text = element_text(size=13,face="bold")) +
  coord_flip() +
  theme(axis.text.y = element_text(size=13,face="bold")) +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5, size = 12)) +
  scale_fill_manual(values=c(` tissue specific`="plum3", `2 tissues`="palegreen2", `3 tissues`="indianred1", `Common`="steelblue2")) +
  ggtitle("Distribution of the tissue-specificity\nof circadian genes in Mouse") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  geom_point(data=Tableau_plot,aes(Tissue,valeur.Ordonnee, size=Number.of.CoreClock.Genes)) +
  scale_size_continuous(range = c(-0.8,7)) 

head(Tableau_plot$`Number of rhythmic genes`)

# Fichier de sortie sans les Dots: 
tiff("~/Documents/RESULTS/Distribution_CircadGenes_Mouse.tiff", width = 10.4, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()
# Fichier de sortie avec les Dots: 
tiff("~/Documents/RESULTS/CoreClock_CircadGenes_Mouse.tiff", width = 11, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()













#### Pareil mais avec les paralogues :
setwd("Documents/DATA/CircaDB/MoGene_Stats")

Listtissues <- list.files(pattern="stats",full.names=FALSE) 
Listtissues_Bis <- gsub("_stats", "", Listtissues)
Listtissues_Bis_Title <- gsub("mogene_"," ", Listtissues_Bis)
Listtissues_Bis_legende <- gsub("_", " ", Listtissues_Bis_Title)
Nbtissue <- length(Listtissues)


List_count_totale <- matrix()
List_count_CoreClock <- matrix()
for(i in 1:Nbtissue)
{ 
  filename <- paste(Listtissues[i],sep="") 
  Mouse_Tissue_1_stats <- read.table(filename, header = F, quote = "",sep = ",", strip.white = TRUE, stringsAsFactors = FALSE)
  Mouse_Tissue_1_signif <- subset(Mouse_Tissue_1_stats,V10<0.05) 
  Mouse_Tissue_1 <- merge(Correspondance_ID, Mouse_Tissue_1_signif, by.x = "ID", by.y = "V1") 
  Mouse_Tissue_Circad <- as.data.frame(unique(Mouse_Tissue_1$Gene))
  Mouse_Tissue_Circad <- unique(Mouse_Tissue_Circad)
  colnames(Mouse_Tissue_Circad) <- "Gene"
  
  # On ne garde que les paralogues circadiens : 
  Mouse_Tissue_Circad_Paralog <- merge(Mouse_Tissue_Circad,Paralogous_List_unique, by="Gene")
  colnames(Mouse_Tissue_Circad_Paralog) <- "Circad.GeneID"
  
  # On Merge avec la count-liste :
  count_List_par_tissue <- merge(Mouse_Tissue_Circad_Paralog,Count_List, by="Circad.GeneID")
  # On change le nom de la colonne freq, car sinon il fait la somme apr??s : 
  colnames(count_List_par_tissue)[2] <- "frequency"
  # Et on fait la proportion de genes commun : 
  Freq_list_par_tissue <- count(count_List_par_tissue,'frequency' )
  colnames(Freq_list_par_tissue)[1] <- "Nb of tissues"
  colnames(Freq_list_par_tissue)[2] <- "Nb of CircadGenes" 
  # On considere si un gene au moins dans 4 tissus alors il est "commun": 
  count_tot_common_genes <-  sum(Freq_list_par_tissue$`Nb of CircadGenes`[4:14],na.rm = TRUE)
  Freq_list_par_tissue <- as.data.frame(Freq_list_par_tissue[1:4,2])
  colnames(Freq_list_par_tissue) <- Listtissues_Bis_legende[i]
  Freq_list_par_tissue[4,] <- count_tot_common_genes
  row.names(Freq_list_par_tissue) <- c(" tissue specific","2 tissues","3 tissues","Common")
  
  # On cr?? un fichier avec juste les "core clock genes" present dans ce tissus et avec le nombre de tissus tot dans lesquels on les retrouve,
  # et on lui fait le meme procede :
  CoreClock_genes_par_tissue <- merge(count_List_par_tissue,Mouse_BioSystem_Circadian_Genes, by.x = "Circad.GeneID", by.y = "Gene" )
  Freq_CoreClock_genes_par_tissue <- count(CoreClock_genes_par_tissue,'frequency' )
  colnames(Freq_CoreClock_genes_par_tissue)[1] <- "Nb of tissues"
  colnames(Freq_CoreClock_genes_par_tissue)[2] <- "Nb of CircadGenes" 
  value1 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(1,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  value2 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(2,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  value3 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(3,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  Tot_Commum <-  sum(Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(4:14,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)],na.rm = TRUE)
  # On fait le tableau final pour les core clock genes:
  a <- c(value1,value2,value3,Tot_Commum)
  row_names <- c(" tissue specific","2 tissues","3 tissues","Common")
  FreqFinal_CoreClock_tissue <- data.frame(a, row.names = row_names)
  FreqFinal_CoreClock_tissue[is.na(FreqFinal_CoreClock_tissue)] <- 0
  colnames(FreqFinal_CoreClock_tissue) <- Listtissues_Bis_legende[i]
  
  # On fait la liste totale : 
  List_count_totale <- cbind(List_count_totale,Freq_list_par_tissue)
  List_count_CoreClock <- cbind(List_count_CoreClock,FreqFinal_CoreClock_tissue)
}


List_count_totale <- List_count_totale[,-1]
Reverse_List_count_totale <- as.data.frame(t(List_count_totale))
# On ajoute une colonne avec le Nb total de genes circadiens par tissus  que l'on r??cup??re du fichier plus haut :
Reverse_List_count_totale <- cbind(Reverse_List_count_totale,Tot_CircadGenes_par_tissue)
# On range dans l'ordre en fonction du Nombre totale de genes circadiens : 
Reverse_List_count_totale_ordered <- Reverse_List_count_totale[order(Reverse_List_count_totale$Tot.CircadGenes,decreasing = FALSE),]
Reverse_List_count_totale_ordered <- Reverse_List_count_totale_ordered[,-5]
# On cre un tableau avec les futurs abscisses ou ordonnees des points :
Dots_Ordonnees <- matrix()
Dots_Ordonnees$tissueSpeOrdonnee <- Reverse_List_count_totale_ordered$`2 tissues` + Reverse_List_count_totale_ordered$`3 tissues` + 
  Reverse_List_count_totale_ordered$Common + Reverse_List_count_totale_ordered$` tissue specific`/2
Dots_Ordonnees$tissues2.Ordonnee <- Reverse_List_count_totale_ordered$`2 tissues`/2 + Reverse_List_count_totale_ordered$`3 tissues` + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$tissues3.Ordonnee <- Reverse_List_count_totale_ordered$`3 tissues`/2 + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$Common.Ordonnee <- Reverse_List_count_totale_ordered$Common/2
Dots_Ordonnees <- as.data.frame(Dots_Ordonnees)
Dots_Ordonnees <- Dots_Ordonnees[,-1]
rownames(Dots_Ordonnees) <- rownames(Reverse_List_count_totale_ordered)
Dots_Ordonnees <- as.data.frame(t(Dots_Ordonnees))
Dots_Ordonnees <- melt(cbind(Dots_Ordonnees, ind = rownames(Dots_Ordonnees)), id.vars = c('ind'))
colnames(Dots_Ordonnees) <- c("Ordonnes", "Tissue", "valeur.Ordonnee")
# Puis seulement maintenant on r??verse puis on transform par melt() :
Reverse_List_count_totale_ordered <- as.data.frame(t(Reverse_List_count_totale_ordered))
Tableau_plot <- melt(cbind(Reverse_List_count_totale_ordered, ind = rownames(Reverse_List_count_totale_ordered)), id.vars = c('ind'))
# On ajoute les valeur des ordonnees pour les futurs points noirs : 
Tableau_plot$valeur.Ordonnee <- Dots_Ordonnees[,3]
colnames(Tableau_plot) <- c("specificity","Tissue","Nb.of.genes","valeur.Ordonnee")

# On fait pareil avec la partie des core clock genes : 
List_count_CoreClock <- List_count_CoreClock[,-1]
Reverse_List_count_CoreClock <- as.data.frame(t(List_count_CoreClock))
Reverse_CoreClock_count_totale <- cbind(Reverse_List_count_CoreClock,Tot_CircadGenes_par_tissue)
Reverse_CoreClock_ordered <- Reverse_CoreClock_count_totale[order(Reverse_CoreClock_count_totale$Tot.CircadGenes,decreasing = FALSE),]
Reverse_CoreClock_ordered <- Reverse_CoreClock_ordered[,-5]
Reverse_CoreClock_ordered <- as.data.frame(t(Reverse_CoreClock_ordered))
Tableau_CoreClock <- melt(cbind(Reverse_CoreClock_ordered, ind = rownames(Reverse_CoreClock_ordered)), id.vars = c('ind'))

# On ajoute alors les valeurs des core clock genes pour les futurs points noirs : 
Tableau_plot$Number.of.CoreClock.Genes <- Tableau_CoreClock[,3]


plot <- 
  ggplot(Tableau_plot,aes(x = Tissue, y = Nb.of.genes,fill = specificity,width=.6)) +
  geom_bar(position = "fill",stat = "identity") +
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Dark2") +
  theme(legend.text = element_text(size=13,face="bold")) +
  coord_flip() +
  theme(axis.text.y = element_text(size=13,face="bold")) +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5, size = 12)) +
  ggtitle("Distribution of the number of circadian PARALOGS in Mouse") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  geom_point(data=Tableau_plot,aes(Tissue,valeur.Ordonnee, size=Number.of.CoreClock.Genes)) +
  scale_size_continuous(range = c(-0.8,7)) 


# Fichier de sortie : 
tiff("~/Documents/RESULTS/Distribution_Circad_Paralog_Mouse.tiff", width = 10.4, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()

# Fichier de sortie : 
tiff("~/Documents/RESULTS/CoreClock_ParalogCircad_Mouse.tiff", width = 11, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()










#### Pareil mais avec les paralogues from recent duplicates (Tina's List):
setwd("Documents/DATA/CircaDB/MoGene_Stats")

Listtissues <- list.files(pattern="stats",full.names=FALSE) 
Listtissues_Bis <- gsub("_stats", "", Listtissues)
Listtissues_Bis_Title <- gsub("mogene_"," ", Listtissues_Bis)
Listtissues_Bis_legende <- gsub("_", " ", Listtissues_Bis_Title)
Nbtissue <- length(Listtissues)

List_count_totale <- matrix()
List_count_CoreClock <- matrix()
for(i in 1:Nbtissue)
{ 
  filename <- paste(Listtissues[i],sep="") 
  Mouse_Tissue_1_stats <- read.table(filename, header = F, quote = "",sep = ",", strip.white = TRUE, stringsAsFactors = FALSE)
  Mouse_Tissue_1_signif <- subset(Mouse_Tissue_1_stats,V10<0.05) 
  Mouse_Tissue_1 <- merge(Correspondance_ID, Mouse_Tissue_1_signif, by.x = "ID", by.y = "V1") 
  Mouse_Tissue_Circad <- as.data.frame(unique(Mouse_Tissue_1$Gene))
  Mouse_Tissue_Circad <- unique(Mouse_Tissue_Circad)
  colnames(Mouse_Tissue_Circad) <- "Gene"
  
  # On ne garde que les paralogues circadiens : 
  Mouse_Tissue_Circad_Paralog <- merge(Mouse_Tissue_Circad,Paralog_Mouse_Recent_Duplicates_ListUnique, by="Gene")
  colnames(Mouse_Tissue_Circad_Paralog) <- "Circad.GeneID"
  
  # On Merge avec la count-liste :
  count_List_par_tissue <- merge(Mouse_Tissue_Circad_Paralog,Count_List, by="Circad.GeneID")
  # On change le nom de la colonne freq, car sinon il fait la somme apr??s : 
  colnames(count_List_par_tissue)[2] <- "frequency"
  # Et on fait la proportion de genes commun : 
  Freq_list_par_tissue <- count(count_List_par_tissue,'frequency' )
  colnames(Freq_list_par_tissue)[1] <- "Nb of tissues"
  colnames(Freq_list_par_tissue)[2] <- "Nb of CircadGenes" 
  # On considere si un gene au moins dans 4 tissus alors il est "commun": 
  count_tot_common_genes <-  sum(Freq_list_par_tissue$`Nb of CircadGenes`[4:14],na.rm = TRUE)
  Freq_list_par_tissue <- as.data.frame(Freq_list_par_tissue[1:4,2])
  colnames(Freq_list_par_tissue) <- Listtissues_Bis_legende[i]
  Freq_list_par_tissue[4,] <- count_tot_common_genes
  row.names(Freq_list_par_tissue) <- c(" tissue specific","2 tissues","3 tissues","Common")
  
  # On cr?? un fichier avec juste les "core clock genes" present dans ce tissus et avec le nombre de tissus tot dans lesquels on les retrouve,
  # et on lui fait le meme procede :
  CoreClock_genes_par_tissue <- merge(count_List_par_tissue,Mouse_BioSystem_Circadian_Genes, by.x = "Circad.GeneID", by.y = "Gene" )
  Freq_CoreClock_genes_par_tissue <- count(CoreClock_genes_par_tissue,'frequency' )
  colnames(Freq_CoreClock_genes_par_tissue)[1] <- "Nb of tissues"
  colnames(Freq_CoreClock_genes_par_tissue)[2] <- "Nb of CircadGenes" 
  value1 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(1,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  value2 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(2,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  value3 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(3,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  Tot_Commum <-  sum(Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(4:14,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)],na.rm = TRUE)
  # On fait le tableau final pour les core clock genes:
  a <- c(value1,value2,value3,Tot_Commum)
  row_names <- c(" tissue specific","2 tissues","3 tissues","Common")
  FreqFinal_CoreClock_tissue <- data.frame(a, row.names = row_names)
  FreqFinal_CoreClock_tissue[is.na(FreqFinal_CoreClock_tissue)] <- 0
  colnames(FreqFinal_CoreClock_tissue) <- Listtissues_Bis_legende[i]
  
  # On fait la liste totale : 
  List_count_totale <- cbind(List_count_totale,Freq_list_par_tissue)
  List_count_CoreClock <- cbind(List_count_CoreClock,FreqFinal_CoreClock_tissue)
}


List_count_totale <- List_count_totale[,-1]
Reverse_List_count_totale <- as.data.frame(t(List_count_totale))
# On ajoute une colonne avec le Nb total de genes circadiens par tissus  que l'on r??cup??re du fichier plus haut :
Reverse_List_count_totale <- cbind(Reverse_List_count_totale,Tot_CircadGenes_par_tissue)
# On range dans l'ordre en fonction du Nombre totale de genes circadiens : 
Reverse_List_count_totale_ordered <- Reverse_List_count_totale[order(Reverse_List_count_totale$Tot.CircadGenes,decreasing = FALSE),]
Reverse_List_count_totale_ordered <- Reverse_List_count_totale_ordered[,-5]
# On cre un tableau avec les futurs abscisses ou ordonnees des points :
Dots_Ordonnees <- matrix()
Dots_Ordonnees$tissueSpeOrdonnee <- Reverse_List_count_totale_ordered$`2 tissues` + Reverse_List_count_totale_ordered$`3 tissues` + 
  Reverse_List_count_totale_ordered$Common + Reverse_List_count_totale_ordered$` tissue specific`/2
Dots_Ordonnees$tissues2.Ordonnee <- Reverse_List_count_totale_ordered$`2 tissues`/2 + Reverse_List_count_totale_ordered$`3 tissues` + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$tissues3.Ordonnee <- Reverse_List_count_totale_ordered$`3 tissues`/2 + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$Common.Ordonnee <- Reverse_List_count_totale_ordered$Common/2
Dots_Ordonnees <- as.data.frame(Dots_Ordonnees)
Dots_Ordonnees <- Dots_Ordonnees[,-1]
rownames(Dots_Ordonnees) <- rownames(Reverse_List_count_totale_ordered)
Dots_Ordonnees <- as.data.frame(t(Dots_Ordonnees))
Dots_Ordonnees <- melt(cbind(Dots_Ordonnees, ind = rownames(Dots_Ordonnees)), id.vars = c('ind'))
colnames(Dots_Ordonnees) <- c("Ordonnes", "Tissue", "valeur.Ordonnee")
# Puis seulement maintenant on r??verse puis on transform par melt() :
Reverse_List_count_totale_ordered <- as.data.frame(t(Reverse_List_count_totale_ordered))
Tableau_plot <- melt(cbind(Reverse_List_count_totale_ordered, ind = rownames(Reverse_List_count_totale_ordered)), id.vars = c('ind'))
# On ajoute les valeur des ordonnees pour les futurs points noirs : 
Tableau_plot$valeur.Ordonnee <- Dots_Ordonnees[,3]
colnames(Tableau_plot) <- c("specificity","Tissue","Nb.of.genes","valeur.Ordonnee")

# On fait pareil avec la partie des core clock genes : 
List_count_CoreClock <- List_count_CoreClock[,-1]
Reverse_List_count_CoreClock <- as.data.frame(t(List_count_CoreClock))
Reverse_CoreClock_count_totale <- cbind(Reverse_List_count_CoreClock,Tot_CircadGenes_par_tissue)
Reverse_CoreClock_ordered <- Reverse_CoreClock_count_totale[order(Reverse_CoreClock_count_totale$Tot.CircadGenes,decreasing = FALSE),]
Reverse_CoreClock_ordered <- Reverse_CoreClock_ordered[,-5]
Reverse_CoreClock_ordered <- as.data.frame(t(Reverse_CoreClock_ordered))
Tableau_CoreClock <- melt(cbind(Reverse_CoreClock_ordered, ind = rownames(Reverse_CoreClock_ordered)), id.vars = c('ind'))

# On ajoute alors les valeurs des core clock genes pour les futurs points noirs : 
Tableau_plot$Number.of.CoreClock.Genes <- Tableau_CoreClock[,3]


plot <- 
  ggplot(Tableau_plot,aes(x = Tissue, y = Nb.of.genes,fill = specificity,width=.6)) +
  geom_bar(position = "fill",stat = "identity") +
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Set3") +
  theme(legend.text = element_text(size=13,face="bold")) +
  coord_flip() +
  theme(axis.text.y = element_text(size=13,face="bold")) +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5, size = 12)) +
  ggtitle("Distribution of the number of circadian PARALOGS \n from RECENT DUPLICATION in Mouse (Tina's List)") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  geom_point(data=Tableau_plot,aes(Tissue,valeur.Ordonnee, size=Number.of.CoreClock.Genes)) +
  scale_size_continuous(range = c(-2,0)) 

# Fichier de sortie : 
tiff("~/Documents/RESULTS/Distribution_Circad_Paralog_RecentDuplic_Mouse.tiff", width = 10.4, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()

# Fichier de sortie : 
tiff("~/Documents/RESULTS/CoreClock_RecentDuplicCircad_Mouse.tiff", width = 11, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()














#### Pareil mais avec les paralogues from RECENT duplicates (Jialin's List):
setwd("Documents/DATA/CircaDB/MoGene_Stats")

Listtissues <- list.files(pattern="stats",full.names=FALSE) 
Listtissues_Bis <- gsub("_stats", "", Listtissues)
Listtissues_Bis_Title <- gsub("mogene_"," ", Listtissues_Bis)
Listtissues_Bis_legende <- gsub("_", " ", Listtissues_Bis_Title)
Nbtissue <- length(Listtissues)

List_count_totale <- matrix()
List_count_CoreClock <- matrix()
for(i in 1:Nbtissue)
{ 
  filename <- paste(Listtissues[i],sep="") 
  Mouse_Tissue_1_stats <- read.table(filename, header = F, quote = "",sep = ",", strip.white = TRUE, stringsAsFactors = FALSE)
  Mouse_Tissue_1_signif <- subset(Mouse_Tissue_1_stats,V10<0.05) 
  Mouse_Tissue_1 <- merge(Correspondance_ID, Mouse_Tissue_1_signif, by.x = "ID", by.y = "V1") 
  Mouse_Tissue_Circad <- as.data.frame(unique(Mouse_Tissue_1$Gene))
  Mouse_Tissue_Circad <- unique(Mouse_Tissue_Circad)
  colnames(Mouse_Tissue_Circad) <- "Gene"
  
  # On ne garde que les paralogues circadiens : 
  Mouse_Tissue_Circad_Paralog <- merge(Mouse_Tissue_Circad,Paralog_RecentDuplic_Jialin_MurinaeSciurognathi, by="Gene")
  colnames(Mouse_Tissue_Circad_Paralog) <- "Circad.GeneID"
  
  # On Merge avec la count-liste :
  count_List_par_tissue <- merge(Mouse_Tissue_Circad_Paralog,Count_List, by="Circad.GeneID")
  # On change le nom de la colonne freq, car sinon il fait la somme apr??s : 
  colnames(count_List_par_tissue)[2] <- "frequency"
  # Et on fait la proportion de genes commun : 
  Freq_list_par_tissue <- count(count_List_par_tissue,'frequency' )
  colnames(Freq_list_par_tissue)[1] <- "Nb of tissues"
  colnames(Freq_list_par_tissue)[2] <- "Nb of CircadGenes" 
  # On considere si un gene au moins dans 4 tissus alors il est "commun": 
  count_tot_common_genes <-  sum(Freq_list_par_tissue$`Nb of CircadGenes`[4:14],na.rm = TRUE)
  Freq_list_par_tissue <- as.data.frame(Freq_list_par_tissue[1:4,2])
  colnames(Freq_list_par_tissue) <- Listtissues_Bis_legende[i]
  Freq_list_par_tissue[4,] <- count_tot_common_genes
  row.names(Freq_list_par_tissue) <- c(" tissue specific","2 tissues","3 tissues","Common")
  
  # On cr?? un fichier avec juste les "core clock genes" present dans ce tissus et avec le nombre de tissus tot dans lesquels on les retrouve,
  # et on lui fait le meme procede :
  CoreClock_genes_par_tissue <- merge(count_List_par_tissue,Mouse_BioSystem_Circadian_Genes, by.x = "Circad.GeneID", by.y = "Gene" )
  Freq_CoreClock_genes_par_tissue <- count(CoreClock_genes_par_tissue,'frequency' )
  colnames(Freq_CoreClock_genes_par_tissue)[1] <- "Nb of tissues"
  colnames(Freq_CoreClock_genes_par_tissue)[2] <- "Nb of CircadGenes" 
  value1 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(1,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  value2 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(2,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  value3 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(3,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  Tot_Commum <-  sum(Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(4:14,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)],na.rm = TRUE)
  # On fait le tableau final pour les core clock genes:
  a <- c(value1,value2,value3,Tot_Commum)
  row_names <- c(" tissue specific","2 tissues","3 tissues","Common")
  FreqFinal_CoreClock_tissue <- data.frame(a, row.names = row_names)
  FreqFinal_CoreClock_tissue[is.na(FreqFinal_CoreClock_tissue)] <- 0
  colnames(FreqFinal_CoreClock_tissue) <- Listtissues_Bis_legende[i]
  
  # On fait la liste totale : 
  List_count_totale <- cbind(List_count_totale,Freq_list_par_tissue)
  List_count_CoreClock <- cbind(List_count_CoreClock,FreqFinal_CoreClock_tissue)
}


List_count_totale <- List_count_totale[,-1]
Reverse_List_count_totale <- as.data.frame(t(List_count_totale))
# On ajoute une colonne avec le Nb total de genes circadiens par tissus  que l'on r??cup??re du fichier plus haut :
Reverse_List_count_totale <- cbind(Reverse_List_count_totale,Tot_CircadGenes_par_tissue)
# On range dans l'ordre en fonction du Nombre totale de genes circadiens : 
Reverse_List_count_totale_ordered <- Reverse_List_count_totale[order(Reverse_List_count_totale$Tot.CircadGenes,decreasing = FALSE),]
Reverse_List_count_totale_ordered <- Reverse_List_count_totale_ordered[,-5]
# On cre un tableau avec les futurs abscisses ou ordonnees des points :
Dots_Ordonnees <- matrix()
Dots_Ordonnees$tissueSpeOrdonnee <- Reverse_List_count_totale_ordered$`2 tissues` + Reverse_List_count_totale_ordered$`3 tissues` + 
  Reverse_List_count_totale_ordered$Common + Reverse_List_count_totale_ordered$` tissue specific`/2
Dots_Ordonnees$tissues2.Ordonnee <- Reverse_List_count_totale_ordered$`2 tissues`/2 + Reverse_List_count_totale_ordered$`3 tissues` + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$tissues3.Ordonnee <- Reverse_List_count_totale_ordered$`3 tissues`/2 + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$Common.Ordonnee <- Reverse_List_count_totale_ordered$Common/2
Dots_Ordonnees <- as.data.frame(Dots_Ordonnees)
Dots_Ordonnees <- Dots_Ordonnees[,-1]
rownames(Dots_Ordonnees) <- rownames(Reverse_List_count_totale_ordered)
Dots_Ordonnees <- as.data.frame(t(Dots_Ordonnees))
Dots_Ordonnees <- melt(cbind(Dots_Ordonnees, ind = rownames(Dots_Ordonnees)), id.vars = c('ind'))
colnames(Dots_Ordonnees) <- c("Ordonnes", "Tissue", "valeur.Ordonnee")
# Puis seulement maintenant on r??verse puis on transform par melt() :
Reverse_List_count_totale_ordered <- as.data.frame(t(Reverse_List_count_totale_ordered))
Tableau_plot <- melt(cbind(Reverse_List_count_totale_ordered, ind = rownames(Reverse_List_count_totale_ordered)), id.vars = c('ind'))
# On ajoute les valeur des ordonnees pour les futurs points noirs : 
Tableau_plot$valeur.Ordonnee <- Dots_Ordonnees[,3]
colnames(Tableau_plot) <- c("specificity","Tissue","Nb.of.genes","valeur.Ordonnee")

# On fait pareil avec la partie des core clock genes : 
List_count_CoreClock <- List_count_CoreClock[,-1]
Reverse_List_count_CoreClock <- as.data.frame(t(List_count_CoreClock))
Reverse_CoreClock_count_totale <- cbind(Reverse_List_count_CoreClock,Tot_CircadGenes_par_tissue)
Reverse_CoreClock_ordered <- Reverse_CoreClock_count_totale[order(Reverse_CoreClock_count_totale$Tot.CircadGenes,decreasing = FALSE),]
Reverse_CoreClock_ordered <- Reverse_CoreClock_ordered[,-5]
Reverse_CoreClock_ordered <- as.data.frame(t(Reverse_CoreClock_ordered))
Tableau_CoreClock <- melt(cbind(Reverse_CoreClock_ordered, ind = rownames(Reverse_CoreClock_ordered)), id.vars = c('ind'))

# On ajoute alors les valeurs des core clock genes pour les futurs points noirs : 
Tableau_plot$Number.of.CoreClock.Genes <- Tableau_CoreClock[,3]

plot <- 
  ggplot(Tableau_plot,aes(x = Tissue, y = Nb.of.genes,fill = specificity,width=.6)) +
  geom_bar(position = "fill",stat = "identity") +
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Set3") +
  theme(legend.text = element_text(size=13,face="bold")) +
  coord_flip() +
  theme(axis.text.y = element_text(size=13,face="bold")) +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5, size = 12)) +
  ggtitle("Distribution of the number of circadian PARALOGS \n from RECENT DUPLICATION (Murinae & Sciurognathi) in Mouse") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  geom_point(data=Tableau_plot,aes(Tissue,valeur.Ordonnee, size=Number.of.CoreClock.Genes)) +
  scale_size_continuous(range = c(-2,0)) 

# Fichier de sortie : 
tiff("~/Documents/RESULTS/Distribution_Circad_Paralog_RecentDuplic_MurinaeSciurognathi_Mouse.tiff", width = 10.4, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()

# Fichier de sortie : 
tiff("~/Documents/RESULTS/CoreClock_RecentDuplicCircad_MurinaeSciurognathi_Mouse.tiff", width = 11, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()















#### Pareil mais avec les paralogues from OLD duplicates (Chordates_Vertebrate) (Jialin's List):
setwd("Documents/DATA/CircaDB/MoGene_Stats")

Listtissues <- list.files(pattern="stats",full.names=FALSE) 
Listtissues_Bis <- gsub("_stats", "", Listtissues)
Listtissues_Bis_Title <- gsub("mogene_"," ", Listtissues_Bis)
Listtissues_Bis_legende <- gsub("_", " ", Listtissues_Bis_Title)
Nbtissue <- length(Listtissues)

List_count_totale <- matrix()
List_count_CoreClock <- matrix()
for(i in 1:Nbtissue)
{ 
  filename <- paste(Listtissues[i],sep="") 
  Mouse_Tissue_1_stats <- read.table(filename, header = F, quote = "",sep = ",", strip.white = TRUE, stringsAsFactors = FALSE)
  Mouse_Tissue_1_signif <- subset(Mouse_Tissue_1_stats,V10<0.05) 
  Mouse_Tissue_1 <- merge(Correspondance_ID, Mouse_Tissue_1_signif, by.x = "ID", by.y = "V1") 
  Mouse_Tissue_Circad <- as.data.frame(unique(Mouse_Tissue_1$Gene))
  Mouse_Tissue_Circad <- unique(Mouse_Tissue_Circad)
  colnames(Mouse_Tissue_Circad) <- "Gene"
  
  # On ne garde que les paralogues circadiens : 
  Mouse_Tissue_Circad_Paralog <- merge(Mouse_Tissue_Circad,Paralog_OldDuplic_UniqList, by="Gene")
  colnames(Mouse_Tissue_Circad_Paralog) <- "Circad.GeneID"
  
  # On Merge avec la count-liste :
  count_List_par_tissue <- merge(Mouse_Tissue_Circad_Paralog,Count_List, by="Circad.GeneID")
  # On change le nom de la colonne freq, car sinon il fait la somme apr??s : 
  colnames(count_List_par_tissue)[2] <- "frequency"
  # Et on fait la proportion de genes commun : 
  Freq_list_par_tissue <- count(count_List_par_tissue,'frequency' )
  colnames(Freq_list_par_tissue)[1] <- "Nb of tissues"
  colnames(Freq_list_par_tissue)[2] <- "Nb of CircadGenes" 
  # On considere si un gene au moins dans 4 tissus alors il est "commun": 
  count_tot_common_genes <-  sum(Freq_list_par_tissue$`Nb of CircadGenes`[4:14],na.rm = TRUE)
  Freq_list_par_tissue <- as.data.frame(Freq_list_par_tissue[1:4,2])
  colnames(Freq_list_par_tissue) <- Listtissues_Bis_legende[i]
  Freq_list_par_tissue[4,] <- count_tot_common_genes
  row.names(Freq_list_par_tissue) <- c(" tissue specific","2 tissues","3 tissues","Common")
  
  # On cr?? un fichier avec juste les "core clock genes" present dans ce tissus et avec le nombre de tissus tot dans lesquels on les retrouve,
  # et on lui fait le meme procede :
  CoreClock_genes_par_tissue <- merge(count_List_par_tissue,Mouse_BioSystem_Circadian_Genes, by.x = "Circad.GeneID", by.y = "Gene" )
  Freq_CoreClock_genes_par_tissue <- count(CoreClock_genes_par_tissue,'frequency' )
  colnames(Freq_CoreClock_genes_par_tissue)[1] <- "Nb of tissues"
  colnames(Freq_CoreClock_genes_par_tissue)[2] <- "Nb of CircadGenes" 
  value1 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(1,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  value2 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(2,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  value3 <- Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(3,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)]
  Tot_Commum <-  sum(Freq_CoreClock_genes_par_tissue$`Nb of CircadGenes`[match(4:14,Freq_CoreClock_genes_par_tissue$`Nb of tissues`)],na.rm = TRUE)
  # On fait le tableau final pour les core clock genes:
  a <- c(value1,value2,value3,Tot_Commum)
  row_names <- c(" tissue specific","2 tissues","3 tissues","Common")
  FreqFinal_CoreClock_tissue <- data.frame(a, row.names = row_names)
  FreqFinal_CoreClock_tissue[is.na(FreqFinal_CoreClock_tissue)] <- 0
  colnames(FreqFinal_CoreClock_tissue) <- Listtissues_Bis_legende[i]
  
  # On fait la liste totale : 
  List_count_totale <- cbind(List_count_totale,Freq_list_par_tissue)
  List_count_CoreClock <- cbind(List_count_CoreClock,FreqFinal_CoreClock_tissue)
}


List_count_totale <- List_count_totale[,-1]
Reverse_List_count_totale <- as.data.frame(t(List_count_totale))
# On ajoute une colonne avec le Nb total de genes circadiens par tissus  que l'on r??cup??re du fichier plus haut :
Reverse_List_count_totale <- cbind(Reverse_List_count_totale,Tot_CircadGenes_par_tissue)
# On range dans l'ordre en fonction du Nombre totale de genes circadiens : 
Reverse_List_count_totale_ordered <- Reverse_List_count_totale[order(Reverse_List_count_totale$Tot.CircadGenes,decreasing = FALSE),]
Reverse_List_count_totale_ordered <- Reverse_List_count_totale_ordered[,-5]
# On cre un tableau avec les futurs abscisses ou ordonnees des points :
Dots_Ordonnees <- matrix()
Dots_Ordonnees$tissueSpeOrdonnee <- Reverse_List_count_totale_ordered$`2 tissues` + Reverse_List_count_totale_ordered$`3 tissues` + 
  Reverse_List_count_totale_ordered$Common + Reverse_List_count_totale_ordered$` tissue specific`/2
Dots_Ordonnees$tissues2.Ordonnee <- Reverse_List_count_totale_ordered$`2 tissues`/2 + Reverse_List_count_totale_ordered$`3 tissues` + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$tissues3.Ordonnee <- Reverse_List_count_totale_ordered$`3 tissues`/2 + Reverse_List_count_totale_ordered$Common
Dots_Ordonnees$Common.Ordonnee <- Reverse_List_count_totale_ordered$Common/2
Dots_Ordonnees <- as.data.frame(Dots_Ordonnees)
Dots_Ordonnees <- Dots_Ordonnees[,-1]
rownames(Dots_Ordonnees) <- rownames(Reverse_List_count_totale_ordered)
Dots_Ordonnees <- as.data.frame(t(Dots_Ordonnees))
Dots_Ordonnees <- melt(cbind(Dots_Ordonnees, ind = rownames(Dots_Ordonnees)), id.vars = c('ind'))
colnames(Dots_Ordonnees) <- c("Ordonnes", "Tissue", "valeur.Ordonnee")
# Puis seulement maintenant on r??verse puis on transform par melt() :
Reverse_List_count_totale_ordered <- as.data.frame(t(Reverse_List_count_totale_ordered))
Tableau_plot <- melt(cbind(Reverse_List_count_totale_ordered, ind = rownames(Reverse_List_count_totale_ordered)), id.vars = c('ind'))
# On ajoute les valeur des ordonnees pour les futurs points noirs : 
Tableau_plot$valeur.Ordonnee <- Dots_Ordonnees[,3]
colnames(Tableau_plot) <- c("specificity","Tissue","Nb.of.genes","valeur.Ordonnee")

# On fait pareil avec la partie des core clock genes : 
List_count_CoreClock <- List_count_CoreClock[,-1]
Reverse_List_count_CoreClock <- as.data.frame(t(List_count_CoreClock))
Reverse_CoreClock_count_totale <- cbind(Reverse_List_count_CoreClock,Tot_CircadGenes_par_tissue)
Reverse_CoreClock_ordered <- Reverse_CoreClock_count_totale[order(Reverse_CoreClock_count_totale$Tot.CircadGenes,decreasing = FALSE),]
Reverse_CoreClock_ordered <- Reverse_CoreClock_ordered[,-5]
Reverse_CoreClock_ordered <- as.data.frame(t(Reverse_CoreClock_ordered))
Tableau_CoreClock <- melt(cbind(Reverse_CoreClock_ordered, ind = rownames(Reverse_CoreClock_ordered)), id.vars = c('ind'))

# On ajoute alors les valeurs des core clock genes pour les futurs points noirs : 
Tableau_plot$Number.of.CoreClock.Genes <- Tableau_CoreClock[,3]

plot <- 
  ggplot(Tableau_plot,aes(x = Tissue, y = Nb.of.genes,fill = specificity,width=.6)) +
  geom_bar(position = "fill",stat = "identity") +
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Set3") +
  theme(legend.text = element_text(size=13,face="bold")) +
  coord_flip() +
  theme(axis.text.y = element_text(size=13,face="bold")) +
  theme(axis.text.x=element_text(angle=0,hjust=1,vjust=0.5, size = 12)) +
  ggtitle("Distribution of the number of circadian PARALOGS \n from OLD DUPLICATION (Chordates_Vertebrate) in Mouse") + 
  theme(plot.title = element_text(lineheight=.8, face="bold")) +
  geom_point(data=Tableau_plot,aes(Tissue,valeur.Ordonnee, size=Number.of.CoreClock.Genes)) +
  scale_size_continuous(range = c(-1,1),breaks = 1) 


# Fichier de sortie : 
tiff("~/Documents/RESULTS/Distribution_Circad_Paralog_OldDuplic_ChordatesVertebrate_Mouse.tiff", width = 10.4, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()

# Fichier de sortie : 
tiff("~/Documents/RESULTS/CoreClock_OldDuplicCircad_Mouse.tiff", width = 11, height = 7, units = 'in', res = 300)
multiplot(plot, cols=1)
dev.off()
