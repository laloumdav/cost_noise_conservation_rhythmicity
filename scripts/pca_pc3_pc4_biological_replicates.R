library(ggplot2)
library(ellipse)
library(ggrepel) 
library(cowplot)
library(officer)
library(magrittr)
### ### ### ### ### ### ### ### ### 
### PCA on biological replicates ##   PC3 & PC4
### ### ### ### ### ### ### ### ### 
list.pca.data <- list()
list.elypses.data <- list()


#1.1# ARABIDOPSIS leaves transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/transcriptome/transcript_level.txt", h=T, fill=T)
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("LD", rep(c("00", "04", "08", "12", "16", "20"), each=3), sep="")
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.data <- na.omit(replicates.data)

## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4])
# Make elypses
centroids <- aggregate(cbind(X,Y)~Sample, pca.data, mean)
row.names(centroids) <- centroids$Sample
elypses.data  <- do.call(rbind,lapply(unique(pca.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(pca.data[pca.data$Sample==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.5),
             stringsAsFactors=FALSE)))

## redo PCA without elypses ##
replicates.data <- replicates.data.raw[, -1]
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.timepoint <- na.omit(replicates.timepoint)
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:3, 6), sep="")

# Store in lists
list.pca.data[[1]] <- cbind(pca.data, data.frame(species=rep("arabidopsis\nLeaves", nrow(pca.data)),
                                                 data=rep("transcripts", nrow(pca.data))))
list.elypses.data[[1]] <- cbind(elypses.data, data.frame(species=rep("arabidopsis\nLeaves", nrow(elypses.data)),
                                                         data=rep("transcripts", nrow(elypses.data))))





#1.2# ARABIDOPSIS leaves proteome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/arabidopsis/leaves/proteome/protein_level.txt", h=T, fill=T)
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("LL", rep(c("12", "16", "20", "00", "04", "08"), 5), sep="")
#length(timepoint.names)
#length(colnames(replicates.data))
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

replicates.data <- na.omit(replicates.data)
## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4])
# Make elypses
centroids <- aggregate(cbind(X,Y)~Sample, pca.data, mean)
row.names(centroids) <- centroids$Sample
elypses.data  <- do.call(rbind,lapply(unique(pca.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(pca.data[pca.data$Sample==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.55),
             stringsAsFactors=FALSE)))

## redo PCA without elypses ##
replicates.data <- replicates.data.raw[, -1]
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

replicates.data <- na.omit(replicates.data)
## PCA ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:5, each=6), sep="")

# Store in lists
list.pca.data[[2]] <- cbind(pca.data, data.frame(species=rep("arabidopsis\nLeaves", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[2]] <- cbind(elypses.data, data.frame(species=rep("arabidopsis\nLeaves", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))






#2.1# CYANOBACTERIA transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/cyanobacteria/transcriptome/transcript_level.txt", h=T, fill=T)
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("LL", rep(c("04", "08", "12", "16", "20", "00"), 4), sep="")[1:ncol(replicates.data)]
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.data <- na.omit(replicates.data)

## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4])
# Make elypses
centroids <- aggregate(cbind(X,Y)~Sample, pca.data, mean)
row.names(centroids) <- centroids$Sample
elypses.data  <- do.call(rbind,lapply(unique(pca.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(pca.data[pca.data$Sample==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.5),
             stringsAsFactors=FALSE)))

## redo PCA without elypses ##
replicates.data <- replicates.data.raw[, -1]
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.timepoint <- na.omit(replicates.timepoint)
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:4, each=6), sep="")[1:ncol(replicates.data)]

# Store in lists
list.pca.data[[3]] <- cbind(pca.data, data.frame(species=rep("cyanobacteria", nrow(pca.data)),
                                                 data=rep("transcripts", nrow(pca.data))))
list.elypses.data[[3]] <- cbind(elypses.data, data.frame(species=rep("cyanobacteria", nrow(elypses.data)),
                                                         data=rep("transcripts", nrow(elypses.data))))




#2.2# CYANOBACTERIA proteome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/cyanobacteria/proteome/protein_level.txt", h=T, fill=T)
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("ZT", rep(c("06.5", "09.5", "11.5", "12.5", "14.5", "17.5", "21.5", "23.5", "00.5", "02.5"), 2), sep="")
#length(timepoint.names)
#length(colnames(replicates.data))
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.data <- na.omit(replicates.data)
## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4])
# Make elypses
centroids <- aggregate(cbind(X,Y)~Sample, pca.data, mean)
row.names(centroids) <- centroids$Sample
elypses.data  <- do.call(rbind,lapply(unique(pca.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(pca.data[pca.data$Sample==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.2),
             stringsAsFactors=FALSE)))

## redo PCA without elypses ##
replicates.data <- replicates.data.raw[, -1]
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.data <- na.omit(replicates.data)
## PCA ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:2, each=10), sep="")

# Store in lists
list.pca.data[[4]] <- cbind(pca.data, data.frame(species=rep("cyanobacteria", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[4]] <- cbind(elypses.data, data.frame(species=rep("cyanobacteria", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))




#3.1# OSTEOCOCCUS transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/ostreococcus/transcriptome/transcript_level.txt", h=T, fill=T)
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("LD", rep(c("03", "06", "09", "12", "15", "18", "21", "24"), each=3), sep="")
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

replicates.data <- na.omit(replicates.data)
## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4])
# Make elypses
centroids <- aggregate(cbind(X,Y)~Sample, pca.data, mean)
row.names(centroids) <- centroids$Sample
elypses.data  <- do.call(rbind,lapply(unique(pca.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(pca.data[pca.data$Sample==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.5),
             stringsAsFactors=FALSE)))

## redo PCA without elypses ##
replicates.data <- replicates.data.raw[, -1]
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

replicates.data <- na.omit(replicates.data)
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:3, 8), sep="")

# Store in lists
list.pca.data[[5]] <- cbind(pca.data, data.frame(species=rep("ostreococcus", nrow(pca.data)),
                                                 data=rep("transcripts", nrow(pca.data))))
list.elypses.data[[5]] <- cbind(elypses.data, data.frame(species=rep("ostreococcus", nrow(elypses.data)),
                                                         data=rep("transcripts", nrow(elypses.data))))





#3.2# OSTEOCOCCUS proteome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/ostreococcus/proteome/protein_level.txt", h=T, fill=T)
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("ZT", rep(c("00", "04", "08", "12", "16", "20"), 5), sep="")
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

replicates.data <- na.omit(replicates.data)
## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4])
# Make elypses
centroids <- aggregate(cbind(X,Y)~Sample, pca.data, mean)
row.names(centroids) <- centroids$Sample
elypses.data  <- do.call(rbind,lapply(unique(pca.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(pca.data[pca.data$Sample==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.6),
             stringsAsFactors=FALSE)))

## redo PCA without elypses ##
replicates.data <- replicates.data.raw[, -1]
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

replicates.data <- na.omit(replicates.data)
## PCA ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:5, each=6), sep="")

# Store in lists
list.pca.data[[6]] <- cbind(pca.data, data.frame(species=rep("ostreococcus", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[6]] <- cbind(elypses.data, data.frame(species=rep("ostreococcus", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))






#4.1# MOUSE liver transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/transcriptome/transcript_level.txt", 
                                  head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
replicates.data <- replicates.data.raw[, 1:4*-1]
timepoint.names <- paste("ZT", rep(c("00", "02", "04", "06", "08", "10", "12", "14", "16", "18", "20", "22"), 2), sep="")
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.data <- na.omit(replicates.data)
## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4])
# Make elypses
centroids <- aggregate(cbind(X,Y)~Sample, pca.data, mean)
row.names(centroids) <- centroids$Sample
elypses.data  <- do.call(rbind,lapply(unique(pca.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(pca.data[pca.data$Sample==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.2),
             stringsAsFactors=FALSE)))

## redo PCA without elypses ##
replicates.data <- replicates.data.raw[, 1:4*-1]
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.data <- na.omit(replicates.data)
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:2, each=12), sep="")

# Store in lists
list.pca.data[[7]] <- cbind(pca.data, data.frame(species=rep("mouse\nLiver", nrow(pca.data)),
                                                 data=rep("transcripts", nrow(pca.data))))
list.elypses.data[[7]] <- cbind(elypses.data, data.frame(species=rep("mouse\nLiver", nrow(elypses.data)),
                                                         data=rep("transcripts", nrow(elypses.data))))





#4.2# Mouse liver proteome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/proteome/protein_level.txt", head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
replicates.data <- replicates.data.raw[, c(-1,-2)]
timepoint.names <- paste("ZT", rep(c("00", "03", "06", "09", "12", "15", "18", "21"), 2), sep="")
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

replicates.data <- na.omit(replicates.data)
## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4])
# Make elypses
centroids <- aggregate(cbind(X,Y)~Sample, pca.data, mean)
row.names(centroids) <- centroids$Sample
elypses.data  <- do.call(rbind,lapply(unique(pca.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(pca.data[pca.data$Sample==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.2),
             stringsAsFactors=FALSE)))

## redo PCA without elypses ##
replicates.data <- replicates.data.raw[, c(-1,-2)]
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

replicates.data <- na.omit(replicates.data)
## PCA ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:2, each=8), sep="")

# Store in lists
list.pca.data[[8]] <- cbind(pca.data, data.frame(species=rep("mouse\nLiver", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[8]] <- cbind(elypses.data, data.frame(species=rep("mouse\nLiver", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))







### store TOT data and make all PCA plots ###
TOT.pca.data <- data.frame(Sample=NA, 
                           X=NA,
                           Y=NA,
                           var.per=NA,
                           time.point=NA,
                           species=NA, 
                           data=NA)
TOT.elypses.data <- data.frame(Sample=NA, 
                               X=NA,
                               Y=NA,
                               species=NA, 
                               data=NA)
for (i in 1:length(list.pca.data)) {
  TOT.pca.data <- rbind(TOT.pca.data, list.pca.data[[i]])
  TOT.elypses.data <- rbind(TOT.elypses.data, list.elypses.data[[i]])
}
TOT.pca.data <- TOT.pca.data[-1, ]
TOT.elypses.data <- TOT.elypses.data[-1, ]

# Function to adjust overall legend size
addSmallLegend <- function(myPlot, pointSize = 2, textSize = 8, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

ggplotList <- list()
for (i in 1:length(list.pca.data)) {
  myPlot <- ggplot(list.pca.data[[i]], aes(x=X, y=Y, colour = Sample)) +
    geom_point(aes(color = time.point)) +
    geom_path(data = list.elypses.data[[i]]) +
    geom_text_repel(aes(label = Sample, color = time.point), 
                    size = 3,
                    segment.color = 'black', segment.size = 0.2) +
    facet_wrap(species ~ data, scales = "free", nrow = 1) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size=0.5))) +
    labs(x = paste("PC3 - ", list.pca.data[[i]]$var.per[3], "%", sep=""), 
         y = paste("PC4 - ", list.pca.data[[i]]$var.per[4], "%", sep=""),
         color = "time-point") +
    theme(axis.title.x = element_text(size=11, face = "bold"),
          axis.title.y = element_text(size=11, face = "bold"),
          strip.text.x = element_text(size = 12))
  
  ggplotList[[i]] <- addSmallLegend(myPlot)
}


PCAbiologicalReplicates <- ggdraw() +
  draw_plot(ggplotList[[1]], 0, .75, .5, .25) +
  draw_plot(ggplotList[[2]], .5, .75, .5, .25) +
  
  draw_plot(ggplotList[[3]], 0, .5, .5, .25) +
  draw_plot(ggplotList[[4]], .5, .5, .5, .25) +
  
  draw_plot(ggplotList[[5]], 0, .25, .5, .25) +
  draw_plot(ggplotList[[6]], .5, .25, .5, .25) +
  
  draw_plot(ggplotList[[7]], 0, 0, .5, .25) +
  draw_plot(ggplotList[[8]], .5, 0, .5, .25) 
#draw_plot_label(c("A", "B", "C"), c(0, 0, 0.5), c(1, 0.5, 0.5), size = 15)

# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/supplementary_information/PC3_PC4.pdf", height = 12, width = 7)
print(PCAbiologicalReplicates)
dev.off()


##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##PCA_pc3_pc4_biological_replicates##") %>%
  body_add_gg(value = PCAbiologicalReplicates, width = 6.5, height = 9) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")






###################################################################
### Optional PCA for Mouse : normalized data VS raw counts data ###
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/liver/proteome/raw_protein_level.txt", head=TRUE, fill=TRUE, sep = "\t", check.names = FALSE)
replicates.data <- replicates.data.raw[, c(-1,-2,-3)]
replicates.data <- replicates.data[apply(replicates.data, 1, function(row) all(row !=0 )), ]
# remove zero variance rows from the dataset, setting variance not equal to zero (for the prcomp function)
for (i in 1:ncol(replicates.data)) {
  replicates.data <- replicates.data[apply(replicates.data, 1, var) != 0,]
}
timepoint.names <- paste("ZT", rep(c("00", "03", "06", "09", "12", "15", "18", "21"), 2), sep="")
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.data <- na.omit(replicates.data)
## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4])
# Make elypses
centroids <- aggregate(cbind(X,Y)~Sample, pca.data, mean)
row.names(centroids) <- centroids$Sample
elypses.data  <- do.call(rbind,lapply(unique(pca.data$Sample),function(t)
  data.frame(Sample=as.character(t),
             ellipse(cov(pca.data[pca.data$Sample==t,2:3]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.2),
             stringsAsFactors=FALSE)))

## redo PCA without elypses ##
replicates.data <- replicates.data.raw[, c(-1,-2,-3)]
# remove zero variance rows from the dataset, setting variance not equal to zero (for the prcomp function)
replicates.data <- replicates.data[apply(replicates.data, 1, function(row) all(row !=0 )), ]
for (i in 1:ncol(replicates.data)) {
  replicates.data <- replicates.data[apply(replicates.data, 1, var) != 0,]
}

#replicates.data <- na.omit(replicates.data)
## PCA ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,3],
                       Y=pca$x[,4],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:2, each=8), sep="")

pca.data <- cbind(pca.data, data.frame(species=rep("mouse", nrow(pca.data)),
                                       data=rep("proteins (raw counts)", nrow(pca.data))))
elypses.data <- cbind(elypses.data, data.frame(species=rep("mouse", nrow(elypses.data)),
                                               data=rep("proteins (raw counts)", nrow(elypses.data))))

ggplot(pca.data, aes(x=X, y=Y, colour = Sample)) +
  geom_point(aes(color = time.point)) +
  geom_path(data = elypses.data) +
  geom_text_repel(aes(label = Sample, color = time.point), size = 3) +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(x = paste("PC3 - ", pca.data$var.per[3], "%", sep=""), 
       y = paste("PC4 - ", pca.data$var.per[4], "%", sep=""),
       color = "time-point") +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        strip.text.x = element_text(size = 12))

