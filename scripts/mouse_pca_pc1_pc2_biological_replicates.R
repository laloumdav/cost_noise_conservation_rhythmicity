library(ggplot2)
library(ellipse)
library(ggrepel) 
library(cowplot)
library(officer)
library(magrittr)
### ### ### ### ### ### ### ### ### 
### PCA on biological replicates ##   PC1 & PC2
### ### ### ### ### ### ### ### ### 
list.pca.data <- list()
list.elypses.data <- list()


#1.1# MOUSE Tendon transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/tendon/proteome/protein_level.txt", h=T, fill=T, sep="\t")
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("CT", rep(c("03", "07", "11", "15", "19", "23"), 2), sep="")
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.data <- na.omit(replicates.data)

## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
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
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:2, each=6), sep="")

# Store in lists
list.pca.data[[1]] <- cbind(pca.data, data.frame(species=rep("Mouse\nTendon", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[1]] <- cbind(elypses.data, data.frame(species=rep("Mouse\nTendon", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))




#1.2# MOUSE Cartilage transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/cartilage/proteome/protein_level.txt", h=T, fill=T, sep="\t")
replicates.data <- replicates.data.raw[, c(-1,-2)]
timepoint.names <- paste("ZT", rep(c("02", "06", "10", "14", "18", "22"), 2), sep="")
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.data <- na.omit(replicates.data)

## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
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
replicates.data <- replicates.data.raw[, c(-1,-2)]
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

#replicates.timepoint <- na.omit(replicates.timepoint)
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:2, each=6), sep="")

# Store in lists
list.pca.data[[2]] <- cbind(pca.data, data.frame(species=rep("Mouse\nCartilage", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[2]] <- cbind(elypses.data, data.frame(species=rep("Mouse\nCartilage", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))





#1.3# MOUSE Forebrain transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/forebrain/proteome/protein_level.txt", h=T, fill=T, sep="\t")
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("ZT", c("20", "00", rep(c("04", "08", "12", "16", "20", "00"), 2)), sep="")
colnames(replicates.data) <- timepoint.names
rownames(replicates.data) <- paste("gene", 1:nrow(replicates.data), sep="")

replicates.data <- na.omit(replicates.data)

## PCA for elypses first ##
pca <- prcomp(t(replicates.data), scale=TRUE) 
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
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
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", c(rep(1:2, each=6), 3, 3), sep="")

# Store in lists
list.pca.data[[3]] <- cbind(pca.data, data.frame(species=rep("Mouse\nForebrain", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[3]] <- cbind(elypses.data, data.frame(species=rep("Mouse\nForebrain", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))





### ### ### ### ### ### ### ### ### 
### PCA on biological replicates ##   PC3 & PC4
### ### ### ### ### ### ### ### ### 


#1.1# MOUSE Tendon transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/tendon/proteome/protein_level.txt", h=T, fill=T, sep="\t")
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("CT", rep(c("03", "07", "11", "15", "19", "23"), 2), sep="")
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
pca.data$Sample <- paste("rep", rep(1:2, each=6), sep="")

# Store in lists
list.pca.data[[4]] <- cbind(pca.data, data.frame(species=rep("Mouse\nTendon", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[4]] <- cbind(elypses.data, data.frame(species=rep("Mouse\nTendon", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))




#1.2# MOUSE Cartilage transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/cartilage/proteome/protein_level.txt", h=T, fill=T, sep="\t")
replicates.data <- replicates.data.raw[, c(-1,-2)]
timepoint.names <- paste("ZT", rep(c("02", "06", "10", "14", "18", "22"), 2), sep="")
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
replicates.data <- replicates.data.raw[, c(-1,-2)]
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
pca.data$Sample <- paste("rep", rep(1:2, each=6), sep="")

# Store in lists
list.pca.data[[5]] <- cbind(pca.data, data.frame(species=rep("Mouse\nCartilage", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[5]] <- cbind(elypses.data, data.frame(species=rep("Mouse\nCartilage", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))





#1.3# MOUSE Forebrain transcriptome
######################################
replicates.data.raw <- read.table("~/Documents/cost_theory_workingspace/DATA/mouse/forebrain/proteome/protein_level.txt", h=T, fill=T, sep="\t")
replicates.data <- replicates.data.raw[, -1]
timepoint.names <- paste("ZT", c("20", "00", rep(c("04", "08", "12", "16", "20", "00"), 2)), sep="")
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
pca.data$Sample <- paste("rep", c(rep(1:2, each=6), 3, 3), sep="")

# Store in lists
list.pca.data[[6]] <- cbind(pca.data, data.frame(species=rep("Mouse\nForebrain", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[6]] <- cbind(elypses.data, data.frame(species=rep("Mouse\nForebrain", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))




### ### ### ### ### ### ### ### ### 
### PCA on biological replicates ##   PC5 & PC6
### ### ### ### ### ### ### ### ### 

#2.1# Mouse liver proteome
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
                       X=pca$x[,5],
                       Y=pca$x[,6])
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
                       X=pca$x[,5],
                       Y=pca$x[,6],
                       var.per=pca.var.per)
pca.data$time.point <- timepoint.names
# ordering time-points 
pca.data$time.point <- factor(pca.data$time.point, level = unique(timepoint.names))
# change names of points
pca.data$Sample <- paste("rep", rep(1:2, each=8), sep="")

# Store in lists
list.pca.data[[7]] <- cbind(pca.data, data.frame(species=rep("mouse\nLiver", nrow(pca.data)),
                                                 data=rep("proteins", nrow(pca.data))))
list.elypses.data[[7]] <- cbind(elypses.data, data.frame(species=rep("mouse\nLiver", nrow(elypses.data)),
                                                         data=rep("proteins", nrow(elypses.data))))



### ### ### ### ### ### ### ### ### 
### PCA on biological replicates ##   PC1 & PC4
### ### ### ### ### ### ### ### ### 

#2.1# Mouse liver proteome
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
                       X=pca$x[,1],
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
                       X=pca$x[,1],
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
for (i in 1:3) {
  myPlot <- ggplot(list.pca.data[[i]], aes(x=X, y=Y, colour = Sample)) +
    geom_point(aes(color = time.point)) +
    geom_path(data = list.elypses.data[[i]]) +
    geom_text_repel(aes(label = Sample, color = time.point), 
                    size = 3,
                    segment.color = 'black', segment.size = 0.2) +
    facet_wrap(species ~ data, scales = "free", nrow = 1) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size=0.5))) +
    labs(x = paste("PC1 - ", list.pca.data[[i]]$var.per[1], "%", sep=""), 
         y = paste("PC2 - ", list.pca.data[[i]]$var.per[2], "%", sep=""),
         color = "time-point") +
    theme(axis.title.x = element_text(size=11, face = "bold"),
          axis.title.y = element_text(size=11, face = "bold"),
          strip.text.x = element_text(size = 12))
  
  ggplotList[[i]] <- addSmallLegend(myPlot)
}

for (i in 4:6) {
  myPlot <- ggplot(list.pca.data[[i]], aes(x=X, y=Y, colour = Sample)) +
    geom_point(aes(color = time.point)) +
    geom_path(data = list.elypses.data[[i]]) +
    geom_text_repel(aes(label = Sample, color = time.point), 
                    size = 3,
                    segment.color = 'black', segment.size = 0.2) +
    facet_wrap(species ~ data, scales = "free", nrow = 1) +
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size=0.5))) +
    labs(x = paste("PC3 - ", list.pca.data[[i]]$var.per[1], "%", sep=""), 
         y = paste("PC4 - ", list.pca.data[[i]]$var.per[2], "%", sep=""),
         color = "time-point") +
    theme(axis.title.x = element_text(size=11, face = "bold"),
          axis.title.y = element_text(size=11, face = "bold"),
          strip.text.x = element_text(size = 12))
  
  ggplotList[[i]] <- addSmallLegend(myPlot)
}

i <- 7
myPlot <- ggplot(list.pca.data[[i]], aes(x=X, y=Y, colour = Sample)) +
  geom_point(aes(color = time.point)) +
  geom_path(data = list.elypses.data[[i]]) +
  geom_text_repel(aes(label = Sample, color = time.point), 
                  size = 3,
                  segment.color = 'black', segment.size = 0.2) +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=0.5))) +
  labs(x = paste("PC5 - ", list.pca.data[[i]]$var.per[1], "%", sep=""), 
       y = paste("PC6 - ", list.pca.data[[i]]$var.per[2], "%", sep=""),
       color = "time-point") +
  theme(axis.title.x = element_text(size=11, face = "bold"),
        axis.title.y = element_text(size=11, face = "bold"),
        strip.text.x = element_text(size = 12))
ggplotList[[i]] <- addSmallLegend(myPlot)

i <- 8
myPlot <- ggplot(list.pca.data[[i]], aes(x=X, y=Y, colour = Sample)) +
  geom_point(aes(color = time.point)) +
  geom_path(data = list.elypses.data[[i]]) +
  geom_text_repel(aes(label = Sample, color = time.point), 
                  size = 3,
                  segment.color = 'black', segment.size = 0.2) +
  facet_wrap(species ~ data, scales = "free", nrow = 1) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=0.5))) +
  labs(x = paste("PC1 - ", list.pca.data[[i]]$var.per[1], "%", sep=""), 
       y = paste("PC4 - ", list.pca.data[[i]]$var.per[2], "%", sep=""),
       color = "time-point") +
  theme(axis.title.x = element_text(size=11, face = "bold"),
        axis.title.y = element_text(size=11, face = "bold"),
        strip.text.x = element_text(size = 12))
ggplotList[[i]] <- addSmallLegend(myPlot)


PCAbiologicalReplicates <- ggdraw() +
  draw_plot(ggplotList[[1]], 0, .75, .5, .25) +
  draw_plot(ggplotList[[4]], .5, .75, .5, .25) +
  
  draw_plot(ggplotList[[2]], 0, .5, .5, .25) +
  draw_plot(ggplotList[[5]], .5, .5, .5, .25) +
  
  draw_plot(ggplotList[[3]], 0, .25, .5, .25) +
  draw_plot(ggplotList[[6]], .5, .25, .5, .25) +
  
  draw_plot(ggplotList[[7]], 0, 0, .5, .25) +
  draw_plot(ggplotList[[8]], .5, 0, .5, .25) 


# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/supplementary_information/PCA_mouse_proteom.pdf", height = 12, width = 7)
print(PCAbiologicalReplicates)
dev.off()

