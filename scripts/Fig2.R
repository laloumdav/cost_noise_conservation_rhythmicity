####################################################
####################################################
#########   Main figures grid arrangement   ########
####################################################
####################################################
library(cowplot)
library(ggpubr)
library(magick)

rds.objects.dir <- "~/Documents/cost_theory_workingspace/rds_objects/"
#list.files(rds.objects.dir)

AAsynthCostComparedPerGroup <- readRDS(paste(rds.objects.dir, "AAsynthCostComparedPerGroup.rds", sep=""))
totalCostComparedPerGroup <- readRDS(paste(rds.objects.dir, "totalCostComparedPerGroup.rds", sep=""))
proteins_maxMeanExprComparedPerGroup <- readRDS(paste(rds.objects.dir, "proteins_maxMeanExprComparedPerGroup.rds", sep=""))
proteins_meanMeanExprComparedPerGroup <- readRDS(paste(rds.objects.dir, "proteins_meanMeanExprComparedPerGroup.rds", sep=""))
proteinLengthComparedPerGroup <- readRDS(paste(rds.objects.dir, "proteinLengthComparedPerGroup.rds", sep=""))
averageAASynthCostsHistograms <- readRDS(paste(rds.objects.dir, "averageAASynthCostsHistograms.rds", sep=""))
#aMolWeightComparedPerGroup <- readRDS(paste(rds.objects.dir, "MolWeightComparedPerGroup.rds", sep=""))

proteins_meanMeanExprComparedPerGroup <- proteins_meanMeanExprComparedPerGroup + rremove("legend")
proteins_maxMeanExprComparedPerGroup <- proteins_maxMeanExprComparedPerGroup + rremove("legend")
AAsynthCostComparedPerGroup <- AAsynthCostComparedPerGroup + rremove("legend")
proteinLengthComparedPerGroup <- proteinLengthComparedPerGroup + rremove("legend")
#MolWeightComparedPerGroup <- MolWeightComparedPerGroup + rremove("legend")


adjusted.first.row_1 <- plot_grid(NULL,
                                  totalCostComparedPerGroup,
                                  #NA,
                                  labels = c("", "a"),
                                  ncol = 1, nrow = 2,
                                  rel_heights = c(.4, .6),
                                  label_size = 14)
adjusted.first.row_2 <- plot_grid(AAsynthCostComparedPerGroup,
                                  proteinLengthComparedPerGroup,
                                  labels = c("b", "c"),
                                  ncol = 1, nrow = 2,
                                  label_size = 14)


logo <- paste(rds.objects.dir, "formula.png", sep="")
row.1 <- plot_grid(adjusted.first.row_1,
          adjusted.first.row_2,
          #last.column.grid,
          labels = c("", ""),
          ncol = 2, nrow = 1,
          label_size = 14)

row.2 <- plot_grid(proteins_meanMeanExprComparedPerGroup,
                   proteins_maxMeanExprComparedPerGroup,
                   labels = c("d", "e"),
                   ncol = 2, nrow = 1,
                   label_size = 14)

plots <- plot_grid(row.1,
                   row.2,
                   ncol = 1, nrow = 2,
                   rel_heights = c(.7, .3))

final.draw.plot <- ggdraw() +
  draw_plot(plots) +
  draw_image(logo, x = .58, y = .72,
             width = 0.67, height = 0.27,
             hjust = 1)

pdf("~/Documents/cost_theory_workingspace/Fig2.pdf", width = 9, height = 9) # or other device
#pdf("~/Documents/cost_theory_workingspace/Fig2_Akashi.pdf", width = 9, height = 9) # or other device
print(final.draw.plot)
dev.off()


