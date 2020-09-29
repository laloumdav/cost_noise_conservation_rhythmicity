library(ggplot2)
library(ggridges)
library(officer)
library(magrittr)

species <- c("arabidopsis", "cyanobacteria", "ostreococcus", "mouse")
tissue <- c("leaves", NA, NA, "liver")

main.dir <- gsub("/NA", "", paste("~/Documents/cost_theory_workingspace/DATA", species, tissue, sep="/"))
proteome.dir <- paste(main.dir, "proteome", sep="/")
transcriptome.dir <- paste(main.dir, "transcriptome", sep="/")

transcriptome.data <- list()
proteome.data <- list()
tot.data <- data.frame(average.AA.synth.cost=NA, 
                       mean.value=NA,
                       species=NA)
for (i in 1:length(species)) {
  proteome.data[[i]] <- read.table(paste(proteome.dir[i], "proteome_data.txt", sep="/"), head=TRUE, fill=TRUE, check.names = FALSE)
  
  tot.data.proteins <- data.frame(average.AA.synth.cost=proteome.data[[i]]$aa.synthesis.average.cost,
                                  mean.value=round(mean(proteome.data[[i]]$aa.synthesis.average.cost, na.rm = TRUE), 2),
                                  species=rep(species[i], nrow(proteome.data[[i]])))
  
  tot.data <- rbind(tot.data, tot.data.proteins)
}

tot.data <- tot.data[-1,]

tot.data <- transform(tot.data, species = factor(species, levels=c("mouse", "arabidopsis", "cyanobacteria", "ostreococcus"),
                                    labels=c(expression("mouse (Ne ~ 2,5-12 ."*10^"6"*")"),
                                             expression("arabidopsis (Ne ~ 127 ."*10^"6"*")"),
                                             expression("cyanobacteria (Ne ~ 142 ."*10^"9"*")"),
                                             expression("ostreococcus (Ne ~ 6-29 ."*10^"6"*")"))))

tot.data.mean <- unique(tot.data[, c("species", "mean.value")])
tot.data.mean$y.axis.text <- c(56, 16, 24, 158)
tot.data.mean$x.axis.text <- c(20.9, 19.9, 20.2, 21.2)

### PLOTS ###
averageAASynthCostsHistograms <- ggplot(tot.data, aes(x=average.AA.synth.cost)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.8, na.rm = TRUE, bins = 100) + 
  geom_vline(aes(xintercept = mean.value), color="blue", linetype="dashed") +
  xlim(16, 24) +
  geom_text(data = tot.data.mean, aes(label=mean.value, x = x.axis.text, y = y.axis.text),
            color="blue", size = 3, 
            check_overlap = TRUE) +
  facet_wrap( ~ species, scales = "free", ncol = 2, labeller = label_parsed) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "averaged AA synthesis costs") 

########
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##average_AA_synth_costs_histograms##") %>%
  body_add_gg(value = averageAASynthCostsHistograms, width = 5) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")


########
# Also save it as an R object:
averageAASynthCostsHistograms <- ggplot(tot.data, aes(x=average.AA.synth.cost)) + 
  geom_histogram(aes(y = ..count..), position = 'identity', alpha=0.8, na.rm = TRUE, bins = 100) + 
  xlim(16, 24) +
  geom_vline(aes(xintercept = mean.value), color="blue", linetype="dashed") +
  geom_text(data = tot.data.mean, aes(label=mean.value, x = x.axis.text, y = y.axis.text),
            color="blue", size = 3, 
            check_overlap = TRUE) +
  facet_wrap( ~ species, scales = "free", ncol = 2, labeller = label_parsed) +
  theme_bw(base_size = 11, base_rect_size = 0.1) +
  scale_color_grey() +
  labs(x = "averaged AA synthesis costs") +
  theme(axis.text.y = element_text(size=5), 
      axis.title.y = element_text(size=7), 
      axis.text.x = element_text(size=5), 
      axis.title.x = element_text(size=7), 
      strip.background = element_blank(),
      strip.text = element_text(size=7, face = "bold", margin = margin(0)),
      axis.ticks.x = element_line(size = 0.1),
      axis.line.x = element_line(size = 0.1),
      axis.ticks.y = element_line(size = 0.1),
      axis.line.y = element_line(size = 0.1))

out.rds.name <- paste("~/Documents/cost_theory_workingspace/rds_objects", "averageAASynthCostsHistograms.rds", sep="/")
saveRDS(averageAASynthCostsHistograms, file = out.rds.name)


