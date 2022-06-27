library(ggplot2)
library(ggrepel)
library(latex2exp)
library(officer)
library(magrittr)

main.dir <- "~/Documents/cost_theory_workingspace"
aa_biosynthesis_cost <- read.csv(paste(main.dir, "DATA/AA_biosynthesis_cost.csv", sep="/"), skip=6)

AAcostsLinearRelationships <- ggplot(aa_biosynthesis_cost, aes(x= cost_Akashi, y= cost_Wagner, colour="indianred2")) +
  geom_abline(intercept = 0, slope = 1, size = .2) +
  geom_point(color="black") +
  xlim(0, 80) + ylim(0, 80) +
  geom_label_repel(aes(label = Amino.Acid),
                   box.padding   = 0.5, 
                   point.padding = 0.3,
                   segment.color = 'indianred2') +
  #theme_classic() +
  labs(y = TeX("$\\textbf{c_{AA}}$ (Wagner, 2005)"), 
       x = TeX("$\\textbf{c_{AA}}$ (Akashi and Gojobori, 2002)")) +
  theme(axis.text.y = element_text(size=10), 
        axis.text.x = element_text(size=10), 
        legend.position = "none")

# Export as pdf image
pdf("~/Documents/cost_theory_workingspace/additional_files/AAcostsLinearRelationships.pdf", height = 5, width = 5)
print(AAcostsLinearRelationships)
dev.off()

##############
# Add the plots to the word document:
read_docx(path = "~/Documents/cost_theory_workingspace/figures_rmd.docx") %>%
  cursor_begin() %>% 
  
  cursor_reach(keyword = "##AA_costs_Linear_Relationships##") %>%
  body_add_gg(value = AAcostsLinearRelationships, width = 4, height = 4) %>%
  cursor_forward() %>%
  body_remove()  %>%
  
  cursor_end() %>%
  print(target = "~/Documents/cost_theory_workingspace/figures_rmd.docx")


