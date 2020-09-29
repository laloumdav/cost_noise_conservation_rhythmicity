##################################################
### Protein biosynthesis cost files preparation ###
##################################################
library(seqinr)
library(DescTools)
library(ggpubr)

main.dir <- "~/Documents/cost_theory_workingspace"

aa_biosynthesis_cost <- read.csv(paste(main.dir, "DATA/protein_AA_composition/AA_biosynthesis_cost.csv", sep="/"), skip=6)
colnames(aa_biosynthesis_cost)[1] <- "AA"

#############
#
#############
is.integer0 <- function(x)
{  is.integer(x) && length(x) == 0L }

fasta.files.list <- list.files(paste(main.dir, "DATA/protein_AA_composition/with_pi_s/", sep="/"))
fasta.files.list <- paste(main.dir, "DATA/protein_AA_composition/with_pi_s", grep(".fasta.gz", fasta.files.list, value = TRUE), sep="/")

for (n in 1:length(fasta.files.list)) {
  fasta.file <- read.fasta(file = fasta.files.list[n], seqtype = "AA")
  
  aa.synthesis.average.cost.per.protein <- data.frame(Protein.ID=NA, protein.length=NA, aa.synthesis.average.cost_Akashi=NA, aa.synthesis.average.cost_Wagner=NA)
  for (i in 1:length(fasta.file)) {
    protein.id <- names(fasta.file[i])
    protein.id <- gsub("\\..*", "", protein.id) # to deselect in case of Transcript.ID instead of Protein.ID
    if( !is.integer0(grep("[|]", protein.id)) ) {
      protein.id <- unlist(strsplit(protein.id, split = "[|]"))[2]
    }
    protein.length <- as.numeric(length(fasta.file[[i]]))
    aa.distribution <- data.frame(AA = Freq(fasta.file[[i]])[, "level"],
                                  freq = as.numeric(Freq(fasta.file[[i]])[, "freq"]))
    aa.distribution <- merge(aa.distribution, aa_biosynthesis_cost, by="AA")
    aa.distribution$average.cost.per.aa_Akashi <- aa.distribution$freq * aa.distribution$cost_Akashi / protein.length
    aa.distribution$average.cost.per.aa_Wagner <- aa.distribution$freq * aa.distribution$cost_Wagner / protein.length
    aa.synthesis.average.cost_Akashi <- round(sum(aa.distribution$average.cost.per.aa_Akashi), 2)
    aa.synthesis.average.cost_Wagner <- round(sum(aa.distribution$average.cost.per.aa_Wagner), 2)
    
    aa.synthesis.average.cost.per.protein[i, ] <- c(protein.id, protein.length, aa.synthesis.average.cost_Akashi, aa.synthesis.average.cost_Wagner)
  }
  # check: Must be TRUE
  isTRUE(length(unique(aa.synthesis.average.cost.per.protein$Protein.ID)) == nrow(aa.synthesis.average.cost.per.protein))
  
  final.file.name <- paste(gsub(".fasta.gz", "", fasta.files.list[n]), "_AA_synthesis_average_cost_per_protein.txt", sep="")
  write.table(aa.synthesis.average.cost.per.protein, final.file.name, quote = FALSE, sep = "\t", row.names = FALSE)
  
}
