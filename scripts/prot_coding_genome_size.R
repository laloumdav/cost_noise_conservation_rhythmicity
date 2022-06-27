##################################################
### Add the size of the genome coding for proteins ###
##################################################

main.dir <- "~/Documents/cost_theory_workingspace"

files.list <- list.files(paste(main.dir, "DATA/protein_AA_composition/with_pi_s/", sep="/"))
files.list <- paste(main.dir, "DATA/protein_AA_composition/with_pi_s", grep("AA_synthesis_average_cost_per_protein.txt", files.list, value = TRUE), sep="/")

for (n in 1:length(files.list)) {
  AA.file <- read.table(files.list[n], sep="\t", head=TRUE, fill=TRUE, check.names = FALSE)
  AA.file$protein.coding.genome.size = sum(AA.file$protein.length)

  write.table(AA.file, files.list[n], quote = FALSE, sep = "\t", row.names = FALSE)
  
}
