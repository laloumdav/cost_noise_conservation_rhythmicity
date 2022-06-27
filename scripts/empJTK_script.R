#################
### empJTK ###
#################

if (dir.exists("/scratch/cluster/monthly/dlaloum/")==TRUE){
  main.dir <- "/scratch/cluster/monthly/dlaloum/Documents/cost_theory_workingspace"
} else {
  main.dir <- "~/Documents/cost_theory_workingspace"
}

file.dir <- paste(paste(main.dir, "DATA/ostreococcus", sep = "/"), "/", sep="")
tissue.file <- "protein_level.txt"
file.name <- paste(file.dir, tissue.file, sep = "")
raw.dataset <- read.table(file.name, head=TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, fill=TRUE)

empJTK.python.directory <- paste(main.dir, "/scripts/empJTK_master/", sep = "")
# check if empJTK python script is present
if (file.exists(paste(empJTK.python.directory, "eJTK-CalcP.py", sep = ""))==FALSE){
  stop("empJTK files are missing. Should be in rhythm_detection_benchmark/scripts/empJTK_master/")
}



### SCRIPT ###
waveforms <- "ref_files/waveform_cosine.txt"
period <- "ref_files/period24.txt"
phases.searched <- "ref_files/phases_00-22_by2.txt"
asymmetries.searched <- "ref_files/asymmetries_02-22_by2.txt"
output <- "empJTK"

### Bash script for empJTK ###
system(paste("chmod 755 ", empJTK.python.directory, "eJTK-CalcP.py", sep=""), wait=TRUE)
system(paste("cd ", empJTK.python.directory, "bin/ ; python setup.py build_ext --inplace", sep=""), wait = TRUE)

empJTK_bash.script <- paste("./eJTK-CalcP.py",
                            "-f", file.name,
                            "-w", waveforms,
                            "-p", period,
                            "-s", phases.searched,
                            "-a", asymmetries.searched,
                            "-x", output,
                            sep = " ")

empJTK_bash.script <- paste("cd ", empJTK.python.directory, " ; ", empJTK_bash.script, sep="")
empJTK_bash.script <- paste("echo waveform searched:", " ; ", "head ", empJTK.python.directory, "ref_files/waveform_cosine.txt", " ; ", empJTK_bash.script, sep="")
empJTK_bash.script <- paste("echo period searched:", " ; ", "head ", empJTK.python.directory, "ref_files/period24.txt", " ; ", empJTK_bash.script, sep="")
system(empJTK_bash.script, wait=TRUE)

# Then we get back the "GAMMA" empJTK files:
empJTK.files <- list.files(file.dir)[grep(tissue, list.files(file.dir))]
empJTK.files <- empJTK.files[grep("empJTK", empJTK.files)]
empJTK.gamma.file <- empJTK.files[grep("Gamma", empJTK.files)]
if (length(empJTK.gamma.file)!=1) { print("ERROR: empJTK_Gamma_output_file doesn't exist or is not unique") }
empJTK.gamma.file <- paste(file.dir, empJTK.gamma.file, sep="")
empJTK <- read.table(empJTK.gamma.file, head=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# We keep only info needed
empJTK <- empJTK[,c("ID", "P", "empP", "Period", "Phase", "Max_Amp")]
colnames(empJTK) <- c("ID", "raw.pvalue", "default.pvalue", "period", "phase", "amplitude")
# We remove empJTK_output_files
empJTK.files <- paste(file.dir, empJTK.files, sep = "")
empJTK.files <- paste(empJTK.files, collapse = " ")
remove.empJTK.files_bash.script <- paste("rm", empJTK.files, sep = " ")
system(remove.empJTK.files_bash.script, wait = TRUE)


# Write the new one:
write.table(empJTK, paste(file.dir, "empJTK.txt", sep = ""), row.names = F, quote = F, sep = "\t")
