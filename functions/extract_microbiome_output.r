library(optparse)
library(tidyverse)
library(dplyr)
library(stringr)

option_list <- list(
  make_option(c("--sample_name"), action = "store", help = "sample name"),
  make_option(c("--output_file"), action = "store", help = "path to kraken output file"),
  make_option(c("--kraken_report"), action = "store", help = "path to kraken report"),
  make_option(c("--mpa_report"), action = "store", help = "path to standard kraken mpa report"),
  make_option(c("--out_path"), action = "store", help = "output path"),
  make_option(c("--keep_original"), action = "store", default = T, help = "delete original fastq file? T/F"),
  make_option(c("--ntaxid"), action = "store", default = 8000, help = "number of taxids to extract at a time")
)
opt <- parse_args(OptionParser(option_list = option_list))

kr <- read.delim(opt$kraken_report, header = FALSE)
kr <- kr[-c(1:2), ]
mpa <- read.delim(opt$mpa_report, header = FALSE)
n <- str_which(mpa$V1, "(?i)Bacteria|Fungi|Viruses")
taxid <- kr$V7[n]
taxid.list <- split(taxid, ceiling(seq_along(taxid) / opt$ntaxid))

file_name <- paste0(opt$sample_name, ".microbiome.output.txt")
if (file.exists(file.path(opt$out_path, file_name))) {
  file.remove(file.path(opt$out_path, file_name))
}

for (i in seq_along(taxid.list)) {
  print(paste("Extracting output data", i, "/", length(taxid.list)))

  taxid <- paste0("(taxid ", taxid.list[[i]], ")", collapse = "\\|")
  taxid <- paste0("'", taxid, "'")
  str <- paste0(
    "grep -w ", taxid, " ", opt$output_file, " >> ",
    file.path(opt$out_path, file_name)
  )
  system(str)
}

if (opt$keep_original == FALSE) {
  file.remove(opt$output_file)
}

print("Done")
