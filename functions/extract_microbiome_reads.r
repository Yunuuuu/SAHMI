library(optparse)
library(tidyverse)
library(dplyr)
library(stringr)

option_list <- list(
  make_option(c("--sample_name"), action = "store", help = "sample name"),
  make_option(c("--fq"), action = "store", help = "path to fastq file"),
  make_option(c("--kraken_report"), action = "store", help = "path to standard kraken report"),
  make_option(c("--mpa_report"), action = "store", help = "path to standard kraken mpa report"),
  make_option(c("--out_path"), action = "store", help = "output path"),
  make_option(c("--keep_original"), action = "store", default = TRUE, help = "delete original fastq file? T/F"),
  make_option(c("--ntaxid"), action = "store", default = 7000, help = "number of taxids to extract at a time")
)
opt <- parse_args(OptionParser(option_list = option_list))

kr <- read.delim(opt$kraken_report, header = FALSE)
kr <- kr[-c(1:2), ]
mpa <- read.delim(opt$mpa_report, header = FALSE)
n <- str_which(mpa$V1, "(?i)Bacteria|Fungi|Viruses")
taxid <- kr$V7[n]
taxid.list <- split(taxid, ceiling(seq_along(taxid) / opt$ntaxid))
file_name <- paste0(opt$sample_name, "_line_numbers.txt")
if (file.exists(file.path(opt$out_path, file_name))) {
  file.remove(file.path(opt$out_path, file_name))
}

for (i in seq_along(taxid.list)) {
  print(paste("Finding reads", i, "/", length(taxid.list)))

  taxid <- paste0("taxid|", taxid.list[[i]], collapse = "\\|")
  taxid <- paste0("'", taxid, "'")
  str <- paste0("grep -wn ", taxid, " ", opt$fq, " | grep -Eo '^[^:]+' >> ", opt$out_path, opt$sample_name, "_line_numbers.txt")
  system(str)
}

h <- read.delim(file.path(opt$out_path, file_name), header = FALSE)
r <- h$V1 + 1
d <- data.frame(h = h$V1, r = r) %>%
  rownames_to_column("n") %>%
  pivot_longer(-n)
write.table(d$value,
  file = file.path(opt$out_path, file_name),
  row.names = FALSE, col.names = FALSE
)

print("Extracting reads")
str <- paste0(
  "awk 'NR==FNR{ a[$1]; next }FNR in a' ",
  file.path(opt$out_path, file_name), " ",
  opt$fq, " > ",
  file.path(opt$out_path, paste0(opt$sample_name, ".fa"))
)
system(str)

str <- paste0(
  "sed -i 's/@/>/g' ",
  file.path(opt$out_path, paste0(opt$sample_name, ".fa"))
)
system(str)

file.remove(file.path(opt$out_path, file_name))

if (!as.logical(opt$keep_original)) {
  file.remove(opt$fq)
}

print("Done")
