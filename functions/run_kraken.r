library(optparse)

option_list = list(
  make_option(c("--sample"), action="store", help = "sample name"),
  make_option(c("--fq1"), action="store", help = "path to fastq 1 file"),
  make_option(c("--fq2"), action="store", help = "path to fastq 2 file"),
  make_option(c("--out_path"), action="store", help = "output directory path"),
  make_option(c("--ncbi_blast_path"), action="store", help = "path to ncbi-blast"),
  make_option(c("--Kraken2Uniq_path"), action="store", help = "path to Kraken2 main 'kraken2' function"),
  make_option(c("--kraken_database_path"), action="store", help = "path to kraken database")
)
opt = parse_args(OptionParser(option_list = option_list))

library(stringr)

# run Kraken
str = paste0(str,
             'export PATH=$PATH:', opt$ncbi_blast_path, ' \n\n',
             opt$Kraken2Uniq_path, ' \\\n',
             '--db ', opt$kraken_database_path, ' \\\n',
             '--threads 24 \\\n',
             '--paired \\\n',
             '--use-names \\\n',
             '--report-minimizer-data \\\n',
             '--classified-out ', paste0(opt$out_path, opt$sample), '#.fq', ' \\\n',
             '--output ', paste0(opt$out_path, opt$sample, '.kraken.output.txt'), ' \\\n',
             '--report ', paste0(opt$out_path, opt$sample, '.kraken.report.txt'), ' \\\n',
             paste0(f1[i], ' \\\n'),
             paste0(f2[i], '\n\n')
)


# create MPA style report and standard Kraken report
str = paste0(str, 
             'cut -f1-3,6-8 ', 
             paste0(opt$out_path, opt$sample, '.kraken.report.txt'),
             ' > ',
             paste0(opt$out_path, opt$sample, '.kraken.report.std.txt'),
             '\n\n',
             '/projects/sd948/bassel/Kraken2Uniq/kraken2-master/kreport2mpa.py \\\n',
             '-r ', paste0(opt$out_path, opt$sample, '.kraken.report.std.txt'), ' \\\n',
             '-o ', paste0(opt$out_path, opt$sample, '.kraken.report.mpa.txt'), ' \\\n',
             '--intermediate-ranks',
             '\n\n'
)

system(str)