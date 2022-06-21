#!/usr/bin/env Rscript

dir = snakemake@params[[1]]
pattern.o  = snakemake@params[[2]]
save_file = snakemake@output[[1]]

#dir = "/Users/pagudeloromero/Downloads/Summary"
#pattern.o="*_nt_summary_mmseqs2.txt"
#save_file = "/Users/pagudeloromero/Downloads/Summary/Summary_aa_mmseqs2.txt"

combine = function(dir,pattern){
  #-- Create list of text files
  print("Create list of text files")
  file_list = list.files(path=dir, pattern=pattern,full.names=T)
  #-- Read the files
  print("Reading files")
  files <- lapply(file_list, function(x) {read.table(file = x, header = T, sep ="\t")})
  # Combine them
  print("Combining files")
  df <- do.call("rbind", lapply(files, as.data.frame))
  return(df)
}

mmseq.db = grepl("*_nt_*", pattern.o)
if (mmseq.db == TRUE) {
  print("-- nt --")
  pattern.n  = "*_nt_summary_mmseqs2.txt"
  combined_df = combine(dir,pattern.n)
} else {
  #-- Create list of text files
  print("-- aa --")
  pattern.a  = "*_aa_summary_mmseqs2.txt"
  combined_df = combine(dir,pattern.a)
}

print("Save file")
write.table(combined_df,file=save_file,col.names=T,row.names=F,sep="\t",quote=F)
