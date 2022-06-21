#!/usr/bin/env Rscript

aa_file = snakemake@input[[1]]
nt_file  = snakemake@input[[2]]
save_file = snakemake@output[[1]]

aa=read.table(file=aa_file,sep="\t",header=T,fill=T)
nt=read.table(file=nt_file,sep="\t",header=T,fill=T)

combine.df = function(df1,df2){
  print("Combine files")
  df=data.frame(rbind(df1,df2))
  print("Sort files by samples name")
  df=df[order(df$sample,decreasing=F),]
  return(df)
}

aa_nt = combine.df(aa,nt)

write.table(aa_nt, save_file, col.names=T,row.names=FALSE,sep="\t",quote=F)
