#!/usr/bin/env Rscript
##!/usr/bin/Rscript
##! /path/to/Rscript

#args = commandArgs(trailingOnly=TRUE)
#lca_file  = args[1]
#aln_file  = args[2]
#name_DB   = args[3]
#save_file = args[4]

lca_file  = snakemake@input[[1]]
aln_file  = snakemake@input[[2]]
name_DB   = snakemake@input[[3]]
save_file = snakemake@output[[1]]

#-- lca
print("reading file - lca")
tax=read.table(file=lca_file,sep="\t",header=T,fill=T)
tax[tax==""] <- "NA"
tax2 = tax$lca_family

tax2_validation = all(is.na(tax2))

if (tax2_validation == TRUE) {
		  print("empty file")
		  df = data.frame(matrix(nrow=0, ncol=0))
write.table(df,file=save_file,col.names=T,row.names=F,sep="\t",quote=F)
} else {
  print("non empty file")
#-- Baltimore
print("Baltimore - bal")
Balti=read.table(file=name_DB,sep="\t",header=T)
family_balti <- function(family) {
		df1=data.frame(matrix(nrow=0, ncol=3))
		colnames(df1)=c("Baltimore_group","Baltimore_Genome_composition","ICTV_Year")

		line = Balti[grep(family, Balti$Family),]
		line = line[c(2:4)]

		if (nrow(line) == 0){
				df1[nrow(df1) + 1,] <- c(NA, NA, NA)
		} else{
				df1[nrow(df1) + 1,] <- line
		}
		return(df1)
}

df = data.frame(matrix(nrow=0, ncol=3))
colnames(df)=c("Baltimore_group","Baltimore_Genome_composition","ICTV_Year")

for (i in tax2) {
  name <- c("Balti_group")
  x = family_balti(i)
  df = data.frame(rbind(df,x))
  assign(name,df)
}

tax_v2 = data.frame(cbind(tax,Balti_group))
tax_v3 = tax_v2[!grepl("^0", tax_v2$lca_taxid),]

#-- Identification
print("Sample identification")
basename = basename(lca_file)
base = gsub("_tax_viral_(.*)_lca_reformated_header.tsv","",as.character(basename))
sample = rep(base, nrow(tax_v3))
mmseq.db = grepl("*_nt_*", basename)

if (mmseq.db == TRUE) {
		  database = rep("nucleotide", nrow(tax_v3))
} else {
		  database = rep("amino acid", nrow(tax_v3))
}
tax_v4 = data.frame(cbind(sample,database,tax_v3))

#-- Aligment
print("Aligment - aln")
aln=read.table(file=aln_file,sep="\t",header=T,fill=T)

alingment <- function(contig) {
		x1 = aln[grepl(contig, aln$aln_query),]
		x2 = x1[order(x1$aln_bits, decreasing = T),]
		x3 = data.frame(x2[1,])
		x4 = subset(x3, select = c("aln_query","aln_target","aln_taxid","aln_taxname","aln_evalue","aln_pident","aln_mismatch","aln_alnlen","aln_bits") )
		return(x4)
}

aln.query = tax_v4$lca_query
aln.df = data.frame(matrix(nrow=0, ncol=9))
colnames(aln.df)=c("aln_query","aln_target","aln_taxid","aln_taxname","aln_evalue","aln_pident","aln_mismatch","aln_alnlen","aln_bits")
for (x in aln.query) {
		name <- c("Aligment")
		x = alingment(x)
		aln.df = data.frame(rbind(aln.df,x))
		assign(name,aln.df)
}

#- Final file
tax_v5 = data.frame(cbind(tax_v4,Aligment))

write.table(tax_v5,file=save_file,col.names=T,row.names=F,sep="\t",quote=F)
}
