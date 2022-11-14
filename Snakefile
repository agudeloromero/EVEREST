""" 
workflow Snakefile
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au

███████╗██╗   ██╗███████╗██████╗ ███████╗███████╗████████╗
██╔════╝██║   ██║██╔════╝██╔══██╗██╔════╝██╔════╝╚══██╔══╝
█████╗  ██║   ██║█████╗  ██████╔╝█████╗  ███████╗   ██║   
██╔══╝  ╚██╗ ██╔╝██╔══╝  ██╔══██╗██╔══╝  ╚════██║   ██║   
███████╗ ╚████╔╝ ███████╗██║  ██║███████╗███████║   ██║   
╚══════╝  ╚═══╝  ╚══════╝╚═╝  ╚═╝╚══════╝╚══════╝   ╚═╝   
"""

import os
DIR = os.getcwd()
configfile: "config/config.yaml"

SAMPLES,  = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))

if config["sequencing"] == 'PE':
		include: "SMK/01_raw_qc_PE.smk"
		include: "SMK/02_trimming_adaptors_PE.smk"
		include: "SMK/03_Host_removal_PE.smk"
		include: "SMK/04_DeNovo_PE.smk"
elif config["sequencing"] == 'SE':
		include: "SMK/01_raw_qc_SE.smk"
		include: "SMK/02_trimming_adaptors_SE.smk"
		include: "SMK/03_Host_removal_SE.smk"
		include: "SMK/04_DeNovo_SE.smk"
else:
		include: "SMK/01_raw_qc_Long.smk"
		include: "SMK/02_trimming_adaptors_Long.smk"
		include: "SMK/03_Host_removal_Long.smk"
		include: "SMK/04_DeNovo_Long.smk" 
include: "SMK/05_Cleaning_contigs.smk"
include: "SMK/06_Taxonomy_nt.smk"
include: "SMK/07_Taxonomy_aa.smk"
include: "SMK/08_Summary_Taxonomy.smk"

#if config["annotation"] == 'RASTt':
#		include: "SMK/07_annotation_RASTt.smk"
#else:
#		include: "SMK/07_annotation_Prokka.smk"

rule all:
	input:
		#rawQC
		os.path.join(config["output_DIR"],"EVEREST/multiQC_rep/raw_qc_multiqc_report.html"),
		#adaptorsQC
		os.path.join(config["output_DIR"],"EVEREST/multiQC_rep/trim_adaptors_multiqc_report.html"),
		#Host_removal
		os.path.join(config["output_DIR"],"EVEREST/multiQC_rep/fastq_before_merge_multiqc_report.html"),
		#DeNovo
#		expand(os.path.join(config["output_DIR"], "EVEREST/SPADES/{sample}/scaffolds.fasta"), sample=SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_all_seqs.fasta"), sample=SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_cluster.tsv"), sample=SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq.fasta"), sample=SAMPLES),
		#Viral contings enrichment
		expand(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"), sample=SAMPLES),
		#Classification
		os.path.join(config["output_DIR"], "EVEREST/Summary/Summary_nt_mmseqs2.txt"),
		os.path.join(config["output_DIR"], "EVEREST/Summary/Summary_aa_mmseqs2.txt"),
		#Summaries
		os.path.join(config["output_DIR"], "EVEREST/Summary/Summary_mmseqs2.txt"),

