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

SAMPLES, = glob_wildcards(os.path.join(config["samples_DIR"],"{sample}_R1.fastq.gz"))

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
include: "SMK/06_Taxonomy_aa.smk"
include: "SMK/06_Taxonomy_nt.smk"
#if config["annotation"] == 'RASTt':
#		include: "SMK/07_annotation_RASTt.smk"
#else:
#		include: "SMK/07_annotation_Prokka.smk"

rule all:
	input:
		#rawQC
		"results/multiQC_rep/raw_qc_multiqc_report.html",
		#adaptorsQC
		"results/multiQC_rep/trim_adaptors_multiqc_report.html",
		#Host_removal
		"index_minimap2/index",
		"results/multiQC_rep/fastq_before_merge_multiqc_report.html",
		#DeNovo
		expand("results/SPADES/{sample}/scaffolds.fasta", sample=SAMPLES),
		#Viral contings enrichment
		expand("results/CHECKV/{sample}/viruses_rename.fna", sample=SAMPLES),
		#Classification
		expand("results/Summary/{sample}_nt_summary_mmseqs2.txt", sample=SAMPLES),
		expand("results/Summary/{sample}_aa_summary_mmseqs2.txt", sample=SAMPLES),
		"results/Summary/Summary_mmseqs2.txt",
		"results/Summary/Summary_aa_mmseqs2.txt",
		"results/Summary/Summary_nt_mmseqs2.txt",

