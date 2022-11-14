"""
workflow 05_Cleaning_contigs.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd() 

#configfile: "config/config.yaml"
#SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))

#rule all:
#	input:
#		expand(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"), sample = SAMPLES),

rule SEQKIT_filter:
	input:
		fasta  = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq.fasta"),
	output:
		filter = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq_FilterLen.fasta"),
	params:
		"seq -m 5000",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S02_SEQKIT_filter_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S02_SEQKIT_filter_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/seqkit.yml"),
	message:
		"remove contigs < 5000 bp",
	shell:
		(" seqkit {params} {input.fasta} -o {output.filter} 2> {log} ")

checkpoint VIRSORTER_detect:
	input:
		fasta = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq_FilterLen.fasta"),
		db = config["VIRSORTER_DB"],
	output:
		dir = directory(os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}")),
		viral = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-combined.fa"),
		score = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-score.tsv"),
		boundary = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-boundary.tsv"),
	params:
		type = "run --keep-original-seq",
		groups = "--include-groups dsDNAphage,NCLDV,RNA,ssDNA",
		others = "--min-length 5000 --min-score 0.5",
		threads = "-j 7 all",
#		dir  = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/"),
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S03_VIRSORTER_detect_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S03_VIRSORTER_detect_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/virsorter2.yml"),
	message:
		"rid off from non-viral contigs",
	shell:
		(" virsorter {params.type} -i {input.fasta} -w {output.dir} --db-dir {input.db} {params.groups} {params.others} {params.threads} 2> {log} ")

def get_VIRSORTER_files(wildcards):
	checkpoint_output = checkpoints.VIRSORTER_detect.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output,"final-viral-combined.fa"))

checkpoint CHECKV_viral_seq:
	input:
		fasta = get_VIRSORTER_files,
		DB = config["CHECKV_DB"],
#		fasta = os.path.join(config["output_DIR"], "EVEREST/VIRSORTER/{sample}/final-viral-combined.fa"),
	output:
		dir  = directory(os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}")),
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses.fna"),
	params:
		type = "end_to_end",
		threads = "-t 7",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S04_CHECKV_viral_seq_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S04_CHECKV_viral_seq_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/checkv.yml"),
	message:
		"rid off non-viral reads",
	shell:
		(" checkv {params.type} -d {input.DB} {input.fasta} {output.dir} 2> {log} ")

def get_CHECKV_files(wildcards):
	checkpoint_output = checkpoints.CHECKV_viral_seq.get(**wildcards).output[0]
	return expand(os.path.join(checkpoint_output,"viruses.fna"))

rule RENAME_viral_seq:
	input:
		get_CHECKV_files,
#		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses.fna"),
	output:
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"),
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S05_RENAME_viral_seq_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S05_RENAME_viral_seq_{sample}.txt"),
	message:
		"rename viral contigs",
	shell:
		(" sed 's/||.*//' {input} > {output.fasta} 2> {log} ")
