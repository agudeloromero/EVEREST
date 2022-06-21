"""
workflow 02_trimming_adaptors.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

rule TRIMM_PE:
	input:
		f1 = os.path.join(config["samples_DIR"], "{sample}_R1.fastq.gz"),
		f2 = os.path.join(config["samples_DIR"], "{sample}_R2.fastq.gz"),
		fasta = config["adaptor"],
	output:
		p1  = temp("results/TRIMM/{sample}_trimm_pair_R1.fastq.gz"),
		up1 = temp("results/TRIMM/{sample}_trimm_unpair_R1.fastq.gz"),
		p2  = temp("results/TRIMM/{sample}_trimm_pair_R2.fastq.gz"),
		up2 = temp("results/TRIMM/{sample}_trimm_unpair_R2.fastq.gz"),
	params:
		threads = "7",
	log:
		"results/logs/02_TRIMM_PE_{sample}.log",
	benchmark:
	        "results/benchmarks/02_TRIMM_PE_{sample}.txt"
	conda:
		 os.path.join(DIR, "envs/QC.yml"),
	message:
		"trim adaptors",
	shell:
		(" trimmomatic PE \
		-threads {params.threads} \
		-phred33 \
		{input.f1} {input.f2} \
		{output.p1} {output.up1} \
		{output.p2} {output.up2} \
		ILLUMINACLIP:{input.fasta}:2:30:10 2> {log} ")

rule CAT_pair_unpair:
	input:
		p1  = "results/TRIMM/{sample}_trimm_pair_R1.fastq.gz",
		up1 = "results/TRIMM/{sample}_trimm_unpair_R1.fastq.gz",
		p2  = "results/TRIMM/{sample}_trimm_pair_R2.fastq.gz",
		up2 = "results/TRIMM/{sample}_trimm_unpair_R2.fastq.gz",
	output:
		c1 = temp("results/TRIMM/{sample}_trimm_cat_R1.fastq.gz"),
		c2 = temp("results/TRIMM/{sample}_trimm_cat_R2.fastq.gz"),
	log:
		"results/logs/02_CAT_pair_{sample}.log",
	benchmark:
	    "results/benchmarks/02_CAT_pair_unpair_{sample}.txt"
	message:
		"concatenate pair and unpair",
	shell:
		(" cat {input.p1} {input.up1} > {output.c1} 2> {log};\
		cat {input.p2} {input.up2} > {output.c2} 2>> {log} ")

rule TRIMM_fastqc:
	input:
		expand(["results/TRIMM/{sample}_trimm_cat_R1.fastq.gz", "results/TRIMM/{sample}_trimm_cat_R2.fastq.gz"], sample = SAMPLES),
	output:
		html = expand(["results/TRIMM/{sample}_trimm_cat_R1_fastqc.html","results/TRIMM/{sample}_trimm_cat_R2_fastqc.html"], sample = SAMPLES),
	params:
		threads = "7",
	log:
		expand(["results/logs/02_TRIMM_fastqc_{sample}.log"], sample = SAMPLES)
	conda:
		 os.path.join(DIR, "envs/QC.yml"),
	message:
		"FASTQC -  trimm fastq files"
	shell:
		("fastqc {input} -t {params.threads} --outdir results/TRIMM/ 2> {log} ")

rule TRIMM_multiQC:
	input:
		expand(["results/TRIMM/{sample}_trimm_cat_R1_fastqc.html","results/TRIMM/{sample}_trimm_cat_R2_fastqc.html"], sample = SAMPLES)
	output:
		"results/multiQC_rep/trim_adaptors_multiqc_report.html"
	log:
		expand(["results/logs/02_TRIMM_multiQC_{sample}.log"], sample = SAMPLES)
	conda:
		 os.path.join(DIR, "envs/QC.yml"),
	message:
		" MultiQC report - after adaptor "
	shell:
		(" multiqc results/TRIMM/ -o results/multiQC_rep/ -i trim_adaptors ")

