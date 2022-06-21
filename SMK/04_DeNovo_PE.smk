"""
workflow 04_DeNovo.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

rule BBMAP_merge:
	input:
		f1  = "results/FASTQ/{sample}_unmapped_cat_dedup_norm_R1.fastq.gz",
		f2  = "results/FASTQ/{sample}_unmapped_cat_dedup_norm_R2.fastq.gz",
	output:
		o   = temp("results/FASTQ/{sample}_unmapped_cat_R1_merge.fastq.gz"),
		ou1 = temp("results/FASTQ/{sample}_unmapped_cat_unmerge_R1.fastq.gz"),
		ou2 = temp("results/FASTQ/{sample}_unmapped_cat_unmerge_R2.fastq.gz"),
	params:
		"minoverlap=20 strict=t"
	log:
		"results/logs/04_BBMAP_merge_{sample}.log",
	benchmark:
		"results/benchmarks/04_BBMAP_merge_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"merge reads R1 + R2",
	shell:
		(" bbmerge.sh -Xmx20000m in1={input.f1} in2={input.f2} \
		out={output.o} outu={output.ou1} outu2={output.ou2} {params} 2> {log} ")

rule TRIMM_unmerge:
	input:
		f1  = "results/FASTQ/{sample}_unmapped_cat_unmerge_R1.fastq.gz",
		f2  = "results/FASTQ/{sample}_unmapped_cat_unmerge_R2.fastq.gz",
		fasta = config["adaptor"],
	output:
		p1  = temp("results/FASTQ/{sample}_unmapped_cat_unmerge_pair_R1.fastq.gz"),
		up1 = temp("results/FASTQ/{sample}_unmapped_cat_unmerge_unpair_R1.fastq.gz"),
		p2  = temp("results/FASTQ/{sample}_unmapped_cat_unmerge_pair_R2.fastq.gz"),
		up2 = temp("results/FASTQ/{sample}_unmapped_cat_unmerge_unpair_R2.fastq.gz"),
	params:
		threads = "5",
		others  = "LEADING:10 TRAILING:10 SLIDINGWINDOW:3:15 MINLEN:50",
	log:
		"results/logs/04_TRIMM_unmerge_{sample}.log",
	benchmark:
		"results/benchmarks/04_TRIMM_unmerge_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/QC.yml"),
	message:
		"trimm unmerge",
	shell:
		(" trimmomatic PE -threads {params.threads} -phred33 {input.f1} {input.f2} \
		{output.p1} {output.up1} {output.p2} {output.up2} \
		ILLUMINACLIP:{input.fasta}:2:30:10 {params.others} 2> {log} ")

rule TRIMM_merge:
	input:
		f1 = "results/FASTQ/{sample}_unmapped_cat_R1_merge.fastq.gz",
		fasta = config["adaptor"],
	output:
		o = temp("results/FASTQ/{sample}_unmapped_cat_R1_merge_trimm.fastq.gz"),
	params:
		threads = "5",
		others  = "LEADING:10 TRAILING:10 SLIDINGWINDOW:3:15 MINLEN:50",
	log:
		"results/logs/04_TRIMM_merge_{sample}.log",
	benchmark:
		"results/benchmarks/04_TRIMM_merge_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/QC.yml"),
	message:
		"trimm merge",
	shell:
		(" trimmomatic SE -threads {params.threads} -phred33 {input.f1} \
		{output.o} ILLUMINACLIP:{input.fasta}:2:30:10 {params.others} 2> {log} ")

rule SPADES_DeNovo:
	input:
		m  = "results/FASTQ/{sample}_unmapped_cat_R1_merge_trimm.fastq.gz",
		f  = "results/FASTQ/{sample}_unmapped_cat_unmerge_pair_R1.fastq.gz",
		r  = "results/FASTQ/{sample}_unmapped_cat_unmerge_pair_R2.fastq.gz",
#		u1 = "results/FASTQ/{sample}_unmapped_cat_unmerge_unpair_R1.fastq.gz",
#		u2 = "results/FASTQ/{sample}_unmapped_cat_unmerge_unpair_R2.fastq.gz",
	output:
		fasta = "results/SPADES/{sample}/scaffolds.fasta",
	params:
		dir = "results/SPADES/{sample}",
		type = "--meta",
		threads = "-t 10",
		extra = "--only-assembler",
	log:
		"results/logs/04_SPADES_DeNovo_{sample}.log",
	benchmark:
		"results/benchmarks/04_SPADES_DeNovo_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/spades.yml"),
	message:
		"De novo assembly",
	shell:
		(" spades.py {params.type} --merged {input.m} \
		-1 {input.f} -2 {input.r} \
		-o {params.dir} \
		{params.threads} {params.extra} &> {log} ")
#--s1 {input.u1} --s2 {input.u2} \
