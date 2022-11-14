"""
workflow 04_DeNovo.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

#configfile: "config/config.yaml"
#SAMPLES, = glob_wildcards(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R1.fastq.gz"))

#rule all:
#	input:
#		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_all_seqs.fasta"), sample=SAMPLES),
#		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_cluster.tsv"), sample=SAMPLES),
#		expand(os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq.fasta"), sample=SAMPLES),
##		expand(os.path.join(config["output_DIR"], "EVEREST/SPADES/{sample}/scaffolds.fasta"), sample = SAMPLES),

rule BBMAP_merge:
	input:
		f1  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R1.fastq.gz"),
		f2  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R2.fastq.gz"),
	output:
		o   = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R1_merge.fastq.gz")),
		ou1 = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_R1.fastq.gz")),
		ou2 = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_R2.fastq.gz")),
	params:
		"minoverlap=20 strict=t"
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R04_S01_BBMAP_merge_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R04_S01_BBMAP_merge_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		" merge reads R1 + R2" ,
	shell:
		(" bbmerge.sh -Xmx20000m in1={input.f1} in2={input.f2} \
		out={output.o} outu={output.ou1} outu2={output.ou2} {params} 2> {log} ")

rule TRIMM_unmerge:
	input:
		f1  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_R1.fastq.gz"),
		f2  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_R2.fastq.gz"),
		fasta = config["adaptor"],
	output:
		p1  = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_pair_R1.fastq.gz")),
		up1 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_unpair_R1.fastq.gz"),
		p2  = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_pair_R2.fastq.gz")),
		up2 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_unpair_R2.fastq.gz"),
	params:
		threads = "5",
		others  = "LEADING:10 TRAILING:10 SLIDINGWINDOW:3:15 MINLEN:50",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R04_S02_TRIMM_unmerge_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R04_S02_TRIMM_unmerge_{sample}.txt"),
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
		f1 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R1_merge.fastq.gz"),
		fasta = config["adaptor"],
	output:
		o = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R1_merge_trimm.fastq.gz")),
	params:
		threads = "5",
		others  = "LEADING:10 TRAILING:10 SLIDINGWINDOW:3:15 MINLEN:50",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R04_S03_TRIMM_merge_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R04_S03_TRIMM_merge_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/QC.yml"),
	message:
		"trimm merge",
	shell:
		(" trimmomatic SE -threads {params.threads} -phred33 {input.f1} \
		{output.o} ILLUMINACLIP:{input.fasta}:2:30:10 {params.others} 2> {log} ")

checkpoint SPADES_DeNovo:
	input:
		m  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R1_merge_trimm.fastq.gz"),
		f  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_pair_R1.fastq.gz"),
		r  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_pair_R2.fastq.gz"),
		u1 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_unpair_R1.fastq.gz"),
		u2 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_unmerge_unpair_R2.fastq.gz"),
	output:
		fasta = os.path.join(config["output_DIR"], "EVEREST/SPADES/{sample}/scaffolds.fasta"),
		dir = directory(os.path.join(config["output_DIR"], "EVEREST/SPADES/{sample}")),
	params:
		type = "--meta",
		threads = "-t 10",
		extra = "--only-assembler",
#		type = "--careful",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R04_S04_SPADES_DeNovo_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R04_S04_SPADES_DeNovo_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/spades.yml"),
	message:
		"De novo assembly",
	shell:
		(" spades.py {params.type} --merge {input.m} \
		-1 {input.f} -2 {input.r} \
		-s {input.u1} -s {input.u2} \
		-o {output.dir} \
		{params.threads} {params.extra} &> {log} ")
# this doesn't work in --meta		--s1 {input.u1} --s2 {input.u2} \

def get_SPADES_files(wildcards):
	checkpoint_output = checkpoints.SPADES_DeNovo.get(**wildcards).output[1]
	print(checkpoint_output)
	return expand(os.path.join(checkpoint_output,"scaffolds.fasta"))

rule MMSEQ2_eLinclust:
	input:
		get_SPADES_files,
#       fasta = os.path.join(config["output_DIR"], "EVEREST/SPADES/{sample}/scaffolds.fasta"),
	output:
		all  = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_all_seqs.fasta"),
		clus = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_cluster.tsv"),
		rep  = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_rep_seq.fasta"),
	params:
		prefix  = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}"),
		tmp_dir = os.path.join(config["output_DIR"], "EVEREST/MMSEQ_eLinclust/{sample}_tmp"),
		threads = "7",
		type    = "--min-seq-id 0.98 --kmer-per-seq-scale 0.3 --sort-results 1 --alignment-mode 3 --cov-mode 1",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R05_S01_MMSEQ2_eLinclust_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R05_S01_MMSEQ2_eLinclust_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/MMSEQS.yml"),
	message:
		"cluster sequences"
	shell:
		(" mmseqs easy-linclust {input} {params.prefix} {params.tmp_dir} --threads {params.threads} {params.type} 2> {log} ")

