"""
workflow 05_Cleaning_contigs.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd() 

rule MMSEQ2_eLinclust:
	input:
		fasta = "results/SPADES/{sample}/scaffolds.fasta",
	output:
		all  = "results/MMSEQ_eLinclust/{sample}_all_seqs.fasta",
		clus = "results/MMSEQ_eLinclust/{sample}_cluster.tsv",
		rep  = "results/MMSEQ_eLinclust/{sample}_rep_seq.fasta",
	params:
		prefix  = "results/MMSEQ_eLinclust/{sample}",
		tmp_dir = "results/MMSEQ_eLinclust/{sample}_tmp",
		threads = "7",
		type    = "--min-seq-id 0.98 --kmer-per-seq-scale 0.3 --sort-results 1 --alignment-mode 3 --cov-mode 1",
	log:
		"results/logs/05_MMSEQ2_eLinclust_{sample}.log",
	benchmark:
		"results/benchmarks/05_MMSEQ2_eLinclust_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/MMSEQS.yml"),
	message:
		"cluster sequences"
	shell:
		(" mmseqs easy-linclust {input.fasta} {params.prefix} {params.tmp_dir} --threads {params.threads} {params.type} 2> {log} ")

rule SEQKIT_filter:
	input:
		fasta  = "results/MMSEQ_eLinclust/{sample}_rep_seq.fasta",
	output:
		filter = "results/MMSEQ_eLinclust/{sample}_rep_seq_FilterLen.fasta",
	params:
		"seq -m 5000",
	log:
		"results/logs/05_SEQKIT_filter_{sample}.log",
	benchmark:
		"results/benchmarks/05_SEQKIT_filter_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/seqkit.yml"),
	message:
		"remove contigs < 5000 bp",
	shell:
		(" seqkit {params} {input.fasta} -o {output.filter} 2> {log} ")

rule VIRSORTER_detect:
	input:
		fasta = "results/MMSEQ_eLinclust/{sample}_rep_seq_FilterLen.fasta",
		db = config["VIRSORTER_DB"],
	output:
		viral = "results/VIRSORTER/pass1/{sample}/final-viral-combined.fa",
		score = "results/VIRSORTER/pass1/{sample}/final-viral-score.tsv",
		boundary = "results/VIRSORTER/pass1/{sample}/final-viral-boundary.tsv",
	params:
		type = "run --keep-original-seq",
		dir  = "results/VIRSORTER/pass1/{sample}/",
		groups = "--include-groups dsDNAphage,NCLDV,RNA,ssDNA",
		others = "--min-length 5000 --min-score 0.5",
		threads = "-j 28 all",
	log:
		"results/logs/05_VIRSORTER_detect_{sample}.log",
	benchmark:
		"results/benchmarks/05_VIRSORTER_detect_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/virsorter2.yml"),
	message:
		"rid off from non-viral contigs",
	shell:
		(" virsorter {params.type} -i {input.fasta} -w {params.dir} --db-dir {input.db} {params.groups} {params.others} {params.threads} 2> {log} ")

rule CHECKV_viral_seq:
	input:
		fasta = "results/VIRSORTER/pass1/{sample}/final-viral-combined.fa",
		DB = config["CHECKV_DB"],
	output:
		fasta1 = "results/CHECKV/{sample}/viruses.fna",
	params:
		type = "end_to_end",
		dir  = "results/CHECKV/{sample}",
		threads = "-t 7",
	log:
		"results/logs/05_CHECKV_viral_seq_{sample}.log",
	benchmark:
		"results/benchmarks/05_CHECKV_viral_seq_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/checkv.yml"),
	message:
		"rid off non-viral reads",
	shell:
		(" checkv {params.type} -d {input.DB} {input.fasta} {params.dir} 2> {log} ")

rule RENAME_viral_seq:
	input:
		fasta = "results/CHECKV/{sample}/viruses.fna",
	output:
		fasta = "results/CHECKV/{sample}/viruses_rename.fna"
	log:
		"results/logs/05_RENAME_viral_seq_{sample}.log",
	benchmark:
		"results/benchmarks/05_RENAME_viral_seq_{sample}.txt"
	message:
		"rename viral contigs",
	shell:
		(" sed 's/||.*//' {input.fasta} > {output.fasta} 2> {log} ")
