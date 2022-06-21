"""
workflow 06_Taxonomy_aa.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd() 

rule MMSEQ_eTaxonomy_aa:
	input:
		fasta = "results/CHECKV/{sample}/viruses_rename.fna",
		DB    = config["MMSEQ_Viral_DB_aa"],
	output:
		lca    = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_lca.tsv",
		report = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_report",
		aln    = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_tophit_aln",
		tophit = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_tophit_report",
	params:
		type = "easy-taxonomy",
		prefix = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa",
		tmp  = "results/MMSEQ_eTaxonomy/{sample}_tmp",
		other= "--min-length 30 -a --tax-lineage 1 --search-type 2 -e 1e-5 --majority 0.5 --vote-mode 1",
		lca  = "--lca-ranks superkingdom,phylum,class,order,family,genus,species",
		sen  = "--start-sens 1 --sens-steps 3 -s 7 --lca-mode 3 --shuffle 0",
		out  = "--format-output query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage",
		threads = "--threads 7",
	log:
		temp("results/logs/06_MMSEQ_eTaxonomy_aa_{sample}.log"),
	benchmark:
		"results/benchmarks/06_MMSEQ_eTaxonomy_aa_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/MMSEQS.yml"),
	message:
		"taxonomy aa",
	shell:
		(" mmseqs {params.type} {input.fasta} {input.DB} {params.prefix} {params.tmp} {params.other} {params.lca} {params.sen} {params.out} {params.threads} 2> {log} ")

rule MMSEQ_eTaxonomy_aa_aln_header:
	input:
		aln1 = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_tophit_aln",
	output:
		aln2 = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_tophit_aln.txt",
	log:
		temp("results/logs/06_MMSEQ_eTaxonomy_aa_aln_header_{sample}.log"),
	benchmark:
		"results/benchmarks/06_MMSEQ_eTaxonomy_aa_aln_header_{sample}.txt"
	message:
		"add header to taxonomy aa"
	shell:
		(" sed '1i aln_query\taln_target\taln_evalue\taln_pident\taln_fident\taln_nident\taln_mismatch\taln_qcov\taln_tcov\taln_qstart\taln_qend\taln_qlen\taln_tstart\taln_tend\taln_tlen\taln_alnlen\taln_bits\taln_qheader\taln_theader\taln_taxid\taln_taxname\taln_taxlineage' {input.aln1} > {output.aln2} 2> {log} ")

rule TAXONKIT_aa_reformated:
	input:
		tsv = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_lca.tsv",
		DB  = config["TAX_aa"],
	output:
		tsv = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_lca_reformated.tsv",
	log:
		temp("results/logs/06_MMSEQ_eTaxonomy_aa_lca_reformated_{sample}.log"),
	benchmark:
		"results/benchmarks/06_MMSEQ_eTaxonomy_aa_lca_reformated_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/taxonkit.yml"),
	message:
		  "taxonomy aa",
	shell:
		(" taxonkit lineage {input.tsv} --data-dir {input.DB} -i 2 | taxonkit reformat --data-dir {input.DB} -i 11 -f '{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}' -F --fill-miss-rank | cut --complement -f5,6,7,8,9,10 > {output.tsv} 2> {log}  ")

rule MMSEQ_eTaxonomy_aa_lca_header:
	input:
		tsv = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_lca_reformated.tsv",
	output:
		tsv = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_lca_reformated_header.tsv",
	log:
		temp("results/logs/06_MMSEQ_eTaxonomy_aa_lca_header_{sample}.log"),
	benchmark:
		"results/benchmarks/06_MMSEQ_eTaxonomy_aa_lca_header_{sample}.txt"
	message:
		"add header to taxonomy aa",
	shell:
		(" sed '1 i\lca_query\tlca_taxid\tlca_taxonomic_rank\tlca_taxonomic_name\tlca_taxlineage\tlca_kingdom\tlca_phylum\tlca_class\tlca_order\tlca_family\tlca_genus\tlca_species' {input.tsv} > {output.tsv} 2> {log} ")

checkpoint R_aa_summary_single:
	input:
		"results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_lca_reformated_header.tsv","results/MMSEQ_eTaxonomy/{sample}_tax_viral_aa_tophit_aln.txt",config["BALTIMORE"]
	output:
		"results/Summary/{sample}_aa_summary_mmseqs2.txt"
	log:
		temp("results/logs/06_R_aa_summary_{sample}.log"),
	benchmark:
		"results/benchmarks/06_R_aa_summary_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/R.yml"),
	message:
		  "compiling summary - aa",
	script:
		os.path.join(DIR, "scr/Summary_script.R")

rule R_aa_summary_all:
	output:
		"results/Summary/Summary_aa_mmseqs2.txt",
	params:
		"results/Summary","_aa_"
	conda:
		os.path.join(DIR, "envs/R.yml"),
	message:
		"Combine Samples for summary - aa"
	script:
		os.path.join(DIR, "scr/Summary_CombineSamples_script.R")

rule R_summary_all:
	input:
		"results/Summary/Summary_aa_mmseqs2.txt","results/Summary/Summary_nt_mmseqs2.txt",
	output:
		"results/Summary/Summary_mmseqs2.txt",
	conda:
		os.path.join(DIR, "envs/R.yml"),
	message:
		"Combine Samples for summary "
	script:
		os.path.join(DIR, "scr/Summary_CombineSamples_all_script.R")

