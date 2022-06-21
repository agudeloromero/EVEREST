"""
workflow 06_Taxonomy_nt.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd() 

rule MMSEQ_eTaxonomy_nt:
	input:
		fasta  = "results/CHECKV/{sample}/viruses_rename.fna",
		DB   = config["MMSEQ_Viral_DB_nt"],
	output:
		lca    = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_lca.tsv",
		report = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_report",
		aln    = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_tophit_aln",
		tophit = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_tophit_report",
	params:
		type = "easy-taxonomy",
		prefix = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt",
		tmp  = "results/MMSEQ_eTaxonomy/{sample}_tmp",
		other= "--min-length 100 -a --tax-lineage 2 --search-type 2 -e 1e-20",
		lca  = "--lca-ranks superkingdom,phylum,class,order,family,genus,species",
		sen  = " --start-sens 2 -s 7 --sens-steps 3",
		out  = "--format-output query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage",
		threads = "--threads 7",
	log:
		temp("results/logs/06_MMSEQ_eTaxonomy_nt_{sample}.log"),
	benchmark:
		"results/benchmarks/06_MMSEQ_eTaxonomy_nt_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/MMSEQS.yml"),
	message:
		"Taxonomy with nt",
	shell:
		(" mmseqs {params.type} {input.fasta} {input.DB} {params.prefix} {params.tmp} {params.other} {params.lca} {params.sen} {params.out} {params.threads} 2> {log} ")

rule MMSEQ_eTaxonomy_nt_aln_header:
	input:
		aln1 = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_tophit_aln",
	output:
		aln2 = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_tophit_aln.txt",
	log:
		temp("results/logs/06_MMSEQ_eTaxonomy_nt_aln_header_{sample}.log"),
	benchmark:
		"results/benchmarks/06_MMSEQ_eTaxonomy_nt_aln_header_{sample}.txt"
	message:
		"add header to taxonomy nt"
	shell:
		(" sed '1i aln_query\taln_target\taln_evalue\taln_pident\taln_fident\taln_nident\taln_mismatch\taln_qcov\taln_tcov\taln_qstart\taln_qend\taln_qlen\taln_tstart\taln_tend\taln_tlen\taln_alnlen\taln_bits\taln_qheader\taln_theader\taln_taxid\taln_taxname\taln_taxlineage' {input.aln1} > {output.aln2} 2> {log} ")

rule TAXONKIT_nt_reformated:
	input:
		tsv = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_lca.tsv",
		DB  = config["TAX_nt"],
	output:
		tsv = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_lca_reformated.tsv",
	log:
		temp("results/logs/06_MMSEQ_eTaxonomy_nt_lca_reformated_{sample}.log"),
	benchmark:
		"results/benchmarks/06_MMSEQ_eTaxonomy_nt_lca_reformated_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/taxonkit.yml"),
	message:
		  "taxonomy nt",
	shell:
		(" taxonkit lineage {input.tsv} --data-dir {input.DB} -i 2 | taxonkit reformat --data-dir {input.DB} -i 7 -f '{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}' -F --fill-miss-rank | cut --complement -f5,6 > {output.tsv} 2> {log}  ")

rule MMSEQ_eTaxonomy_nt_lca_header:
	input:
		tsv = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_lca_reformated.tsv",
	output:
		tsv = "results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_lca_reformated_header.tsv",
	log:
		temp("results/logs/06_MMSEQ_eTaxonomy_nt_lca_header_{sample}.log"),
	benchmark:
		"results/benchmarks/06_MMSEQ_eTaxonomy_nt_lca_header_{sample}.txt"
	message:
		"add header to taxonomy nt",
	shell:
		(" sed '1 i\lca_query\tlca_taxid\tlca_taxonomic_rank\tlca_taxonomic_name\tlca_taxlineage\tlca_kingdom\tlca_phylum\tlca_class\tlca_order\tlca_family\tlca_genus\tlca_species' {input.tsv} > {output.tsv} 2> {log}  ")

checkpoint R_nt_summary_single:
	input:
		"results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_lca_reformated_header.tsv","results/MMSEQ_eTaxonomy/{sample}_tax_viral_nt_tophit_aln.txt",config["BALTIMORE"]
	output:
		"results/Summary/{sample}_nt_summary_mmseqs2.txt"
	log:
		temp("results/logs/06_R_nt_summary_{sample}.log"),
	benchmark:
		"results/benchmarks/06_R_nt_summary_{sample}.txt"
	conda:
		os.path.join(DIR, "envs/R.yml"),
	message:
		  "compiling summary - nt",
	script:
		os.path.join(DIR, "scr/Summary_script.R")

rule R_nt_summary_all:
	output:
		"results/Summary/Summary_nt_mmseqs2.txt",
	params:
		"results/Summary","_nt_"
	conda:
		os.path.join(DIR, "envs/R.yml"),
	message:
		"Combine Samples for summary - nt"
	script:
		os.path.join(DIR, "scr/Summary_CombineSamples_script.R")

