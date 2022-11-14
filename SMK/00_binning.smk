"""
workflow 00_bining.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

configfile: "config/config.yaml"
SAMPLES, = glob_wildcards(os.path.join(config["output_DIR"],"EVEREST/CHECKV/{sample}/viruses_rename.fna"))
#SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))

rule all:
	input:
		expand(os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/fasta"), sample=SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/prodigal/viruses_rename_gene.fna"), sample=SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/prodigal/viruses_rename_protein.faa"), sample=SAMPLES),
		expand(os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/binnig"), sample=SAMPLES),

rule SEQKIT_split_fasta:
	input:
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"),
	output:
		dir = directory(os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/fasta")),
#	log:
#	benckmark:
	conda:
		os.path.join(DIR, "envs/seqkit.yml"),
	message:
		"splitting fasta for binning",
	shell:
		(" seqkit split2 --by-size 1 {input.fasta} -O {output.dir} ; \
		for f in {output.dir}/*.fna ; do mv $f {output.dir}/$( seqkit seq --name --only-id $f).fasta ; done ")

rule PRODIGAL_gene_protein:
	input:
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"),
	output:
		gene = os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/prodigal/viruses_rename_gene.fna"),
		protein = os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/prodigal/viruses_rename_protein.faa"),
#	log:
#	benckmark:
	conda:
		os.path.join(DIR, "envs/vrhyme.yml"),
	message:
		"creating genes and proteins from contigs",
	shell:
		(" prodigal -i {input.fasta} -o {output.gene} -a {output.protein} -p meta ")

rule vRHYME_binning:
	input:
		vRhyme = config["vRhyme"], 
		fastq = expand([os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"), os.path.join(config["input_DIR"],"{sample}_R2.fastq.gz")], sample=SAMPLES),
		fasta = os.path.join(config["output_DIR"], "EVEREST/CHECKV/{sample}/viruses_rename.fna"),
		gene = os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/prodigal/viruses_rename_gene.fna"),
		protein = os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/prodigal/viruses_rename_protein.faa"),
	output:
		dir = directory(os.path.join(config["output_DIR"], "EVEREST/vRhyme/{sample}/binnig")),
	params:
		"-t 4 --method longest",
##	log:
##	benckmark:
	conda:
		os.path.join(DIR, "envs/vrhyme.yml"),
	message:
		"vRhyme - binning  contigs",
	shell:
		(" {input.vRhyme} -i {input.fasta} -g {input.gene} -p {input.protein} -r {input.fastq} \
		-o {output.dir} {params} ")









