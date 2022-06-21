"""
workflow 01_raw_qc.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

#configfile: "config/config.yaml"
#DIRs = config["samples_DIR"]
#SAMPLES, = glob_wildcards(os.path.join(config["samples_DIR"],"{sample}_R1.fastq.gz"))

rule FASTQC_raw:
	input:
		expand([os.path.join(config["samples_DIR"], "{sample}_R1.fastq.gz"), os.path.join(config["samples_DIR"], "{sample}_R2.fastq.gz")], sample = SAMPLES),
	output:
		html = expand(["results/raw_fastqc/{sample}_R1_fastqc.html","results/raw_fastqc/{sample}_R2_fastqc.html"], sample = SAMPLES),
	params:
		threads = "7",
	log:
		expand(["results/logs/01_FASTQC_raw_{sample}.log"], sample = SAMPLES),
	conda:
		os.path.join(DIR, "envs/QC.yml"),
	message:
		"FASTQC - raw fastq files"
	shell:
		(" fastqc {input} -t {params.threads} --outdir results/raw_fastqc/ 2> {log} ")

rule multiQC_raw:
	input:
		expand(["results/raw_fastqc/{sample}_R1_fastqc.html","results/raw_fastqc/{sample}_R2_fastqc.html"], sample = SAMPLES),
	output:
		"results/multiQC_rep/raw_qc_multiqc_report.html"
#	log:
#	    expand(["results/logs/01_multiQC_raw_{sample}.log"], sample = SAMPLES)
	conda:
		os.path.join(DIR, "envs/QC.yml"),	
	message:
		"Raw multiQC"
	shell:
		(" multiqc results/raw_fastqc/ -o results/multiQC_rep -i raw_qc ")
