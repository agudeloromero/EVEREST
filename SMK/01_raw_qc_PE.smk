"""
workflow 01_raw_qc.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

#configfile: "config/config.yaml"
#SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))

#rule all:
#	input:
#		os.path.join(config["output_DIR"],"EVEREST/multiQC_rep/raw_qc_multiqc_report.html"),

rule FASTQC_raw:
	input:
		expand([os.path.join(config["input_DIR"], "{sample}_R1.fastq.gz"), os.path.join(config["input_DIR"], "{sample}_R2.fastq.gz")], sample = SAMPLES),
	output:
		html = expand([os.path.join(config["output_DIR"],"EVEREST/raw_fastqc/{sample}_R1_fastqc.html"),os.path.join(config["output_DIR"],"EVEREST/raw_fastqc/{sample}_R2_fastqc.html")], sample = SAMPLES),
	params:
		threads = "7",
		in_dir = os.path.join(config["output_DIR"], "EVEREST/raw_fastqc"),
	log:
		expand([os.path.join(config["output_DIR"],"EVEREST/logs/R01_S01_FASTQC_raw_{sample}.log")], sample = SAMPLES),
	conda:
		os.path.join(DIR, "envs/QC.yml"),
	message:
		"FASTQC - raw fastq files"
	shell:
		(" fastqc {input} -t {params.threads} --outdir {params.in_dir} 2> {log} ")

rule multiQC_raw:
	input:
		expand([os.path.join(config["output_DIR"],"EVEREST/raw_fastqc/{sample}_R1_fastqc.html"), os.path.join(config["output_DIR"],"EVEREST/raw_fastqc/{sample}_R2_fastqc.html")], sample = SAMPLES),
	output:
		os.path.join(config["output_DIR"],"EVEREST/multiQC_rep/raw_qc_multiqc_report.html"),
	params:
		in_dir = os.path.join(config["output_DIR"],"EVEREST/raw_fastqc"),
	log:
	    os.path.join(config["output_DIR"],"EVEREST/logs/R01_S02_multiQC_raw.log"),
	conda:
		os.path.join(DIR, "envs/QC.yml"),	
	message:
		"Raw multiQC"
	shell:
		(" multiqc {params.in_dir} -n {output} -i raw_qc 2> {log} ")
