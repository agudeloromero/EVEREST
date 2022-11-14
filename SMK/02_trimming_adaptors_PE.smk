"""
workflow 02_trimming_adaptors.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

#configfile: "config/config.yaml"
#SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))

#rule all:
#	input:
#		os.path.join(config["output_DIR"], "EVEREST/multiQC_rep/trim_adaptors_multiqc_report.html"),

rule BBMAP_phix:
	input:
		f1 = os.path.join(config["input_DIR"], "{sample}_R1.fastq.gz"),
		f2 = os.path.join(config["input_DIR"], "{sample}_R2.fastq.gz"),
#		phix = config["PHIX"],
	output:
		clean1 = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_clean_R1.fastq.gz"),
		clean2 = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_clean_R2.fastq.gz"),
		unclean1 = temp(os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_noclean_R1.fastq.gz")),
		unclean2 = temp(os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_noclean_R2.fastq.gz")),
		stats_phix = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_stats_phix.txt"),
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R02_S01_BBMAP_phix_{sample}.log"),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R02_S01_BBMAP_phix_{sample}.txt"),
	params:
		mem = "-Xmx20000m",
		extra = "ordered cardinality k=31 hdist=1",
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"BBMAP PhiX screening and removal",
	shell:
		(" bbduk.sh {params.mem} in1={input.f1} in2={input.f2} out1={output.clean1} out2={output.clean2} \
		outm1={output.unclean1} outm2={output.unclean2} ref=artifacts,phix stats={output.stats_phix} {params.extra} 2> {log} ")

rule TRIMM_PE:
	input:
		f1 = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_clean_R1.fastq.gz"),
		f2 = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_clean_R2.fastq.gz"),
		fasta = config["adaptor"],
	output:
		p1  = temp(os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_pair_R1.fastq.gz")),
		up1 = temp(os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_unpair_R1.fastq.gz")),
		p2  = temp(os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_pair_R2.fastq.gz")),
		up2 = temp(os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_unpair_R2.fastq.gz")),
	params:
		threads = "7",
	log:
		os.path.join(config["output_DIR"],"EVEREST/logs/R02_S02_TRIMM_PE_{sample}.log"),
	benchmark:
	    os.path.join(config["output_DIR"],"EVEREST/benchmarks/R02_S02_TRIMM_PE_{sample}.txt"),
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
		p1  = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_pair_R1.fastq.gz"),
		up1 = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_unpair_R1.fastq.gz"),
		p2  = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_pair_R2.fastq.gz"),
		up2 = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_unpair_R2.fastq.gz"),
	output:
		c1 = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_cat_R1.fastq.gz"),
		c2 = os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_cat_R2.fastq.gz"),
	log:
		os.path.join(config["output_DIR"],"EVEREST/logs/R02_S03_CAT_pair_{sample}.log"),
	benchmark:
	    os.path.join(config["output_DIR"],"EVEREST/benchmarks/R02_S03_CAT_pair_unpair_{sample}.txt"),
	message:
		"concatenate pair and unpair",
	shell:
		(" cat {input.p1} {input.up1} > {output.c1} ;\
		cat {input.p2} {input.up2} > {output.c2} ")

rule FASTQC_trimm:
	input:
		expand([os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_cat_R1.fastq.gz"), os.path.join(config["output_DIR"], "EVEREST/TRIMM/{sample}_trimm_cat_R2.fastq.gz")], sample = SAMPLES),
	output:
		html = expand([os.path.join(config["output_DIR"],"EVEREST/TRIMM/{sample}_trimm_cat_R1_fastqc.html"), os.path.join(config["output_DIR"], "EVEREST/TRIMM/{sample}_trimm_cat_R2_fastqc.html")], sample = SAMPLES),
	params:
		threads = "7",
		in_dir =  os.path.join(config["output_DIR"], "EVEREST/TRIMM"),
		out_dir = os.path.join(config["output_DIR"], "EVEREST/TRIMM"),
	log:
		expand([os.path.join(config["output_DIR"], "EVEREST/logs/R02_S04_TRIMM_fastqc_{sample}.log")], sample = SAMPLES),
	conda:
		 os.path.join(DIR, "envs/QC.yml"),
	message:
		" FASTQC -  trimm fastq files "
	shell:
		(" fastqc {input} -t {params.threads} --outdir {params.in_dir} 2> {log} ")

rule multiQC_trimm:
	input:
		expand([os.path.join(config["output_DIR"], "EVEREST/TRIMM/{sample}_trimm_cat_R1_fastqc.html"), os.path.join(config["output_DIR"], "EVEREST/TRIMM/{sample}_trimm_cat_R2_fastqc.html")], sample = SAMPLES),
	output:
		os.path.join(config["output_DIR"], "EVEREST/multiQC_rep/trim_adaptors_multiqc_report.html"),
	params:
		in_dir = os.path.join(config["output_DIR"], "EVEREST/TRIMM"),
	log:
		os.path.join(config["output_DIR"],"EVEREST/logs/R02_S05_TRIMM_multiQC.log"),
	conda:
		 os.path.join(DIR, "envs/QC.yml"),
	message:
		" MultiQC report - after adaptor "
	shell:
		(" multiqc {params.in_dir} -n {output} -i trim_adaptors 2> {log} ")

