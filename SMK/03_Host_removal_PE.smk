"""
workflow 03_Host_removal.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd()

configfile: "config/config.yaml"
#SAMPLES, = glob_wildcards(os.path.join(config["input_DIR"],"{sample}_R1.fastq.gz"))
#MINIMAP, = glob_wildcards(os.path.join(config["output_DIR"], "EVEREST/TRIMM/{minimap}_trimm_cat_R1.fastq.gz"))
#DUPLI,   = glob_wildcards(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{duplis}_unmapped_cat_dedup.fastq.gz"))
#NORM,   =  glob_wildcards(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{norms}_unmapped_cat_dedup_norm_R1.fastq.gz"))

#rule all:
#	input:


rule MINIMAP2_index:
	input:
		fasta = config["genome"],
	output:
		index = os.path.join(config["output_DIR"],"index_minimap2/index"),
	params:
		threads = "7",
	log:
		os.path.join(config["output_DIR"],"EVEREST/logs/R03_S01_MINIMAP2_index.log"),
	benchmark: 
		os.path.join(config["output_DIR"],"EVEREST/benchmarks/R03_S01_MINIMAP2_index.txt"),
	conda:
		os.path.join(DIR, "envs/minimap2.yml"),
	message:
		"index minimap2",
	shell:
		(" minimap2 -t {params.threads} -d {output.index} {input.fasta} 2> {log} ")

rule MINIMAP2_host_removal:
	input:
		p1  = os.path.join(config["output_DIR"], "EVEREST/TRIMM/{sample}_trimm_cat_R1.fastq.gz"),
		p2  = os.path.join(config["output_DIR"], "EVEREST/TRIMM/{sample}_trimm_cat_R2.fastq.gz"),
		inx = os.path.join(config["output_DIR"], "index_minimap2/index"),
	output:
		o1 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_R1.fastq"),
		o2 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_R2.fastq"),
		s  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_singletons.fastq"),
	params:
		type_m = "-ax sr --secondary=no", 
		threads = "7",
		type_u = "view -f 4 -h",
		type_s = "sort -@ 7",
		type_f = "fastq -NO -@ 7",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S02_MINIMAP2_host_removal_{sample}.log"),
	benchmark: 
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R03_S02_MINIMAP2_host_removal_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/minimap2.yml"),
	message:
		"obtain unmapped from minimap2 and convert to fastq",
	shell:
		(" ( minimap2 {params.type_m} -t {params.threads} {input.inx} {input.p1} {input.p2} \
		| samtools {params.type_u} - | samtools {params.type_s} \
		| samtools {params.type_f} - -1 {output.o1} -2 {output.o2} -s {output.s} ) &> {log} ")  

rule BBMAP_singletons_PE:
	input:
		s  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_singletons.fastq"),
	output:
		s1 = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_singletons_R1.fastq")),
		s2 = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_singletons_R2.fastq")),
	params:
		mem = "-Xmx20000m", 
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S03_BBMAP_singletons_PE_{sample}.log"),
	benchmark: 
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R03_S03_BBMAP_singletons_PE_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"BBMAP for rescue the PE from singletons",
	shell:
		(" reformat.sh {params.mem} in={input.s} out={output.s1} out2={output.s2} 2> {log} ")

rule CAT_PE:
	input:
		o1 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_R1.fastq"),
		o2 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_R2.fastq"),
		s1 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_singletons_R1.fastq"),
		s2 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_singletons_R2.fastq"),
	output:
		f1 = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R1.fastq")),
		f2 = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R2.fastq")),
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S04_CAT_PE_{sample}.log"),
	benchmark: 
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R03_S04_CAT_PE_{sample}.txt"),
	message:
		" Concatenate: pairs and singletons pairs from unmapped reads ",
	shell:
		(" cat {input.o1} {input.s1} > {output.f1}; \
		 cat {input.o2} {input.s2} > {output.f2} ")

rule PIGZ_fastq:
	input:
		expand([os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R1.fastq"),os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R2.fastq")], sample=SAMPLES),
	output:
		temp(expand([os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R1.fastq.gz"), os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R2.fastq.gz")], sample=SAMPLES)),
	params:
		"-p 7 -5",
	log:
		expand([os.path.join(config["output_DIR"], "EVEREST/logs/R03_S05_PIGZ_fastq_{sample}.log")], sample=SAMPLES),
	conda:
		os.path.join(DIR, "envs/minimap2.yml"),
	message:
		"compress fastq files",
	shell:
		(" pigz {params} {input} ")

rule BBMAP_dup:
	input:
		f1 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R1.fastq.gz"),
		f2 = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_R2.fastq.gz"),
	output:
		f  = temp(os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup.fastq.gz")),
	params:
		mem = "-Xmx20000m",
		other = "ac=f s=10 e=7 minidentity=92",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S06_BBMAP_dup_{sample}.log"),
	benchmark: 
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R03_S06_BBMAP_dup_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"BBMAP deduplication",
	shell:
		(" dedupe.sh {params.mem} {params.other} in1={input.f1} in2={input.f2} out={output.f} 2> {log} ")

rule BBMAP_deduped_reformat:
	input:
		f = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup.fastq.gz"),
	output:
		f1  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_R1.fastq.gz"),
		f2  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_R2.fastq.gz"),
	params:
		mem = "-Xmx20000m",
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S07_BBMAP_dup_ref_{sample}.log"),
	benchmark: 
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/S03_S07_BBMAP_dup_ref_{sample}.txt"),
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"BBMAP deduped reformat",
	shell:
		(" reformat.sh {params.mem} in={input.f} out={output.f1} out2={output.f2} 2> {log} ")

rule BBMAP_duduped_normalisation:
	input:
		f1  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_R1.fastq.gz"),
		f2  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_R2.fastq.gz"),
	output:
		f1  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R1.fastq.gz"),
		f2  = os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R2.fastq.gz"),	
	log:
		os.path.join(config["output_DIR"], "EVEREST/logs/R03_S08_BBMAP_dup_norm_{sample}.log"),
	benchmark: 
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R03_S08_BBMAP_dup_norm_{sample}.txt"),
	params:
		mem = "-Xmx20000m",
	conda:
		os.path.join(DIR, "envs/BBMAP.yml"),
	message:
		"BBMAP digital normalisation",
	shell:
		(" bbnorm.sh {params.mem} in={input.f1} in2={input.f2} out={output.f1} out2={output.f2} 2> {log} ")

rule FASTQC_before_merge:
	input:
		expand([os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R1.fastq.gz"), os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R2.fastq.gz")], sample=SAMPLES),
	output:
		html = expand([os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R1_fastqc.html"), os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R2_fastqc.html")], sample=SAMPLES),
	params:
		threads = "7",
		out_dir = os.path.join(config["output_DIR"], "EVEREST/FASTQ"),
	log:
		expand([os.path.join(config["output_DIR"], "EVEREST/logs/R03_S09_FASTQC_before_merge_{sample}.log")], sample=SAMPLES),
	conda:
		os.path.join(DIR, "envs/QC.yml"),
	message:
		"FASTQC - fastq files before merge",
	shell:
		(" fastqc {input} -t {params.threads} --outdir {params.out_dir}  2> {log} ")
	
rule multiQC_before_merge:
	input:
		expand([os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R1_fastqc.html"), os.path.join(config["output_DIR"], "EVEREST/FASTQ/{sample}_unmapped_cat_dedup_norm_R2_fastqc.html")], sample=SAMPLES),
	output:
		os.path.join(config["output_DIR"],"EVEREST/multiQC_rep/fastq_before_merge_multiqc_report.html"),
	params:
		in_dir = os.path.join(config["output_DIR"],"EVEREST/FASTQ"),
	log:
		os.path.join(config["output_DIR"],"EVEREST/logs/R03_S10_multiQC_before_merge.log"),
	conda:
		os.path.join(DIR, "envs/QC.yml"),
	message:
		"MultiQC report - fastq files before merge",
	shell:
		(" multiqc {params.in_dir} -n {output} -i fastq_before_merge 2> {log} ")
