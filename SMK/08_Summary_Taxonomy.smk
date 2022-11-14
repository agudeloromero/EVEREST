"""
workflow 08_Summary_Taxonomy.smk
Author: Patricia Agudelo-Romero
email : Patricia.AgudeloRomero@telethonkids.org.au
"""

import os
DIR = os.getcwd() 

#configfile: "config/config.yaml"
#NTSUMM, = glob_wildcards(os.path.join(config["output_DIR"], "EVEREST/Summary/{ntsumm}_nt_summary_mmseqs2.txt"))
#AASUMM, = glob_wildcards(os.path.join(config["output_DIR"], "EVEREST/Summary/{aasumm}_aa_summary_mmseqs2.txt"))

#rule all:
#	input:
#		os.path.join(config["output_DIR"], "EVEREST/Summary/Summary_mmseqs2.txt"),

rule R_summary_all:
	input:
		aa_dir = os.path.join(config["output_DIR"], "EVEREST/Summary/Summary_aa_mmseqs2.txt"),
		nt_dir = os.path.join(config["output_DIR"], "EVEREST/Summary/Summary_nt_mmseqs2.txt"),
	output:
		os.path.join(config["output_DIR"], "EVEREST/Summary/Summary_mmseqs2.txt"),
	log:
		temp(os.path.join(config["output_DIR"], "EVEREST/logs/R08_S03_R_summary_all.log")),
	benchmark:
		os.path.join(config["output_DIR"], "EVEREST/benchmarks/R08_S03_R_summary_all.txt"),
	conda:
		os.path.join(DIR, "envs/R.yml"),
	message:
		"Combine Samples for summary "
	script:
		os.path.join(DIR, "scr/Summary_CombineSamples_all_script.R")

