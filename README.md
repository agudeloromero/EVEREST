# EVEREST (pipEline for Viral assEmbly and chaRactEriSaTion)

The current version of EVEREST is operational and still under development.

---
**Running EVEREST:**

EVEREST is a snakemake pipeline for virus discovery, its main purpose is to characterise phage genomes isolates but can be also used for all the virome.
1. Clone EVEREST repository.
```
git clone --recursive https://github.com/agudeloromero/EVEREST.git
```

2. Create conda environment using everest.yml file.
```
conda env -n everest create -f everest.yml
```

3. EVEREST databases have to be requesting by email, but they will open soon.

**How to run**

Edit conf/config.yml file to point to the input, output and databases directories. 

Input directory should contain the .fastq.gz files to analyse ( -j option have to be set depending on number of available cores).
```
snakemake --use-conda -k -j 2 --keep-going
```

Files are expected to be named as "name_R1" and "name_R2" plus extension. In case you need to rename then, please see this example (https://github.com/agudeloromero/rename_fastq_files).

