
# **Snakemake Installation**

## **Via Conda/Mamba**

In case you do not have mamba in your Mambaforge
```
$ conda install -n base -c conda-forge mamba -y
```

## **Full Installation**
```
$ conda activate base
$ conda create -c anaconda -n everest python=3.9.12 -y
$ mamba -c conda-forge -c bioconda snakemake=7.3.6 -y
```

## **Activation of the environment**
```
$ conda activate everest
```
