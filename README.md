# MITOAS (Mitochondrial Genome Assembly and Standardization Pipeline)

# Overview

## Novel mitochondrial genome assembly and standardization pipeline for paired-end short DNA reads




<code> <b>data/</b> </code> : Reference sequences for MitoS annotation and a csv file for standardizing strand orientation

<code> <b>bin/</b> </code> : Python scripts used in the snakemake pipeline

<code> <b>pipeline.smk</b> </code> : Snakemake file


# Setup

```
git clone https://github.com/bcgsc/mitogenome_assembly_pipeline.git
```


``` Dependencies ```

Python 3.9.1

Snakemake 5.30.2

Pandas 1.1.5

Numpy 1.19.4

BLAST 2.10.1+

Biopython 1.79

Bbmap 38.18

Seqtk 1.1-r91

Regex 2022.3.15 

AbySS 3.82

ABySS-Sealer (ABySS) 2.2.5

ntJoin 3.82

BWA 0.7.17-r1188

Samtools 1.10

Pilon 1.23

Mitos 2.0.8


---

```Special Installation Instructions for MitoS```

***Please note: mitos 2.0.8 uses python 2.7.18, if downloaded to the same environment with the other dependencies, a conflict will occur.***

Please install mitos in a new environment called "mitos" using the instruction below:

```
conda install mitos -c bioconda -m -n mitos
```


# Running MITOAS


### 1. Set up library folder and softlink paired-end read files

```Mandatory Directory Structure :```  
> <library_folder_name>
 >> r1.fq.gz  
 >> r2.fq.gz

***Note: Paired-end read files must be named "r1.fq.gz" and "r2.fq.gz"***


```
mkdir <library_folder_name>
cd <library_folder_name>

# Softlink Original Read Files to Library Folder
ln -sfT <target_read_1_fastq.gz> r1.fq.gz
ln -sfT <target_read_2_fastq.gz> r2.fq.gz
```

### 2.  Required Parameters for MITOAS

`library`: library folder name 

`k_sweeps`: k-mer sweeps used for ABySS 

`kc`: k-mer minimum coverage multiplicity cutoff for ABySS

`B`: bloom filter size for ABySS

`b`: bloom filter size for abyss-sealer (gap-filling step)

`P`: number of max alternate paths used in abyss-sealer (gap-filling step)

`blastdb_path`: path to the reference data folder for BLAST

`blastdb_name`: name for the reference database

`ref_path`: path to the fasta file containing reference sequences

`bf_end_recovery`: bloom filter size for abyss-sealer (end recovery step)

`P_end_recovery`: number of max alternate paths for abyss-sealer (end recovery step)

`mito_gencode`: csv file that contains the mitochondrial genetic code (e.g., 2, 5) for your target species (see below for the required format)



mito-genetic code can be searched on https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

More information on optimizing k and kc for ABySS: https://github.com/bcgsc/abyss#optimizing-the-parameters-k-and-kc

More information on abyss-sealer parameters (b, P, bf_end_recovery, P_end_recovery): https://github.com/bcgsc/abyss/tree/master/Sealer


`Required csv file format`
```
Scientific Name,Mitochondrial Genetic Code
Anoplopoma fimbria,2
Acipenser fulvescens,2
Sebastes zacentrus,2
```
`Build BLAST database using local mitogenome sequences`
```
makeblastdb -in <ref_path> -out <blastdb_name> -dbtype nucl -parse_seqids
```
`Examples`

Running the pipeline with a single set of parameters

```
snakemake -s pipeline.smk --config library=sea_otter_libraryID001 k_sweeps=96 kc_sweeps=4 B=30G b=30G P=5 blastdb_path=path/to/blastdb_folder blastdb_name=mito_blastdb ref_path=path/to/fasta bf_end_recovery=30G P_end_recovery=5 mito_gencode=mito_gencode.csv --cores 10
```


Running the pipeline for 2 libraries using different k/kc sweeps 

```
snakemake -s pipeline.smk --config library=[sea_otter_libraryID001,sea_otter_libraryID002]  k_sweeps=[64, 80, 96] kc_sweeps=[2,3] B=30G b=30G P=5 blastdb_path=path/to/blastdb_folder blastdb_name=mito_blastdb ref_path=path/to/fasta bf_end_recovery=30G P_end_recovery=5 mito_gencode=mito_gencode.csv --cores 10
```
If you have a config file, run
```
snakemake  -s pipeline.smk --configfile config.yaml --cores 10
```

