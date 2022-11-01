# mtGasp (Mitochondrial Genome Assembly and Standardization Pipeline)



## Novel mitochondrial genome assembly and standardization pipeline for paired-end short DNA reads



# Setup
1. Clone the repository

```
git clone https://github.com/bcgsc/mtGasp.git
```
2. Add the mtGasp directory to your PATH, use the following command to check if it is added correctly

```
echo $PATH
```


3. Install dependencies



``` Dependencies ```

Anaconda 

Python 3.9+

Snakemake 

Pandas 

Numpy 

BLAST 

Biopython 

Seqtk 

AbySS v2.2.0+

ABySS-Sealer 

ntJoin

BWA 

Samtools 

Pilon 

Mitos 2.0.8


---

```Special Installation Instructions for MitoS```

***Please note: mitos 2.0.8 uses python 2.7, if downloaded to the same environment with the other dependencies, a conflict will occur.***

Please install mitos in a new environment called "mitos" using the instruction below:

```
conda install mitos -c bioconda -m -n mitos
```


# Running mtGasp

### Required Parameters 

`-o` or `--out_dir`: output folder name ***(full path or relative path)***

`-r1` or `--read1`: compressed fastq.gz file containing the forward reads from paired-end sequencing ***(MUST be full path)***

`-r2` or `--read2`: compressed fastq.gz file containing the reverse reads from paired-end sequencing ***(MUST be full path)***


`-p` or `--ref_path`: path to the fasta file containing reference sequences that will be used to build blast database ***(MUST be full path)***

(**Please note**: having more fasta sequences in the database will result in increased runtime and memory usage. The best practice is to have a maximum of one or two sequences that belong to the same species group)

`-m` or `--mt_gen`: mitochondrial genetic code (e.g., 2, 5) for your target species

***mito-genetic code can be searched on*** https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

---
###  Optional Parameters for Advanced Users

`-k` or `--kmer`: k-mer size used in the construction of de bruijn graph for ABySS

`-c` or `--kc`: k-mer minimum coverage multiplicity cutoff for ABySS
 
---

### Other Parameters

`-h` or `--help`: show help message and exit

`-t` or `--threads`: number of threads to use

`-n` or `--dry_run`: dry run mtGasp
    


***More information on optimizing k and kc for ABySS***: https://github.com/bcgsc/abyss#optimizing-the-parameters-k-and-kc

---
### Examples

Dry run mtGasp with default parameters

    
```
mtgasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -p /path/to/mito_db/refs.fa -n
```

Run mtGasp with default parameters using 4 threads

```
mtgasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -p /path/to/mito_db/refs.fa -t 4
```

Run mtGasp with custom k-mer size and k-mer coverage cutoff

```
mtgasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -p /path/to/mito_db/refs.fa -k 80 -c 2
```


