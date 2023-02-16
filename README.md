# mtGasp (Mitochondrial Genome Assembly and Standardization Pipeline)
<img src="mtGasp_logo.png" width="200" height="200">


## Mitochondrial genome assembly and standardization pipeline for paired-end short DNA reads


# Credits
Concept: Lauren Coombe and Cecilia Yang

Implementation: Cecilia Yang

# Setup
1. Clone the repository

```
git clone https://github.com/bcgsc/mtGasp.git
```

2. Add the mtGasp directory to your PATH, use the following command to check if it is added correctly




3. Install dependencies


* Conda v4.13.0+ 
* Python v3.9+
* Snakemake 
* Pandas 
* Numpy 
* BLAST 
* Biopython 
* Seqtk 
* ABySS v2.2.0+
* ntJoin
* BWA 
* Samtools 
* Pilon 
* Mitos
* ntCard

---

#### Special Installation Instructions for Mitos


As Mitos uses an older Python version, please install it in a new conda environment called "mitos" using the instructions below:

```
conda create -n mitos
conda activate mitos
conda install python=2.7
conda install 'r-base>=4' r-ggplot2 r-reshape2 openjdk
conda install -c conda-forge biopython
conda install -c bioconda blast=2.9
conda install -c bioconda hmmer=3.2 infernal=1.1 'viennarna<2'
conda install -c bioconda mitos=2.0.8
```





# Running mtGasp

### Required Parameters 

`-o` or `--out_dir=DIR`: output folder name ***(full path or relative path)*** [Required]

`-r1` or `--read1=FILE`: compressed fastq.gz file containing the forward reads from paired-end sequencing ***(MUST be full path)*** [Required]

`-r2` or `--read2=FILE`: compressed fastq.gz file containing the reverse reads from paired-end sequencing ***(MUST be full path)*** [Required]



`-r` or `--ref_path=FILE`: path to the fasta file containing reference sequences that will be used to build blast database ***(MUST be full path)*** [Required]



***How to select sequences for the customized database?***

When reference sequences are available, you can simply store them in the fasta file you intend to use for blast database construction. 

If your target species do not have any reference sequences, you can include scaffolds or reference sequences of closely related species (i.e. the target species and references belong to the same kingdom or phylum) in your fasta file for database construction. 

(**Please note**: having more fasta sequences in the database will result in increased runtime and memory usage. The best practice is to have a maximum of one or two sequences belonging to the same species group (i.e. try to avoid duplicates of the same species))
---


`-m` or `--mt_gen=N`: mitochondrial genetic code (e.g., 2, 5) for your target species [Required]

***mito-genetic code can be searched on*** https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

---
###  Optional Parameters for Advanced Users 

`-k` or `--kmer=N`: k-mer size used in the construction of de bruijn graph for ABySS [96]

`-c` or `--kc=N`: k-mer minimum coverage multiplicity cutoff for ABySS [3]

***More information on optimizing k and kc for ABySS***: https://github.com/bcgsc/abyss#optimizing-the-parameters-k-and-kc

`-a` or `--abyss_fpr=N`: False positive rate for the bloom filter used by abyss during the assembly step [0.005]

`-s` or `--sealer_fpr=N`: False positive rate for the bloom filter used by sealer during the gap filling step [0.01]

`-p` or `--gap_filling_p=N`: Merge at most N alternate paths during sealer gap filling step; use 'nolimit' for no limit [5]

`-b` or `--sealer_k=STRING`: k-mer size(s) used in sealer gap filling ['60,80,100,120']

`-e` or `--end_recov_sealer_fpr=N`: False positive rate for the bloom filter used by sealer during flanking end recovery [0.01]

`-v` or `--end_recov_sealer_k=STRING`: k-mer size used in sealer flanking end recovery ['60,80,100,120']

`-i` or `--end_recov_p=N`: Merge at most N alternate paths during sealer flanking end recovery; use 'nolimit' for no limit [5]

`-t` or `--threads=N`: number of threads to use [3]


 
---

### Other Parameters

`-h` or `--help`: show help message and exit

`-n` or `--dry_run`: dry run mtGasp

`-u` or `--unlock`: Remove a lock implemented by snakemake on the working directory
    






---
### Examples

Dry run mtGasp with default parameters

    
```
mtgasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -r /path/to/mito_db/refs.fa -n
```

Run mtGasp with default parameters using 4 threads

```
mtgasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -r /path/to/mito_db/refs.fa -t 4
```

Run mtGasp with custom k-mer size and k-mer coverage cutoff

```
mtgasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -r /path/to/mito_db/refs.fa -k 80 -c 2
```

Snakemake uses a lock file to prevent other instances of Snakemake from running the same command simultaneously, if your working directory is locked by snakemake, use `-u or --unlock` to unlock the working directory

```
mtgasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -r /path/to/mito_db/refs.fa -u 
```
---
### Output Files

Output files can be found in a subfolder named `final_output` under the user-specified output directory `<out_dir>`. 

The `final_output` contains subfolders with the following format: `<out_dir>_k<k>_kc<kc>`, where k and kc are the k-mer size and k-mer coverage cutoff used for the assembly and out_dir is the user-specified output directory.

In each `<out_dir>_k<k>_kc<kc>` folder, there will be a fasta file named `<out_dir>_k<k>_kc<kc>.final-mtgasp-assembly.fa` storing the standardized mitochondrial sequence(s) and several subfolders and files containing the intermediate files generated by mitos during the mitogenome annotation step. 

Please note: `<out_dir>/final_output/<out_dir>_k<k>_kc<kc>/<out_dir>_k<k>_kc<kc>.final-mtgasp-assembly.fa` is the final output file containing the standardized mitochondrial sequence(s) generated by mtGasp. 

The mitos annotation results can be found in `<out_dir>/final_output/<out_dir>_k<k>_kc<kc>` with a prefix `result.` (e.g. result.fas, result.faa, result.gff)



---
### Clean up mtGasp intermediate files

```
cleanup.py <out_dir>
```
Here, this clean-up script will recursively clean up all intermediate files located in the user-specified mtGasp output directory `<out_dir>`. As a result, only the annotation results and the fasta file storing standardized mitochondrial sequence(s) are kept under `./<out_dir>/standardized_output/`. 



### Summarize mtGasp results

```
summarize.py -i <Input text file> -t <Output tsv file> -l <Output text file>
```
Here, this script will summarize the mtGasp results for all assembly output folders listed in the input text file `<Input text file>`. The output tsv file `<Output tsv file>` will contain the following columns:

`Assembly`: the name of the assembly output folder along with the k and kc values used for the assembly

`Number of sequences`: the number of mitochondrial sequences generated by mtGasp

`Scaffold lengths`: the length(s) of the mitochondrial sequence(s) generated by mtGasp

`Standardization status`: whether the mitochondrial sequence(s) generated by mtGasp is standardized or not

`Number of gaps (pre-GapFilling)`: the number of gaps in the mitochondrial sequence(s) generated by mtGasp before gap filling

`Number of gaps (post-GapFilling)`: the number of gaps in the mitochondrial sequence(s) generated by mtGasp after gap filling


Finally, `<Output text file>` will contain the full path(s) to the assembly output fasta file(s)

---
### How to fine-tune mtGasp parameters




