# mtGrasp (Mitochondrial Reference-Grade Genome Assembly and Standardization Pipeline)



## Mitochondrial genome assembly and standardization pipeline for paired-end short DNA reads


# Credits
Concept: Lauren Coombe and Cecilia Yang

Implementation: Cecilia Yang

# Setup
1. Clone the repository

```
git clone https://github.com/bcgsc/mtGrasp.git
```

2. Add the mtGrasp directory to your PATH, use the following command to check if it is added correctly




3. Install dependencies


* Conda v4.13.0+ 
* Python v3.9+
* Snakemake 
* Pandas 
* Numpy 
* BLAST v2.10+
* Biopython 
* Seqtk 
* ABySS v2.2.0+
* ntJoin
* BWA 
* Samtools 
* Pilon v1.24+
* MITOS
* ntCard

---

#### Special Installation Instructions for MITOS


As MITOS uses an older Python version, please install it in a new conda environment called "mitos" using the instructions below:



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





# Running mtGrasp

### Required Parameters 

`-o` or `--out_dir=DIR`: output folder name ***(full path or relative path)*** [Required]

`-r1` or `--read1=FILE`: compressed fastq.gz file containing the forward reads from paired-end sequencing ***(MUST be full path)*** [Required]

`-r2` or `--read2=FILE`: compressed fastq.gz file containing the reverse reads from paired-end sequencing ***(MUST be full path)*** [Required]



`-r` or `--ref_path=FILE`: path to the fasta file containing reference sequences that will be used to build blast database ***(MUST be full path)*** [Required]



***How to select sequences for the customized database?***

If the reference sequences of the species whose mitogenomes you plan to assemble are available, you can conveniently add them to the fasta file that will later be used by mtGrasp to construct the mitochondrial sequence filtering database.

In case reference sequences are not accessible for your target species, you can consider incorporating scaffolds or reference sequences from closely related species. The ideal scenario would be if both the target and reference species are from the same family or genus, as this would improve the quality of the mitochondrial database by improving the probability of matching the assembled mitochondrial sequence with your mitogenome database. 

However, if such sequences are unavailable, you can move up the taxonomic hierarchy until you find a suitable sequence. Additionally, it is worth trying different sequences as altering the sequences in the database does not impact the final quality of the assembly process, it simply increases the likelihood of successful mitochondrial sequence searching in cases where they are successfully assembled.

(**Please note**: Having more fasta sequences in the database will result in increased runtime and memory usage. The best practice is to have a maximum of one or two sequences belonging to the same species group (i.e. try to avoid duplicates of the same species))
---


`-m` or `--mt_gen=N`: mitochondrial translation table code (e.g., 2, 5, 13) for your target species [Required]

***Mitochondrial translation table code can be searched on*** https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

---
###  Optional Parameters for Advanced Users 

`-k` or `--kmer=N`: k-mer size used in the construction of de bruijn graph for ABySS [91] (Please note: k-mer size must be less than 128)

`-c` or `--kc=N`: k-mer minimum coverage multiplicity cutoff for ABySS [3]

***More information on optimizing -k and -c for ABySS***: https://github.com/bcgsc/abyss#optimizing-the-parameters-k-and-kc

`-sub` or `--subsample=N`: Subsample N read pairs from two paired FASTQ files  [2000000] 

`-nsub` or `--nosubsample`: Run mtGrasp using the entire read dataset without subsampling [False]

`-an` or `--annotate`: Run gene annotation on the final assembly output [False]

`-d` or `--delete`: Delete intermediate subdirectories/files once mtGrasp reaches completion [False]


`-a` or `--abyss_fpr=N`: False positive rate for the bloom filter used by abyss during the assembly step [0.005]

`-s` or `--sealer_fpr=N`: False positive rate for the bloom filter used by sealer during the gap filling step [0.01]

`-p` or `--gap_filling_p=N`: Merge at most N alternate paths during sealer gap filling step; use 'nolimit' for no limit [5]

`-b` or `--sealer_k=STRING`: k-mer size(s) used in sealer gap filling ['60,80,100,120']

`-e` or `--end_recov_sealer_fpr=N`: False positive rate for the bloom filter used by sealer during flanking end recovery [0.01]

`-v` or `--end_recov_sealer_k=STRING`: k-mer size used in sealer flanking end recovery ['60,80,100,120']

`-i` or `--end_recov_p=N`: Merge at most N alternate paths during sealer flanking end recovery; use 'nolimit' for no limit [5]

`-t` or `--threads=N`: number of threads to use [8]

`-ma` or `--mismatch_allowed=N`: Maximum number of mismatches allowed while determining the overlapping region between the two ends of the mitochondrial assembly [1]


 
---

### Other Parameters

`-h` or `--help`: show help message and exit

`-n` or `--dry_run`: dry run mtGrasp

`-u` or `--unlock`: Remove a lock implemented by snakemake on the working directory
    






---
### Examples

Dry run mtGrasp with default parameters

    
```
mtgrasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -r /path/to/mito_db/refs.fa -n
```

Run mtGrasp with default parameters using 4 threads

```
mtgrasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -r /path/to/mito_db/refs.fa -t 4
```

Run mtGrasp with custom k-mer size and k-mer coverage cutoff

```
mtgrasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -r /path/to/mito_db/refs.fa -k 80 -c 2
```

Snakemake uses a lock file to prevent other instances of Snakemake from running the same command simultaneously, if your working directory is locked by snakemake, use `-u or --unlock` to unlock the working directory

```
mtgrasp.py -r1 /path/to/read1.fq.gz -r2 /path/to/read2.fq.gz -o test_out -m 2 -r /path/to/mito_db/refs.fa -u 
```
---
### Where to Look For Output Files

Output files can be found in a subfolder named `final_output` under the user-specified output directory `<out_dir>`. 

The `final_output` contains subfolders with the following format: `<out_dir>_k<k>_kc<kc>`, where k and kc are the k-mer size and k-mer coverage cutoff used for the assembly and out_dir is the user-specified output directory.

In each `<out_dir>_k<k>_kc<kc>` folder, there will be a fasta file named `<out_dir>_k<k>_kc<kc>.final-mtgrasp-assembly.fa` storing the standardized mitochondrial sequence(s) and several subfolders and files containing the intermediate files generated by MITOS during the mitogenome annotation step. 

Please note: `<out_dir>/final_output/<out_dir>_k<k>_kc<kc>/<out_dir>_k<k>_kc<kc>.final-mtgrasp-assembly.fa` is the final output file containing the standardized mitochondrial sequence(s) generated by mtGrasp. 

If the `-an` or `--annotate` argument is provided, mtGrasp will run gene annotation for the final assembly output and the annotation results can be found in `<out_dir>/final_output/annotation_output`.

If you are not interested in the standardized mitogenome sequence(s) or MITOS annotation failed during the standardization, you can find the non-standardized mitogenome sequence output(s) generated by mtGrasp in `<out_dir>/assembly_output/<out_dir>_k<k>_kc<kc>.postsealer.postpilon.fasta`. 



### Summarize mtGrasp results

```
mtgrasp_summarize.py -i <Input text file> -t <Output tsv file> -l <Output text file>
```
Here, this script will summarize the mtGrasp results for all assembly output folders listed in the input text file `<Input text file>`. The output tsv file `<Output tsv file>` will contain the following columns:

`Assembly`: the name of the assembly output folder along with the k and kc values used for the assembly

`Number of sequences`: the number of mitochondrial sequences generated by mtGrasp

`Scaffold lengths`: the length(s) of the mitochondrial sequence(s) generated by mtGrasp

`Standardization status`: whether the mitochondrial sequence(s) generated by mtGrasp is standardized or not

`Number of gaps (pre-GapFilling)`: the number of gaps in the mitochondrial sequence(s) generated by mtGrasp before gap filling

`Number of gaps (post-GapFilling)`: the number of gaps in the mitochondrial sequence(s) generated by mtGrasp after gap filling


Finally, `<Output text file>` will contain the full path(s) to the assembly output fasta file(s)

---
### Standardize your own mitochondrial sequence(s) using mtGrasp
If you have your own mitochondrial sequence(s) and would like to standardize it/them using mtGrasp, you can do so by using mtGrasp's `mtgrasp_standardize.py` script. 

For usage, please run `mtgrasp_standardize.py -h` for help:

`-i` or `--input=FILE`: input fasta file containing the mitochondrial sequence(s) to be standardized [Required]

`-o` or `--out_dir=DIR`: output folder name [Required]

`-c` or `--gencode=N`: mitochondrial translation table code (e.g., 2, 5, 13) for your target species [Required]

***Mitochondrial translation table code can be searched on*** https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

`-p` or `--prefix=STRING`: prefix name for the output files [mtgrasp_standardized]

`-a` or `--annotate`: Run gene annotation on the final assembly output [False]




Example:

```
mtgrasp_standardize.py -i /path/to/your/mito_seq.fa -o /path/to/output_dir -c 2 -p prefix_name -a
```

The final output file containing the standardized mitochondrial sequence(s) can be found in `<out_dir>/<prefix_name>.final-mtgrasp-assembly.fa`


If the `-a` or `--annotate` argument is provided, mtGrasp will run gene annotation for the final assembly output and the annotation results can be found in `<out_dir>/annotation_output`.

The amino acid sequences of the annotated genes can be found in `<out_dir>/annotation_output/result.faa`.

The nucleotide sequences of the annotated genes can be found in `<out_dir>/annotation_output/result.fas`.

The order of the annotated genes can be found in `<out_dir>/annotation_output/result.geneorder`.

Because mtGrasp annotation uses a third-party tool called [MITOS](https://www.sciencedirect.com/science/article/abs/pii/S1055790312003326), any inquires regarding the annotation results should be directed to the [MITOS developers](https://gitlab.com/Bernt/MITOS).

***Please note***: 
- Currently, mtGrasp only supports standardizing animal mitochondrial sequences.
- When you are using `mtgrasp_standardize.py` as a stand-alone tool, the input fasta MUST contain only one mitochondrial sequence, currently, mtGrasp does not support standardizing multiple mitochondrial sequences in one fasta file. 
- The headers of the fasta sequences contain information about whether the sequence is start-site standardized and strand standardized or not. Start-site standardized means the sequence starts with tRNA-Phe and strand standardized means the final sequence is on the positive strand. 
- The fasta headers also contain 'Linear', however, this does not mean the sequence is indeed linear, it simply means the sequence did not go through the circularization step during the mtGrasp pipeline. 





