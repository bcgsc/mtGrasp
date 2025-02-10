[![link](https://img.shields.io/badge/mtGrasp-manuscript-brightgreen)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14506)


# mtGrasp : Mitochondrial Genome Reference-grade Assembly and Standardization Pipeline

Mitochondrial genome assembly and standardization pipeline for paired-end short DNA reads


## Contents

1. [Credit](#credit)
2. [Description](#description)
3. [Installation](#install)
4. [Dependencies](#dependencies)
5. [Demo](#demo)
6. [Usage](#usage)
7. [Examples](#examples)
8. [Output](#output)
9. [MtG standardization](#standardize)
10. [Summarize mtGrasp output](#summary)
11. [Citing](#citing)
12. [License](#license)

## Credit  <a name=credit></a>
Concept: Lauren Coombe and Cecilia Yang
Implementation: Cecilia Yang

## Description <a name=description></a>

The Mitochondrial Genome Reference-grade Assembly and Standardization Pipeline (mtGrasp) is a streamlined, high-throughput mitogenome assembly utility

## Installation <a name=install></a>

### Installation using conda
You can easily install mtGrasp using conda/mamba:
```
conda install -c conda-forge -c bioconda mtgrasp
```

### Installation from the source code
1. Download mtGrasp package

  You can download the latest version of the tarball file [here](https://github.com/bcgsc/mtGrasp/releases/download/v1.1.8/mtGrasp-v1.1.8.tar.gz).


2. Add the mtGrasp directory to your PATH, use the following command to check if it is added correctly
```
echo $PATH
```

3. Install dependencies <a name=dependencies></a>

* Python v3.9+
* Snakemake 
* BLAST v2.9+
* Biopython 
* Seqtk 
* ABySS v2.2.0+
* ntJoin v1.1.3+
* BWA 
* Samtools 
* Pilon v1.24+
* MITOS v2.1.7+
* ntCard

### Installing dependencies using conda:
Recommended (Faster):
```
conda create -n mtgrasp python=3.10 mamba
conda activate mtgrasp
mamba install -c conda-forge -c bioconda snakemake 'blast>=2.9.0' biopython seqtk abyss ntjoin bwa samtools pilon ntcard 'mitos>=2.1.7'
```

Alternative (Slower):
```
conda create -n mtgrasp python=3.10
conda activate mtgrasp
conda install -c conda-forge -c bioconda snakemake 'blast>=2.9.0' biopython seqtk abyss ntjoin bwa samtools pilon ntcard 'mitos>=2.1.7'
```

### Using docker/singularity
Users can also use the docker images located at [Quay.io](https://quay.io/repository/biocontainers/mtgrasp?tab=tags) to run mtGrasp using docker or singularity.

# Test run <a name=demo></a>
### Test-run mtGrasp to ensure all required dependencies are installed properly
The test will take ~5-10min to complete.

```
mtgrasp.py -test
```

If `runmitos.py` is not available on your PATH:
```
mtgrasp.py -test -mp /path/to/mitos_install_dir
```
Note: `/path/to/mitos_install_dir` is the location where the main MITOS script `runmitos.py` is stored
# Running mtGrasp <a name=usage></a>

### Required Parameters 

`-o` or `--out_dir=DIR`: output folder name ***(full path or relative path)***

`-r1` or `--read1=FILE`: compressed fastq.gz file containing the forward reads from paired-end sequencing ***(MUST be full path)***

`-r2` or `--read2=FILE`: compressed fastq.gz file containing the reverse reads from paired-end sequencing ***(MUST be full path)***

`-r` or `--ref_path=FILE`: path to the fasta file containing reference sequences that will be used to build blast database ***(MUST be full path)***


***How to select sequences for the customized database?***

If the reference sequences of the species whose mitogenomes you plan to assemble are available, you can conveniently add them to the fasta file that will later be used by mtGrasp to construct the mitochondrial sequence filtering database.

In case reference sequences are not accessible for your target species, you can consider incorporating scaffolds or reference sequences from closely related species. The ideal scenario would be if both the target and reference species are from the same family or genus, as this would improve the quality of the mitochondrial database by improving the probability of matching the assembled mitochondrial sequence with your mitogenome database. 

However, if such sequences are unavailable, you can move up the taxonomic hierarchy until you find a suitable sequence. Additionally, it is worth trying different sequences as altering the sequences in the database does not impact the final quality of the assembly process, it simply increases the likelihood of successful mitochondrial sequence searching in cases where they are successfully assembled.

(**Please note**: Having more fasta sequences in the database will result in increased runtime and memory usage. (i.e. try to avoid duplicates of the same species))
---


`-m` or `--mt_gen=N`: mitochondrial translation table code (e.g., 2, 5, 13) for your target species

***Mitochondrial translation table code can be searched on*** https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi

---
###  Optional Parameters for Advanced Users 

`-k` or `--kmer=N`: k-mer size used in the construction of de Bruijn graph for ABySS [91]

`-c` or `--kc=N`: k-mer minimum coverage multiplicity cutoff for ABySS [3]

***More information on optimizing -k and -c for ABySS***: https://github.com/bcgsc/abyss#optimizing-the-parameters-k-and-kc

`-sub` or `--subsample=N`: Subsample N read pairs from two paired FASTQ files  [2000000] 

`-nsub` or `--nosubsample`: Run mtGrasp using the entire read dataset without subsampling [False]

`-an` or `--annotate`: Run gene annotation on the final assembly output [False]

`-d` or `--delete`: Delete intermediate subdirectories/files once mtGrasp finishes [False]

`-mp` or `--mitos_path`: Complete path to `runmitos.py` (e.g., `/home/user/path/to/mitos/bin`), this is required if `runmitos.py` is not found on your `PATH` [None]

`-test` or `--test_run`:Test run mtGrasp to ensure all required dependencies are installed [False]

`-p` or `--gap_filling_p=N`: Merge at most N alternate paths during sealer gap filling step; use 'nolimit' for no limit [5]

`-b` or `--sealer_k=STRING`: k-mer size(s) used in Sealer gap filling ['60,80,100,120']

`-sk` or `--end_recov_sealer_k=STRING`: k-mer size used in sealer flanking end recovery ['60,80,100,120']

`-i` or `--end_recov_p=N`: Merge at most N alternate paths during sealer flanking end recovery; use 'nolimit' for no limit [5]

`-t` or `--threads=N`: number of threads to use [8]

`-ma` or `--mismatch_allowed=N`: Maximum number of mismatches allowed while determining the overlapping region between the two ends of the mitochondrial assembly [1]

 
---

### Other Parameters

`-h` or `--help`: show help message and exit

`-n` or `--dry_run`: dry run mtGrasp

`-u` or `--unlock`: Remove a lock implemented by snakemake on the working directory

`-v` or `--version`: Print out the version of mtGrasp and exit






---
### Examples <a name=examples></a>

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
### Where to Look For Output Files <a name=output></a>

Output files can be found in a subfolder named `final_output` under the user-specified output directory `<out_dir>`. 

The `final_output` contains subfolders with the following format: `<out_dir>_k<k>_kc<kc>`, where k and kc are the k-mer size and k-mer coverage cutoff used for the assembly and out_dir is the user-specified output directory.

In each `<out_dir>_k<k>_kc<kc>` folder, there will be a fasta file named `<out_dir>_k<k>_kc<kc>.final-mtgrasp-assembly.fa` storing the standardized mitochondrial sequence(s) and several subfolders and files containing the intermediate files generated by MITOS during the mitogenome annotation step. 

Please note: `<out_dir>/final_output/<out_dir>_k<k>_kc<kc>/<out_dir>_k<k>_kc<kc>.final-mtgrasp-assembly.fa` is the final output file containing the standardized mitochondrial sequence(s) generated by mtGrasp. 

If the `-an` or `--annotate` argument is provided, mtGrasp will run gene annotation for the final assembly output and the annotation results can be found in `<out_dir>/final_output/annotation_output`.

If you are not interested in the standardized mitogenome sequence(s) or MITOS annotation failed during the standardization, you can find the non-standardized mitogenome sequence output(s) generated by mtGrasp in `<out_dir>/assembly_output/<out_dir>_k<k>_kc<kc>.postsealer.postpilon.fasta`. 



---
### Standardize any mitochondrial sequence(s) using mtGrasp <a name=standardize></a>
If you have any mitochondrial sequence(s) and would like to standardize it/them using mtGrasp, you can do so by using mtGrasp's `mtgrasp_standardize.py` script. 

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

The annotated proteins can be found in `<out_dir>/annotation_output/result.faa`.

The annotated transcripts can be found in `<out_dir>/annotation_output/result.fas`.

The order of the annotated genes can be found in `<out_dir>/annotation_output/result.geneorder`.

Because mtGrasp annotation uses a third-party tool called [MITOS](https://www.sciencedirect.com/science/article/abs/pii/S1055790312003326), any inquires regarding the annotation results should be directed to the [MITOS developers](https://gitlab.com/Bernt/MITOS).

***Please note***: 
- Currently, mtGrasp only supports animal mitochondrial sequences.
- When you are using `mtgrasp_standardize.py` as a stand-alone tool, the input fasta MUST contain only one mitochondrial sequence, currently, mtGrasp does not support standardizing multiple mitochondrial sequences in one fasta file. 
- The headers of the fasta sequences contain information about whether the sequence is start-site standardized and strand standardized or not. Start-site standardized means the sequence starts with tRNA-Phe and strand standardized means the final sequence is on the positive strand. 
- The fasta headers also contain 'Linear', however, this does not mean the sequence is indeed linear, it simply means the sequence did not go through the circularization step during the mtGrasp pipeline. 



---

### Summarize mtGrasp results <a name=summary></a>

```
mtgrasp_summarize.py -i <Input text file> -p <Prefix of the summary files>
```

This script will summarize the mtGrasp results for all assembly output folders listed in the input text file `<Input text file>`. The output tsv file `{prefix}_mtgrasp_{mtgrasp_version}_assembly_summary.tsv'` will contain the following columns:

`Assembly`: the name of the assembly output folder along with the k and kc values (number of read pairs if `-sub` is enabled) used for the assembly

`Ns per 1000bp`: the number of Ns per 1000bp in the mitochondrial sequence(s) generated by mtGrasp

`Number of Contigs`: the number of mitochondrial sequences generated by mtGrasp

`Total Number of Base Pairs Per Assembly`: the total number of base pairs in the mitochondrial sequence(s) generated by mtGrasp

`Length of the Longest Contig (bp)`: the length of the longest mitochondrial sequence generated by mtGrasp

`Circular or Linear`: whether the mitochondrial sequence(s) generated by mtGrasp is circular or linear

`Standardization Status`: whether the mitochondrial sequence(s) generated by mtGrasp is start-site standardized and strand standardized


Finally, `{prefix}_mtgrasp_{mtgrasp_version}_path_to_output.txt` will contain the complete path(s) to the assembly output fasta file(s)


## Citing <a name=citing></a>

Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/mtGrasp.svg)](https://github.com/bcgsc/mtGrasp/stargazers) and for using and promoting this free software! We hope that mtGrasp is useful to you and your research.

If you use mtGrasp, please cite:

[mtGrasp: Streamlined reference-grade mitochondrial genome assembly and standardization to enhance metazoan mitogenome resources](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14506)
<pre>
Lopez, M. L. D., Yang, C. L., Coombe, L., Warren, R. L., Allison, M. J., Imbery, J. J., Birol, I., & Helbing, C. C. (2025). mtGrasp: Streamlined reference-grade mitochondrial genome assembly and standardization to enhance metazoan mitogenome resources. Methods in Ecology and Evolution, 00, 1–10. https://doi.org/10.1111/2041-210X.14506
</pre>


## License <a name=license></a>

mtGrasp Copyright (c) 2024-present British Columbia Cancer Agency Branch.  All rights reserved.

mtGrasp is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

For commercial licensing options, please contact
Patrick Rebstein <prebstein@bccancer.bc.ca>
