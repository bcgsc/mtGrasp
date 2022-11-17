import os.path
import shlex
import subprocess
import math

# check if library starts with a slash (full path) or not (relative path)
# if not, add current directory to the path
library = config["out_dir"]
if library.startswith("/"):
    current_dir = ""
else:
    current_dir = os.getcwd() + "/"



# Start of the pipeline
rule all:
     input:
        expand(current_dir + "{library}/standardized_output/{library}_k{k}_kc{kc}/post_standardization.fasta", library = config["out_dir"], k = config["k"], kc = config["kc"])




# Bloom filter size calculation for abyss-pe
def bloom_filter_calc(file):
    with open(file, 'r') as f:
        content = f.readlines()
        # second line: F0
        distinct = int(content[1].strip().split()[1])
        # third line: number of k-mers that only occurred once
        once = int(content[2].strip().split()[1])
        count = distinct - once # the number of k-mers that occurred more than once
    bits = math.ceil(-count/math.log(0.998)) # FPR = 0.2%
    bytes = bits/8
    # convert bytes to B, KB, MB, GB, TB
    power = 10**3
    n = 0
    power_labels = {0 : '', 1: 'K', 2: 'M', 3: 'G', 4: 'T'}
    while bytes > power:
        bytes /= power
        n += 1
    return str(math.ceil(bytes))+power_labels[n]

def bf_abyss(r1, r2, output_dir, k):
    new_dir = f'{output_dir}/bloom_filter_calc/abyss'
    # check if a file exists
    if not os.path.exists(new_dir + '/' + f'freq_k{k}.hist'):
        cmd = f'mkdir -p {new_dir} && cd {new_dir} && ntcard -k{k} -p freq {r1} {r2}'
        os.system(cmd)
        return bloom_filter_calc(new_dir + '/' + f'freq_k{k}.hist')
    else:
        return bloom_filter_calc(new_dir + '/' + f'freq_k{k}.hist')


# Bloom filter size calculation for abyss-sealer
def bf_sealer(r1, r2, output_dir):
    new_dir = f'{output_dir}/bloom_filter_calc/sealer'
    # check if a directory exists
    if not os.path.exists(new_dir):
        cmd = f'mkdir -p {new_dir} && cd {new_dir} && ntcard -k60,80,100,120 -p freq {r1} {r2}' # k-mer sweeps used in sealer: 60, 80, 100, 120
        os.system(cmd)
        kmer = sealer_max_kmer_count(new_dir)
        return bloom_filter_calc(f'{new_dir}/freq_{kmer}.hist')
    else:    
        kmer = sealer_max_kmer_count(new_dir)
        return bloom_filter_calc(f'{new_dir}/freq_{kmer}.hist')

def sealer_max_kmer_count(new_dir): 
    count_list = []
    kmer = []
    for file in os.listdir(new_dir):
          if file.endswith(".hist"):
            with open(f'{new_dir}/{file}', 'r') as f:
                content = f.readlines()
                k = file.split('.')[0].split('_')[1]
                kmer.append(k)
                # second line: F0
                distinct = int(content[1].strip().split()[1])
                # third line: number of k-mers that only occurred once
                once = int(content[2].strip().split()[1])
                # split by tab
                count = distinct - once
                count_list.append(count)
        # get the max count and the kmer 
    return kmer[count_list.index(max(count_list))]


out_dir = current_dir + config["out_dir"]

cmd = f'mkdir -p {out_dir}'
cmd_shlex = shlex.split(cmd)
subprocess.call(cmd_shlex)     
    
#De novo assembly 
rule de_novo_assembly:
     input:
        r1=expand("{read1}", read1=config['r1']),
        r2=expand("{read2}", read2=config['r2'])
     output:
        current_dir + "{library}/abyss/k{k}_kc{kc}-scaffolds.fa"
     benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.de_novo_assembly.benchmark.txt"
     log:
        current_dir + "{library}/abyss/k{k}_kc{kc}.log"
     run:
        bf = bf_abyss(input.r1, input.r2, wildcards.library, wildcards.k)
        shell("abyss-pe --directory {wildcards.library}/abyss v=-v kc={wildcards.kc}  k={wildcards.k}  B={bf}  name=k{wildcards.k}_kc{wildcards.kc} in='{input.r1} {input.r2}' &> {log}")
       

       

#Filtering for mitochondrial sequences
rule select_length:
     input:
        rules.de_novo_assembly.output
     output:
        current_dir + "{library}/blast/k{k}_kc{kc}-scaffolds.1000-20000bp.fa"
     benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.select_length.benchmark.txt"
     shell:
        """awk '!/^>/ {{ next }} {{ getline seq }} length(seq) > 1000 {{ print $0 "\\n" seq }}' {input} | awk '!/^>/ {{ next }} {{ getline seq }} length(seq) < 20000 {{ print $0 "\\n" seq }}' > {output}"""


rule blast:
     input:
         rules.select_length.output
     output:
         current_dir + "{library}/blast/k{k}_kc{kc}-scaffolds.1000-20000bp.nt.blast.tsv"
     benchmark:
         current_dir + "{library}/benchmark/k{k}_kc{kc}.blast.benchmark.txt"
     params:
         ref_fasta=config["ref_path"]
     run:
         # get the directory of the script
         string = subprocess.check_output(['which', 'mtgasp.smk'])
         string = string.decode('utf-8')
         # remove new line character
         string = string.strip()
         # split string by '/'
         script_dir = '/'.join(string.split('/')[0:-1])
         db_name = os.path.splitext(os.path.basename(params.ref_fasta))[0]

         # check if a blast database exists
         if not os.path.exists(f'{script_dir}/blast_db/{db_name}'):
            shell("mkdir -p {script_dir}/blast_db/{db_name} && cd {script_dir}/blast_db/{db_name} && makeblastdb -in {params.ref_fasta} -dbtype nucl -out {db_name} -parse_seqids")
            shell("export BLASTDB={script_dir}/blast_db/{db_name} && blast_best-hit.py {input} {db_name} > {output}")
         else:
            shell("export BLASTDB={script_dir}/blast_db/{db_name} && blast_best-hit.py {input} {db_name} > {output}")

rule create_lists:
     input:
         rules.blast.output
     output:
         ref_list=current_dir + "{library}/blast/k{k}_kc{kc}-query.txt",
         query_list=current_dir + "{library}/blast/k{k}_kc{kc}-ref.txt"
     benchmark:
         current_dir + "{library}/benchmark/k{k}_kc{kc}.create_lists.benchmark.txt"
     shell:
          "extract_tsv_value.py -in {input} -out {output.ref_list} -c ref ; "
          "extract_tsv_value.py -in {input} -out {output.query_list} -c query"
         

rule extract_seq:
     input:
         ref=rules.create_lists.output.ref_list,
         query=rules.create_lists.output.query_list,
         assemblies=rules.select_length.output
     output:
         query_out=current_dir + "{library}/mito_filtering_output/k{k}_kc{kc}-scaffolds.1000-20000bp.blast-mt_db.fa",
         ref_out=current_dir + "{library}/mito_filtering_output/k{k}_kc{kc}-ref.fa"
     benchmark:
         current_dir + "{library}/benchmark/k{k}_kc{kc}.extract_seq.benchmark.txt"
     params:
         ref_fasta=config["ref_path"]
     shell:
         "seqtk subseq {input.assemblies} {input.query} > {output.query_out} ; "
         "seqtk subseq {params.ref_fasta} {input.ref} > {output.ref_out}"


#Prepare for polishing: ntJoin+Sealer or Sealer
rule pre_polishing:
      input:
          target=rules.extract_seq.output.query_out,
          ref=rules.extract_seq.output.ref_out
      output:
          current_dir + "{library}/sealer/k{k}_kc{kc}.postsealer_scaffold.fa"
      params:
        workdir= current_dir + "{library}/mito_filtering_output",
        out=current_dir + "{library}/sealer/k{k}_kc{kc}.postsealer",
        ntjoin_out=current_dir + "{library}/mito_filtering_output/k{k}_kc{kc}-scaffolds.1000-20000bp.blast-mt_db.fa.k32.w500.n1.all.scaffolds.fa",
        r1=config['r1'],
        r2=config['r2']
      benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.pre_polishing.benchmark.txt"
      log:
        sealer=current_dir + "{library}/sealer/k{k}_kc{kc}_sealer.log",
        ntjoin=current_dir + "{library}/mito_filtering_output/k{k}_kc{kc}_ntjoin.log"
      run: 
          import os
          target = os.path.basename(input.target)
          ref = os.path.basename(input.ref)
          log_ntjoin = os.path.basename(log.ntjoin)
          log_sealer = log.sealer
          bf = bf_sealer(params.r1, params.r2, wildcards.library)
          count = sum(1 for line in open(input[0]))
          if count <= 2:
               shell("""abyss-sealer -b{bf} -k 60 -k 80 -k 100 -k 120 -P 5 -o {params.out} -S {input.target} {params.r1} {params.r2} &> {log_sealer}""")
               print("no ntJoin needed")
          else:
               shell("""bash run_ntjoin.sh {params.workdir} {target} {ref} {log_ntjoin} && abyss-sealer -b{bf} -k 60 -k 80 -k 100 -k 120 -P 5 -o {params.out} -S {params.ntjoin_out}  {params.r1} {params.r2} &> {log_sealer}""")
               print("ntJoin needed")



#Polishing
rule bwa_alignment:
      input:
        rules.pre_polishing.output
      output:
        current_dir + "{library}/pilon/k{k}_kc{kc}.sorted.bam"
      params:
        r1=config['r1'],
        r2=config['r2']
      benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.bwa_alignment.benchmark.txt"
      shell:
        "bwa index {input} && "
        "bwa mem {input} {params.r1} {params.r2} | samtools view -b -F 4 | samtools sort -o {output} && "
        "samtools index {output}"

rule polishing:
      input:
        in1=rules.pre_polishing.output,
        in2=rules.bwa_alignment.output
      output:
        current_dir + "{library}/assembly_output/{library}_k{k}_kc{kc}.postsealer.postpilon.fasta"
      params:
        outdir=current_dir + "{library}/pilon/{library}_k{k}_kc{kc}.postsealer.postpilon",
        fasta=current_dir + "{library}/pilon/{library}_k{k}_kc{kc}.postsealer.postpilon.fasta"
      benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.polishing.benchmark.txt"
      shell:
        "pilon --genome {input.in1} --frags {input.in2} --output {params.outdir} --changes --fix all --verbose && "
        "mv {params.fasta} {output}"

# Standardization of strand orientation and start-site
rule end_recovery:
       input:
          rules.polishing.output
       output:
          current_dir + "{library}/end_recovery/{library}_k{k}_kc{kc}/flank_added_assembly.fa"
       params:
          r1=config['r1'],
          r2=config['r2'],
          outdir=current_dir + "{library}/end_recovery/{library}_k{k}_kc{kc}"
       benchmark:
          current_dir + "{library}/benchmark/k{k}_kc{kc}.end_recovery.benchmark.txt"
       run:
          bf = bf_sealer(params.r1, params.r2, wildcards.library)
          shell("end_recover.py {input} {bf} {params.r1} {params.r2} {params.outdir}")

rule standardization:
        input:
            rules.end_recovery.output
        output:
            current_dir + "{library}/standardized_output/{library}_k{k}_kc{kc}/post_standardization.fasta"
        benchmark:
            current_dir + "{library}/benchmark/k{k}_kc{kc}.standardization.benchmark.txt"
        params:
            mito_gencode=config["mt_code"],
            outdir=current_dir + "{library}/standardized_output/{library}_k{k}_kc{kc}"
        shell:
            "mitos_annotation.py {input} {params.mito_gencode} {params.outdir} {wildcards.library}_k{wildcards.k}_kc{wildcards.kc}"
