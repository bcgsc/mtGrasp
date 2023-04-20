import os.path
import shlex
import subprocess
import math
from Bio import SeqIO

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
        expand(current_dir + "{library}/final_output/{library}_k{k}_kc{kc}/{library}_k{k}_kc{kc}.final-mtgasp-assembly.fa", library = config["out_dir"], k = config["k"], kc = config["kc"])




# Bloom filter size calculation for abyss-pe
def bloom_filter_calc(file, FPR):
    with open(file, 'r') as f:
        content = f.readlines()
        # second line: F0
        distinct = int(content[1].strip().split()[1])
        # third line: number of k-mers that only occurred once
        once = int(content[2].strip().split()[1])
        count = distinct - once # the number of k-mers that occurred more than once
    bits = math.ceil(-count/math.log(1-FPR))
    bytes = bits/8
    # convert bytes to B, KB, MB, GB, TB
    power = 10**3
    n = 0
    power_labels = {0 : '', 1: 'K', 2: 'M', 3: 'G', 4: 'T'}
    while bytes > power:
        bytes /= power
        n += 1
    return str(math.ceil(bytes))+power_labels[n]

def bf_abyss(r1, r2, output_dir, k, FPR):
    new_dir = f'{output_dir}/bloom_filter_calc/abyss'
    # check if a file exists
    if not os.path.exists(new_dir + '/' + f'freq_k{k}.hist'):
        cmd = f'mkdir -p {new_dir} && cd {new_dir} && ntcard -k{k} -p freq {r1} {r2}'
        os.system(cmd)
        return bloom_filter_calc(new_dir + '/' + f'freq_k{k}.hist', FPR)
    else:
        return bloom_filter_calc(new_dir + '/' + f'freq_k{k}.hist', FPR)


# Bloom filter size calculation for abyss-sealer
def bf_sealer(r1, r2, output_dir, ntcard_threads, FPR,kmers):
    new_dir = f'{output_dir}/bloom_filter_calc/sealer'
    for k in kmers.split(','):
       hist_file = f'{new_dir}/freq_k{k}.hist'
       if not os.path.exists(hist_file):
            cmd = f'mkdir -p {new_dir} && cd {new_dir} && ntcard -t {ntcard_threads} -k{kmers} -p freq {r1} {r2}' 
            os.system(cmd)
    kmer = sealer_max_kmer_count(new_dir)
    return bloom_filter_calc(f'{new_dir}/freq_{kmer}.hist', FPR)


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

# check if the fasta file has gaps or not
def check_gaps(file):
    gap = 0
    for rec in SeqIO.parse(file, 'fasta'):
       seq = str(rec.seq)
       if len(seq.split('N')) > 1:
          gap +=1 
    return gap


# convert a input string containing k-mer values (e.g., "60,80,90") to the sealer parsable format (e.g., "-k60 -k80 -k90")
def k_string_converter(k_string):
  split_list = k_string.split(',')
  new_k_string = '-k' + ' -k'.join(split_list)
  return new_k_string

def check_blast_tsv(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
            if len(lines) == 1:
                print("No mitochondrial sequence found.")
                exit(0)
            else:
                print("Mitochondrial sequence(s) found")
    except FileNotFoundError:
        print(f"File '{filename}' not found.")
        exit(1)

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
     params:
        abyss_fpr=config['abyss_fpr'],
        threads=config['threads']
     
     run:
        bf = bf_abyss(input.r1, input.r2, wildcards.library, wildcards.k, params.abyss_fpr)
        shell("abyss-pe --directory {wildcards.library}/abyss v=-v kc={wildcards.kc} j={params.threads}  k={wildcards.k}  B={bf}  name=k{wildcards.k}_kc{wildcards.kc} in='{input.r1} {input.r2}' &> {log}")
       

       

#Filtering for mitochondrial sequences
rule select_length:
     input:
        rules.de_novo_assembly.output
     output:
        current_dir + "{library}/blast/k{k}_kc{kc}-scaffolds.fa"
     benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.select_length.benchmark.txt"
     shell:
        """awk '!/^>/ {{ next }} {{ getline seq }} length(seq) > 1000 {{ print $0 "\\n" seq }}' {input} > {output}"""


rule blast:
     input:
         rules.select_length.output
     output:
         current_dir + "{library}/blast/k{k}_kc{kc}-scaffolds.nt.blast.tsv"
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
            shell("mkdir -p {script_dir}/blast_db/{db_name} && cd {script_dir}/blast_db/{db_name} && makeblastdb -in {params.ref_fasta} -dbtype nucl -out {db_name}")
            shell("export BLASTDB={script_dir}/blast_db/{db_name} && blast_best-hit.py {input} {db_name} > {output}")
         else:
            shell("export BLASTDB={script_dir}/blast_db/{db_name} && blast_best-hit.py {input} {db_name} > {output}")
         check_blast_tsv(f'{output}')
         
rule create_lists:
     input:
         rules.blast.output
     output:
         ref_list=current_dir + "{library}/blast/k{k}_kc{kc}-query.txt",
         query_list=current_dir + "{library}/blast/k{k}_kc{kc}-ref.txt"
     benchmark:
         current_dir + "{library}/benchmark/k{k}_kc{kc}.create_lists.benchmark.txt"
     shell:
          "extract_tsv_value.py {input} {output.ref_list} ref ; "
          "extract_tsv_value.py {input} {output.query_list} query"
         

rule extract_seq:
     input:
         ref=rules.create_lists.output.ref_list,
         query=rules.create_lists.output.query_list,
         assemblies=rules.select_length.output
     output:
         query_out=current_dir + "{library}/mito_filtering_output/k{k}_kc{kc}-scaffolds.blast-mt_db.fa",
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
          current_dir + "{library}/prepolishing/k{k}_kc{kc}_scaffold.fa"
      params:
        workdir= current_dir + "{library}/mito_filtering_output",
        out=current_dir + "{library}/prepolishing/k{k}_kc{kc}",
        ntjoin_out=current_dir + "{library}/mito_filtering_output/k{k}_kc{kc}-scaffolds.blast-mt_db.fa.k32.w100.n1.all.scaffolds.fa",
        r1=config['r1'],
        r2=config['r2'],
        sealer_fpr=config['sealer_fpr'],
        threads=config['threads'],
        p=config['p'],
        k=config['sealer_k']
      benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.pre_polishing.benchmark.txt"
      log:
        sealer=current_dir + "{library}/prepolishing/k{k}_kc{kc}_sealer.log",
        ntjoin=current_dir + "{library}/mito_filtering_output/k{k}_kc{kc}_ntjoin.log"
      
      run: 
          import os
          target = os.path.basename(input.target)
          ref = os.path.basename(input.ref)
          log_ntjoin = os.path.basename(log.ntjoin)
          log_sealer = log.sealer
          bf = bf_sealer(params.r1, params.r2, wildcards.library, params.threads, params.sealer_fpr,params.k)
          count = sum(1 for line in open(input[0]))
          k = k_string_converter(params.k)
          if check_gaps(input.target) != 0:
            print("---Start Sealer Gap Filling---")
            if count <= 2:
               shell("""abyss-sealer -b{bf} -j {params.threads} -vv {k} -P {params.p} -o {params.out} -S {input.target} {params.r1} {params.r2} &> {log_sealer}""")
               print("no ntJoin needed")
            else:
               shell("""bash run_ntjoin.sh {params.workdir} {target} {ref} {log_ntjoin} {params.threads} && abyss-sealer -b{bf} -j {params.threads} -vv {k} -P {params.p} -o {params.out} -S {params.ntjoin_out}  {params.r1} {params.r2} &> {log_sealer}""")
               print("ntJoin needed")
          else:
             print("---No Gaps Found, Gap Filling Not Needed---")
             if count <= 2:
               shell("cp {input.target} {output}")
               print("no ntJoin needed")
             else:
                shell("""bash run_ntjoin.sh {params.workdir} {target} {ref} {log_ntjoin} {params.threads} && cp {params.ntjoin_out} {output}""")
                print("ntJoin needed")

               



#Polishing
rule bwa_alignment:
      input:
        rules.pre_polishing.output
      output:
        current_dir + "{library}/pilon/k{k}_kc{kc}.sorted.bam"
      params:
        r1=config['r1'],
        r2=config['r2'],
        threads=config['threads']
      benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.bwa_alignment.benchmark.txt"
      
      shell:
        "bwa index {input} && "
        "bwa mem {input} {params.r1} {params.r2} -t {params.threads} | samtools view -b -F 4 | samtools sort -o {output} -@ {params.threads} && "
        "samtools index {output}"

rule polishing:
      input:
        in1=rules.pre_polishing.output,
        in2=rules.bwa_alignment.output
      output:
        current_dir + "{library}/assembly_output/{library}_k{k}_kc{kc}.postsealer.postpilon.fasta"
      params:
        outdir=current_dir + "{library}/pilon/{library}_k{k}_kc{kc}.postsealer.postpilon",
        fasta=current_dir + "{library}/pilon/{library}_k{k}_kc{kc}.postsealer.postpilon.fasta",
        threads=config['threads'],
        abyss_fpr=config['abyss_fpr'],
        r1=config['r1'],
        r2=config['r2']


      benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.polishing.benchmark.txt"
      
      run:
        shell("pilon -Xmx200g --genome {input.in1} --frags {input.in2} --threads {params.threads} --output {params.outdir} --changes --fix all --verbose && mv {params.fasta} {output}")

# Standardization of strand orientation and start-site
rule end_recovery:
       input:
          rules.polishing.output
       output:
          current_dir + "{library}/end_recovery/{library}_k{k}_kc{kc}/flank_added_assembly.fa"
       params:
          r1=config['r1'],
          r2=config['r2'],
          outdir=current_dir + "{library}/end_recovery/{library}_k{k}_kc{kc}",
          sealer_fpr=config['end_recov_sealer_fpr'],
          threads=config['threads'],
          p=config['end_recov_p'],
          k=config['end_recov_sealer_k'],
          mismatch_allowed=config['mismatch_allowed']

       benchmark:
          current_dir + "{library}/benchmark/k{k}_kc{kc}.end_recovery.benchmark.txt"
       
       run:
          
          bf = bf_sealer(params.r1, params.r2, wildcards.library, params.threads, params.sealer_fpr,params.k)
          k = k_string_converter(params.k)
          shell("end_recover.py {input} {bf} {params.r1} {params.r2} {params.outdir} {params.threads} {params.p} {params.mismatch_allowed} {k}")

rule standardization:
        input:
            rules.end_recovery.output
        output:
            current_dir + "{library}/final_output/{library}_k{k}_kc{kc}/{library}_k{k}_kc{kc}.final-mtgasp-assembly.fa"
        benchmark:
            current_dir + "{library}/benchmark/k{k}_kc{kc}.standardization.benchmark.txt"
        params:
            mito_gencode=config["mt_code"],
            outdir=current_dir + "{library}/final_output/{library}_k{k}_kc{kc}",
            annotate=config["annotate"]
        run:
            if params.annotate=='No':
              shell("mtgasp_standardize.py -i {input} -c {params.mito_gencode} -o {params.outdir} -p {wildcards.library}_k{wildcards.k}_kc{wildcards.kc}")
            else:
               shell("mtgasp_standardize.py -i {input} -c {params.mito_gencode} -o {params.outdir} -p {wildcards.library}_k{wildcards.k}_kc{wildcards.kc} -a")
        
