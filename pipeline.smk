import os.path
current_dir = os.getcwd() + "/"

rule all:
     input:
        expand(current_dir + "{library}/standardized_output/{library}_k{k}_kc{kc}/post_standardization.fasta", library = config["library"], k = config["k_sweeps"], kc = config["kc_sweeps"])

#De novo assembly
rule de_novo_assembly:
     input:
        r1=current_dir + "{library}/r1.fq.gz",
        r2=current_dir + "{library}/r2.fq.gz"
     output:
        current_dir + "{library}/abyss/k{k}_kc{kc}-scaffolds.fa"
     benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.de_novo_assembly.benchmark.txt"
     log:
        current_dir + "{library}/abyss/k{k}_kc{kc}.log"
     params:
        B=config["B"]
     shell:
        "abyss-pe --directory {wildcards.library}/abyss v=-v kc={wildcards.kc}  k={wildcards.k}  B={params.B}  name=k{wildcards.k}_kc{wildcards.kc} in='{input.r1} {input.r2}' &> {log}"
       

       

#Filtering for mitochondrial sequences
rule select_length:
     input:
        rules.de_novo_assembly.output
     output:
        current_dir + "{library}/blast/k{k}_kc{kc}-scaffolds.1000-20000bp.fa"
     benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.select_length.benchmark.txt"
     shell:
        "reformat.sh in={input} out={output} minlength=1000 maxlength=20000"


rule blast:
     input:
         rules.select_length.output
     output:
         current_dir + "{library}/blast/k{k}_kc{kc}-scaffolds.1000-20000bp.nt.blast.tsv"
     benchmark:
         current_dir + "{library}/benchmark/k{k}_kc{kc}.blast.benchmark.txt"
     params:
         db_path=config["blastdb_path"],
         db_name=config["blastdb_name"]
     shell:
         "export BLASTDB={params.db_path} && " 
         "python3 bin/blast_best-hit.py {input} {params.db_name} > {output} "

rule create_lists:
     input:
         rules.blast.output
     output:
         ref_list=current_dir + "{library}/blast/k{k}_kc{kc}-query.txt",
         query_list=current_dir + "{library}/blast/k{k}_kc{kc}-ref.txt"
     benchmark:
         current_dir + "{library}/benchmark/k{k}_kc{kc}.create_lists.benchmark.txt"
     shell:
          "python3 bin/extract_tsv_value.py -in {input} -out {output.ref_list} -c ref ; "
          "python3 bin/extract_tsv_value.py -in {input} -out {output.query_list} -c query"
         

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
          ref=rules.extract_seq.output.ref_out,
          r1=current_dir + "{library}/r1.fq.gz",
          r2=current_dir + "{library}/r2.fq.gz"
      output:
          current_dir + "{library}/sealer/k{k}_kc{kc}.postsealer_scaffold.fa"
      params:
        b=config["b"],
        p=config["P"],
        workdir= current_dir + "{library}/mito_filtering_output",
        out="{library}/sealer/k{k}_kc{kc}.postsealer",
        select_length="{library}/sealer/k{k}_kc{kc}.scaffold_10000-20000.fa"
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

          count = sum(1 for line in open(input[0]))
          if count <= 2:
              shell("reformat.sh in={input.target} out={params.select_length} minlength=10000 maxlength=20000 && abyss-sealer -b{params.b} -k 60 -k 80 -k 100 -k 120 -P {params.p} -o {params.out} -j24 -S {params.select_length} {input.r1} {input.r2} &> {log_sealer}")
              print("no ntJoin needed")
          else:
              shell("bash bin/run_ntjoin.sh {params.workdir} {target} {ref} {log_ntjoin} && reformat.sh in={wildcards.library}/mito_filtering_output/k{wildcards.k}_kc{wildcards.kc}-scaffolds.1000-20000bp.blast-mt_db.fa.k32.w500.n1.all.scaffolds.fa out={params.select_length} minlength=10000 maxlength=20000 && abyss-sealer -b{params.b} -k 60 -k 80 -k 100 -k 120 -P {params.p} -o {params.out} -S {params.select_length}  {input.r1} {input.r2} &> {log_sealer}")
              print("ntJoin needed")



#Polishing
rule bwa_alignment:
      input:
        in1=rules.pre_polishing.output,
        in2=current_dir + "{library}/r1.fq.gz",
        in3=current_dir + "{library}/r2.fq.gz"
      output:
        current_dir + "{library}/pilon/k{k}_kc{kc}.sorted.bam"
      benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.bwa_alignment.benchmark.txt"
      shell:
        "bwa index {input.in1} && "
        "bwa mem {input.in1} {input.in2} {input.in3} | samtools view -b -F 4 | samtools sort -o {output} && "
        "samtools index {output}"

rule polishing:
      input:
        in1=rules.pre_polishing.output,
        in2=rules.bwa_alignment.output
      output:
        current_dir + "{library}/assembly_output/{library}_k{k}_kc{kc}.postsealer.postpilon.fasta"
      benchmark:
        current_dir + "{library}/benchmark/k{k}_kc{kc}.polishing.benchmark.txt"
      shell:
        "pilon --genome {input.in1} --frags {input.in2} --output {wildcards.library}/pilon/{wildcards.library}_k{wildcards.k}_kc{wildcards.kc}.postsealer.postpilon --changes --fix all --verbose && "
        "mv {wildcards.library}/pilon/{wildcards.library}_k{wildcards.k}_kc{wildcards.kc}.postsealer.postpilon.fasta {output}"

rule end_recovery:
       input:
          rules.polishing.output
       output:
          current_dir + "{library}/end_recovery/{library}_k{k}_kc{kc}/flank_added_assembly.fa"
       params:
          bf=config["bf_end_recovery"],
          p=config["P_end_recovery"]
       benchmark:
          current_dir + "{library}/benchmark/k{k}_kc{kc}.end_recovery.benchmark.txt"
       shell:
          "python3 bin/end_recover.py {input} {params.bf} {params.p}"

rule standardization:
        input:
            rules.end_recovery.output
        output:
            current_dir + "{library}/standardized_output/{library}_k{k}_kc{kc}/post_standardization.fasta"
        benchmark:
            current_dir + "{library}/benchmark/k{k}_kc{kc}.standardization.benchmark.txt"
        params:
            mito_gencode=config["mito_gencode"]
        shell:
            "python3 bin/mitos_annotation.py {input} {params.mito_gencode} {wildcards.library}/standardized_output/{wildcards.library}_k{wildcards.k}_kc{wildcards.kc}"
