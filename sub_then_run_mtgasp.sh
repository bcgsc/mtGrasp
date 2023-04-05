#!/bin/bash

out_dir=$1
r1=$2
r2=$3
subsample=$4
read1_base=$5
read2_base=$6
script_dir=$7
threads=$8
mt_gen=$9
kmer=${10}
kc=${11}
ref_path=${12}
abyss_fpr=${13}
sealer_fpr=${14}
p=${15}
sealer_k=${16}
end_recov_sealer_fpr=${17}
end_recov_p=${18}
end_recov_sealer_k=${19}
mismatch_allowed=${20}


current_dir="$(pwd)/"


# check if the subsampled reads are already present
if [ -f ${out_dir}/subsets/${read1_base}_${subsample}.fastq.gz ] && [ -f ${out_dir}/subsets/${read2_base}_${subsample}.fastq.gz ]; then
    echo "Subsampled reads already present in ${out_dir}/subsets directory, skipping subsampling"
else
    echo "Subsampling ${subsample} read pairs and saving to ${out_dir}/subsets directory"
    mkdir -p ${out_dir}/subsets
    seqtk sample -s100 ${r1} ${subsample} | gzip > ${out_dir}/subsets/${read1_base}_${subsample}.fastq.gz &
    seqtk sample -s100 ${r2} ${subsample} | gzip > ${out_dir}/subsets/${read2_base}_${subsample}.fastq.gz &
    wait
fi

snakemake -s ${script_dir}/mtgasp.smk --cores ${threads} -p -k \
--config r1=${current_dir}${out_dir}/subsets/${read1_base}_${subsample}.fastq.gz r2=${current_dir}${out_dir}/subsets/${read2_base}_${subsample}.fastq.gz out_dir=${out_dir} mt_code=${mt_gen} k=${kmer} kc=${kc} ref_path=${ref_path} threads=${threads} \
abyss_fpr=${abyss_fpr} sealer_fpr=${sealer_fpr} p=${p} sealer_k=${sealer_k} end_recov_sealer_fpr=${end_recov_sealer_fpr} \
end_recov_p=${end_recov_p} end_recov_sealer_k=${end_recov_sealer_k} mismatch_allowed=${mismatch_allowed}

rm ${out_dir}/subsets/*${subsample}.fastq.gz
