#!/bin/bash
# Run this script to ensure all required dependencies are installed for mtGrasp
# Make sure to edit the version for new releases


mitos_path=$1

script_directory="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "Running test. Log messages can be found in test.log."
mtgrasp.py -r1 ${script_directory}/SRR20078528_1M_subset_R1.fastq.gz -r2 ${script_directory}/SRR20078528_1M_subset_R2.fastq.gz -nsub -m 2 -o test_out -r ${script_directory}/NC_012920.1.fasta -mp ${mitos_path} &> test.log

# check if final assembly fasta is generated
file=$(find test_out/final_output/test_out_k91_kc3/ -name '*-assembly.fa' -print -quit)
if [[ -n "$file" ]]; then
    echo "mtGrasp executed successfully. All required dependencies are installed."
    echo "The output files can be found in the test_out directory. The final assembly file can be found at: test_out/final_output/test_out_k91_kc3/test_out_k91_kc3.final-mtgrasp_<version>-assembly.fa"
else
    echo "Error: mtGrasp could not run successfully due to missing dependencies. Please check the test.log file and make sure all required dependencies are installed before running the program."
fi
