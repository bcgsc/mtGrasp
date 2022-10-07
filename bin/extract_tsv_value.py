#!/usr/bin/env python3
import argparse
import pandas as pd
         
parser = argparse.ArgumentParser()
parser.add_argument('-in','--input',help='Input TSV file',required=True)
parser.add_argument('-out','--output',help='Output TSV file',required=True)
parser.add_argument('-c','--columnname',help='Column name',required=True)
args = parser.parse_args()

input = args.input
output = args.output
columnname = args.columnname


def extract_tsv_value(input,output,columnname):
    df= pd.read_csv(input, sep='\t')
    list = df[columnname].to_list()
    with open(output, 'w') as f:
        for item in list:
            f.write("%s\n" % item)


if __name__ == '__main__':
    print(extract_tsv_value(input,output,columnname))

