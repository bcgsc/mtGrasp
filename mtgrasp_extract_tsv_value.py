#!/usr/bin/env python3

import pandas as pd
import sys
         


input = sys.argv[1]
output = sys.argv[2]
columnname = sys.argv[3]



df= pd.read_csv(input, sep='\t')
list = df[columnname].to_list()
with open(output, 'w') as f:
    for item in list:
        f.write("%s\n" % item)




