#!/bin/bash
cd $1
ntJoin assemble target=$2 target_weight=1 references=$3  reference_weights='2' k=32 w=100 &> $4
