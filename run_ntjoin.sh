#!/bin/bash
cd $1
ntJoin assemble target=$2 target_weight=1 reference_config=$3  t=$5 k=32 w=100 &> $4
