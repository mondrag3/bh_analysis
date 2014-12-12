#!/bin/bash

dir=/msu/data/t3work5/ivanp/ntuple_analysis/gosam/$1

part=$2

echo -e "\n\E[0;49;93mReweighting $part\E[0m\n"
./bin/reweigh --bh=$dir/bh/"$part".root -o $dir/wt/"$part".root --pdf=CT10nnlo \
  1> $dir/wt/"$part".out 2> $dir/wt/"$part".err
