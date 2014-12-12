#!/bin/bash

dir=/msu/data/t3work5/ivanp/ntuple_analysis/gosam/$1

part=$2

./bin/gosam_hists \
  --bh=$dir/bh/$part.root --sj=$dir/sj/$part.root --wt=$dir/wt/$part.root \
  -o $dir/hist/$part.root \
  -j AntiKt4 -s config/gosam_2j.css # -n 50:100
