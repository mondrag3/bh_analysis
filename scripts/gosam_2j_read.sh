#!/bin/bash

dir=/msu/data/t3work5/ivanp/ntuple_analysis/gosam/2j

part=$1

./bin/gosam_2j \
  --bh=$dir/bh/$part.root --sj=$dir/sj/$part.root --wt=$dir/wt/$part.root \
  -o $dir/$part.root \
  -w Fac0.5Ht_Ren0.5Ht_PDFCT10nnlo_cent \
  -s config/hists.css # -n 50:100
