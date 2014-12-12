#!/bin/bash

dir=/msu/data/t3work5/ivanp/ntuple_analysis/gosam/H2j

part=$1

./bin/gosam_2j \
  --bh=$dir/bh/$part.root --sj=$dir/sj/$part.root --wt=$dir/wt/$part.root \
  -o $dir/hist/$part.root \
  -w Fac0.5Ht_Ren0.5Ht_PDFCT10nnlo_cent -w Fac0.5Ht_Ren0.25Ht_PDFCT10nnlo_cent -w Fac0.5Ht_Ren1Ht_PDFCT10nnlo_cent \
  -w Fac0.25Ht_Ren0.5Ht_PDFCT10nnlo_cent -w Fac0.25Ht_Ren0.25Ht_PDFCT10nnlo_cent \
  -w Fac1Ht_Ren0.5Ht_PDFCT10nnlo_cent -w Fac1Ht_Ren1Ht_PDFCT10nnlo_cent \
  -s config/gosam_2j.css # -n 50:100
