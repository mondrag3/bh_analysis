#!/bin/bash

dir=data

part=$1

./bin/gosam_2j \
  --bh=$dir/bh_"$part"_10k.root \
  --sj=$dir/sj_"$part"_10k.root \
  --wt=$dir/wt_"$part"_10k.root \
  -o out/hist_"$part"_10k.root \
  -w Fac0.5Ht_Ren0.5Ht_PDFCT10nnlo_cent -w Fac0.25Ht_Ren0.25Ht_PDFCT10nnlo_cent \
  -j AntiKt4 -j AntiKt5 -j AntiKt6 \
  -s config/gosam_2j.css # -n 50:100
