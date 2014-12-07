#!/bin/bash

dir=data

part=$1

./bin/read_withFriending \
  --bh=$dir/bh_"$part"_10k.root \
  --sj=$dir/sj_"$part"_10k.root \
  --wt=$dir/wt_"$part"_10k.root \
  -o out/hist_"$part"_10k.root \
  -w Fac0.5Ht_Ren0.5Ht_PDFCT10nnlo_cent \
  -s config/gosam_2j.css # -n 50:100
