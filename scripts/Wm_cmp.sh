#!/bin/bash

#data=~/raida/data
data=/msu/data/t3work5/ivanp/ntuple_analysis/Wm1j7TeV_cp2/data

part=$1

if false; then

echo -e "\n\E[0;49;93mReweighting\E[0m\n"
./bin/reweigh --bh=$data/Wm1j7TeV_"$part"_0.BH.root -o $data/Wm1j7TeV_"$part"_0.weights.root --pdf=CT10 --old-bh

echo -e "\n\E[0;49;93mMaking new histograms\E[0m\n"
./bin/hist_weights $data/Wm1j7TeV_"$part"_0.weights.root config/all_weights.bins $data/Wm1j7TeV_"$part"_0.new_hist_uncut.root

echo -e "\n\E[0;49;93mSelecting old histograms\E[0m\n"
./bin/select_old_weight_hists -i $data/Wm1j7TeV_"$part"_0.HIST.root -o $data/Wm1j7TeV_"$part"_0.old_hist_uncut.root --pdf=CT10 --uncut

fi

echo -e "\n\E[0;49;93mDrawing comparison plots\E[0m\n"
#./bin/draw_together -o out/cmp_"$part"_lin.pdf -i old:$data/Wm1j7TeV_"$part"_0.old_hist_uncut.root -i new:$data/Wm1j7TeV_"$part"_0.new_hist_uncut.root --name-title

#./bin/draw_together -o out/cmp_"$part"_logx.pdf --logx -i old:$data/Wm1j7TeV_"$part"_0.old_hist_uncut.root -i new:$data/Wm1j7TeV_"$part"_0.new_hist_uncut.root --name-title

./bin/draw_together -o $data/../cmp_"$part"_logy.pdf --logy -i old:$data/Wm1j7TeV_"$part"_0.old_hist_uncut.root -i new:$data/Wm1j7TeV_"$part"_0.new_hist_uncut.root --name-title

#./bin/draw_together -o out/cmp_"$part"_logxy.pdf --logx --logy -i old:$data/Wm1j7TeV_"$part"_0.old_hist_uncut.root -i new:$data/Wm1j7TeV_"$part"_0.new_hist_uncut.root --name-title

