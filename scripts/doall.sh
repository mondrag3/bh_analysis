#!/bin/bash

#data=~/raida/data
data=data

part=$1

echo -e "\n\E[0;49;93mReweighting\E[0m\n"
./bin/reweigh --bh=$data/Wm1j7TeV_"$part"_0.BH.root -o $data/Wm1j7TeV_"$part"_0.weights.root --pdf=CT10 --old-bh

#if false; then

echo -e "\n\E[0;49;93mMaking new histograms\E[0m\n"
./bin/hist_weights $data/Wm1j7TeV_"$part"_0.weights.root config/all_weights.bins $data/Wm1j7TeV_"$part"_0.new_hist_uncut.root

echo -e "\n\E[0;49;93mSelecting old histograms\E[0m\n"
./bin/select_old_weight_hists -i $data/Wm1j7TeV_"$part"_0.HIST.root -o $data/Wm1j7TeV_"$part"_0.old_hist_uncut.root --pdf=CT10 --uncut

echo -e "\n\E[0;49;93mDrawing comparison plots\E[0m\n"
./bin/draw_together -o cmp_"$part"_lin.pdf -i old:$data/Wm1j7TeV_"$part"_0.old_hist_uncut.root -i new:$data/Wm1j7TeV_"$part"_0.new_hist_uncut.root --name-title

./bin/draw_together -o cmp_"$part"_logx.pdf --logx -i old:$data/Wm1j7TeV_"$part"_0.old_hist_uncut.root -i new:$data/Wm1j7TeV_"$part"_0.new_hist_uncut.root --name-title

#fi
