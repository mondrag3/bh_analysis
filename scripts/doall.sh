#!/bin/bash

data=~/raida/data

./bin/reweigh --bh=$data/Wm1j7TeV_born_0.BH.root -o $data/Wm1j7TeV_born_0.weights.root --pdf=CT10 --old-bh

./bin/hist_weights $data/Wm1j7TeV_born_0.weights.root config/all_weights.bins $data/Wm1j7TeV_born_0.new_hist_uncut.root

./bin/select_old_weight_hists -i $data/Wm1j7TeV_born_0.HIST.root -o $data/Wm1j7TeV_born_0.old_hist_uncut.root

./bin/draw_together -o cmp_lin.pdf -i old:$data/Wm1j7TeV_born_0.old_hist_uncut.root -i new:$data/Wm1j7TeV_born_0.new_hist_uncut.root

./bin/draw_together -o cmp_logx.pdf --logx -i old:$data/Wm1j7TeV_born_0.old_hist_uncut.root -i new:$data/Wm1j7TeV_born_0.new_hist_uncut.root
