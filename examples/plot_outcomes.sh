#!/bin/sh
plot_outcomes.py --matched_cats=mwacs,vlssr,mrc,sumss,nvss \
	--pref_cats=nvss,sumss,mwacs,mrc,vlssr \
	--input_bayes=bayes_mwacs-v-m-s-n.txt \
	--cat_freqs=180,74,408,843,1400 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:03:00 --split=00:02:00 \
	--num_catalogues=5 --num_combinations=4 \
	--write=10
	
