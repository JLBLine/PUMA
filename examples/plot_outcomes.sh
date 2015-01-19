#!/bin/sh
plot_outcomes.py --matched_cats=mwacs,vlssr,mrc,sumss \
	--pref_cats=sumss,mwacs,mrc,vlssr \
	--input_bayes=bayes_mw-v-m-s.txt \
	--cat_freqs=180,74,408,843 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:03:00 \
	--split=00:01:15 --num_matches=4 #--accept --num_matches=3 #--reject