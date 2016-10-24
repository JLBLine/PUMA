# !/bin/sh
START_TIME=$SECONDS
cross_match.py --separation=180 \
 	--table1=$PUMA_DIR/original_cats/vizier_mwacs.fits \
 	--details1=MWACS,RAJ2000,e_RAJ2000,DEJ2000,e_DEJ2000,180,S180,e_S180,PA,Maj,Min,-,- \
 	--units1=deg,arcmin,deg,arcmin,Jy,Jy,deg,arcmin,arcmin \
 	--ra_lims1=0,100 \
 	--dec_lims1=-55,-20 \
 	--prefix1=mwacs \
 	--table2=$PUMA_DIR/original_cats/vlssr_names.fits \
 	--details2=Name,RA_J2000,RA_err,DEC_J2000,DEC_err,74,Flux,Flux_err,PA,major,minor,-,- \
 	--units2=deg,deg,deg,deg,Jy,Jy,deg,arcsec,arcsec \
 	--ra_lims2=0,140 \
 	--dec_lims2=-20,40 \
 	--prefix2=vlssr
STRING_1="Time to match vlssr: "
TIME_1=$SECONDS
ELAPSED_1=$(($TIME_1 - $START_TIME))
PRINT_1=$STRING_1$ELAPSED_1
echo $PRINT_1
##-------------------------------------------------------------------------------------------------
cross_match.py --separation=180 \
	--table1=$PUMA_DIR/original_cats/vizier_mwacs.fits \
	--details1=MWACS,RAJ2000,e_RAJ2000,DEJ2000,e_DEJ2000,180,S180,e_S180,PA,Maj,Min,-,- \
	--units1=deg,arcmin,deg,arcmin,Jy,Jy,deg,arcmin,arcmin \
	--ra_lims1=0,100 \
	--dec_lims1=-55,-20 \
	--prefix1=mwacs \
	--table2=$PUMA_DIR/original_cats/vizier_mrc.fits \
	--details2=MRC,_RAJ2000,e_RA2000,_DEJ2000,e_DE2000,408,S408,e_S408,-,-,-,-,- \
	--units2=deg,sec,deg,arcsec,Jy,Jy,-,-,- \
	--ra_lims2=150,250 \
	--dec_lims2=-30,10 \
	--prefix2=mrc
STRING_2="Time to match mrc: "
TIME_2=$SECONDS
ELAPSED_2=$(($TIME_2 - $TIME_1))
PRINT_2=$STRING_2$ELAPSED_2
echo $PRINT_2
##-------------------------------------------------------------------------------------------------
cross_match.py --separation=180 \
	--table1=$PUMA_DIR/original_cats/vizier_mwacs.fits \
	--details1=MWACS,RAJ2000,e_RAJ2000,DEJ2000,e_DEJ2000,180,S180,e_S180,PA,Maj,Min,-,- \
	--units1=deg,arcmin,deg,arcmin,Jy,Jy,deg,arcmin,arcmin \
	--ra_lims1=0,100 \
	--dec_lims1=-55,-20 \
	--prefix1=mwacs \
	--table2=$PUMA_DIR/original_cats/sumss_names.fits \
	--details2=name,_RAJ2000,e_RAJ2000,_DEJ2000,e_DEJ2000,843,St,e_St,PA,MajAxis,MinAxis,-,- \
	--units2=deg,arcsec,deg,arcsec,mJy,mJy,deg,arcsec,arcsec \
	--ra_lims2=0,100 \
	--dec_lims2=-60,-30 \
	--prefix2=sumss
STRING_3="Time to match sumss: "
TIME_3=$SECONDS
ELAPSED_3=$(($TIME_3 - $TIME_2))
PRINT_3=$STRING_3$ELAPSED_3
echo $PRINT_3
##-------------------------------------------------------------------------------------------------
cross_match.py --separation=180 \
	--table1=$PUMA_DIR/original_cats/vizier_mwacs.fits \
	--details1=MWACS,RAJ2000,e_RAJ2000,DEJ2000,e_DEJ2000,180,S180,e_S180,PA,Maj,Min,-,- \
	--units1=deg,arcmin,deg,arcmin,Jy,Jy,deg,arcmin,arcmin \
	--ra_lims1=0,100 \
	--dec_lims1=-55,-20 \
	--prefix1=mwacs \
	--table2=$PUMA_DIR/original_cats/vizier_nvss.fits \
	--details2=NVSS,_RAJ2000,e_RAJ2000,_DEJ2000,e_DEJ2000,1400,S1_4,e_S1_4,PA,MajAxis,MinAxis,-,- \
	--units2=deg,sec,deg,arcsec,mJy,mJy,deg,arcsec,arcsec \
	--ra_lims2=100,200 \
	--dec_lims2=-40,0 \
	--prefix2=nvss
STRING_4="Time to match nvss: "
TIME_4=$SECONDS
ELAPSED_4=$(($TIME_4 - $TIME_3))
PRINT_4=$STRING_4$ELAPSED_4
echo $PRINT_4 
##-------------------------------------------------------------------------------------------------
calculate_bayes.py --primary_cat=mwacs --matched_cats=vlssr,mrc,sumss,nvss \
	--primary_freq=180 --matched_freqs=74,408,843,1400 \
	--out_name=bayes_mwacs-v-m-s-n.txt --resolution=00:03:00
STRING_5="Time to calculate posterior probabilities: "
TIME_5=$SECONDS
ELAPSED_5=$(($TIME_5 - $TIME_4))
PRINT_5=$STRING_5$ELAPSED_5
echo $PRINT_5
##-------------------------------------------------------------------------------------------------
make_table.py --matched_cats=mwacs,vlssr,mrc,sumss,nvss \
	--pref_cats=nvss,sumss,mwacs,mrc,vlssr \
	--input_bayes=bayes_mwacs-v-m-s-n.txt \
	--cat_freqs=180,74,408,843,1400 \
	--prob_thresh=0.8,0.95 \
	--epsilon_thresh=0.1 --chi_thresh=10 \
	--resolution=00:03:00 \
	--output_name=puma_mwacs-v-m-s-n_split --verbose --split=00:02:00 \
	--big_cat --prefix=PUMA --format=fits
STRING_6="Time to apply criteria and create final matched table: "
TIME_6=$SECONDS
ELAPSED_6=$(($TIME_6 - $TIME_5))
PRINT_6=$STRING_6$ELAPSED_6
##-------------------------------------------------------------------------------------------------
STRING_7="Total Time: "
TIME_7=$SECONDS
ELAPSED_7=$(($TIME_7 - $START_TIME))
PRINT_7=$STRING_7$ELAPSED_7
echo $PRINT_1
echo $PRINT_2
echo $PRINT_3
echo $PRINT_4
echo $PRINT_5
echo $PRINT_6
echo $PRINT_7
