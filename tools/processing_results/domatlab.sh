#!/bin/bash

# Will run the appropriate matlab scripts to generate the entropy graph
# 
# Vicsek folder location. i.e. ~/path/to/folder/vicsek/ should exist
vicsekDir=~/host/refactor_tests/

# matlab path
matlab=/usr/local/MATLAB/R2014a/bin/matlab
	
# Location to save the output files to. Relative to script
outputdir=results/
mkdir -p $outputdir
# Processing method

for metric in mibin tebin gtebin; do
	# Output figure name
	filename=$metric\_N1000_test

	# Location of the bin files. Note this will be relative to the vicsek folder.
	# i.e ~/path/to/folder/vicsek/xy_run/gtebin
	inputdir=out_test_2/$metric

	matlab -nodesktop -nodisplay  -nosplash -r "setenv DATADIR '$vicsekDir'; [fname, E, R] =  get_files('$inputdir'); display_data(fname, E, R, '$filename');exit"

	mv $filename.png $outputdir/$filename.png
	mv $filename.fig $outputdir/$filename.fig
done