#!/bin/bash

#extract="~/scripts/extract.sh"
#extract="~/host/scripts/extract.sh"
extract="~/host/fixed_scripts/processing_results/extract.sh"

tar -zxvf *.tar.gz
#tar -xvf *.tar
mkdir mibin
mv mibin_* mibin
mkdir tebin
mv tebin_* tebin
mkdir gtebin
mv gtebin_* gtebin
mkdir params
mv params_* params
cd mibin
eval $extract
cd ../tebin
eval $extract
cd ../gtebin
eval $extract
cd ../params
eval $extract
cd ..

