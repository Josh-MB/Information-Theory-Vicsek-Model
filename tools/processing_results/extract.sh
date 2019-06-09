#!/bin/bash

for file in *.log
do
	echo "`grep 'I = ' $file | cut -d ' ' -f3`" > $file.I	
	echo "`grep 'I2 = ' $file | cut -d ' ' -f3`" >> $file.I	
	echo "`grep 'I3 = ' $file | cut -d ' ' -f3`" >> $file.I	
	echo "`grep 'I4 = ' $file | cut -d ' ' -f3`" >> $file.I	
	echo "`grep 'I5 = ' $file | cut -d ' ' -f3`" >> $file.I
done
